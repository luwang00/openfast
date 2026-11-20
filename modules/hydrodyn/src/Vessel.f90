!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Laboratory of the Rockies
!
!    This file is part of HydroDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
!> MMG-type ship maneuvering submodule of HydroDyn (hull, propellers with rudders, speed and steering controllers).
!! One instance per potential-flow body with VesselMod=1. MMG formulas are evaluated in MMG axes (x forward, y starboard,
!! z down), obtained from the OpenFAST body axes by T = diag(1,-1,-1).
!!
!! Continuous-state vector x%x (size p%NStates), in this order:
!!   x(1:NProp)                      thrust T_i of each propeller (N)
!!   x(NProp+1:NProp+NRud)           rudder angle delta (rad, MMG sign) of each propeller with RudFlag=1, in propeller order;
!!                                   p%iRudState(i) is the state index of the rudder of propeller i (p%iRudState(i) = 0 if RudFlag(i) = 0)
!!   x(iXi+1), x(iXi+2), x(iXi+3)    integrators xi_x (m-s), xi_y (m-s), xi_psi (rad-s), with iXi = NProp+NRud
!!   x(iFlt+1:iFlt+6)                only when CtrlWaveFc > 0, with iFlt = iXi+3: wave-filtered controller inputs
!!                                   x error (m), x velocity (m/s), y error (m), y velocity (m/s), heading (rad), yaw rate (rad/s)
module Vessel

   use NWTC_Library
   use Vessel_Types
   use SeaSt_WaveField

   implicit none

   private

   integer(IntKi), parameter :: NFltSig = 6   !< Wave-filtered controller inputs: x error, x velocity, y error, y velocity, heading, yaw rate

   !> Instantaneous kinematics of one Vessel instance
   type :: VslKinType
      real(R8Ki) :: Rg2b(3,3)                 !< Global to body rotation
      real(R8Ki) :: PosO(3)                   !< Global position of the body reference point
      real(R8Ki) :: VelO(3)                   !< Global velocity of the body reference point
      real(R8Ki) :: Omega(3)                  !< Global angular velocity
      real(ReKi) :: u, vm, r                  !< MMG surge, sway (relative to the mean current) and yaw rate at the hull force center
      real(ReKi) :: UTot, Beta, vp, rp        !< Speed, drift angle, v', r'
      real(ReKi) :: uH, vmH                   !< MMG surge and sway at the hull reference point
      real(ReKi) :: Sig(NFltSig)              !< Raw controller inputs
   end type VslKinType

   public :: Vessel_Init
   public :: Vessel_CalcOutput
   public :: Vessel_CalcContStateDeriv
   public :: Vessel_UpdateStates

contains

!----------------------------------------------------------------------------------------------------------------------------------
!> Reads the Vessel input file, sets the parameters, creates the input mesh and initializes the states.
subroutine Vessel_Init( InitInp, u, p, x, m, ErrStat, ErrMsg )
   type(Vessel_InitInputType),       intent(in   ) :: InitInp
   type(Vessel_InputType),           intent(  out) :: u
   type(Vessel_ParameterType),       intent(  out) :: p
   type(Vessel_ContinuousStateType), intent(  out) :: x
   type(Vessel_MiscVarType),         intent(  out) :: m
   integer(IntKi),                   intent(  out) :: ErrStat
   character(*),                     intent(  out) :: ErrMsg

   real(R8Ki)                 :: R0(3,3), PH0(3)
   real(ReKi)                 :: qNd, mxA0, myA0
   logical                    :: mxDef, myDef
   integer(IntKi)             :: i
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   character(*), parameter    :: RoutineName = 'Vessel_Init'

   ErrStat = ErrID_None
   ErrMsg  = ''

   p%BodyID    = InitInp%BodyID
   p%DT        = InitInp%DT
   p%FrameSpd  = InitInp%FrameSpd
   p%WaveField => InitInp%WaveField
   p%RhoW      = p%WaveField%WtrDens
   p%HasWaves  = p%WaveField%WaveMod /= WaveMod_None

   call ReadVesselFile( InitInp, p, mxDef, myDef, ErrStat2, ErrMsg2 ); if (Failed()) return

      ! MMG added masses: DEFAULT from the zero-frequency added mass of this body
   qNd  = 0.5_ReKi*p%RhoW*p%Lpp**2*p%Draft
   mxA0 = InitInp%A0Surge/qNd
   myA0 = InitInp%A0Sway /qNd
   if ( mxDef ) then
      p%mx = InitInp%A0Surge
   else
      if ( abs(p%mx - mxA0) > 0.25_ReKi*abs(mxA0) ) call SetErrStat( ErrID_Warn, 'Body '//trim(Num2LStr(p%BodyID))// &
         ': mx differs by more than 25% from the zero-frequency surge added mass ('//trim(Num2LStr(mxA0))//'); acceptable only if it matches the hull-derivative data set.', ErrStat, ErrMsg, RoutineName )
      p%mx = p%mx*qNd
   end if
   if ( myDef ) then
      p%my = InitInp%A0Sway
   else
      if ( abs(p%my - myA0) > 0.25_ReKi*abs(myA0) ) call SetErrStat( ErrID_Warn, 'Body '//trim(Num2LStr(p%BodyID))// &
         ': my differs by more than 25% from the zero-frequency sway added mass ('//trim(Num2LStr(myA0))//'); acceptable only if it matches the hull-derivative data set.', ErrStat, ErrMsg, RoutineName )
      p%my = p%my*qNd
   end if

   if ( .not. EqualRealNos( InitInp%PtfmRefztRot, 0.0_R8Ki ) ) call SetErrStat( ErrID_Warn, 'Body '//trim(Num2LStr(p%BodyID))// &
      ': the Vessel hull axes are the HydroDyn body axes, but PtfmRefztRot is nonzero.', ErrStat, ErrMsg, RoutineName )

      ! Initial position of the hull reference point (reference of the speed and track controllers)
   R0  = EulerConstructZYX( real(InitInp%PlatformPos(4:6),R8Ki) )
   PH0 = real(InitInp%PlatformPos(1:3),R8Ki) + matmul( transpose(R0), real(InitInp%BodyRefPos + p%HullRef,R8Ki) )
   p%XYH0 = real( PH0(1:2), ReKi )

   call CheckWaveGrid( InitInp, R0, ErrStat2, ErrMsg2 ); if (Failed()) return

      ! States: thrust per propeller, rudder angles, integrators, wave filters
   p%NStates = p%NProp + p%NRud + 3
   p%iXi     = p%NProp + p%NRud
   p%iFlt    = p%iXi + 3
   if ( p%CtrlWaveFc > 0.0_ReKi ) p%NStates = p%NStates + NFltSig
   call AllocAry( x%x, p%NStates, 'x%x', ErrStat2, ErrMsg2 ); if (Failed()) return
   x%x = 0.0_R8Ki
   do i = 1, p%NProp
      x%x(i) = min( max( p%ThrShr(i)*p%T0, p%TMin ), p%TMax )
      if ( p%iRudState(i) > 0 .and. p%AutoPilotMod == 0 ) x%x(p%iRudState(i)) = p%Delta0
   end do
   if ( p%CtrlWaveFc > 0.0_ReKi ) x%x(p%iFlt+5) = atan2( R0(1,2), R0(1,1) )

   call InitMesh( InitInp, u, ErrStat2, ErrMsg2 ); if (Failed()) return
   call InitMisc( p, m, ErrStat2, ErrMsg2 ); if (Failed()) return

   if ( InitInp%UnSum > 0 ) call WriteSummary( InitInp%UnSum, p, mxA0, myA0 )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed

   subroutine CheckWaveGrid( InitInp, R0, ErrStat, ErrMsg )
      type(Vessel_InitInputType), intent(in   ) :: InitInp
      real(R8Ki),                 intent(in   ) :: R0(3,3)
      integer(IntKi),             intent(  out) :: ErrStat
      character(*),               intent(  out) :: ErrMsg
      real(R8Ki)     :: Pos(3)
      integer(IntKi) :: j

      ErrStat = ErrID_None
      ErrMsg  = ''
      if ( .not. p%HasWaves ) return
      do j = 1, p%NProp
         Pos = real(InitInp%PlatformPos(1:3),R8Ki) + matmul( transpose(R0), real(InitInp%BodyRefPos + p%PropPos(:,j),R8Ki) )
         if ( OutsideGrid( Pos ) ) then
            call SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(p%BodyID))//': propeller '//trim(Num2LStr(j))// &
               ' (plus a margin of 0.1 Lpp) lies outside the SeaState wave grid.', ErrStat, ErrMsg, 'CheckWaveGrid' )
            return
         end if
         if ( p%RudFlag(j) == 0 ) cycle
         Pos = real(InitInp%PlatformPos(1:3),R8Ki) + matmul( transpose(R0), real(InitInp%BodyRefPos + p%RudPos(:,j),R8Ki) )
         if ( OutsideGrid( Pos ) ) then
            call SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(p%BodyID))//': rudder '//trim(Num2LStr(j))// &
               ' (plus a margin of 0.1 Lpp) lies outside the SeaState wave grid.', ErrStat, ErrMsg, 'CheckWaveGrid' )
            return
         end if
      end do
   end subroutine CheckWaveGrid

   logical function OutsideGrid( Pos )
      real(R8Ki), intent(in) :: Pos(3)
      real(ReKi) :: XYMin(2), XYMax(2), Margin
      Margin = 0.1_ReKi*p%Lpp
      XYMin  = p%WaveField%VolGridParams%pZero(2:3)
      XYMax  = XYMin + real( p%WaveField%VolGridParams%n(2:3) - 1, ReKi )*p%WaveField%VolGridParams%delta(2:3)
      OutsideGrid = any( real(Pos(1:2),ReKi) - Margin < XYMin ) .or. any( real(Pos(1:2),ReKi) + Margin > XYMax ) .or. &
                    real(Pos(3),ReKi) < -p%WaveField%EffWtrDpth
   end function OutsideGrid
end subroutine Vessel_Init

!----------------------------------------------------------------------------------------------------------------------------------
!> Parses the Vessel input file and its table files into the parameters.
subroutine ReadVesselFile( InitInp, p, mxDef, myDef, ErrStat, ErrMsg )
   type(Vessel_InitInputType), intent(in   ) :: InitInp
   type(Vessel_ParameterType), intent(inout) :: p
   logical,                    intent(  out) :: mxDef, myDef
   integer(IntKi),             intent(  out) :: ErrStat
   character(*),               intent(  out) :: ErrMsg

   type(FileInfoType)          :: FileInfo
   character(1024)             :: Path, ResFile, Line, EchoFile
   character(1024), allocatable:: OWFile(:), RudTblFile(:)
   character(64)               :: FAlphaStr
   logical                     :: Echo, T0Def
   integer(IntKi)              :: CurLine, UnEc, i, ios, WaveLdMod
   real(ReKi)                  :: Row(13), RudDx, RudDz, RudLam, xHnd
   real(ReKi)                  :: Tmp
   integer(IntKi)              :: ErrStat2
   character(ErrMsgLen)        :: ErrMsg2
   character(*), parameter     :: RoutineName = 'ReadVesselFile'

   ErrStat = ErrID_None
   ErrMsg  = ''
   UnEc    = -1
   mxDef   = .false.
   myDef   = .false.

   call ProcessComFile( InitInp%InputFile, FileInfo, ErrStat2, ErrMsg2 ); if (Failed()) return
   call GetPath( InitInp%InputFile, Path )

   CurLine = 3
   call ParseVar( FileInfo, CurLine, 'Echo', Echo, ErrStat2, ErrMsg2 ); if (Failed()) return
   if ( Echo ) then
      EchoFile = trim(InitInp%OutRootName)//'.Vsl'//trim(Num2LStr(p%BodyID))//'.ech'
      call OpenEcho( UnEc, trim(EchoFile), ErrStat2, ErrMsg2 ); if (Failed()) return
      write(UnEc,'(A)') 'Echo file for HydroDyn Vessel input file: '//trim(InitInp%InputFile)
      write(UnEc,'(A)') trim(FileInfo%Lines(1))
      write(UnEc,'(A)') trim(FileInfo%Lines(2))
      CurLine = 3
      call ParseVar( FileInfo, CurLine, 'Echo', Echo, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   end if

   !---------------------- GENERAL
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'HullRefX', p%HullRef(1), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'HullRefY', p%HullRef(2), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'HullRefZ', p%HullRef(3), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   !---------------------- HULL
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'Lpp',     p%Lpp,     ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Breadth', p%Breadth, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Draft',   p%Draft,   ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'zH',      p%zH,      ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'UMin',    p%UMin,    ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'ResMod',  p%ResMod,  ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'R0',      p%R0,      ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'ResFile', ResFile,   ErrStat2, ErrMsg2, UnEc, IsPath=.true. ); if (Failed()) return
   if ( PathIsRelative(ResFile) ) ResFile = trim(Path)//trim(ResFile)

   !---------------------- HULL DERIVATIVES
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'Xvv',   p%XDer(1), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Xvr',   p%XDer(2), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Xrr',   p%XDer(3), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Xvvvv', p%XDer(4), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Yv',    p%YDer(1), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Yr',    p%YDer(2), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Yvvv',  p%YDer(3), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Yvvr',  p%YDer(4), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Yvrr',  p%YDer(5), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Yrrr',  p%YDer(6), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Nv',    p%NDer(1), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Nr',    p%NDer(2), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Nvvv',  p%NDer(3), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Nvvr',  p%NDer(4), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Nvrr',  p%NDer(5), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Nrrr',  p%NDer(6), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   !---------------------- ADDED MASS
   call SkipLine()
   call ParseRealWDefault( 'mx', p%mx, mxDef ); if (ErrStat >= AbortErrLev) return
   call ParseRealWDefault( 'my', p%my, myDef ); if (ErrStat >= AbortErrLev) return

   !---------------------- PROPELLERS
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'NProp', p%NProp, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   if ( p%NProp < 0 ) then
      call SetErrStat( ErrID_Fatal, 'NProp must not be negative.', ErrStat, ErrMsg, RoutineName ); call Cleanup(); return
   end if
   call AllocArrays(); if (ErrStat >= AbortErrLev) return
   call SkipLine(); call SkipLine()
   do i = 1, p%NProp
      call NextRow(); if (ErrStat >= AbortErrLev) return
      read(Line,*,iostat=ios) Row(1:10)
      if ( RowFailed('propeller', i) ) return
      p%PropPos(:,i) = Row(1:3)
      p%PropD(i)     = Row(4)
      p%tP(i)        = Row(5)
      p%wP0(i)       = Row(6)
      p%C1(i)        = Row(7)
      p%C2Pos(i)     = Row(8)
      p%C2Neg(i)     = Row(9)
      p%ThrShr(i)    = Row(10)
   end do

   !---------------------- PROPELLER OPEN-WATER DATA
   call SkipLine(); call SkipLine(); call SkipLine()
   do i = 1, p%NProp
      call NextRow(); if (ErrStat >= AbortErrLev) return
      read(Line,*,iostat=ios) p%PropOWMod(i), p%KTPoly(:,i), p%KQPoly(:,i), OWFile(i), p%RotDir(i), p%EtaR(i)
      if ( RowFailed('propeller open-water', i) ) return
      if ( PathIsRelative(OWFile(i)) ) OWFile(i) = trim(Path)//trim(OWFile(i))
   end do

   !---------------------- RUDDERS
   call SkipLine(); call SkipLine(); call SkipLine()
   do i = 1, p%NProp
      call NextRow(); if (ErrStat >= AbortErrLev) return
      read(Line,*,iostat=ios) p%RudFlag(i), Row(1:13), p%RudLiftMod(i), FAlphaStr, p%StallAng(i), RudTblFile(i)
      if ( RowFailed('rudder', i) ) return
      RudDx        = Row(1)
      RudDz        = Row(2)
      p%RudAR(i)   = Row(3)
      p%RudHR(i)   = Row(4)
      RudLam       = Row(5)
      p%tR(i)      = Row(6)
      p%aH(i)      = Row(7)
      xHnd         = Row(8)
      p%GamRPos(i) = Row(9)
      p%GamRNeg(i) = Row(10)
      p%lRnd(i)    = Row(11)
      p%Eps(i)     = Row(12)
      p%Kappa(i)   = Row(13)
      p%RudPos(:,i)  = p%PropPos(:,i) + (/ -RudDx, 0.0_ReKi, RudDz /)
      p%RudHPos(:,i) = (/ p%HullRef(1) + xHnd*p%Lpp, p%HullRef(2), p%RudPos(3,i) /)
      p%StallAng(i)  = p%StallAng(i)*D2R
      call Conv2UC( FAlphaStr )
      if ( index( adjustl(FAlphaStr), 'DEFAULT' ) == 1 ) then
         p%FAlpha(i) = 6.13_ReKi*RudLam/( RudLam + 2.25_ReKi )
         if ( p%RudFlag(i) == 1 .and. p%RudLiftMod(i) == 1 .and. RudLam <= 0.0_ReKi ) then
            call SetErrStat( ErrID_Fatal, 'RudLam must be positive for the default rudder lift gradient (rudder '//trim(Num2LStr(i))//').', ErrStat, ErrMsg, RoutineName )
            call Cleanup(); return
         end if
      else
         read(FAlphaStr,*,iostat=ios) p%FAlpha(i)
         if ( RowFailed('rudder FAlpha', i) ) return
      end if
      if ( PathIsRelative(RudTblFile(i)) ) RudTblFile(i) = trim(Path)//trim(RudTblFile(i))
   end do

   !---------------------- SPEED CONTROL
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'SpdCtrlMod', p%SpdCtrlMod, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseRealWDefault( 'T0', p%T0, T0Def ); if (ErrStat >= AbortErrLev) return
   call ParseVar( FileInfo, CurLine, 'KPx',  p%KPx,  ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KIx',  p%KIx,  ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KDx',  p%KDx,  ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'TauT', p%TauT, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'TMin', p%TMin, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'TMax', p%TMax, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   !---------------------- STEERING
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'AutoPilotMod', p%AutoPilotMod, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'Delta0',       p%Delta0,       ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'PsiRef',       p%PsiRef,       ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KPy',          p%KPy,          ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KIy',          p%KIy,          ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KDy',          p%KDy,          ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'PsiCmdMax',    p%PsiCmdMax,    ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KPpsi',        p%KPpsi,        ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KIpsi',        p%KIpsi,        ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'KDpsi',        p%KDpsi,        ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'TauDelta',     p%TauDelta,     ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'DeltaRateMax', p%DeltaRateMax, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'DeltaMax',     p%DeltaMax,     ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   p%Delta0       = p%Delta0      *D2R
   p%PsiRef       = p%PsiRef      *D2R
   p%PsiCmdMax    = p%PsiCmdMax   *D2R
   p%DeltaRateMax = p%DeltaRateMax*D2R
   p%DeltaMax     = p%DeltaMax    *D2R

   !---------------------- POWER AND CONTROLLER WAVE FILTER
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'EtaD',       p%EtaD,       ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
   call ParseVar( FileInfo, CurLine, 'CtrlWaveFc', p%CtrlWaveFc, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   !---------------------- WAVE MEAN LOADS (reserved)
   call SkipLine()
   call ParseVar( FileInfo, CurLine, 'WaveLdMod', WaveLdMod, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return

   !---------------------- Tables, derived data and checks
   if ( p%ResMod == 2 ) then
      call ReadResTable( ResFile, p, ErrStat2, ErrMsg2 ); if (Failed()) return
   end if
   do i = 1, p%NProp
      if ( p%PropOWMod(i) == 2 ) then
         call ReadOWTable( OWFile(i), p%OWTbl(i), ErrStat2, ErrMsg2 ); if (Failed()) return
      end if
      if ( p%RudFlag(i) == 1 .and. p%RudLiftMod(i) == 2 ) then
         call ReadRudTable( RudTblFile(i), p%RudTbl(i), ErrStat2, ErrMsg2 ); if (Failed()) return
      end if
   end do

   p%NRud = 0
   do i = 1, p%NProp
      if ( p%RudFlag(i) == 1 ) then
         p%NRud = p%NRud + 1
         p%iRudState(i) = p%NProp + p%NRud
      else
         p%iRudState(i) = 0
      end if
   end do

   if ( T0Def ) then
      Tmp = 0.0_ReKi
      do i = 1, p%NProp
         Tmp = Tmp + p%ThrShr(i)*( 1.0_ReKi - p%tP(i) )
      end do
      if ( Tmp > 0.0_ReKi ) then
         p%T0 = Resistance( p, abs(p%FrameSpd) )/Tmp
      else
         p%T0 = 0.0_ReKi
      end if
   end if

   call CheckParams( WaveLdMod ); if (ErrStat >= AbortErrLev) return

   call Cleanup()

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed

   subroutine Cleanup()
      if ( UnEc > 0 ) close(UnEc)
      UnEc = -1
      call NWTC_Library_DestroyFileInfoType( FileInfo, ErrStat2, ErrMsg2 )
   end subroutine Cleanup

   subroutine SkipLine()
      if ( UnEc > 0 .and. CurLine <= size(FileInfo%Lines) ) write(UnEc,'(A)') trim(FileInfo%Lines(CurLine))
      CurLine = CurLine + 1
   end subroutine SkipLine

   subroutine NextRow()
      if ( CurLine > size(FileInfo%Lines) ) then
         call SetErrStat( ErrID_Fatal, 'Vessel input file "'//trim(InitInp%InputFile)//'" is too short.', ErrStat, ErrMsg, RoutineName )
         call Cleanup()
         return
      end if
      Line = FileInfo%Lines(CurLine)
      if ( UnEc > 0 ) write(UnEc,'(A)') trim(Line)
      CurLine = CurLine + 1
   end subroutine NextRow

   logical function RowFailed( TableName, iRow )
      character(*),   intent(in) :: TableName
      integer(IntKi), intent(in) :: iRow
      RowFailed = ios /= 0
      if ( RowFailed ) then
         call SetErrStat( ErrID_Fatal, 'Error reading row '//trim(Num2LStr(iRow))//' of the '//TableName//' table in "'// &
                          trim(InitInp%InputFile)//'":'//NewLine//'   '//trim(Line), ErrStat, ErrMsg, RoutineName )
         call Cleanup()
      end if
   end function RowFailed

   subroutine ParseRealWDefault( Name, Var, IsDefault )
      character(*), intent(in   ) :: Name
      real(ReKi),   intent(  out) :: Var
      logical,      intent(  out) :: IsDefault
      character(64) :: Str
      call ParseVar( FileInfo, CurLine, Name, Str, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      call Conv2UC( Str )
      IsDefault = index( adjustl(Str), 'DEFAULT' ) == 1
      Var = 0.0_ReKi
      if ( .not. IsDefault ) then
         read(Str,*,iostat=ios) Var
         if ( ios /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Invalid value of '//Name//' in "'//trim(InitInp%InputFile)//'": "'//trim(Str)// &
                             '" (number or DEFAULT expected).', ErrStat, ErrMsg, RoutineName )
            call Cleanup()
         end if
      end if
   end subroutine ParseRealWDefault

   subroutine AllocArrays()
      integer(IntKi) :: n
      n = p%NProp
      call AllocAry( p%PropPos,   3, n, 'PropPos',   ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%PropD,        n, 'PropD',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%tP,           n, 'tP',        ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%wP0,          n, 'wP0',       ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%C1,           n, 'C1',        ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%C2Pos,        n, 'C2Pos',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%C2Neg,        n, 'C2Neg',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%ThrShr,       n, 'ThrShr',    ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%PropOWMod,    n, 'PropOWMod', ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%KTPoly,    3, n, 'KTPoly',    ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%KQPoly,    3, n, 'KQPoly',    ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%HasKQ,        n, 'HasKQ',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%J0,           n, 'J0',        ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RotDir,       n, 'RotDir',    ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%EtaR,         n, 'EtaR',      ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RudFlag,      n, 'RudFlag',   ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%iRudState,    n, 'iRudState', ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RudPos,    3, n, 'RudPos',    ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RudHPos,   3, n, 'RudHPos',   ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RudAR,        n, 'RudAR',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RudHR,        n, 'RudHR',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%tR,           n, 'tR',        ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%aH,           n, 'aH',        ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%GamRPos,      n, 'GamRPos',   ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%GamRNeg,      n, 'GamRNeg',   ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%lRnd,         n, 'lRnd',      ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%Eps,          n, 'Eps',       ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%Kappa,        n, 'Kappa',     ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%RudLiftMod,   n, 'RudLiftMod',ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%FAlpha,       n, 'FAlpha',    ErrStat2, ErrMsg2 ); if (Failed()) return
      call AllocAry( p%StallAng,     n, 'StallAng',  ErrStat2, ErrMsg2 ); if (Failed()) return
      allocate( p%OWTbl(n), p%RudTbl(n), OWFile(n), RudTblFile(n), stat=ErrStat2 )
      if ( ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error allocating the propeller and rudder tables.', ErrStat, ErrMsg, RoutineName )
         call Cleanup()
      end if
   end subroutine AllocArrays

   subroutine CheckParams( WaveLdMod )
      integer(IntKi), intent(in) :: WaveLdMod
      integer(IntKi) :: j
      real(ReKi)     :: k0, k1, k2, Disc

      if ( p%Lpp <= 0.0_ReKi .or. p%Draft <= 0.0_ReKi ) call Fatal( 'Lpp and Draft must be positive.' )
      if ( p%UMin <= 0.0_ReKi ) call Fatal( 'UMin must be positive.' )
      if ( p%ResMod /= 1 .and. p%ResMod /= 2 ) call Fatal( 'ResMod must be 1 or 2.' )
      if ( WaveLdMod /= 0 ) call Fatal( 'WaveLdMod must be 0 (mean wave loads are not yet available).' )
      if ( p%SpdCtrlMod /= 0 .and. p%SpdCtrlMod /= 1 ) call Fatal( 'SpdCtrlMod must be 0 or 1.' )
      if ( p%AutoPilotMod < 0 .or. p%AutoPilotMod > 2 ) call Fatal( 'AutoPilotMod must be 0, 1, or 2.' )
      if ( p%CtrlWaveFc < 0.0_ReKi ) call Fatal( 'CtrlWaveFc must not be negative.' )
      if ( p%NProp > 0 ) then
         if ( abs( sum(p%ThrShr) - 1.0_ReKi ) > 1.0e-3_ReKi ) call Fatal( 'The thrust shares ThrShr must add up to 1.' )
         if ( p%TMax < p%TMin ) call Fatal( 'TMax must not be less than TMin.' )
         call CheckTau( p%TauT, 'TauT' )
      end if
      if ( p%NRud > 0 ) then
         call CheckTau( p%TauDelta, 'TauDelta' )
         if ( p%DeltaRateMax <= 0.0_ReKi .or. p%DeltaMax <= 0.0_ReKi ) call Fatal( 'DeltaRateMax and DeltaMax must be positive.' )
      end if
      if ( p%AutoPilotMod == 2 .and. p%PsiCmdMax <= 0.0_ReKi ) call Fatal( 'PsiCmdMax must be positive.' )
      if ( ErrStat >= AbortErrLev ) return

      do j = 1, p%NProp
         if ( p%PropD(j) <= 0.0_ReKi ) call Fatal( 'PropD must be positive (propeller '//trim(Num2LStr(j))//').' )
         if ( p%RudFlag(j) /= 0 .and. p%RudFlag(j) /= 1 ) call Fatal( 'RudFlag must be 0 or 1 (propeller '//trim(Num2LStr(j))//').' )
         if ( p%EtaR(j) <= 0.0_ReKi ) call Fatal( 'EtaR must be positive (propeller '//trim(Num2LStr(j))//').' )
         if ( .not. EqualRealNos( abs(p%RotDir(j)), 1.0_ReKi ) ) call Fatal( 'RotDir must be +1 or -1 (propeller '//trim(Num2LStr(j))//').' )
         if ( ErrStat >= AbortErrLev ) return

         select case ( p%PropOWMod(j) )
         case (0)
            p%HasKQ(j) = .false.
            p%J0(j)    = 0.0_ReKi
         case (1)
            k0 = p%KTPoly(1,j); k1 = p%KTPoly(2,j); k2 = p%KTPoly(3,j)
            Disc = k1**2 - 4.0_ReKi*k2*k0
            if ( k0 <= 0.0_ReKi .or. Disc < 0.0_ReKi .or. ( -k1 + sqrt(max(Disc,0.0_ReKi)) ) <= 0.0_ReKi ) then
               call Fatal( 'The quadratic K_T of propeller '//trim(Num2LStr(j))//' needs KT0 > 0 and a positive zero-thrust advance ratio.' )
               return
            end if
            p%J0(j) = 2.0_ReKi*k0/( -k1 + sqrt(Disc) )
            if ( k1 > 0.0_ReKi .or. k1 + 2.0_ReKi*k2*p%J0(j) > 0.0_ReKi ) &
               call Fatal( 'The quadratic K_T of propeller '//trim(Num2LStr(j))//' must decrease between J = 0 and the zero-thrust advance ratio.' )
            p%HasKQ(j) = any( abs( p%KQPoly(:,j) ) > 0.0_ReKi )
         case (2)
            p%J0(j)    = TableJ0( p%OWTbl(j) )
            p%HasKQ(j) = .true.
         case default
            call Fatal( 'PropOWMod must be 0, 1, or 2 (propeller '//trim(Num2LStr(j))//').' )
         end select

         if ( p%RudFlag(j) == 1 ) then
            if ( p%RudLiftMod(j) /= 1 .and. p%RudLiftMod(j) /= 2 ) call Fatal( 'RudLiftMod must be 1 or 2 (rudder '//trim(Num2LStr(j))//').' )
            if ( p%RudAR(j) <= 0.0_ReKi .or. p%RudHR(j) <= 0.0_ReKi ) call Fatal( 'RudAR and RudHR must be positive (rudder '//trim(Num2LStr(j))//').' )
         end if
         if ( ErrStat >= AbortErrLev ) return
      end do
      if ( p%NProp > 0 .and. p%EtaD <= 0.0_ReKi .and. .not. all(p%HasKQ) ) call Fatal( 'EtaD must be positive.' )
   end subroutine CheckParams

   subroutine CheckTau( Tau, Name )
      real(ReKi),   intent(in) :: Tau
      character(*), intent(in) :: Name
      if ( Tau <= real(p%DT,ReKi)/2.78_ReKi ) then
         call Fatal( Name//' must exceed DT/2.78 of HydroDyn for a stable explicit integration.' )
      else if ( Tau < 2.0_ReKi*real(p%DT,ReKi) ) then
         call SetErrStat( ErrID_Warn, 'Body '//trim(Num2LStr(p%BodyID))//': '//Name//' is less than twice the HydroDyn time step.', ErrStat, ErrMsg, RoutineName )
      end if
   end subroutine CheckTau

   real(ReKi) function TableJ0( Tbl )
      type(OWTableType), intent(in) :: Tbl
      integer(IntKi) :: j
      TableJ0 = 0.0_ReKi
      do j = 2, size(Tbl%J)
         if ( Tbl%KT(j) <= 0.0_ReKi ) then
            TableJ0 = Tbl%J(j-1) + ( Tbl%J(j) - Tbl%J(j-1) )*Tbl%KT(j-1)/( Tbl%KT(j-1) - Tbl%KT(j) )
            return
         end if
      end do
   end function TableJ0

   subroutine Fatal( Msg )
      character(*), intent(in) :: Msg
      call SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(p%BodyID))//': '//Msg, ErrStat, ErrMsg, RoutineName )
   end subroutine Fatal
end subroutine ReadVesselFile

!----------------------------------------------------------------------------------------------------------------------------------
!> Reads the resistance table: header, title, NRows, two column-header lines, rows of U (m/s), R (N).
subroutine ReadResTable( FileName, p, ErrStat, ErrMsg )
   character(*),               intent(in   ) :: FileName
   type(Vessel_ParameterType), intent(inout) :: p
   integer(IntKi),             intent(  out) :: ErrStat
   character(*),               intent(  out) :: ErrMsg

   type(FileInfoType)      :: TblInfo
   integer(IntKi)          :: nRows, iLine, j
   real(ReKi)              :: Vals(2)
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'ReadResTable'

   ErrStat = ErrID_None
   ErrMsg  = ''
   call ProcessComFile( FileName, TblInfo, ErrStat2, ErrMsg2 ); if (Failed()) return
   iLine = 3
   call ParseVar( TblInfo, iLine, 'NRows', nRows, ErrStat2, ErrMsg2 ); if (Failed()) return
   iLine = iLine + 2
   if ( nRows < 2 ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The resistance table "'//trim(FileName)//'" needs at least two rows.'
      if (Failed()) return
   end if
   call AllocAry( p%ResU, nRows, 'ResU', ErrStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry( p%ResR, nRows, 'ResR', ErrStat2, ErrMsg2 ); if (Failed()) return
   do j = 1, nRows
      call ParseAry( TblInfo, iLine, 'resistance table row '//trim(Num2LStr(j)), Vals, 2, ErrStat2, ErrMsg2 ); if (Failed()) return
      p%ResU(j) = Vals(1)
      p%ResR(j) = Vals(2)
   end do
   ErrStat2 = ErrID_None
   if ( .not. EqualRealNos( p%ResU(1), 0.0_ReKi ) .or. .not. EqualRealNos( p%ResR(1), 0.0_ReKi ) ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The first row of the resistance table "'//trim(FileName)//'" must be U = 0, R = 0.'
   else if ( any( p%ResU(2:nRows) <= p%ResU(1:nRows-1) ) ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The speeds of the resistance table "'//trim(FileName)//'" must increase strictly.'
   end if
   if (Failed()) return
   call NWTC_Library_DestroyFileInfoType( TblInfo, ErrStat2, ErrMsg2 )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call NWTC_Library_DestroyFileInfoType( TblInfo, ErrStat2, ErrMsg2 )
   end function Failed
end subroutine ReadResTable

!----------------------------------------------------------------------------------------------------------------------------------
!> Reads a propeller open-water table: header, title, NRows, two column-header lines, rows of J, K_T, K_Q.
subroutine ReadOWTable( FileName, Tbl, ErrStat, ErrMsg )
   character(*),      intent(in   ) :: FileName
   type(OWTableType), intent(inout) :: Tbl
   integer(IntKi),    intent(  out) :: ErrStat
   character(*),      intent(  out) :: ErrMsg

   type(FileInfoType)      :: TblInfo
   integer(IntKi)          :: nRows, iLine, j, jZero
   real(ReKi)              :: Vals(3)
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'ReadOWTable'

   ErrStat = ErrID_None
   ErrMsg  = ''
   call ProcessComFile( FileName, TblInfo, ErrStat2, ErrMsg2 ); if (Failed()) return
   iLine = 3
   call ParseVar( TblInfo, iLine, 'NRows', nRows, ErrStat2, ErrMsg2 ); if (Failed()) return
   iLine = iLine + 2
   if ( nRows < 2 ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The open-water table "'//trim(FileName)//'" needs at least two rows.'
      if (Failed()) return
   end if
   call AllocAry( Tbl%J,  nRows, 'J',  ErrStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry( Tbl%KT, nRows, 'KT', ErrStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry( Tbl%KQ, nRows, 'KQ', ErrStat2, ErrMsg2 ); if (Failed()) return
   do j = 1, nRows
      call ParseAry( TblInfo, iLine, 'open-water table row '//trim(Num2LStr(j)), Vals, 3, ErrStat2, ErrMsg2 ); if (Failed()) return
      Tbl%J(j)  = Vals(1)
      Tbl%KT(j) = Vals(2)
      Tbl%KQ(j) = Vals(3)
   end do

      ! J must start at 0 and increase; K_T must be positive at J = 0, reach zero, and decrease until it does
   jZero = 0
   do j = 2, nRows
      if ( Tbl%KT(j) <= 0.0_ReKi ) then
         jZero = j
         exit
      end if
   end do
   ErrStat2 = ErrID_None
   if ( .not. EqualRealNos( Tbl%J(1), 0.0_ReKi ) ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The open-water table "'//trim(FileName)//'" must start at J = 0.'
   else if ( any( Tbl%J(2:nRows) <= Tbl%J(1:nRows-1) ) ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The advance ratios of the open-water table "'//trim(FileName)//'" must increase strictly.'
   else if ( Tbl%KT(1) <= 0.0_ReKi .or. jZero == 0 ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The open-water table "'//trim(FileName)//'" must start with K_T > 0 and reach K_T <= 0.'
   else if ( any( Tbl%KT(2:jZero) > Tbl%KT(1:jZero-1) ) ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'K_T of the open-water table "'//trim(FileName)//'" must decrease up to the zero-thrust advance ratio.'
   end if
   if (Failed()) return
   call NWTC_Library_DestroyFileInfoType( TblInfo, ErrStat2, ErrMsg2 )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call NWTC_Library_DestroyFileInfoType( TblInfo, ErrStat2, ErrMsg2 )
   end function Failed
end subroutine ReadOWTable

!----------------------------------------------------------------------------------------------------------------------------------
!> Reads a rudder lift/drag table: header, title, NRows, two column-header lines, rows of alpha (deg), C_L, C_D.
subroutine ReadRudTable( FileName, Tbl, ErrStat, ErrMsg )
   character(*),       intent(in   ) :: FileName
   type(RudTableType), intent(inout) :: Tbl
   integer(IntKi),     intent(  out) :: ErrStat
   character(*),       intent(  out) :: ErrMsg

   type(FileInfoType)      :: TblInfo
   integer(IntKi)          :: nRows, iLine, j
   real(ReKi)              :: Vals(3)
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'ReadRudTable'

   ErrStat = ErrID_None
   ErrMsg  = ''
   call ProcessComFile( FileName, TblInfo, ErrStat2, ErrMsg2 ); if (Failed()) return
   iLine = 3
   call ParseVar( TblInfo, iLine, 'NRows', nRows, ErrStat2, ErrMsg2 ); if (Failed()) return
   iLine = iLine + 2
   if ( nRows < 2 ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The rudder table "'//trim(FileName)//'" needs at least two rows.'
      if (Failed()) return
   end if
   call AllocAry( Tbl%Alpha, nRows, 'Alpha', ErrStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry( Tbl%Cl,    nRows, 'Cl',    ErrStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry( Tbl%Cd,    nRows, 'Cd',    ErrStat2, ErrMsg2 ); if (Failed()) return
   do j = 1, nRows
      call ParseAry( TblInfo, iLine, 'rudder table row '//trim(Num2LStr(j)), Vals, 3, ErrStat2, ErrMsg2 ); if (Failed()) return
      Tbl%Alpha(j) = Vals(1)*D2R
      Tbl%Cl(j)    = Vals(2)
      Tbl%Cd(j)    = Vals(3)
   end do
   if ( any( Tbl%Alpha(2:nRows) <= Tbl%Alpha(1:nRows-1) ) ) then
      ErrStat2 = ErrID_Fatal; ErrMsg2 = 'The angles of attack of the rudder table "'//trim(FileName)//'" must increase strictly.'
      if (Failed()) return
   end if
   call NWTC_Library_DestroyFileInfoType( TblInfo, ErrStat2, ErrMsg2 )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call NWTC_Library_DestroyFileInfoType( TblInfo, ErrStat2, ErrMsg2 )
   end function Failed
end subroutine ReadRudTable

!----------------------------------------------------------------------------------------------------------------------------------
!> Creates the one-node input mesh at the undisplaced body reference point.
subroutine InitMesh( InitInp, u, ErrStat, ErrMsg )
   type(Vessel_InitInputType), intent(in   ) :: InitInp
   type(Vessel_InputType),     intent(inout) :: u
   integer(IntKi),             intent(  out) :: ErrStat
   character(*),               intent(  out) :: ErrMsg

   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'InitMesh'

   ErrStat = ErrID_None
   ErrMsg  = ''
   call MeshCreate( BlankMesh = u%Mesh, IOS = COMPONENT_INPUT, Nnodes = 1, ErrStat = ErrStat2, ErrMess = ErrMsg2, &
                    TranslationDisp = .true., Orientation = .true., TranslationVel = .true., RotationVel = .true. )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MeshPositionNode( u%Mesh, 1, InitInp%BodyRefPos, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MeshConstructElement( u%Mesh, ELEMENT_POINT, ErrStat2, ErrMsg2, 1 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call MeshCommit( u%Mesh, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
end subroutine InitMesh

!----------------------------------------------------------------------------------------------------------------------------------
subroutine InitMisc( p, m, ErrStat, ErrMsg )
   type(Vessel_ParameterType), intent(in   ) :: p
   type(Vessel_MiscVarType),   intent(inout) :: m
   integer(IntKi),             intent(  out) :: ErrStat
   character(*),               intent(  out) :: ErrMsg

   integer(IntKi)          :: n
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'InitMisc'

   ErrStat = ErrID_None
   ErrMsg  = ''
   n = p%NProp
   call AllocAry( m%PrpThr,    n, 'PrpThr',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpTCmd,   n, 'PrpTCmd',   ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpUA,     n, 'PrpUA',     ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpWP,     n, 'PrpWP',     ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpPwrT,   n, 'PrpPwrT',   ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpPwrIdl, n, 'PrpPwrIdl', ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpImm,    n, 'PrpImm',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpJ,      n, 'PrpJ',      ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpRPM,    n, 'PrpRPM',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpTrq,    n, 'PrpTrq',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%PrpPwrD,   n, 'PrpPwrD',   ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudAng,    n, 'RudAng',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudCmd,    n, 'RudCmd',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudUR,     n, 'RudUR',     ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudVR,     n, 'RudVR',     ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudAoA,    n, 'RudAoA',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudFN,     n, 'RudFN',     ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudFT,     n, 'RudFT',     ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call AllocAry( m%RudImm,    n, 'RudImm',    ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if ( ErrStat >= AbortErrLev ) return
   m%PrpThr = 0.0_ReKi;  m%PrpTCmd = 0.0_ReKi;  m%PrpUA  = 0.0_ReKi;  m%PrpWP   = 0.0_ReKi
   m%PrpPwrT = 0.0_ReKi; m%PrpPwrIdl = 0.0_ReKi; m%PrpImm = 0.0_ReKi;  m%PrpJ    = 0.0_ReKi
   m%PrpRPM = 0.0_ReKi;  m%PrpTrq  = 0.0_ReKi;  m%PrpPwrD = 0.0_ReKi
   m%RudAng = 0.0_ReKi;  m%RudCmd  = 0.0_ReKi;  m%RudUR  = 0.0_ReKi;  m%RudVR   = 0.0_ReKi
   m%RudAoA = 0.0_ReKi;  m%RudFN   = 0.0_ReKi;  m%RudFT  = 0.0_ReKi;  m%RudImm  = 0.0_ReKi
   m%FHull  = 0.0_ReKi;  m%FCor    = 0.0_ReKi;  m%FTot   = 0.0_ReKi
end subroutine InitMisc

!----------------------------------------------------------------------------------------------------------------------------------
subroutine WriteSummary( Un, p, mxA0, myA0 )
   integer(IntKi),             intent(in) :: Un
   type(Vessel_ParameterType), intent(in) :: p
   real(ReKi),                 intent(in) :: mxA0, myA0

   real(ReKi)     :: qNd, J, KT, KQ, Eta0, Eta0Max, JMax
   integer(IntKi) :: i, k

   qNd = 0.5_ReKi*p%RhoW*p%Lpp**2*p%Draft
   write(Un,'(//A,I5)') 'Vessel (MMG) Model - Body', p%BodyID
   write(Un,'(A)')      '------------------------------'
   write(Un,'(A42,2X,ES15.6)') "MMG m'_x used                          :", p%mx/qNd
   write(Un,'(A42,2X,ES15.6)') "MMG m'_y used                          :", p%my/qNd
   write(Un,'(A42,2X,ES15.6)') "A0_11/(0.5 rho Lpp^2 d)                :", mxA0
   write(Un,'(A42,2X,ES15.6)') "A0_22/(0.5 rho Lpp^2 d)                :", myA0
   write(Un,'(A42,2X,ES15.6)') 'Resistance at FrameSpd (N)             :', Resistance( p, abs(p%FrameSpd) )
   write(Un,'(A42,2X,ES15.6)') 'Feed-forward total thrust T0 (N)       :', p%T0
   write(Un,'(A42,2X,I15)')    'Number of propellers                   :', p%NProp
   write(Un,'(A42,2X,I15)')    'Number of rudders                      :', p%NRud
   write(Un,'(A42,2X,I15)')    'Number of continuous states            :', p%NStates
   do i = 1, p%NProp
      if ( p%PropOWMod(i) == 0 ) cycle
      Eta0Max = 0.0_ReKi
      JMax    = 0.0_ReKi
      do k = 1, 100
         J = p%J0(i)*real(k,ReKi)/101.0_ReKi
         call OpenWater( p, i, J, KT, KQ )
         if ( KQ > 0.0_ReKi ) then
            Eta0 = J*KT/( TwoPi*KQ )
            if ( Eta0 > Eta0Max ) then
               Eta0Max = Eta0
               JMax    = J
            end if
         end if
      end do
      write(Un,'(A,I2,A,ES12.4,A,ES12.4,A,ES12.4)') 'Propeller ', i, ': zero-thrust J0 = ', p%J0(i), &
         ', max. open-water efficiency eta0 = ', Eta0Max, ' at J = ', JMax
   end do
end subroutine WriteSummary

!----------------------------------------------------------------------------------------------------------------------------------
!> Calm-water resistance R(U).
real(ReKi) function Resistance( p, U )
   type(Vessel_ParameterType), intent(in) :: p
   real(ReKi),                 intent(in) :: U

   integer(IntKi) :: j, n

   if ( p%ResMod == 1 ) then
      Resistance = 0.5_ReKi*p%RhoW*p%Lpp*p%Draft*U**2*p%R0
   else
      n = size(p%ResU)
      if ( U >= p%ResU(n) ) then
         Resistance = p%ResR(n)*( U/p%ResU(n) )**2   ! quadratic beyond the table
      else
         do j = 2, n
            if ( U <= p%ResU(j) ) exit
         end do
         Resistance = p%ResR(j-1) + ( p%ResR(j) - p%ResR(j-1) )*( U - p%ResU(j-1) )/( p%ResU(j) - p%ResU(j-1) )
      end if
   end if
end function Resistance

!----------------------------------------------------------------------------------------------------------------------------------
!> Open-water K_T and K_Q of propeller i at advance ratio J.
subroutine OpenWater( p, i, J, KT, KQ )
   type(Vessel_ParameterType), intent(in   ) :: p
   integer(IntKi),             intent(in   ) :: i
   real(ReKi),                 intent(in   ) :: J
   real(ReKi),                 intent(  out) :: KT, KQ

   integer(IntKi) :: k, n
   real(ReKi)     :: f

   if ( p%PropOWMod(i) == 1 ) then
      KT = p%KTPoly(1,i) + p%KTPoly(2,i)*J + p%KTPoly(3,i)*J**2
      KQ = p%KQPoly(1,i) + p%KQPoly(2,i)*J + p%KQPoly(3,i)*J**2
   else
      n = size(p%OWTbl(i)%J)
      do k = 2, n-1
         if ( J <= p%OWTbl(i)%J(k) ) exit
      end do
      k = min( k, n )
      f  = ( J - p%OWTbl(i)%J(k-1) )/( p%OWTbl(i)%J(k) - p%OWTbl(i)%J(k-1) )
      KT = p%OWTbl(i)%KT(k-1) + f*( p%OWTbl(i)%KT(k) - p%OWTbl(i)%KT(k-1) )
      KQ = p%OWTbl(i)%KQ(k-1) + f*( p%OWTbl(i)%KQ(k) - p%OWTbl(i)%KQ(k-1) )
   end if
end subroutine OpenWater

!----------------------------------------------------------------------------------------------------------------------------------
!> Propeller operating point by thrust identity (Theory 7.4); T and uA must be non-negative.
subroutine PropOperatingPoint( p, i, T, uA, J, n, Q )
   type(Vessel_ParameterType), intent(in   ) :: p
   integer(IntKi),             intent(in   ) :: i
   real(ReKi),                 intent(in   ) :: T, uA
   real(ReKi),                 intent(  out) :: J, n, Q

   real(ReKi)     :: D, c, k0, k1, k2, Disc, KT, KQ, fa, fb, b, B2, Den
   integer(IntKi) :: s, nRows

   J = 0.0_ReKi
   n = 0.0_ReKi
   Q = 0.0_ReKi
   if ( p%PropOWMod(i) == 0 ) return
   if ( T <= 0.0_ReKi .and. uA <= 0.0_ReKi ) return
   D = p%PropD(i)

   if ( p%PropOWMod(i) == 1 ) then
      k0 = p%KTPoly(1,i); k1 = p%KTPoly(2,i); k2 = p%KTPoly(3,i)
      Disc = k1**2*uA**2 - 4.0_ReKi*k0*( k2*uA**2 - T/( p%RhoW*D**2 ) )
      J = 2.0_ReKi*k0*uA/( -k1*uA + sqrt( max( Disc, 0.0_ReKi ) ) )
   else
         ! f(J) = c K_T(J) - T J^2 decreases; find the segment containing its root and solve the segment quadratic
      c     = p%RhoW*D**2*uA**2
      nRows = size(p%OWTbl(i)%J)
      do s = 1, nRows-1
         fb = c*p%OWTbl(i)%KT(s+1) - T*p%OWTbl(i)%J(s+1)**2
         if ( fb <= 0.0_ReKi ) exit
      end do
      s  = min( s, nRows-1 )
      fa = c*p%OWTbl(i)%KT(s) - T*p%OWTbl(i)%J(s)**2
      if ( fa <= 0.0_ReKi ) then
         J = p%OWTbl(i)%J(s)
      else
         b   = ( p%OWTbl(i)%KT(s+1) - p%OWTbl(i)%KT(s) )/( p%OWTbl(i)%J(s+1) - p%OWTbl(i)%J(s) )
         B2  = 2.0_ReKi*T*p%OWTbl(i)%J(s) - c*b
         Den = B2 + sqrt( B2**2 + 4.0_ReKi*T*fa )
         if ( Den > 0.0_ReKi ) then
            J = p%OWTbl(i)%J(s) + 2.0_ReKi*fa/Den
         else
            J = p%OWTbl(i)%J(s+1)
         end if
      end if
   end if
   J = min( max( J, 0.0_ReKi ), p%J0(i) )

   call OpenWater( p, i, J, KT, KQ )
   if ( J <= 0.5_ReKi*p%J0(i) ) then
      n = sqrt( T/( p%RhoW*D**4*KT ) )
   else
      n = uA/( J*D )
   end if
   if ( p%HasKQ(i) ) Q = KQ/p%EtaR(i)*p%RhoW*n**2*D**5
end subroutine PropOperatingPoint

!----------------------------------------------------------------------------------------------------------------------------------
!> Rudder normal and tangential forces per unit dynamic pressure area (C_N, C_T) at angle of attack Alpha.
subroutine RudderCoefs( p, i, Alpha, CN, CT )
   type(Vessel_ParameterType), intent(in   ) :: p
   integer(IntKi),             intent(in   ) :: i
   real(ReKi),                 intent(in   ) :: Alpha
   real(ReKi),                 intent(  out) :: CN, CT

   real(ReKi)     :: a, Cl, Cd, f
   integer(IntKi) :: k, n

   if ( p%RudLiftMod(i) == 1 ) then
      a = Alpha
      if ( p%StallAng(i) > 0.0_ReKi ) a = sign( min( abs(Alpha), p%StallAng(i) ), Alpha )
      CN = p%FAlpha(i)*sin(a)
      CT = 0.0_ReKi
   else
      n = size(p%RudTbl(i)%Alpha)
      if ( Alpha <= p%RudTbl(i)%Alpha(1) ) then
         Cl = p%RudTbl(i)%Cl(1); Cd = p%RudTbl(i)%Cd(1)
      else if ( Alpha >= p%RudTbl(i)%Alpha(n) ) then
         Cl = p%RudTbl(i)%Cl(n); Cd = p%RudTbl(i)%Cd(n)
      else
         do k = 2, n
            if ( Alpha <= p%RudTbl(i)%Alpha(k) ) exit
         end do
         f  = ( Alpha - p%RudTbl(i)%Alpha(k-1) )/( p%RudTbl(i)%Alpha(k) - p%RudTbl(i)%Alpha(k-1) )
         Cl = p%RudTbl(i)%Cl(k-1) + f*( p%RudTbl(i)%Cl(k) - p%RudTbl(i)%Cl(k-1) )
         Cd = p%RudTbl(i)%Cd(k-1) + f*( p%RudTbl(i)%Cd(k) - p%RudTbl(i)%Cd(k-1) )
      end if
      CN = Cl*cos(Alpha) + Cd*sin(Alpha)
      CT = Cd*cos(Alpha) - Cl*sin(Alpha)
   end if
end subroutine RudderCoefs

!----------------------------------------------------------------------------------------------------------------------------------
!> Body-axis (OpenFAST) velocity of body point Rho relative to the mean current.
function RelVelBody( p, k, Rho ) result( V )
   type(Vessel_ParameterType), intent(in) :: p
   type(VslKinType),           intent(in) :: k
   real(ReKi),                 intent(in) :: Rho(3)
   real(R8Ki)                             :: V(3)

   V = k%VelO + cross_product( k%Omega, matmul( transpose(k%Rg2b), real(Rho,R8Ki) ) )
   V(1) = V(1) + real(p%FrameSpd,R8Ki)
   V = matmul( k%Rg2b, V )
end function RelVelBody

!----------------------------------------------------------------------------------------------------------------------------------
!> Kinematics of the instance from the input mesh (Theory 3.2, 3.3, 10).
subroutine CalcKinematics( p, u, k )
   type(Vessel_ParameterType), intent(in   ) :: p
   type(Vessel_InputType),     intent(in   ) :: u
   type(VslKinType),           intent(  out) :: k

   real(R8Ki) :: V(3), wb(3), RhoG(3)
   real(ReKi) :: Ureg

   k%Rg2b  = u%Mesh%Orientation(:,:,1)
   k%PosO  = real( u%Mesh%Position(:,1), R8Ki ) + u%Mesh%TranslationDisp(:,1)
   k%VelO  = real( u%Mesh%TranslationVel(:,1), R8Ki )
   k%Omega = real( u%Mesh%RotationVel(:,1), R8Ki )
   wb      = matmul( k%Rg2b, k%Omega )

      ! Hull force center: hull model, drift angle and nondimensional variables
   V      = RelVelBody( p, k, p%HullRef + (/ 0.0_ReKi, 0.0_ReKi, p%zH /) )
   k%u    = real(  V(1), ReKi )
   k%vm   = real( -V(2), ReKi )
   k%r    = real( -wb(3), ReKi )
   k%UTot = sqrt( k%u**2 + k%vm**2 )
   if ( k%UTot > 0.0_ReKi ) then
      k%Beta = atan2( -k%vm, k%u )
   else
      k%Beta = 0.0_ReKi
   end if
   Ureg = max( k%UTot, p%UMin )
   k%vp = k%vm/Ureg
   k%rp = k%r*p%Lpp/Ureg

      ! Hull reference point: MMG Coriolis term
   V     = RelVelBody( p, k, p%HullRef )
   k%uH  = real(  V(1), ReKi )
   k%vmH = real( -V(2), ReKi )

      ! Controller inputs: position and velocity of the hull reference point in the frame, heading, body yaw rate
   RhoG = matmul( transpose(k%Rg2b), real(p%HullRef,R8Ki) )
   V    = k%VelO + cross_product( k%Omega, RhoG )
   k%Sig(1) = real( k%PosO(1) + RhoG(1), ReKi ) - p%XYH0(1)
   k%Sig(2) = real( V(1), ReKi )
   k%Sig(3) = real( k%PosO(2) + RhoG(2), ReKi ) - p%XYH0(2)
   k%Sig(4) = real( V(2), ReKi )
   k%Sig(5) = real( atan2( k%Rg2b(1,2), k%Rg2b(1,1) ), ReKi )
   k%Sig(6) = real( wb(3), ReKi )
end subroutine CalcKinematics

!----------------------------------------------------------------------------------------------------------------------------------
!> Speed and steering controllers (Theory 10): saturated commands and integrator derivatives with anti-windup.
subroutine Controller( p, k, x, TCmd, DeltaCmd, PsiCmd, dXi )
   type(Vessel_ParameterType), intent(in   ) :: p
   type(VslKinType),           intent(in   ) :: k
   real(R8Ki),                 intent(in   ) :: x(:)
   real(ReKi),                 intent(  out) :: TCmd(:)
   real(ReKi),                 intent(  out) :: DeltaCmd, PsiCmd
   real(ReKi),                 intent(  out) :: dXi(3)

   real(ReKi)     :: s(NFltSig), Xi(3), Tc, Tu, dc, PsiU, ePsi, DeltaU
   integer(IntKi) :: i

   if ( p%CtrlWaveFc > 0.0_ReKi ) then
      s = real( x(p%iFlt+1:p%iFlt+NFltSig), ReKi )
   else
      s = k%Sig
   end if
   Xi  = real( x(p%iXi+1:p%iXi+3), ReKi )
   dXi = 0.0_ReKi

      ! Speed
   if ( p%SpdCtrlMod == 1 ) then
      Tc     = p%T0 - p%KPx*s(1) - p%KDx*s(2) - p%KIx*Xi(1)
      dXi(1) = s(1)
   else
      Tc     = p%T0
   end if
   do i = 1, p%NProp
      Tu      = p%ThrShr(i)*Tc
      TCmd(i) = min( max( Tu, p%TMin ), p%TMax )
      if ( p%SpdCtrlMod == 1 ) then
         dc = -p%KIx*p%ThrShr(i)*s(1)
         if ( ( Tu > p%TMax .and. dc > 0.0_ReKi ) .or. ( Tu < p%TMin .and. dc < 0.0_ReKi ) ) dXi(1) = 0.0_ReKi
      end if
   end do

      ! Steering
   select case ( p%AutoPilotMod )
   case (0)
      PsiCmd   = 0.0_ReKi
      DeltaCmd = min( max( p%Delta0, -p%DeltaMax ), p%DeltaMax )
   case default
      if ( p%AutoPilotMod == 2 ) then
         PsiU   = -( p%KPy*s(3) + p%KIy*Xi(2) + p%KDy*s(4) )
         PsiCmd = min( max( PsiU, -p%PsiCmdMax ), p%PsiCmdMax )
         dXi(2) = s(3)
         dc     = -p%KIy*s(3)
         if ( ( PsiU > p%PsiCmdMax .and. dc > 0.0_ReKi ) .or. ( PsiU < -p%PsiCmdMax .and. dc < 0.0_ReKi ) ) dXi(2) = 0.0_ReKi
      else
         PsiCmd = p%PsiRef
      end if
      ePsi     = WrapAngle( s(5) - PsiCmd )
      DeltaU   = p%KPpsi*ePsi + p%KIpsi*Xi(3) + p%KDpsi*s(6)
      DeltaCmd = min( max( DeltaU, -p%DeltaMax ), p%DeltaMax )
      dXi(3)   = ePsi
      dc       = p%KIpsi*ePsi
      if ( ( DeltaU > p%DeltaMax .and. dc > 0.0_ReKi ) .or. ( DeltaU < -p%DeltaMax .and. dc < 0.0_ReKi ) ) dXi(3) = 0.0_ReKi
   end select
end subroutine Controller

!----------------------------------------------------------------------------------------------------------------------------------
real(ReKi) function WrapAngle( a )
   real(ReKi), intent(in) :: a
   WrapAngle = atan2( sin(a), cos(a) )
end function WrapAngle

!----------------------------------------------------------------------------------------------------------------------------------
!> State derivatives from the kinematics (Theory 7.2, 8.5, 10, 15.4).
subroutine StateDeriv( p, k, x, dxdt )
   type(Vessel_ParameterType), intent(in   ) :: p
   type(VslKinType),           intent(in   ) :: k
   real(R8Ki),                 intent(in   ) :: x(:)
   real(R8Ki),                 intent(  out) :: dxdt(:)

   real(ReKi)     :: TCmd(p%NProp), DeltaCmd, PsiCmd, dXi(3), Rate
   integer(IntKi) :: i, j

   call Controller( p, k, x, TCmd, DeltaCmd, PsiCmd, dXi )
   dxdt = 0.0_R8Ki
   do i = 1, p%NProp
      dxdt(i) = ( real(TCmd(i),R8Ki) - x(i) )/real(p%TauT,R8Ki)
      j = p%iRudState(i)
      if ( j > 0 ) then
         Rate    = ( DeltaCmd - real(x(j),ReKi) )/p%TauDelta
         dxdt(j) = real( min( max( Rate, -p%DeltaRateMax ), p%DeltaRateMax ), R8Ki )
      end if
   end do
   dxdt(p%iXi+1:p%iXi+3) = real( dXi, R8Ki )
   if ( p%CtrlWaveFc > 0.0_ReKi ) then
      do i = 1, NFltSig
         if ( i == 5 ) then
            dxdt(p%iFlt+i) = real( p%CtrlWaveFc*WrapAngle( k%Sig(i) - real(x(p%iFlt+i),ReKi) ), R8Ki )
         else
            dxdt(p%iFlt+i) = real( p%CtrlWaveFc*( k%Sig(i) - real(x(p%iFlt+i),ReKi) ), R8Ki )
         end if
      end do
   end if
end subroutine StateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
!> Continuous-state derivatives (used for linearization).
subroutine Vessel_CalcContStateDeriv( u, p, x, dxdt, ErrStat, ErrMsg )
   type(Vessel_InputType),           intent(in   ) :: u
   type(Vessel_ParameterType),       intent(in   ) :: p
   type(Vessel_ContinuousStateType), intent(in   ) :: x
   type(Vessel_ContinuousStateType), intent(inout) :: dxdt
   integer(IntKi),                   intent(  out) :: ErrStat
   character(*),                     intent(  out) :: ErrMsg

   type(VslKinType) :: k

   ErrStat = ErrID_None
   ErrMsg  = ''
   call CalcKinematics( p, u, k )
   call StateDeriv( p, k, x%x, dxdt%x )
end subroutine Vessel_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
!> Advances the states from t to t + DT with RK4, using inputs interpolated from Inputs at InputTimes. Re-entrant (Theory 11.1).
subroutine Vessel_UpdateStates( t, Inputs, InputTimes, p, x, ErrStat, ErrMsg )
   real(DbKi),                       intent(in   ) :: t
   type(Vessel_InputType),           intent(inout) :: Inputs(:)
   real(DbKi),                       intent(in   ) :: InputTimes(:)
   type(Vessel_ParameterType),       intent(in   ) :: p
   type(Vessel_ContinuousStateType), intent(inout) :: x
   integer(IntKi),                   intent(  out) :: ErrStat
   character(*),                     intent(  out) :: ErrMsg

   type(Vessel_InputType)  :: uInterp
   type(VslKinType)        :: k
   real(R8Ki)              :: x0(p%NStates), k1(p%NStates), k2(p%NStates), k3(p%NStates), k4(p%NStates), dt
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'Vessel_UpdateStates'

   ErrStat = ErrID_None
   ErrMsg  = ''
   dt = real( p%DT, R8Ki )
   x0 = x%x

   call Vessel_CopyInput( Inputs(1), uInterp, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); if (Failed()) return

   call Vessel_Input_ExtrapInterp( Inputs, InputTimes, uInterp, t, ErrStat2, ErrMsg2 ); if (Failed()) return
   call CalcKinematics( p, uInterp, k )
   call StateDeriv( p, k, x0, k1 )

   call Vessel_Input_ExtrapInterp( Inputs, InputTimes, uInterp, t + 0.5_DbKi*p%DT, ErrStat2, ErrMsg2 ); if (Failed()) return
   call CalcKinematics( p, uInterp, k )
   call StateDeriv( p, k, x0 + 0.5_R8Ki*dt*k1, k2 )
   call StateDeriv( p, k, x0 + 0.5_R8Ki*dt*k2, k3 )

   call Vessel_Input_ExtrapInterp( Inputs, InputTimes, uInterp, t + p%DT, ErrStat2, ErrMsg2 ); if (Failed()) return
   call CalcKinematics( p, uInterp, k )
   call StateDeriv( p, k, x0 + dt*k3, k4 )

   x%x = x0 + dt/6.0_R8Ki*( k1 + 2.0_R8Ki*k2 + 2.0_R8Ki*k3 + k4 )

   call Vessel_DestroyInput( uInterp, ErrStat2, ErrMsg2 )

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Vessel_DestroyInput( uInterp, ErrStat2, ErrMsg2 )
   end function Failed
end subroutine Vessel_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Vessel loads about the body reference point in global axes (Theory 5-9), and the quantities kept for output channels.
subroutine Vessel_CalcOutput( Time, u, p, x, m, Force, Moment, ErrStat, ErrMsg )
   real(DbKi),                       intent(in   ) :: Time
   type(Vessel_InputType),           intent(in   ) :: u
   type(Vessel_ParameterType),       intent(in   ) :: p
   type(Vessel_ContinuousStateType), intent(in   ) :: x
   type(Vessel_MiscVarType),         intent(inout) :: m
   real(ReKi),                       intent(  out) :: Force(3)
   real(ReKi),                       intent(  out) :: Moment(3)
   integer(IntKi),                   intent(  out) :: ErrStat
   character(*),                     intent(  out) :: ErrMsg

   type(VslKinType)        :: k
   real(ReKi)              :: TCmd(p%NProp), DeltaCmd, PsiCmd, dXi(3)
   real(ReKi)              :: FB(3), MB(3), FHB(3), MHB(3), FCB(3), MCB(3)
   real(ReKi)              :: q2, XH, YH, NH, vp, rp
   real(ReKi)              :: o(3), Elev, uPc, xPnd, BetaP, C2, OneMinusW, uA, uAp, T, Tp, J, n, Q, PwrE
   real(ReKi)              :: Delta, BetaR, Gam, Eta, Wash, uR, vR, Alpha, UR2, CN, CT, FN, FT, cd, sd
   real(R8Ki)              :: V(3)
   integer(IntKi)          :: i
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'Vessel_CalcOutput'

   ErrStat = ErrID_None
   ErrMsg  = ''

   call CalcKinematics( p, u, k )
   call Controller( p, k, x%x, TCmd, DeltaCmd, PsiCmd, dXi )

   FB = 0.0_ReKi; MB = 0.0_ReKi

      ! Hull (Theory 5): X_H at the hull reference point, Y_H at the hull force center, N_H pure (MMG -> body: y, z, N flip)
   vp = k%vp
   rp = k%rp
   q2 = 0.5_ReKi*p%RhoW*p%Lpp*p%Draft*k%UTot**2
   XH = -Resistance( p, k%UTot ) + q2*( p%XDer(1)*vp**2 + p%XDer(2)*vp*rp + p%XDer(3)*rp**2 + p%XDer(4)*vp**4 )
   YH = q2*( p%YDer(1)*vp + p%YDer(2)*rp + p%YDer(3)*vp**3 + p%YDer(4)*vp**2*rp + p%YDer(5)*vp*rp**2 + p%YDer(6)*rp**3 )
   NH = q2*p%Lpp*( p%NDer(1)*vp + p%NDer(2)*rp + p%NDer(3)*vp**3 + p%NDer(4)*vp**2*rp + p%NDer(5)*vp*rp**2 + p%NDer(6)*rp**3 )
   FHB = 0.0_ReKi; MHB = 0.0_ReKi
   call AddPointLoad( p%HullRef, (/ XH, 0.0_ReKi, 0.0_ReKi /), FHB, MHB )
   call AddPointLoad( p%HullRef + (/ 0.0_ReKi, 0.0_ReKi, p%zH /), (/ 0.0_ReKi, -YH, 0.0_ReKi /), FHB, MHB )
   MHB(3) = MHB(3) - NH

      ! MMG added-mass Coriolis at the hull reference point (Theory 6.2): (m_y v r, -m_x u r) in MMG axes
   FCB = 0.0_ReKi; MCB = 0.0_ReKi
   call AddPointLoad( p%HullRef, (/ p%my*k%vmH*k%r, p%mx*k%uH*k%r, 0.0_ReKi /), FCB, MCB )

   FB = FHB + FCB
   MB = MHB + MCB

      ! Propellers (Theory 7)
   m%PwrE = 0.0_ReKi
   m%PwrD = 0.0_ReKi
   do i = 1, p%NProp
      call Orbital( p%PropPos(:,i), o, Elev ); if (ErrStat >= AbortErrLev) return
      V         = RelVelBody( p, k, p%PropPos(:,i) )
      uPc       = real( V(1), ReKi )
      xPnd      = ( p%PropPos(1,i) - p%HullRef(1) )/p%Lpp
      BetaP     = k%Beta - xPnd*rp
      if ( BetaP > 0.0_ReKi ) then
         C2 = p%C2Pos(i)
      else
         C2 = p%C2Neg(i)
      end if
      OneMinusW = ( 1.0_ReKi - p%wP0(i) )*( 1.0_ReKi + ( 1.0_ReKi - exp( -p%C1(i)*abs(BetaP) ) )*( C2 - 1.0_ReKi ) )
      uA        = OneMinusW*uPc - o(1)
      T         = real( x%x(i), ReKi )
      uAp       = max( uA, 0.0_ReKi )
      Tp        = max( T, 0.0_ReKi )
      if ( uA < 0.0_ReKi .and. m%FirstWarn_Astern ) then
         call SetErrStat( ErrID_Warn, 'Body '//trim(Num2LStr(p%BodyID))//': negative propeller advance speed; the propeller and rudder '// &
                          'models are clamped to u_A = 0. Further warnings are suppressed.', ErrStat, ErrMsg, RoutineName )
         m%FirstWarn_Astern = .false.
      end if

      call AddPointLoad( p%PropPos(:,i), (/ ( 1.0_ReKi - p%tP(i) )*T, 0.0_ReKi, 0.0_ReKi /), FB, MB )

      call PropOperatingPoint( p, i, Tp, uAp, J, n, Q )
      MB(1) = MB(1) - p%RotDir(i)*Q

      PwrE = k%uH*( 1.0_ReKi - p%tP(i) )*T
      m%PrpThr(i)    = T
      m%PrpTCmd(i)   = TCmd(i)
      m%PrpUA(i)     = uA
      m%PrpWP(i)     = 1.0_ReKi - OneMinusW
      m%PrpPwrT(i)   = T*uA
      m%PrpPwrIdl(i) = 0.5_ReKi*Tp*( uAp + sqrt( uAp**2 + 8.0_ReKi*Tp/( p%RhoW*Pi*p%PropD(i)**2 ) ) )
      m%PrpImm(i)    = Elev - ( PointZ( p%PropPos(:,i) ) + 0.5_ReKi*p%PropD(i) )
      m%PrpJ(i)      = J
      m%PrpRPM(i)    = 60.0_ReKi*n
      m%PrpTrq(i)    = Q
      if ( p%HasKQ(i) ) then
         m%PrpPwrD(i) = TwoPi*n*Q
      else
         m%PrpPwrD(i) = PwrE/p%EtaD
      end if
      m%PwrE = m%PwrE + PwrE
      m%PwrD = m%PwrD + m%PrpPwrD(i)
      call CheckEmerged( m%PrpImm(i), 'propeller', i )

         ! Rudder behind this propeller (Theory 8)
      if ( p%iRudState(i) == 0 ) cycle
      Delta = real( x%x(p%iRudState(i)), ReKi )
      call Orbital( p%RudPos(:,i), o, Elev ); if (ErrStat >= AbortErrLev) return
      BetaR = k%Beta - p%lRnd(i)*rp
      if ( BetaR > 0.0_ReKi ) then
         Gam = p%GamRPos(i)
      else
         Gam = p%GamRNeg(i)
      end if
      Eta   = p%PropD(i)/p%RudHR(i)
      Wash  = uAp + p%Kappa(i)*( sqrt( uAp**2 + 8.0_ReKi*Tp/( p%RhoW*Pi*p%PropD(i)**2 ) ) - uAp )
      uR    = p%Eps(i)*sqrt( Eta*Wash**2 + ( 1.0_ReKi - Eta )*uAp**2 )
      vR    = Gam*k%UTot*BetaR + o(2)
      Alpha = Delta - atan2( vR, uR )
      UR2   = uR**2 + vR**2
      call RudderCoefs( p, i, Alpha, CN, CT )
      FN = 0.5_ReKi*p%RhoW*p%RudAR(i)*UR2*CN
      FT = 0.5_ReKi*p%RhoW*p%RudAR(i)*UR2*CT
      cd = cos(Delta)
      sd = sin(Delta)
      call AddPointLoad( p%RudPos(:,i),  (/ -( 1.0_ReKi - p%tR(i) )*( FN*sd + FT*cd ), FN*cd - FT*sd, 0.0_ReKi /), FB, MB )
      call AddPointLoad( p%RudHPos(:,i), (/ 0.0_ReKi, p%aH(i)*FN*cd, 0.0_ReKi /), FB, MB )

      m%RudAng(i) = Delta
      m%RudCmd(i) = DeltaCmd
      m%RudUR(i)  = uR
      m%RudVR(i)  = vR
      m%RudAoA(i) = Alpha
      m%RudFN(i)  = FN
      m%RudFT(i)  = FT
      m%RudImm(i) = Elev - ( PointZ( p%RudPos(:,i) ) + 0.5_ReKi*p%RudHR(i) )
      call CheckEmerged( m%RudImm(i), 'rudder', i )
   end do

      ! Global axes, about the body reference point (Theory 9)
   Force  = real( matmul( transpose(k%Rg2b), real(FB,R8Ki) ), ReKi )
   Moment = real( matmul( transpose(k%Rg2b), real(MB,R8Ki) ), ReKi )
   m%FHull(1:3) = real( matmul( transpose(k%Rg2b), real(FHB,R8Ki) ), ReKi )
   m%FHull(4:6) = real( matmul( transpose(k%Rg2b), real(MHB,R8Ki) ), ReKi )
   m%FCor(1:3)  = real( matmul( transpose(k%Rg2b), real(FCB,R8Ki) ), ReKi )
   m%FCor(4:6)  = real( matmul( transpose(k%Rg2b), real(MCB,R8Ki) ), ReKi )
   m%FTot(1:3)  = Force
   m%FTot(4:6)  = Moment

   m%U      = k%u
   m%Vm     = k%vm
   m%R      = k%r
   m%Beta   = k%Beta
   m%Hdg    = k%Sig(5)
   m%XErr   = k%Sig(1)
   m%YErr   = k%Sig(3)
   m%PsiCmd = PsiCmd

contains
   !> Adds a force given in OpenFAST body axes at body point Rho to the body-axis force and moment about the body reference point.
   subroutine AddPointLoad( Rho, F, FSum, MSum )
      real(ReKi), intent(in   ) :: Rho(3), F(3)
      real(ReKi), intent(inout) :: FSum(3), MSum(3)
      FSum = FSum + F
      MSum = MSum + cross_product( Rho, F )
   end subroutine AddPointLoad

   real(ReKi) function PointZ( Rho )
      real(ReKi), intent(in) :: Rho(3)
      real(R8Ki) :: Pos(3)
      Pos    = k%PosO + matmul( transpose(k%Rg2b), real(Rho,R8Ki) )
      PointZ = real( Pos(3), ReKi )
   end function PointZ

   !> Local wave orbital velocity in MMG axes and wave elevation at body point Rho (Theory 3.2).
   subroutine Orbital( Rho, o, Elev )
      real(ReKi), intent(in   ) :: Rho(3)
      real(ReKi), intent(  out) :: o(3), Elev
      real(R8Ki)     :: Pos(3), Vw(3)
      real(ReKi)     :: PosR(3)
      real(SiKi)     :: FV(3)
      integer(IntKi) :: InWater

      o    = 0.0_ReKi
      Elev = 0.0_ReKi
      if ( .not. p%HasWaves ) return
      Pos  = k%PosO + matmul( transpose(k%Rg2b), real(Rho,R8Ki) )
      PosR = real( Pos, ReKi )
      Elev = real( WaveField_GetNodeTotalWaveElev( p%WaveField, m%WaveField_m, Time, PosR, ErrStat2, ErrMsg2 ), ReKi )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ); if (ErrStat >= AbortErrLev) return
      call WaveField_GetNodeWaveVel( p%WaveField, m%WaveField_m, Time, PosR, .false., .false., InWater, FV, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ); if (ErrStat >= AbortErrLev) return
      if ( InWater == 0 ) return
      Vw    = real( FV, R8Ki )
      Vw(1) = Vw(1) + real( p%FrameSpd, R8Ki )      ! remove the mean current -FrameSpd
      Vw    = matmul( k%Rg2b, Vw )
      o     = real( (/ Vw(1), -Vw(2), -Vw(3) /), ReKi )
   end subroutine Orbital

   subroutine CheckEmerged( Imm, Name, iProp )
      real(ReKi),     intent(in) :: Imm
      character(*),   intent(in) :: Name
      integer(IntKi), intent(in) :: iProp
      if ( Imm < 0.0_ReKi .and. m%FirstWarn_Emerged ) then
         call SetErrStat( ErrID_Warn, 'Body '//trim(Num2LStr(p%BodyID))//': '//Name//' '//trim(Num2LStr(iProp))// &
                          ' emerges from the water (no ventilation model). Further warnings are suppressed.', ErrStat, ErrMsg, RoutineName )
         m%FirstWarn_Emerged = .false.
      end if
   end subroutine CheckEmerged
end subroutine Vessel_CalcOutput

end module Vessel
