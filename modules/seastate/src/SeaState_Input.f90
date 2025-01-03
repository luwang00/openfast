!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2021  National Renewable Energy Laboratory
!
!    This file is part of SeaState.
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
module SeaState_Input
   use                              NWTC_Library
   use                              SeaState_Types
   use                              SeaState_Output
   use                              Waves
   use                              NWTC_RandomNumber ! for parameters pRNG_INTRINSIC and pRNG_RANLUX

   implicit                         none

   contains

!====================================================================================================
subroutine SeaSt_ParseInput( InputFileName, OutRootName, defWtrDens, defWtrDpth, defMSL2SWL, FileInfo_In, InputFileData, ErrStat, ErrMsg )
!     This public subroutine reads the input required for SeaState from the file whose name is an
!     input parameter.
!----------------------------------------------------------------------------------------------------

      ! Passed variables
   character(*),                  intent(in   ) :: InputFileName        !< The name of the input file, for putting in echo file.
   character(*),                  intent(in   ) :: OutRootName          !< The rootname of the echo file, possibly opened in this routine
   real(ReKi),                    intent(in   ) :: defWtrDens           !< default value for water density
   real(ReKi),                    intent(in   ) :: defWtrDpth           !< default value for water depth
   real(ReKi),                    intent(in   ) :: defMSL2SWL           !< default value for mean sea level to still water level
   type(FileInfoType),            INTENT(IN   ) :: FileInfo_In          !< The derived type for holding the file information
   type(SeaSt_InputFile),         INTENT(INOUT) :: InputFileData        ! the SeaState input file data
   integer,                       INTENT(  OUT) :: ErrStat              ! returns a non-zero value when an error occurs
   character(*),                  INTENT(  OUT) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   integer                                      :: UnEc                 ! The local unit number for this module's echo file
   character(1024)                              :: EchoFile             ! Name of SeaState echo file
   character(MaxFileInfoLineLen)                :: Line                 ! String to temporarially hold value of read line
   real(ReKi), allocatable                      :: tmpVec1(:), tmpVec2(:) ! Temporary arrays for WAMIT data
   integer, allocatable                         :: tmpArray(:)          ! Temporary array storage of the joint output list
   real(ReKi), allocatable                      :: tmpReArray(:)        ! Temporary array storage of the joint output list
   character(1)                                 :: Line1                ! The first character of an input line
   integer(IntKi)                               :: CurLine              !< Current entry in FileInfo_In%Lines array
   integer(IntKi)                               :: IOS
   integer(IntKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*),  parameter                     :: RoutineName = 'SeaSt_ParseInput'

      ! Initialize local data
   UnEc     = -1
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   InputFileData%Echo = .FALSE.  ! initialize for error handling (cleanup() routine)


   !-------------------------------------------------------------------------------------------------
   ! General settings
   !-------------------------------------------------------------------------------------------------

   CurLine = 3    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2 )
         if (Failed()) return;

   if ( InputFileData%Echo ) then
      EchoFile = trim(OutRootName)//'.ech'
      call OpenEcho ( UnEc, trim(EchoFile), ErrStat2, ErrMsg2 )
         if (Failed())  return;
      write(UnEc, '(A)') 'Echo file for SeaState primary input file: '//trim(InputFileName)
      ! Write the first three lines into the echo file
      write(UnEc, '(A)') trim(FileInfo_In%Lines(1))
      write(UnEc, '(A)') trim(FileInfo_In%Lines(2))

      CurLine = 3
      call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
   endif


   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! WtrDens - Water density.
   call ParseVarWDefault ( FileInfo_In, CurLine, 'WtrDens', InputFileData%WtrDens, defWtrDens, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WtrDpth - Water depth
   call ParseVarWDefault ( FileInfo_In, CurLine, 'WtrDpth', InputFileData%WtrDpth, defWtrDpth, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! MSL2SWL
   call ParseVarWDefault ( FileInfo_In, CurLine, 'MSL2SWL', InputFileData%MSL2SWL, defMSL2SWL, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for Wave Kinematics data grid spatial discretization
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! X_HalfWidth - Half-width of the domain in the X direction.
   call ParseVar( FileInfo_In, CurLine, 'X_HalfWidth', InputFileData%X_HalfWidth, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Y_HalfWidth - Half-width of the domain in the Y direction.
   call ParseVar( FileInfo_In, CurLine, 'Y_HalfWidth', InputFileData%Y_HalfWidth, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Z_Depth - Depth of the domain the Z direction.
   call ParseVarWDefault ( FileInfo_In, CurLine, 'Z_Depth', InputFileData%Z_Depth, InputFileData%WtrDpth+InputFileData%MSL2SWL, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NX - Number of nodes in half of the X-direction domain.
   call ParseVar( FileInfo_In, CurLine, 'NX', InputFileData%NX, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NY - Number of nodes in half of the Y-direction domain.
   call ParseVar( FileInfo_In, CurLine, 'NY', InputFileData%NY, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NZ - Number of nodes in the Z-direction domain.
   call ParseVar( FileInfo_In, CurLine, 'NZ', InputFileData%NZ, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for waves
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! WaveMod - Wave kinematics model switch. and WavePhase (as appropriate)
   InputFileData%Waves%WavePhase = 0.0
   call ParseVar( FileInfo_In, CurLine, 'WaveMod', InputFileData%WaveMod, ErrStat2, ErrMsg2, UnEc )
   if ( ErrStat2 >= AbortErrLev ) then
      ! try to read the line that just failed, as a string this time to see if it's "1P"
      call ParseVar( FileInfo_In, CurLine, 'WaveMod', Line, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return
      
      call Conv2UC( Line )    ! Convert Line to upper case.
      if ( Line(1:2) == '1P' )  then                     ! The user wants to specify the phase in place of a random phase

         InputFileData%WaveMod   = WaveMod_RegularUsrPh                          ! Internally define WaveMod = 10 to mean regular waves with a specified (nonrandom) phase
         
         read (Line(3:),*,IOSTAT=IOS )  InputFileData%Waves%WavePhase
            call CheckIOS ( IOS, "", 'WavePhase', NumType, ErrStat2, ErrMsg2 )
            if (Failed())  return

         InputFileData%Waves%WavePhase = InputFileData%Waves%WavePhase*D2R       ! Convert the phase from degrees to radians

      else                                               ! The user must have specified WaveMod incorrectly.
         ErrStat2 = ErrID_Fatal
         ErrMsg2 = 'WaveMod incorrectly specified in SeaState input file.'
         if (Failed())  return
      end if

   end if

      ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface.
   call ParseVar( FileInfo_In, CurLine, 'WaveStMod', InputFileData%WaveStMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveTMax - Analysis time for incident wave calculations.
   call ParseVar( FileInfo_In, CurLine, 'WaveTMax', InputFileData%Waves%WaveTMax, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDT - Time step for incident wave calculations
   call ParseVar( FileInfo_In, CurLine, 'WaveDT', InputFileData%Waves%WaveDT, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveHs - Significant wave height
   call ParseVar( FileInfo_In, CurLine, 'WaveHs', InputFileData%Waves%WaveHs, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveTp - Peak spectral period.
   call ParseVar( FileInfo_In, CurLine, 'WaveTp', InputFileData%Waves%WaveTp, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WavePkShp - Peak shape parameter.
   call ParseVarWDefault(FileInfo_In, CurLine, 'WavePkShp', InputFileData%Waves%WavePkShp, &
                          WavePkShpDefault( InputFileData%WaveMod, InputFileData%Waves%WaveHs, InputFileData%Waves%WaveTp), ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


      ! WvLowCOff - Low Cut-off frequency or lower frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).
   call ParseVar( FileInfo_In, CurLine, 'WvLowCOff', InputFileData%WvLowCOff, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

     ! WvHiCOff - High Cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s).
   call ParseVar( FileInfo_In, CurLine, 'WvHiCOff', InputFileData%WvHiCOff, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDir - Mean wave heading direction.
   call ParseVar( FileInfo_In, CurLine, 'WaveDir', InputFileData%WaveDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDirMod -  Directional spreading function {0: None, 1: COS2S}       (-) [Used only if WaveMod=2]
   call ParseVar( FileInfo_In, CurLine, 'WaveDirMod', InputFileData%WaveDirMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDirSpread -  Spreading coefficient [only used if WaveMod=2 and WaveDirMod=1]
   call ParseVar( FileInfo_In, CurLine, 'WaveDirSpread', InputFileData%Waves%WaveDirSpread, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveNDir -  The number of wave directions to calculate [must be odd; only used if WaveDirMod=1]
   call ParseVar( FileInfo_In, CurLine, 'WaveNDir', InputFileData%Waves%WaveNDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveDirRange - Full range of the wave directions from WaveDir - WaveDirRange/2 to WaveDir + WaveDirRange/2 (only used if WaveMod=2 and WaveDirMod=1)
   call ParseVar( FileInfo_In, CurLine, 'WaveDirRange', InputFileData%Waves%WaveDirRange, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! Negative values should be treated as positive.
   InputFileData%Waves%WaveDirRange =  abs( InputFileData%Waves%WaveDirRange )


      ! WaveSeed(1)
   call ParseVar( FileInfo_In, CurLine, 'WaveSeed(1)', InputFileData%Waves%RNG%RandSeed(1), ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      !WaveSeed(2)
   call ParseVar( FileInfo_In, CurLine, 'WaveSeed(2)', Line, ErrStat2, ErrMsg2, UnEc )    ! Read into a string and then parse
      if (Failed())  return;

   read (Line,*,IOSTAT=ErrStat2) Line1  ! check the first character to make sure we don't have T/F, which can be interpreted as 1/-1 or 0 in Fortran
   call Conv2UC( Line1 )
   if ( (Line1 == 'T') .OR. (Line1 == 'F') ) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = ' WaveSeed(2): Invalid RNG type.'
      if (Failed())  return;
   endif

   read (Line,*,IOSTAT=ErrStat2) InputFileData%Waves%RNG%RandSeed(2)

   if (ErrStat2 == 0) then ! the user entered a number
      InputFileData%Waves%RNG%RNG_type = "NORMAL"
      InputFileData%Waves%RNG%pRNG = pRNG_INTRINSIC
   else
      InputFileData%Waves%RNG%RandSeed(2) = 0

      InputFileData%Waves%RNG%RNG_type = adjustl( Line )
      call Conv2UC( InputFileData%Waves%RNG%RNG_type )

      if ( InputFileData%Waves%RNG%RNG_type == "RANLUX") then
         InputFileData%Waves%RNG%pRNG = pRNG_RANLUX
      else
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = ' WaveSeed(2): Invalid alternative random number generator.'
         if (Failed())  return;
      endif

   endif


      ! WaveNDAmp - Flag for normally distributed amplitudes.
   call ParseVar( FileInfo_In, CurLine, 'WaveNDAmp', InputFileData%Waves%WaveNDAmp, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvKinFile
   call ParseVar( FileInfo_In, CurLine, 'WvKinFile', InputFileData%Waves%WvKinFile, ErrStat2, ErrMsg2, UnEc, IsPath=.true. )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for 2nd Order Waves
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! WvDiffQTFF     - Second order waves -- difference forces
   call ParseVar( FileInfo_In, CurLine, 'WvDiffQTF', InputFileData%Waves2%WvDiffQTFF, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvSumQTFF      - Second order waves -- sum forces
   call ParseVar( FileInfo_In, CurLine, 'WvSumQTF', InputFileData%Waves2%WvSumQTFF, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvLowCOffD   -- Minimum frequency used in the difference methods (rad/s)              [Only used if DiffQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvLowCOffD', InputFileData%WvLowCOffD, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvHiCOffD   -- Maximum frequency used in the difference methods  (rad/s)              [Only used if DiffQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvHiCOffD', InputFileData%WvHiCOffD, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvLowCOffS   -- Minimum frequency used in the        sum-QTF     (rad/s)              [Only used if  SumQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvLowCOffS', InputFileData%WvLowCOffS, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WvHiCOffS   -- Maximum frequency used in the        sum-QTF      (rad/s)              [Only used if  SumQTF /= 0]
   call ParseVar( FileInfo_In, CurLine, 'WvHiCOffS', InputFileData%WvHiCOffS, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for constrained wave
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

   ! ConstWaveMod - Constrained wave model switch.
   call ParseVar( FileInfo_In, CurLine, 'ConstWaveMod', InputFileData%Waves%ConstWaveMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


      ! CrestHmax - Crest height
   call ParseVar( FileInfo_In, CurLine, 'CrestHmax', InputFileData%Waves%CrestHmax, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CrestTime -Time of the crest
   call ParseVar( FileInfo_In, CurLine, 'CrestTime', InputFileData%Waves%CrestTime, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CrestXi - X-position of the crest
   call ParseVar( FileInfo_In, CurLine, 'CrestXi', InputFileData%Waves%CrestXi, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CrestYi - Y-position of the crest
   call ParseVar( FileInfo_In, CurLine, 'CrestYi', InputFileData%Waves%CrestYi, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for current
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! CurrMod - Current profile model switch
   call ParseVar( FileInfo_In, CurLine, 'CurrMod', InputFileData%Current%CurrMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrSSV0 - Sub-surface current velocity at still water level
   call ParseVar( FileInfo_In, CurLine, 'CurrSSV0', InputFileData%Current%CurrSSV0, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


      ! CurrSSDirChr - Sub-surface current heading direction
   call ParseVar( FileInfo_In, CurLine, 'CurrSSDir', InputFileData%Current%CurrSSDirChr, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   call Conv2UC( InputFileData%Current%CurrSSDirChr )    ! Convert Line to upper case.


      ! CurrNSRef - Near-surface current reference depth.
   call ParseVar( FileInfo_In, CurLine, 'CurrNSRef', InputFileData%Current%CurrNSRef, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrNSV0 - Near-surface current velocity at still water level.
   call ParseVar( FileInfo_In, CurLine, 'CurrNSV0', InputFileData%Current%CurrNSV0, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrNSDir - Near-surface current heading direction.
   call ParseVar( FileInfo_In, CurLine, 'CurrNSDir', InputFileData%Current%CurrNSDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrDIV - Depth-independent current velocity.
   call ParseVar( FileInfo_In, CurLine, 'CurrDIV', InputFileData%Current%CurrDIV, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! CurrDIDir - Depth-independent current heading direction.
   call ParseVar( FileInfo_In, CurLine, 'CurrDIDir', InputFileData%Current%CurrDIDir, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for the MacCamy-Fuchs diffraction model
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! MacCamy-Fuchs member radius
   call ParseVar( FileInfo_In, CurLine, 'MCFD', InputFileData%MCFD, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------
   if ( InputFileData%Echo )   write(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

         ! SeaSum - Whether or not to generate a summary file
   call ParseVar( FileInfo_In, CurLine,  'SeaStSum', InputFileData%SeaStSum, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

         ! OutSwtch - Specify how to write to an output file
   call ParseVar( FileInfo_In, CurLine, 'OutSwtch', InputFileData%OutSwtch, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

        ! OutFmt - Format for numerical outputs
   call ParseVar( FileInfo_In, CurLine, 'OutFmt', InputFileData%OutFmt, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

         ! OutSFmt - Format for output column headers
   call ParseVar( FileInfo_In, CurLine, 'OutSFmt', InputFileData%OutSFmt, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

            ! NWaveElev - Number of Wave elevations to output
   call ParseVar( FileInfo_In, CurLine, 'NWaveElev', InputFileData%NWaveElev, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


      ! This check is needed here instead of being located in SeaStateInput_ProcessInputData() because
      ! we need to allocate arrays.  If _GetInput() was skipped, then these array would already have
      ! been allocated and populated.

   if ( InputFileData%NWaveElev < 0 .OR. InputFileData%NWaveElev > 9 ) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'NWaveElev must be greater than or equal to zero and less than 10.'
      if (Failed())  return;
   end if

      ! allocate space for the output location arrays:
   call AllocAry( InputFileData%WaveElevxi, InputFileData%NWaveElev, 'WaveElevxi' , ErrStat2, ErrMsg2);  if (Failed())  return;
   call AllocAry( InputFileData%WaveElevyi, InputFileData%NWaveElev, 'WaveElevyi' , ErrStat2, ErrMsg2);  if (Failed())  return;

      ! WaveElevxi
   call ParseAry ( FileInfo_In, CurLine, 'WaveElevxi.', InputFileData%WaveElevxi, InputFileData%NWaveElev, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveElevyi
   call ParseAry ( FileInfo_In, CurLine, 'WaveElevyi.', InputFileData%WaveElevyi, InputFileData%NWaveElev, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! NWaveKin
   call ParseVar( FileInfo_In, CurLine, 'NWaveKin', InputFileData%NWaveKin, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;


      ! This check is needed here instead of being located in SeaStateInput_ProcessInputData() because
      ! we need to allocate arrays.  If _GetInput() was skipped, then these array would already have
      ! been allocated and populated.

   if ( InputFileData%NWaveKin < 0 .OR. InputFileData%NWaveKin > 9 ) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'NWaveKin must be greater than or equal to zero and less than 10.'
      if (Failed())  return;
   end if

      ! allocate space for the output location arrays:
   call AllocAry( InputFileData%WaveKinxi, InputFileData%NWaveKin, 'WaveKinxi' , ErrStat2, ErrMsg2);  if (Failed())  return;
   call AllocAry( InputFileData%WaveKinyi, InputFileData%NWaveKin, 'WaveKinyi' , ErrStat2, ErrMsg2);  if (Failed())  return;
   call AllocAry( InputFileData%WaveKinzi, InputFileData%NWaveKin, 'WaveKinzi' , ErrStat2, ErrMsg2);  if (Failed())  return;

      ! WaveKinxi
   call ParseAry ( FileInfo_In, CurLine, 'WaveKinxi.', InputFileData%WaveKinxi, InputFileData%NWaveKin, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveKinyi
   call ParseAry ( FileInfo_In, CurLine, 'WaveKinyi.', InputFileData%WaveKinyi, InputFileData%NWaveKin, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

      ! WaveKinzi
   call ParseAry ( FileInfo_In, CurLine, 'WaveKinzi.', InputFileData%WaveKinzi, InputFileData%NWaveKin, ErrStat2, ErrMsg2, UnEc )
      if (Failed())  return;

   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT CHANNELS
   !-------------------------------------------------------------------------------------------------

   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(CurLine))    ! Write section break to echo
   CurLine = CurLine + 1

      ! OutList - list of requested parameters to output to a file
   call AllocAry( InputFileData%OutList, MaxOutPts, 'InputFileData%OutList', ErrStat2, ErrMsg2 )
      if (Failed())  return;

   call ReadOutputListFromFileInfo( FileInfo_In, CurLine, InputFileData%OutList, InputFileData%NumOuts, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return;

contains
   !..............................
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      if (allocated(tmpArray  )) deallocate(tmpArray  )
      if (allocated(tmpReArray)) deallocate(tmpReArray)
      if (allocated(tmpVec1   )) deallocate(tmpVec1   )
      if (allocated(tmpVec2   )) deallocate(tmpVec2   )
         ! Cleanup the Echo file and global variables
      if (UnEc > 0)  close ( UnEc )
   end subroutine Cleanup

end subroutine SeaSt_ParseInput

!====================================================================================================
subroutine SeaStateInput_ProcessInitData( InitInp, p, InputFileData, ErrStat, ErrMsg )
!     This private subroutine verifies the input required for HydroDyn is correctly specified.
!----------------------------------------------------------------------------------------------------


      ! Passed variables

   type(SeaSt_InitInputType),     intent( in    )   :: InitInp              ! the SeaState data
   type(SeaSt_ParameterType),     intent( inout )   :: p                    ! the SeaState parameter data
   type(SeaSt_InputFile),         intent( inout )   :: InputFileData        ! the SeaState input file data
   integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs
   character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   integer                                          :: I, count             ! Generic loop counter index
   integer                                          :: J                    ! Generic loop counter index
   integer                                          :: K                    ! Generic loop counter index
   character(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   real(ReKi)                                       :: xpos, ypos, zpos
   real(SiKi)                                       :: TmpFreq

   integer(IntKi)                                   :: ErrStat2, IOS
   character(ErrMsgLen)                             :: ErrMsg2
   character(*), parameter                          :: RoutineName = 'SeaStateInput_ProcessInitData'
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !-------------------------------------------------------------------------
   ! Check environmental conditions
   !-------------------------------------------------------------------------


      ! WtrDens - Water density.
   if ( InputFileData%WtrDens < 0.0 )  then
      call SetErrStat( ErrID_Fatal,'WtrDens must not be negative.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! WtrDpth - Water depth
   if ( InputFileData%WtrDpth + InputFileData%MSL2SWL <= 0.0 )  then
      call SetErrStat( ErrID_Fatal,'WtrDpth + MSL2SWL must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! X_HalfWidth - Half-width of the domain in the X direction (m)
   if ( InputFileData%X_HalfWidth <= 0.0_ReKi ) then
      call SetErrStat( ErrID_Fatal,'X_HalfWidth must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if

       ! Y_HalfWidth - Half-width of the domain in the Y direction (m)
   if ( InputFileData%Y_HalfWidth <= 0.0_ReKi ) then
      call SetErrStat( ErrID_Fatal,'Y_HalfWidth must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if

       ! Z_Depth - Depth of the domain the Z direction (m)
   if ( ( InputFileData%Z_Depth <= 0.0_ReKi ) .or. ( InputFileData%Z_Depth > InputFileData%WtrDpth + InputFileData%MSL2SWL ) ) then
      call SetErrStat( ErrID_Fatal,'Z_Depth must be greater than zero and less than or equal to the WtrDpth + MSL2SWL.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! NX - Number of nodes in half of the X-direction domain
   if ( InputFileData%NX < 2 ) then
      call SetErrStat( ErrID_Fatal,'NX must be greater than or equal to 2.',ErrStat,ErrMsg,RoutineName)
      return
   end if

         ! NY - Number of nodes in half of the Y-direction domain
   if ( InputFileData%NY < 2 ) then
      call SetErrStat( ErrID_Fatal,'NY must be greater than or equal to 2.',ErrStat,ErrMsg,RoutineName)
      return
   end if

         ! NZ - Number of nodes in the Z-direction domain
   if ( InputFileData%NZ < 2 ) then
      call SetErrStat( ErrID_Fatal,'NZ must be greater than or equal to 2.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! WaveMod - Wave kinematics model switch.

   SELECT CASE(InputFileData%WaveMod)
      CASE(WaveMod_None)
      CASE(WaveMod_Regular)
      CASE(WaveMod_RegularUsrPh)
      CASE(WaveMod_JONSWAP)
      CASE(WaveMod_WhiteNoise)
      CASE(WaveMod_UserSpctrm)
      CASE(WaveMod_ExtElev)
      CASE(WaveMod_ExtFull)
      CASE(WaveMod_UserFreq)
      CASE DEFAULT
         call SetErrStat( ErrID_Fatal,'WaveMod must be 0, 1, 1P#, 2, 3, 4, 5, 6, or 7',ErrStat,ErrMsg,RoutineName)
         return
   END SELECT



   ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface.
   IF ( InputFileData%WaveMod == WaveMod_None ) THEN
      InputFileData%WaveStMod = 0_IntKi
   ELSEIF ( InputFileData%WaveMod == WaveMod_ExtFull ) THEN
      IF ( (InputFileData%WaveStMod /= 0) .AND. (InputFileData%WaveStMod /= 1) .AND. &
                                                (InputFileData%WaveStMod /= 3) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WaveStMod must be 0, 1, or 3 when WaveMod = 6.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   ELSE
      IF ( (InputFileData%WaveStMod /= 0) .AND. (InputFileData%WaveStMod /= 1) .AND. &
           (InputFileData%WaveStMod /= 2) .AND. (InputFileData%WaveStMod /= 3) ) THEN
         CALL SetErrStat( ErrID_Fatal,'WaveStMod must be 0, 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF
   

      ! WaveTMax - Analysis time for incident wave calculations.

   if ( InputFileData%WaveMod == WaveMod_None )  then   ! .TRUE if we have incident waves.

      ! TODO: Issue warning if WaveTMax was not already 0.0 in this case.
      ! Setting WaveTMax = 0 breaks interpolation. Should probably set it to just TMax instead.
      if ( .NOT. EqualRealNos(InputFileData%Waves%WaveTMax, InitInp%TMax) ) then
         call WrScr( '  Setting WaveTMax to TMax since WaveMod = 0' )
         InputFileData%Waves%WaveTMax = InitInp%TMax
      end if
      if ( .NOT. EqualRealNos(InputFileData%WaveDir, 0.0_SiKi) ) then
         call WrScr( '  Setting WaveDir to 0.0 since WaveMod = 0' )
         InputFileData%WaveDir = 0.0
      end if
   elseif ( InputFileData%WaveMod == WaveMod_ExtElev ) then   ! User wave elevation file reading in
      if (InitInp%TMax > InputFileData%Waves%WaveTMax ) then
         call SetErrstat( ErrID_Fatal, '  WaveTMax must be larger than the simulation time for user wave elevations (WaveMod == 5).',ErrStat,ErrMsg,RoutineName)
         return
      end if
   else
      if (InitInp%TMax > InputFileData%Waves%WaveTMax ) then
         call WrScr( '  WaveTMax is less then the simulation time.  Wave data will repeat every WaveTMax seconds.')
      end if
   end if


      ! WaveDT - Time step for incident wave calculations

   if ( InputFileData%WaveMod /= WaveMod_None )  then   ! .TRUE if we have incident waves.

      if ( InputFileData%Waves%WaveDT <= 0.0 )  then
         call SetErrStat( ErrID_Fatal,'WaveDT must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      ! When waveMod = 0, should also set WaveDT to InitInp%TMax to keep interpolation working.
      ! Essentially just two time steps, t=0 and t=TMax
      !InputFileData%Waves%WaveDT = 0.0
      InputFileData%Waves%WaveDT  = InitInp%TMax

   end if


       ! WaveHs - Significant wave height
   if ( InputFileData%WaveMod == WaveMod_Regular      .OR. &
        InputFileData%WaveMod == WaveMod_RegularUsrPh .OR. &
        InputFileData%WaveMod == WaveMod_JONSWAP      .OR. &
        InputFileData%WaveMod == WaveMod_WhiteNoise    ) then

      if ( InputFileData%Waves%WaveHs <= 0.0 )  then
         call SetErrStat( ErrID_Fatal,'WaveHs must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         return
      end if
      
   end if


      ! WaveTp - Peak spectral period.
   if ( InputFileData%Waves%WaveTp <= 0.0 )  then
      call SetErrStat( ErrID_Fatal,'WaveTp must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if



       ! WavePkShp - Peak shape parameter
   if ( ( InputFileData%Waves%WavePkShp < 1.0 ) .OR. ( InputFileData%Waves%WavePkShp > 7.0 ) )  then
      call SetErrStat( ErrID_Fatal,'WavePkShp must be greater than or equal to 1 and less than or equal to 7.',ErrStat,ErrMsg,RoutineName)
      return
   end if


      ! WvLowCOff and WvHiCOff - Wave Cut-off frequency

   if ( InputFileData%WvLowCOff < 0 ) then
      call SetErrStat( ErrID_Fatal,'WvLowCOff must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! Threshold upper cut-off based on sampling rate
   if ( EqualRealNos(InputFileData%Waves%WaveDT, 0.0_DbKi) ) then
      InputFileData%WvHiCOff = 10000.0;  ! This is not going to be used because WaveDT is zero.
   else
      TmpFreq = REAL( Pi/InputFileData%Waves%WaveDT,SiKi)
      if ( InputFileData%WvHiCOff > TmpFreq ) then
         InputFileData%WvHiCOff =  TmpFreq
         call SetErrStat( ErrID_Info,'WvHiCOff adjusted to '//trim(num2lstr(TmpFreq))//' rad/s, based on WaveDT.',ErrStat,ErrMsg,RoutineName)
      end if
   end if

   if (InputFileData%WaveMod == WaveMod_JONSWAP    .or. &
       InputFileData%WaveMod == WaveMod_WhiteNoise .or. &
       InputFileData%WaveMod == WaveMod_UserSpctrm .or. &
       InputFileData%WaveMod == WaveMod_ExtElev    .or. &
       InputFileData%WaveMod == WaveMod_UserFreq   ) then
       
      if ( InputFileData%WvLowCOff >= InputFileData%WvHiCOff ) then
         call SetErrSTat( ErrID_Fatal,'WvLowCOff must be less than WvHiCOff.',ErrStat,ErrMsg,RoutineName)
         return
      end if
   else
      ! overwrite these so that ALL frequencies are allowed (otherwise we might exclude frequencies with WaveMod = WaveMod_Regular or WaveMod_RegularUsrPh)
      InputFileData%WvLowCOff = -HUGE(InputFileData%WvLowCOff)
      InputFileData%WvHiCOff  =  HUGE(InputFileData%WvHiCOff )
   end if
   
      ! WaveDir - Wave heading direction.

   if ( ( InputFileData%WaveMod /= WaveMod_None ) .AND. ( InputFileData%WaveMod /= WaveMod_ExtFull ) )  then   ! .TRUE if we have incident waves, but not user input wave data.

      if ( ( InputFileData%WaveDir <= -180.0 ) .OR. ( InputFileData%WaveDir > 180.0 ) )  then
         call SetErrStat( ErrID_Fatal,'WaveDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%WaveDir = 0.0

   end if


      ! Multi-directional waves

      ! Check the WaveDirMod value
   if ( InputFileData%WaveDirMod /= WaveDirMod_None .AND. InputFileData%WaveDirMod /= WaveDirMod_COS2S ) then
      call SetErrStat( ErrID_Fatal,'WaveDirMod must be either 0 (No spreading) or 1 (COS2S spreading function)',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! Check if we are doing multidirectional waves or not.
      ! We can only use multi directional waves on WaveMod=2,3,4
   InputFileData%WaveMultiDir = .FALSE.         ! Set flag to false to start
   IF (InputFileData%WaveDirMod == WaveDirMod_COS2S ) THEN
      if ( InputFileData%WaveMod == WaveMod_JONSWAP .OR. InputFileData%WaveMod == WaveMod_WhiteNoise .OR. InputFileData%WaveMod == WaveMod_UserSpctrm ) then
         InputFileData%WaveMultiDir = .TRUE.
      else
         call SetErrStat( ErrID_Warn,'WaveDirMod unused unless WaveMod == 2, 3, or 4.  Ignoring WaveDirMod.',ErrStat,ErrMsg,RoutineName)
      end if
   ENDIF


      !  Check to see if the for some reason the wave direction spreading range is set to zero.  If it is,
      !  we don't have any spreading, so we will turn off the multidirectional waves.
   if ( InputFileData%WaveMultiDir .AND. EqualRealNos( InputFileData%Waves%WaveDirRange, 0.0_SiKi ) ) then
      call SetErrStat( ErrID_Warn,' WaveDirRange set to zero, so multidirectional waves are turned off.',ErrStat,ErrMsg,RoutineName)
      InputFileData%WaveMultiDir = .FALSE.
   ENDIF



      ! We check the following only if we set WaveMultiDir to true, otherwise ignore them and set them to zero
   if ( InputFileData%WaveMultiDir ) then

         ! Check WaveDirSpread
      if ( InputFileData%Waves%WaveDirSpread <= 0.0 ) then

         call SetErrStat( ErrID_Fatal,'WaveDirSpread cannot negative or zero.',ErrStat,ErrMsg,RoutineName)
         return

      ENDIF


         ! Check that the number of wave directions is a positive odd number.
         !     -> If it is less than 0, error out.
         !     -> If it is even, we will increment it by 1.
      if ( InputFileData%Waves%WaveNDir <= 0_IntKi ) then
         call SetErrStat( ErrID_Fatal,' WaveNDir must be an odd number greater than 0.',ErrStat,ErrMsg,RoutineName)
         return
      ENDIF

         ! Check that the value for WaveNDir is odd
      if ( MODULO( InputFileData%Waves%WaveNDir, 2_IntKi) == 0_IntKi ) then
         InputFileData%Waves%WaveNDir  = InputFileData%Waves%WaveNDir + 1
         call SetErrStat( ErrID_Warn,'WaveNDir must be odd.  Changing the value to '//Num2LStr(InputFileData%Waves%WaveNDir),ErrStat,ErrMsg,RoutineName)
      ENDIF

         ! Now check that the WaveDirRange is less than 360 degrees (not sure why we would want that)
      if ( InputFileData%Waves%WaveDirRange > 360.0_ReKi ) then
         call SetErrStat( ErrID_Fatal,' WaveDirRange should be less than a full circle.',ErrStat,ErrMsg,RoutineName)
      ENDIF

   else  ! Set everything to zero if we aren't going to use it

      InputFileData%Waves%WaveNDir        = 1         ! Only one direction set -- this shouldn't get used later anyhow
      InputFileData%Waves%WaveDirRange    = PiBy2     ! This is so that the constant C=1 in the COS2S function (it shouldn't get called, but in case it does)
      InputFileData%Waves%WaveDirSpread   = 0.0

   end if


      ! WvKinFile

   if ( InputFileData%WaveMod == WaveMod_ExtElev .OR. InputFileData%WaveMod == WaveMod_ExtFull .OR. InputFileData%WaveMod == WaveMod_UserFreq) then      ! .TRUE if we are to read user-supplied wave elevation or wave kinematics file(s).

      if ( LEN_TRIM( InputFileData%Waves%WvKinFile ) == 0 )  then
         call SetErrStat( ErrID_Fatal,'WvKinFile must not be an empty string.',ErrStat,ErrMsg,RoutineName)
         return
      end if

      if ( PathIsRelative( InputFileData%Waves%WvKinFile ) ) then
         call GetPath( TRIM(InitInp%InputFile), TmpPath )
         InputFileData%Waves%WvKinFile    = TRIM(TmpPath)//TRIM(InputFileData%Waves%WvKinFile)
      end if
   
   end if


      !-------------------------------------------------------------------------
      ! Check 2nd Order Waves section
      !-------------------------------------------------------------------------


      ! Difference frequency cutoffs

      ! WvLowCOffD and WvHiCOffD - Wave Cut-off frequency
   if ( InputFileData%WvLowCOffD < 0 ) then
      call SetErrStat( ErrID_Fatal,'WvLowCOffD must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! Check that the order given makes sense.
   if ( InputFileData%WvLowCOffD >= InputFileData%WvHiCOffD ) then
      call SetErrStat( ErrID_Fatal,'WvLowCOffD must be less than WvHiCOffD.',ErrStat,ErrMsg,RoutineName)
      return
   end if


      ! Sum frequency cutoffs

      ! WvLowCOffS and WvHiCOffD - Wave Cut-off frequency
   if ( InputFileData%WvLowCOffS < 0 ) then
      call SetErrStat( ErrID_Fatal,'WvLowCOffS must be greater than or equal to zero.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      ! Check that the order given makes sense.
   if ( InputFileData%WvLowCOffS >= InputFileData%WvHiCOffS ) then
      call SetErrStat( ErrID_Fatal,'WvLowCOffS must be less than WvHiCOffS.',ErrStat,ErrMsg,RoutineName)
      return
   end if

      !-------------------------------------------------------------------------
      ! Check Constrained Waves section
      !-------------------------------------------------------------------------

      ! ConstWaveMod
   select case(InputFileData%Waves%ConstWaveMod)
      case(ConstWaveMod_None)          ! 0
      case(ConstWaveMod_CrestElev)     ! 1
      case(ConstWaveMod_Peak2Trough)   ! 2
      case default
         call SetErrStat( ErrID_Fatal,'ConstWaveMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName)
         return
   end select
   
      ! CrestHmax
   IF ( ( InputFileData%WaveMod == WaveMod_JONSWAP ) .AND. ( InputFileData%Waves%ConstWaveMod /= ConstWaveMod_None ) .AND. &
        ( InputFileData%Waves%CrestHmax < InputFileData%Waves%WaveHs ) ) THEN
      call SetErrStat( ErrID_Fatal,'CrestHmax must be larger than WaveHs.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   
      !-------------------------------------------------------------------------
      ! Check Current section
      !-------------------------------------------------------------------------


      ! CurrMod - Current profile model switch

   if ( ( InputFileData%Current%CurrMod /= 0 ) .AND. ( InputFileData%Current%CurrMod /= 1 ) .AND. ( InputFileData%Current%CurrMod /= 2 ) )  then
      call SetErrStat( ErrID_Fatal,'CurrMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName)
      return
   end if

   if ( ( InputFileData%Current%CurrMod /= 0 ) .AND. ( InputFileData%WaveMod == WaveMod_ExtFull ) )  then
      call SetErrStat( ErrID_Fatal,'CurrMod must be set to 0 when WaveMod is set to 6: user-input wave data.',ErrStat,ErrMsg,RoutineName)
      return
   end if


      ! CurrSSV0 - Sub-surface current velocity at still water level

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.

      if ( InputFileData%Current%CurrSSV0 < 0.0 )  then
         call SetErrStat( ErrID_Fatal,'CurrSSV0 must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%Current%CurrSSV0 = 0.0

   end if


      ! CurrSSDirChr - Sub-surface current heading direction

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.


      if ( TRIM(InputFileData%Current%CurrSSDirChr) == 'DEFAULT' )  then   ! .TRUE. when one wants to use the default value of codirectionality between sub-surface current and incident wave propogation heading directions.

         if ( InputFileData%WaveMod == WaveMod_None ) then
            call SetErrStat( ErrID_Fatal,'CurrSSDir must not be set to ''DEFAULT'' when WaveMod is set to 0.',ErrStat,ErrMsg,RoutineName)
            return
         end if

         InputFileData%Current%CurrSSDir = InputFileData%WaveDir

      else                                   ! The input must have been specified numerically.

         read (InputFileData%Current%CurrSSDirChr,*,IOSTAT=IOS)  InputFileData%Current%CurrSSDir
            call CheckIOS ( IOS, "", 'CurrSSDir', NumType, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if ( ErrStat >= AbortErrLev ) return

         if ( ( InputFileData%Current%CurrSSDir <= -180.0 ) .OR. ( InputFileData%Current%CurrSSDir > 180.0 ) )  then
            call SetErrStat( ErrID_Fatal,'CurrSSDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
            return
         end if

      end if


   else

      InputFileData%Current%CurrSSDir = 0.0

   end if


      ! CurrNSRef - Near-surface current reference depth.

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.

      if ( InputFileData%Current%CurrNSRef <= 0.0 ) then
         call SetErrStat( ErrID_Fatal,'CurrNSRef must be greater than zero.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%Current%CurrNSRef = 0.0

   end if



        ! CurrNSV0 - Near-surface current velocity at still water level.

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.

      if ( InputFileData%Current%CurrNSV0 < 0.0 ) then
         call SetErrStat( ErrID_Fatal,'CurrNSV0 must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%Current%CurrNSV0 = 0.0

   end if


      ! CurrNSDir - Near-surface current heading direction.

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.

      if ( ( InputFileData%Current%CurrNSDir <= -180.0 ) .OR. ( InputFileData%Current%CurrNSDir > 180.0 ) )  then
         call SetErrStat( ErrID_Fatal,'CurrNSDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%Current%CurrNSDir = 0.0

   end if


      ! CurrDIV - Depth-independent current velocity.

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.

      if ( InputFileData%Current%CurrDIV < 0.0 ) then
         call SetErrStat( ErrID_Fatal,'CurrDIV must not be less than zero.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%Current%CurrDIV = 0.0

   end if


      ! CurrDIDir - Depth-independent current heading direction.

   if ( InputFileData%Current%CurrMod == 1 )  then  ! .TRUE if we have standard current.

      if ( ( InputFileData%Current%CurrDIDir <= -180.0 ) .OR. ( InputFileData%Current%CurrDIDir > 180.0 ) ) then
         call SetErrStat( ErrID_Fatal,'CurrDIDir must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
         return
      end if

   else

      InputFileData%Current%CurrDIDir = 0.0

   end if

   !-------------------------------------------------------------------------------------------------
   ! Data section for MacCamy-Fuchs diffraction model
   !-------------------------------------------------------------------------------------------------
   IF ( InputFileData%WaveMod == WaveMod_None .OR. InputFileData%WaveMod == WaveMod_ExtFull ) THEN
      IF ( InputFileData%MCFD > 0.0_SiKi ) THEN
         CALL SetErrStat( ErrID_Fatal,' The MacCamy-Fuchs diffraction model is not compatible with WaveMod = 0 or 6. Need to set MCFD to 0.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for OUTPUT
   !-------------------------------------------------------------------------------------------------


      ! OutSwtch - output file switch

   if ( InputFileData%OutSwtch /= 1 .AND. InputFileData%OutSwtch /= 2 .AND. InputFileData%OutSwtch /= 3 ) then
      call SetErrStat( ErrID_Fatal,'OutSwitch must be set to 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
      return
   end if

   !InputFileData%OutFmt
   !InputFileData%OutSFmt

   ! Shift from MSL to SWL coordinate system
   InputFileData%WaveKinzi(:) = InputFileData%WaveKinzi(:) - InputFileData%MSL2SWL


   !----------------------------------------------------------
   ! Populate data in sub-types from parent or other module types
   !----------------------------------------------------------

      ! Current
         ! For wave kinematic calculations, the effective water depth is the user input water depth (positive valued) + MSL2SWL (positive when SWL is above MSL).
      InputFileData%Current%EffWtrDpth  = InputFileData%WtrDpth + InputFileData%MSL2SWL ! adjusted for the MSL2SWL.


      ! Waves
      InputFileData%Waves%Gravity      = InitInp%Gravity
         ! For wave kinematic calculations, the effective water depth is the user input water depth (positive valued) + MSL2SWL (positive when SWL is above MSL).




!TODO: This is now set with the grid points? GJH 7/11/21
      
      p%NGrid(1) = InputFileData%NX*2-1
      p%NGrid(2) = InputFileData%NY*2-1
      p%NGrid(3) = InputFileData%NZ
      p%NGridPts = p%NGrid(1) * p%NGrid(2) * p%NGrid(3)
      InputFileData%Waves%NGrid = p%NGrid
      InputFileData%Current%NGridPts = p%NGridPts

      call AllocAry( InputFileData%Current%WaveKinGridzi, p%NGridPts, 'WaveKinGridzi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if ( ErrStat >= AbortErrLev ) return


         ! Establish the number and locations where the wave kinematics will be computed
      InputFileData%Waves%NWaveKinGrid   = p%NGridPts                          ! Number of grid points where the incident wave kinematics will be computed (-)
      InputFileData%Waves%NWaveElevGrid  = p%NGrid(1)*p%NGrid(2)               ! Number of XY grid points where the wave elevations are computed
      
      if ( InputFileData%Waves%NWaveElevGrid < 0 ) then
         call SetErrStat( ErrID_Fatal,'Number of nodes in the spatial discretization ('//trim(num2lstr(InputFileData%Waves%NWaveElevGrid))//') must not be negative.',ErrStat,ErrMsg,RoutineName)
         return
      end if
      
      call AllocAry( InputFileData%Waves%WaveKinGridxi, p%NGridPts, 'WaveKinGridxi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call AllocAry( InputFileData%Waves%WaveKinGridyi, p%NGridPts, 'WaveKinGridyi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call AllocAry( InputFileData%Waves%WaveKinGridzi, p%NGridPts, 'WaveKinGridzi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if ( ErrStat >= AbortErrLev ) return
      
      ! Generate grid points
      p%deltaGrid(1) = InputFileData%X_HalfWidth/(InputFileData%NX-1)
      p%deltaGrid(2)= InputFileData%Y_HalfWidth/(InputFileData%NY-1)
      p%deltaGrid(3) = PI / ( 2*(InputFileData%NZ-1) )
      count = 1
      do k = 0, p%NGrid(3) - 1
         zpos = - ( 1.0 - cos( real((p%NGrid(3) - 1) - k, ReKi) * p%deltaGrid(3) )  ) * InputFileData%Z_Depth
         do j = 0, p%NGrid(2)-1
            ypos = -InputFileData%Y_HalfWidth + p%deltaGrid(2)*j
            do i= 0, p%NGrid(1)-1
               xpos = -InputFileData%X_HalfWidth + p%deltaGrid(1)*i
               InputFileData%Waves%WaveKinGridxi(count)      = xpos   ! xi-coordinates for points where the incident wave kinematics will be computed;
               InputFileData%Waves%WaveKinGridyi(count)      = ypos   ! yi-coordinates for points where the incident wave kinematics will be computed;

               InputFileData%Waves%WaveKinGridzi(count)      = zpos   ! zi-coordinates for points where the incident wave kinematics will be computed;
               InputFileData%Current%WaveKinGridzi(count) = InputFileData%Waves%WaveKinGridzi(count)

               !if ( k == 0 ) then
               !   InputFileData%Waves%WaveElevGridxi(count)      = xpos   ! xi-coordinates for points where the incident wave kinematics will be computed;
               !   InputFileData%Waves%WaveElevGridyi(count)      = ypos   ! yi-coordinates for points where the incident wave kinematics will be computed;
               !end if
               count = count + 1
            end do
         end do
      end do

      ! Waves2

            ! If we are using the Waves module, the node information must be copied over.
      InputFileData%Waves2%NWaveKinGrid   = InputFileData%Waves%NWaveKinGrid                          ! Number of points where the incident wave kinematics will be computed (-)
      if ( InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF ) then
         InputFileData%Waves2%Gravity       = InitInp%Gravity
         InputFileData%Waves2%NGrid         = p%NGrid
         InputFileData%Waves2%NWaveElevGrid = InputFileData%Waves%NWaveElevGrid

         call AllocAry( InputFileData%Waves2%WaveKinGridxi, p%NGridPts, 'WaveKinGridxi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         call AllocAry( InputFileData%Waves2%WaveKinGridyi, p%NGridPts, 'WaveKinGridyi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         call AllocAry( InputFileData%Waves2%WaveKinGridzi, p%NGridPts, 'WaveKinGridzi' , ErrStat2, ErrMsg2);  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if ( ErrStat >= AbortErrLev ) return

         InputFileData%Waves2%WaveKinGridxi  = InputFileData%Waves%WaveKinGridxi
         InputFileData%Waves2%WaveKinGridyi  = InputFileData%Waves%WaveKinGridyi
         InputFileData%Waves2%WaveKinGridzi  = InputFileData%Waves%WaveKinGridzi
      ENDIF

      
   !------------------------------------------------------------
   ! Allocate the WaveFieldType to store wave field information
   !------------------------------------------------------------
   ALLOCATE(p%WaveField, STAT=ErrStat2)
   IF (ErrStat2 /=0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating WaveField.",ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
         
   p%WaveField%WtrDpth      = InputFileData%WtrDpth
   p%WaveField%MSL2SWL      = InputFileData%MSL2SWL
   p%WaveField%EffWtrDpth   = InputFileData%WtrDpth + InputFileData%MSL2SWL
   
   p%WaveField%WaveMod      = InputFileData%WaveMod
   p%WaveField%WaveStMod    = InputFileData%WaveStMod
   p%WaveField%WtrDens      = InputFileData%WtrDens     ! may have overwritten default InitInp
   p%WaveField%RhoXg        = p%WaveField%WtrDens*InitInp%Gravity               ! For WAMIT and WAMIT2
   p%WaveField%WaveDir      = InputFileData%WaveDir
   p%WaveField%WaveMultiDir = InputFileData%WaveMultiDir
   p%WaveField%MCFD         = InputFileData%MCFD

   p%WaveField%WvLowCOff    =  InputFileData%WvLowCOff
   p%WaveField%WvHiCOff     =  InputFileData%WvHiCOff
   p%WaveField%WvLowCOffD   =  InputFileData%WvLowCOffD
   p%WaveField%WvHiCOffD    =  InputFileData%WvHiCOffD
   p%WaveField%WvLowCOffS   =  InputFileData%WvLowCOffS
   p%WaveField%WvHiCOffS    =  InputFileData%WvHiCOffS
   p%WaveField%WaveDOmega   =  InputFileData%WaveDOmega          ! For WAMIT and WAMIT2, FIT
   
      
end subroutine SeaStateInput_ProcessInitData

end module SeaState_Input
