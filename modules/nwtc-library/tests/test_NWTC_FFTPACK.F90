module test_NWTC_FFTPACK

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_FFTPACK
use NWTC_Library
use nwtc_library_test_tools

implicit none

private
public :: test_NWTC_FFTPACK_suite

real(SiKi), parameter :: tol = 1.0e-4_SiKi

contains

subroutine test_NWTC_FFTPACK_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      new_unittest("FFT_roundtrip", test_fft_roundtrip), &
      new_unittest("FFT_forward_known", test_fft_forward_known), &
      new_unittest("FFT_cx_vs_real", test_fft_cx_vs_real), &
      new_unittest("CFFT_forward_known", test_cfft_forward_known), &
      new_unittest("COST_known_values", test_cost_known_values), &
      new_unittest("SINT_known_values", test_sint_known_values), &
      new_unittest("FFT_normalize", test_fft_normalize) &
   ]
end subroutine

! Forward then backward (no normalization) gives x*N
subroutine test_fft_roundtrip(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 16
   real(SiKi) :: x(N), x_orig(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   do i = 1, N
      x(i) = sin(2.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi)) &
            + 0.5_SiKi * cos(4.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi))
   end do
   x_orig = x

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Un-normalized roundtrip: forward then backward = x * N
   do i = 1, N
      call check(error, real(x(i), kind=8), real(N, kind=8) * real(x_orig(i), kind=8), thr=real(tol, kind=8))
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! Forward FFT of known signal should produce expected coefficients
subroutine test_fft_forward_known(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat

   ! Constant signal: all 3.0 => DC = N*3 = 24, all others zero
   x = 3.0_SiKi

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! DC component R(1) = sum of all values = N*3 = 24
   call check(error, real(x(1), kind=8), 24.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   ! All other coefficients should be zero
   call check(error, real(x(2), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return
   call check(error, real(x(3), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return
   call check(error, real(x(N), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   call ExitFFT(fft, ErrStat)
end subroutine

! ApplyFFT_cx should produce same result as ApplyFFT for equivalent input
subroutine test_fft_cx_vs_real(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), x_real(N), x_cx(N)
   complex(SiKi) :: H(N/2+1)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   do i = 1, N
      x(i) = sin(2.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi)) + 2.0_SiKi
   end do

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Get spectral coefficients
   x_real = x
   call ApplyFFT_f(x_real, fft, ErrStat)

   ! Build complex array from real interleaved format
   H(1) = CMPLX(x_real(1), 0.0_SiKi, SiKi)
   do i = 2, N/2
      H(i) = CMPLX(x_real(2*i-2), x_real(2*i-1), SiKi)
   end do
   H(N/2+1) = CMPLX(x_real(N), 0.0_SiKi, SiKi)

   ! Backward via complex interface
   call ApplyFFT_cx(x_cx, H, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Backward via real interface (reuse x_real which has spectral data)
   call ApplyFFT(x_real, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Both should produce identical results
   do i = 1, N
      call check(error, real(x_cx(i), kind=8), real(x_real(i), kind=8), thr=real(tol, kind=8))
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! Complex FFT forward of a pure cosine should produce energy at correct bin
subroutine test_cfft_forward_known(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 16
   complex(SiKi) :: c(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   ! Pure cosine at frequency bin 2: cos(2*pi*2*t/N)
   do i = 1, N
      c(i) = CMPLX(cos(4.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi)), 0.0_SiKi, SiKi)
   end do

   call InitCFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyCFFT_f(c, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! DC should be zero
   call check(error, real(abs(c(1)), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   ! Bin 2 (index 3) should have energy; check it's nonzero
   call check(error, real(abs(c(3)), kind=8) > real(tol, kind=8), .true.)
   if (allocated(error)) return

   ! Bins 1 (index 2) and 3 (index 4) should be zero
   call check(error, real(abs(c(2)), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return
   call check(error, real(abs(c(4)), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   call ExitCFFT(fft, ErrStat)
end subroutine

! Cosine transform of known input vs analytical formula
subroutine test_cost_known_values(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 9  ! must be odd
   real(SiKi) :: x(N), y(N), y_expected(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i, j
   real(SiKi) :: pi_val

   pi_val = real(Pi_D, SiKi)

   ! Input: x(j) = j for j=1..N
   do i = 1, N
      x(i) = real(i, SiKi)
   end do
   y = x

   ! Analytical formula for un-normalized DCT-I:
   ! y(J) = X(1) + (-1)^(J-1)*X(N) + sum_{K=2}^{N-1} 2*X(K)*cos((K-1)*(J-1)*pi/(N-1))
   do j = 1, N
      y_expected(j) = x(1) + ((-1.0_SiKi)**(j-1)) * x(N)
      do i = 2, N-1
         y_expected(j) = y_expected(j) + 2.0_SiKi * x(i) * &
            cos(real(i-1, SiKi) * real(j-1, SiKi) * pi_val / real(N-1, SiKi))
      end do
   end do

   call InitCOST(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyCOST(y, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(y(i), kind=8), real(y_expected(i), kind=8), thr=0.01d0)
      if (allocated(error)) return
   end do

   call ExitCOST(fft, ErrStat)
end subroutine

! Sine transform of known input vs analytical formula
subroutine test_sint_known_values(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 9  ! must be odd
   real(SiKi) :: x(N), y(N), y_expected(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i, j
   real(SiKi) :: pi_val

   pi_val = real(Pi_D, SiKi)

   ! Input: endpoints zero, interior = sin(k*pi/(N-1))
   x(1) = 0.0_SiKi
   do i = 2, N-1
      x(i) = sin(real(i-1, SiKi) * pi_val / real(N-1, SiKi))
   end do
   x(N) = 0.0_SiKi
   y = x

   ! Analytical formula for un-normalized DST-I:
   ! y(J) = sum_{K=2}^{N-1} 2*X(K)*sin((K-1)*(J-1)*pi/(N-1))
   y_expected(1) = 0.0_SiKi
   y_expected(N) = 0.0_SiKi
   do j = 2, N-1
      y_expected(j) = 0.0_SiKi
      do i = 2, N-1
         y_expected(j) = y_expected(j) + 2.0_SiKi * x(i) * &
            sin(real(i-1, SiKi) * real(j-1, SiKi) * pi_val / real(N-1, SiKi))
      end do
   end do

   call InitSINT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplySINT(y, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(y(i), kind=8), real(y_expected(i), kind=8), thr=0.01d0)
      if (allocated(error)) return
   end do

   call ExitSINT(fft, ErrStat)
end subroutine

! FFT with normalization on both: forward(/N) then backward(/N) gives x/N
subroutine test_fft_normalize(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), x_orig(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   do i = 1, N
      x(i) = real(i, SiKi)
   end do
   x_orig = x

   call InitFFT(N, fft, NormalizeIn=.TRUE., ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call ApplyFFT(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Both forward and backward divide by N, so result = x*N / N^2 = x/N
   do i = 1, N
      call check(error, real(x(i), kind=8), real(x_orig(i) / real(N, SiKi), kind=8), thr=real(tol, kind=8))
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

end module
