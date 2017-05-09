program fftw_test_1d
 use, intrinsic :: iso_c_binding 
implicit none
include 'fftw3.f03'
integer, parameter :: FO = 12345
integer, parameter :: NI = 100
real(8), parameter   :: PI     = 4.0E0_8 * atan(1.0E0_8)
real(8), parameter   :: N_PERIOD    = 1.00E0_8
real(8), parameter   :: AMPL        = 5.00E0_8

real(8), parameter   :: N_PERIOD2   = 3.00E0_8
real(8), parameter   :: AMPL2       = 2.00E0_8

real(8), parameter   :: SIN_FKT     = 2.0E0_8 * PI / dble(NI) * N_PERIOD
real(8), parameter   :: COS_FKT     = 2.0E0_8 * PI / dble(NI) * N_PERIOD2

character(len=*), parameter :: FILE_OUT_NAME = "spektrum.dat"

type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), allocatable :: out(:)
real(8), allocatable :: u(:)

integer :: i
write(*,*) "========== FFT 1D TEST =============="

allocate(out(ni))
allocate(u(ni))

do i = 1, NI
   u(i) = AMPL * sin(dble(i) * SIN_FKT)
!  u(i) = AMPL * sin(dble(i) * SIN_FKT) + AMPL2 * cos(dble(i) * COS_FKT)
end do

open(unit = fo , file = "solution.dat")

do i = 1, NI
   write(fo,*) u(i)
end do

close(fo)
plan = fftw_plan_dft_r2c_1d(ni           &
                           ,u            &
                           ,out          &
                           ,FFTW_ESTIMATE)
call fftw_execute_dft_r2c(plan, u , out)
call fftw_destroy_plan(plan)

open(unit = fo , file = file_out_name)

do i = 1, NI
   write(fo,*) real(out(i)), aimag(out(I)), sqrt(real(out(i))*real(out(i)) + aimag(out(i))*aimag(out(i)))/NI*2.0E0_8
end do

close(fo)

write(*,*) "========== done ========="
end program fftw_test_1d


