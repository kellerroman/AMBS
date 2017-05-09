program fftw_test_1d
 use, intrinsic :: iso_c_binding 
implicit none
include 'fftw3.f03'
integer, parameter :: FO = 12345
integer, parameter :: NI = 50
integer, parameter :: NJ = 50

real(8), parameter   :: PI     = 4.0E0_8 * atan(1.0E0_8)

real(8), parameter   :: N_PERIOD    = 7.00E0_8
real(8), parameter   :: AMPL        = 2.00E0_8

real(8), parameter   :: N_PERIOD2   = 1.00E0_8
real(8), parameter   :: AMPL2       = 1.00E0_8

real(8), parameter   :: SIN_FKT     = 2.0E0_8 * PI / dble(NI) * N_PERIOD
real(8), parameter   :: COS_FKT     = 2.0E0_8 * PI / dble(NI) * N_PERIOD2

character(len=*), parameter :: FILE_OUT_NAME = "spektrum.dat"

type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), allocatable :: out(:,:)
real(8), allocatable :: u(:,:)

integer :: i,j
real(8) :: fj
write(*,*) "========== FFT 1D TEST =============="

allocate(out(NI,NJ))
allocate(u(NI,NJ))
do i = 1, NI
   do j = 1, NJ
      fj = 1.0E0_8 / sqrt(sqrt(dble(j)))
      fj = 1.0E0_8 / (1+1*dble(j-1)/dble(nj-1))
!     fj = 1.0E0_8
      u(i,j) = AMPL * sin(dble(j) * SIN_FKT)
!     u(i,j) = AMPL * fj * sin(dble(i) * SIN_FKT * fj )
!     u(i,j) = AMPL      * sin(dble(i) * SIN_FKT * fj )
!     u(i,j) = AMPL * fj * sin(dble(i+j) * SIN_FKT * fj )
!     u(i,j) = AMPL * sin(dble(i) * SIN_FKT) + AMPL2 * sin(dble(j) * COS_FKT)
!     if (abs(i-(NI/2)) < 5 .and. abs(j-(NJ/2)) < 2) then
!        u(i,j) = 1.0E0_8
!     else
!        u(i,j) = 0.0E0_8
!     end if
   end do
end do
open(unit = fo , file = "solution.dat")

do j = 1, NJ
   do i = 1, NI
      write(fo,*) i,j,u(i,j)
   end do
end do
close(fo)

plan = fftw_plan_dft_r2c_2d(NI,NJ        &
                           ,u            &
                           ,out          &
                           ,FFTW_ESTIMATE)
call fftw_execute_dft_r2c(plan, u , out)
call fftw_destroy_plan(plan)

open(unit = fo , file = file_out_name)
do j = 1, NJ/2
   do i = 1, NI
      write(fo,*) i,j,sqrt(real(out(i,j))*real(out(i,j)) + aimag(out(i,j))*aimag(out(i,j)))/(NJ*NI)*2.0E0_8 &
                     ,real(out(i,j)), aimag(out(i,j))
   end do
end do
close(fo)

write(*,*) "========== done ========="
end program fftw_test_1d


