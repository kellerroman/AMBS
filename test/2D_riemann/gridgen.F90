program gridgen

use mod_gridgen
use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none

integer :: imax = 100
integer :: jmax = 100
integer :: kmax = 2
integer, parameter :: nVar   = 5

real(kind=8), parameter :: a2d = 180.0D0 / 3.1415927D0
real(kind=8), parameter :: length = 1.0D0
real(kind=8) :: winkel = 0.0D0

real(kind=8) :: tu,tv,trho,tp
character(len = 100) :: arg

real(kind=8) :: mat(2,2)
real(kind=8) :: temp(3)


integer :: i,j,k


write(*,'(A)') "SIMPLE GRID GEN"
i = 1
DO
   CALL get_command_argument(i, arg)
   IF (LEN_TRIM(arg) == 0) EXIT
   if( i == 1) read(arg,*) imax
   if( i == 2) read(arg,*) winkel
   i = i+1
END DO
write(*,'(A,1X,I0,1X,F5.2)') "AXIAL GRID RESOLUTION:",imax,winkel
winkel = winkel / a2d

mat(1,1) = + cos(winkel)
mat(2,1) = - sin(winkel)
mat(1,2) = + sin(winkel)
mat(2,2) = + cos(winkel)

call add_block(imax-1,jmax-1,kmax-1)
call allocate_blocks(nVar)

do k = 1,kmax
   do j = 1,jmax
      do i = 1,imax
         blocks(1) % xyzs(i,j,k,1) = length/dble(imax-1) * dble(i-1)
   
         blocks(1) % xyzs(i,j,k,2) = length/dble(imax-1) * dble(j-1)
         blocks(1) % xyzs(i,j,k,3) = length/dble(imax-1) * dble(k-1)
         temp = blocks(1) % xyzs(i,j,k,:)
         blocks(1) % xyzs(i,j,k,1) = mat(1,1) * temp(1) + mat(2,1) * temp(2)
         blocks(1) % xyzs(i,j,k,2) = mat(1,2) * temp(1) + mat(2,2) * temp(2)
   
      end do
   end do
end do
do i = 1,imax-1
   do j = 1,jmax-1
     
      if (i < imax/2 .and. j < jmax/2) then
         tu   = 0.8939d0
         tv   = 0.8939d0
         trho = 1.1d0
         tp   = 1.1d0
      else if (i < imax/2 .and. j >= jmax/2) then
         tu   = 0.8939d0
         tv   = 0.d0
         trho = 0.5065d0
         tp   = 0.35d0
      else if (i >= imax/2 .and. j < jmax/2) then
         tu   = 0.d0
         tv   = 0.8939d0
         trho = 0.5065d0
         tp   = 0.35d0
      else if (i >= imax/2 .and. j >= jmax/2) then
         tu   = 0.d0
         tv   = 0.d0
         trho = 1.1d0
         tp   = 1.1d0
      end if
      blocks(1) % vars(i,j,:,1) = trho
      blocks(1) % vars(i,j,:,2) = tu
      blocks(1) % vars(i,j,:,3) = tv
      blocks(1) % vars(i,j,:,4) = 0.0D0
      blocks(1) % vars(i,j,:,5) = tp/(gamma-1.0D0)+ &
                     0.5d0 * (tu*tu+tv*tv)*trho
   end do

end do
blocks(1) % boundary_condition(DIR_WEST) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_EAST) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_SOUTH) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_NORTH) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_FRONT) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_BACK) = BC_SYMMETRY

call write_grid()
call write_xdmf()
write(*,'(A)') "done"


end program
