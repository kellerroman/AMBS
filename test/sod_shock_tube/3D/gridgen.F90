program gridgen

use cgns

implicit none

integer :: imax = 100
integer :: jmax = 100
integer :: kmax = 100
integer, parameter :: Dimen = 3
integer, parameter :: nVar   = 4

real(kind=8), parameter :: a2d = 180.0D0 / 3.1415927D0
real(kind=8), parameter :: length = 1.0D0
real(kind=8) :: winkel = 0.0D0


integer :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol,cgns_bc
character(len = 100) :: arg

real(kind=8),allocatable :: xyz (:,:,:,:)
real(kind=8),allocatable :: vec (:,:,:,:)
real(kind=8) :: mat(2,2)
real(kind=8) :: temp(dimen)


integer :: i,j,k

integer(kind=CGSIZE_T) :: isize (dimen,3)
integer(kind=CGSIZE_T) :: ibc   (dimen,2,6)

write(*,'(A)') "GRID GEN for 3D Sod Shock Tube"
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


allocate(xyz (imax,jmax,kmax,Dimen))
allocate(vec (imax-1,jmax-1,kmax-1,nVar ))
do k = 1,kmax
   do j = 1,jmax
      do i = 1,imax
         xyz(i,j,k,1) = length/dble(imax-1) * dble(i-1)
   
         xyz(i,j,k,2) = length/dble(imax-1) * dble(j-1)
         xyz(i,j,k,3) = length/dble(imax-1) * dble(k-1)
         temp = xyz(i,j,k,:)
         xyz(i,j,k,1) = mat(1,1) * temp(1) + mat(2,1) * temp(2)
         xyz(i,j,k,2) = mat(1,2) * temp(1) + mat(2,2) * temp(2)
   
      end do
   end do
end do
do k = 1, kmax-1
   do j = 1,jmax-1
      do i = 1,imax-1
         vec(i,j,k,1  ) = 1.0D0
         vec(i,j,k,2:3) = 0.0D0
         vec(i,j,k,4  ) = 1.0D0 / 0.4D0
         if (i+j+k > 1.5*jmax) then
            vec(i,j,k,1) = 0.125D0
            vec(i,j,k,4) = 0.1D0 / 0.4D0
         end if
      end do
   end do
end do

call cg_open_f("data_in.cgns",CG_MODE_WRITE,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_base_write_f(cgns_file,"Grid",Dimen,Dimen,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

isize(1,1) = imax
isize(2,1) = jmax
isize(3,1) = kmax

isize(1,2) = isize(1,1) - 1
isize(2,2) = isize(2,1) - 1
isize(3,2) = isize(3,1) - 1

isize(1,3) = 0
isize(2,3) = 0
isize(3,3) = 0

call cg_zone_write_f(cgns_file,cgns_base,"Zone 1",isize,Structured,cgns_zone,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!            write grid                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateX",xyz(:,:,:,1),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateY",xyz(:,:,:,2),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateZ",xyz(:,:,:,3),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!            write solution                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
         ,"Density",vec(:,:,:,1),cgns_var,ierror)
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
         ,"Geschw_U",vec(:,:,:,2),cgns_var,ierror)

call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
         ,"Geschw_V",vec(:,:,:,3),cgns_var,ierror)

call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
         ,"Geschw_W",vec(:,:,:,3),cgns_var,ierror)
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
         ,"Energie",vec(:,:,:,4),cgns_var,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!            write boundary                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! WEST INFLOW
ibc(1,1,1)  = 1
ibc(2,1,1)  = 1
ibc(3,1,1)  = 1
ibc(1,2,1)  = 1
ibc(2,2,1)  = jmax - 1
ibc(3,2,1)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'WEST',BCInflow,PointRange,2,ibc(:,:,1),cgns_bc,ierror)

!!! EAST OUTFLOW
ibc(1,1,2)  = imax - 1
ibc(2,1,2)  = 1
ibc(3,1,2)  = 1
ibc(1,2,2)  = imax - 1
ibc(2,2,2)  = jmax - 1
ibc(3,2,2)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'EAST',BCOutflow,PointRange,2,ibc(:,:,2),cgns_bc,ierror)

!!! SOUTH WALL
ibc(1,1,3)  = 1
ibc(2,1,3)  = 1
ibc(3,1,3)  = 1
ibc(1,2,3)  = imax - 1
ibc(2,2,3)  = 1
ibc(3,2,3)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'SOUTH',BCSymmetryPlane,PointRange,2,ibc(:,:,3),cgns_bc,ierror)

!!! WEST INFLOW
ibc(1,1,4)  = 1
ibc(2,1,4)  = jmax - 1
ibc(3,1,4)  = 1
ibc(1,2,4)  = imax - 1
ibc(2,2,4)  = jmax - 1
ibc(3,2,4)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'NORTH',BCSymmetryPlane,PointRange,2,ibc(:,:,4),cgns_bc,ierror)

!!! WEST INFLOW
ibc(1,1,5)  = 1
ibc(2,1,5)  = 1
ibc(3,1,5)  = 1
ibc(1,2,5)  = imax - 1
ibc(2,2,5)  = jmax - 1
ibc(3,2,5)  = 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'FRONT',BCSymmetryPlane,PointRange,2,ibc(:,:,5),cgns_bc,ierror)

!!! WEST INFLOW
ibc(1,1,6)  = 1
ibc(2,1,6)  = 1
ibc(3,1,6)  = kmax - 1
ibc(1,2,6)  = imax - 1
ibc(2,2,6)  = jmax - 1
ibc(3,2,6)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'BACK',BCSymmetryPlane,PointRange,2,ibc(:,:,6),cgns_bc,ierror)

call cg_close_f(cgns_file,ierror)

if (ierror /= CG_OK) call cg_error_exit_f()

write(*,'(A)') "done"

end program
