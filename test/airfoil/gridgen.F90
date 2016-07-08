program gridgen
use cgns

implicit none

real(kind=8),parameter :: RGas = 287.102D0

integer, parameter :: imax = 71
integer, parameter :: jmax = 48
integer, parameter :: kmax = 2

integer, parameter :: Dimen = 3
integer, parameter :: nVar   = 4

real(kind=8) :: xyz (imax,jmax,kmax,Dimen)
real(kind=8),allocatable :: vec (:,:,:,:)
real(kind=8) :: p,gamma,rho,u,v,gm1,T

integer :: i,j

integer :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol,cgns_bc

integer(kind=CGSIZE_T) :: isize (dimen,3)
integer(kind=CGSIZE_T) :: ibc   (dimen,2,6)

write(*,*) "AIRFOIL converter"

open(unit=2,file="fort.1.fine")
do i = 1,imax
    do j = 1, jmax
        read(2,1000) xyz(i,j,1,1:2)
        xyz(i,j,1,3) = 0.0D0
        xyz(i,j,2,:) = xyz(i,j,1,:)
        xyz(i,j,2,3) = xyz(i,j,1,3) + 2.0D-2
    end do
    read(2,2000)
end do
close(2)
allocate(vec (imax-1,jmax-1,kmax-1,nVar ))
Gamma  = 1.4D0
gm1 = gamma - 1.0D0
rho = 0.01d0
u = 694.44d0
v = 0.0D0
T = 300.0D0
p = rho * T * RGas
vec(:,:,:,1) = rho
vec(:,:,:,2) = u 
vec(:,:,:,3) = v
vec(:,:,:,4) = p / gm1 + 0.5D0*rho*(u*u+v*v)

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
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
   ,"Geschw_U",vec(:,:,:,2),cgns_var,ierror) 
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
   ,"Geschw_V",vec(:,:,:,3),cgns_var,ierror) 
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
   ,"Geschw_W",vec(:,:,:,3),cgns_var,ierror) 
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
   ,"Energie",vec(:,:,:,4),cgns_var,ierror) 
if (ierror /= CG_OK) call cg_error_exit_f() 

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
if (ierror /= CG_OK) call cg_error_exit_f()

!!! EAST OUTFLOW
ibc(1,1,2)  = imax - 1
ibc(2,1,2)  = 1
ibc(3,1,2)  = 1
ibc(1,2,2)  = imax - 1
ibc(2,2,2)  = jmax - 1
ibc(3,2,2)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'EAST',BCOutflow,PointRange,2,ibc(:,:,2),cgns_bc,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!!! SOUTH WALL
ibc(1,1,3)  = 1
ibc(2,1,3)  = 1
ibc(3,1,3)  = 1
ibc(1,2,3)  = imax - 1
ibc(2,2,3)  = 1
ibc(3,2,3)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'SOUTH',BCWall,PointRange,2,ibc(:,:,3),cgns_bc,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!!! WEST INFLOW
ibc(1,1,4)  = 1
ibc(2,1,4)  = jmax - 1
ibc(3,1,4)  = 1
ibc(1,2,4)  = imax - 1
ibc(2,2,4)  = jmax - 1
ibc(3,2,4)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'NORTH',BCOutflow,PointRange,2,ibc(:,:,4),cgns_bc,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!!! WEST INFLOW
ibc(1,1,5)  = 1
ibc(2,1,5)  = 1
ibc(3,1,5)  = 1
ibc(1,2,5)  = imax - 1
ibc(2,2,5)  = jmax - 1
ibc(3,2,5)  = 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'FRONT',BCSymmetryPlane,PointRange,2,ibc(:,:,5),cgns_bc,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!!! WEST INFLOW
ibc(1,1,6)  = 1
ibc(2,1,6)  = 1
ibc(3,1,6)  = kmax - 1
ibc(1,2,6)  = imax - 1
ibc(2,2,6)  = jmax - 1
ibc(3,2,6)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'BACK',BCSymmetryPlane,PointRange,2,ibc(:,:,6),cgns_bc,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_close_f(cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

write(*,*) "done"

1000 format(1x,e11.4,5x,e11.4)
2000 format (1x)
end program
