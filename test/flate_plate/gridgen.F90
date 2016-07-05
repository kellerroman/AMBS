program gridgen
use cgns

implicit none

integer, parameter :: imax = 100
integer, parameter :: jmax = 30
integer, parameter :: kmax = 2

integer, parameter :: Version = 1000
INTEGER, PARAMETER :: ioout = 10

integer, parameter :: Dimen = 3
integer, parameter :: nBlock = 1
integer, parameter :: nVar   = 4

real(kind=8),parameter :: RGas = 287.102D0

integer, parameter :: bc(4) = (/-2,-3,-1,-4/)

real(kind=8) :: xyz (imax,jmax,kmax,Dimen)
real(kind=8),allocatable :: vec (:,:,:,:)
real(kind=8) :: p,gamma,rho,u,v,gm1,t,dy,ma

integer :: i,j

integer(kind = CGSIZE_T) :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol,cgns_bc

integer(kind=CGSIZE_T) :: isize (dimen,3)
integer(kind=CGSIZE_T) :: ibc   (dimen,2,6)

write(*,*) "AIRFOIL converter"
dy = 1.0D-4
do i = 1,imax
    dy = 5.0D-5
    do j = 1, jmax
      xyz(i,j,:,1) = 1.0D0 * dble(i-21) / dble(imax-1)
      if (j == 1) then
         xyz(i,j,:,2) =  0.0D0
      else
         xyz(i,j,:,2) = xyz(i,j-1,:,2) + dy
         dy = dy * 1.1D0
      end if
        xyz(i,j,1,3) = 0.0D0
        xyz(i,j,2,3) = xyz(i,j,1,3) + 2.0D-3
    end do
end do
allocate(vec (imax-1,jmax-1,kmax-1,nVar+2 ))
gamma  = 1.4D0
gm1 = gamma - 1.0D0
t = 388.889D0
ma = 0.1D0
p = 41368.5D0
rho = p / (RGAS * t)
u = ma * sqrt( gamma * RGAS * t) 
v = 0.0D0
write(*,*) rho,p,t,u

vec(:,:,:,1) = rho
vec(:,:,:,2) = u
vec(:,:,:,3) = v
vec(:,:,:,4) = p / gm1 + 0.5D0*rho*(u*u+v*v)
vec(:,:,:,5) = p
vec(:,:,:,6) = t
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
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'SOUTH',BCWall,PointRange,2,ibc(:,:,3),cgns_bc,ierror)

!!! WEST INFLOW
ibc(1,1,4)  = 1
ibc(2,1,4)  = jmax - 1
ibc(3,1,4)  = 1
ibc(1,2,4)  = imax - 1
ibc(2,2,4)  = jmax - 1
ibc(3,2,4)  = kmax - 1
call cg_boco_write_f(cgns_file,cgns_base,cgns_zone,'NORTH',BCOutflow,PointRange,2,ibc(:,:,4),cgns_bc,ierror)

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



open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")

write (ioout) Version,Dimen,nBlock

i = 3*4
write(ioout) i

write (ioout)  bc
close (ioout)
write(*,*) "done"

end program
