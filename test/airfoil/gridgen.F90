program gridgen
use cgns

implicit none

real(kind=8),parameter :: RGas = 287.102D0

integer, parameter :: imax = 71
integer, parameter :: jmax = 48
integer, parameter :: kmax = 2
integer, parameter :: Version = 1000
INTEGER, PARAMETER :: ioout = 10

integer, parameter :: Dimen = 3
integer, parameter :: nBlock = 1
integer, parameter :: nVar   = 4

integer, parameter :: bc(4) = (/-2,-3,-1,-4/)

real(kind=8) :: xyz (imax,jmax,kmax,Dimen)
real(kind=8),allocatable :: vec (:,:,:,:)
real(kind=8) :: p,gamma,rho,u,v,gm1,T

integer :: i,j

integer :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol

integer(kind=CGSIZE_T) :: isize (dimen,3)

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
call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateX",xyz(:,:,:,1),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateY",xyz(:,:,:,2),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateZ",xyz(:,:,:,3),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

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


call cg_close_f(cgns_file,ierror)

if (ierror /= CG_OK) call cg_error_exit_f()



open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")

write (ioout) Version,Dimen,nBlock

i = 3*4
write(ioout) i

write (ioout)  bc
close (ioout)
write(*,*) "done"

1000 format(1x,e11.4,5x,e11.4)
2000 format (1x)
end program
