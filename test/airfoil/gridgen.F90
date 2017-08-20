program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC
implicit none

integer, parameter :: imax = 71
integer, parameter :: jmax = 48
integer, parameter :: kmax = 3

integer, parameter :: nVar   = 5

real(kind=8) :: p,rho,u,v,T

integer :: i,j

write(*,*) "AIRFOIL converter"
debug = .true.

call add_block(imax-1,jmax-1,kmax-1)
call add_block(imax-1,jmax-1,kmax-1)

call allocate_blocks(nVar)

open(unit=2,file="fort.1.fine")
do i = 1,imax
    do j = 1, jmax
       read(2,1000) blocks(1) % xyzs(i,j,1,1:2)
       blocks(1) % xyzs(i,j,1,3) = 0.0D0
       blocks(1) % xyzs(i,j,2,:) = blocks(1) % xyzs(i,j,1,:)
       blocks(1) % xyzs(i,j,2,3) = blocks(1) % xyzs(i,j,1,3) + 2.0D-2
       blocks(1) % xyzs(i,j,3,:) = blocks(1) % xyzs(i,j,1,:)
       blocks(1) % xyzs(i,j,3,3) = blocks(1) % xyzs(i,j,2,3) + 2.0D-2
       blocks(2) % xyzs(i,j,1,1:2) =  blocks(1) % xyzs(i,j,1,1:2)
       blocks(2) % xyzs(i,j,1,3) = 6.0D-2
       blocks(2) % xyzs(i,j,2,:) = blocks(2) % xyzs(i,j,1,:)
       blocks(2) % xyzs(i,j,2,3) = blocks(2) % xyzs(i,j,1,3) + 2.0D-2
       blocks(2) % xyzs(i,j,3,:) = blocks(2) % xyzs(i,j,1,:)
       blocks(2) % xyzs(i,j,3,3) = blocks(2) % xyzs(i,j,2,3) + 2.0D-2
    end do
    read(2,2000)
end do
close(2)
rho = 0.01d0
u = 694.44d0
v = 0.0D0
T = 300.0D0
p = rho * T * RGas
blocks(1) % vars(:,:,:,1) = rho
blocks(1) % vars(:,:,:,2) = u 
blocks(1) % vars(:,:,:,3) = v
blocks(1) % vars(:,:,:,4) = v
blocks(1) % vars(:,:,:,5) = p / gm1 + 0.5D0*rho*(u*u+v*v)

blocks(1) % boundary_condition(DIR_WEST) = BC_INFLOW
blocks(1) % boundary_condition(DIR_EAST) = BC_OUTFLOW
blocks(1) % boundary_condition(DIR_SOUTH) = BC_WALL
blocks(1) % boundary_condition(DIR_NORTH) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_FRONT) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_BACK) = BC_SYMMETRY

call write_grid()
call write_xdmf()
write(*,*) "done"

1000 format(1x,e11.4,5x,e11.4)
2000 format (1x)
end program
