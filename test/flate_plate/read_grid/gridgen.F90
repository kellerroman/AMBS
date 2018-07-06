program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none

integer            :: imax = 100
integer            :: jmax = 30
integer, parameter :: kmax = 1

integer, parameter :: nVar   = 5

real(kind=8) :: p,rho,u,v,t,ma

integer :: i,j,k,b
integer :: io
real(kind=8), allocatable :: coords(:,:,:)

open(newunit=io,file="fplate.grd")

read(io,*)
read(io,*)
read(io,*)
read(io,*) imax,jmax
read(io,*)
write(*,*) "Gitterdimensions from File:", imax,jmax
allocate(coords(imax+1,jmax+1,2))
do j = 1, jmax + 1
   do i = 1, imax + 1
      read(io,*) coords(i,j,:)
   end do
end do

close(io)

write(*,*) "Flate PLate Grid generator"
call add_block(imax, jmax,kmax)
call allocate_blocks(nVar)
io = 0
b = 1
do i = 1,blocks(b) % npkts(1)
   do j = 1, blocks(b) % npkts(2)
      do k = 1, blocks(b) % npkts(3)
       blocks(b) % xyzs(i,j,k,1:2) = coords(i+io,j,:)
       blocks(b) % xyzs(i,j,k,3) = dble(k-1) !*  2.0D-3
       end do
    end do
end do
io = io + blocks(b) % nCells(1)
deallocate(coords)
t = 288.15
ma = 0.5D0
p = 1E5
rho = p / (RGAS * t)
u = ma * sqrt( gamma * RGAS * t) 
v = 0.0D0

write(*,'(A20,1X,ES10.3)') "Density:",rho
write(*,'(A20,1X,ES10.3)') "Velocity:",u
write(*,'(A20,1X,ES10.3)') "Mach-Number:",ma
write(*,'(A20,1X,ES10.3)') "Temperature:",t
write(*,'(A20,1X,ES10.3)') "Pressure:",p
!write(*,*) rho,p,t,u

blocks(b) % vars(:,:,:,1) = rho
blocks(b) % vars(:,:,:,2) = u
blocks(b) % vars(:,:,:,3) = v
blocks(b) % vars(:,:,:,4) = v
blocks(b) % vars(:,:,:,5) = p / gm1 + 0.5D0*rho*(u*u+v*v)

blocks(b) % boundary_condition(DIR_WEST)  % bc_type = BC_INFLOW_SUB
blocks(b) % boundary_condition(DIR_WEST)  % velocity = u
blocks(b) % boundary_condition(DIR_EAST)  % bc_type = BC_OUTFLOW
blocks(1) % boundary_condition(DIR_EAST)  % pressure = p
blocks(b) % boundary_condition(DIR_SOUTH) % bc_type = BC_WALL
blocks(b) % boundary_condition(DIR_NORTH) % bc_type = BC_OUTFLOW
blocks(1) % boundary_condition(DIR_NORTH)  % pressure = p
blocks(b) % boundary_condition(DIR_FRONT) % bc_type = BC_SYMMETRY
blocks(b) % boundary_condition(DIR_BACK)  % bc_type = BC_SYMMETRY

call write_grid()
call write_xdmf()
write(*,'(A)') "done"

end program
