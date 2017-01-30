program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none

integer, parameter :: imax1 = 10
integer, parameter :: imax = 100
integer, parameter :: jmax = 30
integer, parameter :: kmax = 2

integer, parameter :: nVar   = 5

real(kind=8) :: p,rho,u,v,gm1,t,dy,ma,vorz

integer :: i,j,b

write(*,*) "Flate PLate Grid geneerator"
call add_block(imax1-1,jmax-1,kmax-1)
call add_block(imax-1, jmax-1,kmax-1)
call allocate_blocks(nVar)
dy = 1.0D-4
vorz = -1.0D0 * (imax1-1)/(imax-1)
do b = 1, 2
   do i = 1,blocks(b) % npkts(1)
       dy = 5.0D-4
       blocks(b) % xyzs(i,:,:,1) = vorz + 1.0D0 * dble(i-1) / dble(imax-1)
       do j = 1, blocks(b) % npkts(2)
         if (j == 1) then
            blocks(b) % xyzs(i,j,:,2) =  0.0D0
         else
            blocks(b) % xyzs(i,j,:,2) = blocks(b) % xyzs(i,j-1,:,2) + dy
            dy = dy * 1.1D0
         end if
           blocks(b) % xyzs(i,j,1,3) = 0.0D0
           blocks(b) % xyzs(i,j,2,3) = blocks(b) % xyzs(i,j,1,3) + 2.0D-3
       end do
   end do
   vorz = 0.0D0
end do
gm1 = gamma - 1.0D0
t = 388.889D0
ma = 0.1D0
p = 41368.5D0
rho = p / (RGAS * t)
u = ma * sqrt( gamma * RGAS * t) 
v = 0.0D0
!write(*,*) rho,p,t,u

do b = 1, 2
   blocks(b) % vars(:,:,:,1) = rho
   blocks(b) % vars(:,:,:,2) = u
   blocks(b) % vars(:,:,:,3) = v
   blocks(b) % vars(:,:,:,4) = v
   blocks(b) % vars(:,:,:,5) = p / gm1 + 0.5D0*rho*(u*u+v*v)
end do
b = 1
blocks(b) % boundary_condition(DIR_WEST) = BC_INFLOW
blocks(b) % boundary_condition(DIR_EAST) = 2
blocks(b) % boundary_condition(DIR_SOUTH) = BC_SYMMETRY
b = 2
blocks(b) % vars(:,:,:,2) = 0.0d0
blocks(b) % vars(:,:,:,5) = p / gm1 + 0.5D0*rho*(u*u+v*v)
blocks(b) % boundary_condition(DIR_WEST) = 1
blocks(b) % boundary_condition(DIR_EAST) = BC_OUTFLOW
blocks(b) % boundary_condition(DIR_SOUTH) = BC_WALL
do b = 1, 2
   !blocks(1) % boundary_condition(DIR_NORTH) = BC_SYMMETRY
   blocks(b) % boundary_condition(DIR_NORTH) = BC_OUTFLOW
   blocks(b) % boundary_condition(DIR_FRONT) = BC_SYMMETRY
   blocks(b) % boundary_condition(DIR_BACK) = BC_SYMMETRY
end do

call write_grid()
call write_xdmf()
write(*,'(A)') "done"

end program
