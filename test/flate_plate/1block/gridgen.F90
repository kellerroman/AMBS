program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none

integer, parameter :: imax = 100
integer, parameter :: jmax = 30
integer, parameter :: kmax = 2

integer, parameter :: nVar   = 5

real(kind=8) :: p,rho,u,v,t,dy,ma

integer :: i,j

write(*,*) "Flate PLate Grid geneerator"
call add_block(imax-1,jmax-1,kmax-1)
call allocate_blocks(nVar)
dy = 1.0D-4
do i = 1,imax
    dy = 5.0D-4
    do j = 1, jmax
      blocks(1) % xyzs(i,j,:,1) = 1.0D0 * dble(i-21) / dble(imax-1)
      if (j == 1) then
         blocks(1) % xyzs(i,j,:,2) =  0.0D0
      else
         blocks(1) % xyzs(i,j,:,2) = blocks(1) % xyzs(i,j-1,:,2) + dy
         dy = dy * 1.1D0
      end if
        blocks(1) % xyzs(i,j,1,3) = 0.0D0
        blocks(1) % xyzs(i,j,2,3) = blocks(1) % xyzs(i,j,1,3) + 2.0D-3
    end do
end do
t = 388.889D0
ma = 0.1D0
p = 41368.5D0
rho = p / (RGAS * t)
u = ma * sqrt( gamma * RGAS * t) 
v = 0.0D0
write(*,*) rho,p,t,u

blocks(1) % vars(:,:,:,1) = rho
blocks(1) % vars(:,:,:,2) = u
blocks(1) % vars(:,:,:,3) = v
blocks(1) % vars(:,:,:,4) = v
blocks(1) % vars(:,:,:,5) = p / gm1 + 0.5D0*rho*(u*u+v*v)

blocks(1) % boundary_condition(DIR_WEST)  % bc_type = BC_INFLOW_SUB
blocks(1) % boundary_condition(DIR_WEST)  % velocity = u
blocks(1) % boundary_condition(DIR_EAST)  % bc_type = BC_OUTFLOW
blocks(1) % boundary_condition(DIR_EAST)  % pressure = p
blocks(1) % boundary_condition(DIR_SOUTH) % bc_type = BC_WALL
!blocks(1) % boundary_condition(DIR_NORTH) % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_NORTH) % bc_type = BC_OUTFLOW
blocks(1) % boundary_condition(DIR_NORTH) % pressure = p
blocks(1) % boundary_condition(DIR_FRONT) % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_BACK)  % bc_type = BC_SYMMETRY

call write_grid()
call write_xdmf()
write(*,'(A)') "done"

end program
