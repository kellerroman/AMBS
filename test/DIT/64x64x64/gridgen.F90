program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none 
integer       , parameter   :: NVAR   = 5
real(kind = 8), parameter   :: PI     = 3.183E-00 
!!!!! GIT DIMENSIONS
integer       , parameter   :: IMAX   = 65
integer       , parameter   :: JMAX   = 65
integer       , parameter   :: KMAX   = 65
real(kind = 8), parameter   :: LENGTH = 1.D-02 * PI

integer                     :: fu
integer                     :: i,j,k
real(kind = 8)              :: rho, T, p, E

write(*,'(A)') "GRID GEN for 64^3 DIT"
call add_block(imax-1,jmax-1,kmax-1)
call allocate_blocks(nVar)
do k = 1,kmax
   do j = 1,jmax
      do i = 1,imax
      blocks(1) % xyzs(i,j,k,1) = length/dble(imax-1) * dble(i-1)
      blocks(1) % xyzs(i,j,k,2) = length/dble(imax-1) * dble(j-1)
      blocks(1) % xyzs(i,j,k,3) = length/dble(imax-1) * dble(k-1)
   
      end do
   end do
end do
open(newunit=fu,file="inco64.in")
p = 1.0D+5
T = 300.0d0
rho  = p / ( T * RGAS )
do k = 1, kmax-1
   do j = 1,jmax-1
      do i = 1,imax-1
         blocks(1) % vars(i,j,k,1  ) = rho
         read(fu,*) blocks(1) % vars(i,j,k,2:4)
         blocks(1) % vars(i,j,k,2:4) =  blocks(1) % vars(i,j,k,2:4) * 60.0E0_REAL_KIND
         E = p / GM1 + 0.5D0 * rho * &
                     ( blocks(1) % vars(i,j,k,2) * blocks(1) % vars(i,j,k,2) &
                     + blocks(1) % vars(i,j,k,3) * blocks(1) % vars(i,j,k,3) &
                     + blocks(1) % vars(i,j,k,4) * blocks(1) % vars(i,j,k,4) )
         blocks(1) % vars(i,j,k,5  ) = E
      end do
   end do
end do
close(fu)

blocks(1) % boundary_condition(DIR_WEST)  % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_EAST)  % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_SOUTH) % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_NORTH) % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_FRONT) % bc_type = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_BACK)  % bc_type = BC_SYMMETRY

call write_grid()
call write_xdmf()
write(*,'(A)') "done"

end program
