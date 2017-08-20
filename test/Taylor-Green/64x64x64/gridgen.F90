program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none 
integer       , parameter   :: NVAR   = 5
real(kind = 8), parameter   :: PI     = 4 * atan(1.0E0_8)
!!!!! GIT DIMENSIONS
integer       , parameter   :: IMAX   = 65
integer       , parameter   :: JMAX   = 65
integer       , parameter   :: KMAX   = 65
!integer       , parameter   :: IMAX   = 33
!integer       , parameter   :: JMAX   = 33
!integer       , parameter   :: KMAX   = 33
real(kind = 8), parameter   :: LENGTH = 1.D-02
real(kind = 8), parameter   :: DOMAIN_LENGTH = 2.D-00 * PI * LENGTH

integer                     :: i,j,k,b
real(kind = 8)              :: T, p, E
real(kind = 8)              :: x,y,z
real(kind = 8)              :: rho0, p0, V0

write(*,'(A)') "GRID GEN for 64^3 DIT"
call add_block(imax-1,jmax-1,kmax-1)
call allocate_blocks(nVar)
b = 1
do k = 1,kmax
   do j = 1,jmax
      do i = 1,imax
         blocks(1) % xyzs(i,j,k,1) = DOMAIN_LENGTH/dble(IMAX-1) * dble(i-1)
         blocks(1) % xyzs(i,j,k,2) = DOMAIN_LENGTH/dble(JMAX-1) * dble(j-1)
         blocks(1) % xyzs(i,j,k,3) = DOMAIN_LENGTH/dble(KMAX-1) * dble(k-1)
      end do
   end do
end do
V0 = 30.0D0
T = 300.0d0
p0 = 1.0D5
rho0 =  p0 / ( T * RGAS )
do k = 1, kmax-1
   do j = 1,jmax-1
      do i = 1,imax-1
         x = 0.125D0 * ( blocks(b) % xyzs(i  ,j  ,k  ,1) &
                       + blocks(b) % xyzs(i+1,j  ,k  ,1) & 
                       + blocks(b) % xyzs(i  ,j+1,k  ,1) & 
                       + blocks(b) % xyzs(i  ,j  ,k+1,1) & 
                       + blocks(b) % xyzs(i+1,j+1,k  ,1) & 
                       + blocks(b) % xyzs(i+1,j  ,k+1,1) & 
                       + blocks(b) % xyzs(i  ,j+1,k+1,1) & 
                       + blocks(b) % xyzs(i+1,j+1,k+1,1) )
         y = 0.125D0 * ( blocks(b) % xyzs(i  ,j  ,k  ,2) &
                       + blocks(b) % xyzs(i+1,j  ,k  ,2) & 
                       + blocks(b) % xyzs(i  ,j+1,k  ,2) & 
                       + blocks(b) % xyzs(i  ,j  ,k+1,2) & 
                       + blocks(b) % xyzs(i+1,j+1,k  ,2) & 
                       + blocks(b) % xyzs(i+1,j  ,k+1,2) & 
                       + blocks(b) % xyzs(i  ,j+1,k+1,2) & 
                       + blocks(b) % xyzs(i+1,j+1,k+1,2) )
         z = 0.125D0 * ( blocks(b) % xyzs(i  ,j  ,k  ,3) &
                       + blocks(b) % xyzs(i+1,j  ,k  ,3) & 
                       + blocks(b) % xyzs(i  ,j+1,k  ,3) & 
                       + blocks(b) % xyzs(i  ,j  ,k+1,3) & 
                       + blocks(b) % xyzs(i+1,j+1,k  ,3) & 
                       + blocks(b) % xyzs(i+1,j  ,k+1,3) & 
                       + blocks(b) % xyzs(i  ,j+1,k+1,3) & 
                       + blocks(b) % xyzs(i+1,j+1,k+1,3) )
         blocks(1) % vars(i,j,k,2) =  V0 * sin (x / LENGTH) * cos(y / LENGTH) * cos(z / LENGTH) 
         blocks(1) % vars(i,j,k,3) = -V0 * cos (x / LENGTH) * sin(y / LENGTH) * cos(z / LENGTH) 
         blocks(1) % vars(i,j,k,4) = 0.0d0
         p = p0 + rho0*V0*V0 / 16.0D0 * (cos(2.0D0 * x / LENGTH) + cos(2.0D0 * y / LENGTH) ) & 
                                      * (cos(2.0D0 * z / LENGTH) + 2.0D0)
         blocks(1) % vars(i,j,k,1) =  p / ( T * RGAS )
         E = p / GM1 + 0.5D0 * rho0 * &
                     ( blocks(1) % vars(i,j,k,2) * blocks(1) % vars(i,j,k,2) &
                     + blocks(1) % vars(i,j,k,3) * blocks(1) % vars(i,j,k,3) &
                     + blocks(1) % vars(i,j,k,4) * blocks(1) % vars(i,j,k,4) )
         blocks(1) % vars(i,j,k,5) = E
      end do
   end do
end do

blocks(1) % boundary_condition(DIR_WEST)  = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_EAST)  = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_SOUTH) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_NORTH) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_FRONT) = BC_SYMMETRY
blocks(1) % boundary_condition(DIR_BACK)  = BC_SYMMETRY

call write_grid()
call write_xdmf()
write(*,'(A)') "done"

end program
