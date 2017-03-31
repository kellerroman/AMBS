program gridgen

   use hdf5 ! this module contains all necessary modules
   use mod_gridgen
   use const_mod!, only: BC_OUTFLOW, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC

implicit none

integer            :: imax1 = 10
integer            :: imax = 100
integer            :: jmax = 30
integer, parameter :: kmax = 1

integer, parameter :: nVar   = 5

real(kind=8) :: p,rho,u,v,t,ma

integer :: i,j,k,b
integer :: io
integer :: ni_git
real(kind=8), allocatable :: coords(:,:,:)
logical :: blockgrenze_found

open(newunit=io,file="fplate.grd")

read(io,*)
read(io,*)
read(io,*)
read(io,*) ni_git,jmax
read(io,*)
write(*,*) "Gitterdimensions from File:", ni_git,jmax
allocate(coords(ni_git+1,jmax+1,2))
blockgrenze_found = .false.
do j = 1, jmax + 1
   do i = 1, ni_git + 1
      read(io,*) coords(i,j,:)
      if (.not. blockgrenze_found) then
         if (abs(coords(i,j,1)) <= 1E-8) then
            blockgrenze_found = .true.
            imax1 = i - 1
            imax = ni_git - imax1
            write(*,*) "Blockgrenze bei",i
            write(*,*) "Blockdimensionen:",imax1,imax,jmax
            write(*,*) coords(i,1,1)
         end if
      end if
   end do
end do

close(io)

write(*,*) "Flate PLate Grid generator"
call add_block(imax1,jmax,kmax)
call add_block(imax, jmax,kmax)
call allocate_blocks(nVar)
io = 0
do b = 1, 2
   do i = 1,blocks(b) % npkts(1)
      do j = 1, blocks(b) % npkts(2)
         do k = 1, blocks(b) % npkts(3)
          blocks(b) % xyzs(i,j,k,1:2) = coords(i+io,j,:)
          blocks(b) % xyzs(i,j,k,3) = dble(k-1) !*  2.0D-3
          end do
       end do
   end do
   io = io + blocks(b) % nCells(1)
end do
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
!blocks(b) % vars(:,:,:,2) = 0.0d0
!blocks(b) % vars(:,:,:,5) = p / gm1 + 0.5D0*rho*(u*u+v*v)
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
