program flate_plate_extract
use mod_post
use const_mod, only: REAL_KIND
implicit none
integer, parameter :: Y_SPU = 3
integer :: fo
integer :: b, i, j, k, s, var

integer :: i_start,i_end
integer :: j_start,j_end

character(len=*), parameter :: file_out = "data_1d.csv"
real(REAL_KIND), allocatable :: grenzschicht_dicke(:)
real(REAL_KIND) :: avg_vel,y1,y2,v1,v2
!< Welche Richtung soll variert werden

write(*,*) "========== EXTRACT 1D SOLUTION =============="
debug = .true.
debug = .false.
call read_solution()
write(*,*) "Solution with ",nBlock," Blocks, ",nVar, " Variables  and ",nSol," Solutions read"
do i = 1, nVar
   write(*,*) i, varnames(i)
end do
b = 1
i = 1
k = 1
avg_vel = 0.0E0_REAL_KIND
do j = 1, blocks(b) % nCells(2)
   avg_vel = avg_vel + blocks(b) % solutions(nSol) % vars(i,j,k,Y_SPU)
end do
avg_vel = avg_vel / blocks(b) % nCells(2) 
!write(*,*) avg_vel
allocate( grenzschicht_dicke(blocks(b) % nCells(1)))
do i = 1, blocks(b) % nCells(1)
!   avg_vel = 0.0E0_REAL_KIND
!   do j = 1, blocks(b) % nCells(2)
!      avg_vel = avg_vel + blocks(b) % vars(i,j,k,Y_SPU)
!   end do
!   avg_vel = avg_vel / blocks(b) % nCells(2) 
   do j = 1, blocks(b) % nCells(2)
      if  ( blocks(b) % solutions(nSol) % vars(i,j,k,Y_SPU) >= 0.99E0_REAL_KIND * avg_vel) then
         if (j > 1) then
            y1 = 0.5E0_REAL_KIND * ( blocks(b) % coords(i,j,k,2) + blocks(b) % coords(i,j-1,k,2))
            y2 = 0.5E0_REAL_KIND * ( blocks(b) % coords(i,j,k,2) + blocks(b) % coords(i,j+1,k,2))
            v1 = blocks(b) % solutions(nSol) % vars(i,j-1,k,Y_SPU)
            v2 = blocks(b) % solutions(nSol) % vars(i,j  ,k,Y_SPU)
            grenzschicht_dicke(i) = y1 + (y2-y1) / (v2-v1) * (0.99E0_REAL_KIND * avg_vel-v1)
         else
            grenzschicht_dicke(i) = 0.0E0_REAL_KIND 
            !* ( blocks(b) % coords(i,j,k,Y_SPU) + blocks(b) % coords(k,j+1,k,Y_SPU))
         end if
         exit
      end if
   end do
   !write(*,*) i,j, grenzschicht_dicke(i), blocks(b) % solutions(1) % vars(i,j-1:j+1,k,Y_SPU),avg_vel
end do

open(newunit = fo , file = trim(file_out))
   i_start = 1
   i_end   = blocks(b) % nCells(1)
   j_start = 1 !max(1,blocks(b) % nCells(2) / 2)
   j_end   = 1 !max(1,blocks(b) % nCells(2) / 2)
write(fo,'(A20)',advance="no") "COORD"
do var = 1,nVar
   write(fo,'(",",A20)',advance="no") trim(varnames(var))
end do
write(fo,*)
do i = i_start, i_end
   do j = j_start, j_end
      write(fo ,'(F20.13)',advance = "no") (blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i+1,j,k,1)) * 0.5d0
      do s = 1, nSol
         do var = 1,nVar
            write(fo,'(",",F20.13)',advance = "no") blocks(b) % solutions(s) % vars(i,j,k,var)
         end do
      end do
      write(fo,*)
   end do
end do
write (fo,*) 
close(fo)

open(newunit = fo , file = "grenzschicht_dicke.csv")
j = 1
k = 1
do i = 1, blocks(b) % nCells(1) 
   write(fo,*) (blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i+1,j,k,1)) * 0.5d0, grenzschicht_dicke(i) 
end do
deallocate ( grenzschicht_dicke)
close(fo)

open(newunit = fo , file = "profile.csv")
j = 1
i = 1
k = 1
do j = blocks(b) % nCells(2),1,-1 
write(fo,*) (blocks(b) % coords(i,j,k,2)+blocks(b) % coords(i,j+1,k,2)) * 0.5d0,blocks(b) % solutions(nSol) % vars(:,j,k,Y_SPU) 
end do
close(fo)
write(*,*) "========== EXTRACT 1D SOLUTION done ========="
end program flate_plate_extract
