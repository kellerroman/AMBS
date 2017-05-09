program sol2vtk
use mod_post
implicit none
integer :: fo
integer :: b, i ,j ,k ,d, n

character(len = 100) :: arg
character(len=100) :: filename_sol

integer :: s,var
integer :: var_dir = 1
CHARACTER(LEN=5) :: BlockNr

!< Welche Richtung soll variert werden
filename_sol = "data_out.h5"
write(*,*) "========== SOLUTION TO VTK =============="

call read_solution(trim(filename_sol))

write(*,*) "Solution with ",nBlock," Blocks, ",nVar, " Variables  and ",nSol," Solutions read"

!do b = 1,nBlock
b = 1
do s = 1,nSol
   write(*,'(I0,1X)',ADVANCE="No") s
   !write(*,666) s,nSol
   write(blocknr,'(i5.5)') s
   open(NEWUNIT=fo,file="sol"//BlockNr//".vtk")
   write(fo,'(a)') '# vtk DataFile Version 3.1'
   write(fo,'(a)') 'Structured Grid file from AMBS'
   write(fo,'(a)') 'ASCII'
   !write(fo,*)
   write(fo,'(a)') 'DATASET STRUCTURED_GRID'
   write(fo,'(a11,3I4)') 'DIMENSIONS ', blocks(b) % nPkts(:)
   write(fo,*) 'POINTS ',product(blocks(b) % nPkts(:)),' double'
   do k = 1,blocks(b) % nPkts(3)
      do j = 1,blocks(b) % nPkts(2)
         do i = 1,blocks(b) % nPkts(1)
            write(fo,*) (blocks(b) % coords(i,j,k,n),n=1,3)
         end do
      end do
   end do
   write(fo,*)
   write(fo,"(A,1X,I0)") 'CELL_DATA',product(blocks(b) % nCells(:))
   do d = 1, nVar
      write(fo,"(A)") 'SCALARS '//TRIM(varnames(d))//' double'
      write(fo,"(A)") 'LOOKUP_TABLE default'
      do k = 1,blocks(b) % nCells(3)
          do j = 1,blocks(b) % nCells(2)
              do i = 1,blocks(b) % nCells(1)
                  write(fo,*) blocks(b) % solutions(s) % vars(i,j,k,d)
              end do
          end do
      end do
   end do
   close(fo)
end do

write(*,*)
write(*,*) "========== SOLUTION TO VTK done ========="
!666 FORMAT('+',T2,'CASE ',I3,' OF ',I3,' COMPUTED')
end program sol2vtk
