program extract_1D
use mod_post
use file_io_mod, only: VARNAME_SPU,VARNAME_SPV,VARNAME_SPW
 use, intrinsic :: iso_c_binding 
implicit none
include 'fftw3.f03'
integer, parameter :: fo = 12345
character(len=*), parameter :: file_data_in = "data_out.h5"
character(len=*), parameter :: file_out_name = "velocities_"
character(len=*), parameter :: file_out_ending = ".dat"

integer :: b, i ,j ,k ,s
integer :: index_spu,index_spv,index_spw
integer(INT_KIND):: nCell
type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), allocatable :: out(:,:,:,:)

real(REAL_KIND) :: energy

write(*,*) "========== EXTRACT VELOCITIES FROM DIT SOLUTION =============="

call read_solution(trim(file_data_in))

if (nBlock > 1) then
   write(*,*) "Extract 1D nur möglich mit einem Block:",nBlock
   stop 1
end if
Write(*,'(A4,3(1X,A4))') "#Bl","NI","NJ","NK"
do b = 1,nBlock
   write(*,'(I4,3(1X,I4))') b, blocks(b) % nCells
   nCell = nCell + product(blocks(b) % nCells)
end do
write(*,*) "Overall Cell Number:", nCell
   write(*,*) "Number of Variables on File",nVar
b = 1
s = 1
write(*,*) varnames
index_spu = -1
index_spv = -1
index_spw = -1

do s = 1, nVar
   if (varnames(s) == VARNAME_SPU) then
      index_spu = s
   else if (varnames(s) == VARNAME_SPV) then
      index_spv = s
   else if (varnames(s) == VARNAME_SPW) then
      index_spw = s
   end if
end do
if (index_spu == -1 ) then
   write(*,*) VARNAME_SPU//" not found"
   stop 1
end if
if (index_spv == -1 ) then
   write(*,*) VARNAME_SPV//" not found"
   stop 1
end if
if (index_spw == -1 ) then
   write(*,*) VARNAME_SPW//" not found"
   stop 1
end if
allocate(out(blocks(b) % nCells(1),blocks(b) % nCells(2),blocks(b) % nCells(3),3))
plan = fftw_plan_dft_r2c_3d(blocks(b) % nCells(1) &
                           ,blocks(b) % nCells(2) &
                           ,blocks(b) % nCells(3) &
                           ,blocks(b) % solutions(1) % vars(:,:,:,index_spu) &
                           ,out(:,:,:,1) &
                           ,FFTW_ESTIMATE)
open(unit = fo , file = file_out_name//trim(solnames(s))//file_out_ending)
do s = 1, nSol
   
   call fftw_execute_dft_r2c(plan, blocks(b) % solutions(s) % vars(:,:,:,index_spu), out(:,:,:,1))
   call fftw_execute_dft_r2c(plan, blocks(b) % solutions(s) % vars(:,:,:,index_spv), out(:,:,:,2))
   call fftw_execute_dft_r2c(plan, blocks(b) % solutions(s) % vars(:,:,:,index_spw), out(:,:,:,3))
   do k = 1, blocks(b) % nCells(3)
      do j = 1, blocks(b) % nCells(2)
         do i = 1, blocks(b) % nCells(1)
!            write(fo,*) blocks(b) % solutions(s) % vars(i,j,k,index_spu) &
!                      , blocks(b) % solutions(s) % vars(i,j,k,index_spu) &
!                      , blocks(b) % solutions(s) % vars(i,j,k,index_spu)
            energy = 0.5E0_REAL_KIND * (sum(out(i,j,k,:)*out(i,j,k,:)))
         end do
      end do
   end do
end do
call fftw_destroy_plan(plan)
close(fo)
write(*,*) "========== EXTRACT VELOCITIES FROM DIT SOLUTION ========="
end program extract_1D


