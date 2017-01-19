program extract_1D
!use const_mod
!use data_mod
!use hdf5
!use file_io_mod
use mod_post
implicit none
integer :: fo
integer :: b, i ,j ,k ,d

integer :: i_start,i_end
integer :: j_start,j_end

logical :: fexists

character(len = 100) :: arg
character(len=100) :: filename_sol
character(len=*), parameter :: file_out = "data_1d.csv"

!real(REAL_KIND), allocatable :: temp_coord(:,:,:)
!   
!   integer     ::   error ! Error flag
!   integer(hid_t) :: file_id       ! file identifier
!   integer(hid_t) :: dset_id       ! dataset identifier
!   integer(hid_t) :: group_id      ! dataset identifier
!   integer(hid_t) :: group_id1     ! dataset identifier
!   integer(hid_t) :: group_id2     ! dataset identifier
!   integer(hid_t) :: dspace_id     ! dataspace identifier
!   integer        :: solution_type ! dataspace identifier
!   integer        :: var_type ! dataspace identifier
!   integer(HSIZE_T),allocatable :: dims(:,:)
!   integer(HSIZE_T),allocatable :: dims_data(:,:)
!   integer(HSIZE_T) :: maxdims(3)
!   integer     :: var
!   integer     ::   hdf5_nSol, nVar_in
!   character(len=len(GROUP_BLOCK)+2) :: block_group
!   character(len=20),allocatable :: solution_name(:)
!   character(len=30),allocatable :: varName_in(:)
!   integer :: fu

   !integer :: nd,sol
!   integer:: rank
integer :: s,var
integer :: var_dir = 1
!< Welche Richtung soll variert werden
filename_sol = "data_out.h5"
write(*,*) "========== EXTRACT 1D SOLUTION =============="
i = 1
DO
   CALL get_command_argument(i, arg)
   IF (LEN_TRIM(arg) == 0) EXIT
   write(*,*) arg
   if (trim(arg) == "x") then
       var_dir = 1
   else if (trim(arg) == "y") then
       var_dir = 2
   else if (trim(arg) == "xy") then
       var_dir = 3
   else if (i == 1) then
      filename_sol = trim(arg)
   end if
   i = i+1
END DO

call read_solution(trim(filename_sol))
write(*,*) "Solution with ",nBlock," Blocks, ",nVar, " Variables  and ",nSol," Solutions read"
b = 1
open(newunit = fo , file = trim(file_out))
if (var_dir == 1) then
   i_start = 1
   i_end   = blocks(b) % nCells(1)
   j_start = max(1,blocks(b) % nCells(2) / 2)
   j_end   = max(1,blocks(b) % nCells(2) / 2)
   k = 1
else if (var_dir == 2) then
   i_start = max(1,blocks(b) % nCells(1) / 2)
   i_end   = max(1,blocks(b) % nCells(1) / 2)
   j_start = 1
   j_end   = blocks(b) % nCells(2)
   k = 1
else if (var_dir == 3) then
   i_start = 1
   i_end   = blocks(b) % nCells(1)
   j_start = max(1,blocks(b) % nCells(2) / 2)
   j_end   = max(1,blocks(b) % nCells(2) / 2)
   k = 1
end if 
write(fo,'(A20)',advance="no") "COORD"
do var = 1,nVar
   write(fo,'(",",A20)',advance="no") trim(varnames(var))
end do
write(fo,*)
do i = i_start, i_end
   if (var_dir == 3) then
      j_start = i
      j_end   = i
   end if
   do j = j_start, j_end
      if (var_dir == 1) then
         write(fo ,'(F20.13)',advance = "no") (blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i+1,j,k,1)) * 0.5d0
      else if (var_dir == 2) then
         write(fo ,'(F20.13)',advance = "no") (blocks(b) % coords(i,j,k,2)+blocks(b) % coords(i,j+1,k,2)) * 0.5d0
      else if (var_dir == 3) then
         write(fo ,'(F20.13)',advance = "no") sqrt((blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i,j+1,k,1)) * &
                                                   (blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i,j+1,k,1)) * 0.125D0&
                                                 + (blocks(b) % coords(i,j,k,2)+blocks(b) % coords(i,j+1,k,2)) * &
                                                   (blocks(b) % coords(i,j,k,2)+blocks(b) % coords(i,j+1,k,2)) * 0.125D0)
      end if

      do s = 1, nSol
         do var = 1,nVar!,nVar
            write(fo,'(",",F20.13)',advance = "no") blocks(b) % solutions(s) % vars(i,j,k,var)
         end do
      end do
      write(fo,*)
   end do
end do
write (fo,*) 
close(fo)
write(*,*) "========== EXTRACT 1D SOLUTION done ========="
end program extract_1D


