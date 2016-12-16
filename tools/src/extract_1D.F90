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
character(len=*), parameter :: filename_sol = "data_out.h5"
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
   end if
   i = i+1
END DO

call read_solution(trim(filename_sol))
write(*,*) "Solution with ",nBlock," Blocks, ",nVar, " Variables  and ",nSol," Solutions read"
!inquire(file=trim(filename_sol),exist=fexists)
!
!if(.not. fexists) then
!  write(*,*) "Data Input File konnte nicht gefunden werden: "//TRIM(filename_sol),__FILE__,__LINE__
!end if
!
!! Initialize FORTRAN interface.
!call h5open_f(error)
!
!! Open an existing file.
!call h5fopen_f (filename_sol, h5f_acc_rdwr_f, file_id, error)
!
!call h5gopen_f(file_id,GROUP_GRID,group_id,error)
!
!call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)
!
!allocate(dims(3,nBlock))
!allocate(dims_data(3,nBlock))
!
!
!if (nBlock > 1) then
!   write(*,*) "Extract 1D nur möglich mit einem Block:",nBlock
!   stop 1
!end if
!allocate(blocks(nBlock))
!write(*,*) "Number of Blocks:", nBlock
!Write(*,'(A4,3(1X,A4))') "#Bl","NI","NJ","NK"
!do b = 1, nBlock
!
!   write(block_group,'(A,I0)') GROUP_BLOCK,b
!
!   call h5gopen_f(group_id,block_group,group_id2,error)
!   
!   ! Open an existing dataset.
!   call h5dopen_f(group_id2, COORD_NAME(1), dset_id, error)
!   call h5dget_space_f(dset_id,dspace_id,error)
!   
!   call h5sget_simple_extent_ndims_f(dspace_id,Dimen,error)
!   call h5sget_simple_extent_dims_f(dspace_id,dims(:,b),maxdims,error)
!   call h5dclose_f(dset_id, error)
!
!   nFaces = Dimen * 2
!   nCorners = 2**Dimen
!   blocks(b) % nPkts = 1
!   blocks(b) % nCells = 1
!   blocks(b) % nPkts(1:Dimen) = int(dims(1:Dimen,b))
!   write(*,'(I4,3(1X,I4))') b, blocks(b) % nCells
!   allocate (blocks(b) % coords(1 : blocks(b)%nPkts(1) &
!                               ,1 : blocks(b)%nPkts(2) &
!                               ,1 : blocks(b)%nPkts(3) &
!                               ,Dimen))
!   allocate (temp_coord(blocks(b)%nPkts(1),blocks(b)%nPkts(2),blocks(b)%nPkts(3)))
!   do d = 1,Dimen
!      call h5dopen_f(group_id2, COORD_NAME(d), dset_id, error)
!      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, temp_coord, dims(:,b), error)
!      blocks(b) % coords(1:blocks(b)%nPkts(1) &
!                        ,1:blocks(b)%nPkts(2) &
!                        ,1:blocks(b)%nPkts(3),d) = temp_coord
!      call h5dclose_f(dset_id, error)
!   end do
!   call h5gclose_f(group_id2, error)
!   deallocate (temp_coord)
!end do
!call h5gclose_f(group_id, error)
!
!call h5gopen_f(file_id,GROUP_DATA,group_id,error)
!call h5gn_members_f(file_id, GROUP_DATA, hdf5_nSol, error)
!write(*,*) "Anzahl der Bloecke:" ,nBlock
!write(*,*) "Anzahl der Loesungen:" ,hdf5_nSol
!allocate(solution_name(hdf5_nSol))
!do sol = 1, hdf5_nSol
!   call h5gget_obj_info_idx_f(file_id, GROUP_DATA, sol-1,solution_name(sol), solution_type, error)
!!   write(*,*) solution_name(sol)
!   call h5gopen_f(group_id,solution_name(sol),group_id1,error) ! OPEN TIMESTEP GROUP
!   do b = 1,nBlock
!      write(block_group,'(A,I0)') GROUP_BLOCK,b
!      call h5gopen_f(group_id1,block_group,group_id2,error)   ! open block Group
!      call h5gn_members_f(group_id1, block_group, nVar_in, error)
!
!      if (sol == 1 .and. b == 1) &
!         write(*,*) "Number of Variables on File",nVar_in
!      if (sol == 1 .and. b == 1) &
!         allocate (varname_in(nVar_in))
!      do var = 1, nVar_in
!         call h5gget_obj_info_idx_f(group_id1, block_group, var-1,varName_in(var), var_type, error)
!         call h5dopen_f(group_id2, varName_in(var), dset_id, error)
!         if (sol == 1 .and. var == 1) then
!            call h5dget_space_f(dset_id,dspace_id,error)
!            
!            call h5sget_simple_extent_dims_f(dspace_id,dims_data(:,b),maxdims,error)
!            blocks(b) % nCells(1:Dimen) = int(dims_data(1:Dimen,b))
!            nCell = nCell + product(blocks(b) % nCells)
!            allocate ( blocks(b) % vars      (blocks(b)%nCells(1) &
!                                             ,blocks(b)%nCells(2) &
!                                             ,blocks(b)%nCells(3) &
!                                             ,nVar_in*hdf5_nSol))
!         end if
!         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, blocks(b)%vars(:,:,:,var+(sol-1)*nVar_in), dims_data(:,b), error)
!
!         call h5dclose_f(dset_id, error)
!      end do
!      call h5gclose_f(group_id2, error)
!   end do
!   call h5gclose_f(group_id1, error) !CLOSE TIMESTEP GROUP
!end do
!call h5gclose_f(group_id, error)  !CLOSE DATA GROUP
!! close the file.
!call h5fclose_f(file_id, error)
!write(*,*) "Overall Cell Number:", nCell
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
         do var = 1,nVar,nVar
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


