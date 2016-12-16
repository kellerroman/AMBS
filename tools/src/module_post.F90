module mod_post
use hdf5
use file_io_mod, only: SOL_NAME_LENGTH
use screen_io_mod, only: error_wr
use const_mod, only: REAL_KIND, INT_KIND
implicit none
private
integer         , parameter :: VARNAME_LENGTH        = 20
character(len=*), parameter :: GROUP_GRID            = "grid"
character(len=*), parameter :: GROUP_DATA            = "data"
character(len=*), parameter :: GROUP_BLOCK           = "block"

character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

integer, parameter     ::   RANK = 3                        ! dataset rank

integer :: nBlock
integer :: nVar
integer :: nSol
integer :: Dimen
logical :: debug = .false.

character(len=VARNAME_LENGTH),  allocatable :: varnames(:)
character(len=SOL_NAME_LENGTH), allocatable :: solnames(:)

type :: tsol
   real(REAL_KIND), allocatable :: vars(:,:,:,:)
end type

type :: tblock
   integer :: ncells(3)
   integer :: npkts(3)
   real(REAL_KIND), allocatable :: coords(:,:,:,:)
   type(tsol), allocatable :: solutions(:)
end type

type(tblock), allocatable :: blocks(:)

interface read_solution
   module procedure read_solution_gf ! given Filename
   module procedure read_solution_sf ! stadnard filename
end interface read_solution

public:: blocks, debug, nblock, nVar, nSol, read_solution, varnames, solnames

contains
   subroutine read_solution_sf()
      implicit none
      character(len=*), parameter :: FILENAME              = "data_out.h5"  ! file name HDF5
      call read_solution_gf(FILENAME)
   end subroutine


   subroutine read_solution_gf(filename)
   implicit none

   character(len=*), intent(in) :: filename

   integer(HSIZE_T), dimension(RANK) :: dims
   integer(HSIZE_T), dimension(RANK) :: dims2
   integer(HSIZE_T), dimension(RANK) :: maxdims
   integer     ::   error ! error flag

   integer(hid_t) :: file_id       ! file identifier
   integer(hid_t) :: dset_id       ! dataset identifier
   integer(hid_t) :: group_id      ! dataset identifier
   integer(hid_t) :: group_id1     ! dataset identifier
   integer(hid_t) :: group_id2     ! dataset identifier
   integer(hid_t) :: dspace_id     ! dataspace identifier
      
   character(len=7) :: block_group

   integer :: b, s, d, var
   integer :: solution_type, var_type ! dataspace identifier

   logical :: fexists

   inquire(file=trim(filename),exist=fexists)
   
   if(.not. fexists) then
     call error_wr("Data Input File konnte nicht gefunden werden: "//TRIM(filename),__FILE__,__LINE__)
   end if
   
   if (debug) &
      write(*,*) "write_grid: Writing HDF5 output File"

   ! Initialize FORTRAN interface.
   call h5open_f(error)

   ! Open an existing file.
   call h5fopen_f (filename, h5f_acc_rdwr_f, file_id, error)

   call h5gopen_f(file_id,GROUP_GRID,group_id,error)

   call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)

   allocate(blocks(nBlock))
   do b = 1, nBlock

      write(block_group,'(A,I0)') GROUP_BLOCK, b
      call h5gopen_f(group_id,block_group,group_id2,error)
      
      ! Open an existing dataset.
      call h5dopen_f(group_id2, COORD_NAME(1), dset_id, error)
      call h5dget_space_f(dset_id,dspace_id,error)
      
      call h5sget_simple_extent_ndims_f(dspace_id,dimen,error)
      call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

      call h5dclose_f(dset_id, error)
      blocks(b) % nPkts = INT(dims,INT_KIND)

      blocks(b) % nCells = max(1,blocks(b) % nPkts - 1)

      allocate (blocks(b) % coords(blocks(b)%nPkts(1) &
                                ,blocks(b)%nPkts(2) &
                                ,blocks(b)%nPkts(3) &
                                ,Dimen))

      do d = 1,dimen
         call h5dopen_f(group_id2, COORD_NAME(d), dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, blocks(b) % coords (:,:,:,d), dims, error)
         call h5dclose_f(dset_id, error)
      end do
      call h5gclose_f(group_id2, error)
   end do
   call h5gclose_f(group_id, error) ! CLOSE GRID GROUP

   call h5gopen_f(file_id,GROUP_DATA,group_id,error)

   call h5gn_members_f(file_id, GROUP_DATA, nSol, error)
   if (debug) write(*,*) "Number of Solutions", nSol
   
   allocate(solnames(nSol))
   do b = 1, nBlock
      allocate(blocks(b) % solutions(nSol))
   end do
   do s = 1, nSol

      call h5gget_obj_info_idx_f(file_id, GROUP_DATA, s - 1,solnames(s), solution_type, error)
      if (debug) write(*,*) "Opening group ", solnames(s)
      call h5gopen_f(group_id,solnames(s),group_id1,error) ! OPEN TIMESTEP GROUP
      if (error /= 0) then
         call error_wr("HDF5 Error",__FILE__,__LINE__)
      end if
      do b = 1, nBlock
         write(block_group,'(A,I0)') GROUP_BLOCK, b
         !write(*,*) "opening", block_group
         if (debug) write(*,*) "Opening group ", block_group
         call h5gopen_f(group_id1,block_group,group_id2,error)
         if ( s == 1 .and. b == 1) then
            call h5gn_members_f(group_id1, block_group, nVar, error)
            allocate( varnames(nVar))
            if (debug) &
               write(*,*) "Anzahl der Variablen:", nVar
         end if
         allocate(blocks(b) % solutions(s) % vars(blocks(b) % nCells(1) &
                                                 ,blocks(b) % nCells(2) &
                                                 ,blocks(b) % nCells(3) &
                                                 ,nVar))
         do var = 1, nVar
            call h5gget_obj_info_idx_f(group_id1, block_group, var-1,varnames(var), var_type, error)

            call h5dopen_f(group_id2, varnames(var), dset_id, error)

            call h5dget_space_f(dset_id,dspace_id,error)
            
            call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, blocks(b)  % solutions(s) % vars(:,:,:,var), dims, error)
            call h5dclose_f(dset_id, error)
         end do
         call h5gclose_f(group_id2, error)
      end do
      
      call h5gclose_f(group_id1, error) !CLOSE TIMESTEP GROUP
   end do
   call h5gclose_f(group_id, error)  !CLOSE DATA GROUP
   ! close the file.
   call h5fclose_f(file_id, error)
   
   ! close fortran interface.
   call h5close_f(error)


   end subroutine

end module mod_post
