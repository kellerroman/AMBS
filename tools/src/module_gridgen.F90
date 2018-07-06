module mod_gridgen
use hdf5
use const_mod
implicit none
private
integer         , parameter :: VARNAME_LENGTH        = 20
character(len=*), parameter :: FILENAME              = "data_in.h5"  ! file name HDF5
character(len=*), parameter :: FILENAME_XDMF         = "grid.xdmf"   ! filename XDMF
character(len=*), parameter :: FILENAME_BC           = "bc.bin"  ! file name HDF5
character(len=*), parameter :: GROUP_GRID            = "grid"
character(len=*), parameter :: GROUP_DATA            = "data"
character(len=*), parameter :: GROUP_BLOCK           = "block"

character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

character(len=VARNAME_LENGTH), parameter :: VARNAME_RHO   = "Density"     ! dataset name
character(len=VARNAME_LENGTH), parameter :: VARNAME_SPU   = "Geschw_U"    ! dataset name
character(len=VARNAME_LENGTH), parameter :: VARNAME_SPV   = "Geschw_V"    ! dataset name
character(len=VARNAME_LENGTH), parameter :: VARNAME_SPW   = "Geschw_W"    ! dataset name
character(len=VARNAME_LENGTH), parameter :: VARNAME_ENE   = "Energie"     ! dataset name
integer, parameter     ::   RANK = 3                        ! dataset rank
integer :: nblock = 0
integer :: nVar = 4
integer :: ncells_temp(10,3)
logical :: blocks_allocated = .false.
logical :: debug = .false.

type :: tboundary_condition
   integer           :: bc_type
   real(kind = 8)    :: pressure
   real(kind = 8)    :: velocity
end type tboundary_condition

type :: tblock
   integer                     :: ncells(3)
   integer                     :: npkts(3)
   real(kind = 8), allocatable :: xyzs(:,:,:,:)
   real(kind = 8), allocatable :: vars(:,:,:,:)
   type(tboundary_condition)   :: boundary_condition(6)
end type

type(tblock), allocatable :: blocks(:)
interface add_block
   module procedure add_block_scalar
   module procedure add_block_array
end interface add_block

public:: blocks,  add_block, allocate_blocks, debug, write_grid, write_xdmf

contains
   subroutine add_block_array(ncells)
      implicit none
      integer, intent(in) :: ncells(3)
      call add_block(ncells(1),ncells(2),ncells(3))
   end subroutine

   subroutine add_block_scalar(ni,nj,nk)
      implicit none
      integer, intent(in) :: ni,nj,nk
      if (blocks_allocated) then
         write(*,*) "add_block: Cannot add Block"
         write(*,*) "Blocks are allready allocated"
         stop 1
      end if
      nblock = nblock + 1
      ncells_temp(nblock,1)  = ni
      ncells_temp(nblock,2)  = nj
      ncells_temp(nblock,3)  = nk
      if (debug) &
         write(*,*) "Adding Block: NB:",nblock,"(ni,nj,nk)",ni,nj,nk
   end subroutine

   subroutine allocate_blocks(anzahl_var)
      implicit none
      integer, intent(in) :: anzahl_var
      integer :: i
      if (blocks_allocated) then
         write(*,*) "allocate_blocks: Canno allocate Blocks"
         write(*,*) "Blocks are allready allocated"
         stop 1
      end if
      if (nblock > 0) then
         blocks_allocated = .true.
         if (anzahl_var < 5) then
            write(*,*) "allocate_blocks: Cannot allocate Blocks"
            write(*,*) "Passed nVar is to small",anzahl_var
         else
            nvar = anzahl_var
         end if

         if (debug) &
            write(*,*) "Allocate Blocks",nblock,"nVar:",nVar

         allocate( blocks(nblock))
         do i = 1,nblock
            blocks(i) % ncells = ncells_temp(i,:)
            blocks(i) % npkts = ncells_temp(i,:) + 1
            if (debug) then
               write(*,*) i, blocks(i)%ncells, blocks(i) % npkts
            end if
            allocate(blocks(i) % xyzs(blocks(i) %  npkts(1),blocks(i) %  npkts(2),blocks(i) %  npkts(3),3))
            allocate(blocks(i) % vars(blocks(i) % ncells(1),blocks(i) % ncells(2),blocks(i) % ncells(3),nvar)) 
            blocks(i) % boundary_condition(:) % pressure = -1.0d0
            blocks(i) % boundary_condition(:) % velocity = -1.0d0
         end do
      else
         write(*,*) "allocate_blocks: Cannot allocate Blocks"
         write(*,*) "no Blocks defined"
         stop 1
      end if
   end subroutine
   subroutine write_grid()
      implicit none

   integer(hsize_t), dimension(RANK) :: DIMS
   integer(hsize_t), dimension(RANK) :: DIMS2
   integer     ::   error ! error flag

   integer(hid_t) :: file_id       ! file identifier
   integer(hid_t) :: dset_id       ! dataset identifier
   integer(hid_t) :: group_id      ! dataset identifier
   integer(hid_t) :: group_id1     ! dataset identifier
   integer(hid_t) :: group_id2     ! dataset identifier
   integer(hid_t) :: dspace_id     ! dataspace identifier
      
   character(len=7) :: block_group

   character(len=VARNAME_LENGTH), allocatable :: varname_out(:)

   integer :: nb 
   integer :: nd
   integer                         :: file_unit

   allocate(varname_out(nvar))
   
   varname_out = [VARNAME_RHO,VARNAME_SPU,VARNAME_SPV,VARNAME_SPW,VARNAME_ENE]
   if (debug) &
      write(*,*) "write_grid: Writing HDF5 output File"

   ! initialize fortran interface.
   call h5open_f(error)
   
   ! create a new file using default properties.
   call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)

   ! Create a group grid in the file.
   call h5gcreate_f(file_id,  GROUP_GRID, group_id,  error)

   !   Write Coordinates into grid/blockN/
   do nb = 1, nBlock
      write(block_group,'(A,I0)') GROUP_BLOCK,nb
      dims  = blocks(nb) % npkts 

      ! create the dataspace.
      call h5screate_simple_f(rank, dims, dspace_id, error)
   
      ! Create a group named for block1 in the file.
      call h5gcreate_f(group_id, block_group, group_id2, error)
      
      do nd = 1, RANK
         call h5dcreate_f(group_id2, COORD_NAME(nd), h5t_native_double, dspace_id, &
                         dset_id, error)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, blocks(nb) % xyzs(:,:,:,nd), dims, error)
         call h5dclose_f(dset_id, error)
      end do

      ! terminate access to the data space of the grid and create a new one for the data
      call h5sclose_f(dspace_id, error)
      ! Close the group.
      call h5gclose_f(group_id2, error)
   end do

   ! Close the group.
   call h5gclose_f(group_id, error)


   ! Create a group data in the file.
   call h5gcreate_f(file_id,  GROUP_DATA, group_id,  error)
   call h5gcreate_f(group_id,  "0", group_id1,  error)

   !   Write Coordinates into data/blockN/
   do nb = 1, nBlock
      write(block_group,'(A,I0)') GROUP_BLOCK,nb
      dims2 = blocks(nb) % ncells 
      call h5screate_simple_f(rank, dims2, dspace_id, error)

      ! Create a group named for block1 in the file.
      call h5gcreate_f(group_id1, block_group, group_id2, error)
       
      do nd = 1, nVar
      
         call h5dcreate_f(group_id2, varname_out(nd), h5t_native_double, dspace_id, &
              dset_id, error)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, blocks(nb) % vars(:,:,:,nd), dims2, error)
         call h5dclose_f(dset_id, error)
      end do

      ! terminate access to the data space.
      call h5sclose_f(dspace_id, error)
      ! Close the group.
      call h5gclose_f(group_id2, error)
   end do

   call h5gclose_f(group_id1, error)
   ! Close the group.
   call h5gclose_f(group_id, error)

   
   ! close the file.
   call h5fclose_f(file_id, error)
   
   ! close fortran interface.
   call h5close_f(error)
   open(newunit = file_unit, file=trim(FILENAME_BC),form="unformatted",access="stream")
   do nb = 1, nBlock
      do nd = 1,6
         write(file_unit) blocks(nb) % boundary_condition(nd) % bc_type

         if  (blocks(nb) % boundary_condition(nd) % bc_type == BC_OUTFLOW) then
            if (blocks(nb) % boundary_condition(nd) % pressure < 0.0D0) then
               write(*,*) "Pressure for Outplow Boundary Condition not set"
               write(*,*) "BLOCK",nb,"DIRECTION: ",DIR_NAMES(nd)
               stop 1
            end if
            write(file_unit) blocks(nb) % boundary_condition(nd) % pressure

         else if  (blocks(nb) % boundary_condition(nd) % bc_type == BC_INFLOW_SUB) then
            if (blocks(nb) % boundary_condition(nd) % velocity < 0.0D0) then
               write(*,*) "Velocity for Inflow SUBSONIC Boundary Condition not set"
               write(*,*) "BLOCK",nb,"DIRECTION: ",DIR_NAMES(nd)
               stop 1
            end if
            write(file_unit) blocks(nb) % boundary_condition(nd) % velocity
         end if
      end do
   end do
   close(file_unit)

   end subroutine

   subroutine write_xdmf()
   implicit none
   integer :: fu

   character(len=VARNAME_LENGTH), allocatable :: varname_out(:)

   character(len=7) :: block_group

   integer :: nb 
   integer :: nd
   integer :: d

   allocate(varname_out(nvar))
   
   varname_out = [VARNAME_RHO,VARNAME_SPU,VARNAME_SPV,VARNAME_SPW,VARNAME_ENE]

   if (debug) &
      write(*,*) "write_xdmf: Writing XDMF output File"
      
   open(newunit=fu,file=FILENAME_XDMF)

   write(fu,'(A)') '<?xml version="1.0" ?>'
   write(fu,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
   write(fu,'(A)') '<Xdmf Version="2.0">'
   write(fu,'(A)') ' <Domain>'
   do nb = 1, nBlock

      write(block_group,'(A,I0)') GROUP_BLOCK,nb

      write(fu,'(A,I0,A)') '   <Grid Name="Block',nb,'" GridType="Uniform">'
      write(fu,'(A,3(I0,1X),A)') '     <Topology TopologyType="3DSMesh" NumberOfElements="'&
                                ,(blocks(nb) % npkts(d),d=RANK,1,-1) &
                                ,'"/>'
      write(fu,'(A)') '     <Geometry GeometryType="X_Y_Z">'
      do nd = 1, RANK
         write(fu,'(A,3(I0,1X),A)') '       <DataItem Dimensions="'&
                                   ,(blocks(nb) % npkts(d),d=RANK,1,-1) &
                                   ,'" NumberType="Float" Precision="7" Format="HDF">'
         write(fu,'(8X,A,":",3("/",A))') FILENAME,GROUP_GRID,trim(block_group),COORD_NAME(nd)
         write(fu,'(A)') '       </DataItem>'
      end do
      write(fu,'(A)') '     </Geometry>'
      do nd = 1, nVar
         write(fu,'(3A)') '     <Attribute Name="',trim(varname_out(nd)),'" Center="Cell">'
         write(fu,'(A,3(I0,1X),A)') '       <DataItem Dimensions="'&
                                   ,(blocks(nb) % npkts(d),d=1,RANK) &
                                   ,'" NumberType="Float" Precision="7" Format="HDF">'
         write(fu,'(8X,A,":",4("/",A))') FILENAME,GROUP_DATA,"0"&
                                        ,trim(block_group),trim(varname_out(nd))
         write(fu,'(A)') '       </DataItem>'
         write(fu,'(A)') '     </Attribute>'
      end do
      write(fu,'(A)') '   </Grid>'
   end do
   write(fu,'(A)') ' </Domain>'
   write(fu,'(A)') '</Xdmf>'

   close(fu)
   end subroutine
end module
