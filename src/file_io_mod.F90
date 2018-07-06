module file_io_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: module for file output
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 04.03.2016
! 
! LAST CHANGE: 04.03.2016
!
! CHANGELOG:
! 04.03.2016,RK: Start of Coding
!
! SUBROUTINES:
!     datin_control : Read control file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use const_mod
   use hdf5
   use screen_io_mod
implicit none

   character(len=*), parameter :: VARNAME_RHO           = "Density"  
   character(len=*), parameter :: VARNAME_SPU           = "Geschw_U" 
   character(len=*), parameter :: VARNAME_SPV           = "Geschw_V" 
   character(len=*), parameter :: VARNAME_SPW           = "Geschw_W" 
   character(len=*), parameter :: VARNAME_ENE           = "Energie"  
   character(len=*), parameter :: VARNAME_PRE           = "Druck"
   character(len=*), parameter :: VARNAME_TEMP          = "Temperatur"
   character(len=*), parameter :: VARNAME_MACH          = "Mach_Number"

   character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

   character(len=*)  , parameter   :: FOLDER_INPUT       = "./"
   character(len=*)  , parameter   :: FOLDER_OUTPUT      = "./"
   character(len=100)              :: filename_control   = FOLDER_INPUT //"config.bin"
   character(len=100)              :: filename_data_in   = FOLDER_INPUT //"data_in.h5"
   character(len=100)              :: filename_bc_in     = FOLDER_INPUT //"bc.bin"
   character(len=100)              :: filename_data_out  = FOLDER_OUTPUT//"data_out.h5"
   
   character(len=*), parameter :: GROUP_GRID            = "grid"
   character(len=*), parameter :: GROUP_DATA            = "data"
   character(len=*), parameter :: GROUP_BLOCK           = "block"

   integer                           , parameter :: CONFIG_FILE_VERSION = 2 

   integer         , parameter :: SOL_NAME_LENGTH = 11
   logical :: override_sol = .false.
   logical :: write_sol_header = .true.
contains
   subroutine datin_control()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: routine reads in binary control file which was parsed with preprocessor
!          python script from ASCII format to BINARY format
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 04.03.2016
! 
! LAST CHANGE: 04.03.2016
!
! CHANGELOG:
! 04.03.2016,RK: Start of Coding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use control_mod, only: max_iteration   &
                        , solution_out    &
                        , residual_out    &
                        , equation        &
                        , turbulence      &
                        , space_order     &
                        , riemann_solver  &
                        , timestep_method &
                        , time_order      &
                        , cfl             &
                        , timestep        &
                        , c_les_sgs
   implicit none 
      character(len=*), parameter :: SOLVERS(2) = ["Roe          ","Lax-Friedrich"]
      character(len=*), parameter :: STR_EQUATION(2) = ["Euler        " ,"Navier-Stokes"]
      integer                         :: file_unit
      integer                         :: cfv_version
      logical :: fexists
      inquire(file=trim(filename_control),exist=fexists)
      if (.not. fexists) then
         call error_wr("Config Datei: '"//trim(filename_control)//"' nicht gefunden!",__FILE__,__LINE__)
      end if
      open(newunit = file_unit, file=trim(filename_control),form="unformatted",access="stream")
      read(file_unit) cfv_version
      if (cfv_version /= CONFIG_FILE_VERSION) then
         close(file_unit,status='delete')
         call error_wr("Config Datei Version stimmt nicht",__FILE__,__LINE__)
      end if
      read(file_unit) max_iteration   
      read(file_unit) solution_out    
      read(file_unit) residual_out    
      read(file_unit) equation        
      read(file_unit) turbulence      
      read(file_unit) space_order     
      read(file_unit) riemann_solver  
      read(file_unit) timestep_method 
      read(file_unit) time_order      
      read(file_unit) cfl             
      read(file_unit) timestep
      read(file_unit) c_les_sgs
      close(file_unit,status='delete')
      
      write(*,'(A30,2X,I0)')     "Number of Iterations",max_iteration
      write(*,'(A30,2X,I0)')     "Solution Output every",solution_out
      write(*,'(A30,2X,I0)')     "Residual Output every",residual_out
      write(*,'(A30,2X,I0)')     "Spatial Order",space_order
      write(*,'(A30,2X,A)')      "Equation",STR_EQUATION(equation)
      write(*,'(A30,2X,ES10.4)') "Timestep",Timestep    
      write(*,'(A30,2X,A)')      "Riemann Solver",trim(SOLVERS(riemann_solver))
      write(*,'(A30,2X,ES10.4)') "LES SGS Konstant", c_les_sgs

   end subroutine datin_control
   
   subroutine datin_sol()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: routine reads in Solution 
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 15.03.2016
! 
! LAST CHANGE: 15.03.2016
!
! CHANGELOG:
! 15.03.2016,RK: Start of Coding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use data_mod
   use hdf5
   implicit none
      !integer :: ib, d,ibb
      logical :: fexists
      !integer                :: normalindex(3)
      !integer                :: normallist,ndataset,normaldatatype,iptset,ibocotype
      !integer                :: dir_of_boundary
      
      real(REAL_KIND), allocatable :: data_in(:,:,:)
      integer     ::   error ! Error flag
      integer(hid_t) :: file_id       ! file identifier
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: group_id      ! dataset identifier
      integer(hid_t) :: group_id1     ! dataset identifier
      integer(hid_t) :: group_id2     ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer        :: solution_type ! dataspace identifier
      integer        :: var_type ! dataspace identifier
      integer(HSIZE_T) :: dims(3)
      integer(HSIZE_T) :: maxdims(3)
      integer     ::   ib, d, var
      integer     ::   hdf5_nSol, nVar_in
      character(len=len(GROUP_BLOCK)+2) :: block_group
      character(len=10) :: solution_name
      character(len=30) :: varName_in
      nCell = 0
      
      inquire(file=trim(filename_data_in),exist=fexists)
      
      if(.not. fexists) then
        call error_wr("Data Input File konnte nicht gefunden werden: " &
                      //TRIM(filename_data_in),__FILE__,__LINE__)
      end if
      
      
      ! Initialize FORTRAN interface.
      call h5open_f(error)

      ! Open an existing file.
      call h5fopen_f (filename_data_in, h5f_acc_rdwr_f, file_id, error)

      call h5gopen_f(file_id,GROUP_GRID,group_id,error)

      call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
!      if (nBlock > 1) then
!         write(*,*) "Only one Block supported at the Moment"
!         write(*,*) "Number of Blocks in File:",nBlock
!         nBlock = 1
!      end if
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
      allocate(blocks(nBlock))
      do ib = 1, nBlock


         write(block_group,'(A,I0)') GROUP_BLOCK, ib
         call h5gopen_f(group_id,block_group,group_id2,error)
         
         ! Open an existing dataset.
         call h5dopen_f(group_id2, COORD_NAME(1), dset_id, error)
         call h5dget_space_f(dset_id,dspace_id,error)
         
         call h5sget_simple_extent_ndims_f(dspace_id,dimen,error)
         call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

         call h5dclose_f(dset_id, error)
         blocks(ib) % nPkts = INT(dims,INT_KIND)

         blocks(ib) % nCells = max(1,blocks(ib) % nPkts - 1)
         nCell = nCell + product(blocks(ib) % nCells)

!         write(*,*) dimen,dims
         
         nFaces = dimen * 2
         nCorners = 2 ** dimen
         !!!! MUSS VERSCHOBEN WERDEN FUER MULTIBLOCK
         call allocate_vars(ib)

         allocate (data_in(blocks(ib)%nPkts(1),blocks(ib)%nPkts(2),blocks(ib)%nPkts(3)))
         do d = 1,dimen
            call h5dopen_f(group_id2, COORD_NAME(d), dset_id, error)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_in, dims, error)
            blocks(ib) % coords (1:blocks(ib)%nPkts(1) &
                               ,1:blocks(ib)%nPkts(2) &
                               ,1:blocks(ib)%nPkts(3),d) = data_in
            call h5dclose_f(dset_id, error)
         end do
         deallocate(data_in)
         call h5gclose_f(group_id2, error)
         ! terminate access to the data space.
         !call h5sclose_f(dspace_id, error)
      end do
      call h5gclose_f(group_id, error) ! CLOSE GRID GROUP
      call h5gopen_f(file_id,GROUP_DATA,group_id,error)

      call h5gn_members_f(file_id, GROUP_DATA, hdf5_nSol, error)
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
      if (hdf5_nSol > 1) then
         call error_wr("More than One Solution in Restart-File",__FILE__,__LINE__)
      end if
      call h5gget_obj_info_idx_f(file_id, GROUP_DATA, 0,solution_name, solution_type, error)
      !write(*,*) solution_name
      call h5gopen_f(group_id,solution_name,group_id1,error) ! OPEN TIMESTEP GROUP
      do ib = 1, nBlock
         write(block_group,'(A,I0)') GROUP_BLOCK, ib
         !write(*,*) "opening", block_group
         call h5gopen_f(group_id1,block_group,group_id2,error)
         call h5gn_members_f(group_id1, block_group, nVar_in, error)
         if (ib == 1) &
         write(*,'(A,I0,A)',ADVANCE="NO") "Variablen(", nVar_in,"):"
         allocate (data_in(blocks(ib)%nCells(1),blocks(ib)%nCells(2),blocks(ib)%nCells(3)))
         do var = 1, nVar_in
            call h5gget_obj_info_idx_f(group_id1, block_group, var-1,varName_in, var_type, error)
         if (ib == 1) &
            write(*,'(1X,A)',ADVANCE="NO") trim(varName_in)
            call h5dopen_f(group_id2, varName_in, dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            
            call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_in, dims, error)
            call h5dclose_f(dset_id, error)
            select case(VarName_in)
               case(VarName_Rho)
                  blocks(ib) % vars(1:blocks(ib) % nCells(1)            &
                                   ,1:blocks(ib) % nCells(2)            &
                                   ,1:blocks(ib) % nCells(3),1) = data_in
               case(VarName_SpU)
                  blocks(ib) % vars(1:blocks(ib) % nCells(1)            &
                                   ,1:blocks(ib) % nCells(2)            &
                                   ,1:blocks(ib) % nCells(3),2) = data_in
               case(VarName_SpV)
                  blocks(ib) % vars(1:blocks(ib) % nCells(1)            &
                                   ,1:blocks(ib) % nCells(2)            &
                                   ,1:blocks(ib) % nCells(3),3) = data_in
               case(VarName_SpW)
                  blocks(ib) % vars(1:blocks(ib) % nCells(1)            &
                                   ,1:blocks(ib) % nCells(2)            &
                                   ,1:blocks(ib) % nCells(3),4) = data_in
               case(VarName_Ene)
                  blocks(ib) % vars(1:blocks(ib) % nCells(1)            &
                                   ,1:blocks(ib) % nCells(2)            &
                                   ,1:blocks(ib) % nCells(3),5) = data_in
               case default
                  write(*,*) varname_in," nicht erkannt"
             end select
         end do
         if (ib == 1) write(*,*)
         deallocate(data_in)
         call h5gclose_f(group_id2, error)
      end do
      
      call h5gclose_f(group_id1, error) !CLOSE TIMESTEP GROUP
      call h5gclose_f(group_id, error)  !CLOSE DATA GROUP
   ! close the file.
   call h5fclose_f(file_id, error)
   
   ! close fortran interface.
   call h5close_f(error)
   end subroutine datin_sol

   subroutine datin_bc()
      use data_mod, only: nBlock, blocks
   implicit none
      integer                         :: file_unit
      logical                         :: fexists
      integer                         :: ib,bc


      inquire(file=trim(filename_bc_in),exist=fexists)
      if (.not. fexists) then
         call error_wr("Config Datei: '"//trim(filename_bc_in)//"' nicht gefunden!",__FILE__,__LINE__)
      end if
      open(newunit = file_unit, file=trim(filename_bc_in),form="unformatted",access="stream")
      write(*,'(88("="))')
      write(*,'(A4,3(1X,A5),6(1X,A10))') "B#","NI","NJ","NK","WEST","EAST","SOUTH","NORTH","FRONT","BACK"
      write(*,'(88("-"))')
      do ib = 1, nBlock
      write(*,'(I4,3(1X,I5))',ADVANCE="NO") ib, blocks(ib) % NCells
         !write(*,'("========== BLOCK",I2.2,"==========")') ib
         do bc = 1,6
            read(file_unit) blocks(ib) % boundary(bc) % bc_type
            if (blocks(ib) % boundary(bc) % bc_type == BC_OUTFLOW) then
               read(file_unit) blocks(ib) % boundary(bc) % pressure
            else if (blocks(ib) % boundary(bc) % bc_type == BC_INFLOW_SUB) then
               read(file_unit) blocks(ib) % boundary(bc) % pressure
            end if
            write(*,'(1X,A10)',ADVANCE="NO") bc_names(blocks(ib) % boundary(bc) % bc_type)
         end do
         write(*,*)
      end do
      write(*,'(88("="))')
   end subroutine datin_bc
   subroutine datout_sol()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: routine writes out Solution 
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 14.04.2016
! 
! LAST CHANGE: 14.04.2016
!
! CHANGELOG:
! 14.04.2016,RK: Start of Coding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use data_mod
   use control_mod, only: current_iteration
   implicit none
   integer :: b, d, var
   character(len=SOL_NAME_LENGTH)  :: solname
   
   integer(hsize_t) :: dims(3)
!   real(kind=dp) :: T
   real(REAL_KIND), allocatable :: data_out(:,:,:)

   character(len=7) :: block_group
   integer(hid_t) :: file_id       ! file identifier
   integer(hid_t) :: dset_id       ! dataset identifier
   integer(hid_t) :: group_id      ! dataset identifier
   integer(hid_t) :: group_id1     ! dataset identifier
   integer(hid_t) :: group_id2     ! dataset identifier
   integer(hid_t) :: dspace_id     ! dataspace identifier
   integer     ::   error ! error flag
   character(len=20) :: varname_out

   CALL h5open_f(error)
   if (write_sol_header) then
      !   call wr("Writing New Solution to File "//trim(filename_data_out),4)
      if (.not. override_sol) write_sol_header = .false.

      call h5fcreate_f(filename_data_out, h5f_acc_trunc_f, file_id, error)

      call h5gcreate_f(file_id,  GROUP_GRID, group_id,  error)

      do b = 1, nBlock
         dims(:) = blocks(b) % nPkts  (1:Dimen)
         !isize(:) = blocks(b) % nCells (1:Dimen)
         call h5screate_simple_f(Dimen, dims, dspace_id, error)
         write(block_group,'(A,I0)') GROUP_BLOCK,b
         call h5gcreate_f(group_id, block_group, group_id2, error)
         allocate (data_out ( blocks(b) % nPkts(1), blocks(b) % nPkts(2), blocks(b) % nPkts(3) ))
         do d = 1, Dimen
            data_out = blocks(b) % coords(1:blocks(b) % nPkts(1)            &
                                        ,1:blocks(b) % nPkts(2)            &
                                        ,1:blocks(b) % nPkts(3),d)
            call h5dcreate_f(group_id2, COORD_NAME(d), h5t_native_double, dspace_id, &
                            dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dims, error)
            call h5dclose_f(dset_id, error)
         end do
         ! Close the group.
         call h5gclose_f(group_id2, error)
         deallocate(data_out)
      end do
      ! Close the group.
      call h5gclose_f(group_id, error)
      call h5gcreate_f(file_id,  GROUP_DATA, group_id,  error)
   else
      !write(*,*) "Open excisting file"
      CALL h5fopen_f (trim(filename_data_out), H5F_ACC_RDWR_F, file_id, error)
      call h5gopen_f(file_id,GROUP_DATA,group_id,error)
   end if

   write(solname,'(I10.10)') current_iteration
   call screen_wr("Writing Solution@"//trim(solname)//" to File "//trim(filename_data_out),1)
   ! Create a group data in the file.
   call h5gcreate_f(group_id,  trim(solname), group_id1,  error)
   do b = 1,nBlock
      write(block_group,'(A,I0)') GROUP_BLOCK,b
      dims = blocks(b) % ncells 
      call h5screate_simple_f(Dimen, dims, dspace_id, error)

      ! Create a group named for block1 in the file.
      call h5gcreate_f(group_id1, block_group, group_id2, error)

      allocate (data_out(blocks(b) % nCells(1),blocks(b) % nCells(2),blocks(b) % nCells(3))) 

      do var = 1, 12!sol_out_nVar
         select case(var)!VarName(var))
         case(1)!VarName_Rho)
            varname_out = VarName_Rho
            data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_RHO)
         case(2)!VarName_SpU)
            varname_out = VarName_SpU
            data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_SPU) 
         case(3)!VarName_SpV)
            varname_out = VarName_SpV
            data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_SPV)
         case(4)!VarName_SpW)
            varname_out = VarName_SpW
            data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_SPW)
         case(5)!VarName_Ene)
            varname_out = VarName_Ene
            data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_ENE)
         case(7)!VarName_Pre)
            varname_out = "Druck"
            data_out = blocks(b) % pressures(1:blocks(b) % nCells(1)  &
                                            ,1:blocks(b) % nCells(2)  &
                                            ,1:blocks(b) % nCells(3))
         case(8)!VarName_Temp)
            varname_out = "Temperatur"
            data_out = blocks(b) % temperatures(1:blocks(b) % nCells(1)  &
                                               ,1:blocks(b) % nCells(2)  &
                                               ,1:blocks(b) % nCells(3))
         case(6)
            varname_out = "Mach Number"
            data_out = sqrt( &
                       blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_SPU) ** 2 &
                     + blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_SPV) ** 2 &
                     + blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3),VEC_SPW) ** 2)&
                     / blocks(b) % schallges(1:blocks(b) % nCells(1)  &
                                       ,1:blocks(b) % nCells(2)  &
                                       ,1:blocks(b) % nCells(3))
         case(9)!VarName_Temp)
            varname_out = "timestep"
            data_out = blocks(b) % timesteps(1:blocks(b) % nCells(1)  &
                                            ,1:blocks(b) % nCells(2)  &
                                            ,1:blocks(b) % nCells(3))
!         case(9)!VarName_Temp)
!            varname_out = "dUdyJ"
!            data_out = blocks(b) % dudnJ(1:blocks(b) % nCells(1)  &
!                                        ,1:blocks(b) % nCells(2)  &
!                                        ,1:blocks(b) % nCells(3),1,2)
         case(10)!VarName_Temp)
            varname_out = "visFluxesJ"
            data_out = blocks(b) % visFluxesJ(1:blocks(b) % nCells(1)  &
                                        ,1:blocks(b) % nCells(2)  &
                                        ,1:blocks(b) % nCells(3),VEC_SPU)
         case(11)!VarName_Temp)
            varname_out = "FluxesJ"
            data_out = blocks(b) % FluxesJ(1:blocks(b) % nCells(1)  &
                                        ,1:blocks(b) % nCells(2)  &
                                        ,1:blocks(b) % nCells(3),VEC_SPU)
         case(12)!VarName_Temp)
            varname_out = "Residual"
            data_out = blocks(b) % residuals(1:blocks(b) % nCells(1)  &
                                            ,1:blocks(b) % nCells(2)  &
                                            ,1:blocks(b) % nCells(3),VEC_SPU)
         end select
         
         call h5dcreate_f(group_id2, varname_out, h5t_native_double, dspace_id, &
              dset_id, error)
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dims, error)
         call h5dclose_f(dset_id, error)
      end do
      deallocate(data_out)
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
   end subroutine datout_sol
end module file_io_mod
