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
   use cgns
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
   character(len=100)              :: filename_data_out  = FOLDER_OUTPUT//"data_out.h5"
   
   character(len=*), parameter :: GROUP_GRID            = "grid"
   character(len=*), parameter :: GROUP_DATA            = "data"
   character(len=*), parameter :: GROUP_BLOCK           = "block"

   integer                           , parameter :: CONFIG_FILE_VERSION = 2 
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
      
      write(*,'(A30,2X,I0)') "Number of Iterations",max_iteration
      write(*,'(A30,2X,I0)') "Solution Output every",solution_out
      write(*,'(A30,2X,I0)') "Residual Output every",residual_out
      write(*,'(A30,2X,I0)') "Spatial Order",space_order
      write(*,'(A30,2X,ES10.4)') "Timestep",Timestep    
      write(*,'(A30,2X,A)') "Riemann Solver",trim(SOLVERS(riemann_solver))
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
      integer(hid_t) :: group_id2     ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer(HSIZE_T) :: dims(3)
      integer(HSIZE_T) :: maxdims(3)
      integer     ::   b, d
   character(len=len(GROUP_BLOCK)+2) :: block_group
      nCell = 0
      
      inquire(file=trim(filename_data_in),exist=fexists)
      
      if(.not. fexists) then
        call error_wr("Data Input File konnte nicht gefunden werden: "//TRIM(filename_data_in),__FILE__,__LINE__)
      end if
      
      
      ! Initialize FORTRAN interface.
      call h5open_f(error)

      ! Open an existing file.
      call h5fopen_f (filename_data_in, h5f_acc_rdwr_f, file_id, error)

      call h5gopen_f(file_id,GROUP_GRID,group_id,error)

      call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
      if (nBlock > 1) then
         write(*,*) "Only one Block supported at the Moment"
         write(*,*) "Number of Blocks in File:",nBlock
         nBlock = 1
      end if
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
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

      blocks(b) % nCells = min(1,blocks(b) % nPkts - 1)
      nCell = nCell + product(blocks(b) % nCells)

      write(*,*) dimen,dims
      
      nFaces = dimen * 2
      nCorners = 2 ** dimen
      call allocate_vars()

      allocate (data_in(blocks(b)%nPkts(1),blocks(b)%nPkts(2),blocks(b)%nPkts(3)))
      do d = 1,dimen
         call h5dopen_f(group_id2, COORD_NAME(d), dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_in, dims, error)
         blocks(b) % coords (1:blocks(b)%nPkts(1) &
                            ,1:blocks(b)%nPkts(2) &
                            ,1:blocks(b)%nPkts(3),d) = data_in
         call h5dclose_f(dset_id, error)
      end do
      deallocate(data_in)
      end do
      stop 0 
!      
!      do ib = 1,nBlock
!      associate (b => blocks(ib))
!         cgns_zone = ib
!!         call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,cgns_zonename,isize,ierror)
!         if (ierror /= CG_OK) call cg_error_exit_f()
!      
!!         call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
!      
!         allocate (data_in(b%nPkts(1),b%nPkts(2),b%nPkts(3)))
!      
!         do d = 1,dimen
!!            call cg_coord_read_f(cgns_file,cgns_base,cgns_zone,COORD_NAME(d),RealDouble &
!!                                ,istart,isize(:,1),data_in,ierror)
!            if (ierror /= CG_OK) call cg_error_exit_f()
!            b % coords (1:b%nPkts(1) &
!                       ,1:b%nPkts(2) &
!                       ,1:b%nPkts(3),d) = data_in
!         end do
!         deallocate (data_in)
!          !!!! CHECKING if more than one solution.
!!          call cg_nsols_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,ierror)
!          if (ierror /= CG_OK) call cg_error_exit_f()
!          if (cgns_nSol /= 1) then
!             !call error_wr("More than One Solution in Restart-File",__FILE__,__LINE__)
!          end if
!          !!!! Checking if Cell-Centered
!!          call cg_sol_info_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,cgns_solname,data_location,ierror)
!          if (data_location .ne. CellCenter) then
!             call error_wr("Not Cell-Centered Data",__FILE__,__LINE__)
!          end if
!!          call cg_nfields_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,nVar_in,ierror)
!          allocate (data_in(b % nCells(1),b % nCells(2),b % nCells(3)))
!          do var = 1, nVar_in
!!             call cg_field_info_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,var,datatype,varname_in,ierror)
!             if (ierror /= CG_OK) call cg_error_exit_f()
!             call cg_field_read_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,VarName_in,datatype       &
!                                 ,istart,isize(:,2),data_in,ierror)
!             if (ierror /= CG_OK) call cg_error_exit_f()
!      
!              select case(VarName_in)
!                 case(VarName_Rho)
!                    b % vars(1:b % nCells(1)            &
!                            ,1:b % nCells(2)            &
!                            ,1:b % nCells(3),1) = data_in
!                 case(VarName_SpU)
!                    b % vars(1:b % nCells(1)            &
!                            ,1:b % nCells(2)            &
!                            ,1:b % nCells(3),2) = data_in
!                 case(VarName_SpV)
!                    b % vars(1:b % nCells(1)            &
!                            ,1:b % nCells(2)            &
!                            ,1:b % nCells(3),3) = data_in
!                 case(VarName_SpW)
!                    b % vars(1:b % nCells(1)            &
!                            ,1:b % nCells(2)            &
!                            ,1:b % nCells(3),4) = data_in
!                 case(VarName_Ene)
!                    b % vars(1:b % nCells(1)            &
!                            ,1:b % nCells(2)            &
!                            ,1:b % nCells(3),5) = data_in
!                 case default
!                    write(*,*) varname_in," nicht erkannt"
!               end select
!          end do
!          deallocate(data_in)
!         !  find out number of BCs that exist under this zone
!          call cg_nbocos_f(cgns_file,cgns_base,cgns_zone,cgns_nBoco,ierror)
!          if (cgns_nBoco /= 6) then
!            call error_wr(" Not enough Boundary Conditions in CGNS File",__FILE__,__LINE__) 
!          end if 
!         !  do loop over the total number of BCs
!          do ibb = 1, cgns_nBoco
!         !  get BC info
!!            call cg_boco_info_f(cgns_file,cgns_base,cgns_zone,ibb,               &
!!                 cgns_boconame,ibocotype,iptset,npts,normalindex,normallistflag,        &
!!                 normaldatatype,ndataset,ierror)
!            if (iptset .ne. PointRange) then
!              call error_wr("Error.  For this program, BCs must be set"// &
!                            " up as PointRange type"//&
!                            trim(PointSetTypeName(iptset)),__FILE__,__LINE__)
!            end if
!            !  read point range in here
!!            call cg_boco_read_f(cgns_file,cgns_base,cgns_zone,ibb,               &
!!                 ipnts,normallist,ierror)
!            dir_of_boundary = -1
!            if      (ipnts(1,1) == 1 .and. ipnts(1,2) == 1 ) then
!               dir_of_boundary = DIR_WEST
!            else if (ipnts(1,1) == b % nCells(1) .and. ipnts(1,2) == b % nCells(1) ) then
!               dir_of_boundary = DIR_EAST
!            else if (ipnts(2,1) == 1 .and. ipnts(2,2) == 1 ) then
!               dir_of_boundary = DIR_SOUTH
!            else if (ipnts(2,1) == b % nCells(2) .and. ipnts(2,2) == b % nCells(2) ) then
!               dir_of_boundary = DIR_NORTH
!            else if (ipnts(3,1) == 1 .and. ipnts(3,2) == 1 ) then
!               dir_of_boundary = DIR_FRONT
!            else if (ipnts(3,1) == b % nCells(3) .and. ipnts(3,2) == b % nCells(3) ) then
!               dir_of_boundary = DIR_BACK
!            else 
!               call error_wr("Boundary not in any direcion",__FILE__,__LINE__)
!            end if
!            write(6,'(a6,": ",a15," range=",3(1x,2i5))') &
!                   DIR_NAMES(dir_of_boundary), trim(BCTypeName(ibocotype)) &
!                  ,ipnts(1,1),ipnts(1,2) &
!                  ,ipnts(2,1),ipnts(2,2) &
!                  ,ipnts(3,1),ipnts(3,2)
!            select case (ibocotype) 
!            case (BCInflow)
!               b % boundary(dir_of_boundary) % bc_type = BC_INFLOW
!            case (BCOutflow)
!               b % boundary(dir_of_boundary) % bc_type = BC_OUTFLOW
!            case (BCWall)
!               b % boundary(dir_of_boundary) % bc_type = BC_WALL
!            case (BCSymmetryPlane)
!               b % boundary(dir_of_boundary) % bc_type = BC_SYMMETRY
!            case (BCGeneral)
!               b % boundary(dir_of_boundary) % bc_type = BC_PERIODIC
!            case default
!               call error_wr("BYType unknown: "//trim(BCTypeName(ibocotype)),__FILE__,__LINE__) 
!            end select
!          enddo
!          if ( b % nCells(3) == 1) then ! 2D SIMULATION
!            b % boundary(DIR_FRONT) % bc_type = BC_SYMMETRY
!            b % boundary(DIR_BACK)  % bc_type = BC_SYMMETRY
!          end if
!       end associate
!       end do
!      if (ierror /= CG_OK) call cg_error_exit_f()
   end subroutine datin_sol
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
      character(len=32)  :: solname
      integer :: cgns_file,cgns_coord,cgns_sol,cgns_var,ierror
      integer :: cgns_base,cgns_zone
!      integer(kind=CGSIZE_T) :: nVar_in
      
      integer(kind=CGSIZE_T),allocatable :: isize(:,:)
      character(len=32)  :: cgns_zonename
      character(len=32)  :: cgns_varname
      character(len=32),parameter  :: cgns_basename = "SOLUTION"
!      real(kind=dp) :: T
      real(REAL_KIND), allocatable :: data_out(:,:,:)

   write(solname,'(I0)') current_iteration
   if (write_sol_header) then
      !   call wr("Writing New Solution to File "//trim(file_data_out),4)
      if (.not. override_sol) write_sol_header = .false.

      call cg_open_f(trim(filename_data_out),CG_MODE_WRITE,cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      call cg_base_write_f(cgns_file,cgns_basename,dimen,dimen,cgns_base,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      allocate(isize(Dimen,3))
      isize = 0

      do b = 1, nBlock

         isize(:,1) = blocks(b) % nPkts  (1:Dimen)
         isize(:,2) = blocks(b) % nCells (1:Dimen)
         write(cgns_zonename,'(A,I3.3)') "BLOCK",b
         call cg_zone_write_f(cgns_file,cgns_base,cgns_zonename,isize,Structured,cgns_zone,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         allocate (data_out ( blocks(b) % nPkts(1), blocks(b) % nPkts(2), blocks(b) % nPkts(3) ))
         do d = 1, Dimen

            data_out = blocks(b) % coords(1:blocks(b) % nPkts(1)            &
                                        ,1:blocks(b) % nPkts(2)            &
                                        ,1:blocks(b) % nPkts(3),d)
            call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,coord_name(d) &
                                 ,data_out,cgns_coord,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()
         end do
         deallocate(data_out)

      end do
      call cg_biter_write_f(cgns_file,cgns_base,"TimeIterValues",1,ierror)


   else
!     call wr("Appending Solution to existing File: "//trim(file_data_out),4)
      call cg_open_f(trim(filename_data_out),CG_MODE_MODIFY,cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
!     call wr("reading nbases",4)
      call cg_nbases_f(cgns_file,cgns_base,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

!      call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
!      if (ierror /= CG_OK) call cg_error_exit_f()
!      if (zonetype /= Structured) then
!         call error_out("Only Structured Grid supported."//TRIM(file_git_in),__FILE__,__LINE__)
!      end if

!      call cg_biter_read_f(cgns_file,cgns_base,linkname,iter_in_file,ierror)
!      iter_in_file = iter_in_file + 1
!     write(*,*) "Iterations on File:",iter_in_file
!      call cg_biter_write_f(cgns_file,cgns_base,linkname,iter_in_file,ierror)

!      if (iter_in_file == 2) then
!
!      else
!         call cg_goto_f(cgns_file,cgns_base,ierror,'BaseIterativeData_t',1,'end')
!         allocate(iter_time(iter_in_file))
!
!         iter_time = dble(iteration)
!         call cg_array_write_f('TimeValues',RealDouble,1,2,iter_time,ierror)
!
!      end if
!      if (iter_in_file == 2) then
!         do b = 1,nBlock
!            cgns_zone = b
!            call cg_ziter_write_f(cgns_file,cgns_base,cgns_zone,'ZoneIterativeData',ierror)
!            call cg_goto_f(cgns_file,cgns_base,ierror,'Zone_t',cgns_zone,'ZoneIterativeData_t',1,'end')
!            call cg_array_write_f('FlowSolutionPointers',Character,1,32,solname,ierror)
!         end do
!      else
!      end if
   end if

   do b = 1,nBlock
      cgns_zone = b
      allocate (data_out(blocks(b) % nCells(1),blocks(b) % nCells(2),blocks(b) % nCells(3))) 
      call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,solname,CellCenter,cgns_sol,ierror)

      if (ierror /= CG_OK) call cg_error_exit_f()
      do var = 1, 6!sol_out_nVar
         select case(var)!VarName(var))
            case(1)!VarName_Rho)
               cgns_varname = VarName_Rho
               data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),1)
            case(2)!VarName_SpU)
               cgns_varname = VarName_SpU
               data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),2) 
            case(3)!VarName_SpV)
               cgns_varname = VarName_SpV
               data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),3)
            case(4)!VarName_SpW)
               cgns_varname = VarName_SpW
               data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),4)
            case(5)!VarName_Ene)
               cgns_varname = VarName_Ene
               data_out = blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),5)
!            case(VarName_Pre)
!               data_out = blocks(b) % P(1:block(b) % nCell(1)  &
!                                      ,1:block(b) % nCell(2)  &
!                                      ,1:block(b) % nCell(3))
!            case(VarName_Temp)
!               data_out = blocks(b) % T(1:block(b) % nCell(1)  &
!                                      ,1:block(b) % nCell(2)  &
!                                      ,1:block(b) % nCell(3))
            case(6)
               cgns_varname = "Mach Number"
               data_out = sqrt( &
                          blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),2) ** 2 &
                        + blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),3) ** 2 &
                        + blocks(b) % vars(1:blocks(b) % nCells(1)  &
                                          ,1:blocks(b) % nCells(2)  &
                                          ,1:blocks(b) % nCells(3),4) ** 2)&
                        / sqrt( gamma * RGas                      & 
                              * blocks(b) % temperatures(1:blocks(b) % nCells(1)  &
                                            ,1:blocks(b) % nCells(2)  &
                                            ,1:blocks(b) % nCells(3)))
            end select
            
            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble       &
                                    ,cgns_VarName,data_out,cgns_var,ierror)

       end do
       deallocate(data_out)
   end do
   call cg_close_f(cgns_file,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()
   end subroutine datout_sol
end module file_io_mod
