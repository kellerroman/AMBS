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
   use cgnslib
   use cgns_types, only: CGSIZE_T
   use screen_io_mod
implicit none
   integer, parameter               :: io_len_VarName = 20

   character ( len = io_len_Varname ), parameter :: VarName_Rho        = "Density"
   character ( len = io_len_Varname ), parameter :: VarName_SpU        = "Geschw_U"
   character ( len = io_len_Varname ), parameter :: VarName_SpV        = "Geschw_V"
   character ( len = io_len_Varname ), parameter :: VarName_SpW        = "Geschw_W"
   character ( len = io_len_Varname ), parameter :: VarName_Ene        = "Energie"
   character ( len = io_len_Varname ), parameter :: VarName_Pre        = "Druck"
   character ( len = io_len_Varname ), parameter :: VarName_Temp       = "Temperatur"
   character ( len = io_len_Varname ), parameter :: VarName_Mach       = "Mach_Number"

   character(len=32) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]

   character(len=*)  , parameter   :: FOLDER_INPUT       = "./"
   character(len=*)  , parameter   :: FOLDER_OUTPUT      = "./"
   character(len=100)              :: filename_control   = FOLDER_INPUT //"config.bin"
   character(len=100)              :: filename_data_in   = FOLDER_INPUT //"data_in.cgns"
   character(len=100)              :: filename_data_out  = FOLDER_OUTPUT//"data_out.cgns"
   
   integer                           , parameter :: CONFIG_FILE_VERSION = 1 
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
                        , timestep
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
      close(file_unit,status='delete')
      
      write(*,'(A30,2X,I0)') "Number of Iterations",max_iteration
      write(*,'(A30,2X,I0)') "Solution Output every",solution_out
      write(*,'(A30,2X,I0)') "Residual Output every",residual_out
      write(*,'(A30,2X,I0)') "Spatial Order",space_order
      write(*,'(A30,2X,ES10.4)') "Timestep",Timestep    
      write(*,'(A30,2X,A)') "Riemann Solver",trim(SOLVERS(riemann_solver))

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
   implicit none
      integer :: ib, d
      logical :: fexists
      integer(kind=CGSIZE_T) :: cgns_file,ierror,cgns_base,cgns_physDim,cgns_zone,zonetype
      integer(kind=CGSIZE_T) :: cgns_nSol,nVar_in,datatype,data_location,var,cgns_dimen
      integer(kind=CGSIZE_T) :: cgns_nBlock
      
      integer(kind=CGSIZE_T),allocatable :: isize(:,:),istart(:)
      character(len=32)  :: cgns_solname,varname_in
      character(len=32)  :: cgns_zonename
      character(len=32)  :: cgns_basename
!      real(kind=dp) :: T
      real(REAL_KIND), allocatable :: data_in(:,:,:)
      nCell = 0
      
      inquire(file=trim(filename_data_in),exist=fexists)
      
      if(.not. fexists) then
        call error_wr("Data Input File konnte nicht gefunden werden: "//TRIM(filename_data_in),__FILE__,__LINE__)
      end if
      
      
      call cg_open_f(trim(filename_data_in),CG_MODE_READ,cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
      
      call cg_nbases_f(cgns_file,cgns_base,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
     
      if (cgns_base /= 1) then
         call error_wr("Data Input File has more than one base",__FILE__,__LINE__)
      end if
      
      call cg_base_read_f(cgns_file,cgns_base,cgns_basename,cgns_dimen,cgns_physDim,ierror)
      dimen = int(cgns_dimen,INT_KIND)
      allocate(isize(dimen,3))
      allocate(istart(dimen))
      
      nFaces = dimen * 2
      nCorners = 2**dimen
      istart = 1
      call cg_nzones_f(cgns_file,cgns_base,cgns_nBlock,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
      nBlock = int(cgns_nBlock,INT_KIND)
      allocate(blocks(nBlock))
      
      do ib = 1,nBlock
      associate (b => blocks(ib))
         cgns_zone = ib
         call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,cgns_zonename,isize,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
      
         call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         if (zonetype /= Structured) then
            call error_wr("Only Structured Grid supported."//TRIM(filename_data_in),__FILE__,__LINE__)
         end if
         b % nPkts = 1
         b % nCells = 1
         b % nPkts (1:dimen) = int(isize(1:dimen,1),INT_KIND)
         b % nCells(1:dimen) = int(isize(1:dimen,2),INT_KIND)
         nCell = nCell + product(b % nCells)
      end associate
      end do
      
      call allocate_vars()
      
      do ib = 1,nBlock
      associate (b => blocks(ib))
         cgns_zone = ib
         call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,cgns_zonename,isize,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
      
         call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
      
         allocate (data_in(b%nPkts(1),b%nPkts(2),b%nPkts(3)))
      
         do d = 1,dimen
            call cg_coord_read_f(cgns_file,cgns_base,cgns_zone,COORD_NAME(d),RealDouble &
                                ,istart,isize(:,1),data_in,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()
            b % coords (1:b%nPkts(1) &
                       ,1:b%nPkts(2) &
                       ,1:b%nPkts(3),d) = data_in
         end do
         deallocate (data_in)
          !!!! CHECKING if more than one solution.
          call cg_nsols_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,ierror)
          if (ierror /= CG_OK) call cg_error_exit_f()
          if (cgns_nSol /= 1) then
             !call error_wr("More than One Solution in Restart-File",__FILE__,__LINE__)
          end if
          !!!! Checking if Cell-Centered
          call cg_sol_info_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,cgns_solname,data_location,ierror)
          if (data_location .ne. CellCenter) then
             call error_wr("Not Cell-Centered Data",__FILE__,__LINE__)
          end if
          call cg_nfields_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,nVar_in,ierror)
          allocate (data_in(b % nCells(1),b % nCells(2),b % nCells(3)))
          do var = 1, nVar_in
             call cg_field_info_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,var,datatype,varname_in,ierror)
             if (ierror /= CG_OK) call cg_error_exit_f()
             call cg_field_read_f(cgns_file,cgns_base,cgns_zone,cgns_nSol,VarName_in,datatype       &
                                 ,istart,isize(:,2),data_in,ierror)
             if (ierror /= CG_OK) call cg_error_exit_f()
      
              select case(VarName_in)
                 case(VarName_Rho)
                    b % vars(1:b % nCells(1)            &
                            ,1:b % nCells(2)            &
                            ,1:b % nCells(3),1) = data_in
                 case(VarName_SpU)
                    b % vars(1:b % nCells(1)            &
                            ,1:b % nCells(2)            &
                            ,1:b % nCells(3),2) = data_in
                 case(VarName_SpV)
                    b % vars(1:b % nCells(1)            &
                            ,1:b % nCells(2)            &
                            ,1:b % nCells(3),3) = data_in
                 case(VarName_SpW)
                    b % vars(1:b % nCells(1)            &
                            ,1:b % nCells(2)            &
                            ,1:b % nCells(3),4) = data_in
                 case(VarName_Ene)
                    b % vars(1:b % nCells(1)            &
                            ,1:b % nCells(2)            &
                            ,1:b % nCells(3),5) = data_in
                 case default
                    write(*,*) varname_in," nicht erkannt"
               end select
          end do
          deallocate(data_in)
      end associate
       end do
      call cg_close_f(cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
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
      integer(kind=CGSIZE_T) :: cgns_file,cgns_coord,cgns_sol,cgns_var,ierror
      integer(kind=CGSIZE_T) :: cgns_base,cgns_zone
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
      do var = 1, 5!sol_out_nVar
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
!            case(VarName_Mach)
!               data_out = blocks(b) % Q(1:block(b) % nCell(1)  &
!                                      ,1:block(b) % nCell(2)  &
!                                      ,1:block(b) % nCell(3),2) &
!                        / sqrt( gamma * RGas                      & 
!                              * block(b) % T(1:block(b) % nCell(1)  &
!                                            ,1:block(b) % nCell(2)  &
!                                            ,1:block(b) % nCell(3)))
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
