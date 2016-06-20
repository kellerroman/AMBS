program extract_1D
use const_mod
use data_mod
!use control
use cgns
use file_io_mod
implicit none
integer, parameter :: fo = 12345
integer :: b, i ,j ,k ,d

integer :: i_start,i_end
integer :: j_start,j_end

logical :: fexists
integer(kind=CGSIZE_T) :: cgns_file,ierror,cgns_base,PhysDim,cgns_zone,zonetype
integer(kind=CGSIZE_T) :: nSol,nVar_in,datatype,data_location,var,sol
character(len = 100) :: arg
integer(kind=CGSIZE_T),allocatable :: isize(:,:),istart(:)
character(len=32)  :: solname,cgns_git_basename
character(len=*), parameter :: file_data_in = "data_out.cgns"
character(len=32),allocatable  :: varname_in(:)
character(len=*), parameter :: file_out = "data.vtk"
real(REAL_KIND), allocatable :: temp_coord(:,:,:)
integer :: var_dir = 1
!< Welche Richtung soll variert werden

write(*,*) "========== SOLUTION 2 VTK  =================="
nCell = 0
i = 1

inquire(file=trim(file_data_in),exist=fexists)

if(.not. fexists) then
  write(*,*) "Data Input File konnte nicht gefunden werden: "//TRIM(file_data_in),__FILE__,__LINE__
end if


call cg_open_f(trim(file_data_in),CG_MODE_READ,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_nbases_f(cgns_file,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

if (cgns_base /= 1) then
   write(*,*) "Data Input File has more than one base",__FILE__,__LINE__
end if

call cg_base_read_f(cgns_file,cgns_base,cgns_git_basename,Dimen,PhysDim,ierror)

allocate(isize(Dimen,3))
allocate(istart(Dimen))

nFaces = Dimen * 2
nCorners = 2**Dimen
istart = 1

call cg_nzones_f(cgns_file,cgns_base,nBlock,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
allocate(blocks(nBlock))
write(*,*) "Number of Blocks:", nBlock
Write(*,'(A4,3(1X,A4))') "#Bl","NI","NJ","NK"
do b = 1,nBlock
   cgns_zone = b
   call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,solname,isize,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()

   call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()
   if (zonetype /= Structured) then
      write(*,*) "Only Structured Grid supported."//TRIM(file_data_in),__FILE__,__LINE__
   end if
   blocks(b) % nPkts = 1
   blocks(b) % nCells = 1
   blocks(b) % nPkts(1:Dimen) = int(isize(1:Dimen,1),4)
   blocks(b) % nCells(1:Dimen) = int(isize(1:Dimen,2),4)
   write(*,'(I4,3(1X,I4))') b, blocks(b) % nCells
   nCell = nCell + product(blocks(b) % nCells)
   allocate (blocks(b) % coords(1 : blocks(b)%nPkts(1) &
                           ,1 : blocks(b)%nPkts(2) &
                           ,1 : blocks(b)%nPkts(3) &
                           ,Dimen))

end do
write(*,*) "Overall Cell Number:", nCell
do b = 1,nBlock
   cgns_zone = b
   call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,solname,isize,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()

   call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)

   allocate (temp_coord(blocks(b)%nPkts(1),blocks(b)%nPkts(2),blocks(b)%nPkts(3)))

   do d = 1,Dimen
      call cg_coord_read_f(cgns_file,cgns_base,cgns_zone,coord_name(d),RealDouble &
                          ,istart,isize(:,1),temp_coord,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
      blocks(b) % coords(1:blocks(b)%nPkts(1) &
                        ,1:blocks(b)%nPkts(2) &
                        ,1:blocks(b)%nPkts(3),d) = temp_coord
   end do
   deallocate (temp_coord)

   !!!! CHECKING if more than one solution.
   call cg_nsols_f(cgns_file,cgns_base,cgns_zone,nSol,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()
   write(*,*) "Number of Solutions on File:", nSol
   !!!! Checking if Cell-Centered
   call cg_sol_info_f(cgns_file,cgns_base,cgns_zone,1,solname,data_location,ierror)
   if (data_location .ne. CellCenter) then
      write(*,*) "Not Cell-Centered Data",__FILE__,__LINE__
      stop 1
   end if
   call cg_nfields_f(cgns_file,cgns_base,cgns_zone,1,nVar_in,ierror)
   write(*,*) "Number of Variables on File",nVar_in
   allocate ( blocks(b) % vars      (blocks(b)%nCells(1) &
                                    ,blocks(b)%nCells(2) &
                                    ,blocks(b)%nCells(3) &
                                    ,nVar_in*nSol))
   allocate (varname_in(nVar_in))
   do sol = 1, nSol
      do var = 1, nVar_in
         call cg_field_info_f(cgns_file,cgns_base,cgns_zone,sol,var,datatype,varname_in(var),ierror)
         call cg_field_read_f(cgns_file,cgns_base,cgns_zone,sol,VarName_in(var),datatype       &
                             ,istart,isize(:,2) &
                             ,blocks(b)%vars(:,:,:,var+(sol-1)*nVar_in) &
                             ,ierror)
      end do
   end do
end do
call cg_close_f(cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

write(*,*) "========== SOLUTION 2 VTK  =================="
contains
    function vtk_blk_data(block,filename) result(err)
    !---------------------------------------------------------------------------------------------------------------------------
    ! Function for writing block data.
    !---------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(block_type), intent(IN):: block    ! Block-level data.
    character(*),      intent(IN)::    filename ! File name of the output block file.
    integer(I_P)::                     err      ! Error trapping flag: 0 no errors, >0 error occurs.
    integer(I_P)::         ni1,ni2,nj1,nj2,nk1,nk2                         ! Bounds of dimensions of node-centered data.
    integer(I_P)::         ci1,ci2,cj1,cj2,ck1,ck2                         ! Bounds of dimensions of cell-centered data.
    integer(I_P)::         nnode,ncell                                     ! Number of nodes and cells.
    integer(I_P)::         i,j,k,s                                         ! Counters.
    !---------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------
    if (pp_format%node) then
      ! tri-linear interpolation of cell-centered values at nodes
      call interpolate_primitive(block=block,P=P)
    endif
    ! initialize the block dimensions
    call compute_dimensions(block=block,                                     &
                            ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2, &
                            ci1=ci1,ci2=ci2,cj1=cj1,cj2=cj2,ck1=ck1,ck2=ck2)
    nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1)
    ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1)
    ! storing species densities into an array for avoiding problems with nonzero rank pointer
    ! initializing VTK file
    err = VTK_INI_XML(output_format = 'binary',         &
                      filename      = trim(filename),   &
                      mesh_topology = 'StructuredGrid', &
                      nx1 = ni1, nx2 = ni2,             &
                      ny1 = nj1, ny2 = nj2,             &
                      nz1 = nk1, nz2 = nk2)
    ! saving auxiliary data (time and time step)
    err = VTK_FLD_XML(fld_action='open')
    err = VTK_FLD_XML(fld=block%global%t,fname='TIME')
    err = VTK_FLD_XML(fld=block%global%n,fname='CYCLE')
    err = VTK_FLD_XML(fld_action='close')
    ! saving the geometry
    err = VTK_GEO_XML(nx1 = ni1, nx2 = ni2,                    &
                      ny1 = nj1, ny2 = nj2,                    &
                      nz1 = nk1, nz2 = nk2,                    &
                      NN = nnode,                              &
                      X=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%x, &
                      Y=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%y, &
                      Z=block%node(ni1:ni2,nj1:nj2,nk1:nk2)%z)
    if (.not.meshonly) then
      ! saving dependent variables
      err = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
      do s=1,global%Ns
        err=VTK_VAR_XML(NC_NN=ncell,varname='r('//trim(str(.true.,s))//')',var=r(s,ci1:ci2,cj1:cj2,ck1:ck2))
      enddo
      err=VTK_VAR_XML(NC_NN=ncell,varname='r',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%d  )
      err=VTK_VAR_XML(NC_NN=ncell,varname='u',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%x)
      err=VTK_VAR_XML(NC_NN=ncell,varname='v',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%y)
      err=VTK_VAR_XML(NC_NN=ncell,varname='w',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%v%z)
      err=VTK_VAR_XML(NC_NN=ncell,varname='p',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%p  )
      err=VTK_VAR_XML(NC_NN=ncell,varname='g',var=block%C(ci1:ci2,cj1:cj2,ck1:ck2)%P%g  )
      err=VTK_DAT_XML(var_location ='cell',var_block_action = 'close')
    endif
    ! closing VTK file
    err = VTK_GEO_XML()
    err = VTK_END_XML()
    return
    !---------------------------------------------------------------------------------------------------------------------------
    endfunction vtk_blk_data
end program extract_1D


