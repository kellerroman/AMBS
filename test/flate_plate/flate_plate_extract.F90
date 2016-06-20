program flate_plate_extract
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
integer(kind=CGSIZE_T),allocatable :: isize(:,:),istart(:)
character(len=32)  :: solname,cgns_git_basename
character(len=*), parameter :: file_data_in = "data_out.cgns"
character(len=32),allocatable  :: varname_in(:)
character(len=*), parameter :: file_out = "data_1d.csv"
real(REAL_KIND), allocatable :: temp_coord(:,:,:)
real(REAL_KIND), allocatable :: grenzschicht_dicke(:)
real(REAL_KIND) :: avg_vel,y1,y2,v1,v2
!< Welche Richtung soll variert werden

write(*,*) "========== EXTRACT 1D SOLUTION =============="
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
if (nBlock > 1) then
   write(*,*) "Extract 1D nur möglich mit einem Block:",nBlock
   stop 1
end if
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

b = 1
i = 1
k = 1
avg_vel = 0.0E0_REAL_KIND
do j = 1, blocks(b) % nCells(2)
   avg_vel = avg_vel + blocks(b) % vars(i,j,k,2)
end do
avg_vel = avg_vel / blocks(b) % nCells(2) 
write(*,*) avg_vel
allocate( grenzschicht_dicke(blocks(b) % nCells(1)))
do i = 1, blocks(b) % nCells(1)
!   avg_vel = 0.0E0_REAL_KIND
!   do j = 1, blocks(b) % nCells(2)
!      avg_vel = avg_vel + blocks(b) % vars(i,j,k,2)
!   end do
!   avg_vel = avg_vel / blocks(b) % nCells(2) 
   do j = 1, blocks(b) % nCells(2)
      if  ( blocks(b) % vars(i,j,k,2) >= 0.99E0_REAL_KIND * avg_vel) then
         if (j > 1) then
            y1 = 0.5E0_REAL_KIND * ( blocks(b) % coords(i,j,k,2) + blocks(b) % coords(i,j-1,k,2))
            y2 = 0.5E0_REAL_KIND * ( blocks(b) % coords(i,j,k,2) + blocks(b) % coords(i,j+1,k,2))
            v1 = blocks(b) % vars(i,j-1,k,2)
            v2 = blocks(b) % vars(i,j  ,k,2)
            grenzschicht_dicke(i) = y1 + (y2-y1) / (v2-v1) * (avg_vel-v1)
         else
            grenzschicht_dicke(i) = 0.0E0_REAL_KIND 
            !* ( blocks(b) % coords(i,j,k,2) + blocks(b) % coords(k,j+1,k,2))
         end if
         exit
      end if
   end do
   write(*,*) i,j, grenzschicht_dicke(i), blocks(b) % vars(i,j-1:j+1,k,2),avg_vel
end do

open(unit = fo , file = trim(file_out))
   i_start = 1
   i_end   = blocks(b) % nCells(1)
   j_start = 1 !max(1,blocks(b) % nCells(2) / 2)
   j_end   = 1 !max(1,blocks(b) % nCells(2) / 2)
write(fo,'(A20)',advance="no") "COORD"
do var = 1,nVar_in
   write(fo,'(",",A20)',advance="no") trim(VarName_in(var))
end do
write(fo,*)
do i = i_start, i_end
   do j = j_start, j_end
      write(fo ,'(F20.13)',advance = "no") (blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i+1,j,k,1)) * 0.5d0
      do var = 1,nVar_in*nSol
         write(fo,'(",",F20.13)',advance = "no") blocks(b) % vars(i,j,k,var)
      end do
      write(fo,*)
   end do
end do
write (fo,*) 
close(fo)
open(unit = fo , file = "grenzschicht_dicke.csv")
j = 1
k = 1
do i = 1, blocks(b) % nCells(1) 
   write(fo,*) (blocks(b) % coords(i,j,k,1)+blocks(b) % coords(i+1,j,k,1)) * 0.5d0, grenzschicht_dicke(i) 
end do
deallocate ( grenzschicht_dicke)
close(fo)




open(unit = fo , file = "profile.csv")
j = 1
i = 1
k = 1
do j = blocks(b) % nCells(2),1,-1 
write(fo,*) (blocks(b) % coords(i,j,k,2)+blocks(b) % coords(i,j+1,k,2)) * 0.5d0,blocks(b) % vars(:,j,k,2) 
end do
close(fo)
write(*,*) "========== EXTRACT 1D SOLUTION done ========="
end program flate_plate_extract
