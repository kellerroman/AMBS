!==================================================================================================!
!===     Ascii-Config Parser and Binary Config Writer                                           ===!
!===                                                                                            ===!
!===     author: Roman Keller                                                                   ===!
!===                                                                                            ===!
!===     date: 30.01.2017                                                                       ===!
!===                                                                                            ===!
!===                                                                                            ===!
!==================================================================================================!
module datin_para
implicit none
   private
   integer, parameter         :: CONFIG_FILE_VERSION     =  2 
   integer, parameter, public :: VARNAME_LENGTH          = 20
   integer, parameter, public :: RETURN_CODE_GOOD        =  0
   integer, parameter, public :: RETURN_CODE_NOT_IN_LIST = -1
   integer, parameter, public :: RETURN_CODE_NOT_A_INT   = -2
   integer, parameter, public :: RETURN_CODE_NOT_A_REAL  = -3
   integer, parameter         :: DATATYPE_INT            =  1
   integer, parameter         :: DATATYPE_REAL           =  2
   type :: t_datin_para
      character(len=VARNAME_LENGTH) :: varname
      logical :: required
      logical :: is_set
      integer :: datatype
      integer, pointer :: int_var
      real(kind=8), pointer :: real_var
      type(t_datin_para), pointer :: next => null()
   end type t_datin_para
   type(t_datin_para), pointer :: first => null()

   integer, protected :: n_datin_para = 0

   public :: add_parameter, list_parameter, set_para, get_unset_para, get_unset_paras, free_all, write_binary_config_file
   interface add_parameter
      module procedure add_int
      module procedure add_real
   end interface add_parameter
   contains
      subroutine add_int(str,var_pointer,init_value)
      implicit none
      character(len=*), intent(in) :: str
      integer, target, intent(in) :: var_pointer
      integer, intent(in), optional :: init_value

      type(t_datin_para), pointer :: new, tmp
      ! Speicher reservieren
      allocate(new)
      ! Werte setzen
      new % varname = lower_case(str)
      new % datatype = DATATYPE_INT
      new % int_var => var_pointer
      if ( present(init_value) ) then
         new % int_var = init_value
         new % is_set   = .true.
         new % required = .false.
      else
         new % required = .true.
         new % is_set   = .false.
      end if
      ! Am Beginn der Liste einf?gen
      if (.not. associated(first)) then
         first => new
      else
         tmp => first
         first => new
         first % next => tmp
       end if
       n_datin_para = n_datin_para + 1
      end subroutine add_int

      subroutine add_real(str,var_pointer,init_value)
      implicit none
      character(len=*), intent(in) :: str
      real(kind=8), target, intent(in) :: var_pointer
      real(kind=8), intent(in), optional :: init_value

      type(t_datin_para), pointer :: new, tmp
      ! Speicher reservieren
      allocate(new)
      ! Werte setzen
      new % varname = lower_case(str)
      new % datatype = DATATYPE_REAL
      new % real_var => var_pointer
      if ( present(init_value) ) then
         new % real_var = init_value
         new % is_set   = .true.
         new % required = .false.
      else
         new % required = .true.
         new % is_set   = .false.
      end if
      ! Am Beginn der Liste einf?gen
      if (.not. associated(first)) then
         first => new
      else
         tmp => first
         first => new
         first % next => tmp
       end if
       n_datin_para = n_datin_para + 1
      end subroutine add_real

      subroutine list_parameter()
      implicit none
      integer :: npara
      type(t_datin_para), pointer :: tmp
      npara = 0
      tmp => first
      do 
         if (.not. associated(tmp) ) exit
         npara = npara + 1
         if ( tmp % datatype == DATATYPE_INT) then
            if (tmp % is_set .or. .not. tmp % required) then 
               write(*,'(A20," = ",I0)') tmp % varname, tmp % int_var
            else
               write(*,'(A20      )') tmp % varname
            end if
         else if ( tmp % datatype == DATATYPE_REAL) then
            if (tmp % is_set .or. .not. tmp % required) then 
               write(*,'(A20," = ",ES10.4)') tmp % varname, tmp % real_var
            else
               write(*,'(A20,1A       )') tmp % varname
            end if
         end if 
         tmp =>tmp % next
      end do

      end subroutine list_parameter


      subroutine write_binary_config_file()
      implicit none
      integer :: npara
      integer :: file_unit
      type(t_datin_para), pointer :: tmp
      npara = 0
      open(newunit = file_unit, file="config.bin",form="unformatted",access="stream")
      write(file_unit) CONFIG_FILE_VERSION
      tmp => first
      do 
         if (.not. associated(tmp) ) exit
         npara = npara + 1
         if ( tmp % datatype == DATATYPE_INT) then
            write(file_unit) tmp % int_var
         else if ( tmp % datatype == DATATYPE_REAL) then
            write(file_unit) tmp % real_var
         end if 
         tmp =>tmp % next
      end do
      close(file_unit)

      end subroutine write_binary_config_file


      function set_para(varname,varvalue) result (return_value)
      implicit none
      character(len=VARNAME_LENGTH), intent(in) :: varname
      character(len=*), intent(in) :: varvalue
      integer :: return_value
      type(t_datin_para), pointer :: tmp

      return_value = RETURN_CODE_GOOD
      tmp => first
      do 
         if (.not. associated(tmp) ) then
            return_value = RETURN_CODE_NOT_IN_LIST
            exit
         end if
   
         if (tmp % varname == lower_case(varname)) then
            if (tmp % datatype == DATATYPE_INT) then
               if (is_integer(varvalue)) then
                  !write(*,*) "Changing ",trim(varname)," to value ",varvalue
                  read(varvalue,*) tmp % int_var
                  tmp % is_set = .true.
               else
                  return_value = RETURN_CODE_NOT_A_INT
                  exit
               end if
            else if (tmp % datatype == DATATYPE_REAL) then
               if (is_real(varvalue)) then
                  !write(*,*) "Changing ",trim(varname)," to value ",varvalue
                  read(varvalue,*) tmp % real_var
                  tmp % is_set = .true.
               else
                  return_value = RETURN_CODE_NOT_A_REAL
                  exit
               end if
            end if
            exit
         end if
         
         tmp =>tmp % next
      end do

      
      end function set_para

      function get_unset_para(last_para) result (para_name)
      implicit none
      character(len=VARNAME_LENGTH), intent(in) :: last_para
      character(len=VARNAME_LENGTH) :: para_name
      type(t_datin_para), pointer :: tmp
      tmp => first
      para_name = ""
      do 
         if (.not. associated(tmp) ) exit
         if (    tmp % required              .and. &
            .not.tmp % is_set                .and. &
                 tmp % varname /= last_para) then
            para_name = tmp % varname
            exit
         end if
         tmp =>tmp % next
      end do

      end function get_unset_para

      function get_unset_paras() result (para_name)
      implicit none
      character(len=VARNAME_LENGTH),allocatable :: para_name(:)
      type(t_datin_para), pointer :: tmp
      integer :: nunset
      tmp => first
      nunset = 0
      do 
         if (.not. associated(tmp) ) exit
         if (    tmp % required              .and. &
            .not.tmp % is_set             ) then
            nunset = nunset + 1
         end if
         tmp =>tmp % next
      end do
      allocate(para_name(nunset))
      nunset = 0
      tmp => first
      do 
         if (.not. associated(tmp) ) exit
         if (    tmp % required              .and. &
            .not.tmp % is_set             ) then
            nunset = nunset + 1
            para_name(nunset) = tmp % varname
         end if
         tmp =>tmp % next
      end do
      end function get_unset_paras
      subroutine free_all()
      implicit none
      type(t_datin_para), pointer :: tmp
      do        
         tmp => first
         if (.not. associated(tmp)) exit
         first => first%next
         deallocate(tmp)
         end do                     
      end subroutine free_all

      function lower_case( input_string ) result ( output_string )
         ! -- argument and result
         implicit none
         character( * ), intent( in )       :: input_string
         character( len( input_string ) )   :: output_string
         integer                            :: ii,ic,nlen
         nlen = len(input_string)
         do ii=1,nlen
            ic = ichar(input_string(ii:ii))
            if (ic >= 65 .and. ic <= 90) then
               output_string(ii:ii) = char(ic+32)
            else
               output_string(ii:ii) = char(ic)
            end if
         end do
      end function lower_case
      logical function is_integer(string)
         implicit none
            integer :: ipos
            character (len = *) :: string
            is_integer =  .true.
            do ipos = 1,len_trim(string)
               if ((ichar(string(ipos:ipos)) >57 .or. ichar(string(ipos:ipos)) <48 )&
                  .and.string(ipos:ipos) /= "-" .and.string(ipos:ipos) /= "+") then
      !            write(*,*) trim(string) , "ist keine integer-zahl"
                  is_integer =  .false.
                  return
               end if
            end do

         end function

      logical function is_real(string)
      implicit none
         integer :: ipos
         integer :: elem_pos
         character (len = *) :: string
         is_real =  .true.
         elem_pos = 1
         !!!! aufteilung eines realen zahl   - 0 . 0 d + 1
         ! vorzeichen                  elem_pos = 1
         ! zahl vor dem komma          elem_pos = 2
         ! trennzeichen                elem_pos = 3
         ! zahl nach dem komma         elem_pos = 4
         ! exponent                    elem_pos = 5
         ! exponenten-vorzeichen       elem_pos = 6
         ! exponent                    elem_pos = 7
         do ipos = 1,len_trim(string)
            if (ichar(string(ipos:ipos)) <=57 .and. ichar(string(ipos:ipos)) >=48) then
   !            write(*,*) trim(string) , "ist keine integer-zahl"
               select case(elem_pos)
                  case (2)
                  case (4)
                  case (5)
                     elem_pos = 7
                  case (7)
                  case default
                     elem_pos = elem_pos + 1
               end select
            else if (string(ipos:ipos) == "-" .or. string(ipos:ipos) == "+") then
               select case(elem_pos)
                  case (2:4)
                     is_real = .false.
                     return
                  case (6:7)
                     is_real = .false.
                     return
                  case default
                     elem_pos = elem_pos + 1
               end select
            else if (string(ipos:ipos) == ".") then
               select case(elem_pos)
                  case (1:2)
                     elem_pos = elem_pos + 1
                  case default
                     is_real = .false.
                     return
               end select
            else if (string(ipos:ipos) == "d" .or. string(ipos:ipos) == "e" .or.  &
                     string(ipos:ipos) == "D" .or. string(ipos:ipos) == "E") then
               select case(elem_pos)
                  case (2:4)
                     elem_pos = 5
                  case default
                     is_real = .false.
                     return
               end select
            end if
         end do
         if (elem_pos <= 2 .or. &
             elem_pos == 5 .or. &
             elem_pos == 6 ) then
            is_real = .false.
            return
         end if
      end function
end module datin_para
program datin
   use datin_para
implicit none
integer :: fu
integer :: stat
integer :: pos
character(len=100) :: line
character(len=VARNAME_LENGTH) :: varname
character(len=90) :: varvalue

integer :: iterations, sol_out_screen_int, res_out_screen_int, res_nr_out
integer :: equation, turbulence, space_order
integer :: riemann_solver, timestep_method,time_order
real(kind = 8) :: c_les_sgs, cfl, timestep
real(kind = 8) :: DT
character(len=VARNAME_LENGTH),allocatable :: unset_paras(:)
integer :: i

write(*,*) "START DATIN"
!!!!! LIST WIRD IN UMGEKEHRTER REIHENFOLGE GESCHRIEBEN!!
!call add_parameter("RES_NR_OUT",res_nr_out)
call add_parameter("C_LES_SGS",c_les_sgs)
call add_parameter("TIMESTEP",timestep)
call add_parameter("CFL",cfl)
call add_parameter("TIME_ORDER",time_order)
call add_parameter("TIMESTEP_METHOD",timestep_method)
call add_parameter("RIEMANN_SOLVER",riemann_solver)
call add_parameter("SPACE_ORDER",space_order)
call add_parameter("TURBULENCE",turbulence)
call add_parameter("EQUATION",equation)
call add_parameter("RES_OUT_SCREEN_INT",res_out_screen_int)
call add_parameter("SOL_OUT_SCREEN_INT",sol_out_screen_int)
call add_parameter("ITERATIONS",iterations)

open(newunit=fu,file="config.cfg")
do 
   read(fu,'(A)',iostat=stat) line

   if (stat < 0) then
      !write(*,*) "End of File"
      exit
   end if
   line = trim(adjustl(line))
   pos = index(line,"!")
   if (pos > 0 ) then
      line = line(1:pos-1)
   end if
   if (len_trim(line) == 0) cycle
   !write(*,*) line
   pos = index(line,"=")

   if (pos < 1) then
      write(*,'("konnte line nicht verarbeiten, kein = gefunden:",a)') trim(line)
      stop 1
   end if
!   write(*,*) trim(line(1:pos-1)),trim(line(pos+1:))
   varname = trim(line(1:pos-1))
   varvalue = trim(adjustl(line(pos+1:)))
   select case(set_para(varname,varvalue))
   case(RETURN_CODE_NOT_IN_LIST)
      !write(*,*) "'"//trim(varname)//"' is not in list"
      !stop 1
   case(RETURN_CODE_NOT_A_INT)
      write(*,*) "value for '"//trim(varname)//"' is not an INTEGER ",varvalue
      stop 1
   case(RETURN_CODE_NOT_A_REAL)
      write(*,*) "value for '"//trim(varname)//"' is not an REAL ",varvalue
      stop 1
   end select
   
end do
close(fu)
varname = ""
!do
!   varname = get_unset_para(varname)
!   if (varname == "") then
!      exit
!   else
!      write(*,*) varname, " must be set"
!      stop 1
!   end if
!end do
unset_paras =  get_unset_paras()
if (ubound(unset_paras,1) > 0 ) then
   write(*,*) "There are unset Parameters which require Values"
   do i = 1,ubound(unset_paras,1)
      write(*,*) i, unset_paras(i)
   end do
   stop 1
end if
call list_parameter()
call write_binary_config_file()
call free_all()
write(*,*) "ENDE  DATIN"

end program
