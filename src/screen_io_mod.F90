module screen_io_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: module for screen output
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 04.03.2016
! 
! LAST CHANGE: 17.10.2016
!
! CHANGELOG:
! 04.03.2016,RK: Start of Coding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use, intrinsic :: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT
   use const_mod
implicit none
   integer, parameter         :: SCREEN_WIDTH         = 100 
   character(len=*),parameter :: FORMAT_SEPLINE       = '(100("="),/,100("="))'
   character(len=*),parameter :: FORMAT_LINE          = '(3("="),1X,A92,1X,3("="))'
   character(len=*),parameter :: RED_START            = achar(27)//"[31m"
   character(len=*),parameter :: RED_END              = achar(27)//"[0m"                    
   character(len=*),parameter :: GREEN_START          = achar(27)//"[32m"
   character(len=*),parameter :: GREEN_END            = achar(27)//"[0m"                    
contains
subroutine screen_wr_start()
implicit none
   integer :: start_time(8)
   call date_and_time(values = start_time)
   write(stdout,'(A)') RED_START
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) "AMBS 2016                                                   "
   write(stdout,FORMAT_LINE) ""
   write(stdout,'(3("="),23X,"Date:",4X,I2.2,"-",I2.2,"-",I4.4,52X,3("="))') start_time(3),start_time(2),start_time(1)
   write(stdout,'(3("="),23X,"Time:",4X,I2.2,"-",I2.2,"-",I2.2,54X,3("="))') start_time(5),start_time(6),start_time(7)
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,'(A)') RED_END
end subroutine
subroutine screen_wr_start_calc()
implicit none
   write(stdout,'(A)') GREEN_START
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,FORMAT_LINE) "Start of Calculation                                       "
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,'(A)') GREEN_END
end subroutine
subroutine screen_wr_end()
   use control_mod, only: get_runtime
implicit none
   integer :: date_time(8)
   call date_and_time(values = date_time)
   write(stdout,'(A)') RED_START
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) "AMBS done!                                                 "
   write(stdout,FORMAT_LINE) ""
   write(stdout,'(3("="),23X,"Date:",4X,I2.2,"-",I2.2,"-",I4.4,52X,3("="))') date_time(3),date_time(2),date_time(1)
   write(stdout,'(3("="),23X,"Time:",4X,I2.2,"-",I2.2,"-",I2.2,54X,3("="))') date_time(5),date_time(6),date_time(7)
   write(stdout,'(3("="),23X,"Runtime:",1X,F10.2,"s",51X,3("="))') get_runtime()
                                                                  
                                                                  
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,'(A)') RED_END
end subroutine
subroutine screen_wr(text,level)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: standard screen output with optional level and Coloring
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
   implicit none
      character(len=*),intent(in)           :: text
      integer         ,intent(in), optional :: level
      character(len=*), parameter           :: PRE_TEXT  = "=================" 

      if (present(level)) then
         if (level == 1) then
            write(stdout,'(A)')  RED_START     // &
                                      PRE_TEXT      // &
                                      " "//text//" "// &
                                      PRE_TEXT      // &
                                      RED_END
         end if
      else
         write(stdout,'(A)') text
      end if
   end subroutine screen_wr

   subroutine error_wr(text,errorfile,errorline)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: standard screen output with optional level and Coloring
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
   implicit none
      character(len=*),intent(in)           :: text
      character(len=*),intent(in)           :: errorfile
      integer         ,intent(in)           :: errorline
      character(len=*), parameter           :: PRE_TEXT  = "=================" 

      write(stderr,'(A,I0,A)')  RED_START         // &
                                PRE_TEXT            // &
                               "IN "//errorfile     // &
                               "@ ",errorline        , &
                                PRE_TEXT            // &
                                RED_END
      write(stderr,'(A)')  RED_START            // &
                                PRE_TEXT            // &
                                " "//text//" "      // &
                                PRE_TEXT            // &
                                RED_END
      stop 1
   end subroutine error_wr

   subroutine screen_residual()
      use control_mod ,only: current_iteration,res_avg,res_max,solution_time, get_runtime, fio
   implicit none
      real(kind = REAL_KIND) :: runtime
      runtime = get_runtime()
      write(stdout,'(I10,6(1X,ES10.4))') current_iteration & 
                                           , res_max(1) &
                                           , res_avg(1) &
                                           , res_max(2) &
                                           , res_avg(2) &
                                           , solution_time &
                                           , runtime
      write(fio,'(I10,6(1X,ES10.4))') current_iteration & 
                                           , res_max(1) &
                                           , res_avg(1) &
                                           , res_max(2) &
                                           , res_avg(2) &
                                           , solution_time &
                                           , runtime
   end  subroutine screen_residual
end module screen_io_mod
