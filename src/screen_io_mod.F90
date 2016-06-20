module screen_io_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: module for screen output
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 04.03.2016
! 
! LAST CHANGE: 04.03.2016
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
   write(stdout,'(A)') RED_START
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) "AMBS by ROMAN KELLER                                        "
   write(stdout,FORMAT_LINE) ""
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
implicit none
   write(stdout,'(A)') RED_START
   write(stdout,FORMAT_SEPLINE) 
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) ""
   write(stdout,FORMAT_LINE) "AMBS done!                                                 "
   write(stdout,FORMAT_LINE) ""
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
      use control_mod ,only: current_iteration,res_avg,res_max,solution_time
   implicit none
      write(stdout,'(I10,3(1X,ES10.4))') current_iteration & 
                                           , res_max &
                                           , res_avg &
                                           , solution_time
   end  subroutine screen_residual
end module screen_io_mod
