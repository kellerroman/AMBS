module screen_io_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use ISO_FORTRAN_ENV, only: OUTPUT_UNIT, ERROR_UNIT
   use const_mod
implicit none
   character(len=*),parameter :: RED_START = achar(27)//"[31m"
   character(len=*),parameter :: RED_END   = achar(27)//"[0m"
contains
   subroutine screen_wr(text,level)
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
      integer         ,intent(in), optional :: level
      character(len=*), parameter           :: PRE_TEXT  = "=================" 

      if (present(level)) then
         if (level == 1) then
            write(OUTPUT_UNIT,'(A)')  RED_START     // &
                                      PRE_TEXT      // &
                                      " "//text//" "// &
                                      PRE_TEXT      // &
                                      RED_END
         end if
      else
         write(OUTPUT_UNIT,'(A)') text
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

      write(ERROR_UNIT,'(A,I0,A)')  RED_START         // &
                                PRE_TEXT            // &
                               "IN "//errorfile     // &
                               "@ ",errorline        , &
                                PRE_TEXT            // &
                                RED_END
      write(ERROR_UNIT,'(A)')  RED_START            // &
                                PRE_TEXT            // &
                                " "//text//" "      // &
                                PRE_TEXT            // &
                                RED_END
      stop 1
   end subroutine error_wr

   subroutine screen_residual()
      use control_mod ,only: current_iteration,res_avg,res_max,solution_time
   implicit none
      write(OUTPUT_UNIT,'(I10,3(1X,ES10.4))') current_iteration & 
                                           , res_max &
                                           , res_avg &
                                           , solution_time
   end  subroutine screen_residual

end module screen_io_mod
