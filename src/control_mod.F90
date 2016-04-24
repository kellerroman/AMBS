module control_mod
   use const_mod
implicit none

logical :: end_main_loop
logical :: residual_on_screen
logical :: solution_to_file
integer(INT_KIND) :: max_iteration = 10000
integer(INT_KIND) :: current_iteration = 0
integer(INT_KIND) :: nBoundaryCells = 2
integer(INT_KIND) :: space_order = 1
real(REAL_KIND) :: res_max
real(REAL_KIND) :: res_avg
real(REAL_KIND) :: solution_time = 0.0E+0_REAL_KIND
real(REAL_KIND) :: timestep      = 1.0E-5_REAL_KIND

contains
   subroutine main_loop_control
   implicit none
      current_iteration = current_iteration + 1
      solution_time = solution_time + timestep

      residual_on_screen = .false.
      solution_to_file   = .false.

      res_max = 0.0E0_REAL_KIND
      res_avg = 0.0E0_REAL_KIND

!      if (mod(current_iteration,10) == 0) then
!         solution_to_file = .true.
!      end if
      if (current_iteration >= max_iteration) then

         end_main_loop = .true.
         solution_to_file = .true.
         residual_on_screen = .true.
      end if

      if ( current_iteration <= 50 .or. &
         & mod(current_iteration,100) == 0) then
         residual_on_screen = .true.
      end if
!         residual_on_screen = .true.
   end subroutine main_loop_control
end module control_mod
