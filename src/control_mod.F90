module control_mod
   use const_mod
implicit none

logical :: end_main_loop
logical :: residual_on_screen
logical :: solution_to_file
integer(INT_KIND) :: max_iteration
integer(INT_KIND) :: current_iteration = 0
integer(INT_KIND) :: nBoundaryCells = 2
integer(INT_KIND) :: solution_out 
integer(INT_KIND) :: residual_out
integer(INT_KIND) :: equation
integer(INT_KIND) :: turbulence
integer(INT_KIND) :: space_order
integer(INT_KIND) :: limiter
integer(INT_KIND) :: riemann_solver
integer(INT_KIND) :: timestep_method            ! Wie wird der Timestep berechnet / gewaehlt
                                                ! 1: constant aus config
                                                ! 2: global constant min berechnet
                                                ! 3: local  unterschiedlich berechnet
integer(INT_KIND) :: time_order

integer(INT_KIND) :: res_out1 = 1      ! Residuum 1 welches gemittelt und ausgegben wird
integer(INT_KIND) :: res_out2 = 2      ! Residuum 2 welches gemittelt und ausgegben wird

integer(INT_KIND),private :: start_time(8)

real(REAL_KIND) :: res_max(2)
real(REAL_KIND) :: res_avg(2)
real(REAL_KIND) :: solution_time = 0.0E+0_REAL_KIND
real(REAL_KIND) :: cfl
real(REAL_KIND) :: timestep

real(REAL_KIND) :: pressure_out = 1.0E+5_REAL_KIND

real(REAL_KIND) :: c_les_sgs
!< SUBGRID CONSTANT FOR LES LENGTHSCALE

integer :: fio
contains

   subroutine main_loop_control
   implicit none
      current_iteration = current_iteration + 1
      solution_time = solution_time + timestep

      residual_on_screen = .false.
      solution_to_file   = .false.

      res_max = 0.0E0_REAL_KIND
      res_avg = 0.0E0_REAL_KIND

      if (mod(current_iteration,solution_out) == 0) then
         solution_to_file = .true.
      end if
      if (current_iteration >= max_iteration) then

         end_main_loop = .true.
         solution_to_file = .true.
         residual_on_screen = .true.
      end if

      if ( current_iteration <= 50 .or. &
         & mod(current_iteration,residual_out) == 0) then
         residual_on_screen = .true.
      end if
      if (timestep_method == 2) then
         timestep = 1.0E20_REAL_KIND
      end if

!         residual_on_screen = .true.
   end subroutine main_loop_control

   function get_runtime() result (runtime)
   real(REAL_KIND) :: runtime
   integer(INT_KIND) :: current_time(8)

   call date_and_time(values = current_time)

   runtime =  real(current_time(7) - start_time(7)) &           !sekonds 
             +real(current_time(8) - start_time(8))/ 1.0d3 &    ! miliseconds
             +real(current_time(6) - start_time(6))* 6.0d1 &    ! minutes
             +real(current_time(5) - start_time(5))* 3.6d3      ! hours

   end function

   subroutine set_start_time()
   call date_and_time(values = start_time)
   end subroutine

end module control_mod
