program AMBS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: main program routine 
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 04.03.2016
! 
! LAST CHANGE: 08.09.2016
!
! CHANGELOG:
! 04.03.2016,RK: Start of Coding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use boundary_mod, only: init_boundary, set_boundary
   use const_mod
   use control_mod
   use data_mod
   use face_values_mod, only: face_values
   use screen_io_mod
   use time_disc_mod, only: update_residual & 
                          , update_sol
   use file_io_mod
   use inv_fluxes_mod, only: calc_fluxes
   use viscous_fluxes_mod, only: calc_face_gradients,calc_viscous_fluxes
   use thermprop_mod, only: update_thermprop
   use turb_mod, only: dudn_cell_center, smagorinsky_SGS
implicit none
   integer :: b
call screen_wr_start()

call datin_control()
call datin_sol()
call datin_bc()
call calc_grid()
call init_sol()

do b = 1, nBlock
   call update_thermprop(blocks(b) ) 
   call init_boundary(blocks(b),nBoundaryCells)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!                     START OF CALCULATION                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call set_start_time()
call screen_wr_start_calc()
end_main_loop = .false.

current_iteration = 0
main_loop: do while (.not. end_main_loop)
   call main_loop_control()
   block_loop: do b = 1, nBlock
      call set_boundary(blocks(b),nBoundaryCells)
      call face_values ( blocks(b) % vars     &
                       , blocks(b) % faceVarsLeftI & 
                       , blocks(b) % faceVarsRightI & 
                       , blocks(b) % faceVarsLeftJ & 
                       , blocks(b) % faceVarsRightJ & 
                       , blocks(b) % faceVarsLeftK & 
                       , blocks(b) % faceVarsRightK & 
                       , nBoundaryCells)
      call calc_fluxes ( blocks(b) )
      call calc_face_gradients(blocks(b))
      call dudn_cell_center ( blocks(b) % vars              &
                            , blocks(b) % dnw               &
                            , blocks(b) % dns               &
                            , blocks(b) % dnb               &
                            , blocks(b) % volumes           &
                            , blocks(b) % dUdN              &
                            , nBoundaryCells                )
!      call smagorinsky_SGS( blocks(b) % viscosities         &
!                          , blocks(b) % dUdN                &
!                          , blocks(b) % vars                &
!                          , blocks(b) % lles                &
!                          , nBoundaryCells                  ) 
      if (equation == EQU_TYP_NS) then
         call calc_viscous_fluxes(blocks(b)) 
      end if
      call update_residual(blocks(b) % residuals,blocks(b),res_max,res_avg)
   end do block_loop

   if(residual_on_screen) call screen_residual()

   block_loop_update: do b = 1, nBlock
      call update_sol(blocks(b))
      call update_thermprop(blocks(b) ) 
   end do block_loop_update
   if (solution_to_file) call datout_sol()
end do main_loop
call screen_wr_end()
end program AMBS
