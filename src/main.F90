program AMBS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PURPOSE: main program routine 
!
! AUTHOR: Roman Keller(RK)
!
! START DATE: 04.03.2016
! 
! LAST CHANGE: 12.04.2016
!
! CHANGELOG:
! 04.03.2016,RK: Start of Coding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
implicit none
   integer :: b
call screen_wr("AMBS",1)

call datin_control()
call datin_sol()
call calc_grid()
call init_sol()
do b = 1, nBlock
   call init_boundary(blocks(b),nBoundaryCells)
end do
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
   call update_residual(blocks(b),res_max,res_avg)
!   write(*,'(A3,5(1X,F10.4,1X,10X))') "rho",blocks(b) % vars(:,1,1,1)
!   write(*,'(A3,5(1X,F10.4,1X,10X))') "SpU",blocks(b) % vars(:,1,1,2)
!   write(*,'(A3,5(1X,F10.4,1X,10X))') "SpV",blocks(b) % vars(:,1,1,3)
!   write(*,'(A3,5(1X,F10.4,1X,10X))') "SpW",blocks(b) % vars(:,1,1,4)
!   write(*,'(A3,5(1X,F10.4,1X,10X))') "Ene",blocks(b) % vars(:,1,1,5)
!   write(*,'(A3,22X,3(1X,ES10.3,1X,10X))') "Rs1",blocks(b) % residuals(:,1,1,1)
!   write(*,'(A3,22X,3(1X,ES10.3,1X,10X))') "Rs2",blocks(b) % residuals(:,1,1,2)
!   write(*,'(A3,22X,3(1X,ES10.3,1X,10X))') "Rs3",blocks(b) % residuals(:,1,1,3)
!   write(*,'(A3,22X,3(1X,ES10.3,1X,10X))') "Rs4",blocks(b) % residuals(:,1,1,4)
!   write(*,'(A3,22X,3(1X,ES10.3,1X,10X))') "Rs5",blocks(b) % residuals(:,1,1,5)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "Fl1",blocks(b) % fluxesI(:,1,1,1)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "Fl2",blocks(b) % fluxesI(:,1,1,2)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "Fl3",blocks(b) % fluxesI(:,1,1,3)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "Fl4",blocks(b) % fluxesI(:,1,1,4)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "Fl5",blocks(b) % fluxesI(:,1,1,5)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "FAr",blocks(b) % CellFaceAreasI(:,1,1)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "VeX",blocks(b) % CellFaceVecsI(1,:,1,1)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "VeY",blocks(b) % CellFaceVecsI(2,:,1,1)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "VeZ",blocks(b) % CellFaceVecsI(3,:,1,1)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "rhL",blocks(b) % faceVarsLeftI(:,1,1,2)
!   write(*,'(A3,5(1X,10X,1X,ES10.3))') "rhR",blocks(b) % faceVarsRightI(:,1,1,2)
   end do block_loop

   if(residual_on_screen) call screen_residual()

   block_loop_update: do b = 1, nBlock
      call update_sol(blocks(b))
   end do block_loop_update
   if (solution_to_file) call datout_sol()
end do main_loop
call screen_wr("AMBS done",1)
end program AMBS
