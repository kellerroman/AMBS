program AMBS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use boundary_mod, only: init_boundary, set_boundary
   use const_mod
   use control_mod
   use data_mod
   use face_values_mod, only: face_values
   use screen_io_mod
   use time_disc_mod, only: update_residual & 
                          , update_sol      &
                          , calc_timestep
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
   write(*,*) b, blocks(b) % vars(0,1,1,2)
end do

open(newunit = fio, file="residual.dat")
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
      call set_boundary(blocks,b,nBoundaryCells)
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

!      call dudn_cell_center ( blocks(b) % vars              &
!                            , blocks(b) % dnw               &
!                            , blocks(b) % dns               &
!                            , blocks(b) % dnb               &
!                            , blocks(b) % volumes           &
!                            , blocks(b) % dUdN              &
!                            , nBoundaryCells                )

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

   if (residual_on_screen) call screen_residual()

   block_loop_update: do b = 1, nBlock
      call calc_timestep(blocks(b))
      call update_sol(blocks(b))
      call update_thermprop(blocks(b) ) 
   end do block_loop_update
   if (solution_to_file) call datout_sol()
!   write(*,*) " vfi  ", blocks(1) % CellFaceAreasI(150,1,1)*blocks(1) % visFluxesI(150,1,1,2:3) &
!                       ,blocks(1) % CellFaceAreasI(150,1,1)*blocks(1) % visFluxesI(150,1,1,5)
!   write(*,*) " vfi+1", blocks(1) % CellFaceAreasI(151,1,1)*blocks(1) % visFluxesI(151,1,1,2:3) &
!                       ,blocks(1) % CellFaceAreasI(151,1,1)*blocks(1) % visFluxesI(151,1,1,5)
!   write(*,*) " vfj  ", blocks(1) % CellFaceAreasJ(150,1,1)*blocks(1) % visFluxesJ(150,1,1,2:3) &
!                       ,blocks(1) % CellFaceAreasJ(150,1,1)*blocks(1) % visFluxesJ(150,1,1,5)
!   write(*,*) " vfj+1", blocks(1) % CellFaceAreasJ(150,2,1)*blocks(1) % visFluxesJ(150,2,1,2:3) &
!                       ,blocks(1) % CellFaceAreasJ(150,2,1)*blocks(1) % visFluxesJ(150,2,1,5)
!   write(*,*) "  fi  ", blocks(1) % CellFaceAreasI(150,1,1)*blocks(1) %    FluxesI(150,1,1,2:3) &
!                       ,blocks(1) % CellFaceAreasI(150,1,1)*blocks(1) %    FluxesI(150,1,1,5)
!   write(*,*) "  fi+1", blocks(1) % CellFaceAreasI(151,1,1)*blocks(1) %    FluxesI(151,1,1,2:3) &
!                       ,blocks(1) % CellFaceAreasI(151,1,1)*blocks(1) %    FluxesI(151,1,1,5)
!   write(*,*) "  fj  ", blocks(1) % CellFaceAreasJ(150,1,1)*blocks(1) %    Fluxesj(150,1,1,2:3) &
!                       ,blocks(1) % CellFaceAreasJ(150,1,1)*blocks(1) %    Fluxesj(150,1,1,5)
!   write(*,*) "  fj+1", blocks(1) % CellFaceAreasJ(150,2,1)*blocks(1) %    Fluxesj(150,2,1,2:3) &
!                       ,blocks(1) % CellFaceAreasJ(150,2,1)*blocks(1) %    Fluxesj(150,2,1,5)
!   write(*,*) "  U   ",blocks(1) % vars (150,1:9,1,VEC_SPU)
!   write(*,*) "  dudy",blocks(1) % dudnJ(150,1:9,1,GRAD_SPU,GRAD_DY)
!   write(*,*) "  res ",blocks(1) % Residuals(150,1:9,1,2) * blocks(1) % volumes(150,1:9,1) * timestep
!   write(*,*) "  dt  ",blocks(1) % timesteps(150,1:9,1) 
!   write(*,*) "  vol ",(blocks(1) % volumes(150,1,1)/blocks(1) % volumes(150,1:9,1))**(1.0D0/4.0D0)
!  b = 0
!  write(*,*) blocks(b) % vars(0,1,1,2), blocks(b) % vars(0,1,1,5), blocks(b) % pressures(0,1,1), blocks(b) % temperatures(0,1,1)
!  write(*,*) blocks(b) % vars(1,1,1,2), blocks(b) % vars(1,1,1,5), blocks(b) % pressures(1,1,1), blocks(b) % temperatures(1,1,1)
   !write(*,*) minval(blocks(1) % timesteps),maxval(blocks(1) % timesteps)
end do main_loop
!b = 1
!write(*,*) b, blocks(b) % vars(blocks(b) % nCells(1)+1,1,1,:)
!write(*,*) b, blocks(b) % vars(blocks(b) % nCells(1)  ,1,1,:)
!b = 2
!write(*,*) b, blocks(b) % vars(                      0,1,1,:)
!write(*,*) b, blocks(b) % vars(                      1,1,1,:)
close(fio)
call screen_wr_end()
end program AMBS
