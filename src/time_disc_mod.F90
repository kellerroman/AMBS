module time_disc_mod   
   use const_mod
   use data_mod, only: block_type
implicit none
contains
   subroutine update_residual(residuals,block,res_max,res_avg)
      use control_mod ,only: res_out1, res_out2
   implicit none
   real(REAL_KIND), intent(out) :: residuals(:,:,:,:)
   type(block_type), intent(inout) :: block
   real(REAL_KIND) , intent(inout) :: res_max(2)
   real(REAL_KIND) , intent(inout) :: res_avg(2)
   integer :: i,j,k
   real(REAL_KIND) :: re(2)
   do k = 1, ubound(residuals,3)
      do j = 1, ubound(residuals,2)
         do i = 1, ubound(residuals,1)
            residuals(i,j,k,:) = &
               +  block % CellFaceAreasI (i  ,j  ,k  ) &
               *  ( block % fluxesI      (i  ,j  ,k  ,:) & 
                  + block % visFluxesI   (i  ,j  ,k  ,:))&  
               -  block % CellFaceAreasI (i+1,j  ,k  ) &
               *  ( block % fluxesI      (i+1,j  ,k  ,:) & 
                  + block % visFluxesI   (i+1,j  ,k  ,:))& 
               +  block % CellFaceAreasJ (i  ,j  ,k  ) &
               *  ( block % fluxesJ      (i  ,j  ,k  ,:) & 
                  + block % visFluxesJ   (i  ,j  ,k  ,:))& 
               -  block % CellFaceAreasJ (i  ,j+1,k  ) &
               *  ( block % fluxesJ      (i  ,j+1,k  ,:) & 
                  + block % visFluxesJ   (i  ,j+1,k  ,:))& 
               +  block % CellFaceAreasK (i  ,j  ,k  ) &
               *  ( block % fluxesK      (i  ,j  ,k  ,:) & 
                  + block % visFluxesK   (i  ,j  ,k  ,:))& 
               -  block % CellFaceAreasK (i  ,j  ,k+1) &
               *  ( block % fluxesK      (i  ,j  ,k+1,:)&
                  + block % visFluxesK   (i  ,j  ,k+1,:))
            re(1) = abs(residuals(i,j,k,res_out1))
            re(2) = abs(residuals(i,j,k,res_out2))
            res_avg = res_avg + re
            res_max(1) = max(res_max(1),re(1))
            res_max(2) = max(res_max(2),re(2))
            if (ubound(residuals,1) >= 100) then
               if (i == 10 .and. j == 1)  write(6661,*) block % vars(i,j:j+3,k,VEC_SPU)
               if (i == 65 .and. j == 1)  write(6662,*) block % vars(i,j:j+3,k,VEC_SPU)
               if (i == 120 .and. j == 1)  write(6663,*) block % vars(i,j:j+3,k,VEC_SPU)
               if (i == 128 .and. j == 1)  write(6664,*) block % vars(i,j:j+3,k,VEC_SPU)
            end if
!            if ( i == 10 .and. j <= 2) then
!               write(*,*) i,j,k   
!               write(*,*) "RHO",block % vars (i,j-1:j+1,k,1)
!               write(*,*) "SPU",block % vars (i,j-1:j+1,k,2)
!               write(*,*) "i  ",block % fluxesI        (i  ,j  ,k  ,1:3) 
!               write(*,*) "i+1",block % fluxesI        (i+1,j  ,k  ,1:3) 
!               write(*,*) "i+2",block % fluxesI        (i+2,j  ,k  ,1:3) 
!               write(*,*) "j  ",block % fluxesJ        (i  ,j  ,k  ,1:3)  
!               write(*,*) "j+1",block % fluxesJ        (i  ,j+1,k  ,1:3)  
!               write(*,*) "Res",block % residuals(i,j,k,1:3)
!               write(*,*) block % visFluxesI   (i  ,j  ,k  ,:)  
!               write(*,*) block % visFluxesI   (i+1,j  ,k  ,:) 
!               write(*,*) block % visFluxesJ   (i  ,j  ,k  ,:) 
!               write(*,*) block % visFluxesJ   (i  ,j+1,k  ,:) 
!               write(*,*) block % visFluxesK   (i  ,j  ,k  ,:) 
!               write(*,*) block % visFluxesK   (i  ,j  ,k+1,:)
!               write(*,*) block % CellFaceAreasJ (i,j,k)
!               !stop
!            end if
         end do
      end do
   end do
   end subroutine update_residual
   subroutine update_sol(block)
      use control_mod, only: timestep, timestep_method
   implicit none
   type(block_type), intent(inout) :: block

   integer :: i,j,k
   real(REAL_KIND) :: tt

   real(REAL_KIND) :: rho
   do k = 1, ubound(block % residuals,3)
      do j = 1, ubound(block % residuals,2)
         do i = 1, ubound(block % residuals,1)
            if (timestep_method == 3) then
               tt = block % timesteps(i,j,k) * 0.9
            else
               tt = timestep
            end if
            block % cons_vars(i,j,k,:) = block % cons_vars (i,j,k,:) &
                                       + block % residuals (i,j,k,:) &
                                       * block % volumes   (i,j,k  ) & 
                                       * tt
                                       !!! VERSUCH EINER LOCALEN ZEITSCHRITT STEUERUNG
                                       !* min((block % volumes   (i,1,k  ) & 
                                       !/ block % volumes   (i,j,k  ) )**(1.0d0/4.0D0) & 
                                       !, 2.0D0)
            rho = block % cons_vars(i,j,k,VEC_RHO)
            block % vars(i,j,k,VEC_RHO) = rho
            block % vars(i,j,k,VEC_SPU:VEC_SPW) = block % cons_vars(i,j,k,VEC_SPU:VEC_SPW) / rho
            block % vars(i,j,k,VEC_ENE) = block % cons_vars(i,j,k,VEC_ENE)
         end do
      end do
   end do
   end subroutine update_sol
   subroutine calc_timestep(block)
      use control_mod, only: timestep, timestep_method
   implicit none
   type(block_type), intent(inout) :: block
   !real(REAL_KIND) :: inv_rho
   !real(REAL_KIND) :: f1,f2,fac
   integer :: i,j,k
   real(REAL_KIND) :: re_x, re_y,re_z,re

   do k = 1, ubound(block % timesteps,3)
      do j = 1, ubound(block % timesteps,2)
         do i = 1, ubound(block % timesteps,1)
         block % timesteps(i,j,k) = 1.0D0 / ( &
                                    abs(block % vars(i,j,k,VEC_SPU)) / block % cellsizes(i,j,k,DIR_EASTWEST)     &
                                  + abs(block % vars(i,j,k,VEC_SPV)) / block % cellsizes(i,j,k,DIR_SOUTHNORTH)   &
                                  + abs(block % vars(i,j,k,VEC_SPW)) / block % cellsizes(i,j,k,DIR_FRONTBACK)    &
                                  + block % schallges(i,j,k) * sqrt(                                        &
                                                      1.0D0 / block % cellsizes(i,j,k,DIR_EASTWEST  )**2    &
                                                    + 1.0D0 / block % cellsizes(i,j,k,DIR_SOUTHNORTH)**2    &
                                                    + 1.0D0 / block % cellsizes(i,j,k,DIR_FRONTBACK )**2   ))
         re_x = abs(block % vars(i,j,k,VEC_SPU)) * block % cellsizes(i,j,k,DIR_EASTWEST)
         re_y = abs(block % vars(i,j,k,VEC_SPV)) * block % cellsizes(i,j,k,DIR_SOUTHNORTH) 
         re_z = abs(block % vars(i,j,k,VEC_SPW)) * block % cellsizes(i,j,k,DIR_FRONTBACK) 
         re = max(min (re_x,re_y,re_z) &
                  * block % vars(i,j,k,VEC_RHO) / block % viscosities(i,j,k) &
                 ,1.0E0_REAL_KIND)
         block % timesteps(i,j,k) = block % timesteps(i,j,k) / (1.0E0_REAL_KIND + 2.0E0_REAL_KIND / re)
         end do
      end do
   end do

   if (timestep_method == 2) then
      do k = 1, ubound(block % timesteps,3)
         do j = 1, ubound(block % timesteps,2)
            do i = 1, ubound(block % timesteps,1)
               timestep = min(timestep, block % timesteps(i,j,k))
            end do
         end do
      end do
   end if
   end subroutine calc_timestep
end module time_disc_mod
