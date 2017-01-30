module time_disc_mod   
   use const_mod
   use data_mod, only: block_type
implicit none
contains
   subroutine update_residual(residuals,block,res_max,res_avg)
   implicit none
   real(REAL_KIND), intent(out) :: residuals(:,:,:,:)
   type(block_type), intent(inout) :: block
   real(REAL_KIND) , intent(inout) :: res_max
   real(REAL_KIND) , intent(inout) :: res_avg
   integer :: i,j,k
   real(REAL_KIND) :: re
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
            re = abs(residuals(i,j,k,1))
            res_avg = res_avg + re
            res_max = max(res_max,re)
!            if ( i == 50 .and. j == 1) then
!            write(*,*) i,j,k   
!            write(*,*) "RHO",block % vars (i:i+1,j,k,1)
!            write(*,*) "SPU",block % vars (i:i+1,j,k,2)
!            write(*,*) "i  ",block % fluxesI        (i  ,j  ,k  ,1:3) 
!            write(*,*) "i+1",block % fluxesI        (i+1,j  ,k  ,1:3) 
!            write(*,*) "i+2",block % fluxesI        (i+2,j  ,k  ,1:3) 
!            write(*,*) "j  ",block % fluxesJ        (i  ,j  ,k  ,1:3)  
!            write(*,*) "j+1",block % fluxesJ        (i  ,j+1,k  ,1:3)  
!            write(*,*) "Res",block % residuals(i,j,k,1:3)
!            !stop
!         end if
         end do
      end do
   end do
   end subroutine update_residual
   subroutine update_sol(block)
      use control_mod, only: timestep
   implicit none
   type(block_type), intent(inout) :: block

   integer :: i,j,k

   real(REAL_KIND) :: rho
   do k = 1, ubound(block % residuals,3)
      do j = 1, ubound(block % residuals,2)
         do i = 1, ubound(block % residuals,1)
            block % cons_vars(i,j,k,:) = block % cons_vars(i,j,k,:) &
                                       + block % residuals (i,j,k,:) &
                                       * block % volumes(i,j,k) & 
                                       * timestep
            rho = block % cons_vars(i,j,k,VEC_RHO)
            block % vars(i,j,k,VEC_RHO) = rho
            block % vars(i,j,k,VEC_SPU:VEC_SPW) = block % cons_vars(i,j,k,VEC_SPU:VEC_SPW) / rho
            block % vars(i,j,k,VEC_ENE) = block % cons_vars(i,j,k,VEC_ENE)
            block % pressures(i,j,k) = (GAMMA - 1.0E0_REAL_KIND)         &
                                     * (                                 &
                                         block % vars(i,j,k,VEC_ENE)     &
                                       - 0.5E0_REAL_KIND                 &
                                       * rho                             &
                                       * ( block % vars (i,j,k,VEC_SPU)  &
                                         * block % vars (i,j,k,VEC_SPU)  & 
                                         + block % vars (i,j,k,VEC_SPV)  &
                                         * block % vars (i,j,k,VEC_SPV)  &
                                         + block % vars (i,j,k,VEC_SPW)  &
                                         * block % vars (i,j,k,VEC_SPW)  &
                                         )                               &
                                       )
            block % temperatures (i,j,k) = block % pressures (i,j,k) & 
                                         / ( RGas * rho )
         end do
      end do
   end do
   end subroutine update_sol
end module time_disc_mod
