   use const_mod
   use data_mod, only: block_type
implicit none
contains
   subroutine update_residual(block,res_max,res_avg)
   implicit none
   type(block_type), intent(inout) :: block
   real(REAL_KIND) , intent(inout) :: res_max
   real(REAL_KIND) , intent(inout) :: res_avg
   integer :: i,j,k
   real(REAL_KIND) :: res
   do k = 1, ubound(block % residuals,3)
      do j = 1, ubound(block % residuals,2)
         do i = 1, ubound(block % residuals,1)
            block % residuals(i,j,k,:) = &
               +  block % CellFaceAreasI (i  ,j  ,k  ) &
               *  block % fluxesI        (i  ,j  ,k  ,:) & 
               -  block % CellFaceAreasI (i+1,j  ,k  ) &
               *  block % fluxesI        (i+1,j  ,k  ,:) & 
               +  block % CellFaceAreasJ (i  ,j  ,k  ) &
               *  block % fluxesJ        (i  ,j  ,k  ,:) & 
               -  block % CellFaceAreasJ (i  ,j+1,k  ) &
               *  block % fluxesJ        (i  ,j+1,k  ,:) & 
               +  block % CellFaceAreasK (i  ,j  ,k  ) &
               *  block % fluxesK        (i  ,j  ,k  ,:) & 
               -  block % CellFaceAreasK (i  ,j  ,k+1) &
               *  block % fluxesK        (i  ,j  ,k+1,:)
            res = abs(block % residuals(i,j,k,1))
            res_avg = res_avg + res
            res_max = max(res_max,res)
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
            rho = block % cons_vars(i,j,k,1)
            block % vars(i,j,k,1) = rho
            block % vars(i,j,k,2:4) = block % cons_vars(i,j,k,2:4) / rho
            block % vars(i,j,k,5) = block % cons_vars(i,j,k,5)
         end do
      end do
   end do
   end subroutine update_sol
end module time_disc_mod
