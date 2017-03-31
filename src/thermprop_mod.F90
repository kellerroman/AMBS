!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!     THERMPROP 

module thermprop_mod
   use const_mod
implicit none
   real(REAL_KIND)                     :: Pr    = 0.7E0_REAL_KIND
   real(REAL_KIND)                     :: Sc    = 0.7E0_REAL_KIND
   real(REAL_KIND),parameter           :: C1    = 1.458E-6_REAL_KIND
   real(REAL_KIND),parameter           :: C2    = 110.4E0_REAL_KIND
   real(REAL_KIND)                     :: Cp    = GAMMA*RGAS/ (GAMMA-1.0E0_REAL_KIND)
contains
   subroutine update_thermprop(block)
      use data_mod, only: block_type
   implicit none
   integer :: i,j,k
   real(REAL_KIND) :: T,rho
   type(block_type), intent(inout) :: block

   do k = 1, block % nCells(3)
      do j = 1, block % nCells(2)
         do i = 1, block % nCells(1)
               rho = block % vars (i,j,k,VEC_RHO)
               block % pressures(i,j,k)    = GM1                               &
                                           * ( block % vars(i,j,k,VEC_ENE)     &
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
               T                            = block % pressures (i,j,k) & 
                                            / ( RGas * rho )
               block % temperatures (i,j,k) = T
               block % viscosities (i,j,k)  = C1 * (T**1.5E0_REAL_KIND) / (T  + C2)
               block % heatKoeffs  (i,j,k)  = block % viscosities(i,j,k) * Cp / Pr
               block % schallges   (i,j,k) = sqrt( gamma * RGas * block % temperatures(i,j,k))
           end do
         end do
      end do
   end subroutine update_thermprop
end module
