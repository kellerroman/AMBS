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
   real(REAL_KIND) :: T
   type(block_type), intent(inout) :: block

   do k = 1, block % nCells(3)
      do j = 1, block % nCells(2)
         do i = 1, block % nCells(1)
               T = block % temperatures (i,j,k)
               block % viscosities (i,j,k) = C1 * (T**1.5E0_REAL_KIND) / (T  + C2)
               block % heatKoeffs  (i,j,k) = block % viscosities(i,j,k) * Cp / Pr
            end do
         end do
      end do
   end subroutine update_thermprop
end module
