module turb_mod
   use const_mod
implicit none
contains
subroutine dudn_cell_center(cellvars &
                           ,dnw,dns,dnb & 
                           ,vol         &
                           ,dudn     &
                           ,nBoundaryCells)
   implicit none
   integer, intent(in) :: nBoundaryCells
   real(REAL_KIND), intent(in)  :: cellvars(1-nBoundaryCells:,1-nBoundaryCells:,1-nBoundaryCells:,:)
   real(REAL_KIND), intent(in)  :: dnw(:,:,:,:)
   real(REAL_KIND), intent(in)  :: dns(:,:,:,:)
   real(REAL_KIND), intent(in)  :: dnb(:,:,:,:)
   real(REAL_KIND), intent(in)  :: vol(:,:,:)
   
   real(REAL_KIND), intent(out) :: dudn(:,:,:,:,:)

   integer :: i,j,k
   real(REAL_KIND) :: temp

   do k = 1, ubound(cellvars,3) - nBoundaryCells
      do j = 1, ubound(cellvars,2) - nBoundaryCells
         do i = 1, ubound(cellvars,1) - nBoundaryCells
            temp = 0.5d0 * vol(i,j,k)
            dudn(i,j,k,GRAD_SPU,GRAD_DX) =(- cellvars(i+1,j  ,k  ,VEC_SPU) * dnw(1,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPU) * dns(1,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPU) * dnb(1,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPU) * dnw(1,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPU) * dns(1,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPU) * dnb(1,i  ,j  ,k  ) &
                                          ) * temp
            dudn(i,j,k,GRAD_SPU,GRAD_DY) =(- cellvars(i+1,j  ,k  ,VEC_SPU) * dnw(2,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPU) * dns(2,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPU) * dnb(2,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPU) * dnw(2,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPU) * dns(2,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPU) * dnb(2,i  ,j  ,k  ) &
                                          ) * temp
            dudn(i,j,k,GRAD_SPU,GRAD_DZ) =(- cellvars(i+1,j  ,k  ,VEC_SPU) * dnw(3,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPU) * dns(3,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPU) * dnb(3,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPU) * dnw(3,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPU) * dns(3,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPU) * dnb(3,i  ,j  ,k  ) &
                                          ) * temp

            dudn(i,j,k,GRAD_SPV,GRAD_DX) =(- cellvars(i+1,j  ,k  ,VEC_SPV) * dnw(1,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPV) * dns(1,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPV) * dnb(1,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPV) * dnw(1,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPV) * dns(1,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPV) * dnb(1,i  ,j  ,k  ) &
                                          ) * temp
            dudn(i,j,k,GRAD_SPV,GRAD_DY) =(- cellvars(i+1,j  ,k  ,VEC_SPV) * dnw(2,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPV) * dns(2,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPV) * dnb(2,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPV) * dnw(2,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPV) * dns(2,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPV) * dnb(2,i  ,j  ,k  ) &
                                          ) * temp
            dudn(i,j,k,GRAD_SPV,GRAD_DZ) =(- cellvars(i+1,j  ,k  ,VEC_SPV) * dnw(3,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPV) * dns(3,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPV) * dnb(3,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPV) * dnw(3,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPV) * dns(3,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPV) * dnb(3,i  ,j  ,k  ) &
                                          ) * temp

            dudn(i,j,k,GRAD_SPW,GRAD_DX) =(- cellvars(i+1,j  ,k  ,VEC_SPW) * dnw(1,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPW) * dns(1,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPW) * dnb(1,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPW) * dnw(1,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPW) * dns(1,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPW) * dnb(1,i  ,j  ,k  ) &
                                          ) * temp
            dudn(i,j,k,GRAD_SPW,GRAD_DY) =(- cellvars(i+1,j  ,k  ,VEC_SPW) * dnw(2,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPW) * dns(2,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPW) * dnb(2,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPW) * dnw(2,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPW) * dns(2,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPW) * dnb(2,i  ,j  ,k  ) &
                                          ) * temp
            dudn(i,j,k,GRAD_SPW,GRAD_DZ) =(- cellvars(i+1,j  ,k  ,VEC_SPW) * dnw(3,i+1,j  ,k  ) &
                                           - cellvars(i  ,j+1,k  ,VEC_SPW) * dns(3,i  ,j+1,k  ) &
                                           - cellvars(i  ,j  ,k+1,VEC_SPW) * dnb(3,i  ,j  ,k+1) &
                                           + cellvars(i-1,j  ,k  ,VEC_SPW) * dnw(3,i  ,j  ,k  ) &
                                           + cellvars(i  ,j-1,k  ,VEC_SPW) * dns(3,i  ,j  ,k  ) &
                                           + cellvars(i  ,j  ,k-1,VEC_SPW) * dnb(3,i  ,j  ,k  ) &
                                          ) * temp

!            write(*,*)  " == X == "
!            write(*,*)                       cellvars(i+1,j  ,k  ,VEC_SPU) , "dnw1", dnw(1,i+1,j  ,k  ) 
!            write(*,*)                       cellvars(i  ,j+1,k  ,VEC_SPU) , "dns1", dns(1,i  ,j+1,k  )
!            write(*,*)                       cellvars(i  ,j  ,k+1,VEC_SPU) , "dnb1", dnb(1,i  ,j  ,k+1)
!            write(*,*)                       cellvars(i-1,j  ,k  ,VEC_SPU) , "dnw1", dnw(1,i  ,j  ,k  )
!            write(*,*)                       cellvars(i  ,j-1,k  ,VEC_SPU) , "dns1", dns(1,i  ,j  ,k  )
!            write(*,*)                       cellvars(i  ,j  ,k-1,VEC_SPU) , "dnb1", dnb(1,i  ,j  ,k  )
!            write(*,*)  " == Y == "
!            write(*,*)                       cellvars(i+1,j  ,k  ,VEC_SPV) , "dnw2", dnw(2,i+1,j  ,k  ) 
!            write(*,*)                       cellvars(i  ,j+1,k  ,VEC_SPV) , "dns2", dns(2,i  ,j+1,k  )
!            write(*,*)                       cellvars(i  ,j  ,k+1,VEC_SPV) , "dnb2", dnb(2,i  ,j  ,k+1)
!            write(*,*)                       cellvars(i-1,j  ,k  ,VEC_SPV) , "dnw2", dnw(2,i  ,j  ,k  )
!            write(*,*)                       cellvars(i  ,j-1,k  ,VEC_SPV) , "dns2", dns(2,i  ,j  ,k  )
!            write(*,*)                       cellvars(i  ,j  ,k-1,VEC_SPV) , "dnb2", dnb(2,i  ,j  ,k  )
!            write(*,*)  " == Z == "
!            write(*,*)                       cellvars(i+1,j  ,k  ,VEC_SPW) , "dnw3", dnw(3,i+1,j  ,k  ) 
!            write(*,*)                       cellvars(i  ,j+1,k  ,VEC_SPW) , "dns3", dns(3,i  ,j+1,k  )
!            write(*,*)                       cellvars(i  ,j  ,k+1,VEC_SPW) , "dnb3", dnb(3,i  ,j  ,k+1)
!            write(*,*)                       cellvars(i-1,j  ,k  ,VEC_SPW) , "dnw3", dnw(3,i  ,j  ,k  )
!            write(*,*)                       cellvars(i  ,j-1,k  ,VEC_SPW) , "dns3", dns(3,i  ,j  ,k  )
!            write(*,*)                       cellvars(i  ,j  ,k-1,VEC_SPW) , "dnb3", dnb(3,i  ,j  ,k  )
!
!            write(*,*) dudn(i,j,k,:,1)
!            write(*,*) dudn(i,j,k,:,2)
!            write(*,*) dudn(i,j,k,:,3)
!
!            stop

         end do
      end do
   end do 
end subroutine dudn_cell_center

subroutine smagorinsky_SGS (viskosity        &
                           ,dudn             &
                           ,cellvars         &
                           ,lles             &
                           ,nBoundaryCells   )
   implicit none
   integer, intent(in) :: nBoundaryCells
   real(REAL_KIND), intent(out) :: viskosity (1-nBoundaryCells:,1-nBoundaryCells:,1-nBoundaryCells:)
   real(REAL_KIND), intent(in)  :: cellvars  (1-nBoundaryCells:,1-nBoundaryCells:,1-nBoundaryCells:,:)
   real(REAL_KIND), intent(in)  :: lles      (:,:,:)
   real(REAL_KIND), intent(in)  :: dudn      (:,:,:,:,:)

   integer :: i,j,k!,dir1,dir2
   real(REAL_KIND) :: temp,t12,t13,t23

   do k = 1, ubound(cellvars,3) - nBoundaryCells
      do j = 1, ubound(cellvars,2) - nBoundaryCells
         do i = 1, ubound(cellvars,1) - nBoundaryCells
            temp = 0.0d0
!            do dir2 = 1,3
!               do dir1 = 1,3
!                  temp = temp + dudn(i,j,k,dir1,dir2) * dudn(i,j,k,dir1,dir2)
!               end do
!            end do
!            temp = 2.0D0 * temp
            t12 = dudn(i,j,k,1,2) + dudn(i,j,k,2,1)
            t13 = dudn(i,j,k,1,3) + dudn(i,j,k,3,1)
            t23 = dudn(i,j,k,2,3) + dudn(i,j,k,3,2)
            temp = 2.0d0 * ( &
                             dudn(i,j,k,1,1) * dudn(i,j,k,1,1) &
                           + dudn(i,j,k,2,2) * dudn(i,j,k,2,2) &
                           + dudn(i,j,k,3,3) * dudn(i,j,k,3,3) &
                           ) &
                 + t12 * t12 + t13 * t13 + t23 * t23

!            write(*,*) viskosity(i,j,k) , cellvars(i,j,k,VEC_RHO)  &
!                                                * lles (i,j,k)             &
!                                                * sqrt(temp)
!            write(*,*) dudn(i,j,k,:,:)
!            stop

            viskosity(i,j,k) = viskosity(i,j,k) + cellvars(i,j,k,VEC_RHO)  &
                                                * lles (i,j,k)             &
                                                * sqrt(temp)
         end do
      end do
   end do
end subroutine smagorinsky_SGS
end module
