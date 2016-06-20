module face_values_mod
   use const_mod
   use control_mod, only: space_order

implicit none
contains
   subroutine face_values( cellvars &
                         , faceVarsLeftI &
                         , faceVarsRightI &
                         , faceVarsLeftJ &
                         , faceVarsRightJ &
                         , faceVarsLeftK &
                         , faceVarsRightK &
                         , nBoundaryCells)
   implicit none
   integer, intent(in) :: nBoundaryCells
   real(REAL_KIND), intent(in)  :: cellvars(1-nBoundaryCells:,1-nBoundaryCells:,1-nBoundaryCells:,:)
   real(REAL_KIND), intent(out) :: faceVarsLeftI(:,:,:,:)
   real(REAL_KIND), intent(out) :: faceVarsRightI(:,:,:,:)
   real(REAL_KIND), intent(out) :: faceVarsLeftJ(:,:,:,:)
   real(REAL_KIND), intent(out) :: faceVarsRightJ(:,:,:,:)
   real(REAL_KIND), intent(out) :: faceVarsLeftK(:,:,:,:)
   real(REAL_KIND), intent(out) :: faceVarsRightK(:,:,:,:)

   integer :: i,j,k
   if (space_order == 1) then
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells + 1
              faceVarsLeftI  (i,j,k,:) =  cellvars (i-1,j  ,k  ,:)
              faceVarsRightI (i,j,k,:) =  cellvars (i  ,j  ,k  ,:)
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells + 1
            do i = 1, ubound(cellvars,1) - nBoundaryCells 
              faceVarsLeftJ  (i,j,k,:) =  cellvars (i  ,j-1,k  ,:)
              faceVarsRightJ (i,j,k,:) =  cellvars (i  ,j  ,k  ,:)
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells + 1
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells
              faceVarsLeftK  (i,j,k,:) =  cellvars (i  ,j  ,k-1,:)
              faceVarsRightK (i,j,k,:) =  cellvars (i  ,j  ,k  ,:)
            end do
         end do
      end do
   else if (space_order == 2) then
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells + 1
              faceVarsLeftI (i,j,k,:) = ( cellvars(i  ,j,k,:) - cellvars(i-1,j,k,:) )  &
                                      / ( cellvars(i-1,j,k,:) - cellvars(i-2,j,k,:) + EPSI)
              faceVarsRightI(i,j,k,:) = ( cellvars(i  ,j,k,:) - cellvars(i-1,j,k,:) )  &
                                      / ( cellvars(i+1,j,k,:) - cellvars(i  ,j,k,:) + EPSI)
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells + 1
              call minmod(faceVarsLeftI (i,j,k,:))
              call minmod(faceVarsRightI(i,j,k,:))
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells + 1
              faceVarsLeftI (i,j,k,:) =   cellvars(i-1,j,k,:) - faceVarsLeftI (i,j,k,:) * HALF  &
                                      * ( cellvars(i-1,j,k,:) - cellvars(i-2,j,k,:) )
              faceVarsRightI(i,j,k,:) =   cellvars(i  ,j,k,:) - faceVarsRightI(i,j,k,:) * HALF  &
                                      * ( cellvars(i+1,j,k,:) - cellvars(i  ,j,k,:) )
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells + 1
            do i = 1, ubound(cellvars,1) - nBoundaryCells 
              faceVarsLeftJ (i,j,k,:) = ( cellvars(i,j  ,k,:) - cellvars(i,j-1,k,:) )  &
                                      / ( cellvars(i,j-1,k,:) - cellvars(i,j-2,k,:) + EPSI)
              faceVarsRightJ(i,j,k,:) = ( cellvars(i,j  ,k,:) - cellvars(i,j-1,k,:) )  &
                                      / ( cellvars(i,j+1,k,:) - cellvars(i,j  ,k,:) + EPSI)
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells + 1
            do i = 1, ubound(cellvars,1) - nBoundaryCells
              call minmod(faceVarsLeftJ (i,j,k,:))
              call minmod(faceVarsRightJ(i,j,k,:))
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells
         do j = 1, ubound(cellvars,2) - nBoundaryCells + 1
            do i = 1, ubound(cellvars,1) - nBoundaryCells
              faceVarsLeftJ (i,j,k,:) =   cellvars(i,j-1,k,:) - faceVarsLeftJ (i,j,k,:) * HALF  &
                                      * ( cellvars(i,j-1,k,:) - cellvars(i,j-2,k,:) )
              faceVarsRightJ(i,j,k,:) =   cellvars(i,j  ,k,:) - faceVarsRightJ(i,j,k,:) * HALF  &
                                      * ( cellvars(i,j+1,k,:) - cellvars(i,j  ,k,:) )
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells + 1
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells 
              faceVarsLeftK (i,j,k,:) = ( cellvars(i,j,k  ,:) - cellvars(i,j,k-1,:) )  &
                                      / ( cellvars(i,j,k-1,:) - cellvars(i,j,k-2,:) + EPSI)
              faceVarsRightK(i,j,k,:) = ( cellvars(i,j,k  ,:) - cellvars(i,j,k-1,:) )  &
                                      / ( cellvars(i,j,k+1,:) - cellvars(i,j,k  ,:) + EPSI)
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells + 1
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells
              call minmod(faceVarsLeftK (i,j,k,:))
              call minmod(faceVarsRightK(i,j,k,:))
            end do
         end do
      end do
      do k = 1, ubound(cellvars,3) - nBoundaryCells + 1 
         do j = 1, ubound(cellvars,2) - nBoundaryCells
            do i = 1, ubound(cellvars,1) - nBoundaryCells
              faceVarsLeftK (i,j,k,:) =   cellvars(i,j,k-1,:) - faceVarsLeftK (i,j,k,:) * HALF  &
                                      * ( cellvars(i,j,k-1,:) - cellvars(i,j,k-2,:) )
              faceVarsRightK(i,j,k,:) =   cellvars(i,j,k  ,:) - faceVarsRightK(i,j,k,:) * HALF  &
                                      * ( cellvars(i,j,k+1,:) - cellvars(i,j,k  ,:) )
            end do
         end do
      end do
   end if
   end subroutine face_values
   subroutine minmod ( a )
      implicit none
      real(REAL_KIND), intent(inout) ::a(:)
      integer :: n
      do n = 1, ubound(a,1)
            a(n) = max(0.0E0_REAL_KIND,min(1.0E0_REAL_KIND,a(n)))
      end do
   end subroutine minmod
end module face_values_mod
