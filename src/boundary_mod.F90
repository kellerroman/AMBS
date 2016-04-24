module boundary_mod
   use data_mod, only: block_type

implicit none
contains
   subroutine init_boundary(block,nBoundaryCells)
   implicit none
   type(block_type), intent(inout) :: block
   integer, intent(in) :: nBoundaryCells
   
   integer :: i

   do i = 0, - nBoundaryCells + 1, -1
   block % vars(i,:,:,:) = block % vars(1,:,:,:)
   block % vars(block % nPkts(1)-i,:,:,:) = block % vars(block % nCells(1),:,:,:)   
   block % vars(:,i,:,:) = block % vars(:,1,:,:)
   block % vars(:,block % nPkts(2)-i,:,:) = block % vars(:,block % nCells(2),:,:)   
   block % vars(:,:,i,:) = block % vars(:,:,1,:)
   block % vars(:,:,block % nPkts(3)-i,:) = block % vars(:,:,block % nCells(3),:)   
   end do
   end subroutine init_boundary
   subroutine set_boundary(block)
   implicit none
   type(block_type), intent(inout) :: block
   
   block % vars(0,:,:,:) = block % vars(1,:,:,:)
   block % vars(block % nPkts(1),:,:,:) = block % vars(block % nCells(1),:,:,:)   
   block % vars(:,0,:,:) = block % vars(:,1,:,:)
   block % vars(:,block % nPkts(2),:,:) = block % vars(:,block % nCells(2),:,:)   
   block % vars(:,:,0,:) = block % vars(:,:,1,:)
   block % vars(:,:,block % nPkts(3),:) = block % vars(:,:,block % nCells(3),:)   
   end subroutine set_boundary
end module boundary_mod
