module boundary_mod
   use const_mod
   use data_mod, only: block_type

implicit none
contains
   subroutine init_boundary(block,nBoundaryCells)
   implicit none
      type(block_type), intent(inout) :: block
      integer, intent(in) :: nBoundaryCells
      
      integer :: i
   
      do i = 0, - nBoundaryCells + 1, -1
         block % vars         (i,:,:,:)                  = block % vars         (1,:,:,:)
         block % vars         (block % nPkts(1)-i,:,:,:) = block % vars         (block % nCells(1),:,:,:)   
         block % vars         (:,i,:,:)                  = block % vars         (:,1,:,:)
         block % vars         (:,block % nPkts(2)-i,:,:) = block % vars         (:,block % nCells(2),:,:)   
         block % vars         (:,:,i,:)                  = block % vars         (:,:,1,:)
         block % vars         (:,:,block % nPkts(3)-i,:) = block % vars         (:,:,block % nCells(3),:)   
         block % temperatures (i,:,:)                    = block % temperatures (1,:,:)
         block % temperatures (block % nPkts(1)-i,:,:)   = block % temperatures (block % nCells(1),:,:)
         block % temperatures (:,i,:)                    = block % temperatures (:,1,:)
         block % temperatures (:,block % nPkts(2)-i,:)   = block % temperatures (:,block % nCells(2),:)
         block % temperatures (:,:,i)                    = block % temperatures (:,:,1)
         block % temperatures (:,:,block % nPkts(3)-i)   = block % temperatures (:,:,block % nCells(3))
         block % viscosities  (i,:,:)                    = block % viscosities  (1,:,:)
         block % viscosities  (block % nPkts(1)-i,:,:)   = block % viscosities  (block % nCells(1),:,:)
         block % viscosities  (:,i,:)                    = block % viscosities  (:,1,:)
         block % viscosities  (:,block % nPkts(2)-i,:)   = block % viscosities  (:,block % nCells(2),:)
         block % viscosities  (:,:,i)                    = block % viscosities  (:,:,1)
         block % viscosities  (:,:,block % nPkts(3)-i)   = block % viscosities  (:,:,block % nCells(3))
         block % heatKoeffs   (i,:,:)                    = block % heatKoeffs   (1,:,:)
         block % heatKoeffs   (block % nPkts(1)-i,:,:)   = block % heatKoeffs   (block % nCells(1),:,:)
         block % heatKoeffs   (:,i,:)                    = block % heatKoeffs   (:,1,:)
         block % heatKoeffs   (:,block % nPkts(2)-i,:)   = block % heatKoeffs   (:,block % nCells(2),:)
         block % heatKoeffs   (:,:,i)                    = block % heatKoeffs   (:,:,1)
         block % heatKoeffs   (:,:,block % nPkts(3)-i)   = block % heatKoeffs   (:,:,block % nCells(3))
      end do
!      block % boundary(DIR_EAST ) % bc_type = BC_OUTFLOW
!      block % boundary(DIR_WEST ) % bc_type = BC_INFLOW
!      block % boundary(DIR_SOUTH) % bc_type = BC_WALL
!      block % boundary(DIR_NORTH) % bc_type = BC_SYMMETRY
!      block % boundary(DIR_FRONT) % bc_type = BC_SYMMETRY
!      block % boundary(DIR_BACK ) % bc_type = BC_SYMMETRY
   end subroutine init_boundary
   subroutine set_boundary(block,nBoundaryCells)
   implicit none
      type(block_type), intent(inout) :: block
      integer, intent(in) :: nBoundaryCells
      integer :: i,j,k
      integer :: ib,i1,ig
      integer :: jb,j1,jg
      real(REAL_KIND) :: un
      !< U_Normal Mormalengeschwindigkeit


      !!!! OST   
      select case (block % boundary(DIR_EAST) % bc_type)
         case(BC_OUTFLOW)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(block % nPkts(1)-i,:,:,:) = block % vars(block % nCells(1),:,:,:)   
            end do
         case(BC_INFLOW) 
         case(BC_WALL)
         case(BC_SYMMETRY)
            do k = 1, block % nCells(3)
               do j = 1, block % nCells(2)
                  ig = block % nPkts(1) !!! GEMETRIE INDEX
                  do i = 0, - nBoundaryCells + 1, -1
                     i1 = block % nCells(1) + i !!! FELDZELLE INDEX
                     ib = block % nPkts(1)  - i !!! BOUNDARYZELLE INDEX

                     un =                                    &
                        ( block % vars (i1, j ,k ,2)         &
                        * block % cellFaceVecsI(1, ig, j, k) &
                        + block % vars (i1, j, k, 3)         &
                        * block % cellFaceVecsI(2, ig, j, k) &
                        + block % vars (i1, j, k, 4)         &
                        * block % cellFaceVecsI(3, ig, j, k) )

                     block % vars(ib,j,k,:) = block % vars(i1,j,k,:)   
                     block % vars(ib,j,k,2) = block % vars(i1,j,k,2) &
                                            - un * block % cellFaceVecsI(1,ig,j,k)
                     block % vars(ib,j,k,3) = block % vars(i1,j,k,3) &
                                            - un * block % cellFaceVecsI(2,ig,j,k)
                     block % vars(ib,j,k,4) = block % vars(i1,j,k,4) &
                                            - un * block % cellFaceVecsI(3,ig,j,k)
                  end do
               end do
            end do
      end select

      !!!! SOUTH 
      select case (block % boundary(DIR_SOUTH) % bc_type)
         case(BC_OUTFLOW)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,i,:,:) = block % vars(:,1,:,:)   
            end do
         case(BC_INFLOW) 
         case(BC_WALL)
            jg = 1          !!!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 1, block % nCells(3)
               do j = 0, - nBoundaryCells + 1, -1
                  j1 = 1 - j !!! FELDZELLE INDEX
                  jb = 0 + j !!! BOUNDARYZELLE INDEX
                  do i = 1, 20
                     block % vars(i,jb,k,VEC_RHO) =   block % vars(i,j1,k,VEC_RHO)
                     block % vars(i,jb,k,VEC_SPU) =   block % vars(i,j1,k,VEC_SPU)
                     block % vars(i,jb,k,VEC_SPV) = - block % vars(i,j1,k,VEC_SPV)
                     block % vars(i,jb,k,VEC_SPW) =   block % vars(i,j1,k,VEC_SPW)
                     block % vars(i,jb,k,VEC_ENE) =   block % vars(i,j1,k,VEC_ENE)
                     block % temperatures(i,jb,k) =   block % temperatures(i,j1,k)
                     block % viscosities (i,jb,k) =   block % viscosities (i,j1,k)
                     block % heatKoeffs  (i,jb,k) =   block % heatKoeffs  (i,j1,k)
                  end do
                  do i = 20, block % nCells(1)
                     block % vars(i,jb,k,VEC_RHO) =   block % vars(i,j1,k,VEC_RHO)
                     block % vars(i,jb,k,VEC_SPU) = - block % vars(i,j1,k,VEC_SPU)
                     block % vars(i,jb,k,VEC_SPV) = - block % vars(i,j1,k,VEC_SPV)
                     block % vars(i,jb,k,VEC_SPW) = - block % vars(i,j1,k,VEC_SPW)
                     block % vars(i,jb,k,VEC_ENE) =   block % vars(i,j1,k,VEC_ENE)
                     block % temperatures(i,jb,k) =   block % temperatures(i,j1,k)
                     block % viscosities (i,jb,k) =   block % viscosities (i,j1,k)
                     block % heatKoeffs  (i,jb,k) =   block % heatKoeffs  (i,j1,k)
                  end do
               end do
            end do
   
         case(BC_SYMMETRY)
            jg = 1          !!!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 1, block % nCells(3)
               do j = 0, - nBoundaryCells + 1, -1
                  j1 = 1 - j !!! FELDZELLE INDEX
                  jb = 0 + j !!! BOUNDARYZELLE INDEX
                  do i = 1, block % nCells(1)
                     un = block % vars (i, j1, k ,2)         &
                        * block % cellFaceVecsJ(1, i, jg, k) &
                        + block % vars (i, j1, k, 3)         &
                        * block % cellFaceVecsJ(2, i, jg, k) &
                        + block % vars (i, j1, k, 4)         &
                        * block % cellFaceVecsJ(3, i, jg, k) 
!                     write(*,*) i,j,k,un
                     block % vars(i,jb,k,:) = block % vars(i,j1,k,:)   
                     block % vars(i,jb,k,2) = block % vars(i,j1,k,2) &
                                            - un * block % cellFaceVecsJ(1,i,jg,k)
                     block % vars(i,jb,k,3) = block % vars(i,j1,k,3) &
                                            - un * block % cellFaceVecsJ(2,i,jg,k)
                     block % vars(i,jb,k,4) = block % vars(i,j1,k,4) &
                                            - un * block % cellFaceVecsJ(3,i,jg,k)
                     if (i == 25 .and. jb == 0) then
!                        write(*,*) block % vars(i,0:1,k,2:3) , un
!                        write(*,*) block % cellFaceVecsJ(1:3,i,jg,k)

                     end if
                  end do
               end do
            end do
            i = 25
            k = 1
      end select

      !!!! NORTH 
      select case (block % boundary(DIR_NORTH) % bc_type)
         case(BC_OUTFLOW)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,block % nPkts(2)-i,:,:) = block % vars(:,block % nCells(2),:,:)   
            end do
         case(BC_INFLOW) 
         case(BC_WALL)
         case(BC_SYMMETRY)
            jg = block % nPkts (2)        !!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 1, block % nCells(3)
               do j = 0, - nBoundaryCells + 1, -1
               j1 = block % nCells(2) + j !!! FELDZELLE INDEX
               jb = block % nPkts (2) - j !!! BOUNDARYZELLE INDEX
                  do i = 1, block % nCells(1)
                     un = block % vars (i, j1, k ,2)         &
                        * block % cellFaceVecsJ(1, i, jg, k) &
                        + block % vars (i, j1, k, 3)         &
                        * block % cellFaceVecsJ(2, i, jg, k) &
                        + block % vars (i, j1, k, 4)         &
                        * block % cellFaceVecsJ(3, i, jg, k) 
!                     write(*,*) i,j,k,un
                     block % vars(i,jb,k,:) = block % vars(i,j1,k,:)   
                     block % vars(i,jb,k,2) = block % vars(i,j1,k,2) &
                                            - un * block % cellFaceVecsJ(1,i,jg,k)
                     block % vars(i,jb,k,3) = block % vars(i,j1,k,3) &
                                            - un * block % cellFaceVecsJ(2,i,jg,k)
                     block % vars(i,jb,k,4) = block % vars(i,j1,k,4) &
                                            - un * block % cellFaceVecsJ(3,i,jg,k)
                  end do
               end do
            end do
      end select
      do i = 0, - nBoundaryCells + 1, -1
!         block % vars(i,:,:,:) = block % vars(1,:,:,:)                                       ! WEST
!         block % vars(block % nPkts(1)-i,:,:,:) = block % vars(block % nCells(1),:,:,:)      ! EAST
      end do
      do i = 0, - nBoundaryCells + 1, -1
!         block % vars(:,i,:,:) = block % vars(:,1,:,:)                                       ! SOUTH
!         block % vars(:,block % nPkts(2)-i,:,:) = block % vars(:,block % nCells(2),:,:)      ! NORTH
      end do
      do i = 0, - nBoundaryCells + 1, -1
         block % vars(:,:,i,:) = block % vars(:,:,1,:)                                       ! FRONT
         block % vars(:,:,block % nPkts(3)-i,:) = block % vars(:,:,block % nCells(3),:)      ! BACK
      end do
   end subroutine set_boundary
end module boundary_mod
