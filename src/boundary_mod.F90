module boundary_mod
   use const_mod
   use data_mod, only: block_type

implicit none
   real(kind = REAL_KIND) :: inflow_vel
contains
   subroutine init_boundary(block,nBoundaryCells)
      use screen_io_mod ,only : error_wr
   implicit none
      type(block_type), intent(inout) :: block
      integer, intent(in) :: nBoundaryCells
      
      integer :: i,j,k,bc_dir
   
      do i = 0, - nBoundaryCells + 1, -1
         block % vars         (i,:,:,:)                  = block % vars         (1,:,:,:)
         block % vars         (block % nPkts(1)-i,:,:,:) = block % vars         (block % nCells(1),:,:,:)   
         block % vars         (:,i,:,:)                  = block % vars         (:,1,:,:)
         block % vars         (:,block % nPkts(2)-i,:,:) = block % vars         (:,block % nCells(2),:,:)   
         block % vars         (:,:,i,:)                  = block % vars         (:,:,1,:)
         block % vars         (:,:,block % nPkts(3)-i,:) = block % vars         (:,:,block % nCells(3),:)   
         block % pressures    (i,:,:)                    = block % pressures    (1,:,:)
         block % pressures    (block % nPkts(1)-i,:,:)   = block % pressures    (block % nCells(1),:,:)
         block % pressures    (:,i,:)                    = block % pressures    (:,1,:)
         block % pressures    (:,block % nPkts(2)-i,:)   = block % pressures    (:,block % nCells(2),:)
         block % pressures    (:,:,i)                    = block % pressures    (:,:,1)
         block % pressures    (:,:,block % nPkts(3)-i)   = block % pressures    (:,:,block % nCells(3))
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

      do bc_dir = 1, 6
         if (block % boundary(bc_dir) % bc_type == BC_INFLOW) then
            if (bc_dir == DIR_WEST) then
               inflow_vel = block % vars (1,1,1,VEC_SPU)
               do k = 1, block % nCells(3) 
                  do j = 1, block % nCells(2)
                     if ( abs (block % vars (1,j,k,VEC_SPU) - inflow_vel) > 1E-10) then
                        call error_wr("Inflow Profile is not constant",__FILE__,__LINE__)
                     end if
                  end do
               end do
               if (abs(block % vars(1,1,1,VEC_SPV)) > 1.0D-10) then
                  call error_wr("SPU V not Zero at Inflow",__FILE__,__LINE__)
               end if
               if (abs(block % vars(1,1,1,VEC_SPW)) > 1.0D-10) then
                  call error_wr("SPU W not Zero at Inflow",__FILE__,__LINE__)
               end if
               write(*,*) "Inflow-Velocity at WEST",inflow_vel
            else
               call error_wr("Inflow Boundary only possible at EAST Border",__FILE__,__LINE__)
            end if
         end if
      end do
         
!      block % boundary(DIR_EAST ) % bc_type = BC_OUTFLOW
!      block % boundary(DIR_WEST ) % bc_type = BC_INFLOW
!      block % boundary(DIR_SOUTH) % bc_type = BC_WALL
!      block % boundary(DIR_NORTH) % bc_type = BC_SYMMETRY
!      block % boundary(DIR_FRONT) % bc_type = BC_SYMMETRY
!      block % boundary(DIR_BACK ) % bc_type = BC_SYMMETRY
   end subroutine init_boundary

   subroutine set_boundary(blocks,b,nBoundaryCells)
      use control_mod ,only: pressure_out
   implicit none
      type(block_type), intent(inout) :: blocks(:)
      integer, intent(in) :: b
      integer, intent(in) :: nBoundaryCells
      integer :: i,j,k
      integer :: ib,i1,ig
      integer :: jb,j1,jg
      integer :: kb,k1,kg
      !integer :: ob !other block
      real(REAL_KIND) :: un
      real(REAL_KIND) :: deltap
      real(REAL_KIND) :: sof ! Speed of Sound
      real(REAL_KIND) :: rhosof ! 1 / (RHO *  Speed of Sound)
      !< U_Normal Mormalengeschwindigkeit

      associate (block => blocks(b))
! ****************************************************************************************************
! ****************************************************************************************************
! ************************       WEST      ***********************************************************
! ****************************************************************************************************
! ****************************************************************************************************
      select case (block % boundary(DIR_WEST) % bc_type)
         case(BC_PERIODIC)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(i,:,:,:) = block % vars(block % nCells(1)+i,:,:,:)   
            end do
         case(BC_INFLOW) 
            ig = 1 !!! GEMETRIE INDEX
            i1 = 1
            do k = 1, block % nCells(3)
               do j = 1, block % nCells(2)
                  do i  = 1, nBoundaryCells
                     ib = 1-i
                     sof = sqrt(GAMMA * RGas * block % temperatures(i1,j,k) )

                     rhosof = block % vars(i1,j,k,VEC_RHO) * sof

                     block % pressures(ib,j,k)    = 0.5E0_REAL_KIND * ( pressure_out + block % pressures(i1,j,k)  &
                           - rhosof * block % abscellFaceVecsI(1,ig,j,k) * ( inflow_vel - block % vars(i1,j,k,VEC_SPU) )) 

                     deltap = block % pressures(i1,j,k) - pressure_out

                     sof = 1.0E0_REAL_KIND  / sof

                     rhosof = 1.0E0_REAL_KIND / rhosof

                     block % vars(ib,j,k,VEC_RHO) = block % vars (i1,j,k,VEC_RHO) - deltap * sof * sof

                     block % vars(ib,j,k,VEC_SPU) = block % vars (i1,j,k,VEC_SPU) + deltap * rhosof &
                                                  * block % abscellFaceVecsI(1,ig,j,k)

                     block % vars(ib,j,k,VEC_SPV) = block % vars (i1,j,k,VEC_SPV) + deltap * rhosof &
                                                  * block % abscellFaceVecsI(2,ig,j,k)

                     block % vars(ib,j,k,VEC_SPW) = block % vars (i1,j,k,VEC_SPW) + deltap * rhosof &
                                                  * block % abscellFaceVecsI(3,ig,j,k)

                     block % pressures(ib,j,k)    = pressure_out

                     block % vars(ib,j,k,VEC_ENE) = block % pressures(ib,j,k) / GM1 &
                                                  + 0.5E0_REAL_KIND                     &
                                                  * block % vars(ib,j,k,VEC_RHO)         &
                                                  * ( block % vars (ib,j,k,VEC_SPU)  &
                                                    * block % vars (ib,j,k,VEC_SPU)  & 
                                                    + block % vars (ib,j,k,VEC_SPV)  &
                                                    * block % vars (ib,j,k,VEC_SPV)  &
                                                    + block % vars (ib,j,k,VEC_SPW)  &
                                                    * block % vars (ib,j,k,VEC_SPW)  &
                                                    )
                  end do
               end do
            end do
!            do i = 0, - nBoundaryCells + 1, -1
!               block % pressures   (i,:,:) = block % pressures(block % nCells(1) + i,:,:)
!               block % vars(i,:,:,VEC_ENE) = block % pressures(i,:,:) / (GM1)    &
!                                           + 0.5E0_REAL_KIND                     &
!                                           * block % vars(i,:,:,VEC_RHO)         &
!                                           * ( block % vars (i,:,:,VEC_SPU)  &
!                                             * block % vars (i,:,:,VEC_SPU)  & 
!                                             + block % vars (i,:,:,VEC_SPV)  &
!                                             * block % vars (i,:,:,VEC_SPV)  &
!                                             + block % vars (i,:,:,VEC_SPW)  &
!                                             * block % vars (i,:,:,VEC_SPW)  &
!                                             )
!            end do
         case(BC_SYMMETRY)
            ig = 1          !!!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 1, block % nCells(3)
               do j = 0, block % nCells(2)
                  do i = 1,  nBoundaryCells
                     i1 = i     !!! FELDZELLE INDEX
                     ib = 1 - i !!! BOUNDARYZELLE INDEX
                     un = block % vars (i1, j, k ,2)         &
                        * block % cellFaceVecsI(1, ig, j, k) &
                        + block % vars (i1, j, k, 3)         &
                        * block % cellFaceVecsI(2, ig, j, k) &
                        + block % vars (i1, j, k, 4)         &
                        * block % cellFaceVecsI(3, ig, j, k) 
                     block % vars(ib,j,k,:) = block % vars(i1,j,k,:)   
                     block % vars(ib,j,k,2) = block % vars(i1,j,k,2)             &
                                            - un * block % cellFaceVecsI(1,ig,j,k)
                     block % vars(ib,j,k,3) = block % vars(i1,j,k,3)             &
                                            - un * block % cellFaceVecsI(2,ig,j,k)
                     block % vars(ib,j,k,4) = block % vars(i1,j,k,4)             &
                                            - un * block % cellFaceVecsI(3,ig,j,k)
                  end do
               end do
            end do
         case(1:)
            associate (ob => blocks(block % boundary(DIR_WEST) % bc_type))
            do i = 0, - nBoundaryCells + 1, -1
               block % vars        (i,:,:,:) = ob % vars(ob % nCells(1) + i,:,:,:)
               block % pressures   (i,:,:)   = ob % pressures   (ob % nCells(1) + i,:,:)
               block % temperatures(i,:,:)   = ob % temperatures(ob % nCells(1) + i,:,:)
               block % viscosities (i,:,:)   = ob % viscosities(ob % nCells(1) + i,:,:)
               block % heatKoeffs  (i,:,:)   = ob % heatKoeffs(ob % nCells(1) + i,:,:)
            end do
            end associate

      end select

! ****************************************************************************************************
! ****************************************************************************************************
! ************************       OST       ***********************************************************
! ****************************************************************************************************
! ****************************************************************************************************
      select case (block % boundary(DIR_EAST) % bc_type)
         case(BC_OUTFLOW)
            ig = block % nPkts(1) !!! GEMETRIE INDEX
            i1 = block % nCells(1)
            do k = 1, block % nCells(3)
               do j = 1, block % nCells(2)
                  do i = 1, nBoundaryCells
                     ib = block % nCells(1)+i
                     deltap = block % pressures(i1,j,k) - pressure_out

                     sof = 1.0E0_REAL_KIND  / sqrt(GAMMA * RGas * block % temperatures(i1,j,k) )

                     rhosof = 1.0E0_REAL_KIND /  block % vars(i1,j,k,VEC_RHO) * sof

                     block % vars(ib,j,k,VEC_RHO) = block % vars (i1,j,k,VEC_RHO) - deltap * sof * sof

                     block % vars(ib,j,k,VEC_SPU) = block % vars (i1,j,k,VEC_SPU) + deltap * rhosof &
                                                  * block % abscellFaceVecsI(1,ig,j,k)

                     block % vars(ib,j,k,VEC_SPV) = block % vars (i1,j,k,VEC_SPV) + deltap * rhosof &
                                                  * block % abscellFaceVecsI(2,ig,j,k)

                     block % vars(ib,j,k,VEC_SPW) = block % vars (i1,j,k,VEC_SPW) + deltap * rhosof &
                                                  * block % abscellFaceVecsI(3,ig,j,k)

                     block % pressures(ib,j,k)    = pressure_out

                     block % vars(ib,j,k,VEC_ENE) = block % pressures(ib,j,k) / GM1 &
                                                  + 0.5E0_REAL_KIND                     &
                                                  * block % vars(ib,j,k,VEC_RHO)         &
                                                  * ( block % vars (ib,j,k,VEC_SPU)  &
                                                    * block % vars (ib,j,k,VEC_SPU)  & 
                                                    + block % vars (ib,j,k,VEC_SPV)  &
                                                    * block % vars (ib,j,k,VEC_SPV)  &
                                                    + block % vars (ib,j,k,VEC_SPW)  &
                                                    * block % vars (ib,j,k,VEC_SPW)  &
                                                    )
                  end do
               end do
            end do
         case(BC_INFLOW) 
         case(BC_WALL)
         case(BC_SYMMETRY)
            ig = block % nPkts(1) !!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 1, block % nCells(3)
               do j = 0, block % nCells(2)
                  do i = 1,  nBoundaryCells
                     i1 = block % nPkts(1) - i !!! FELDZELLE INDEX
                     ib = block % nCells(1) + i !!! BOUNDARYZELLE INDEX
                     un = block % vars (i1, j, k ,2)         &
                        * block % cellFaceVecsI(1, ig, j, k) &
                        + block % vars (i1, j, k, 3)         &
                        * block % cellFaceVecsI(2, ig, j, k) &
                        + block % vars (i1, j, k, 4)         &
                        * block % cellFaceVecsI(3, ig, j, k) 
                     block % vars(ib,j,k,:) = block % vars(i1,j,k,:)   
                     block % vars(ib,j,k,2) = block % vars(i1,j,k,2)             &
                                            - un * block % cellFaceVecsI(1,ig,j,k)
                     block % vars(ib,j,k,3) = block % vars(i1,j,k,3)             &
                                            - un * block % cellFaceVecsI(2,ig,j,k)
                     block % vars(ib,j,k,4) = block % vars(i1,j,k,4)             &
                                            - un * block % cellFaceVecsI(3,ig,j,k)
                  end do
               end do
            end do
         case(BC_PERIODIC)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(block % nPkts(1)-i,:,:,:) = block % vars(1-i,:,:,:)   
            end do
         case(1:)
            associate (ob => blocks(block % boundary(DIR_EAST) % bc_type))
            do i = 0, - nBoundaryCells + 1, -1
               block % vars        (block % nPkts(1)-i,:,:,:) = ob % vars        (1-i,:,:,:)
               block % pressures   (block % nPkts(1)-i,:,:)   = ob % pressures   (1-i,:,:)
               block % temperatures(block % nPkts(1)-i,:,:)   = ob % temperatures(1-i,:,:)
               block % viscosities (block % nPkts(1)-i,:,:)   = ob % viscosities (1-i,:,:)
               block % heatKoeffs  (block % nPkts(1)-i,:,:)   = ob % heatKoeffs  (1-i,:,:)
            end do
            end associate
      end select

! ****************************************************************************************************
! ****************************************************************************************************
! ************************       SOUTH     ***********************************************************
! ****************************************************************************************************
! ****************************************************************************************************
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
                  do i = 1, 32
                     block % vars(i,jb,k,VEC_RHO) =   block % vars(i,j1,k,VEC_RHO)
                     block % vars(i,jb,k,VEC_SPU) =   block % vars(i,j1,k,VEC_SPU)
                     block % vars(i,jb,k,VEC_SPV) = - block % vars(i,j1,k,VEC_SPV)
                     block % vars(i,jb,k,VEC_SPW) =   block % vars(i,j1,k,VEC_SPW)
                     block % vars(i,jb,k,VEC_ENE) =   block % vars(i,j1,k,VEC_ENE)
                     block % temperatures(i,jb,k) =   block % temperatures(i,j1,k)
                     block % viscosities (i,jb,k) =   block % viscosities (i,j1,k)
                     block % heatKoeffs  (i,jb,k) =   block % heatKoeffs  (i,j1,k)
                  end do
                  do i = 33, block % nCells(1)
                  !do i = 1, block % nCells(1)
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
                  end do
               end do
            end do
            i = 25
            k = 1
         case(BC_PERIODIC)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,i,:,:) = block % vars(:,block % nCells(2)+i,:,:)   
            end do
      end select

! ****************************************************************************************************
! ****************************************************************************************************
! ************************       NORTH     ***********************************************************
! ****************************************************************************************************
! ****************************************************************************************************
      select case (block % boundary(DIR_NORTH) % bc_type)
         case(BC_OUTFLOW)
            jg = block % nPkts(2) !!! GEMETRIE INDEX
            j1 = block % nCells(2)
            do k = 1, block % nCells(3)
               do j = 1, nBoundaryCells
                  jb = block % nCells(2)+j
                  do i = 1, block % nCells(1)
                     deltap = block % pressures(i,j1,k) - pressure_out

                     sof = 1.0E0_REAL_KIND  / sqrt(GAMMA * RGas * block % temperatures(i,j1,k) )

                     rhosof = 1.0E0_REAL_KIND /  block % vars(i,j1,k,VEC_RHO) * sof

                     block % vars(i,jb,k,VEC_RHO) = block % vars (i,j1,k,VEC_RHO) - deltap * sof * sof
                                      
                     block % vars(i,jb,k,VEC_SPU) = block % vars (i,j1,k,VEC_SPU) + deltap * rhosof &
                                                  * block % abscellFaceVecsJ(1,i,jg,k)
                                      
                     block % vars(i,jb,k,VEC_SPV) = block % vars (i,j1,k,VEC_SPV) + deltap * rhosof &
                                                  * block % abscellFaceVecsJ(2,i,jg,k)
                                      
                     block % vars(i,jb,k,VEC_SPW) = block % vars (i,j1,k,VEC_SPW) + deltap * rhosof &
                                                  * block % abscellFaceVecsJ(3,i,jg,k)

                     block % pressures(i,jb,k)    = pressure_out

                     block % vars(i,jb,k,VEC_ENE) = block % pressures(i,jb,k) / GM1 &
                                                  + 0.5E0_REAL_KIND                     &
                                                  * block % vars(i,jb,k,VEC_RHO)         &
                                                  * ( block % vars (i,jb,k,VEC_SPU)  &
                                                    * block % vars (i,jb,k,VEC_SPU)  & 
                                                    + block % vars (i,jb,k,VEC_SPV)  &
                                                    * block % vars (i,jb,k,VEC_SPV)  &
                                                    + block % vars (i,jb,k,VEC_SPW)  &
                                                    * block % vars (i,jb,k,VEC_SPW)  &
                                                    )
                  end do
               end do
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
         case(BC_PERIODIC)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,block % nPkts(2)-i,:,:) = block % vars(:,1-i,:,:)   
            end do
      end select

! ****************************************************************************************************
! ****************************************************************************************************
! ************************       FRONT     ***********************************************************
! ****************************************************************************************************
! ****************************************************************************************************
      select case (block % boundary(DIR_FRONT) % bc_type)
         case(BC_PERIODIC)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,:,i,:) = block % vars(:,:,block % nCells(3)+i,:)   
            end do
         case(BC_SYMMETRY)
            kg = 1        !!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 0, - nBoundaryCells + 1, -1
               k1 = 1 - k !!! FELDZELLE INDEX
               kb = 0 + k !!! BOUNDARYZELLE INDEX
               do j = 1, block % nCells(2)
                  do i = 1, block % nCells(1)
                     un = block % vars (i, j, k1 ,2)         &
                        * block % cellFaceVecsK(1, i, j, kg) &
                        + block % vars (i, j, k1, 3)         &
                        * block % cellFaceVecsK(2, i, j, kg) &
                        + block % vars (i, j, k1, 4)         &
                        * block % cellFaceVecsK(3, i, j, kg) 
                     block % vars(i,j,kb,:) = block % vars(i,j,k1,:)   
                     block % vars(i,j,kb,2) = block % vars(i,j,k1,2) &
                                            - un * block % cellFaceVecsK(1,i,j,kg)
                     block % vars(i,j,kb,3) = block % vars(i,j,k1,3) &
                                            - un * block % cellFaceVecsK(2,i,j,kg)
                     block % vars(i,j,kb,4) = block % vars(i,j,k1,4) &
                                            - un * block % cellFaceVecsK(3,i,j,kg)
                  end do
               end do
            end do
         case  default
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,:,i,:) = block % vars(:,:,1,:)
            end do
      end select

! ****************************************************************************************************
! ****************************************************************************************************
! ************************       BACK      ***********************************************************
! ****************************************************************************************************
! ****************************************************************************************************
      select case (block % boundary(DIR_BACK) % bc_type)
         case(BC_PERIODIC)
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,:,block % nPkts(3)-i,:) = block % vars(:,:,1-i,:)   
            end do
         case(BC_SYMMETRY)
            kg = block % nPkts (3)        !!! GEOMETRIE ZELLE (aendert sich nicht für due unterschiedlichen BOUNDARY ZELLEN 
            do k = 0, - nBoundaryCells + 1, -1
            k1 = block % nCells(3) + k !!! FELDZELLE INDEX
            kb = block % nPkts (3) - k !!! BOUNDARYZELLE INDEX
               do j = 1, block % nCells(2)
                  do i = 1, block % nCells(1)
                     un = block % vars (i, j, k1 ,2)         &
                        * block % cellFaceVecsK(1, i, j, kg) &
                        + block % vars (i, j, k1, 3)         &
                        * block % cellFaceVecsK(2, i, j, kg) &
                        + block % vars (i, j, k1, 4)         &
                        * block % cellFaceVecsK(3, i, j, kg) 
                     block % vars(i,j,kb,:) = block % vars(i,j,k1,:)   
                     block % vars(i,j,kb,2) = block % vars(i,j,k1,2) &
                                            - un * block % cellFaceVecsK(1,i,j,kg)
                     block % vars(i,j,kb,3) = block % vars(i,j,k1,3) &
                                            - un * block % cellFaceVecsK(2,i,j,kg)
                     block % vars(i,j,kb,4) = block % vars(i,j,k1,4) &
                                            - un * block % cellFaceVecsK(3,i,j,kg)
                  end do
               end do
            end do
         case  default
            do i = 0, - nBoundaryCells + 1, -1
               block % vars(:,:,block % nPkts(3)-i,:) = block % vars(:,:,block % nCells(3),:)
            end do
      end select

! ================================================================================
! ================================================================================
! ================================================================================
! ================================================================================
! ================================================================================
      do i = 0, - nBoundaryCells + 1, -1
!         block % vars(i,:,:,:) = block % vars(1,:,:,:)                                       ! WEST
!         block % vars(block % nPkts(1)-i,:,:,:) = block % vars(block % nCells(1),:,:,:)      ! EAST
      end do
      do i = 0, - nBoundaryCells + 1, -1
!         block % vars(:,i,:,:) = block % vars(:,1,:,:)                                       ! SOUTH
!         block % vars(:,block % nPkts(2)-i,:,:) = block % vars(:,block % nCells(2),:,:)      ! NORTH
      end do
      end associate
   end subroutine set_boundary
end module boundary_mod
