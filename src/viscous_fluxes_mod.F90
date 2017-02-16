module viscous_fluxes_mod
implicit none
contains
   subroutine calc_face_gradients(block)
      use const_mod
      use data_mod, only: block_type,Dimen
   implicit none 

   type(block_type), intent(inout) :: block

   integer :: i,j,k,v,dir 
   do dir = 1, Dimen  !!! schleife über die Richtungen dx/dy/dz
      do v  = 1, Dimen  !!! schleife über die Geschwindigkeiten U = (u,v,w)
         do k = 1, block % nCells(3)
            do j = 1, block % nCells(2)
               do i = 1, block % nPkts(1)
                  block % dUdnI ( i,j,k,v,dir) = (block % vars (i,j,k,v+1) - block % vars(i-1,j,k,v+1) ) &
                                               * block % swpDistVecsI(dir,i,j,k)
               end do
            end do
         end do
      end do
      do k = 1, block % nCells(3)
         do j = 1, block % nCells(2)
            do i = 1, block % nPkts(1)
               block % dUdnI ( i,j,k,GRAD_TEMP,dir) = (block % temperatures (i,j,k) - block % temperatures(i-1,j,k) ) &
                                                    * block % swpDistVecsI(dir,i,j,k)
            end do
         end do
      end do
   end do
   do dir = 1, Dimen  !!! schleife über die Richtungen dx/dy/dz
      do v = 1, Dimen  !!! schleife über die Geschwindigkeiten U = (u,v,w)
         do k = 1, block % nCells(3)
            do j = 1, block % nPkts(2)
               do i = 1, block % nCells(1)
                  block % dUdnJ ( i,j,k,v,dir) = (block % vars (i,j,k,v+1) - block % vars(i,j-1,k,v+1) ) &
                                               * block % swpDistVecsJ(dir,i,j,k)
!                  if (i == 10 .and. j == 1) then
!                     write(*,*) dir,v,block % dUdnJ (i,j,k,v,dir),block % vars (i,j,k,v+1) &
!                              ,block % vars(i,j-1,k,v+1)  , block % swpDistVecsJ(dir,i,j,k)
!                  end if
               end do
            end do
         end do
      end do
      do k = 1, block % nCells(3)
         do j = 1, block % nPkts(2)
            do i = 1, block % nCells(1)
               block % dUdnJ ( i,j,k,GRAD_TEMP,dir) = (block % temperatures (i,j,k) - block % temperatures(i,j-1,k) ) &
                                                    * block % swpDistVecsJ(dir,i,j,k)
            end do
         end do
      end do
   end do
   do dir = 1, Dimen  !!! schleife über die Richtungen dx/dy/dz
      do v = 1, Dimen  !!! schleife über die Geschwindigkeiten U = (u,v,w)
         do k = 1, block % nPkts(3)
            do j = 1, block % nCells(2)
               do i = 1, block % nCells(1)
                  block % dUdnK ( i,j,k,v,dir) = (block % vars (i,j,k,v+1) - block % vars(i,j,k-1,v+1) ) &
                                               * block % swpDistVecsK(dir,i,j,k)
               end do
            end do
         end do
      end do
      do k = 1, block % nPkts(3)
         do j = 1, block % nCells(2)
            do i = 1, block % nCells(1)
               block % dUdnK ( i,j,k,GRAD_TEMP,dir) = (block % temperatures (i,j,k) - block % temperatures(i,j,k-1) ) &
                                                    * block % swpDistVecsK(dir,i,j,k)
            end do
         end do
      end do
   end do
   end subroutine calc_face_gradients
   subroutine calc_viscous_fluxes(block)
      use const_mod
      use data_mod, only: block_type,dimen
   implicit none 

   type(block_type), intent(inout) :: block

   integer :: i,j,k
   real(REAL_KIND) :: mu_local
   real(REAL_KIND) :: la_local
   real(REAL_KIND) :: sp_local(dimen)

   do k = 1, block % nCells(3)
      do j = 1, block % nCells(2)
         do i = 1, block % nPkts(1)
            mu_local = (block % viscosities (i,j,k) + block % viscosities (i-1,j,k) ) * HALF
            la_local = (block % heatKoeffs  (i,j,k) + block % heatKoeffs  (i-1,j,k) ) * HALF
            sp_local = (block % vars  (i,j,k,2:4) + block % vars  (i-1,j,k,2:4) ) * HALF
            block % tauI (i,j,k,TAU_TXX) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnI ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnI ( i,j,k,GRAD_SPV,GRAD_DY) &
                                           -       block % dUdnI ( i,j,k,GRAD_SPW,GRAD_DZ) )
            block % tauI (i,j,k,TAU_TYY) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnI ( i,j,k,GRAD_SPV,GRAD_DY) &
                                           -       block % dUdnI ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnI ( i,j,k,GRAD_SPW,GRAD_DZ) )
            block % tauI (i,j,k,TAU_TZZ) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnI ( i,j,k,GRAD_SPW,GRAD_DZ) &
                                           -       block % dUdnI ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnI ( i,j,k,GRAD_SPV,GRAD_DY) )
            block % tauI (i,j,k,TAU_TXY) = mu_local &
                                         * (       block % dUdnI ( i,j,k,GRAD_SPU,GRAD_DY) &
                                           +       block % dUdnI ( i,j,k,GRAD_SPV,GRAD_DX) )
            block % tauI (i,j,k,TAU_TXZ) = mu_local &
                                         * (       block % dUdnI ( i,j,k,GRAD_SPU,GRAD_DZ) &
                                           +       block % dUdnI ( i,j,k,GRAD_SPW,GRAD_DX) )
            block % tauI (i,j,k,TAU_TYZ) = mu_local &
                                         * (       block % dUdnI ( i,j,k,GRAD_SPV,GRAD_DZ) &
                                           +       block % dUdnI ( i,j,k,GRAD_SPW,GRAD_DY) )
            block % heatfluxI (i,j,k,1:dimen) = - la_local * block % dudnI(i,j,k,GRAD_TEMP,1:dimen)

            block % visFluxesI(i,j,k,2) = - ( block % abscellFaceVecsI(1,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TXX)    &
                                            + block % abscellFaceVecsI(2,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TXY)    &
                                            + block % abscellFaceVecsI(3,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TXZ)    )
            block % visFluxesI(i,j,k,3) = - ( block % abscellFaceVecsI(1,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TXY)    &
                                            + block % abscellFaceVecsI(2,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TYY)    &
                                            + block % abscellFaceVecsI(3,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TYZ)    )
            block % visFluxesI(i,j,k,4) = - ( block % abscellFaceVecsI(1,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TXZ)    &
                                            + block % abscellFaceVecsI(2,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TYZ)    &
                                            + block % abscellFaceVecsI(3,i,j,k) &
                                            * block % tauI(i,j,k,TAU_TZZ)    )
            block % visFluxesI(i,j,k,5) = - ( block % abscellFaceVecsI(1,i,j,k) &
                                            * ( block % tauI(i,j,k,TAU_TXX) * sp_local(1) &
                                              + block % tauI(i,j,k,TAU_TXY) * sp_local(2) &
                                              + block % tauI(i,j,k,TAU_TXZ) * sp_local(3) &
                                              - block % heatfluxI(i,j,k,GRAD_DX) )      &
                                            + block % abscellFaceVecsI(2,i,j,k) &
                                            * ( block % tauI(i,j,k,TAU_TXY) * sp_local(1) &
                                              + block % tauI(i,j,k,TAU_TYY) * sp_local(2) &
                                              + block % tauI(i,j,k,TAU_TYZ) * sp_local(3) &
                                              - block % heatfluxI(i,j,k,GRAD_DY) )      &
                                            + block % abscellFaceVecsI(3,i,j,k) &
                                            * ( block % tauI(i,j,k,TAU_TXZ) * sp_local(1) &
                                              + block % tauI(i,j,k,TAU_TYZ) * sp_local(2) &
                                              + block % tauI(i,j,k,TAU_TZZ) * sp_local(3) &
                                              - block % heatfluxI(i,j,k,GRAD_DZ) ) )
         end do
      end do
   end do
   do k = 1, block % nCells(3)
      do j = 1, block % nPkts(2)
         do i = 1, block % nCells(1)
            mu_local = (block % viscosities(i,j,k) + block % viscosities (i,j-1,k) ) * HALF
            la_local = (block % heatKoeffs  (i,j,k) + block % heatKoeffs  (i,j-1,k) ) * HALF
            sp_local = (block % vars  (i,j,k,2:4) + block % vars  (i,j-1,k,2:4) ) * HALF
            block % tauJ (i,j,k,TAU_TXX) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnJ ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnJ ( i,j,k,GRAD_SPV,GRAD_DY) &
                                           -       block % dUdnJ ( i,j,k,GRAD_SPW,GRAD_DZ) )
            block % tauJ (i,j,k,TAU_TYY) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnJ ( i,j,k,GRAD_SPV,GRAD_DY) &
                                           -       block % dUdnJ ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnJ ( i,j,k,GRAD_SPW,GRAD_DZ) )
            block % tauJ (i,j,k,TAU_TZZ) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnJ ( i,j,k,GRAD_SPW,GRAD_DZ) &
                                           -       block % dUdnJ ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnJ ( i,j,k,GRAD_SPV,GRAD_DY) )
            block % tauJ (i,j,k,TAU_TXY) = mu_local &
                                         * (       block % dUdnJ ( i,j,k,GRAD_SPU,GRAD_DY) &
                                           +       block % dUdnJ ( i,j,k,GRAD_SPV,GRAD_DX) )
            block % tauJ (i,j,k,TAU_TXZ) = mu_local &
                                         * (       block % dUdnJ ( i,j,k,GRAD_SPU,GRAD_DZ) &
                                           +       block % dUdnJ ( i,j,k,GRAD_SPW,GRAD_DX) )
            block % tauJ (i,j,k,TAU_TYZ) = mu_local &
                                         * (       block % dUdnJ ( i,j,k,GRAD_SPV,GRAD_DZ) &
                                           +       block % dUdnJ ( i,j,k,GRAD_SPW,GRAD_DY) )
            block % heatfluxJ (i,j,k,1:dimen) = - la_local * block % dudnJ(i,j,k,GRAD_TEMP,1:dimen)

            block % visFluxesJ(i,j,k,2) = - ( block % absCellFaceVecsJ(1,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TXX)    &
                                            + block % absCellFaceVecsJ(2,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TXY)    &
                                            + block % absCellFaceVecsJ(3,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TXZ)    )
!            if (i == 10 .and. j <= 2) then
!               write(*,*) block%tauJ(i,j,k,TAU_TXY), block % absCellFaceVecsJ(2,i,j,k),block % visFluxesJ(i,j,k,2), mu_local &
!                                                  ,block % dUdnJ ( i,j,k,GRAD_SPU,GRAD_DY) &
!                                                  ,block % dUdnJ ( i,j,k,GRAD_SPV,GRAD_DX)
!            end if

            block % visFluxesJ(i,j,k,3) = - ( block % absCellFaceVecsJ(1,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TXY)    &
                                            + block % absCellFaceVecsJ(2,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TYY)    &
                                            + block % absCellFaceVecsJ(3,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TYZ)    )

            block % visFluxesJ(i,j,k,4) = - ( block % absCellFaceVecsJ(1,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TXZ)    &
                                            + block % absCellFaceVecsJ(2,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TYZ)    &
                                            + block % absCellFaceVecsJ(3,i,j,k) &
                                            * block % tauJ(i,j,k,TAU_TZZ)    )

            block % visFluxesJ(i,j,k,5) = - ( block % absCellFaceVecsJ(1,i,j,k) &
                                            * ( block % tauJ(i,j,k,TAU_TXX) * sp_local(1) &
                                              + block % tauJ(i,j,k,TAU_TXY) * sp_local(2) &
                                              + block % tauJ(i,j,k,TAU_TXZ) * sp_local(3) &
                                              - block % heatfluxJ(i,j,k,GRAD_DX) )      &
                                            + block % absCellFaceVecsJ(2,i,j,k) &
                                            * ( block % tauJ(i,j,k,TAU_TXY) * sp_local(1) &
                                              + block % tauJ(i,j,k,TAU_TYY) * sp_local(2) &
                                              + block % tauJ(i,j,k,TAU_TYZ) * sp_local(3) &
                                              - block % heatfluxJ(i,j,k,GRAD_DY) )      &
                                            + block % absCellFaceVecsJ(3,i,j,k) &
                                            * ( block % tauJ(i,j,k,TAU_TXZ) * sp_local(1) &
                                              + block % tauJ(i,j,k,TAU_TYZ) * sp_local(2) &
                                              + block % tauJ(i,j,k,TAU_TZZ) * sp_local(3) &
                                              - block % heatfluxJ(i,j,k,GRAD_DZ) ) )
         end do
      end do
   end do
   do k = 1, block % nPkts(3)
      do j = 1, block % nCells(2)
         do i = 1, block % nCells(1)
            mu_local = (block % viscosities(i,j,k) + block % viscosities (i,j,k-1) ) * HALF
            la_local = (block % heatKoeffs  (i,j,k) + block % heatKoeffs  (i,j,k-1) ) * HALF
            sp_local = (block % vars  (i,j,k,2:4) + block % vars  (i,j,k-1,2:4) ) * HALF
            block % tauK (i,j,k,TAU_TXX) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnK ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnK ( i,j,k,GRAD_SPV,GRAD_DY) &
                                           -       block % dUdnK ( i,j,k,GRAD_SPW,GRAD_DZ) )
            block % tauK (i,j,k,TAU_TYY) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnK ( i,j,k,GRAD_SPV,GRAD_DY) &
                                           -       block % dUdnK ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnK ( i,j,k,GRAD_SPW,GRAD_DZ) )
            block % tauK (i,j,k,TAU_TZZ) = TWOTHIRD * mu_local &
                                         * ( TWO * block % dUdnK ( i,j,k,GRAD_SPW,GRAD_DZ) &
                                           -       block % dUdnK ( i,j,k,GRAD_SPU,GRAD_DX) &
                                           -       block % dUdnK ( i,j,k,GRAD_SPV,GRAD_DY) )
            block % tauK (i,j,k,TAU_TXY) = mu_local &
                                         * (       block % dUdnK ( i,j,k,GRAD_SPU,GRAD_DY) &
                                           +       block % dUdnK ( i,j,k,GRAD_SPV,GRAD_DX) )
            block % tauK (i,j,k,TAU_TXZ) = mu_local &
                                         * (       block % dUdnK ( i,j,k,GRAD_SPU,GRAD_DZ) &
                                           +       block % dUdnK ( i,j,k,GRAD_SPW,GRAD_DX) )
            block % tauK (i,j,k,TAU_TYZ) = mu_local &
                                         * (       block % dUdnK ( i,j,k,GRAD_SPV,GRAD_DZ) &
                                           +       block % dUdnK ( i,j,k,GRAD_SPW,GRAD_DY) )
            block % heatfluxK (i,j,k,1:dimen) = - la_local * block % dudnK(i,j,k,GRAD_TEMP,1:dimen)

            block % visFluxesK(i,j,k,2) = - ( block % abscellFaceVecsK(1,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TXX)    &
                                            + block % abscellFaceVecsK(2,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TXY)    &
                                            + block % abscellFaceVecsK(3,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TXZ)    )
            block % visFluxesK(i,j,k,3) = - ( block % abscellFaceVecsK(1,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TXY)    &
                                            + block % abscellFaceVecsK(2,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TYY)    &
                                            + block % abscellFaceVecsK(3,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TYZ)    )
            block % visFluxesK(i,j,k,4) = - ( block % abscellFaceVecsK(1,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TXZ)    &
                                            + block % abscellFaceVecsK(2,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TYZ)    &
                                            + block % abscellFaceVecsK(3,i,j,k) &
                                            * block % tauK(i,j,k,TAU_TZZ)    )
            block % visFluxesK(i,j,k,5) = - ( block % abscellFaceVecsK(1,i,j,k) &
                                            * ( block % tauK(i,j,k,TAU_TXX) * sp_local(1) &
                                              + block % tauK(i,j,k,TAU_TXY) * sp_local(2) &
                                              + block % tauK(i,j,k,TAU_TXZ) * sp_local(3) &
                                              - block % heatfluxK(i,j,k,GRAD_DX) )      &
                                            + block % abscellFaceVecsK(2,i,j,k) &
                                            * ( block % tauK(i,j,k,TAU_TXY) * sp_local(1) &
                                              + block % tauK(i,j,k,TAU_TYY) * sp_local(2) &
                                              + block % tauK(i,j,k,TAU_TYZ) * sp_local(3) &
                                              - block % heatfluxK(i,j,k,GRAD_DY) )      &
                                            + block % abscellFaceVecsK(3,i,j,k) &
                                            * ( block % tauK(i,j,k,TAU_TXZ) * sp_local(1) &
                                              + block % tauK(i,j,k,TAU_TYZ) * sp_local(2) &
                                              + block % tauK(i,j,k,TAU_TZZ) * sp_local(3) &
                                              - block % heatfluxK(i,j,k,GRAD_DZ) ) )
         end do
      end do
   end do

   end subroutine calc_viscous_fluxes
end module viscous_fluxes_mod
