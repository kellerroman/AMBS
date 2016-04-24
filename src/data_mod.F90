module data_mod
use const_mod
implicit none
type :: block_type

   integer(INT_KIND) :: nPkts(3)

   integer(INT_KIND) :: nCells(3)

   real(REAL_KIND), allocatable :: coords         (:,:,:,:)

   real(REAL_KIND), allocatable :: vars           (:,:,:,:)

   real(REAL_KIND), allocatable :: cons_vars      (:,:,:,:)

   real(REAL_KIND), allocatable :: temperatures   (:,:,:)

   real(REAL_KIND), allocatable :: pressures      (:,:,:)

   real(REAL_KIND), allocatable :: viscosities    (:,:,:)

   real(REAL_KIND), allocatable :: heatKoeffs     (:,:,:)

   real(REAL_KIND), allocatable :: residuals      (:,:,:,:)

   real(REAL_KIND), allocatable :: visFluxes      (:,:,:,:)

   real(REAL_KIND), allocatable :: volumes        (:,:,:)

   real(REAL_KIND), allocatable :: faceVarsLeftI  (:,:,:,:)
   !< Variable auf der Zellfläche in I-Richtung von der Linken Seite (i-1)

   real(REAL_KIND), allocatable :: faceVarsRightI (:,:,:,:)
   !< Variable auf der Zellfläche in I-Richtung von der Rechte nSeite (i)

   real(REAL_KIND), allocatable :: faceVarsLeftJ  (:,:,:,:)
   !< Variable auf der Zellfläche in J-Richtung von der Linken Seite (j-1)

   real(REAL_KIND), allocatable :: faceVarsRightJ(:,:,:,:)
   !< Variable auf der Zellfläche in J-Richtung von der Rechte nSeite (j)

   real(REAL_KIND), allocatable :: faceVarsLeftK  (:,:,:,:)
   !< Variable auf der Zellfläche in K-Richtung von der Linken Seite (k-1)

   real(REAL_KIND), allocatable :: faceVarsRightK(:,:,:,:)
   !< Variable auf der Zellfläche in K-Richtung von der Rechte nSeite (k)

   real(REAL_KIND), allocatable :: fluxesI        (:,:,:,:)
   real(REAL_KIND), allocatable :: fluxesJ        (:,:,:,:)
   real(REAL_KIND), allocatable :: fluxesK        (:,:,:,:)

   real(REAL_KIND), allocatable :: cellFaceAreasI (:,:,:)
   !< Fläche der Zellflächen(3D) & Zellkanten (2D) in I Richtung (EAST & WEST)
   !< NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable :: cellFaceAreasJ (:,:,:)
   !< Fläche der Zellflächen(3D) & Zellkanten (2D) in J Richtung (NORTH & SOUTH)
   !< NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable :: cellFaceAreasK (:,:,:)
   !< Fläche der Zellflächen(3D) & Zellkanten (2D) in K Richtung (FROMT & BACK)
   !< NCell(1), NCell(2), NPkt(3)

   real(REAL_KIND), allocatable :: cellFaceVecsI (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in I Richtung (EAST & WEST)
   !< Dimen, NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable :: cellFaceVecsJ (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in J Richtung (NORTH & SOUTH)
   !< Dimen, NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable :: cellFaceVecsK (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in K Richtung (FROMT & BACK)
   !< Dimen, NPkt(1), NCell(2), NPkt(3)

   real(REAL_KIND), allocatable :: swpDistVecsI (:,:,:,:)
   !< Vektor mit Schwerpunkt-Abständen für die Berechnung der Ableitungen in I Richtung (EAST & WEST)
   !< Dimen, NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable :: swpDistVecsJ (:,:,:,:)
   !< Vektor mit Schwerpunkt-Abständen für die Berechnung der Ableitungen in J Richtung (NORTH & SOUTH)
   !< Dimen, NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable :: swpDistVecsK (:,:,:,:)
   !< Vektor mit Schwerpunkt-Abständen für die Berechnung der Ableitungen in K Richtung (FROMT & BACK)
   !< Dimen, NPkt(1), NCell(2), NPkt(3)

end type block_type

type(block_type), allocatable   :: blocks      (:)
integer(INT_KIND)               :: nBlock
integer(INT_KIND)               :: dimen
integer(INT_KIND)               :: nFaces
integer(INT_KIND)               :: nCorners
integer(INT_KIND)               :: nCell
integer(INT_KIND)               :: nVar

contains
   subroutine allocate_vars()
      use control_mod, only: nBoundaryCells
   implicit none

   integer :: ib
   integer :: ci,cj,ck
   ! Number of Cells
   integer :: pi,pj,pk
   ! NUmber of Grid Points
   integer :: cis,cjs,cks
   ! Cellsindex Start including boundary cells
   integer :: cie,cje,cke
   ! Cellsindex End including boundary cells
   
      nvar = 5
      block_loop: do ib = 1, nBlock
      associate (b => blocks(ib))
         ci = b % nCells(1)
         cj = b % nCells(2)
         ck = b % nCells(3)
   
         pi = b % nPkts(1)
         pj = b % nPkts(2)
         pk = b % nPkts(3)

         cis = 1 - nBoundaryCells
         cjs = 1 - nBoundaryCells
         cks = 1 - nBoundaryCells

         cie = ci + nBoundaryCells
         cje = cj + nBoundaryCells
         cke = ck + nBoundaryCells

         allocate(b % coords           (       pi,pj,pk,dimen ) )
         allocate(b % vars             (cis:cie,cjs:cje,cks:cke,nvar  ) )
         allocate(b % temperatures     (cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % pressures        (cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % viscosities      (cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % heatKoeffs       (cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % cons_vars        (       ci,cj,ck,nvar  ) )
         allocate(b % residuals        (       ci,cj,ck,nvar  ) )
         allocate(b % visFluxes        (       pi,pj,pk,nvar  ) )

         allocate(b % faceVarsLeftI    (       pi,cj,ck,nvar  ) )
         allocate(b % faceVarsRightI   (       pi,cj,ck,nvar  ) )

         allocate(b % faceVarsLeftJ    (       ci,pj,ck,nvar  ) )
         allocate(b % faceVarsRightJ   (       ci,pj,ck,nvar  ) )

         allocate(b % faceVarsLeftK    (       ci,cj,pk,nvar  ) )
         allocate(b % faceVarsRightK   (       ci,cj,pk,nvar  ) )

         allocate(b % fluxesI          (       pi,cj,ck,nvar  ) )
         allocate(b % fluxesJ          (       ci,pj,ck,nvar  ) )
         allocate(b % fluxesK          (       ci,cj,pk,nvar  ) )

         allocate(b % volumes          (       ci,cj,ck       ) )
         allocate(b % CellFaceAreasI    (       pi,cj,ck       ) )
         allocate(b % CellFaceAreasJ    (       ci,pj,ck       ) )
         allocate(b % CellFaceAreasK    (       ci,cj,pk       ) )
         allocate(b % CellFaceVecsI     ( dimen,pi,cj,ck       ) )
         allocate(b % CellFaceVecsJ     ( dimen,ci,pj,ck       ) )
         allocate(b % CellFaceVecsK     ( dimen,ci,cj,pk       ) )
   
         allocate(b % swpDistVecsI      ( dimen,pi,cj,ck       ) )
         allocate(b % swpDistVecsJ      ( dimen,ci,pj,ck       ) )
         allocate(b % swpDistVecsK      ( dimen,ci,cj,pk       ) )
   
   
      end associate
   
      end do block_loop


   end subroutine allocate_vars
   subroutine calc_grid()

   implicit none
      integer     :: ib, i,j,k
      real(REAL_KIND), allocatable :: swps(:,:,:,:)
      real(REAL_KIND), allocatable :: dnw(:,:,:,:)
      real(REAL_KIND), allocatable :: dns(:,:,:,:)
      real(REAL_KIND), allocatable :: dnb(:,:,:,:)
      real(REAL_KIND), dimension(dimen) :: vec1, vec2

      do ib = 1,nBlock
      associate (b  => blocks(ib) ,&
                 ci => blocks(ib) % nCells(1) ,&
                 cj => blocks(ib) % nCells(2) ,&
                 ck => blocks(ib) % nCells(3) ,&
                 pi => blocks(ib) % nPkts(1) ,&
                 pj => blocks(ib) % nPkts(2) ,&
                 pk => blocks(ib) % nPkts(3) )
         allocate (swps(dimen,ci,cj,ck))
         allocate (dnw (dimen,pi,cj,ck))
         allocate (dns (dimen,ci,pj,ck))
         allocate (dnb (dimen,ci,cj,pk))
         !!
         !! CELL CENTER POINTS
         !!
         do k = 1, ck
            do j = 1, cj
               do i = 1, ci
                  swps(:,i,j,k) = 0.125E0_REAL_KIND * ( &
                                + b % coords(i  ,j  ,k  ,:) & 
                                + b % coords(i+1,j  ,k  ,:) &
                                + b % coords(i  ,j+1,k  ,:) &
                                + b % coords(i+1,j+1,k  ,:) &
                                + b % coords(i  ,j  ,k+1,:) &
                                + b % coords(i+1,j  ,k+1,:) &
                                + b % coords(i  ,j+1,k+1,:) &
                                + b % coords(i+1,j+1,k+1,:) )
               end do
            end do
         end do
         !!
         !! CELL VOLUME
         !!
         DO  K = 1,ck 
            DO J = 1,cj 
               DO I = 1,pi 
                  dnw(1,I,J,K) =  0.5D+0 * ( (b % coords(I,J,K,2) - b % coords(I,J+1,K+1,2))&
                                         *(b % coords(I,J,K+1,3)-b % coords(I,J+1,K,3))&
                                         -(b % coords(I,J,K,3)-b % coords(I,J+1,K+1,3))&
                                         *(b % coords(I,J,K+1,2)-b % coords(I,J+1,K,2)) )
                  dnw(2,I,J,K) =  - 0.5D+0 * ( (b % coords(I,J,K,3)-b % coords(I,J+1,K+1,3))&
                                         *(b % coords(I,J,K+1,1)-b % coords(I,J+1,K,1))&
                                         -(b % coords(I,J,K,1)-b % coords(I,J+1,K+1,1))&
                                         *(b % coords(I,J,K+1,3)-b % coords(I,J+1,K,3)) )
                  dnw(3,I,J,K) =  0.5D+0 * ( (b % coords(I,J,K,1)-b % coords(I,J+1,K+1,1))&
                                         *(b % coords(I,J,K+1,2)-b % coords(I,J+1,K,2))&
                                         -(b % coords(I,J,K,2)-b % coords(I,J+1,K+1,2))&
                                         *(b % coords(I,J,K+1,1)-b % coords(I,J+1,K,1)) )
               end do
            end do
         end do
     
         DO K = 1,ck 
            DO J = 1,pj 
               DO I = 1,ci 
                  dns(1,I,J,K) =  0.5D+0 * ( (b % coords(I,J,K,2)-b % coords(I+1,J,K+1,2))&
                                         *(b % coords(I+1,J,K,3)-b % coords(I,J,K+1,3))&
                                         -(b % coords(I,J,K,3)-b % coords(I+1,J,K+1,3))&
                                         *(b % coords(I+1,J,K,2)-b % coords(I,J,K+1,2)) )
                  dns(2,I,J,K) =  - 0.5D+0 * ( (b % coords(I,J,K,3)-b % coords(I+1,J,K+1,3))&
                                         *(b % coords(I+1,J,K,1)-b % coords(I,J,K+1,1))&
                                         -(b % coords(I,J,K,1)-b % coords(I+1,J,K+1,1))&
                                         *(b % coords(I+1,J,K,3)-b % coords(I,J,K+1,3)) )
                  dns(3,I,J,K) =  0.5D+0 * ( (b % coords(I,J,K,1)-b % coords(I+1,J,K+1,1))&
                                         *(b % coords(I+1,J,K,2)-b % coords(I,J,K+1,2))&
                                         -(b % coords(I,J,K,2)-b % coords(I+1,J,K+1,2))&
                                         *(b % coords(I+1,J,K,1)-b % coords(I,J,K+1,1)) )
               end do
            end do
         end do
     
         DO K = 1,pk 
            DO J = 1,cj  
               DO I = 1,ci 
                  dnb(1,I,J,K) =  0.5D+0 * ( (b % coords(I,J,K,2)-b % coords(I+1,J+1,K,2))&
                                         *(b % coords(I,J+1,K,3)-b % coords(I+1,J,K,3))&
                                         -(b % coords(I,J,K,3)-b % coords(I+1,J+1,K,3))&
                                         *(b % coords(I,J+1,K,2)-b % coords(I+1,J,K,2)) )
                  dnb(2,I,J,K) =  - 0.5D+0 * ( (b % coords(I,J,K,3)-b % coords(I+1,J+1,K,3))&
                                         *(b % coords(I,J+1,K,1)-b % coords(I+1,J,K,1))&
                                         -(b % coords(I,J,K,1)-b % coords(I+1,J+1,K,1))&
                                         *(b % coords(I,J+1,K,3)-b % coords(I+1,J,K,3)) )
                  dnb(3,I,J,K) =  0.5D+0 * ( (b % coords(I,J,K,1)-b % coords(I+1,J+1,K,1))&
                                         *(b % coords(I,J+1,K,2)-b % coords(I+1,J,K,2))&
                                         -(b % coords(I,J,K,2)-b % coords(I+1,J+1,K,2))&
                                         *(b % coords(I,J+1,K,1)-b % coords(I+1,J,K,1)) )
               end do
            end do
         end do
                 
         DO K = 1,ck 
            DO J = 1,cj 
               DO I = 1,ci 
                  b % volumes(I,J,K) = 0.D+0
     !                 - W + E -
                  vec1(1) = b % coords(I,J,K,1)  +b % coords(I,J+1,K,1)  +&
                           b % coords(I,J,K+1,1)  +b % coords(I,J+1,K+1,1)-&
                        b % coords(I+1,J,K,1)-b % coords(I+1,J+1,K,1)-&
                        b % coords(I+1,J,K+1,1)-b % coords(I+1,J+1,K+1,1)
                  vec1(2) = b % coords(I,J,K,2)  +b % coords(I,J+1,K,2)  +&
                        b % coords(I,J,K+1,2)  +b % coords(I,J+1,K+1,2)-&
                        b % coords(I+1,J,K,2)-b % coords(I+1,J+1,K,2)-&
                        b % coords(I+1,J,K+1,2)-b % coords(I+1,J+1,K+1,2)
                  vec1(3) = b % coords(I,J,K,3)  +b % coords(I,J+1,K,3)  +&
                        b % coords(I,J,K+1,3)  +b % coords(I,J+1,K+1,3)-&
                        b % coords(I+1,J,K,3)-b % coords(I+1,J+1,K,3)-&
                        b % coords(I+1,J,K+1,3)-b % coords(I+1,J+1,K+1,3)
                  b % volumes(I,J,K) = b % volumes(I,J,K)&
                       + vec1(1) * ( + DnW(2,I,J,K) + DnW(2,I+1,J,K) )&
                       + vec1(2) * ( - DnW(1,I,J,K) - DnW(1,I+1,J,K) )&
                       + vec1(3) * ( + DnW(3,I,J,K) + DnW(3,I+1,J,K) )
     
     !                 - S + N -
                  vec1(1) = b % coords(I,J,K,1)  +b % coords(I,J,K+1,1)  +&
                        b % coords(I+1,J,K,1)  +b % coords(I+1,J,K+1,1)-&
                        b % coords(I,J+1,K,1)-b % coords(I,J+1,K+1,1)-&
                        b % coords(I+1,J+1,K,1)-b % coords(I+1,J+1,K+1,1)
                  vec1(2) = b % coords(I,J,K,2)  +b % coords(I,J,K+1,2)  +&
                        b % coords(I+1,J,K,2)  +b % coords(I+1,J,K+1,2)-&
                        b % coords(I,J+1,K,2)-b % coords(I,J+1,K+1,2)-&
                        b % coords(I+1,J+1,K,2)-b % coords(I+1,J+1,K+1,2)
                  vec1(3) = b % coords(I,J,K,3)  +b % coords(I,J,K+1,3)  +&
                        b % coords(I+1,J,K,3)  +b % coords(I+1,J,K+1,3)-&
                        b % coords(I,J+1,K,3)-b % coords(I,J+1,K+1,3)-&
                        b % coords(I+1,J+1,K,3)-b % coords(I+1,J+1,K+1,3)
                  b % volumes(I,J,K) = b % volumes(I,J,K)&
                       + vec1(1) * ( + dns(2,I,J,K) + dns(2,I,J+1,K) )&
                       + vec1(2) * ( - dns(1,I,J,K) - dns(1,I,J+1,K) )&
                       + vec1(3) * ( + dns(3,I,J,K) + dns(3,I,J+1,K) )
     
     !                 - B + F -
                  vec1(1) = b % coords(I,J,K,1)  +b % coords(I+1,J,K,1)  +&
                        b % coords(I,J+1,K,1)  +b % coords(I+1,J+1,K,1)-&
                        b % coords(I,J,K+1,1)-b % coords(I+1,J,K+1,1)-&
                       b % coords(I,J+1,K+1,1)-b % coords(I+1,J+1,K+1,1)
                  vec1(2) = b % coords(I,J,K,2)  +b % coords(I+1,J,K,2)  +&
                        b % coords(I,J+1,K,2)  +b % coords(I+1,J+1,K,2)-&
                        b % coords(I,J,K+1,2)-b % coords(I+1,J,K+1,2)-&
                        b % coords(I,J+1,K+1,2)-b % coords(I+1,J+1,K+1,2)
                  vec1(3) = b % coords(I,J,K,3)  +b % coords(I+1,J,K,3)  +&
                        b % coords(I,J+1,K,3)  +b % coords(I+1,J+1,K,3)-&
                        b % coords(I,J,K+1,3)-b % coords(I+1,J,K+1,3)-&
                        b % coords(I,J+1,K+1,3)-b % coords(I+1,J+1,K+1,3)
                  b % volumes(I,J,K) = b % volumes(I,J,K)&
                       + vec1(1) * ( + dnb(2,I,J,K) + dnb(2,I,J,K+1) )&
                       + vec1(2) * ( - dnb(1,I,J,K) - dnb(1,I,J,K+1) )&
                       + vec1(3) * ( + dnb(3,I,J,K) + dnb(3,I,J,K+1) )
                  b % volumes(I,J,K) = 8.0E0_REAL_KIND / b % volumes(I,J,K) 
               end do
            end do
         end do
         !!
         !!    Cell Center Distanes for gradients
         !!
         do k = 1, ck
            do j = 1, cj
               do i = 2, ci
                  b % swpDistVecsI(:,i,j,k) = swps(:,i,j,k) - swps(:,i-1,j,k)
               end do
            end do
         end do
         do k = 1, ck
            do j = 2, cj
               do i = 1, ci
                  b % swpDistVecsJ(:,i,j,k) = swps(:,i,j,k) - swps(:,i,j-1,k)
               end do
            end do
         end do
         do k = 2, ck
            do j = 1, cj
               do i = 1, ci
                  b % swpDistVecsK(:,i,j,k) = swps(:,i,j,k) - swps(:,i,j,k-1)
               end do
            end do
         end do
         !!
         !!    non block boundary:
         !!
         do k = 1, ck
            do j = 1, cj
               i = 1 
               b % swpDistVecsI(:,i ,j,k) = 0.125E0_REAL_KIND * ( &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j  ,k  ,:) )& 
                             + b % coords(i+1,j  ,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j+1,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j+1,k  ,:) )& 
                             + b % coords(i+1,j+1,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j  ,k+1,:) )& 
                             + b % coords(i+1,j  ,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j+1,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j+1,k+1,:) )& 
                             + b % coords(i+1,j+1,k+1,:) )
               i = pi 
               b % swpDistVecsI(:,i ,j,k) = 0.125E0_REAL_KIND * ( &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i-1,j  ,k  ,:) )& 
                             + b % coords(i-1,j  ,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j+1,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i-1,j+1,k  ,:) )& 
                             + b % coords(i-1,j+1,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i-1,j  ,k+1,:) )& 
                             + b % coords(i-1,j  ,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j+1,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i-1,j+1,k+1,:) )& 
                             + b % coords(i-1,j+1,k+1,:) )
            end do
         end do

         do k = 1, ck
            do i = 1, ci
               j = 1
               b % swpDistVecsJ(:,i ,j,k) = 0.125E0_REAL_KIND * ( &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j+1,k  ,:) )& 
                             + b % coords(i  ,j+1,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j+1,k  ,:) )& 
                             + b % coords(i+1,j+1,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j+1,k+1,:) )& 
                             + b % coords(i  ,j+1,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j  ,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j+1,k+1,:) )& 
                             + b % coords(i+1,j+1,k+1,:) )
               j = pj
               b % swpDistVecsJ(:,i ,j,k) = 0.125E0_REAL_KIND * ( &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j-1,k  ,:) )& 
                             + b % coords(i  ,j-1,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j-1,k  ,:) )& 
                             + b % coords(i+1,j-1,k  ,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j-1,k+1,:) )& 
                             + b % coords(i  ,j-1,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j  ,k+1,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j-1,k+1,:) )& 
                             + b % coords(i+1,j-1,k+1,:) )
            end do
         end do

         do j = 1, cj
            do i = 1, ci
               k = 1
               b % swpDistVecsK(:,i ,j,k) = 0.125E0_REAL_KIND * ( &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j  ,k+1,:) )& 
                             + b % coords(i  ,j  ,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j  ,k+1,:) )& 
                             + b % coords(i+1,j  ,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j+1,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j+1,k+1,:) )& 
                             + b % coords(i  ,j+1,k+1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j+1,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j+1,k+1,:) )& 
                             + b % coords(i+1,j+1,k+1,:) )
               k = pk
               b % swpDistVecsK(:,i ,j,k) = 0.125E0_REAL_KIND * ( &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j  ,k-1,:) )& 
                             + b % coords(i  ,j  ,k-1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j  ,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j  ,k-1,:) )& 
                             + b % coords(i+1,j  ,k-1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i  ,j+1,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i  ,j+1,k-1,:) )& 
                             + b % coords(i  ,j+1,k-1,:) &
                             - ( 2.0E0_REAL_KIND * b % coords(i+1,j+1,k  ,:) & 
                               - 1.0E0_REAL_KIND * b % coords(i+1,j+1,k-1,:) )& 
                             + b % coords(i+1,j+1,k-1,:) )
               b % swpDistVecsK(:,i,j,1 ) = b % swpDistVecsK(:,i,j,2 )
               b % swpDistVecsK(:,i,j,pk) = b % swpDistVecsK(:,i,j,ck)
            end do
         end do
         !! 
         !!   Cell Faces
         !!
         do k = 1, ck
            do j = 1, cj
               do i = 1, pi
                  vec1 = b % coords(i,j+1,k+1,:) - b % coords(i,j  ,k  ,:)
                  vec2 = b % coords(i,j+1,k  ,:) - b % coords(i,j  ,k+1,:)
                  call cross(vec1,vec2,b % cellFaceVecsI(:,i,j,k)) 
                  call vec_length(b % cellFaceVecsI(:,i,j,k) &
                                 ,b % cellFaceAreasI(i,j,k) )
                  b % cellFaceVecsI(:,i,j,k)  = abs(b % cellFaceVecsI  (:,i,j,k)) &
                                              / b % cellFaceAreasI (i,j,k) 
               end do
            end do
         end do
         do k = 1, ck
            do j = 1, pj
               do i = 1, ci
                  vec1 = b % coords(i+1,j,k+1,:) - b % coords(i  ,j,k  ,:)
                  vec2 = b % coords(i  ,j,k+1,:) - b % coords(i+1,j,k  ,:)
                  call cross(vec1,vec2,b % cellFaceVecsJ(:,i,j,k)) 
                  call vec_length(b % cellFaceVecsJ(:,i,j,k) &
                                 ,b % cellFaceAreasJ(i,j,k) )
                  b % cellFaceVecsJ(:,i,j,k)  = abs(b % cellFaceVecsJ  (:,i,j,k)) &
                                              / b % cellFaceAreasJ (i,j,k) 
               end do
            end do
         end do
         do k = 1, pk
            do j = 1, cj
               do i = 1, ci
                  vec1 = b % coords(i+1,j+1,k,:) - b % coords(i  ,j  ,k,:)
                  vec2 = b % coords(i+1,j  ,k,:) - b % coords(i  ,j+1,k,:)
                  call cross(vec1,vec2,b % cellFaceVecsK(:,i,j,k)) 
                  call vec_length(b % cellFaceVecsK(:,i,j,k) &
                                 ,b % cellFaceAreasK(i,j,k) )
                  b % cellFaceVecsK(:,i,j,k)  = abs(b % cellFaceVecsK  (:,i,j,k)) &
                                              / b % cellFaceAreasK (i,j,k) 
               end do
            end do
         end do
      end associate
      deallocate (swps)
      deallocate (dnw )
      deallocate (dns )
      deallocate (dnb )
   end do
   contains
      subroutine cross(a, b, vec_out)
         implicit none
         real(REAL_KIND), dimension(dimen), intent(out) :: vec_out 
         real(REAL_KIND), dimension(dimen), intent(in)  :: a, b
         vec_out(1) = a(2) * b(3) - a(3) * b(2)
         vec_out(2) = a(3) * b(1) - a(1) * b(3)
         vec_out(3) = a(1) * b(2) - a(2) * b(1)
      end subroutine cross
   
      subroutine vec_length(a, length)
         implicit none
         real(REAL_KIND),                   intent(out) :: length  
         real(REAL_KIND), dimension(dimen), intent(in)  :: a
         length = sqrt(a(1) * a(1) + a(2) * a(2) + a(3) * a(3) )
      end subroutine vec_length
      subroutine abs_cross(a, b, vec)
         implicit none
         real(REAL_KIND), dimension(3), intent(inout) :: vec
         real(REAL_KIND), dimension(3), intent(in)  :: a, b
         vec(1) = vec(1) + abs(a(2) * b(3) - a(3) * b(2))
         vec(2) = vec(2) + abs(a(3) * b(1) - a(1) * b(3))
         vec(3) = vec(3) + abs(a(1) * b(2) - a(2) * b(1))
      end subroutine abs_cross
   end subroutine
   subroutine init_sol()
   implicit none
   integer i,j,k,b
   real(REAL_KIND) :: rho
   do b = 1, nBlock
      do k = 1, blocks(b) % nCells(3)
         do j = 1, blocks(b) % nCells(2)
            do i = 1, blocks(b) % nCells(1)
               rho = blocks(b) % vars (i,j,k,1)
               blocks(b) % cons_vars (i,j,k,1) = rho
               blocks(b) % cons_vars (i,j,k,2:4) = rho &
                     * blocks(b) % vars(i,j,k,2:4)
               blocks(b) % cons_vars (i,j,k,5) = blocks(b) % vars (i,j,k,5)
            end do
         end do
      end do
   end do
   end subroutine
end module data_mod
