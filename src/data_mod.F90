module data_mod
use const_mod
implicit none
type :: boundary_type
   integer :: bc_type
end type boundary_type
type :: block_type
   type(boundary_type)           :: boundary(6)
   integer(INT_KIND)             :: nPkts(3)

   integer(INT_KIND)             :: nCells(3)

   real(REAL_KIND), allocatable  :: coords         (:,:,:,:)
   !< Koordinaten der Gitterpunkte (nPkts;dimen)
   real(REAL_KIND), allocatable  :: vars           (:,:,:,:)
   !< primitiver Variablenvektor in den Zellmittelpunkten (nCells, nVar)
   real(REAL_KIND), allocatable  :: cons_vars      (:,:,:,:)
   !< Konservativer Variablenvektor in den Zellmittelpunkten (nCells, nVar)
   !< wird in der Routine UPDATE_SOL in TIME_DISC_MOD verwendet 
   real(REAL_KIND), allocatable  :: temperatures   (:,:,:)
   !< Temperaturen in den Zellmittelpunkten (nCells)
   real(REAL_KIND), allocatable  :: pressures      (:,:,:)
   !< Druck in den Zellmittelpunkten (nCells)
   real(REAL_KIND), allocatable  :: viscosities    (:,:,:)
   !< viskositaet im Zellmittelpunkt
   real(REAL_KIND), allocatable  :: heatKoeffs     (:,:,:)
   !< Waermeleitungskoeffizient
   real(REAL_KIND), allocatable  :: residuals      (:,:,:,:)
   !< Residuum

   real(REAL_KIND), allocatable  :: volumes        (:,:,:)
   !< Inerses Zellvolumen der jeweiligen Zelle
   real(REAL_KIND), allocatable  :: faceVarsLeftI  (:,:,:,:)
   !< Variable auf der Zellfläche in I-Richtung von der Linken Seite (i-1)

   real(REAL_KIND), allocatable  :: faceVarsRightI (:,:,:,:)
   !< Variable auf der Zellfläche in I-Richtung von der Rechte nSeite (i)

   real(REAL_KIND), allocatable  :: faceVarsLeftJ  (:,:,:,:)
   !< Variable auf der Zellfläche in J-Richtung von der Linken Seite (j-1)

   real(REAL_KIND), allocatable  :: faceVarsRightJ(:,:,:,:)
   !< Variable auf der Zellfläche in J-Richtung von der Rechte nSeite (j)

   real(REAL_KIND), allocatable  :: faceVarsLeftK  (:,:,:,:)
   !< Variable auf der Zellfläche in K-Richtung von der Linken Seite (k-1)

   real(REAL_KIND), allocatable  :: faceVarsRightK(:,:,:,:)
   !< Variable auf der Zellfläche in K-Richtung von der Rechte nSeite (k)

   real(REAL_KIND), allocatable  :: fluxesI        (:,:,:,:)
   real(REAL_KIND), allocatable  :: fluxesJ        (:,:,:,:)
   real(REAL_KIND), allocatable  :: fluxesK        (:,:,:,:)

   real(REAL_KIND), allocatable  :: visFluxesI     (:,:,:,:)
   real(REAL_KIND), allocatable  :: visFluxesJ     (:,:,:,:)
   real(REAL_KIND), allocatable  :: visFluxesK     (:,:,:,:)

   real(REAL_KIND), allocatable  :: tauI           (:,:,:,:)
   real(REAL_KIND), allocatable  :: tauJ           (:,:,:,:)
   real(REAL_KIND), allocatable  :: tauK           (:,:,:,:)

   real(REAL_KIND), allocatable  :: heatfluxI      (:,:,:,:)
   real(REAL_KIND), allocatable  :: heatfluxJ      (:,:,:,:)
   real(REAL_KIND), allocatable  :: heatfluxK      (:,:,:,:)

   !> DuDn at the faces
   real(REAL_KIND), allocatable  :: dUdnI          (:,:,:,:,:)
   real(REAL_KIND), allocatable  :: dUdnJ          (:,:,:,:,:)
   real(REAL_KIND), allocatable  :: dUdnK          (:,:,:,:,:)

   real(REAL_KIND), allocatable  :: dUdn           (:,:,:,:,:)
   !<DuDn at the Cell center, for turbulence Model
   !< Dimension: (i,j,k,Komponente,Space-Dir)

   real(REAL_KIND), allocatable  :: cellFaceAreasI (:,:,:)
   !< Fläche der Zellflächen(3D) & Zellkanten (2D) in I Richtung (EAST & WEST)
   !< NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable  :: cellFaceAreasJ (:,:,:)
   !< Fläche der Zellflächen(3D) & Zellkanten (2D) in J Richtung (NORTH & SOUTH)
   !< NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable  :: cellFaceAreasK (:,:,:)
   !< Fläche der Zellflächen(3D) & Zellkanten (2D) in K Richtung (FROMT & BACK)
   !< NCell(1), NCell(2), NPkt(3)

   real(REAL_KIND), allocatable  :: cellFaceVecsI (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in I Richtung (EAST & WEST)
   !< Dimen, NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable  :: cellFaceVecsJ (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in J Richtung (NORTH & SOUTH)
   !< Dimen, NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable  :: cellFaceVecsK (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in K Richtung (FROMT & BACK)
   !< Dimen, NPkt(1), NCell(2), NPkt(3)

   real(REAL_KIND), allocatable  :: abscellFaceVecsI (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in I Richtung (EAST & WEST)
   !< Dimen, NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable  :: abscellFaceVecsJ (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in J Richtung (NORTH & SOUTH)
   !< Dimen, NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable  :: abscellFaceVecsK (:,:,:,:)
   !< normierter Normalenvektor der Zellflächen(3D) & Zellkanten (2D) in K Richtung (FROMT & BACK)
   !< Dimen, NPkt(1), NCell(2), NPkt(3)

   real(REAL_KIND), allocatable  :: swpDistVecsI (:,:,:,:)
   !< Vektor mit Schwerpunkt-Abständen für die Berechnung der Ableitungen in I Richtung (EAST & WEST)
   !< Dimen, NPkt(1), NCell(2), NCell(3)

   real(REAL_KIND), allocatable  :: swpDistVecsJ (:,:,:,:)
   !< Vektor mit Schwerpunkt-Abständen für die Berechnung der Ableitungen in J Richtung (NORTH & SOUTH)
   !< Dimen, NCell(1), NPkt(2), NCell(3)

   real(REAL_KIND), allocatable  :: swpDistVecsK (:,:,:,:)
   !< Vektor mit Schwerpunkt-Abständen für die Berechnung der Ableitungen in K Richtung (FROMT & BACK)
   !< Dimen, NPkt(1), NCell(2), NPkt(3)

   real(REAL_KIND), allocatable :: dnw(:,:,:,:)
   !< Normalenvektoren der Zellseite Westseite
   real(REAL_KIND), allocatable :: dns(:,:,:,:)
   !< Normalenvektoren der Zellseite Suedseite
   real(REAL_KIND), allocatable :: dnb(:,:,:,:)
   !< Normalenvektoren der Zellseite Backseite

   real(REAL_KIND), allocatable :: cellsizes (:,:,:,:) 
   !< Zelllaengen
   !< Dimensionen: (i,j,k,Richtung)
   real(REAL_KIND), allocatable :: lles(:,:,:)
end type block_type

type(block_type), allocatable    :: blocks      (:)
integer(INT_KIND)                :: nBlock
integer(INT_KIND)                :: dimen
integer(INT_KIND)                :: nFaces
integer(INT_KIND)                :: nCorners
integer(INT_KIND)                :: nCell
integer(INT_KIND)                :: nVar

contains
   subroutine allocate_vars(ib)
      use control_mod, only: nBoundaryCells
   implicit none

   integer, intent(in) :: ib
   !< Blockindex
   integer :: ci,cj,ck
   ! Number of Cells
   integer :: pi,pj,pk
   ! NUmber of Grid Points
   integer :: cis,cjs,cks
   ! Cellsindex Start including boundary cells
   integer :: cie,cje,cke
   ! Cellsindex End including boundary cells
   
      nvar = 5
!      block_loop: do ib = 1, nBlock
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

         allocate(b % coords            (       pi     ,pj     ,pk     ,dimen ) )
         allocate(b % vars              (       cis:cie,cjs:cje,cks:cke,nvar  ) )
         allocate(b % temperatures      (       cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % pressures         (       cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % viscosities       (       cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % heatKoeffs        (       cis:cie,cjs:cje,cks:cke       ) )
         allocate(b % cons_vars         (       ci     ,cj     ,ck     ,nvar  ) )
         allocate(b % residuals         (       ci     ,cj     ,ck     ,nvar  ) )

         allocate(b % faceVarsLeftI     (       pi     ,cj     ,ck     ,nvar  ) )
         allocate(b % faceVarsRightI    (       pi     ,cj     ,ck     ,nvar  ) )

         allocate(b % faceVarsLeftJ     (       ci     ,pj     ,ck     ,nvar  ) )
         allocate(b % faceVarsRightJ    (       ci     ,pj     ,ck     ,nvar  ) )

         allocate(b % faceVarsLeftK     (       ci     ,cj     ,pk     ,nvar  ) )
         allocate(b % faceVarsRightK    (       ci     ,cj     ,pk     ,nvar  ) )

         allocate(b % fluxesI           (       pi     ,cj     ,ck     ,nvar  ) )
         allocate(b % fluxesJ           (       ci     ,pj     ,ck     ,nvar  ) )
         allocate(b % fluxesK           (       ci     ,cj     ,pk     ,nvar  ) )

         allocate(b % visFluxesI        (       pi     ,cj     ,ck     ,nvar  ) )
         allocate(b % visFluxesJ        (       ci     ,pj     ,ck     ,nvar  ) )
         allocate(b % visFluxesK        (       ci     ,cj     ,pk     ,nvar  ) )

         allocate(b % tauI              (       pi     ,cj     ,ck     ,6     ) )
         allocate(b % tauJ              (       ci     ,pj     ,ck     ,6     ) )
         allocate(b % tauK              (       ci     ,cj     ,pk     ,6     ) )

         allocate(b % heatfluxI         (       pi     ,cj     ,ck     ,dimen ) )
         allocate(b % heatfluxJ         (       ci     ,pj     ,ck     ,dimen ) )
         allocate(b % heatfluxK         (       ci     ,cj     ,pk     ,dimen ) )

         allocate(b % volumes           (       ci     ,cj     ,ck            ) )
         allocate(b % CellFaceAreasI    (       pi     ,cj     ,ck            ) )
         allocate(b % CellFaceAreasJ    (       ci     ,pj     ,ck            ) )
         allocate(b % CellFaceAreasK    (       ci     ,cj     ,pk            ) )
         allocate(b % CellFaceVecsI     ( dimen,pi     ,cj     ,ck            ) )
         allocate(b % CellFaceVecsJ     ( dimen,ci     ,pj     ,ck            ) )
         allocate(b % CellFaceVecsK     ( dimen,ci     ,cj     ,pk            ) )
         allocate(b % absCellFaceVecsI  ( dimen,pi     ,cj     ,ck            ) )
         allocate(b % absCellFaceVecsJ  ( dimen,ci     ,pj     ,ck            ) )
         allocate(b % absCellFaceVecsK  ( dimen,ci     ,cj     ,pk            ) )
   
         allocate(b % swpDistVecsI      ( dimen,pi     ,cj     ,ck            ) )
         allocate(b % swpDistVecsJ      ( dimen,ci     ,pj     ,ck            ) )
         allocate(b % swpDistVecsK      ( dimen,ci     ,cj     ,pk            ) )
   
         allocate(b % dUdnI             (       pi     ,cj     ,ck     ,dimen+1,dimen) )
         allocate(b % dUdnJ             (       ci     ,pj     ,ck     ,dimen+1,dimen) )
         allocate(b % dUdnK             (       ci     ,cj     ,pk     ,dimen+1,dimen) )

         allocate(b % dUdn              (       ci     ,cj     ,ck     ,dimen  ,dimen) )
         allocate(b % dnw               ( dimen,pi     ,cj     ,ck            ) )
         allocate(b % dns               ( dimen,ci     ,pj     ,ck            ) )
         allocate(b % dnb               ( dimen,ci     ,cj     ,pk            ) )

         allocate(b % cellsizes         (       ci     ,cj     ,ck     ,dimen ) ) 

         allocate(b % lles              (       ci     ,cj     ,ck            ) ) 

         b % visFluxesI = 0.0E0_REAL_KIND
         b % visFluxesJ = 0.0E0_REAL_KIND
         b % visFluxesK = 0.0E0_REAL_KIND
      end associate
   
!      end do block_loop


   end subroutine allocate_vars
   subroutine calc_grid()
      use control_mod, only: c_les_sgs

   implicit none
      integer     :: ib, i,j,k
      real(REAL_KIND), allocatable :: swps(:,:,:,:)
      real(REAL_KIND), dimension(dimen) :: vec1, vec2
      real(REAL_KIND) :: len
      do ib = 1,nBlock
      associate (b  => blocks(ib) ,&
                 ci => blocks(ib) % nCells(1) ,&
                 cj => blocks(ib) % nCells(2) ,&
                 ck => blocks(ib) % nCells(3) ,&
                 pi => blocks(ib) % nPkts(1) ,&
                 pj => blocks(ib) % nPkts(2) ,&
                 pk => blocks(ib) % nPkts(3) )
         allocate (swps(dimen,ci,cj,ck))
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
                  b % dnw(1,I,J,K) =  0.5D+0 * ( (b % coords(I  ,J  ,K  ,2) - b % coords(I  ,J+1,K+1,2)) &
                                             * (  b % coords(I  ,J  ,K+1,3) - b % coords(I  ,J+1,K  ,3)) &
                                             - (  b % coords(I  ,J  ,K  ,3) - b % coords(I  ,J+1,K+1,3)) &
                                             * (  b % coords(I  ,J  ,K+1,2) - b % coords(I  ,J+1,K  ,2)) )
                  b % dnw(2,I,J,K) = -0.5D+0 * ( (b % coords(I  ,J  ,K  ,3) - b % coords(I  ,J+1,K+1,3)) &
                                             * (  b % coords(I  ,J  ,K+1,1) - b % coords(I  ,J+1,K  ,1)) &
                                             - (  b % coords(I  ,J  ,K  ,1) - b % coords(I  ,J+1,K+1,1)) &
                                             * (  b % coords(I  ,J  ,K+1,3) - b % coords(I  ,J+1,K  ,3)) )
                  b % dnw(3,I,J,K) =  0.5D+0 * ( (b % coords(I  ,J  ,K  ,1) - b % coords(I  ,J+1,K+1,1)) &
                                             * (  b % coords(I  ,J  ,K+1,2) - b % coords(I  ,J+1,K  ,2)) &
                                             - (  b % coords(I  ,J  ,K  ,2) - b % coords(I  ,J+1,K+1,2)) &
                                             * (  b % coords(I  ,J  ,K+1,1) - b % coords(I  ,J+1,K  ,1)) )
                  if (i == 30 .and. j == 30 .and. k == 30 ) then
                  write(*,*) "======", i,j,k, "======"
                  write(*,*) "dnw1", b % dnw(1,i,j,k)
                  write(*,*) b % coords(I  ,J  ,K  ,2) , b % coords(I  ,J+1,K+1,2)
                  write(*,*) b % coords(I  ,J  ,K+1,3) , b % coords(I  ,J+1,K  ,3)
                  write(*,*) b % coords(I  ,J  ,K  ,3) , b % coords(I  ,J+1,K+1,3)
                  write(*,*) b % coords(I  ,J  ,K+1,2) , b % coords(I  ,J+1,K  ,2)
                  write(*,*) "dnw2", b % dnw(2,i,j,k)
                  write(*,*) b % coords(I  ,J  ,K  ,3) , b % coords(I  ,J+1,K+1,3)
                  write(*,*) b % coords(I  ,J  ,K+1,1) , b % coords(I  ,J+1,K  ,1)
                  write(*,*) b % coords(I  ,J  ,K  ,1) , b % coords(I  ,J+1,K+1,1)
                  write(*,*) b % coords(I  ,J  ,K+1,3) , b % coords(I  ,J+1,K  ,3)
                  write(*,*) "dnw3", b % dnw(3,i,j,k)
                  write(*,*) b % coords(I  ,J  ,K  ,1) , b % coords(I  ,J+1,K+1,1)
                  write(*,*) b % coords(I  ,J  ,K+1,2) , b % coords(I  ,J+1,K  ,2)
                  write(*,*) b % coords(I  ,J  ,K  ,2) , b % coords(I  ,J+1,K+1,2)
                  write(*,*) b % coords(I  ,J  ,K+1,1) , b % coords(I  ,J+1,K  ,1)
               !stop
            end if
               end do
            end do
         end do
     
         DO K = 1,ck 
            DO J = 1,pj 
               DO I = 1,ci 
                  b % dns(1,I,J,K) =  0.5D+0 * ( (b % coords(I  ,J,K,2)-b % coords(I+1,J,K+1,2)) &
                                             * (  b % coords(I+1,J,K,3)-b % coords(I  ,J,K+1,3)) &
                                             - (  b % coords(I  ,J,K,3)-b % coords(I+1,J,K+1,3)) &
                                             * (  b % coords(I+1,J,K,2)-b % coords(I  ,J,K+1,2)) )
                  b % dns(2,I,J,K) = -0.5D+0 * ( (b % coords(I  ,J,K,3)-b % coords(I+1,J,K+1,3)) &
                                             * (  b % coords(I+1,J,K,1)-b % coords(I  ,J,K+1,1)) &
                                             - (  b % coords(I  ,J,K,1)-b % coords(I+1,J,K+1,1)) &
                                             * (  b % coords(I+1,J,K,3)-b % coords(I  ,J,K+1,3)) )
                  b % dns(3,I,J,K) =  0.5D+0 * ( (b % coords(I  ,J,K,1)-b % coords(I+1,J,K+1,1)) &
                                             * (  b % coords(I+1,J,K,2)-b % coords(I  ,J,K+1,2)) &
                                             - (  b % coords(I  ,J,K,2)-b % coords(I+1,J,K+1,2)) &
                                             * (  b % coords(I+1,J,K,1)-b % coords(I  ,J,K+1,1)) )
                  if (i == 30 .and. j == 30 .and. k == 30 ) then
                  write(*,*) "======", i,j,k, "======"
                  write(*,*) "dns1", b % dns(1,I,J,K)
                  write(*,*) b % coords(I,J,K,2),b % coords(I+1,J,K+1,2)
                  write(*,*) b % coords(I+1,J,K,3),b % coords(I,J,K+1,3)
                  write(*,*) b % coords(I,J,K,3),b % coords(I+1,J,K+1,3)
                  write(*,*) b % coords(I+1,J,K,2),b % coords(I,J,K+1,2)
                  write(*,*) "dns2", b % dns(2,I,J,K) 
                  write(*,*) b % coords(I,J,K,3),b % coords(I+1,J,K+1,3)
                  write(*,*) b % coords(I+1,J,K,1),b % coords(I,J,K+1,1)
                  write(*,*) b % coords(I,J,K,1),b % coords(I+1,J,K+1,1)
                  write(*,*) b % coords(I+1,J,K,3),b % coords(I,J,K+1,3)
                  write(*,*) "dns3", b % dns(3,I,J,K)
                  write(*,*) b % coords(I,J,K,1),b % coords(I+1,J,K+1,1)
                  write(*,*) b % coords(I+1,J,K,2),b % coords(I,J,K+1,2)
                  write(*,*) b % coords(I,J,K,2),b % coords(I+1,J,K+1,2)
                  write(*,*) b % coords(I+1,J,K,1),b % coords(I,J,K+1,1)
               !stop
            end if
               end do
            end do
         end do
     
         DO K = 1,pk 
            DO J = 1,cj  
               DO I = 1,ci 
                  b % dnb(1,I,J,K) =  0.5D+0 * ( (b % coords(I,J  ,K,2)-b % coords(I+1,J+1,K,2)) &
                                             * (  b % coords(I,J+1,K,3)-b % coords(I+1,J  ,K,3)) &
                                             - (  b % coords(I,J  ,K,3)-b % coords(I+1,J+1,K,3)) &
                                             * (  b % coords(I,J+1,K,2)-b % coords(I+1,J  ,K,2)) )
                  b % dnb(2,I,J,K) = -0.5D+0 * ( (b % coords(I,J  ,K,3)-b % coords(I+1,J+1,K,3)) &
                                             * (  b % coords(I,J+1,K,1)-b % coords(I+1,J  ,K,1)) &
                                             - (  b % coords(I,J  ,K,1)-b % coords(I+1,J+1,K,1)) &
                                             * (  b % coords(I,J+1,K,3)-b % coords(I+1,J  ,K,3)) )
                  b % dnb(3,I,J,K) =  0.5D+0 * ( (b % coords(I,J  ,K,1)-b % coords(I+1,J+1,K,1)) &
                                             * (  b % coords(I,J+1,K,2)-b % coords(I+1,J  ,K,2)) &
                                             - (  b % coords(I,J  ,K,2)-b % coords(I+1,J+1,K,2)) &
                                             * (  b % coords(I,J+1,K,1)-b % coords(I+1,J  ,K,1)) )
                  if (i == 30 .and. j == 30 .and. k == 30 ) then
                  write(*,*) "======", i,j,k, "======"
                  write(*,*)  "dnb1", b % dnb(1,I,J,K)
                  write(*,*)  b % coords(I,J  ,K,2),b % coords(I+1,J+1,K,2)
                  write(*,*)  b % coords(I,J+1,K,3),b % coords(I+1,J  ,K,3)
                  write(*,*)  b % coords(I,J  ,K,3),b % coords(I+1,J+1,K,3)
                  write(*,*)  b % coords(I,J+1,K,2),b % coords(I+1,J  ,K,2)
                  write(*,*)  "dnb2", b % dnb(2,I,J,K)
                  write(*,*)  b % coords(I,J  ,K,3),b % coords(I+1,J+1,K,3)
                  write(*,*)  b % coords(I,J+1,K,1),b % coords(I+1,J  ,K,1)
                  write(*,*)  b % coords(I,J  ,K,1),b % coords(I+1,J+1,K,1)
                  write(*,*)  b % coords(I,J+1,K,3),b % coords(I+1,J  ,K,3)
                  write(*,*)  "dnb3", b % dnb(3,I,J,K)
                  write(*,*)  b % coords(I,J  ,K,1),b % coords(I+1,J+1,K,1)
                  write(*,*)  b % coords(I,J+1,K,2),b % coords(I+1,J  ,K,2)
                  write(*,*)  b % coords(I,J  ,K,2),b % coords(I+1,J+1,K,2)
                  write(*,*)  b % coords(I,J+1,K,1),b % coords(I+1,J  ,K,1)
               !stop
            end if
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
                       + vec1(1) * ( + b % DnW(2,I,J,K) + b % DnW(2,I+1,J,K) )&
                       + vec1(2) * ( - b % DnW(1,I,J,K) - b % DnW(1,I+1,J,K) )&
                       + vec1(3) * ( + b % DnW(3,I,J,K) + b % DnW(3,I+1,J,K) )
     
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
                       + vec1(1) * ( + b % dns(2,I,J,K) + b % dns(2,I,J+1,K) )&
                       + vec1(2) * ( - b % dns(1,I,J,K) - b % dns(1,I,J+1,K) )&
                       + vec1(3) * ( + b % dns(3,I,J,K) + b % dns(3,I,J+1,K) )
     
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
                       + vec1(1) * ( + b % dnb(2,I,J,K) + b % dnb(2,I,J,K+1) )&
                       + vec1(2) * ( - b % dnb(1,I,J,K) - b % dnb(1,I,J,K+1) )&
                       + vec1(3) * ( + b % dnb(3,I,J,K) + b % dnb(3,I,J,K+1) )
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
         !!  damit nur eine Multiplication mit dem Delta der Zellwerte notwendig ist, wird der Vektor
         !!  normalisiert (Anteil) und druch die Länge geteilt (Gradient)
         !!  dudx = (u[i] - u[i-1]) / len * dx_norm
         do k = 1, ck
            do j = 1, cj
               do i = 1, pi
                  call vec_length(b % swpDistVecsI(:,i,j,k),len)
                  b % swpDistVecsI(:,i,j,k) = b % swpDistVecsI(:,i,j,k) / (len * len)
               end do
            end do
         end do
         do k = 1, ck
            do j = 1, pj
               do i = 1, ci
                  call vec_length(b % swpDistVecsJ(:,i,j,k),len)
                  b % swpDistVecsJ(:,i,j,k) = b % swpDistVecsJ(:,i,j,k) / (len * len)
               end do
            end do
         end do
         do k = 1, pk
            do j = 1, cj
               do i = 1, ci
                  call vec_length(b % swpDistVecsK(:,i,j,k),len)
                  b % swpDistVecsK(:,i,j,k) = b % swpDistVecsK(:,i,j,k) / (len * len)
               end do
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
                  b % cellFaceVecsI(:,i,j,k)  = b % cellFaceVecsI  (:,i,j,k) &
                                              / b % cellFaceAreasI (i,j,k) 
                  b % abscellFaceVecsI(:,i,j,k)  = abs(b % cellFaceVecsI  (:,i,j,k))
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
                  b % cellFaceVecsJ(:,i,j,k)  = b % cellFaceVecsJ  (:,i,j,k) &
                                              / b % cellFaceAreasJ (i,j,k) 
                  b % abscellFaceVecsJ(:,i,j,k)  = abs(b % cellFaceVecsJ  (:,i,j,k))
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
                  b % cellFaceVecsK(:,i,j,k)  = b % cellFaceVecsK  (:,i,j,k) &
                                              / b % cellFaceAreasK (i,j,k) 
                  b % abscellFaceVecsK(:,i,j,k)  = abs(b % cellFaceVecsK  (:,i,j,k))
               end do
            end do
         end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!    Zelllaengen (abstand der Zellseitenschwerpunkte   !!!!!!!!!!!
!!!!  heatflux array i sused as temp (due to their dimension)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!! OST/WEST-Richtung
         do k = 1, ck
            do j = 1, cj
               do i = 1, pi
                  b % heatfluxI(i,j,k,:) = 0.25d0 *                              &
                                                  ( b % coords ( i  ,j  ,k  ,:)  &
                                                  + b % coords ( i  ,j+1,k  ,:)  &
                                                  + b % coords ( i  ,j  ,k+1,:)  &
                                                  + b % coords ( i  ,j+1,k+1,:)  )
                  
               end do
               do i = 1, ci
                  b % cellsizes(i,j,k,DIR_EASTWEST) = sqrt  &
                                                      ( (b % heatfluxI (i+1,j,k,1) - b % heatfluxI (i  ,j,k,1)) **2 &
                                                      + (b % heatfluxI (i+1,j,k,2) - b % heatfluxI (i  ,j,k,2)) **2 &
                                                      + (b % heatfluxI (i+1,j,k,3) - b % heatfluxI (i  ,j,k,3)) **2 )
               end do
            end do
         end do
         
!!!!!!!!! SUED/NORD-Richtung
         do k = 1, ck
            do j = 1, pj
               do i = 1, ci
                  b % heatfluxJ(i,j,k,:) = 0.25d0 *                              &
                                                  ( b % coords ( i  ,j  ,k  ,:)  &
                                                  + b % coords ( i+1,j  ,k  ,:)  &
                                                  + b % coords ( i  ,j  ,k+1,:)  &
                                                  + b % coords ( i+1,j  ,k+1,:)  )
                  
               end do
            end do
            do j = 1, cj
               do i = 1, ci
                  b % cellsizes(i,j,k,DIR_SOUTHNORTH) = sqrt  &
                                                      ( (b % heatfluxJ (i,j+1,k,1) - b % heatfluxJ (i,j,k,1)) **2 &
                                                      + (b % heatfluxJ (i,j+1,k,2) - b % heatfluxJ (i,j,k,2)) **2 &
                                                      + (b % heatfluxJ (i,j+1,k,3) - b % heatfluxJ (i,j,k,3)) **2 )
               end do
            end do
         end do
!!!!!!!!! FRONT/BACK-Richtung
         do k = 1, pk
            do j = 1, cj
               do i = 1, ci
                  b % heatfluxK(i,j,k,:) = 0.25d0 *                              &
                                                  ( b % coords ( i  ,j  ,k  ,:)  &
                                                  + b % coords ( i+1,j  ,k  ,:)  &
                                                  + b % coords ( i  ,j+1,k  ,:)  &
                                                  + b % coords ( i+1,j+1,k  ,:)  )
                  
               end do
            end do
         end do
         do k = 1, ck
            do j = 1, cj
               do i = 1, ci
                  b % cellsizes(i,j,k,DIR_FRONTBACK) = sqrt  &
                                                      ( (b % heatfluxK (i,j,k+1,1) - b % heatfluxK (i,j,k,1)) **2 &
                                                      + (b % heatfluxK (i,j,k+1,2) - b % heatfluxK (i,j,k,2)) **2 &
                                                      + (b % heatfluxK (i,j,k+1,3) - b % heatfluxK (i,j,k,3)) **2 )
               end do
            end do
         end do
         do k = 1, ck
            do j = 1, cj
               do i = 1, ci
                  b % lles(i,j,k) = (c_les_sgs * maxval(b % cellsizes(i,j,k,:) )) ** 2
               end do
            end do
         end do
      end associate
      deallocate (swps)
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
!      subroutine abs_cross(a, b, vec)
!         implicit none
!         real(REAL_KIND), dimension(3), intent(inout) :: vec
!         real(REAL_KIND), dimension(3), intent(in)  :: a, b
!         vec(1) = vec(1) + abs(a(2) * b(3) - a(3) * b(2))
!         vec(2) = vec(2) + abs(a(3) * b(1) - a(1) * b(3))
!         vec(3) = vec(3) + abs(a(1) * b(2) - a(2) * b(1))
!      end subroutine abs_cross
   end subroutine
   subroutine init_sol()
   implicit none
   integer i,j,k,b
   real(REAL_KIND) :: rho
   do b = 1, nBlock
      do k = 1, blocks(b) % nCells(3)
         do j = 1, blocks(b) % nCells(2)
            do i = 1, blocks(b) % nCells(1)
               rho = blocks(b) % vars (i,j,k,VEC_RHO)
               blocks(b) % cons_vars (i,j,k,VEC_RHO) = rho
               blocks(b) % cons_vars (i,j,k,VEC_SPU:VEC_SPW) = rho &
                     * blocks(b) % vars(i,j,k,VEC_SPU:VEC_SPW)
               blocks(b) % cons_vars (i,j,k,VEC_ENE) = blocks(b) % vars (i,j,k,VEC_ENE)
               blocks(b) % pressures(i,j,k) = (GAMMA - 1.0E0_REAL_KIND)             &
                                            * (                                     &
                                                blocks(b) % vars(i,j,k,VEC_ENE)     &
                                              - 0.5E0_REAL_KIND                     &
                                              * rho                                 &
                                              * ( blocks(b) % vars (i,j,k,VEC_SPU)  &
                                                * blocks(b) % vars (i,j,k,VEC_SPU)  & 
                                                + blocks(b) % vars (i,j,k,VEC_SPV)  &
                                                * blocks(b) % vars (i,j,k,VEC_SPV)  &
                                                + blocks(b) % vars (i,j,k,VEC_SPW)  &
                                                * blocks(b) % vars (i,j,k,VEC_SPW)  &
                                                )                                   &
                                              )
               blocks(b) % temperatures (i,j,k) = blocks(b) % pressures (i,j,k) & 
                     / ( RGas * rho)
            end do
         end do
      end do
   end do
   end subroutine
end module data_mod
