      SUBROUTINE INITCO (Q,QC,P,T,IB1,JB1,KB1,BN,HW,QLI,SUY2,WLIM, &
                         IV,APDF,X,Y,Z,XCC,YCC,ZCC,SLIM)
!----------------------------------------------------------------------*
!     INITIAL CONDITIONS FOR - DENSITY       RHO      : Q(I,J,1)       *
!                            - SPEED         U        : Q(I,J,2)       *
!                            - SPEED         V        : Q(I,J,3)       *
!                            - RHO * TOTAL                             *
!                              ENERGY        RHO*E_tot: Q(I,J,IEN)     *
!                            - PRESSURE      P        : P(I,J)         *
!                            - TEMPERATURE   T        : T(I,J)         *
!----------------------------------------------------------------------*

 USE CALPRO
 USE MPIVAR
 USE PROPM,ONLY: MW,F1SST,PKT,BF
 USE REFVAL
 USE SDATA
 USE SYLIMI
 USE TURBC
 USE CONTRI
 USE BLOCKE, ONLY : NB
 USE SUBSONIC_BC
#ifdef USE_MPI
 USE UFOUT, ONLY: MPIIO_BLOCK_INFO
#endif
#ifdef USE_STATISTIC
 USE STATISTIC
#endif

 IMPLICIT NONE

!!------------------------------------!!
!!------------PARAMETER---------------!!
!!------------------------------------!!

 INTEGER :: IB1
 INTEGER :: JB1
 INTEGER :: KB1
 INTEGER :: BN
 INTEGER :: IV
 INTEGER :: APDF
 REAL(KIND = 8) :: Q(IB1,JB1,KB1,IV+1)
 REAL(KIND = 8) :: QC(IB1,JB1,KB1,IV+1)
 REAL(KIND = 8) :: P(IB1,JB1,KB1)
 REAL(KIND = 8) :: T(IB1,JB1,KB1)
 REAL(KIND = 8) :: HW(JB1,KB1,IV-IFT)
 REAL(KIND = 8) :: QLI(IB1,JB1,KB1,IV+1)
 REAL(KIND = 8) :: SUY2(IB1,JB1,KB1)
 REAL(KIND = 8) :: WINI(IB1,JB1,KB1)
 REAL(KIND = 8) :: X(IB1+1,JB1+1,KB1+1)
 REAL(KIND = 8) :: Y(IB1+1,JB1+1,KB1+1)
 REAL(KIND = 8) :: Z(IB1+1,JB1+1,KB1+1)
 REAL(KIND = 8) :: XCC(IB1,JB1,KB1)
 REAL(KIND = 8) :: YCC(IB1,JB1,KB1)
 REAL(KIND = 8) :: ZCC(IB1,JB1,KB1)
 REAL(KIND = 8) :: SLIM(IB1,JB1,KB1)
 REAL(KIND = 8) :: WLIM(IB1,JB1,KB1)
!!------------------------------------!!
!!---------LOCAL VARIABLES------------!!
!!------------------------------------!!
 INTEGER :: I
 INTEGER :: ISP
 INTEGER :: J
 INTEGER :: K
 INTEGER :: L
 INTEGER :: IG1
 REAL(KIND = 8) :: EPSI
 REAL(KIND = 8) :: KI
 REAL(KIND = 8) :: TI
 REAL(KIND = 8) :: MI
 REAL(KIND = 8) :: MWMIN
 REAL(KIND = 8) :: GASC
 REAL(KIND = 8) :: DUMMY(IB1,JB1,KB1)
 REAL(KIND = 8) :: F1
 REAL(KIND = 8) :: F12,TUS,TUT


 INTEGER :: NJ,NK,NB_IN,NVAR_IN
 
 REAL(KIND = 8),allocatable :: initial_vel(:,:,:,:)
IG1 = IB1*JB1*KB1

#ifdef COMMENTS
 WRITE(*,*) 'SUBROUTINE INITCO'
#endif

!     ===============================
!     array allocation (local values)
!     ===============================
!     setzen start

#ifdef USE_STATISTIC
!     *********************************
!     initialization of Statistic-Start
!     *********************************
 QAMW  = 0.D+00
 TAMW  = 0.D+00
 PAMW  = 0.D+00
 UUAMW = 0.D+00
 VVAMW = 0.D+00
 WWAMW = 0.D+00
 UVAMW = 0.D+00
 UWAMW = 0.D+00
 VWAMW = 0.D+00
 MYLAMW = 0.D+00

 QSTD = 0.D+00
 TSTD = 0.D+00
 PSTD = 0.D+00
 UVSTD = 0.D+00
 UWSTD = 0.D+00
 VWSTD = 0.D+00
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!      SET BOUNDARY CONDITONS
!!
!!      DIESER BEREICH KOMMT IN DIE INPUT.DAT DATEI
!!      MOMENTAN NOCH IN DER INITCO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (.NOT.allocated(inflow_group) ) THEN
   Inflow_Subsonic         = .FALSE. !
   Inflow_NonRefl          = .FALSE. ! .TRUE. !
   Inflow_MassFlow_Control = .FALSE. !

   Outflow_Subsonic        = .FALSE. ! .TRUE.  !
   Outflow_NonRefl         = .FALSE. ! .TRUE.  !

   DICHTE_FLUSS_WHERE      = -1

   PRESSURE_RECONSTRUCTION = 1 !2

   n_MF_BC                 = 1

   ALLOCATE (inflow_group(n_MF_BC))

      POUT = 18.D+0 * 1.0D+5
   !< VORDEFINIERTER AUSGANGSDRUCK für UNTERSCHALL AUSSTRÖMUNG

   inflow_group(1) % MF_soll = 45.0D-3 ! 11.15D-3
   !< Vordefinierter Massenstrom durch den Einlass 1

   SUB_BC_BLOCK(1) % INFLOW_NUMBER = 1
   !< Block 1 gehört zu Einlass 1

   BC_lx  = 2.0D-1
   !< Länge des Rechengebiets

END IF

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! SET IF INFLOW PARAMETER ARE READ FROM FILE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   KI   = 1.D-4
   EPSI = 100.D+0
   KMIN = KI   * 0.005D+0
   EMIN = EPSI * 0.100D+0
   SLIM = 1.D+00
   
   OPEN(10,file='input/inco.dat')
!   MPIIO_BLOCK_INFO(10) !I
!   MPIIO_BLOCK_INFO(11) !J
!   MPIIO_BLOCK_INFO(12) !K
#ifdef USE_MPI
   allocate(initial_vel(MPIIO_BLOCK_INFO(4),MPIIO_BLOCK_INFO(5),MPIIO_BLOCK_INFO(6),3))
!  write(*,*) RANK,MPIIO_BLOCK_INFO(4:6), MPIIO_BLOCK_INFO(10),MPIIO_BLOCK_INFO(11),MPIIO_BLOCK_INFO(12) !J
   DO K = 1, MPIIO_BLOCK_INFO(6)
      DO J = 1,MPIIO_BLOCK_INFO(5)
         DO I = 1,MPIIO_BLOCK_INFO(4)
#else
   allocate(initial_vel(IB1,JB1,KB1,3))
   DO K = 1, KB1 
      DO J = 1, JB1
         DO I = 1, IB1
#endif 
            READ(10,*) initial_vel(i,j,k,:)
         END DO
      END DO
   END DO
   CLOSE(10)
   DO K = 1,KB1
      DO J = 1,JB1
         DO I = 1,IB1
#ifdef USE_MPI
            Q(I,J,K,2)   = initial_vel(i+MPIIO_BLOCK_INFO(10),j+MPIIO_BLOCK_INFO(11),k+MPIIO_BLOCK_INFO(12),1) * 60.0D0
            Q(I,J,K,3)   = initial_vel(i+MPIIO_BLOCK_INFO(10),j+MPIIO_BLOCK_INFO(11),k+MPIIO_BLOCK_INFO(12),2) * 60.0D0
            Q(I,J,K,ILV) = initial_vel(i+MPIIO_BLOCK_INFO(10),j+MPIIO_BLOCK_INFO(11),k+MPIIO_BLOCK_INFO(12),3) * 60.0D0
#else
            Q(I,J,K,2)   = initial_vel(i,j,k,1) * 60.0D0
            Q(I,J,K,3)   = initial_vel(i,j,k,2) * 60.0D0
            Q(I,J,K,ILV) = initial_vel(i,j,k,3) * 60.0D0
#endif 
            Q(I,J,K,IFT) = 1.D-20
            Q(I,J,K,ILT) = 1.D+6
            Q(I,J,K,IFS) = 0.233600D+00 
            Q(I,J,K,ILS) = 0.766400D+00
            P(I,J,K)     = 1.D+5
            T(I,J,K)     = 300.D-00
         END DO
      END DO
   END DO
   deallocate(initial_vel)
!     ============================
!     initialisation k-omega-model
!     ============================
 IF (TMOD.EQ.1) THEN
    DO K = 1,KB1
       DO J = 1,JB1
          DO I = 1,IB1
             WINI(I,J,K)   = Q(I,J,K,ILT)
          END DO
       END DO
    END DO
 END IF

!     ==============================================
!     reference values for limiters (maximum values)
!     normalization values used in AUSM
!     ==============================================
 QREF = 0.2D+0
 QREF(1)   = 2.0D+1
 QREF(2)   = 150.D+0
 QREF(3)   = 0.D+0
 IF (AXSYM.EQ.2) THEN
   QREF(ILV) = 150.D+0
 ENDIF
 QREF(IEN) = 20000.D+0
 QREF(IFT) = 1000.D+0
 QREF(ILT) = 1.0D+4
 IF(APDF.GE.1)THEN
   QREF(IFP)   = 50.
   QREF(ILP)   = 1.D-2
 END IF
 ATREF(1)  = 340.D+0
 ASREF(1)  = 340.D+0
 HTREF(1)  = 20000.D+0
 PREF(1)   = 20000.D+0
 EKREF(1)  = 0.04D+0
 IF (AXSYM.EQ.2) THEN
    QREF(ILV) = 150.D+0
 ENDIF
 ATREF(1)  = 340.D+0
 ASREF(1)  = 340.D+0
 HTREF(1)  = 20000.D+0
 PREF(1)   = 5.2D+6
 EKREF(1)  = 1.4D+0

!     ==============================================================
!     definition for data output in case of time accurate simulation
!     ==============================================================
 IF ((IBOUT.GT.0).AND.(TAI.EQ.1)) THEN
    DO I = 1,IBOUT
       IDOUT(I) = 0.0005D+0 * DBLE(I)
    END DO
 END IF
!     setzen ende

!     =======================
!     deallocate local arrays
!     =======================
 RETURN
      END SUBROUTINE INITCO
       SUBROUTINE WALLDIST (IB1,JB1,KB1,IB2,JB2,KB2,X,Y,Z,XCC,YCC,ZCC,YZ,BN)
!----------------------------------------------------------------------!
!        BERECHUNG DES WANDABSTANDES BEI PROGRAMMINITIALISIERUNG.
!----------------------------------------------------------------------!
      USE MPIVAR

      IMPLICIT NONE

!%$ ------------------------------ $%!
!%$ --------PARAMETER------------- $%!
!%$ ------------------------------ $%!
      INTEGER :: IB1
      INTEGER :: JB1
      INTEGER :: KB1
      INTEGER :: IB2
      INTEGER :: JB2
      INTEGER :: KB2
      INTEGER :: BN
      REAL(KIND = 8) :: X(IB2,JB2,KB2)
      REAL(KIND = 8) :: Y(IB2,JB2,KB2)
      REAL(KIND = 8) :: Z(IB2,JB2,KB2)
      REAL(KIND = 8) :: XCC(IB1,JB1,KB1)
      REAL(KIND = 8) :: YCC(IB1,JB1,KB1)
      REAL(KIND = 8) :: ZCC(IB1,JB1,KB1)
      REAL(KIND = 8) :: YZ(IB1,JB1,KB1)
!%$ ------------------------------ $%!
!%$ ------LOKALE VARIABLEN-------- $%!
!%$ ------------------------------ $%!
      INTEGER :: I
      INTEGER :: J
      INTEGER :: K
      REAL(KIND = 8) :: ZP
      REAL(KIND = 8) :: YP
      REAL(KIND = 8) :: XP
      REAL(KIND = 8) :: YO
      REAL(KIND = 8) :: RADIUS,YWAND
!     ========================
!     distance to nearest wall
!     ========================
!     setzen start
 DO K=1,KB1
    DO J = 1,JB1
        DO I = 1,IB1
           YZ(I,J,K) = 1.0D20
        END DO
    END DO
 END DO
!     setzen ende
      END SUBROUTINE WALLDIST
