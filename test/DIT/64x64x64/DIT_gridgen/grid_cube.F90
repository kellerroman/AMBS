PROGRAM CUBE
!AUTHOR:       MARKUS KINDLER, ROMAN KELLER
!START DATE:   22.02.2016
!LAST CHANGES: 22.02.2016
!PURPOSE:      ERSTELLT EIN LES GITTER FUER DIE DIT
! CHANGELOG:
! 22.02.2016, RK: START AUF BASIS DER VERSION VONMARKUS KINDLER
! 22.02.2016, RK: CHANGED TO BINARY GIT FILE OUTPUT
  IMPLICIT NONE

  integer :: i,j,k
  real(kind = 8), allocatable, dimension(:,:,:) :: x,y,z 

  integer               , parameter :: AXSYM = 2
  integer               , parameter :: NUM_OF_BLOCKS = 1                                            ! no. of blocks

  integer               , parameter :: NCELL = 64
  integer               , parameter :: IM = NCELL + 1
  integer               , parameter :: JM = NCELL + 1
  integer               , parameter :: KM = NCELL + 1
  !real(kind = 8 )       , parameter :: PI = 3.14159265359D+0
  real(kind = 8 )       , parameter :: PI =3.183E-00 
  logical               , parameter :: TEC = .FALSE.

  open(10,file='git.bin',form="UNFORMATTED",access="STREAM",status="REPLACE")

  ALLOCATE(X(IM,JM,KM),Y(IM,JM,KM),Z(IM,JM,KM))

  DO K = 1,KM
     DO J = 1,JM
        DO I = 1,IM
           X(I,J,K) = REAL(I-1)/REAL(IM-1)*1.D-02*PI
           Y(I,J,K) = REAL(J-1)/REAL(JM-1)*1.D-02*PI
           Z(I,J,K) = REAL(K-1)/REAL(KM-1)*1.D-02*PI
        END DO
     END DO
  END DO
  write(10) axsym
  write(10) num_of_blocks 
  WRITE(10) IM,JM,KM
  DO K = 1,KM
     DO J = 1,JM
        DO I = 1,IM
           WRITE(10) X(I,J,K),Y(I,J,K),Z(I,J,K)
        END DO
     END DO
  END DO
  
      ! ZS.DAT SCHREIBEN

      OPEN(20,file='zs.dat')

      WRITE(20,10000) 
      WRITE(20,10001)
      WRITE(20,10000) 
      WRITE(20,*) ' '
      WRITE(20,10002) 1
      WRITE(20,*) ' '

      WRITE(20,10003) 1      
      WRITE(20,10004) IM-1    
      WRITE(20,10005) JM-1 
      WRITE(20,10006) KM-1    
      WRITE(20,10007) 1 
      WRITE(20,10008) 1
      WRITE(20,10009) 1
      WRITE(20,10010) 1
      WRITE(20,10011) 1
      WRITE(20,10012) IM-1
      WRITE(20,10013) 1
      WRITE(20,10014) 1
      WRITE(20,10015) 1
      WRITE(20,10016) IM-1
      WRITE(20,10017) 1
      WRITE(20,10018) 1   
      WRITE(20,*) ' '

10000 FORMAT('---------------------------------')
10001 FORMAT('Blockstruktur und Randbedingungen')
10002 FORMAT(('number of blocks on CPU .............  NB  ='),2X,I3) 
10003 FORMAT(('B L O C K'),1X,I3)
10004 FORMAT(('NUMBER OF CELLS IN I - DIRECTION ....  IM  ='),2X,I3) 
10005 FORMAT(('NUMBER OF CELLS IN J - DIRECTION ....  JM  ='),2X,I3)
10006 FORMAT(('NUMBER OF CELLS IN K - DIRECTION ....  KM  ='),2X,I3)
10007 FORMAT(('Connection West Side ................  NCW ='),2X,I3)
10008 FORMAT(('Connection East Side ................  NCE ='),2X,I3)  
10009 FORMAT(('NUMBER OF SECTIONS NORTH SIDE .......  NSN ='),2X,I3)    
10010 FORMAT(('   CONNECTION NORTH SIDE ............  NCN ='),2X,I3)
10011 FORMAT(('   Start at connected Block..........  KNW ='),2X,I3)
10012 FORMAT(('   End at connected Block............  KNE ='),2X,I3)
10013 FORMAT(('NUMBER OF SECTIONS SOUTH SIDE .......  NSS ='),2X,I3)
10014 FORMAT(('   CONNECTION SOUTH SIDE ............  NCS ='),2X,I3)
10015 FORMAT(('   Start at connected Block..........  KSW ='),2X,I3)
10016 FORMAT(('   End at connected Block............  KSE ='),2X,I3)
10017 FORMAT(('Connection Front Side ...............  NCF ='),2X,I3)
10018 FORMAT(('Connection Back Side ................  NCB ='),2X,I3)  
10019 FORMAT(('   thermal BC .......................  IAS ='),2X,I3)
10020 FORMAT(('   constant flux / vector ...........  IBW ='),2X,I3)
10021 FORMAT(('   wall temperature..................  TW  ='),F5.0)
      DEALLOCATE(X,Y,Z)
      CLOSE(10)
      CLOSE(20)


END PROGRAM CUBE
