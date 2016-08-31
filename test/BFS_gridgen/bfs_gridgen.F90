program BFS_GRIDGEN
!AUTHOR:       ROMAN KELLER
!START DATE:   22.04.2016
!LAST CHANGES: 22.04.2016
!PURPOSE:      ERSTELLT EIN LES GITTER FUER DIE Backward Facing Step 
! CHANGELOG:
! 22.04.2016, RK: START AUF BASIS DES DIT GRIDGEN VON MARKUS KINDLER/ROMAN KELLER

!         __________________________________________
!        |                   |                      |
!        |                   |   /                  |
!  dy1   | nj1    1          |dz/       3           |
!        |       ni1         | / nk                 |
!    dw1 |___________________|/_____________________|
!   dinc1                    |                      |
!               dx1          |                      |
!                      dy2   |  nj2     2           |
!                            |         ni2          |
!                            |______________________|
!                           dw2          
!                          dinc2       dx2
implicit none

   integer               , parameter :: AXSYM        = 2
   integer               , parameter :: NUM_OF_BLOCK = 3                                            ! no. of blocks
   integer               , parameter :: DIM          = 3                                            ! no. of blocks
 
   integer               , parameter :: NI1 = 128
   integer               , parameter :: NI2 = 5 * NI1
   integer               , parameter :: NJ2 = 40 
   integer               , parameter :: NJ1 = 4 * NJ2
   integer               , parameter :: NK1 = 32

   integer               , parameter :: POS_MIN = 1
   integer               , parameter :: POS_MAX = 2
 
   integer               , parameter :: POS_X = 1
   integer               , parameter :: POS_Y = 2
   integer               , parameter :: POS_Z = 3

   real(kind = 8 )       , parameter :: H   =  3.80E-02 
   real(kind = 8 )       , parameter :: DX1 =  4.00E-00 * H
   real(kind = 8 )       , parameter :: DX2 = 20.00E-00 * H
   real(kind = 8 )       , parameter :: DY1 =  4.00E-00 * H 
   real(kind = 8 )       , parameter :: DY2 =  1.00E-00 * H
   real(kind = 8 )       , parameter :: DK  =  2.00E-00 * H
 
   type :: block
      integer :: ncoords(DIM)
      real(kind = 8), allocatable :: coords(:,:,:,:)
      ! dim,i,j,k
      real(kind = 8) :: coords_min_max(2,DIM)
      ! min/max,x/y/k
   end type block
 
 
   type(block) :: blocks(NUM_OF_BLOCK)
   integer :: b,i,j,k,n,lv
   integer :: ni,nj,nk
   real(kind = 8 ), allocatable :: step(:)
   real(kind = 8 ) :: ci



   write(*,*) "Grid Generator for Backward Facing Step"
   b = 1
   blocks(b) % ncoords = [NI1,NJ1,NK1]
   blocks(b) % coords_min_max(POS_MIN,POS_X) = -DX1
   blocks(b) % coords_min_max(POS_MAX,POS_X) = 0.0d0
   blocks(b) % coords_min_max(POS_MIN,POS_Y) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_Y) = DY1 
   blocks(b) % coords_min_max(POS_MIN,POS_Z) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_Z) = DK

   b = 2
   blocks(b) % ncoords = [NI2,NJ2,NK1]
   blocks(b) % coords_min_max(POS_MIN,POS_X) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_X) = DX2  
   blocks(b) % coords_min_max(POS_MIN,POS_Y) = -DY2
   blocks(b) % coords_min_max(POS_MAX,POS_Y) = 0.0d0
   blocks(b) % coords_min_max(POS_MIN,POS_Z) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_Z) = DK

   b = 3
   blocks(b) % ncoords = [NI2,NJ1,NK1]
   blocks(b) % coords_min_max(POS_MIN,POS_X) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_X) = DX2  
   blocks(b) % coords_min_max(POS_MIN,POS_Y) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_Y) = DY1
   blocks(b) % coords_min_max(POS_MIN,POS_Z) = 0.0d0
   blocks(b) % coords_min_max(POS_MAX,POS_Z) = DK

   do b = 1, NUM_OF_BLOCK
      ni = blocks(b) % ncoords(1)+1
      nj = blocks(b) % ncoords(2)+1
      nk = blocks(b) % ncoords(3)+1
      allocate ( blocks(b) % coords( DIM, ni, nj, nk)  )
      do n = 1,DIM
         allocate ( step(blocks(b) % ncoords(n) + 1))
         call constant_cell_size(step                                   &
                                ,blocks(b) % coords_min_max(POS_MIN,n)  &
                                ,blocks(b) % coords_min_max(POS_MAX,n)  &
                                ,blocks(b) % ncoords(n))
         if ((n == 1 .and. b > 1 ) &
         .or.(n==2)) then
            call variableCellSizeStartDefined(step&
                                ,blocks(b) % coords_min_max(POS_MIN,n)  &
                                ,blocks(b) % coords_min_max(POS_MAX,n)  &
                                ,1.0d-6,1.3D0&
                                ,blocks(b) % ncoords(n),ci)
         end if
         do k = 1,nk
            do j = 1,nj
               do i = 1, ni
                  if (n == 1) then
                     lv = i
                  else if (n == 2) then
                     lv = j
                  else 
                     lv = k
                  end if
                  blocks(b) % coords(n,i,j,k) = step(lv)
               end do
            end do
         end do
         deallocate (step)
      end do
   end do
 
 
 
   open(10,file='git.bin',form="unformatted",access="stream",status="replace")
   write(10) AXSYM
   write(10) NUM_OF_BLOCK 
   do b = 1, NUM_OF_BLOCK
      ni = blocks(b) % ncoords(1)+1
      nj = blocks(b) % ncoords(2)+1
      nk = blocks(b) % ncoords(3)+1
      write(10) ni,nj,nk 
   end do
   do b = 1, NUM_OF_BLOCK
      ni = blocks(b) % ncoords(1)+1
      nj = blocks(b) % ncoords(2)+1
      nk = blocks(b) % ncoords(3)+1
      do k = 1,nk
         do j = 1,nj
            do i = 1,ni
               write(10) blocks(b) % coords(:,i,j,k) 
            end do
         end do
      end do
   end do 
   close(10)
   write(*,*) "   DONE  "
end program BFS_GRIDGEN
subroutine interpolate(x,y,n)
! this routine interpolates along the x-Axis the y-Values
implicit none
real( kind = 8), intent(in) :: x(n+1)
real( kind = 8), intent(inout) :: y(n+1)
integer , intent(in) :: n
! number of Cells, dimension of x and y are n+1 because those are grid-points

real( kind = 8) :: l

integer :: i

l = (y(n+1) - y(1)) / (dble(x(n+1)) - dble(x(1)) )

do i = 2, n
   y(i) = y(1) + l * (dble(x(i)) - dble(x(1)))
enddo

end subroutine
subroutine constant_cell_size(x,x_start,x_end,ni)
implicit none
real( kind = 8 ), intent(inout) :: x(ni+1)
real( kind = 8 ), intent(in) :: x_start, x_end
integer :: ni
integer :: i,nstep
real( kind = 8) :: stepsize
nstep = ni
stepsize = (x_end-x_start) / dble(nstep)
x(1) = x_start
do i = 1, nstep
   x(i+1) = x(i) + stepsize
end do
end subroutine constant_cell_size



subroutine variableCellSizeEndDefined(x,x_start,x_end,dn,ci_max,n,ci)
   implicit none
   integer, intent(in) :: n
   real( kind = 8), intent(in) :: x_start,x_end,dn,ci_max
   real( kind = 8), intent(out) :: ci,x(n+1)

   integer :: num_iter , i, n1,n2
   real( kind = 8) :: f,l,dn2,l1,l2
   write(*,*) "================================== START variableCellSizeEndDefined ===================================="
   l = x_end - x_start
   num_iter = 0
   ci = ci_max
   n2 = n
   l2 = l
   if (dn*(1.0d0-ci**n)/(1.0d0-ci) < l) then
      write(*,*) "ERROR: variableCellSizeEndDefine"
      write(*,*) "laenge:",l, "ist mit dn:",dn,"ci:",ci,"n:",n,"nicht darstellbar"
      write(*,*) "maximale Laenge:",dn*(1.0d0-ci**n)/(1.0d0-ci)
      stop 1    
   end if
!   write(*,*) "==========",l,n,dn,ci_max

   do
      num_iter = num_iter + 1
      dn2 = l2 / n2
      n1 = 1 + log(dn2/dn)/log(ci)
      n2 = n - n1
      l1 = dn * (1.0d0-ci**(n1))/(1.0d0-ci)
      l2 = l- l1
      f = n2 * dn2 + l1
      write(*,*) num_iter,n1,n2,l1,l2,dn2,f
      
      if (abs(f-l) <= 1e-8) then
         exit
      end if

      if (num_iter >= 100) then
         write(*,*) "Not Converged: variableCellSizeEndDefine"
         stop 1
      end if
   end do
   x(1) = x_start
   dn2 = l2 / n2

   do i = 1,n2
      x(i+1) = x(i) + dn2
   end do
   do i = n2+1,n
      x(i+1) = x(i) + dn * ci ** (n-i)
   end do
   write(*,*) dn2,dn * ci ** (n-n2-1)
   x(n+1) = x_end
   write(*,*) "======================= DONE ===== variableCellSizeEndDefined ===================================="
end subroutine variableCellSizeEndDefined

subroutine variableCellSizeStartDefined(x,x_start,x_end,dn,ci_max,n,ci)
   implicit none
   integer, intent(in) :: n
   real( kind = 8), intent(in) :: x_start,x_end,dn,ci_max
   real( kind = 8), intent(out) :: ci,x(n+1)

   integer :: num_iter , i, n1,n2
   real( kind = 8) :: f,l,dn2,l1,l2

   l = x_end - x_start
   num_iter = 0
   ci = ci_max
   n2 = n
   l2 = l

   if (dn*(1.0d0-ci**n)/(1.0d0-ci) < l) then
      write(*,*) "ERROR: variableCellSizeStartDefine"
      write(*,*) "laenge:",l, "ist mit dn:",dn,"ci:",ci,"n:",n,"nicht darstellbar"
      write(*,*) "maximale lã„nge:",dn*(1.0d0-ci**n)/(1.0d0-ci)
      stop 1
   end if

   !write(*,*) "==========",l,n,dn,ci_max
   do
      num_iter = num_iter + 1
      dn2 = l2 / n2
      n1 = 1 + log(dn2/dn)/log(ci)
      n2 = n - n1
      l1 = dn * (1.0d0-ci**n1)/(1.0d0-ci)
      l2 = l- l1
      f = n2*dn2+l1
   !   write(*,*) n1,n2,l1,l2,dn2,f
      
      if (abs(f-l) <= 1e-8) then
         exit
      end if

      if (num_iter >= 100) then
         write(*,*) "Not Converged: variableCellSizeStartDefine"
         stop 1
      end if
   end do
   x(1) = x_start

   do i = 1,n1
      x(i+1) = x(i) + dn * ci ** (i-1)
   end do
   do i = n1+1,n
      x(i+1) = x(i) + dn2
   end do
   x(n+1) = x_end
end subroutine variableCellSizeStartDefined
