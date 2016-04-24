module test_mod
implicit none
   integer, parameter :: REAL_KIND = selected_real_kind(15 , 307)
   integer, parameter :: INT_KIND  = selected_int_kind(8)
contains
   subroutine test_sub
   implicit none
      real(REAL_KIND) :: a(3) = [1,2,3]
      real(REAL_KIND) :: b(3) = [3,21,1]
!     real(REAL_KIND) :: c(3) = [1,2,3]

      write(*,*) test_func(a,b)
   end subroutine

   function test_func(a,b) result ( vec )
   implicit none


      real(REAL_KIND), dimension(3), intent(in) :: a,b
      real(REAL_KIND), dimension(3)             :: vec 

      vec = a + b

   end function


end module



program func_test
use test_mod
implicit none


call test_sub


      



end program

