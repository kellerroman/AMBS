module const_mod
implicit none
public
   integer, parameter :: REAL_KIND = selected_real_kind(15 , 307)
   integer, parameter :: INT_KIND  = selected_int_kind(8)


   integer, parameter :: DIR_EASTWEST        = 1
   integer, parameter :: DIR_SOUTHNORTH      = 2
   integer, parameter :: DIR_FRONTBACK       = 3
   
   integer, parameter :: DIR_EAST            = 1
   integer, parameter :: DIR_WEST            = 2
   integer, parameter :: DIR_SOUTH           = 3
   integer, parameter :: DIR_NORTH           = 4
   integer, parameter :: DIR_FRONT           = 5
   integer, parameter :: DIR_BACK            = 6

   real(REAL_KIND),parameter :: EPSI    = 1.00E-20_REAL_KIND
   real(REAL_KIND),parameter :: ZERO    = 0.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: ONE     = 1.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: TWO     = 2.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: HALF    = 0.50E+00_REAL_KIND
   real(REAL_KIND),parameter :: QUARTER = 0.25E+00_REAL_KIND
   real(REAL_KIND),parameter :: FIFTH   = 0.20E+00_REAL_KIND
   real(REAL_KIND),parameter :: GAMMA   = 1.40E+00_REAL_KIND               ! Ratio of specific heats
end module const_mod
