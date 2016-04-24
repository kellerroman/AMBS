module const_mod
implicit none
public
   integer, parameter :: REAL_KIND = selected_real_kind(15 , 307)
   integer, parameter :: INT_KIND  = selected_int_kind(8)


   integer, parameter :: DIR_EASTWEST        = 1
   integer, parameter :: DIR_SOUTHNORTH      = 1
   integer, parameter :: DIR_FRONTBACK       = 1
   
   integer, parameter :: DIR_EAST            = 1
   integer, parameter :: DIR_WEST            = 2
   integer, parameter :: DIR_SOUTH           = 3
   integer, parameter :: DIR_NORTH           = 4
   integer, parameter :: DIR_FRONT           = 5
   integer, parameter :: DIR_BACK            = 6

   real(REAL_KIND),parameter :: EPSI = 1.0E-10_REAL_KIND
   real(REAL_KIND),parameter :: HALF = 5.0E-1_REAL_KIND
end module const_mod
