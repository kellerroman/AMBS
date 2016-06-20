module const_mod
implicit none
public
   integer, parameter :: REAL_KIND = selected_real_kind(15 , 307)
   integer, parameter :: INT_KIND  = selected_int_kind(8)


   integer, parameter :: DIR_EASTWEST        = 1
   integer, parameter :: DIR_SOUTHNORTH      = 2
   integer, parameter :: DIR_FRONTBACK       = 3
   
   integer, parameter :: DIR_WEST            = 1
   integer, parameter :: DIR_EAST            = 2
   integer, parameter :: DIR_SOUTH           = 3
   integer, parameter :: DIR_NORTH           = 4
   integer, parameter :: DIR_FRONT           = 5
   integer, parameter :: DIR_BACK            = 6

   integer, parameter :: VEC_RHO             = 1
   integer, parameter :: VEC_SPU             = 2
   integer, parameter :: VEC_SPV             = 3
   integer, parameter :: VEC_SPW             = 4
   integer, parameter :: VEC_ENE             = 5

   integer, parameter :: GRAD_DX             = 1
   integer, parameter :: GRAD_DY             = 2
   integer, parameter :: GRAD_DZ             = 3

   integer, parameter :: GRAD_SPU            = 1
   integer, parameter :: GRAD_SPV            = 2
   integer, parameter :: GRAD_SPW            = 3
   integer, parameter :: GRAD_TEMP           = 4

   integer, parameter :: TAU_TXX             = 1
   integer, parameter :: TAU_TYY             = 2
   integer, parameter :: TAU_TZZ             = 3
   integer, parameter :: TAU_TXY             = 4
   integer, parameter :: TAU_TXZ             = 5
   integer, parameter :: TAU_TYZ             = 6

   integer, parameter :: BC_OUTFLOW          = 1
   integer, parameter :: BC_INFLOW           = 2
   integer, parameter :: BC_WALL             = 3
   integer, parameter :: BC_SYMMETRY         = 4

   real(REAL_KIND),parameter :: EPSI    = 1.00E-20_REAL_KIND
   real(REAL_KIND),parameter :: ZERO    = 0.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: ONE     = 1.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: TWO     = 2.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: HALF    = 0.50E+00_REAL_KIND
   real(REAL_KIND),parameter :: QUARTER = 0.25E+00_REAL_KIND
   real(REAL_KIND),parameter :: FIFTH   = 0.20E+00_REAL_KIND
   real(REAL_KIND),parameter :: TWOTHIRD= 2.00E+00_REAL_KIND / 3.00E+00_REAL_KIND

   real(REAL_KIND),parameter :: GAMMA   = 1.40E+00_REAL_KIND   ! Ratio of specific heats
   real(REAL_KIND),parameter :: RGAS    = 287.102E0_REAL_KIND
end module const_mod
