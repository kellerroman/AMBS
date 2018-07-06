module const_mod
implicit none
public
   integer, parameter :: REAL_KIND = selected_real_kind(15 , 307)
   integer, parameter :: INT_KIND  = selected_int_kind(8)

   enum, bind(C)
      enumerator :: DIR_EASTWEST = 1, DIR_SOUTHNORTH, DIR_FRONTBACK
   
      enumerator :: DIR_WEST= 1, DIR_EAST, DIR_SOUTH, DIR_NORTH, DIR_FRONT, DIR_BACK

      enumerator :: VEC_RHO = 1, VEC_SPU, VEC_SPV, VEC_SPW, VEC_ENE

      enumerator :: GRAD_DX = 1, GRAD_DY, GRAD_DZ

      enumerator :: GRAD_SPU = 1, GRAD_SPV, GRAD_SPW, GRAD_TEMP

      enumerator :: TAU_TXX = 1, TAU_TYY, TAU_TZZ, TAU_TXY, TAU_TXZ, TAU_TYZ

      enumerator :: BC_OUTFLOW = -6, BC_INFLOW_SUB, BC_INFLOW, BC_WALL, BC_SYMMETRY, BC_PERIODIC, BC_FARFIELD

      enumerator :: EQU_TYPE_EULER = 1, EQU_TYP_NS

   end enum

   real(REAL_KIND),parameter :: EPSI    = 1.00E-20_REAL_KIND
   real(REAL_KIND),parameter :: ZERO    = 0.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: ONE     = 1.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: TWO     = 2.00E+00_REAL_KIND
   real(REAL_KIND),parameter :: HALF    = 0.50E+00_REAL_KIND
   real(REAL_KIND),parameter :: QUARTER = 0.25E+00_REAL_KIND
   real(REAL_KIND),parameter :: FIFTH   = 0.20E+00_REAL_KIND
   real(REAL_KIND),parameter :: TWOTHIRD= 2.00E+00_REAL_KIND / 3.00E+00_REAL_KIND

   real(REAL_KIND),parameter :: GAMMA   = 1.40E+00_REAL_KIND   ! Ratio of specific heats
   real(REAL_KIND),parameter :: GM1     = 0.40E+00_REAL_KIND   ! Ratio of specific heats
   real(REAL_KIND),parameter :: RGAS    = 287.102E0_REAL_KIND

   character(len = 5), parameter  :: DIR_NAMES(6) = ["WEST " &
                                                    ,"EAST " &
                                                    ,"SOUTH" &
                                                    ,"NORTH" &
                                                    ,"FRONT" &
                                                    ,"BACK " ]
   contains
   function bc_names(i) result(string)
      implicit none
      character(len = :), allocatable :: string
      integer, intent(in) :: i
      character(len = 2) :: tmp
      select case(i)
         case(1:)
            write(tmp,'(I0)') i
            string = "BloCo "//trim(tmp)
         case(BC_OUTFLOW)
            string = "OUTFLOW"
         case(BC_INFLOW_SUB)
            string = "INFLOW_SUB"
         case(BC_INFLOW)
            string = "INFLOW"
         case(BC_WALL)
            string = "WALL"
         case(BC_SYMMETRY)
            string = "SYMMETRY"
         case(BC_PERIODIC)
            string = "PERIODIC"
         case(BC_FARFIELD)
            string = "FARFIELD"
         end select
   end function
!   character(len = 8), parameter  :: BC_NAMES(-4:0) =  ["OUTFLOW " &
!                                                    ,"INFLOW  " &
!                                                    ,"WALL    " &
!                                                    ,"SYMMETRY" &
!                                                    ,"PERIODIC" ]
end module const_mod
