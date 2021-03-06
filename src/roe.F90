!********************************************************************************
!* -- 3D Roe's Flux Function with an entropy fix and without tangent vectors --
!*
!* NOTE: This version does not use any tangent vector.
!*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
!*
!* This subroutine computes the Roe flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Roe flux is implemented in the following form:
!*
!*   Numerical flux = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ], 
!*
!*  where
!*
!*    An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.
!*
!* The dissipation term, |An|dU, is actually computed as
!*
!*     sum_{k=1,5} |lambda_k| * (LdU)_k * r_k,
!*
!* where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
!* and r_k is the k-th right-eigenvector.
!*
!* ------------------------------------------------------------------------------
!*  Input: primL(1:5) =  Left state (rhoL, uL, vL, wR, pL)
!*         primR(1:5) = Right state (rhoR, uR, vR, wR, pR)
!*           njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!*
!* Output:  flux(1:5) = The Roe flux with an entropy fix
!* ------------------------------------------------------------------------------
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, April 2012. http://www.cfdbooks.com
!********************************************************************************
 subroutine inviscid_roe_n(primL, primR, njk,  num_flux)

 implicit none
 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 real(p2), intent( in) :: primL(5), primR(5) ! Input: primitive variables
 real(p2), intent( in) :: njk(3)             ! Input: face normal vector

!Output
 real(p2), intent(out) :: num_flux(5)        ! Output: numerical flux

!Some constants
 real(p2) ::  zero = 0.0_p2
 real(p2) ::   one = 1.0_p2
 real(p2) ::   two = 2.0_p2
 real(p2) ::  half = 0.5_p2
 real(p2) :: fifth = 0.2_p2

!Local variables
 real(p2) :: nx, ny, nz                   ! Normal vector
 real(p2) :: uL, uR, vL, vR, wL, wR       ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR           ! Primitive variables.
 real(p2) :: qnL, qnR                     ! Normal velocities
 real(p2) :: aL, aR, HL, HR               ! Speed of sound, Total enthalpy
 real(p2) :: RT,rho,u,v,w,H,a,qn          ! Roe-averages
 real(p2) :: drho,dqn,dp,LdU(4)           ! Wave strengths
 real(p2) :: du, dv, dw                   ! Velocity differences
 real(p2) :: ws(4), R(5,4)                ! Wave speeds and right-eigenvectors
 real(p2) :: dws(4)                       ! Width of a parabolic fit for entropy fix
 real(p2) :: fL(5), fR(5), diss(5)        ! Fluxes ad dissipation term
 real(p2) :: gamma = 1.4_p2               ! Ratio of specific heats

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Primitive and other variables.

!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL                                        !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity

!Wave Strengths

   drho = rhoR - rhoL !Density difference
     dp =   pR - pL   !Pressure difference
    dqn =  qnR - qnL  !Normal velocity difference

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho                         !Shear wave strength (not really, just a factor)

!Absolute values of the wave Speeds

  ws(1) = abs(qn-a) !Left-moving acoustic wave
  ws(2) = abs(qn)   !Entropy wave
  ws(3) = abs(qn+a) !Right-moving acoustic wave
  ws(4) = abs(qn)   !Shear waves

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Right Eigenvectors
!Note: Two shear wave components are combined into one, so that tangent vectors
!      are not required. And that's why there are only 4 vectors here.
!      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx
  R(3,1) = v - a*ny
  R(4,1) = w - a*nz
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx
  R(3,3) = v + a*ny
  R(4,3) = w + a*nz
  R(5,3) = H + a*qn

! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  R(1,4) = zero
  R(2,4) = du - dqn*nx
  R(3,4) = dv - dqn*ny
  R(4,4) = dw - dqn*nz
  R(5,4) = u*du + v*dv + w*dw - qn*dqn

!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
         + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

  num_flux = half * (fL + fR - diss)

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a

 end subroutine inviscid_roe_n
!--------------------------------------------------------------------------------
