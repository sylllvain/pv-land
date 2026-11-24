!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE compute_U_plus              ######
!######                                                      ######
!##################################################################
!##################################################################

subroutine compute_U_plus(z, theta, a, beta_deg, L, clearance, theta_0, Lo, Uplus)

!============================================================
! Purpose:
! --------
! Compute the dimensionless velocity U+ at height z, accounting
! for whether the point lies above or within the PV canopy, and
! for the atmospheric stability above the canopy.
!
! Inputs:
! -------
! z         : vertical coordinate above ground [m]
! theta     : wind direction angle [rad]
! a         : row spacing between PV panels [m]
! beta_deg  : PV panel tilt angle [degrees]
! L         : PV panel length (panel width) [m]
! clearance : clearance height below the panel [m]
! theta_0   : park orientation angle [rad]
! Lo        : Obukhov length [m]. 
!
! Output:
! -------
! Uplus     : dimensionless velocity U+ [-]
!
! Note:
! -----
! Stability is introduced in the above-canopy profile following
! Harman & Finnigan (2007), Eqs. 14-15, through two
! contributions:
!   (i)  the inertial-sublayer MOST momentum function psi_m, and
!   (ii) the stability dependence of the roughness-sublayer function
!        psi_hat_m, i.e. phi_m(z'/Lo) inside the RSL integral (Eq. 15).
!============================================================

  implicit none

  !-------------------------
  ! Input variables
  !-------------------------
  real, intent(in)  :: z          ! Height above ground [m]
  real, intent(in)  :: theta      ! Wind direction [rad]
  real, intent(in)  :: a          ! Row spacing [m]
  real, intent(in)  :: beta_deg   ! Panel tilt angle [deg]
  real, intent(in)  :: L          ! Panel length (width) [m]
  real, intent(in)  :: clearance  ! Ground clearance [m]
  real, intent(in)  :: theta_0    ! Park orientation angle [rad]
  real, intent(in)  :: Lo         ! Obukhov length [m] (<0 unstable, >0 stable)

  !-------------------------
  ! Output variable
  !-------------------------
  real, intent(out) :: Uplus      ! Dimensionless velocity [-]

  !-------------------------
  ! Local constants
  !-------------------------
  real, parameter :: pi = 3.141592653589793  ! Pi [-]

  !-------------------------
  ! Local variables
  !-------------------------
  real :: beta       ! Panel tilt angle [rad]
  real :: h_PV       ! Effective canopy height [m]
  real :: lambda_f   ! Frontal density parameter [-]

  ! Convert tilt angle from degrees to radians
  beta = beta_deg * pi / 180.0

  ! Compute effective PV height
  h_PV = L * sin(beta) + clearance

  ! Compute frontal density parameter
  call compute_lambda_f(a, L, beta, lambda_f)

  ! Select velocity profile depending on height
  if (z >= h_PV) then
     call U_above(z, theta, lambda_f, h_PV, beta_deg, theta_0, Lo, Uplus)
  else
     call U_within(z, theta, lambda_f, h_PV, beta_deg, theta_0, Lo, Uplus)
  endif

end subroutine compute_U_plus


!============================================================
! SUBROUTINE: compute_lambda_f
!
! Purpose:
! --------
! Compute the frontal density parameter lambda_f.
!
! Inputs:
! -------
! a     : row spacing [m]
! L     : panel length [m]
! beta  : panel tilt angle [rad]
!
! Output:
! -------
! lambda_f : frontal density parameter [-]
!============================================================
subroutine compute_lambda_f(a, L, beta, lambda_f)
  implicit none

  real, intent(in)  :: a        ! Row spacing [m]
  real, intent(in)  :: L        ! Panel length [m]
  real, intent(in)  :: beta     ! Tilt angle [rad]
  real, intent(out) :: lambda_f ! Frontal density [-]

  lambda_f = (L * sin(beta)) / (a + L * cos(beta))
end subroutine compute_lambda_f


!============================================================
! SUBROUTINE: alpha1_lambda
!
! Purpose:
! --------
! Compute coefficient alpha1 as a function of lambda_f.
!
! Input:
! ------
! lambda_f : frontal density [-]
!
! Output:
! -------
! alpha1   : displacement scaling coefficient [-]
!============================================================
subroutine alpha1_lambda(lambda_f, alpha1)
  implicit none

  real, intent(in)  :: lambda_f ! Frontal density [-]
  real, intent(out) :: alpha1   ! Coefficient [-]

  real, parameter :: alpha_inf = 0.95 / (3.38 + 1.0) ! Asymptotic value [-]
  real, parameter :: lambda_c  = 0.48                ! Transition density [-]

  alpha1 = (alpha_inf - 0.15) * lambda_f**2 / &
           (lambda_f**2 + lambda_c**2) + 0.15
end subroutine alpha1_lambda


!============================================================
! SUBROUTINE: A1_gamma
!
! Purpose:
! --------
! Compute the asymmetric gamma-like coefficient A1(lambda_f).
!
! Input:
! ------
! lambda_f : frontal density [-]
!
! Output:
! -------
! A1       : directional amplitude coefficient [-]
!============================================================
subroutine A1_gamma(lambda_f, A1)
  implicit none

  real, intent(in)  :: lambda_f ! Frontal density [-]
  real, intent(out) :: A1       ! Amplitude coefficient [-]

  real, parameter :: lambda_peak = 0.2   ! Peak density [-]
  real, parameter :: p = 0.3             ! Shape parameter [-]
  real, parameter :: xpt = 0.204         ! Reference x-point [-]
  real, parameter :: ypt = 0.138         ! Reference y-point [-]

  real :: q       ! Exponential decay factor [-]
  real :: h_xpt   ! Normalization helper [-]
  real :: A1_max  ! Maximum amplitude [-]

  if (lambda_f <= 0.0) then
     A1 = 0.0
     return
  endif

  q = p / lambda_peak
  h_xpt = (xpt / lambda_peak)**p * exp(-q * (xpt - lambda_peak))
  A1_max = ypt / h_xpt

  A1 = A1_max * (lambda_f / lambda_peak)**p * &
       exp(-q * (lambda_f - lambda_peak))
end subroutine A1_gamma


!============================================================
! SUBROUTINE: d_theta
!
! Purpose:
! --------
! Compute the displacement height d as a function of wind
! direction and canopy geometry.
!
! Inputs:
! -------
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
! h_PV      : canopy height [m]
!
! Output:
! -------
! d         : displacement height [m]
!============================================================
subroutine d_theta(theta, lambda_f, h_PV, theta_0, d)
  implicit none

  real, intent(in)  :: theta     ! Wind direction [rad]
  real, intent(in)  :: lambda_f  ! Frontal density [-]
  real, intent(in)  :: h_PV      ! Canopy height [m]
  real, intent(in)  :: theta_0   ! Park orientation angle [rad]
  real, intent(out) :: d         ! Displacement height [m]

  real, parameter :: A = 3.38 ! Empirical constant [-]

  real :: alpha1 ! Scaling coefficient [-]

  call alpha1_lambda(lambda_f, alpha1)

  d = alpha1 * (A * cos(theta - theta_0)**2 + 1.0) * h_PV
end subroutine d_theta


!============================================================
! SUBROUTINE: r_beta
!
! Purpose:
! --------
! Compute the asymmetry ratio r as a function of panel tilt.
!
! Input:
! ------
! beta_deg : panel tilt angle [deg]
!
! Output:
! -------
! r        : asymmetry ratio [-]
!============================================================
subroutine r_beta(beta_deg, r)
  implicit none

  real, intent(in)  :: beta_deg ! Tilt angle [deg]
  real, intent(out) :: r        ! Ratio [-]

  r = 1.70195e-7 * beta_deg**4 - 3.06351e-5 * beta_deg**3 + &
      2.19236e-3 * beta_deg**2 - 7.32403e-2 * beta_deg + 1.0
end subroutine r_beta


!============================================================
! SUBROUTINE: beta_theta
!
! Purpose:
! --------
! Compute the directional beta coefficient used in roughness
! length estimation.
!
! Inputs:
! -------
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
! beta_deg  : panel tilt angle [deg]
!
! Output:
! -------
! beta_t    : directional beta coefficient [-]
!============================================================
subroutine beta_theta(theta, lambda_f, beta_deg, theta_0, beta_t)
  implicit none

  real, intent(in)  :: theta, lambda_f, beta_deg, theta_0
  real, intent(out) :: beta_t

  real, parameter :: C = 8.72e-2
  real, parameter :: pi = 3.141592653589793

  real :: A1, r
  real :: term1, term2

  call A1_gamma(lambda_f, A1)
  call r_beta(beta_deg, r)

  term1 = max(cos(theta - theta_0), 0.0)**2
  term2 = max(cos(theta - theta_0 - pi), 0.0)**2

  beta_t = A1 * term1 + A1 * r * term2 + C
end subroutine beta_theta


!============================================================
! SUBROUTINE: z0_theta
!
! Purpose:
! --------
! Compute the aerodynamic roughness length z0.
!
! Inputs:
! -------
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
! h_PV      : canopy height [m]
! beta_deg  : panel tilt angle [deg]
!
! Output:
! -------
! z0        : roughness length [m]
!============================================================
subroutine z0_theta(theta, lambda_f, h_PV, beta_deg, theta_0, z0)
  implicit none

  real, intent(in)  :: theta, lambda_f, h_PV, beta_deg, theta_0
  real, intent(out) :: z0

  real, parameter :: kappa = 0.387 ! von Kármán constant [-]

  real :: z_m     ! Mixing height [m]
  real :: d       ! Displacement height [m]
  real :: beta_t  ! Directional beta [-]

  z_m = 3.26 * h_PV

  call d_theta(theta, lambda_f, h_PV, theta_0, d)
  call beta_theta(theta, lambda_f, beta_deg, theta_0, beta_t)

  z0 = (z_m - d) * exp(-kappa / beta_t)
end subroutine z0_theta


!============================================================
! SUBROUTINE: c1_theta
!
! Purpose:
! --------
! Compute coefficient c1 controlling wake effects above canopy.
!
! Input:
! ------
! theta : wind direction [rad]
!
! Output:
! -------
! c1    : wake coefficient [-]
!============================================================
subroutine c1_theta(theta, theta_0, c1)
  implicit none

  real, intent(in)  :: theta
  real, intent(in)  :: theta_0
  real, intent(out) :: c1

  real :: val

  val = 0.28 * exp(1.51 * cos(theta - theta_0)) + &
        0.18 * exp(1.51 * cos(theta - theta_0 - 3.141592653589793)) - 0.61

  c1 = max(val, 0.0)
end subroutine c1_theta


!============================================================
! SUBROUTINE: c3_theta
!
! Purpose:
! --------
! Compute coefficient c3 controlling exponential decay inside
! the canopy.
!
! Inputs:
! -------
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
!
! Output:
! -------
! c3        : decay coefficient [-]
!============================================================
subroutine c3_theta(theta, lambda_f, theta_0, c3)
  implicit none

  real, intent(in)  :: theta, lambda_f, theta_0
  real, intent(out) :: c3

  c3 = (5.0 * lambda_f) / (0.48 * cos(theta - theta_0) + 1.0)
end subroutine c3_theta


!============================================================
! SUBROUTINE: alpha_w
!
! Purpose:
! --------
! Compute weighting factor between logarithmic and exponential
! profiles inside the canopy.
!
! Input:
! ------
! lambda_f : frontal density [-]
!
! Output:
! -------
! w        : weighting factor [0Ã¢â‚¬â€œ1]
!============================================================
subroutine alpha_w(lambda_f, w)
  implicit none

  real, intent(in)  :: lambda_f
  real, intent(out) :: w

  w = lambda_f**3 / (lambda_f**3 + 0.04**3)
end subroutine alpha_w


!============================================================
! SUBROUTINE: expint_E1
!
! Purpose:
! --------
! Compute the exponential integral E1(x).
!
! Input:
! ------
! x   : argument (x > 0) [-]
!
! Output:
! -------
! res : exponential integral E1(x) [-]
!============================================================
subroutine expint_E1(x, res)
  implicit none

  real, intent(in)  :: x
  real, intent(out) :: res

  real, parameter :: gamma = 0.5772156649015328606
  real :: term, sum
  integer :: k

  if (x <= 0.0) then
     res = huge(res)
     return
  endif

  if (x <= 1.0) then
     ! Convergent power series (accurate for small x); UNCHANGED so the
     ! neutral profile is bit-for-bit preserved (its arguments are < 1).
     sum = 0.0
     term = 1.0

     do k = 1, 100
        term = -term * x / real(k)
        sum  = sum + term / real(k)
        if (abs(term) < 1.0e-7) exit
     enddo

     res = -log(x) - gamma - sum
  else
     ! Continued fraction (modified Lentz), convergent for all x >= 1.
     block
       !real, parameter :: eps = 1.0e-12, fpmin = 1.0e-300
       real, parameter :: eps = 1.0e-12, fpmin = 1.0e-20 !modif sylvain
       real :: b, c, dd, h, an, del
       integer :: i
       b  = x + 1.0
       c  = 1.0 / fpmin
       dd = 1.0 / b
       h  = dd
       do i = 1, 200
          an = -real(i) * real(i)
          b  = b + 2.0
          dd = 1.0 / (an * dd + b)
          c  = b + an / c
          del = c * dd
          h  = h * del
          if (abs(del - 1.0) < eps) exit
       enddo
       res = h * exp(-x)
     end block
  endif
end subroutine expint_E1


!============================================================
! SUBROUTINE: psi_m_MOST
!
! Purpose:
! --------
! Integrated Monin-Obukhov similarity (MOST) stability function
! for momentum, psi_m(zeta), obtained by integrating the flux-
! gradient relation phi_m used by Harman & Finnigan (2007, Eq. 10):
!
!     phi_m(zeta) = (1 - 16 zeta)^(-1/4)   for zeta <  0  (unstable)
!     phi_m(zeta) = 1 + 5 zeta             for zeta >= 0  (stable)
!
! The integrated forms are the classical Businger-Dyer / Paulson
! (1970) expressions:
!
!   unstable:  psi_m = 2 ln[(1+x)/2] + ln[(1+x^2)/2]
!                        - 2 atan(x) + pi/2,    x = (1-16 zeta)^(1/4)
!   stable:    psi_m = -5 zeta
!
! psi_m(0) = 0 and the two branches join continuously at zeta = 0.
!
! Input:
! ------
! zeta : stability parameter (z-z_d)/Lo [-]
!
! Output:
! -------
! psi  : integrated MOST momentum function psi_m(zeta) [-]
!
! Note:
! -----
! In strongly stable conditions the linear form -5 zeta grows without
! bound (and would drive U+ to unphysical values), because the linear
! Businger-Dyer relation phi_m = 1 + 5 zeta is only valid for weak-to-
! moderate stability. Field data show a departure from linearity and a
! transition to "z-less" stratification for zeta greater than about 1
! (Grachev et al. 2007; Grachev et al. 2013). We therefore cap the
! stable-side argument at zeta_max (parameter below), i.e. we restrict
! the linear function to its range of validity. Unstable side is left
! unchanged (phi_m stays bounded and no runaway occurs there).
!
! NB: zeta_max must be kept identical to the one in psihat_stab_corr.
!============================================================
subroutine psi_m_MOST(zeta, psi)
  implicit none

  real, intent(in)  :: zeta   ! Stability parameter [-]
  real, intent(out) :: psi    ! Integrated MOST function [-]

  real, parameter :: pi = 3.141592653589793
  ! Upper bound on the stable stability parameter (validity limit of the
  ! linear stable flux-profile relation; z-less regime beyond).
  real, parameter :: zeta_max = 1.0

  real :: x      ! (1 - 16 zeta)^(1/4)
  real :: zeta_s ! capped stability parameter (stable side)

  if (zeta < 0.0) then
     ! Unstable stratification
     x = (1.0 - 16.0 * zeta)**0.25
     psi = 2.0 * log(0.5 * (1.0 + x))      &
         +       log(0.5 * (1.0 + x * x))  &
         - 2.0 * atan(x)                   &
         + 0.5 * pi
  else
     ! Stable (and neutral, zeta = 0) stratification.
     ! Cap the argument at zeta_max: for zeta <= zeta_max the validated
     ! linear form is used; beyond, psi_m saturates at -5 zeta_max.
     zeta_s = min(zeta, zeta_max)
     psi = -5.0 * zeta_s
  endif
end subroutine psi_m_MOST


!============================================================
! SUBROUTINE: psihat_stab_corr
!
! Purpose:
! --------
! Stability correction to the roughness-sublayer (RSL) function,
! following the FULL Harman & Finnigan (2007) formulation, Eq. 15.
! With their RSL shape (Eq. 19), 1 - phi_hat_m(z') = c1' exp(-a z'),
! the neutral limit (phi_m = 1) gives the analytic exponential
! integral c1' E1(a (z-z_d)) already used in U_above. Here we return
! ONLY the stability CORRECTION (per unit of the prefactor c1'):
!
!   dpsi = INT_{z-z_d}^{inf} [phi_m(z'/Lo) - 1] exp(-a z') / z' dz'
!
! so that psi_hat_m = c1' * ( E1(a(z-z_d)) + dpsi ). dpsi -> 0 in
! neutral conditions, hence the neutral profile is preserved exactly.
!
! The MOST similarity function phi_m (HF2007 Eq. 10):
!   phi_m(zeta) = (1 - 16 zeta)^(-1/4)   if zeta < 0 (unstable)
!               = 1 + 5 zeta             if zeta >= 0 (stable)
!
! Evaluation:
!   * Stable (Lo > 0): the stable phi_m is used only within its validity
!     range by capping its argument at zeta_max (see psi_m_MOST). The
!     correction stays analytic. Beyond z*, phi_m is held at 1 + 5 zeta_max.
!   * Unstable (Lo < 0): no elementary primitive; evaluated by
!     Gauss-Laguerre quadrature after the substitution
!     z' = (z-z_d) + s/a, which maps the integral to INT_0^inf e^{-s} g(s) ds:
!        dpsi = exp(-a(z-z_d)) * SUM_i w_i [phi_m((..)/Lo)-1]/(a(z-z_d)+s_i).
!     32 nodes give < 5e-4 absolute error down to Lo = -5 m. No cap is
!     needed on the unstable side (phi_m stays bounded).
!
! Inputs:
! -------
! zmd    : z - z_d, height above the virtual origin [m] (> 0)
! arate  : RSL decay rate a = c2 / (4 (h_PV - z_d)^2), as in U_above
! Lo     : Obukhov length [m]
!
! Output:
! -------
! dpsi   : stability correction to psi_hat_m, per unit prefactor [-]
!============================================================
subroutine psihat_stab_corr(zmd, arate, Lo, dpsi)
  implicit none

  real, intent(in)  :: zmd     ! z - z_d [m]
  real, intent(in)  :: arate   ! RSL decay rate a [1/m in the code's convention]
  real, intent(in)  :: Lo      ! Obukhov length [m]
  real, intent(out) :: dpsi

  ! Above this magnitude the Obukhov length is treated as neutral.
  real, parameter :: Lo_neutral = 1.0e5
  ! Upper bound on the stable stability parameter (see psi_m_MOST).
  ! MUST be kept identical to the zeta_max in psi_m_MOST.
  real, parameter :: zeta_max = 1.0

  integer, parameter :: NGL = 32
  real :: xgl(NGL), wgl(NGL)

  ! Gauss-Laguerre nodes (abscissae) for INT_0^inf e^{-s} f(s) ds
  data xgl / &
    4.44893658332674d-02, 2.34526109519618d-01, &
    5.76884629301886d-01, 1.07244875381782d+00, &
    1.72240877644465d+00, 2.52833670642579d+00, &
    3.49221327302199d+00, 4.61645676974977d+00, &
    5.90395850417424d+00, 7.35812673318624d+00, &
    8.98294092421260d+00, 1.07830186325400d+01, &
    1.27636979867427d+01, 1.49311397555226d+01, &
    1.72924543367153d+01, 1.98558609403361d+01, &
    2.26308890131968d+01, 2.56286360224592d+01, &
    2.88621018163235d+01, 3.23466291539647d+01, &
    3.61004948057520d+01, 4.01457197715394d+01, &
    4.45092079957549d+01, 4.92243949873086d+01, &
    5.43337213333969d+01, 5.98925091621340d+01, &
    6.59753772879350d+01, 7.26876280906627d+01, &
    8.01874469779135d+01, 8.87353404178924d+01, &
    9.88295428682840d+01, 1.11751398097938d+02 /

  ! Gauss-Laguerre weights
  data wgl / &
    1.09218341952425d-01, 2.10443107938792d-01, &
    2.35213229669843d-01, 1.95903335972874d-01, &
    1.29983786286068d-01, 7.05786238657156d-02, &
    3.17609125091740d-02, 1.19182148348382d-02, &
    3.73881629461142d-03, 9.80803306614925d-04, &
    2.14864918801358d-04, 3.92034196798780d-05, &
    5.93454161286845d-06, 7.41640457866734d-07, &
    7.60456787912055d-08, 6.35060222662562d-09, &
    4.28138297104074d-10, 2.30589949189124d-11, &
    9.79937928872677d-13, 3.23780165772916d-14, &
    8.17182344342055d-16, 1.54213383339381d-17, &
    2.11979229016355d-19, 2.05442967378800d-21, &
    1.34698258663736d-23, 5.66129413039746d-26, &
    1.41856054546297d-28, 1.91337549445399d-31, &
    1.19224876009820d-34, 2.67151121923966d-38, &
    1.33861694210642d-42, 4.51053619389878d-48 /

  integer :: i
  real :: zeta, phim, s, acc
  real :: zstar, E1_zmd, E1_zstar   ! for the capped stable branch

  if (Lo == 0.0 .or. abs(Lo) >= Lo_neutral) then
     dpsi = 0.0
     return
  endif

  if (Lo > 0.0) then
     ! Stable. The RSL correction is analytic. To respect the validity
     ! range of the linear stable phi_m, the argument z'/Lo is capped at
     ! zeta_max, i.e. phi_m = 1 + 5 min(z'/Lo, zeta_max). Let z* = zeta_max*Lo
     ! be the height where z'/Lo reaches zeta_max. Splitting the integral:
     !   for z1 < z* :  dpsi = (5/(a Lo))[e^{-a z1} - e^{-a z*}]
     !                          + 5 zeta_max E1(a z*)
     !   for z1 >= z*:  dpsi = 5 zeta_max E1(a z1)
     ! (the two pieces join continuously at z1 = z*). z1 = zmd.
     zstar = zeta_max * Lo
     if (zmd < zstar) then
        call expint_E1(arate * zstar, E1_zstar)
        dpsi = (5.0 / (arate * Lo)) * ( exp(-arate * zmd) - exp(-arate * zstar) ) &
             + 5.0 * zeta_max * E1_zstar
     else
        call expint_E1(arate * zmd, E1_zmd)
        dpsi = 5.0 * zeta_max * E1_zmd
     endif
  else
     ! Unstable: Gauss-Laguerre quadrature of the correction integral.
     ! (No cap: the unstable phi_m stays bounded and does not run away.)
     acc = 0.0
     do i = 1, NGL
        s    = xgl(i)
        zeta = (zmd + s / arate) / Lo             ! < 0 (Lo < 0)
        phim = (1.0 - 16.0 * zeta)**(-0.25)
        acc  = acc + wgl(i) * (phim - 1.0) / (arate * zmd + s)
     end do
     dpsi = exp(-arate * zmd) * acc
  endif
end subroutine psihat_stab_corr


!============================================================
! SUBROUTINE: U_above
!
! Purpose:
! --------
! Compute the dimensionless velocity profile above the canopy.
!
! Inputs:
! -------
! z         : height above ground [m]
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
! h_PV      : canopy height [m]
! beta_deg  : panel tilt angle [deg]
! Lo        : Obukhov length [m] (<0 unstable, >0 stable;
!             0 or |Lo| very large -> neutral)
!
! Output:
! -------
! U         : dimensionless velocity [-]
!
! Formulation (Harman & Finnigan 2007, Eqs. 14-15), normalized by u_*:
!
!   U+(z) = (1/kappa) * [ ln((z-z_d)/z0)
!                         - psi_m((z-z_d)/Lo) + psi_m(z0/Lo)   <- inertial sublayer
!                         + psi_hat_m(z, Lo) ]                  <- roughness sublayer
!
! Two distinct stability contributions are included:
!
!  (i) the inertial-sublayer MOST correction, the two psi_m terms;
!
! (ii) the roughness-sublayer function psi_hat_m, which is itself
!      stability dependent (HF2007 Eq. 15): phi_m(z'/Lo) sits INSIDE
!      the RSL integral. We write it as
!         psi_hat_m = c1 exp(-c2/2) * [ E1(a(z-z_d)) + dpsi(z,Lo) ]
!      where E1 is the neutral part and dpsi (from psihat_stab_corr)
!      is the stability correction. dpsi -> 0 in neutral conditions.
!
! Both contributions vanish in neutral stratification
!============================================================
subroutine U_above(z, theta, lambda_f, h_PV, beta_deg, theta_0, Lo, U)
  implicit none

  real, intent(in)  :: z, theta, lambda_f, h_PV, beta_deg, theta_0
  real, intent(in)  :: Lo
  real, intent(out) :: U

  real, parameter :: kappa = 0.387
  real, parameter :: c2 = 0.17
  ! Above this magnitude the Obukhov length is treated as neutral.
  real, parameter :: Lo_neutral = 1.0e5

  real :: d, z0, c1
  real :: arate, arg, E1, dpsihat
  real :: psi_z, psi_0   ! MOST stability terms at (z-d) and z0

  call d_theta(theta, lambda_f, h_PV, theta_0, d)
  call z0_theta(theta, lambda_f, h_PV, beta_deg, theta_0, z0)
  call c1_theta(theta, theta_0, c1)

  ! Roughness-sublayer function.
  ! Neutral part: c1 exp(-c2/2) E1(a (z-z_d)), with a the RSL decay rate.
  arate = c2 / (4.0 * (h_PV - d)**2)
  arg   = arate * (z - d)
  call expint_E1(arg, E1)

  ! Inertial-sublayer MOST correction (HF2007 Eq. 14). Lo = 0 is used as a
  ! "neutral / no data" sentinel; |Lo| >= Lo_neutral is also treated as
  ! neutral. In both cases the psi_m terms vanish.
  if (Lo == 0.0 .or. abs(Lo) >= Lo_neutral) then
     psi_z = 0.0
     psi_0 = 0.0
  else
     call psi_m_MOST((z - d) / Lo, psi_z)
     call psi_m_MOST(z0 / Lo, psi_0)
  endif

  ! Roughness-sublayer stability correction (HF2007 Eq. 15): phi_m(z'/Lo)
  ! inside the RSL integral. Returns 0 in neutral conditions.
  call psihat_stab_corr(z - d, arate, Lo, dpsihat)

  U = (1.0 / kappa) * ( log((z - d) / z0)          &
                        - psi_z + psi_0             &
                        + c1 * exp(-c2 / 2.0) * (E1 + dpsihat) )
end subroutine U_above


!============================================================
! SUBROUTINE: U_within
!
! Purpose:
! --------
! Compute the dimensionless velocity profile within the canopy.
!
! Inputs:
! -------
! z         : height above ground [m]
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
! h_PV      : canopy height [m]
! beta_deg  : panel tilt angle [deg]
! Lo        : Obukhov length [m] (<0 unstable, >0 stable;
!             0 or |Lo| very large -> neutral)
!
! Output:
! -------
! U         : dimensionless velocity [-]
!
! Note:
! -----
! Stability enters only through the canopy-top velocity U+(h_PV),
! computed with the stability-corrected above-canopy profile. The
! log-shape and exponential parts within the canopy are kept in
! their neutral form, consistent with the approach of applying the
! MOST correction above the canopy only (Harman & Finnigan 2007).
!============================================================
subroutine U_within(z, theta, lambda_f, h_PV, beta_deg, theta_0, Lo, U)
  implicit none

  real, intent(in)  :: z, theta, lambda_f, h_PV, beta_deg, theta_0
  real, intent(in)  :: Lo
  real, intent(out) :: U

  real :: U_h      ! Velocity at canopy top [-]
  real :: w        ! Weighting factor [-]
  real :: d_par    ! Displacement height (perpendicular) [m]
  real :: z0_par   ! Roughness length (perpendicular) [m]
  real :: log_part ! Logarithmic contribution [-]
  real :: exp_part ! Exponential contribution [-]
  real :: c3       ! Decay coefficient [-]

  ! Canopy-top velocity: carries the stability correction through Lo
  call U_above(h_PV, theta, lambda_f, h_PV, beta_deg, theta_0, Lo, U_h)
  call alpha_w(lambda_f, w)

  w = w * cos(theta - theta_0)**2

  call d_theta(theta_0 + 3.141592653589793 / 2.0, lambda_f, h_PV, theta_0, d_par)
  call z0_theta(theta_0 + 3.141592653589793 / 2.0, lambda_f, h_PV, beta_deg, theta_0, z0_par)
  call c3_theta(theta, lambda_f, theta_0, c3)

  log_part = log((z - d_par) / z0_par) / log((h_PV - d_par) / z0_par)
  exp_part = exp(-c3 * (1.0 - z / h_PV))

  U = (1.0 - w) * U_h * log_part + w * U_h * exp_part
end subroutine U_within

!============================================================
! SUBROUTINE: Rib_to_Lobk
!
! Purpose:
! --------
! Calculate L_obk.
!
! Inputs:
! -------
! z         : height above ground [m]
! theta     : wind direction [rad]
! lambda_f  : frontal density [-]
! h_PV      : canopy height [m]
! beta_deg  : panel tilt angle [deg]
! Lo        : Obukhov length [m] (<0 unstable, >0 stable;
!             0 or |Lo| very large -> neutral)
!
! Output:
! -------
! U         : dimensionless velocity [-]
!
! Note:
! -----
! Stability enters only through the canopy-top velocity U+(h_PV),
! computed with the stability-corrected above-canopy profile. The
! log-shape and exponential parts within the canopy are kept in
! their neutral form, consistent with the approach of applying the
! MOST correction above the canopy only (Harman & Finnigan 2007).
!============================================================
subroutine Rib_to_Lobk(nx,ny,nz,ibgn,iend,jbgn,jend,zp,roufns,wspd,ptsfc,       &
		       pt1, Lo)
		       
	  IMPLICIT NONE

	  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

	  INTEGER :: ibgn              ! i-index where evaluation begins.
	  INTEGER :: iend              ! i-index where evaluation ends.
	  INTEGER :: jbgn              ! j-index where evaluation begins.
	  INTEGER :: jend              ! j-index where evaluation ends.

	  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate
		                           ! defined at w-point of staggered grid.
	  REAL :: roufns(nx,ny)        ! Surface roughness length

	  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

	  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground
		                           ! level (K)
	  REAL :: pt1   (nx,ny)        ! Potential temperature at the
		                           ! 1st scalar
		                           ! point above ground level, (K)
	  REAL :: Lo   (nx,ny)         ! Obukhov length scale (m)

	!
	!-----------------------------------------------------------------------
	!
	!  Misc. local variables:
	!
	!-----------------------------------------------------------------------
	!
	  INTEGER :: i, j

	  REAL :: z1                   ! The height of 1st scalar point above
		                           ! ground level (m)
	  REAL :: dz0                  ! z1-roufns, where roufns is defined as
		                           ! surface roughness length
	  REAL :: bulkri               ! Bulk Richardson number

	  REAL :: stabp                ! Monin-Obukhov STABility Parameter
		                           ! (zeta)
	  REAL :: x1,x0,psim           ! Intermediate parameters needed
	  REAL :: z1drou,qb3pb2
	  REAL :: c7,c8
	  REAL :: sb,qb,pb,thetab,tb   ! During  computations
	  REAL :: a,b,c,d
	  REAL :: tempan,sqrtqb

	  REAL :: zt                   ! Thermal roughness length
	  REAL :: dzt
	  REAL :: ztdrou

	  REAL :: z1droup, ztdroup

	  REAL, PARAMETER :: eps_tol = 1.0d-6
	  
	!
	!-----------------------------------------------------------------------
	!
	!  Include files:
	!
	!-----------------------------------------------------------------------
	!
	  INCLUDE 'phycst.inc'
	  INCLUDE 'globcst.inc'
	  INCLUDE 'sfcphycst.inc'
	  
	!
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
	!  Beginning of executable code...
	!
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
	!
	!-----------------------------------------------------------------------
	!
	!  Following Byun (1990).
	!
	!-----------------------------------------------------------------------
	!
	  DO j=jbgn,jend
		DO i=ibgn,iend

		  z1  = 0.5d0*(zp(i,j,3)-zp(i,j,2))
		  z1  = MIN(z1,z1limit)

		  zt = ztz0*roufns(i,j)           ! Original Zt formulation

		  dzt = z1-zt
		  ztdrou = z1/zt
		  ztdroup = (z1+zt)/zt

		  dz0 = z1-roufns(i,j)
		  z1drou = z1/roufns(i,j)
		  z1droup = (z1+roufns(i,j))/roufns(i,j)

		  bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
		           (wspd(i,j)*wspd(i,j))

		  IF (bulkri <= 0.0d0) THEN
	!
	!-----------------------------------------------------------------------
	!
	!  Unstable case: See equations (28)-(34) in Byun (1990).
	!
	!-----------------------------------------------------------------------
	!
		    bulkri = MAX (bulkri,-10.0d0)

		    sb =bulkri/prantl0l

		    qb=oned9*(c1l+c2l*sb*sb)
		    pb=oned54*(c3l+c4l*sb*sb)

		    qb3pb2=qb**3-pb*pb
		    c7 = (z1*dzt*LOG(z1drou)*LOG(z1drou))/(dz0*dz0*LOG(ztdroup))

		    IF( qb3pb2 >= 0.0d0 ) THEN

		      sqrtqb = SQRT(qb)
		      tempan = MAX( -1.0d0, MIN( 1.0d0, pb/(sqrtqb**3) ) )

		      thetab=ACOS(tempan)
		      stabp =c7*(-2.0d0*sqrtqb*COS(thetab/3.0d0)+c5l)

		    ELSE

		      tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
		      stabp =c7*(-(tb+qb/tb)+c5l)

		    END IF
		    
		  ELSE
	!
	!-----------------------------------------------------------------------
	!
	!  Stable case: See equation (19) in Byun (1990).
	!
	!-----------------------------------------------------------------------
	!
		  sb =bulkri/prantl0l

          c7 = (z1*dzt*LOG(z1drou)*LOG(z1drou))/(dz0*dz0*LOG(ztdroup))  
          stabp = c7/(2.*betahl*(betaml*bulkri-1.)) * &
               (-(2.*betahl*bulkri-1.)+sqrt(1.+4.*(betahl-betaml)*sb))
          
		  END IF

		  Lo(i,j) = z1/stabp

		END DO
	  END DO

	  RETURN
		    
end subroutine Rib_to_Lobk
