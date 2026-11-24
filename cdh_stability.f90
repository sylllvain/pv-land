!##################################################################
!######                                                      ######
!######            SUBROUTINE CDH_STABILITY                  ######
!######                                                      ######
!######   Compute stability-dependent drag coefficients      ######
!######   for heat and momentum using PV-LAND conventions    ######
!######   (kappa, Pr as input parameters).                   ######
!######                                                      ######
!######   Unstable: Businger-Dyer / Paulson (1970)           ######
!######   Stable:   Louis et al. (1981)                      ######
!######                                                      ######
!##################################################################

SUBROUTINE cdh_stability(z, z0, T_sfc, T_air, U_air, &
                         kappa, Pr, cdh, cdm)
!
!-----------------------------------------------------------------------
!  PURPOSE:
!    Compute Cdh = Cu * Ctheta and Cdm = Cu * Cu with Monin-Obukhov
!    stability corrections, using the caller's kappa and Pr values.
!
!  INPUT:
!    z       : interface height (m), e.g. z_can
!    z0      : roughness length (m)
!    T_sfc   : surface temperature (K), e.g. T_nat
!    T_air   : air temperature at height z (K), e.g. T_can
!    U_air   : wind speed at height z (m/s), e.g. U_can
!    kappa   : von Karman constant (0.387 for PV-LAND)
!    Pr      : turbulent Prandtl number (0.71 for PV-LAND)
!
!  OUTPUT:
!    cdh     : drag coefficient for heat = Cu * Ctheta
!    cdm     : drag coefficient for momentum = Cu * Cu
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL, INTENT(IN)  :: z, z0, T_sfc, T_air, U_air
  REAL, INTENT(IN)  :: kappa, Pr
  REAL, INTENT(OUT) :: cdh, cdm

  ! Local
  REAL :: ln_z_z0        ! ln(z/z0)
  REAL :: Cu_n, Cth_n    ! Neutral Cu and Ctheta
  REAL :: Ri_b           ! Bulk Richardson number
  REAL :: U_eff          ! Effective wind speed (with minimum)
  REAL :: psim, psih     ! Stability correction functions
  REAL :: Cu, Ctheta     ! Stability-corrected values
  REAL :: x              ! Intermediate for unstable case
  REAL :: f_m, f_h       ! Stability functions for stable case

  REAL, PARAMETER :: g = 9.81
  REAL, PARAMETER :: gamma_m = 16.0   ! Businger-Dyer constant for momentum
  REAL, PARAMETER :: gamma_h = 16.0   ! Businger-Dyer constant for heat
  REAL, PARAMETER :: Ri_b_min = -2.0  ! Lower bound for Ri_b (strongly unstable)
  REAL, PARAMETER :: U_min = 0.1      ! Minimum wind speed (m/s)
  REAL, PARAMETER :: Cu_max_factor = 4   ! Max Cu/Cu_neutral (unstable limit)
  REAL, PARAMETER :: Cth_max_factor = 6  ! Max Ctheta/Ctheta_neutral (unstable limit)

  ! Neutral values
  ln_z_z0 = LOG(z / z0)
  Cu_n   = kappa / ln_z_z0
  Cth_n  = kappa / (Pr * ln_z_z0)

  ! Effective wind speed
  U_eff = MAX(U_air, U_min)

  ! Bulk Richardson number: Ri_b = (g/T_air) * (T_air - T_sfc) * z / U^2
  ! Note: T_air - T_sfc so that Ri_b < 0 when surface is warmer (unstable)
  Ri_b = (g / T_air) * (T_air - T_sfc) * (z - z0) / (U_eff * U_eff)

  IF (Ri_b < 0.0) THEN
    !-----------------------------------------------------------------
    ! UNSTABLE CASE (T_sfc > T_air): Businger-Dyer + Paulson (1970)
    !
    ! zeta = z/L is related to Ri_b. For the flux-profile relations:
    !   phi_m = (1 - gamma_m * zeta)^(-1/4)
    !   phi_h = (1 - gamma_h * zeta)^(-1/2)
    !
    ! Integrated psi functions (Paulson 1970):
    !   psi_m = 2*ln((1+x)/2) + ln((1+x^2)/2) - 2*atan(x) + pi/2
    !   psi_h = 2*ln((1+y)/2)
    ! where x = (1-gamma_m*zeta)^(1/4), y = (1-gamma_h*zeta)^(1/2)
    !
    ! For simplicity, we use the approximation zeta ≈ Ri_b for 
    ! weakly unstable conditions, valid for |Ri_b| < ~0.5.
    !-----------------------------------------------------------------

    Ri_b = MAX(Ri_b, Ri_b_min)

    ! Momentum correction
    x = (1.0 - gamma_m * Ri_b)**0.25
    psim = 2.0*LOG((1.0+x)/2.0) + LOG((1.0+x*x)/2.0) - 2.0*ATAN(x) + 1.5707963

    ! Heat correction 
    x = (1.0 - gamma_h * Ri_b)**0.5
    psih = 2.0*LOG((1.0+x)/2.0)

    ! Corrected coefficients
    Cu     = kappa / (ln_z_z0 - psim)
    Ctheta = kappa / (Pr * (ln_z_z0 - psih))

    ! Apply limiters (similar to ARPS: max 2x for Cu, 3.33x for Ctheta)
    Cu     = MIN(Cu,     Cu_max_factor  * Cu_n)
    Ctheta = MIN(Ctheta, Cth_max_factor * Cth_n)

    ! Ensure positive
    Cu     = MAX(Cu,     0.1 * Cu_n)
    Ctheta = MAX(Ctheta, 0.1 * Cth_n)

  ELSE
    !-----------------------------------------------------------------
    ! STABLE CASE (T_sfc < T_air): Louis et al. (1981)
    !
    !   Cu     = Cu_n / sqrt(1 + 2*b*Ri_b / sqrt(1+d*Ri_b))
    !   Ctheta = Cth_n * sqrt(1 + 2*b*Ri_b / sqrt(1+d*Ri_b))
    !            / (Pr * (1 + 3*b*Ri_b * sqrt(1+d*Ri_b)))
    !
    ! Simplified: use common form with b=d=5
    !-----------------------------------------------------------------

    f_m = 1.0 / SQRT(1.0 + 2.0*5.0*Ri_b / SQRT(1.0+5.0*Ri_b))
    f_h = SQRT(1.0 + 2.0*5.0*Ri_b / SQRT(1.0+5.0*Ri_b)) &
        / (1.0 + 3.0*5.0*Ri_b * SQRT(1.0+5.0*Ri_b))

    Cu     = Cu_n * f_m
    Ctheta = Cth_n * f_h

    ! Apply lower limit (blimit ~ 0.1 in ARPS)
    Cu     = MAX(Cu,     0.1 * Cu_n)
    Ctheta = MAX(Ctheta, 0.1 * Cth_n)

  END IF

  ! Final drag coefficients
  cdm = Cu * Cu
  cdh = Cu * Ctheta

END SUBROUTINE cdh_stability
