!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE OBUKHOV_UPDATE               ######
!######                                                      ######
!######   Diagnose the Obukhov length from the surface       ######
!######   fluxes of the PREVIOUS small timestep and relax    ######
!######   the model state towards it.                        ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE obukhov_update(ustar, H_sfc, LE_sfc,                       &
                          T_air, qv_air, psfc, rho_air,               &
                          cp_in, lv_in, p0_in, rddcp_in,              &
                          zmd_ref, dt, tau_L, L_floor, L_obk)

!============================================================
! Purpose:
! --------
! Compute the Monin-Obukhov length from the park-scale friction
! velocity and the surface buoyancy flux (L is a model state,
! not a forcing).
!
!     L = - u*^3 * theta_v / ( kappa * g * (w'theta_v')_s )
!
! with the surface virtual heat flux written from the sensible
! and latent heat fluxes:
!
!     (w'theta_v')_s = H/(rho*cp) + 0.61 * theta * LE/(rho*Lv)
!
! Sign convention (identical to profil.f90):
!     L < 0  unstable   (upward buoyancy flux, H_v > 0)
!     L > 0  stable
!     L = 0  sentinel meaning "neutral / no stability correction"
!
! The routine is meant to be called ONCE PER SMALL TIMESTEP with
! the fluxes just computed, i.e. L is lagged by one small step
! (dtsfc = 2 s here). This is the usual explicit treatment in
! land-surface schemes: u* depends on L through the wind profile
! and L depends on u* through the fluxes, so the pair is solved
! by lagging rather than by an inner iteration. With dtsfc << the
! adjustment time of the surface layer the lag is immaterial.
!
! Numerical safeguards (all of them matter in practice):
!   1. the state that is relaxed is 1/L, NOT L: the Obukhov length
!      passes through +-infinity at neutral, so relaxing L itself
!      would be meaningless near sunrise/sunset;
!   2. weak wind / vanishing buoyancy flux fall back to neutral;
!   3. a floor on |L| that is a pure numerical guard, NOT a physical
!      limiter. The validity range of the stability functions is
!      enforced where those functions are evaluated (parameter
!      zeta_max in psi_m_MOST and psihat_stab_corr, profil.f90), not
!      here, so that there is a single place in the code that decides
!      how far the similarity theory is trusted.
!
! Inputs:
! -------
! ustar    : friction velocity at the reference level [m/s]
! H_sfc    : surface sensible heat flux, POSITIVE UPWARD [W/m2]
! LE_sfc   : surface latent heat flux, POSITIVE UPWARD [W/m2]
! T_air    : air temperature at the reference level [K]
! qv_air   : specific humidity at the reference level [kg/kg]
! psfc     : air pressure [Pa]
! rho_air  : air density [kg/m3]
! cp_in    : specific heat of dry air at constant pressure [J/kg/K]
! lv_in    : latent heat of vaporisation [J/kg]
! p0_in    : reference pressure for the potential temperature [Pa]
! rddcp_in : Rd/cp [-]
! zmd_ref  : z_ref - d, height above the virtual origin [m]. Kept in the
!            interface for diagnostics (zeta = zmd_ref/L) but not used to
!            limit L; see the note on safeguard 3 above.
! dt       : small timestep [s]
! tau_L    : relaxation time constant of 1/L [s] (0 => no relaxation)
!
! In/out:
! -------
! L_obk    : Obukhov length [m]; on entry the previous value
!            (0 = neutral), on exit the updated one.
!============================================================

  implicit none

  real, intent(in)    :: ustar, H_sfc, LE_sfc
  real, intent(in)    :: T_air, qv_air, psfc, rho_air
  real, intent(in)    :: cp_in, lv_in, p0_in, rddcp_in
  real, intent(in)    :: zmd_ref, dt, tau_L
  real, intent(in)    :: L_floor   ! numerical floor on |L| [m] (namelist pv_L_floor)
  real, intent(inout) :: L_obk

  !-------------------------
  ! Constants
  !-------------------------
  ! von Karman constant: MUST be the same as the one used in the wind
  ! profile (profil.f90, subroutine U_above), otherwise u* and L are not
  ! consistent with each other.
  real, parameter :: kappa      = 0.387
  real, parameter :: grav       = 9.81      ! [m/s2]
  ! Above this magnitude L is treated as neutral (same value as profil.f90).
  real, parameter :: Lo_neutral = 1.0e5     ! [m]
  ! Numerical floor on |L| [m]. This is only there to keep the arguments of
  ! the stability functions finite; it is not a validity limit (that is
  ! zeta_max in profil.f90). It binds only at very low wind speed.
  ! L_floor is an argument (namelist pv_L_floor, default 0.5 m).  [SYLVAIN]
  ! Below this friction velocity the surface layer is not in a MOST regime.
  real, parameter :: ustar_min  = 0.02      ! [m/s]
  ! Below this buoyancy flux the stratification is treated as neutral.
  real, parameter :: wtv_min    = 1.0e-5    ! [K m/s]

  !-------------------------
  ! Local variables
  !-------------------------
  real :: theta      ! potential temperature at the reference level [K]
  real :: thetav     ! virtual potential temperature [K]
  real :: wtv        ! surface virtual heat flux (w'theta_v')_s [K m/s]
  real :: invL_old   ! 1/L at the previous small timestep [1/m]
  real :: invL_new   ! 1/L diagnosed from the current fluxes [1/m]
  real :: invL       ! relaxed 1/L [1/m]
  real :: invL_cap   ! numerical cap on |1/L| [1/m]
  real :: wrel       ! relaxation weight [-]

  !------------------------------------------------------------
  ! 1. Previous state, as an inverse length (0 = neutral)
  !------------------------------------------------------------
  if (L_obk == 0.0 .or. abs(L_obk) >= Lo_neutral) then
     invL_old = 0.0
  else
     invL_old = 1.0 / L_obk
  endif

  !------------------------------------------------------------
  ! 2. Diagnosed 1/L from the fluxes of the previous small step
  !------------------------------------------------------------
  theta  = T_air * (p0_in / psfc)**rddcp_in
  thetav = theta * (1.0 + 0.61 * qv_air)

  ! Virtual (buoyancy) heat flux, positive upward.
  ! The moisture term is secondary here (about 7% of LE in K m/s units
  ! against 100% of H) but it is what makes L behave correctly over a
  ! well-watered surface at low Bowen ratio.
  wtv = H_sfc / (rho_air * cp_in)                                       &
      + 0.61 * theta * LE_sfc / (rho_air * lv_in)

  if (ustar <= ustar_min .or. abs(wtv) <= wtv_min) then
     ! No usable similarity regime: relax towards neutral rather than
     ! dividing by something that is going to zero.
     invL_new = 0.0
  else
     ! L = -u*^3 thetav / (kappa g wtv)  ->  1/L = -kappa g wtv / (u*^3 thetav)
     invL_new = - kappa * grav * wtv / (ustar**3 * thetav)
  endif

  !------------------------------------------------------------
  ! 3. Relaxation on 1/L
  !------------------------------------------------------------
  if (tau_L > 0.0) then
     wrel = dt / (tau_L + dt)
  else
     wrel = 1.0
  endif
  invL = (1.0 - wrel) * invL_old + wrel * invL_new

  !------------------------------------------------------------
  ! 4. Numerical floor on |L| (see the note in the header)
  !------------------------------------------------------------
  invL_cap = 1.0 / L_floor
  if (invL >  invL_cap) invL =  invL_cap
  if (invL < -invL_cap) invL = -invL_cap

  !------------------------------------------------------------
  ! 5. Back to a length, with the neutral sentinel of profil.f90
  !------------------------------------------------------------
  if (abs(invL) <= 1.0 / Lo_neutral) then
     L_obk = 0.0
  else
     L_obk = 1.0 / invL
  endif

end subroutine obukhov_update
