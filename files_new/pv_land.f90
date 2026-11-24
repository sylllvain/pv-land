!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE PV_LAND                      ######
!######                                                      ######
!######   Main PV-LAND subroutine with RK2 for PV nodes      ######
!######   Uses soilebm_frc, sfcflx, solve_air_node_temps     ######
!######                                                      ######
!######   COMBINED version: outputs both PV diagnostics      ######
!######   AND nat surface irradiances/radiosities            ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE pv_land(nx, ny, nz, nzsoil, nstyps, &
                  L, a, beta_deg, h_PV, h_clear, theta_0, &
                  z0_nat, z0_PV_back, z0_PV_front, albedo_nat, &
				  albedo_PVf, albedo_PVb, &
                  emiss_nat, emiss_PVf, emiss_PVb, &
                  urefi, udiri, trefi, qvrefi, prcrate, L_obk_in, &
                  swouti, swini, lwouti, lwini, prefi, rhorefi, &
                  phii, lambdai, S_diffusei, &
                  T_nat, T_soil, T_can, T_PVa, &
                  T_PVf, T_cell, T_PVb, &
                  qsoil_surf, qsoil_deep, wetcanp, qvsfc, &
                  soiltyp, vegtyp, lai, veg, stypfrct, &
                  H_ref_out, LE_ref_out, Rn_ref_out, G_out, &
                  SW_out_sky, LW_out_sky, H_nat_out, &
                  H_PVa_ref_out, H_can_ref_out, &
                  G_SW_nat_out, J_SW_nat_out, &
                  G_LW_nat_out, J_LW_nat_out, &
                  diag_qLW_PVf, diag_qLW_PVb, &
                  diag_H_PVf, diag_H_PVb, &
                  diag_qcond_PVf, diag_qcond_PVb, &
                  diag_qSW_PVf, diag_qSW_PVb, &
                  diag_qPV, &
                  diag_cdh_can_out, diag_cdh_PVa_out, &
                  diag_T_can_out, diag_T_PVa_out, &
                  diag_U_can_out, diag_U_ref_out, &
                  diag_JLW_PVf_out, diag_JLW_PVb_out)

  IMPLICIT NONE
  
  INCLUDE 'timelvls.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'
  INCLUDE 'soilcst.inc'
  INCLUDE 'grid.inc'
  
  !=====================================================================
  ! INPUT PARAMETERS
  !=====================================================================
  INTEGER, INTENT(IN) :: nx, ny, nz, nzsoil, nstyps
  REAL, INTENT(IN)    :: L, a, beta_deg, h_PV, h_clear, theta_0
  REAL, INTENT(IN)    :: z0_nat, z0_PV_back, z0_PV_front
  REAL, INTENT(IN)    :: albedo_nat, albedo_PVf, albedo_PVb
  REAL, INTENT(IN)    :: emiss_nat, emiss_PVf, emiss_PVb
  
  ! Interpolation arrays (2 time levels)
  REAL, INTENT(IN)    :: urefi(nx,ny,2), udiri(nx,ny,2)
  REAL, INTENT(IN)    :: trefi(nx,ny,2), qvrefi(nx,ny,2)
  REAL, INTENT(IN)    :: prcrate(nx,ny)
  REAL, INTENT(IN)    :: L_obk_in(nx,ny)   ! Obukhov length [m] (L>0 stable, L<0 unstable)
  REAL, INTENT(IN)    :: swouti(nx,ny,2), swini(nx,ny,2)
  REAL, INTENT(IN)    :: lwouti(nx,ny,2), lwini(nx,ny,2)
  REAL, INTENT(IN)    :: prefi(nx,ny,2), rhorefi(nx,ny,2)
  REAL, INTENT(IN)    :: phii(nx,ny,2), lambdai(nx,ny,2), S_diffusei(nx,ny,2)
  
  ! Surface properties
  INTEGER, INTENT(IN) :: soiltyp(nx,ny,nstyps), vegtyp(nx,ny)
  REAL, INTENT(IN)    :: lai(nx,ny), veg(nx,ny)
  REAL, INTENT(IN)    :: stypfrct(nx,ny,nstyps)
  
  !=====================================================================
  ! STATE VARIABLES (INPUT/OUTPUT)
  !=====================================================================
  REAL, INTENT(INOUT) :: T_nat(nx,ny), T_soil(nx,ny)
  REAL, INTENT(INOUT) :: T_can(nx,ny), T_PVa(nx,ny)
  REAL, INTENT(INOUT) :: T_PVf(nx,ny), T_cell(nx,ny), T_PVb(nx,ny)
  REAL, INTENT(INOUT) :: qsoil_surf(nx,ny), qsoil_deep(nx,ny)
  REAL, INTENT(INOUT) :: wetcanp(nx,ny,0:nstyps)
  REAL, INTENT(INOUT) :: qvsfc(nx,ny,0:nstyps)
  
  !=====================================================================
  ! OUTPUT
  !=====================================================================
  REAL, INTENT(OUT)   :: H_ref_out(nx,ny), LE_ref_out(nx,ny)
  REAL, INTENT(OUT)   :: Rn_ref_out(nx,ny), G_out(nx,ny)
  REAL, INTENT(OUT)   :: SW_out_sky(nx,ny), LW_out_sky(nx,ny)
  REAL, INTENT(OUT)   :: H_nat_out(nx,ny)
  
  ! Nat surface diagnostics (from _nat version)
  REAL, INTENT(OUT)   :: H_PVa_ref_out(nx,ny), H_can_ref_out(nx,ny)
  REAL, INTENT(OUT)   :: G_SW_nat_out(nx,ny), J_SW_nat_out(nx,ny)
  REAL, INTENT(OUT)   :: G_LW_nat_out(nx,ny), J_LW_nat_out(nx,ny)
  
  ! PV energy balance diagnostics (from _pv version)
  REAL, INTENT(OUT)   :: diag_qLW_PVf(nx,ny), diag_qLW_PVb(nx,ny)
  REAL, INTENT(OUT)   :: diag_H_PVf(nx,ny), diag_H_PVb(nx,ny)
  REAL, INTENT(OUT)   :: diag_qcond_PVf(nx,ny), diag_qcond_PVb(nx,ny)
  REAL, INTENT(OUT)   :: diag_qSW_PVf(nx,ny), diag_qSW_PVb(nx,ny)
  REAL, INTENT(OUT)   :: diag_qPV(nx,ny)
  REAL, INTENT(OUT)   :: diag_cdh_can_out(nx,ny), diag_cdh_PVa_out(nx,ny)
  REAL, INTENT(OUT)   :: diag_T_can_out(nx,ny), diag_T_PVa_out(nx,ny)
  REAL, INTENT(OUT)   :: diag_U_can_out(nx,ny), diag_U_ref_out(nx,ny)
  REAL, INTENT(OUT)   :: diag_JLW_PVf_out(nx,ny), diag_JLW_PVb_out(nx,ny)
  
  !=====================================================================
  ! LOCAL VARIABLES
  !=====================================================================
  
  ! Time stepping
  INTEGER :: it, i, j
  REAL :: dtsfc2
  REAL, PARAMETER :: pi = 3.14159265359
  
  ! Geometry
  REAL :: beta, w, z0_park, zd_park
  REAL :: lambda_f_park, z0_park_arr(nx,ny)
  
  ! Interpolated forcing
  REAL :: uref(nx,ny), udir(nx,ny), tref(nx,ny), qvref(nx,ny)
  REAL :: swout(nx,ny), swin(nx,ny), lwout(nx,ny), lwin(nx,ny)
  REAL :: psfc(nx,ny), rhosfc(nx,ny)
  REAL :: phi(nx,ny), lambda(nx,ny), S_diffuse(nx,ny)
  
  ! Wind profile
  REAL :: U_can(nx,ny), U_PVa(nx,ny), U_ref(nx,ny)
  REAL :: u_star_ref, Uplus
  REAL :: L_obk            ! Obukhov length used in the wind profile [m]
  REAL :: z_can, z_PVa, z_ref
  
  ! Drag coefficients (calculated with current temperatures)
  REAL :: cdm_nat(nx,ny), cdh_nat(nx,ny), cdq_nat(nx,ny)
  REAL :: cdm_PVb(nx,ny), cdh_PVb(nx,ny), cdq_PVb(nx,ny)
  REAL :: cdm_PVf(nx,ny), cdh_PVf(nx,ny), cdq_PVf(nx,ny)
  REAL :: cdm_can(nx,ny), cdh_can(nx,ny), cdq_can(nx,ny)
  REAL :: cdm_PVa(nx,ny), cdh_PVa(nx,ny), cdq_PVa(nx,ny)
  
  ! PV thermal parameters
  REAL :: delta_PVf, delta_cell, delta_PVb
  REAL :: rho_cp_PVf, rho_cp_cell, rho_cp_PVb
  REAL :: k_glass, k_EVA, k_backsheet
  REAL :: delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet
  REAL :: eta_STC, gamma_PV
  
  ! RK2 for ALL nodes (nat, soil, PV)
  REAL :: frc_nat(nx,ny), frc_soil(nx,ny)
  REAL :: frc_PVf(nx,ny), frc_cell(nx,ny), frc_PVb(nx,ny)
  REAL :: rhs_nat(nx,ny), rhs_soil(nx,ny)
  REAL :: rhs_PVf(nx,ny), rhs_cell(nx,ny), rhs_PVb(nx,ny)
  REAL :: T_nat_half(nx,ny), T_soil_half(nx,ny)
  REAL :: T_PVf_half(nx,ny), T_cell_half(nx,ny), T_PVb_half(nx,ny)
  
  ! RK2 for soil moisture
  REAL :: frc_wsfc(nx,ny), frc_wdp(nx,ny), frc_wcnp(nx,ny)
  REAL :: rhs_wsfc(nx,ny), rhs_wdp(nx,ny), rhs_wcnp(nx,ny)
  REAL :: qsoil_surf_half(nx,ny), qsoil_deep_half(nx,ny)
  REAL :: wetcanp_half(nx,ny,0:nstyps)
  REAL :: wrmax(nx,ny)
  
  ! Diagnostic outputs from compute_all_forcings
  REAL :: H_nat(nx,ny), LE_nat(nx,ny), G_nat(nx,ny)
  REAL :: SW_out_sky_loc(nx,ny), LW_out_sky_loc(nx,ny)
  
  ! Dummy diagnostic arrays for RK2 internal calls
  REAL :: dum_qLW_PVf(nx,ny), dum_qLW_PVb(nx,ny)
  REAL :: dum_H_PVf(nx,ny), dum_H_PVb(nx,ny)
  REAL :: dum_qcond_PVf(nx,ny), dum_qcond_PVb(nx,ny)
  REAL :: dum_qSW_PVf(nx,ny), dum_qSW_PVb(nx,ny)
  REAL :: dum_qPV(nx,ny)
  REAL :: dum_cdh_can(nx,ny), dum_cdh_PVa(nx,ny)
  REAL :: dum_T_can(nx,ny), dum_T_PVa(nx,ny)
  REAL :: dum_JLW_PVf(nx,ny), dum_JLW_PVb(nx,ny)
  REAL :: dum_G_SW_nat(nx,ny), dum_J_SW_nat(nx,ny)
  REAL :: dum_G_LW_nat(nx,ny), dum_J_LW_nat(nx,ny)
  
  ! Atmospheric arrays for dragcoef/sfcflx
  REAL :: zp(nx,ny,nz), zpsoil(nx,ny,nzsoil)
  REAL :: u_wind(nx,ny,nz), v_wind(nx,ny,nz), w_wind(nx,ny,nz)
  REAL :: ptprt(nx,ny,nz), pprt(nx,ny,nz), qv(nx,ny,nz)
  REAL :: ptbar(nx,ny,nz), pbar(nx,ny,nz), qvbar(nx,ny,nz)
  REAL :: rhostr(nx,ny,nz)
  REAL :: tair(nx,ny), qvair(nx,ny)
  
  ! Dummy soil type for dragcoef
  INTEGER :: soiltyp_dummy(nx,ny)
  
  ! Snow depth
  REAL :: snowdpth(nx,ny)
  
  ! For flux calculations
  REAL :: pt_PVa(nx,ny), pt_can(nx,ny), pt_nat(nx,ny)
  REAL :: pt_PVb(nx,ny), pt_PVf(nx,ny), z0_sfc_arr(nx,ny)   ! for dragcoef stability test
  REAL :: wdspsfc(nx,ny)
  REAL, PARAMETER :: Pr = 0.71
  REAL :: l_nat_can       ! Interface length: natural surface -> canopy air (= z_can)
  REAL :: l_PVb_can       ! Interface length: PV back -> canopy air
  REAL :: l_PVf_PVa       ! Interface length: PV front -> air above PV
  
  ! Dedicated zp array for the can->ref, PVa->ref, and nat->can dragcoef calls.
  REAL :: zp_dc(nx,ny,nz)
  REAL :: l_can_ref, l_PVa_ref  ! Interface lengths
  
  REAL :: zp_dc_nat(nx,ny,nz)
  
  REAL, PARAMETER :: U_min_fac = 0.2   ! wind floor (m/s) to prevent divergence
  !===================================================================
  ! PV panel convection
  ! Forced: Sparrow, Ramsey & Mass (1979), prefactor calibrated on site
  ! Free  : Kaplani & Kaplanis (2014)
  ! Mixing: Kaplani Eq.(20), m = 3
  !===================================================================
  ! Forced-convection prefactors, calibrated on the site energy balance (same value both faces).
  REAL, PARAMETER :: C_sparrow_f  = 5
  REAL, PARAMETER :: C_sparrow_b  = 5
  ! Wind exponent, different for each face: the front is a free plate in the flow, the back is
  ! a confined cavity between panel and soil and responds far less to the outside wind.
  REAL, PARAMETER :: p_forced_f = 0.5   ! front face: Sparrow et al. (1979)
  REAL, PARAMETER :: p_forced_b = 0.30  ! back face: measured on site (confined cavity)
  ! Reference Reynolds used to normalise the back-face law, so that both faces share the
  ! same C_sparrow and differ ONLY through the exponent. Re_0 = median back-face Re (on U_can).
  REAL, PARAMETER :: Re_0 = 2.0e5
  REAL :: L_star     ! [m]      Sparrow characteristic length, L* = 4*A/C
  REAL :: Tfilm      ! [K]      film temperature, 0.5*(T_surface + T_air)
  REAL :: kair       ! [W/m/K]  air thermal conductivity at Tfilm
  REAL :: muair      ! [Pa.s]   air dynamic viscosity at Tfilm (Sutherland)
  REAL :: nuair      ! [m2/s]   air kinematic viscosity, mu/rho
  REAL :: alair      ! [m2/s]   air thermal diffusivity, k/(rho*cp)
  REAL :: dTs        ! [K]      signed surface-to-air temperature difference
  REAL :: Re_pv      ! [-]      Reynolds number on L*
  REAL :: Nu_fo      ! [-]      forced-convection Nusselt number
  REAL :: h_cv       ! [W/m2/K] convective heat transfer coefficient
  REAL :: cc_pr      ! [-]      Prandtl function of Kaplani Eq.(8), precomputed once
  REAL :: Tf_na      ! [K]      Kaplani expansion temperature, Ta + 0.25*(Ts-Ta)
  REAL :: Ra_na      ! [-]      Rayleigh number on the physical length L (Kaplani Eq.(6))
  REAL :: theta_na   ! [rad]    plate angle from the vertical (Kaplani convention)
  REAL :: Grc_na     ! [-]      critical Grashof number (Kaplani, laminar-turbulent transition)
  REAL :: Rac_na     ! [-]      critical Rayleigh number, Grc*Pr
  REAL :: Nu_nat     ! [-]      natural-convection Nusselt (Kaplani, on L)
  REAL :: h_free     ! [W/m2/K] free (natural) convection coefficient
  REAL :: h_forced   ! [W/m2/K] forced-convection coefficient (Sparrow, on L_star)
  LOGICAL :: unstab  ! [-]      .TRUE. if the buoyancy arrangement is unstable
  
  REAL :: facteur
  
  !=====================================================================
  ! INITIALIZATION
  !=====================================================================
  
  nsfcst = 1800/dtsfc
  dtsfc2 = dtsfc / 2.0
  
  ! Geometry
  beta = beta_deg * pi / 180.0
  ! Sparrow characteristic length: L* = 4*A/C (plate area / perimeter)
  L_star = 2.0 * L
  ! Prandtl function of the Kaplani Eq.(8) branch
  cc_pr = ( 1.0 + (0.492/Pr)**(9.0/16.0) )**(8.0/27.0)
  w = L * COS(beta)
  
  ! Park-scale roughness: compute lambda_f once (geometry-only)
  CALL compute_lambda_f(a, L, beta, lambda_f_park)
  
  ! Grid setup
  zp(:,:,1) = -2.0 * 3.3 * h_PV
  zp(:,:,2) = 0.0
  zp(:,:,3) = 2.0 * 3.3 * h_PV
  zpsoil(:,:,1) = 0.1
  zpsoil(:,:,2) = 10
  
  ! Heights for wind profile
  z_ref = 0.5 * zp(1,1,3)
  z_can = (h_PV + h_clear) / 4.0

  ! Interface lengths for neutral drag coefficient calculations
  l_nat_can  = z_can       ! natural ground to canopy air
  l_PVb_can  = 0.015        ! PV back to canopy air
  l_PVf_PVa  = 0.015         ! PV front to air above PV

  ! z_PVa derived from l_PVa_ref = z_ref - z_PVa
  ! where z_PVa = (h_PV+h_clear)/2 + l_PVf*cos(beta)
  z_PVa = (h_PV + h_clear) / 2.0 + l_PVf_PVa * COS(beta)
  
  ! Get PV thermal parameters
  CALL pv_thermal_params(delta_PVf, delta_cell, delta_PVb, &
                        rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                        k_glass, k_EVA, k_backsheet, &
                        delta_glass, delta_EVA_top, &
                        delta_EVA_bot, delta_backsheet, &
                        eta_STC, gamma_PV)
  
  ! Initialize arrays
  snowdpth = 0.0
  soiltyp_dummy = 11
  ptprt = 0.0
  pprt = 0.0
  
  !=====================================================================
  ! SMALL TIMESTEP LOOP (30 seconds)
  !=====================================================================
  
  DO it = 1, nsfcst
  
	!print *, "Start small loop, it = ", it
    !-------------------------------------------------------------------
    ! 1. LINEAR INTERPOLATION OF FORCING
    !-------------------------------------------------------------------
    DO j = 1, ny
      DO i = 1, nx
        uref(i,j)   = urefi(i,j,1)   + (urefi(i,j,2)   - urefi(i,j,1))   * REAL(it) / REAL(nsfcst)
        udir(i,j)   = udiri(i,j,1)   + (udiri(i,j,2)   - udiri(i,j,1))   * REAL(it) / REAL(nsfcst)
        tref(i,j)   = trefi(i,j,1)   + (trefi(i,j,2)   - trefi(i,j,1))   * REAL(it) / REAL(nsfcst)
        qvref(i,j)  = qvrefi(i,j,1)  + (qvrefi(i,j,2)  - qvrefi(i,j,1))  * REAL(it) / REAL(nsfcst)
        swout(i,j)  = swouti(i,j,1)  + (swouti(i,j,2)  - swouti(i,j,1))  * REAL(it) / REAL(nsfcst)
        swin(i,j)   = swini(i,j,1)   + (swini(i,j,2)   - swini(i,j,1))   * REAL(it) / REAL(nsfcst)
        lwout(i,j)  = lwouti(i,j,1)  + (lwouti(i,j,2)  - lwouti(i,j,1))  * REAL(it) / REAL(nsfcst)
        lwin(i,j)   = lwini(i,j,1)   + (lwini(i,j,2)   - lwini(i,j,1))   * REAL(it) / REAL(nsfcst)
        psfc(i,j)   = prefi(i,j,1)   + (prefi(i,j,2)   - prefi(i,j,1))   * REAL(it) / REAL(nsfcst)
        rhosfc(i,j) = rhorefi(i,j,1) + (rhorefi(i,j,2) - rhorefi(i,j,1)) * REAL(it) / REAL(nsfcst)
        phi(i,j)    = phii(i,j,1)    + (phii(i,j,2)    - phii(i,j,1))    * REAL(it) / REAL(nsfcst)
        lambda(i,j) = lambdai(i,j,1) + (lambdai(i,j,2) - lambdai(i,j,1)) * REAL(it) / REAL(nsfcst)
        S_diffuse(i,j) = S_diffusei(i,j,1) + (S_diffusei(i,j,2) - S_diffusei(i,j,1)) * REAL(it) / REAL(nsfcst)
      END DO
    END DO
    
    U_ref = uref
    
    !-------------------------------------------------------------------
    ! 1b. COMPUTE PARK-SCALE z0 AND d FROM WIND DIRECTION
    !-------------------------------------------------------------------
    CALL z0_theta(udir(1,1)*pi/180.0, lambda_f_park, h_PV, beta_deg, theta_0, z0_park)
    CALL d_theta(udir(1,1)*pi/180.0, lambda_f_park, h_PV, theta_0, zd_park)
    z0_park_arr(:,:) = z0_park
    
    ! Setup atmospheric state for dragcoef/sfcflx
    tair = tref
    qvair = qvref
    pbar(:,:,:) = psfc(1,1)
    ptbar(:,:,:) = tref(1,1) * (p0/psfc(1,1))**rddcp
    qvbar(:,:,:) = qvref(1,1)
    qv(:,:,:) = qvref(1,1)
    rhostr(:,:,:) = rhosfc(1,1)
    
    u_wind(:,:,:) = uref(1,1) * SIN(udir(1,1) * pi / 180.0)
    v_wind(:,:,:) = uref(1,1) * COS(udir(1,1) * pi / 180.0)
    w_wind(:,:,:) = 0.0
    
    !-------------------------------------------------------------------
    ! 2. WIND PROFILE
    !-------------------------------------------------------------------
    ! Obukhov length for this timestep (from the forcing); neutral if 0 or very large.
    L_obk = L_obk_in(1,1)

    CALL compute_U_plus(z_ref, udir(1,1)*pi/180.0, a, beta_deg, L, h_clear, theta_0, L_obk, Uplus)
    u_star_ref = uref(1,1) / Uplus
    
    CALL compute_U_plus(z_can, udir(1,1)*pi/180.0, a, beta_deg, L, h_clear, theta_0, L_obk, Uplus)
    U_can = u_star_ref * Uplus
    
    CALL compute_U_plus(z_PVa, udir(1,1)*pi/180.0, a, beta_deg, L, h_clear, theta_0, L_obk, Uplus)
    U_PVa = u_star_ref * Uplus

    !-------------------------------------------------------------------
    ! 3. DRAG COEFFICIENTS - computed with current temperatures
    !-------------------------------------------------------------------
    

    zp_dc_nat(:,:,1) = zp(:,:,1)
    zp_dc_nat(:,:,2) = 0.0

    ! --- sol naturel : interface l_nat_can, vent U_can ---
    z0_sfc_arr(:,:)  = z0_nat
    zp_dc_nat(:,:,3) = 2.0*l_nat_can        ! => z1 = l_nat_can
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt_nat(i,j) = T_nat(i,j) * (p0/pbar(i,j,2))**rddcp
        pt_can(i,j) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
      END DO
    END DO
    CALL dragcoef(nx, ny, nz, zp_dc_nat, pt_nat, qvref, soiltyp_dummy, &
                  z0_sfc_arr, max(U_can, U_min_fac), pt_can, cdm_nat, cdh_nat, cdq_nat)

	facteur = 4
	!facteur = 1
	cdm_nat = facteur * cdm_nat
	cdh_nat = facteur * cdh_nat
	cdq_nat = facteur * cdq_nat

    ! --- PV back face: interface l_PVb_can, wind U_can ---
    z0_sfc_arr(:,:)  = z0_PV_back
    zp_dc_nat(:,:,3) = 2.0 * l_PVb_can
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt_PVb(i,j) = T_PVb(i,j) * (p0/pbar(i,j,2))**rddcp
      END DO
    END DO
    CALL dragcoef(nx, ny, nz, zp_dc_nat, pt_PVb, qvref, soiltyp_dummy, &
                  z0_sfc_arr, max(U_can, U_min_fac), pt_can(:,:), cdm_PVb, cdh_PVb, cdq_PVb)

    ! --- cdh_PVb: mixed convection on the BACK face (faces the ground).
    ! dragcoef above still provides cdm_PVb (momentum); only cdh/cdq are replaced here.
    DO j = 1, ny
      DO i = 1, nx
        ! film temperature: property evaluation point for the boundary layer
        Tfilm  = 0.5 * (T_PVb(i,j) + T_can(i,j))
        ! air thermal conductivity k(Tfilm)  [W/m/K]
        kair   = 1.5207e-11*Tfilm**3 - 4.8574e-8*Tfilm**2 + 1.0184e-4*Tfilm - 3.9333e-4
        ! air dynamic viscosity mu(Tfilm), Sutherland's law  [Pa.s]
        muair  = 1.458e-6 * Tfilm**1.5 / (Tfilm + 110.4)
        ! air kinematic viscosity nu = mu/rho  [m2/s]
        nuair  = muair / rhosfc(i,j)
        ! air thermal diffusivity alpha = k/(rho*cp)  [m2/s]
        alair  = kair / (rhosfc(i,j) * cp)
        ! signed surface-to-air temperature difference  [K]
        dTs    = T_PVb(i,j) - T_can(i,j)
        ! buoyancy arrangement: back face points DOWN, so a COLD face (night) is unstable
        unstab = (dTs < 0.0)
        ! Reynolds number on the Sparrow length, using the local canopy wind
        Re_pv  = MAX(U_can(i,j), U_min_fac) * L_star / nuair
        ! forced Nusselt, BACK face: normalised at Re_0 so both faces share C_sparrow
        Nu_fo  = C_sparrow_b * Pr**(1.0/3.0) * SQRT(Re_0) * (Re_pv/Re_0)**p_forced_b
        ! ---- Natural convection: Kaplani & Kaplanis (2014), Sec. 2.1 ----
        ! Rayleigh on the physical length L; Tf after Kaplani.
        Tf_na    = Tfilm - 0.25*dTs
        Ra_na    = 9.81 / Tf_na * MAX(ABS(dTs), 0.05) * L**3 / (nuair*alair)
        theta_na = pi/2.0 - beta         ! angle from the VERTICAL (Kaplani convention)
        IF (unstab) THEN
          ! up-facing hot-plate arrangement (day front / night back): Kaplani Eq.(7)
          IF (beta > pi/6.0) THEN        ! theta < 60 deg
            Grc_na = 1.327e10 * EXP(-3.708*theta_na)
            Rac_na = Grc_na * Pr
            Nu_nat = 0.13*MAX(Ra_na**(1.0/3.0) - Rac_na**(1.0/3.0), 0.0)  &
                   + 0.56*( Rac_na*MAX(SIN(beta),0.0) )**0.25
          ELSE                           ! theta >= 60 deg (near horizontal, turbulent plume)
            Nu_nat = 0.13 * Ra_na**(1.0/3.0)
          END IF
        ELSE
          ! down-facing hot-plate arrangement (day back / night front): Kaplani Eq.(8)-(9)
          IF (beta >= pi/6.0) THEN       ! beta >= 30 deg : Kaplani Eq.(8), all Ra
            Nu_nat = ( 0.825 + 0.387*( Ra_na*MAX(SIN(beta),0.0) )**(1.0/6.0) / cc_pr )**2
          ELSE IF (beta > 2.0*pi/180.0) THEN  ! 2 < beta < 30 deg : Kaplani Eq.(9), inclined
            Nu_nat = 0.56 * ( Ra_na*MAX(SIN(beta),0.0) )**0.25
          ELSE                           ! beta <= 2 deg : Kaplani Eq.(9), horizontal plate
            Nu_nat = 0.58 * Ra_na**0.2
          END IF
        END IF
        Nu_nat = MAX(Nu_nat, 0.0)
        ! Blend on h, not Nu: forced lives on L_star, natural on L (two length scales)
        h_forced = Nu_fo  * kair / L_star
        h_free   = Nu_nat * kair / L
        h_cv     = ( h_forced**3 + h_free**3 )**(1.0/3.0)
        ! convert h to the model's exchange coefficient: H = rho*cp*cdh*U*dT, so U cancels
        cdh_PVb(i,j) = h_cv / ( rhosfc(i,j) * cp * MAX(U_can(i,j), U_min_fac) )
      END DO
    END DO
    ! moisture coefficient follows heat (dry surface: no latent flux anyway)
    cdq_PVb(:,:) = cdh_PVb(:,:)

    ! --- PV front face: interface l_PVf_PVa, wind U_PVa ---
    z0_sfc_arr(:,:)  = z0_PV_front
    zp_dc_nat(:,:,3) = 2.0*l_PVf_PVa
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt_PVf(i,j) = T_PVf(i,j) * (p0/pbar(i,j,2))**rddcp
        pt_PVa(i,j) = T_PVa(i,j) * (p0/pbar(i,j,2))**rddcp
        wdspsfc(i,j)   = U_PVa(i,j)
      END DO
    END DO
    CALL dragcoef(nx, ny, nz, zp_dc_nat, pt_PVf, qvref, soiltyp_dummy, &
                  z0_sfc_arr, max(U_PVa, U_min_fac), pt_PVa(:,:), cdm_PVf, cdh_PVf, cdq_PVf)

    ! --- cdh_PVf : mixed convection on the FRONT face (faces the sky) ---
    DO j = 1, ny
      DO i = 1, nx
        ! film temperature: property evaluation point for the boundary layer
        Tfilm  = 0.5 * (T_PVf(i,j) + T_PVa(i,j))
        ! air thermal conductivity k(Tfilm)  [W/m/K]
        kair   = 1.5207e-11*Tfilm**3 - 4.8574e-8*Tfilm**2 + 1.0184e-4*Tfilm - 3.9333e-4
        ! air dynamic viscosity mu(Tfilm), Sutherland's law  [Pa.s]
        muair  = 1.458e-6 * Tfilm**1.5 / (Tfilm + 110.4)
        ! air kinematic viscosity nu = mu/rho  [m2/s]
        nuair  = muair / rhosfc(i,j)
        ! air thermal diffusivity alpha = k/(rho*cp)  [m2/s]
        alair  = kair / (rhosfc(i,j) * cp)
        ! signed surface-to-air temperature difference  [K]
        dTs    = T_PVf(i,j) - T_PVa(i,j)
        ! buoyancy arrangement: front face points UP, so a HOT face (day) is unstable
        unstab = (dTs > 0.0)
        ! Reynolds number on the Sparrow length, using the local above-panel wind
        Re_pv  = MAX(U_PVa(i,j), U_min_fac) * L_star / nuair
        ! forced Nusselt, FRONT face: free plate, Sparrow et al. (1979).
        ! Re_0-normalised form; with p_forced_f = 0.5 it reduces to C * Pr^(1/3) * Re^0.5.
        Nu_fo  = C_sparrow_f * Pr**(1.0/3.0) * SQRT(Re_0) * (Re_pv/Re_0)**p_forced_f
        ! ---- Natural convection: Kaplani & Kaplanis (2014), Sec. 2.1 ----
        ! Rayleigh on the physical length L; Tf after Kaplani.
        Tf_na    = Tfilm - 0.25*dTs
        Ra_na    = 9.81 / Tf_na * MAX(ABS(dTs), 0.05) * L**3 / (nuair*alair)
        theta_na = pi/2.0 - beta         ! angle from the VERTICAL (Kaplani convention)
        IF (unstab) THEN
          ! up-facing hot-plate arrangement (day front / night back): Kaplani Eq.(7)
          IF (beta > pi/6.0) THEN        ! theta < 60 deg
            Grc_na = 1.327e10 * EXP(-3.708*theta_na)
            Rac_na = Grc_na * Pr
            Nu_nat = 0.13*MAX(Ra_na**(1.0/3.0) - Rac_na**(1.0/3.0), 0.0)  &
                   + 0.56*( Rac_na*MAX(SIN(beta),0.0) )**0.25
          ELSE                           ! theta >= 60 deg (near horizontal, turbulent plume)
            Nu_nat = 0.13 * Ra_na**(1.0/3.0)
          END IF
        ELSE
          ! down-facing hot-plate arrangement (day back / night front): Kaplani Eq.(8)-(9)
          IF (beta >= pi/6.0) THEN       ! beta >= 30 deg : Kaplani Eq.(8), all Ra
            Nu_nat = ( 0.825 + 0.387*( Ra_na*MAX(SIN(beta),0.0) )**(1.0/6.0) / cc_pr )**2
          ELSE IF (beta > 2.0*pi/180.0) THEN  ! 2 < beta < 30 deg : Kaplani Eq.(9), inclined
            Nu_nat = 0.56 * ( Ra_na*MAX(SIN(beta),0.0) )**0.25
          ELSE                           ! beta <= 2 deg : Kaplani Eq.(9), horizontal plate
            Nu_nat = 0.58 * Ra_na**0.2
          END IF
        END IF
        Nu_nat = MAX(Nu_nat, 0.0)
        ! Blend on h, not Nu: forced lives on L_star, natural on L (two length scales)
        h_forced = Nu_fo  * kair / L_star
        h_free   = Nu_nat * kair / L
        h_cv     = ( h_forced**3 + h_free**3 )**(1.0/3.0)
        ! convert h to the model's exchange coefficient: H = rho*cp*cdh*U*dT, so U cancels
        cdh_PVf(i,j) = h_cv / ( rhosfc(i,j) * cp * MAX(U_PVa(i,j), U_min_fac) )
      END DO
    END DO
    ! moisture coefficient follows heat (dry surface: no latent flux anyway)
    cdq_PVf(:,:) = cdh_PVf(:,:)

	    

    ! Canopy air to reference, with stability.
    ! zp_dc(3) is set so that dragcoef's internal z1 = l_can_ref - zd_park.
    
    l_can_ref = z_ref - z_can
    zp_dc(:,:,1) = zp(:,:,1)
    zp_dc(:,:,2) = 0.0
    zp_dc(:,:,3) = 2.0 * (z_ref - max(z_can, zd_park))

    DO j = 1, ny-1
      DO i = 1, nx-1
        pt_can(i,j) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
        wdspsfc(i,j) = U_ref(i,j)
      END DO
    END DO
    
    CALL dragcoef(nx, ny, nz, zp_dc, &
                 pt_can, qvref, soiltyp_dummy, z0_park_arr, &
                 wdspsfc, ptbar(:,:,2), &
                 cdm_can, cdh_can, cdq_can)
                 
    facteur = 10
    !facteur = 1
	cdm_can = facteur * cdm_can
	cdh_can = facteur * cdh_can
	cdq_can = facteur * cdq_can
    
    ! Air above PV to Reference
    l_PVa_ref = z_ref - z_PVa
    zp_dc(:,:,2) = 0.0
    zp_dc(:,:,3) = 2.0 * l_PVa_ref
    
    CALL dragcoef(nx, ny, nz, zp_dc, &
                 pt_PVa, qvref, soiltyp_dummy, z0_park_arr, &
                 wdspsfc, ptbar(:,:,2), &
                 cdm_PVa, cdh_PVa, cdq_PVa)
    
    facteur = 10
    cdm_PVa = facteur * cdm_PVa
	cdh_PVa = facteur * cdh_PVa
	cdq_PVa = facteur * cdq_PVa 
    
    !-------------------------------------------------------------------
    ! 4. SOLVE AIR NODE TEMPERATURES with current T_nat, T_PVb, T_PVf
    !-------------------------------------------------------------------
    CALL solve_air_node_temperatures(nx, ny, &
                                     T_nat, T_PVb, T_PVf, tref, &
                                     U_can, U_PVa, U_ref, &
                                     cdh_nat, cdh_PVb, cdh_can, &
                                     cdh_PVf, cdh_PVa, &
                                     L, a, w, &
                                     T_can, T_PVa)
    
    !-------------------------------------------------------------------
    ! 5. UNIFIED RK2 INTEGRATION FOR ALL NODES
    !-------------------------------------------------------------------
    
    ! STEP 1: Calculate forcing at t^n
    CALL compute_all_forcings(nx, ny, nz, nzsoil, nstyps, &
                             L, a, beta, h_PV, h_clear, w, &
                             z0_nat, albedo_nat, albedo_PVf, albedo_PVb, &
                             emiss_nat, emiss_PVf, emiss_PVb, &
                             soiltyp, vegtyp, lai, veg, &
                             T_nat, T_soil, T_PVf, T_cell, T_PVb, &
                             T_can, T_PVa, &
                             qsoil_surf, qsoil_deep, wetcanp, snowdpth, qvsfc, &
                             uref, tref, qvref, swin, lwin, psfc, rhosfc, prcrate, &
                             U_can, U_PVa, U_ref, &
                             lambda(1,1), phi(1,1), S_diffuse(1,1), &
                             delta_PVf, delta_cell, delta_PVb, &
                             rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                             delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                             k_glass, k_EVA, k_backsheet, &
                             eta_STC, gamma_PV, &
                             zp, u_wind, v_wind, w_wind, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                             cdh_nat, cdh_PVb, cdh_PVf, &
                             cdm_nat, cdm_PVb, cdm_PVf, &
                             cdq_nat, cdq_PVb, cdq_PVf, &
                             frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                             frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                             H_nat, LE_nat, G_nat, Rn_ref_out, &
                             SW_out_sky_loc, LW_out_sky_loc, &
                             dum_qLW_PVf, dum_qLW_PVb, &
                             dum_H_PVf, dum_H_PVb, &
                             dum_qcond_PVf, dum_qcond_PVb, &
                             dum_qSW_PVf, dum_qSW_PVb, &
                             dum_qPV, &
                             dum_cdh_can, dum_cdh_PVa, &
                             dum_T_can, dum_T_PVa, &
                             dum_JLW_PVf, dum_JLW_PVb, &
                             dum_G_SW_nat, dum_J_SW_nat, &
                             dum_G_LW_nat, dum_J_LW_nat)
    
    ! Predictor step: T^(n+1/2) = T^n + (dt/2) * frc^n
    DO j = 1, ny
      DO i = 1, nx
        T_nat_half(i,j)  = T_nat(i,j)  + dtsfc2 * frc_nat(i,j)
        T_soil_half(i,j) = T_soil(i,j) + dtsfc2 * frc_soil(i,j)
        T_PVf_half(i,j)  = T_PVf(i,j)  + dtsfc2 * frc_PVf(i,j)
        T_cell_half(i,j) = T_cell(i,j) + dtsfc2 * frc_cell(i,j)
        T_PVb_half(i,j)  = T_PVb(i,j)  + dtsfc2 * frc_PVb(i,j)
        ! Moisture half-step (clamped to [0, wsat] like ARPS soilebm)
        qsoil_surf_half(i,j) = qsoil_surf(i,j) + dtsfc2 * frc_wsfc(i,j)
        qsoil_surf_half(i,j) = MAX(qsoil_surf_half(i,j), 0.0d0)
        qsoil_surf_half(i,j) = MIN(qsoil_surf_half(i,j), wsat(soiltyp(i,j,1)))
        qsoil_deep_half(i,j) = qsoil_deep(i,j) + dtsfc2 * frc_wdp(i,j)
        qsoil_deep_half(i,j) = MAX(qsoil_deep_half(i,j), 0.0d0)
        qsoil_deep_half(i,j) = MIN(qsoil_deep_half(i,j), wsat(soiltyp(i,j,1)))
      END DO
    END DO
    wetcanp_half = wetcanp
    DO j = 1, ny
      DO i = 1, nx
        wetcanp_half(i,j,1) = wetcanp(i,j,1) + dtsfc2 * frc_wcnp(i,j)
        wetcanp_half(i,j,1) = MAX(wetcanp_half(i,j,1), 0.0d0)
        wetcanp_half(i,j,1) = MIN(wetcanp_half(i,j,1), wrmax(i,j))
      END DO
    END DO
    
    ! Save forcing at t^n
    rhs_nat  = frc_nat
    rhs_soil = frc_soil
    rhs_PVf  = frc_PVf
    rhs_cell = frc_cell
    rhs_PVb  = frc_PVb
    rhs_wsfc = frc_wsfc
    rhs_wdp  = frc_wdp
    rhs_wcnp = frc_wcnp
    
    ! STEP 2: Update air node temps for half-step state, then calc forcing
    CALL solve_air_node_temperatures(nx, ny, &
                                     T_nat_half, T_PVb_half, T_PVf_half, tref, &
                                     U_can, U_PVa, U_ref, &
                                     cdh_nat, cdh_PVb, cdh_can, &
                                     cdh_PVf, cdh_PVa, &
                                     L, a, w, &
                                     T_can, T_PVa)
    
    CALL compute_all_forcings(nx, ny, nz, nzsoil, nstyps, &
                             L, a, beta, h_PV, h_clear, w, &
                             z0_nat, albedo_nat, albedo_PVf, albedo_PVb, &
                             emiss_nat, emiss_PVf, emiss_PVb, &
                             soiltyp, vegtyp, lai, veg, &
                             T_nat_half, T_soil_half, T_PVf_half, T_cell_half, T_PVb_half, &
                             T_can, T_PVa, &
                             qsoil_surf_half, qsoil_deep_half, wetcanp_half, snowdpth, qvsfc, &
                             uref, tref, qvref, swin, lwin, psfc, rhosfc, prcrate, &
                             U_can, U_PVa, U_ref, &
                             lambda(1,1), phi(1,1), S_diffuse(1,1), &
                             delta_PVf, delta_cell, delta_PVb, &
                             rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                             delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                             k_glass, k_EVA, k_backsheet, &
                             eta_STC, gamma_PV, &
                             zp, u_wind, v_wind, w_wind, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                             cdh_nat, cdh_PVb, cdh_PVf, &
                             cdm_nat, cdm_PVb, cdm_PVf, &
                             cdq_nat, cdq_PVb, cdq_PVf, &
                             frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                             frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                             H_nat, LE_nat, G_nat, Rn_ref_out, &
                             SW_out_sky_loc, LW_out_sky_loc, &
                             dum_qLW_PVf, dum_qLW_PVb, &
                             dum_H_PVf, dum_H_PVb, &
                             dum_qcond_PVf, dum_qcond_PVb, &
                             dum_qSW_PVf, dum_qSW_PVb, &
                             dum_qPV, &
                             dum_cdh_can, dum_cdh_PVa, &
                             dum_T_can, dum_T_PVa, &
                             dum_JLW_PVf, dum_JLW_PVb, &
                             dum_G_SW_nat, dum_J_SW_nat, &
                             dum_G_LW_nat, dum_J_LW_nat)
    
    ! Average the two forcings
    DO j = 1, ny
      DO i = 1, nx
        rhs_nat(i,j)  = 0.5 * (rhs_nat(i,j)  + frc_nat(i,j))
        rhs_soil(i,j) = 0.5 * (rhs_soil(i,j) + frc_soil(i,j))
        rhs_PVf(i,j)  = 0.5 * (rhs_PVf(i,j)  + frc_PVf(i,j))
        rhs_cell(i,j) = 0.5 * (rhs_cell(i,j) + frc_cell(i,j))
        rhs_PVb(i,j)  = 0.5 * (rhs_PVb(i,j)  + frc_PVb(i,j))
        rhs_wsfc(i,j) = 0.5 * (rhs_wsfc(i,j) + frc_wsfc(i,j))
        rhs_wdp(i,j)  = 0.5 * (rhs_wdp(i,j)  + frc_wdp(i,j))
        rhs_wcnp(i,j) = 0.5 * (rhs_wcnp(i,j) + frc_wcnp(i,j))
      END DO
    END DO
    
    ! Corrector step: T^(n+1) = T^n + dt * avg(frc)
    DO j = 1, ny
      DO i = 1, nx
        T_nat(i,j)  = T_nat(i,j)  + dtsfc * rhs_nat(i,j)
        T_soil(i,j) = T_soil(i,j) + dtsfc * rhs_soil(i,j)
        T_PVf(i,j)  = T_PVf(i,j)  + dtsfc * rhs_PVf(i,j)
        T_cell(i,j) = T_cell(i,j) + dtsfc * rhs_cell(i,j)
        T_PVb(i,j)  = T_PVb(i,j)  + dtsfc * rhs_PVb(i,j)
        ! Moisture corrector (clamped to [0, wsat] like ARPS soilebm)
        qsoil_surf(i,j) = qsoil_surf(i,j) + dtsfc * rhs_wsfc(i,j)
        qsoil_surf(i,j) = MAX(qsoil_surf(i,j), 0.0d0)
        qsoil_surf(i,j) = MIN(qsoil_surf(i,j), wsat(soiltyp(i,j,1)))
        qsoil_deep(i,j) = qsoil_deep(i,j) + dtsfc * rhs_wdp(i,j)
        qsoil_deep(i,j) = MAX(qsoil_deep(i,j), 0.0d0)
        qsoil_deep(i,j) = MIN(qsoil_deep(i,j), wsat(soiltyp(i,j,1)))
        wetcanp(i,j,1) = wetcanp(i,j,1) + dtsfc * rhs_wcnp(i,j)
        wetcanp(i,j,1) = MAX(wetcanp(i,j,1), 0.0d0)
        wetcanp(i,j,1) = MIN(wetcanp(i,j,1), wrmax(i,j))
      END DO
    END DO
    
  END DO  ! End small timestep loop
  
  !=====================================================================
  ! DIAGNOSTIC: Recompute air node temperatures and fluxes at final state
  !=====================================================================
  
  ! Solve air node temps with final T_nat, T_PVb, T_PVf
  CALL solve_air_node_temperatures(nx, ny, &
                                   T_nat, T_PVb, T_PVf, tref, &
                                   U_can, U_PVa, U_ref, &
                                   cdh_nat, cdh_PVb, cdh_can, &
                                   cdh_PVf, cdh_PVa, &
                                   L, a, w, &
                                   T_can, T_PVa)
  
  ! Recompute all forcings at final state (we only use the flux outputs)
  CALL compute_all_forcings(nx, ny, nz, nzsoil, nstyps, &
                           L, a, beta, h_PV, h_clear, w, &
                           z0_nat, albedo_nat, albedo_PVf, albedo_PVb, &
                           emiss_nat, emiss_PVf, emiss_PVb, &
                           soiltyp, vegtyp, lai, veg, &
                           T_nat, T_soil, T_PVf, T_cell, T_PVb, &
                           T_can, T_PVa, &
                           qsoil_surf, qsoil_deep, wetcanp, snowdpth, qvsfc, &
                           uref, tref, qvref, swin, lwin, psfc, rhosfc, prcrate, &
                           U_can, U_PVa, U_ref, &
                           lambda(1,1), phi(1,1), S_diffuse(1,1), &
                           delta_PVf, delta_cell, delta_PVb, &
                           rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                           delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                           k_glass, k_EVA, k_backsheet, &
                           eta_STC, gamma_PV, &
                           zp, u_wind, v_wind, w_wind, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                           cdh_nat, cdh_PVb, cdh_PVf, &
                           cdm_nat, cdm_PVb, cdm_PVf, &
                           cdq_nat, cdq_PVb, cdq_PVf, &
                           frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                           frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                           H_nat, LE_nat, G_nat, Rn_ref_out, &
                           SW_out_sky, LW_out_sky, &
                           diag_qLW_PVf, diag_qLW_PVb, &
                           diag_H_PVf, diag_H_PVb, &
                           diag_qcond_PVf, diag_qcond_PVb, &
                           diag_qSW_PVf, diag_qSW_PVb, &
                           diag_qPV, &
                           diag_cdh_can_out, diag_cdh_PVa_out, &
                           diag_T_can_out, diag_T_PVa_out, &
                           diag_JLW_PVf_out, diag_JLW_PVb_out, &
                           G_SW_nat_out, J_SW_nat_out, &
                           G_LW_nat_out, J_LW_nat_out)
  ! frc_* are discarded; H_nat, LE_nat, G_nat, Rn_ref_out, SW/LW_out_sky matter
  
  ! Override cdh diagnostics with the TRUE can->ref and PVa->ref values
  diag_cdh_can_out = cdh_can
  diag_cdh_PVa_out = cdh_PVa
  diag_U_can_out = U_can
  diag_U_ref_out = U_ref
  
  
  !=====================================================================
  ! COMPUTE OUTPUT FLUXES (from final-state diagnostics)
  !=====================================================================
  
  ! H_ref = weighted combination of H_PVa_ref and H_can_ref
  DO j = 1, ny-1
    DO i = 1, nx-1
      H_PVa_ref_out(i,j) = rhosfc(i,j) * cp * cdh_PVa(i,j) * U_ref(i,j) * &
                  (T_PVa(i,j) - tref(i,j))
      H_can_ref_out(i,j) = rhosfc(i,j) * cp * cdh_can(i,j) * U_ref(i,j) * &
                  (T_can(i,j) - tref(i,j))
      
      H_ref_out(i,j) = H_PVa_ref_out(i,j) * w / (w + a) + H_can_ref_out(i,j) * a / (w + a)
    END DO
  END DO
  

  LE_ref_out = LE_nat
  G_out = G_nat
  H_nat_out = H_nat

 



END SUBROUTINE pv_land
