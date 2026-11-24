!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE COMPUTE_ALL_FORCINGS            ######
!######                                                      ######
!######   Calculate forcing terms (dT/dt) for ALL nodes      ######
!######   Called TWICE per timestep for RK2 integration       ######
!######   Uses soilebm_frc and sfcflx                        ######
!######                                                      ######
!######   COMBINED version: outputs both PV diagnostics      ######
!######   (diag_qLW_*, diag_H_*, etc.) AND nat surface       ######
!######   irradiances/radiosities (G_SW_nat, J_SW_nat, etc.) ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_all_forcings(nx, ny, nz, nzsoil, nstyps, &
                               L, a, beta, h_PV, h_clear, w, &
                               z0_nat, &
                               albedo_nat, albedo_PVf, albedo_PVb, &
                               emiss_nat, emiss_PVf, emiss_PVb, &
                               soiltyp, vegtyp, lai, veg, &
                               T_nat, T_soil, T_PVf, T_cell, T_PVb, &
                               T_can, T_PVa, &
                               wetsfc, wetdp, wetcanp, snowdpth, qvsfc, &
                               uref, tref, qvref, swin, lwin, psfc, rhosfc, precip, &
                               U_can, U_PVa, U_ref, &
                               lambda, phi, S_diffuse_in, &
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
                               H_nat_out, LE_nat_out, G_nat_out, Rn_ref_out, &
                               SW_out_sky_out, LW_out_sky_out, &
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
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate forcing terms (dT/dt in K/s) for all thermal nodes
!  (natural surface, soil, PV front/cell/back).
!
!  COMBINED version outputting:
!    - PV energy balance diagnostics (from _pv version)
!    - Natural surface irradiances/radiosities (from _nat version)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
  
  ! Grid
  INTEGER, INTENT(IN) :: nx, ny, nz, nzsoil, nstyps
  
  ! Geometry
  REAL, INTENT(IN) :: L, a, beta, h_PV, h_clear, w
  REAL, INTENT(IN) :: z0_nat
  
  ! Optical properties
  REAL, INTENT(IN) :: albedo_nat, albedo_PVf, albedo_PVb
  REAL, INTENT(IN) :: emiss_nat, emiss_PVf, emiss_PVb
  
  ! Soil/vegetation properties
  INTEGER, INTENT(IN) :: soiltyp(nx,ny,nstyps), vegtyp(nx,ny)
  REAL, INTENT(IN) :: lai(nx,ny), veg(nx,ny)
  
  ! Current temperatures (input state)
  REAL, INTENT(IN) :: T_nat(nx,ny), T_soil(nx,ny)
  REAL, INTENT(IN) :: T_PVf(nx,ny), T_cell(nx,ny), T_PVb(nx,ny)
  REAL, INTENT(IN) :: T_can(nx,ny), T_PVa(nx,ny)
  
  ! Soil moisture
  REAL, INTENT(IN) :: wetsfc(nx,ny), wetdp(nx,ny)
  REAL, INTENT(IN) :: wetcanp(nx,ny,0:nstyps), qvsfc(nx,ny,0:nstyps)
  REAL, INTENT(IN) :: snowdpth(nx,ny)
  
  ! Atmospheric forcing
  REAL, INTENT(IN) :: uref(nx,ny), tref(nx,ny), qvref(nx,ny)
  REAL, INTENT(IN) :: swin(nx,ny), lwin(nx,ny)
  REAL, INTENT(IN) :: psfc(nx,ny), rhosfc(nx,ny), precip(nx,ny)
  REAL, INTENT(IN) :: U_can(nx,ny), U_PVa(nx,ny), U_ref(nx,ny)
  
  ! Solar geometry
  REAL, INTENT(IN) :: lambda, phi, S_diffuse_in
  
  ! PV thermal parameters
  REAL, INTENT(IN) :: delta_PVf, delta_cell, delta_PVb
  REAL, INTENT(IN) :: rho_cp_PVf, rho_cp_cell, rho_cp_PVb
  REAL, INTENT(IN) :: delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet
  REAL, INTENT(IN) :: k_glass, k_EVA, k_backsheet
  REAL, INTENT(IN) :: eta_STC, gamma_PV
  
  ! Atmospheric state (for sfcflx)
  REAL, INTENT(IN) :: zp(nx,ny,nz)
  REAL, INTENT(IN) :: u_wind(nx,ny,nz), v_wind(nx,ny,nz), w_wind(nx,ny,nz)
  REAL, INTENT(IN) :: ptprt(nx,ny,nz), pprt(nx,ny,nz), qv(nx,ny,nz)
  REAL, INTENT(IN) :: ptbar(nx,ny,nz), pbar(nx,ny,nz), qvbar(nx,ny,nz)
  
  ! Drag coefficients (computed in pv_land, once per small timestep)
  REAL, INTENT(IN) :: cdh_nat(nx,ny)
  REAL, INTENT(IN) :: cdh_PVb(nx,ny)
  REAL, INTENT(IN) :: cdh_PVf(nx,ny)
  REAL, INTENT(IN) :: cdm_nat(nx,ny)
  REAL, INTENT(IN) :: cdm_PVb(nx,ny)
  REAL, INTENT(IN) :: cdm_PVf(nx,ny)
  REAL, INTENT(IN) :: cdq_nat(nx,ny)
  REAL, INTENT(IN) :: cdq_PVb(nx,ny)
  REAL, INTENT(IN) :: cdq_PVf(nx,ny)
  
  ! OUTPUT: temperature forcing terms (K/s)
  REAL, INTENT(OUT) :: frc_nat(nx,ny), frc_soil(nx,ny)
  REAL, INTENT(OUT) :: frc_PVf(nx,ny), frc_cell(nx,ny), frc_PVb(nx,ny)
  
  ! OUTPUT: soil moisture forcing terms
  REAL, INTENT(OUT) :: frc_wsfc(nx,ny), frc_wdp(nx,ny), frc_wcnp(nx,ny)
  REAL, INTENT(OUT) :: wrmax(nx,ny)
  
  ! Diagnostic outputs
  REAL, INTENT(OUT) :: H_nat_out(nx,ny), LE_nat_out(nx,ny), G_nat_out(nx,ny)
  REAL, INTENT(OUT) :: Rn_ref_out(nx,ny)
  REAL, INTENT(OUT) :: SW_out_sky_out(nx,ny), LW_out_sky_out(nx,ny)
  
  ! PV energy balance diagnostics (W/m2) -- from _pv version
  REAL, INTENT(OUT) :: diag_qLW_PVf(nx,ny), diag_qLW_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_H_PVf(nx,ny), diag_H_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_qcond_PVf(nx,ny), diag_qcond_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_qSW_PVf(nx,ny), diag_qSW_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_qPV(nx,ny)
  REAL, INTENT(OUT) :: diag_cdh_can_out(nx,ny), diag_cdh_PVa_out(nx,ny)
  REAL, INTENT(OUT) :: diag_T_can_out(nx,ny), diag_T_PVa_out(nx,ny)
  REAL, INTENT(OUT) :: diag_JLW_PVf_out(nx,ny), diag_JLW_PVb_out(nx,ny)
  
  ! Natural surface irradiances/radiosities -- from _nat version
  REAL, INTENT(OUT) :: G_SW_nat_out(nx,ny), J_SW_nat_out(nx,ny)
  REAL, INTENT(OUT) :: G_LW_nat_out(nx,ny), J_LW_nat_out(nx,ny)
  
  ! Local - Radiation (8-surface radiosity system)
  REAL :: F_matrix(8,8), F_sky(8)
  REAL :: S_direct(8), S_diffuse, L_down
  REAL :: alpha_surf(8), epsilon_surf(8), T_surf(8)
  REAL :: J_SW(8), G_SW(8), q_SW(8)
  REAL :: J_LW(8), G_LW(8), q_LW(8)
  INTEGER :: info_rad
  
  ! Radiative fluxes per unit area (W/m2)
  REAL :: q_dot_SW_PVf, q_dot_LW_PVf
  REAL :: q_dot_SW_cell
  REAL :: q_dot_SW_PVb, q_dot_LW_PVb
  REAL :: q_dot_SW_nat, q_dot_LW_nat
  REAL :: G_SW_PVf, G_SW_nat
  
  ! Convective fluxes (W/m2) - from sfcflx
  REAL :: ptsflx_nat(nx,ny), qvsflx_nat(nx,ny)
  REAL :: ptsflx_PVb(nx,ny), qvsflx_PVb(nx,ny)
  REAL :: ptsflx_PVf(nx,ny), qvsflx_PVf(nx,ny)
  REAL :: usflx_dummy(nx,ny), vsflx_dummy(nx,ny)
  REAL :: H_PVf, H_PVb, H_nat
  
  ! Local 3D arrays for sfcflx
  REAL :: ptbar_loc(nx,ny,nz), ptprt_loc(nx,ny,nz)
  REAL :: ptsfc_loc(nx,ny)
  
  ! Soil variables (for soilebm_frc)
  REAL :: frc_tsfc(nx,ny), frc_tdeep(nx,ny)
  REAL :: qvsata(nx,ny), qvsat_local(nx,ny), f34(nx,ny)
  REAL :: ct(nx,ny), shflx(nx,ny), lhflx(nx,ny), gflx(nx,ny), rnflx(nx,ny)
  REAL :: evaprg(nx,ny), evaprtr(nx,ny), evaprr(nx,ny)
  REAL :: relief
  REAL :: radsw(nx,ny)
  REAL :: f3
  REAL :: VPD_kPa
  
  ! Conductive fluxes (W/m2)
  REAL :: q_cond_PVf, q_cond_PVb
  REAL :: R_PVf, R_PVb
  
  ! PV electrical power (W/m2)
  REAL :: q_dot_PV, eta_PV, T_cell_C
  
  ! Heat capacities (J/m2/K)
  REAL :: C_PVf, C_cell, C_PVb
  
  INTEGER :: i, j, is
  REAL, PARAMETER :: pi = 3.14159265359
  
  ! Debug call counter
  INTEGER, SAVE :: call_count = 0
  
  !=====================================================================
  ! 1. RADIATION (temperature-dependent)
  !=====================================================================
  
  ! View factors (geometry only)
  CALL compute_view_factor_matrix(L, a, w, h_PV, beta, F_matrix)
  CALL compute_sky_view_factors(F_matrix, F_sky)
  
  ! Surface optical properties (8-surface periodic array layout)
  ! Surfaces: 1=nat_n, 2=PVf_n, 3=PVb_n, 4=PVb_{n-1}, 5=PVf_{n-1},
  !           6=PVb_{n+1}, 7=nat_{n-1}, 8=nat_{n+1}
  alpha_surf(1) = albedo_nat       ! nat,n
  alpha_surf(2) = albedo_PVf       ! PVf,n
  alpha_surf(3) = albedo_PVb       ! PVb,n
  alpha_surf(4) = albedo_PVb       ! PVb,n-1  
  alpha_surf(5) = albedo_PVf       ! PVf,n-1
  alpha_surf(6) = albedo_PVb       ! PVb,n+1
  alpha_surf(7) = albedo_nat       ! nat,n-1
  alpha_surf(8) = albedo_nat       ! nat,n+1  

  epsilon_surf(1) = emiss_nat      ! nat,n
  epsilon_surf(2) = emiss_PVf      ! PVf,n
  epsilon_surf(3) = emiss_PVb      ! PVb,n
  epsilon_surf(4) = emiss_PVb      ! PVb,n-1  
  epsilon_surf(5) = emiss_PVf      ! PVf,n-1
  epsilon_surf(6) = emiss_PVb      ! PVb,n+1
  epsilon_surf(7) = emiss_nat      ! nat,n-1
  epsilon_surf(8) = emiss_nat      ! nat,n+1  

  T_surf(1) = T_nat(1,1)           ! nat,n
  T_surf(2) = T_PVf(1,1)           ! PVf,n
  T_surf(3) = T_PVb(1,1)           ! PVb,n
  T_surf(4) = T_PVb(1,1)           ! PVb,n-1 
  T_surf(5) = T_PVf(1,1)           ! PVf,n-1
  T_surf(6) = T_PVb(1,1)           ! PVb,n+1
  T_surf(7) = T_nat(1,1)           ! nat,n-1
  T_surf(8) = T_nat(1,1)           ! nat,n+1 
  
  ! Direct SW on each surface (diffuse is handled separately in radiosity)
  CALL compute_S_ground(swin(1,1) - S_diffuse_in, L, a, beta, h_PV, lambda, phi, S_direct(1))
  CALL compute_S_PVf(swin(1,1) - S_diffuse_in, L, a, beta, h_PV, lambda, phi, S_direct(2))
  CALL compute_S_PVb(swin(1,1) - S_diffuse_in, L, a, beta, h_PV, lambda, phi, S_direct(3))
  S_direct(4) = S_direct(3)
  S_direct(5) = S_direct(2)
  S_direct(6) = S_direct(3)
  S_direct(7) = S_direct(1)
  S_direct(8) = S_direct(1)
  
  S_diffuse = S_diffuse_in
  L_down = lwin(1,1)
  
  ! SW radiosity
  CALL solve_SW_radiosities(8, alpha_surf, F_matrix, F_sky, S_direct, &
                            S_diffuse, J_SW, info_rad)
  CALL compute_SW_irradiances(8, F_matrix, F_sky, S_direct, S_diffuse, &
                             J_SW, G_SW)
  CALL compute_SW_net_fluxes(8, alpha_surf, G_SW, q_SW)
  
  ! LW radiosity (temperature-dependent)
  CALL solve_LW_radiosities(8, epsilon_surf, T_surf, F_matrix, F_sky, &
                           L_down, J_LW, info_rad)
  CALL compute_LW_irradiances(8, F_matrix, F_sky, J_LW, L_down, G_LW)
  CALL compute_LW_net_fluxes(8, G_LW, J_LW, q_LW)
  
  ! Store LW radiosities for PV diagnostic output (from _pv)
  diag_JLW_PVf_out(1,1) = J_LW(2)
  diag_JLW_PVb_out(1,1) = J_LW(3)
  
  ! Extract per-surface net fluxes
  q_dot_SW_nat = q_SW(1)
  q_dot_LW_nat = q_LW(1)
  q_dot_SW_PVf = q_SW(2)
  q_dot_LW_PVf = q_LW(2)
  q_dot_SW_PVb = q_SW(3)
  q_dot_LW_PVb = q_LW(3)
  G_SW_PVf = G_SW(2)
  G_SW_nat = G_SW(1)
  
  ! Store nat surface irradiances and radiosities for output (from _nat)
  G_SW_nat_out = G_SW(1)   ! SW irradiance received by nat (SW in nat)
  J_SW_nat_out = J_SW(1)   ! SW radiosity leaving nat (SW out nat)
  G_LW_nat_out = G_LW(1)   ! LW irradiance received by nat (LW in nat)
  J_LW_nat_out = J_LW(1)   ! LW radiosity leaving nat (LW out nat)
  
  ! Cell absorbs SW transmitted through front glass
  q_dot_SW_cell = q_dot_SW_PVf
  
  ! Area-weighted net radiation for the park
  Rn_ref_out = (L * (q_dot_SW_PVf + q_dot_LW_PVf) + &
                L * (q_dot_SW_PVb + q_dot_LW_PVb) + &
                (a + w) * (q_dot_SW_nat + q_dot_LW_nat)) / (a + w)
  
  ! Outgoing fluxes to sky (as seen by a horizontal radiometer above the park)
  SW_out_sky_out = (L * F_sky(2) * J_SW(2) + &
                    L * F_sky(3) * J_SW(3) + &
                    (a + w) * F_sky(1) * J_SW(1)) / (a + w)
  
  LW_out_sky_out = (L * F_sky(2) * J_LW(2) + &
                    L * F_sky(3) * J_LW(3) + &
                    (a + w) * F_sky(1) * J_LW(1)) / (a + w)
  
  !=====================================================================
  ! 2. CONVECTION - Using sfcflx with proper 3D arrays
  !=====================================================================
  
  ptprt_loc = 0.0
  
  ! (1) Natural surface -> Canopy air
  DO j = 1, ny-1
    DO i = 1, nx-1
      ptbar_loc(i,j,:) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
      ptsfc_loc(i,j)   = T_nat(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt_loc, pprt, qv, ptbar_loc, pbar, qvbar, &
             U_can, rhosfc, ptsfc_loc, qvref, &
             cdm_nat, cdh_nat, cdq_nat, &
             usflx_dummy, vsflx_dummy, ptsflx_nat, qvsflx_nat)
  
  ! (2) PV back -> Canopy air
  DO j = 1, ny-1
    DO i = 1, nx-1
      ptsfc_loc(i,j) = T_PVb(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt_loc, pprt, qv, ptbar_loc, pbar, qvbar, &
             U_can, rhosfc, ptsfc_loc, qvref, &
             cdm_PVb, cdh_PVb, cdq_PVb, &
             usflx_dummy, vsflx_dummy, ptsflx_PVb, qvsflx_PVb)
  
  ! (3) PV front -> Air above PV
  DO j = 1, ny-1
    DO i = 1, nx-1
      ptbar_loc(i,j,:) = T_PVa(i,j) * (p0/pbar(i,j,2))**rddcp
      ptsfc_loc(i,j)   = T_PVf(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt_loc, pprt, qv, ptbar_loc, pbar, qvbar, &
             U_PVa, rhosfc, ptsfc_loc, qvref, &
             cdm_PVf, cdh_PVf, cdq_PVf, &
             usflx_dummy, vsflx_dummy, ptsflx_PVf, qvsflx_PVf)
  
  !=====================================================================
  ! 3. SOIL FORCINGS (via soilebm_frc)
  !=====================================================================
  
  ! Saturated specific humidity at reference level
  CALL getqvs(nx, ny, 1, 1, nx-1, 1, ny-1, 1, 1, psfc, T_can, qvsata)
  
  DO j = 1, ny-1
    DO i = 1, nx-1
    
      ! f3: stomatal response to VPD
      IF (qvsata(i,j) - qvref(i,j) <= 12.5) THEN
        f3 = 1.0 - 0.06*(qvsata(i,j) - qvref(i,j))
      ELSE
        f3 = 0.25
      END IF
      
      ! f4: stomatal response to temperature
      IF (T_can(i,j) >= 302.15) THEN
        f34(i,j) = 1.0 / (1.0 + EXP(1.18 * (T_can(i,j) - 314.00)))
      ELSE
        f34(i,j) = 1.0 / (1.0 + EXP(-0.41 * (T_can(i,j) - 282.05)))
      END IF

      ! f3: stomatal response to VPD
      VPD_kPa = (qvsata(i,j) - qvref(i,j)) * psfc(i,j) / 622.0

      ! Hyperbolic formulation (Kim & Verma, 1991)
      !f3 = 1.0 / (1.0 + 0.13 * VPD_kPa)

      f34(i,j) = f34(i,j) * f3
      !f34(i,j) = f34(i,j)
    END DO
  END DO
  
  ! Setup radiation inputs for soilebm_frc
  radsw = G_SW_nat  ! Total SW irradiance on natural surface (from radiosity)
  rnflx = q_dot_SW_nat + q_dot_LW_nat  ! Net radiation for natural surface
  relief = 0.0
  
  ! Call soilebm_frc for each soil type
  DO is = nstyps, 1, -1
    CALL soilebm_frc(nx, ny, nz, soiltyp(1,1,is), vegtyp, lai, veg, &
                    T_nat, T_soil, wetsfc, wetdp, wetcanp(1,1,is), snowdpth, &
                    qvsfc(1,1,is), U_can, psfc, rhosfc, precip, &
                    T_can, qvref, cdh_nat, cdq_nat, radsw, rnflx, &
                    shflx, lhflx, gflx, ct, &
                    evaprg, evaprtr, evaprr, qvsat_local, qvsata, f34, &
                    frc_tsfc, frc_tdeep, frc_wsfc, frc_wdp, frc_wcnp, &
                    wrmax, relief)
  END DO
  
  ! Extract soil forcings
  frc_nat = frc_tsfc
  frc_soil = frc_tdeep
  
  ! Extract output fluxes
  H_nat_out(1,1)  = -cp * ptsflx_nat(1,1)
  LE_nat_out(1,1) = lhflx(1,1)
  G_nat_out = gflx
  
  ! --- H_nat COMPARISON: sfcflx vs soilebm_frc ---
  IF (call_count <= 15) THEN
    WRITE(*,'(A)') '--- H_nat COMPARISON ---'
    WRITE(*,'(A,F12.4)') '  sfcflx   H_nat  = ', -cp * ptsflx_nat(1,1)
    WRITE(*,'(A,F12.4)') '  soilebm  shflx  = ', shflx(1,1)
    WRITE(*,'(A,F12.4)') '  diff             = ', -cp * ptsflx_nat(1,1) - shflx(1,1)
    WRITE(*,'(A)') '------------------------'
  END IF
  
  ! --- LE COMPARISON: sfcflx vs soilebm_frc ---
  call_count = call_count + 1
  IF (call_count <= 5) THEN
    WRITE(*,'(A)') '--- LE COMPARISON (nat surface) ---'
    WRITE(*,'(A,F12.4)') '  sfcflx  LE_nat  = ', -lathv * qvsflx_nat(1,1)
    WRITE(*,'(A,F12.4)') '  soilebm lhflx   = ', lhflx(1,1)
    WRITE(*,'(A,F12.4)') '  sfcflx  qvsflx  = ', qvsflx_nat(1,1)
    WRITE(*,'(A,F12.4)') '  soilebm evaprg  = ', evaprg(1,1)
    WRITE(*,'(A,F12.4)') '  soilebm evaprtr = ', evaprtr(1,1)
    WRITE(*,'(A,F12.4)') '  soilebm evaprr  = ', evaprr(1,1)
    WRITE(*,'(A)') '-----------------------------------'
  END IF

  !=====================================================================
  ! 4. PV NODE FORCINGS
  !=====================================================================
  
  DO j = 1, ny-1
    DO i = 1, nx-1
      
      ! Convective fluxes from sfcflx (positive = sfc loses heat)
      H_nat = -cp * ptsflx_nat(i,j)
      H_PVb = -cp * ptsflx_PVb(i,j)
      H_PVf = -cp * ptsflx_PVf(i,j)
      
      ! Conduction through PV layers
      R_PVf = delta_glass / k_glass + delta_EVA_top / k_EVA
      R_PVb = delta_EVA_bot / k_EVA + delta_backsheet / k_backsheet
      
      q_cond_PVf = (T_cell(i,j) - T_PVf(i,j)) / R_PVf
      q_cond_PVb = (T_cell(i,j) - T_PVb(i,j)) / R_PVb
      
      ! PV electrical power output
      T_cell_C = T_cell(i,j) - 273.15
      eta_PV = eta_STC * (1.0 - gamma_PV * (T_cell_C - 25.0))
      q_dot_PV = eta_PV * G_SW_PVf
      
      ! Heat capacities
      C_PVf = rho_cp_PVf * delta_PVf
      C_cell = rho_cp_cell * delta_cell
      C_PVb = rho_cp_PVb * delta_PVb
      
      ! Forcing terms: dT/dt = Q_net / C  (K/s)
      
      ! Front glass: LW radiation + conduction from cell - convective loss
      frc_PVf(i,j) = (q_dot_LW_PVf + q_cond_PVf - H_PVf) / C_PVf
      
      ! Cell: SW absorbed - conduction to front and back - electrical power
      frc_cell(i,j) = (q_dot_SW_cell - q_cond_PVf - q_cond_PVb - q_dot_PV) / C_cell
      
      IF (i==1 .AND. j==1 .AND. call_count <= 4) THEN
        WRITE(*,'(A)') '--- CELL BALANCE ---'
        WRITE(*,'(A,F12.3)') '  q_SW_cell    = ', q_dot_SW_cell
        WRITE(*,'(A,F12.3)') '  q_cond_PVf   = ', q_cond_PVf
        WRITE(*,'(A,F12.3)') '  q_cond_PVb   = ', q_cond_PVb
        WRITE(*,'(A,F12.3)') '  q_dot_PV     = ', q_dot_PV
        WRITE(*,'(A,F12.3)') '  C_cell       = ', C_cell
        WRITE(*,'(A,ES12.4)') '  frc_cell     = ', frc_cell(1,1)
        WRITE(*,'(A,F12.3)') '  T_cell       = ', T_cell(1,1)
        WRITE(*,'(A,F12.3)') '  T_PVf        = ', T_PVf(1,1)
        WRITE(*,'(A,F12.3)') '  T_PVb        = ', T_PVb(1,1)
        WRITE(*,'(A)') '--------------------'
      END IF
      
      ! Back surface: LW radiation + conduction from cell - convective loss
      frc_PVb(i,j) = (q_dot_LW_PVb + q_dot_SW_PVb + q_cond_PVb - H_PVb) / C_PVb
      
      ! Store PV energy balance diagnostics (from _pv)
      diag_qLW_PVf(i,j)  = q_dot_LW_PVf
      diag_qLW_PVb(i,j)  = q_dot_LW_PVb
      diag_H_PVf(i,j)    = H_PVf
      diag_H_PVb(i,j)    = H_PVb
      diag_qcond_PVf(i,j) = q_cond_PVf
      diag_qcond_PVb(i,j) = q_cond_PVb
      diag_qSW_PVf(i,j)  = q_dot_SW_PVf
      diag_qSW_PVb(i,j)  = q_dot_SW_PVb
      diag_qPV(i,j)      = q_dot_PV
      diag_cdh_can_out(i,j)  = cdh_nat(i,j)
      diag_cdh_PVa_out(i,j)  = cdh_PVf(i,j)
      diag_T_can_out(i,j) = T_can(i,j)
      diag_T_PVa_out(i,j) = T_PVa(i,j)
      
    END DO
  END DO

END SUBROUTINE compute_all_forcings
