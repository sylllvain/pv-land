!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE COMPUTE_ALL_FORCINGS            ######
!######                                                      ######
!######   Calculate forcing terms (dT/dt) for ALL nodes      ######
!######   Called TWICE per timestep for RK2 integration      ######
!######   SIMPLIFIED VERSION using soilebm_frc and sfcflx    ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_all_forcings(nx, ny, nz, nzsoil, nstyps, &
                               L, a, beta, h_PV, h_clear, w, &
                               z0_nat, z0_PV, &
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
                               H_nat_out, LE_nat_out, G_nat_out, Rn_ref_out)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate forcing terms (dT/dt in K/s) for ALL nodes (nat, soil, PV).

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
  REAL, INTENT(IN) :: z0_nat, z0_PV
  
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
  
  ! Forcing
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
  
  ! Drag coefficients (from pv_land - calculated once per timestep)
  REAL, INTENT(IN) :: cdh_nat(nx,ny)
  REAL, INTENT(IN) :: cdh_PVb(nx,ny)
  REAL, INTENT(IN) :: cdh_PVf(nx,ny)
  REAL, INTENT(IN) :: cdm_nat(nx,ny)
  REAL, INTENT(IN) :: cdm_PVb(nx,ny)
  REAL, INTENT(IN) :: cdm_PVf(nx,ny)
  REAL, INTENT(IN) :: cdq_nat(nx,ny)
  REAL, INTENT(IN) :: cdq_PVb(nx,ny)
  REAL, INTENT(IN) :: cdq_PVf(nx,ny)
  
  ! OUTPUT: forcing terms (K/s)
  REAL, INTENT(OUT) :: frc_nat(nx,ny), frc_soil(nx,ny)
  REAL, INTENT(OUT) :: frc_PVf(nx,ny), frc_cell(nx,ny), frc_PVb(nx,ny)
  
  ! Diagnostic outputs
  REAL, INTENT(OUT) :: H_nat_out(nx,ny), LE_nat_out(nx,ny), G_nat_out(nx,ny)
  REAL, INTENT(OUT) :: Rn_ref_out(nx,ny)
  
  ! Local - Radiation
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
  
  ! Soil variables (for soilebm_frc)
  REAL :: frc_tsfc(nx,ny), frc_tdeep(nx,ny), frc_wsfc(nx,ny)
  REAL :: frc_wdp(nx,ny), frc_wcnp(nx,ny)
  REAL :: wrmax(nx,ny), qvsata(nx,ny), qvsat_local(nx,ny), f34(nx,ny)
  REAL :: ct(nx,ny), shflx(nx,ny), lhflx(nx,ny), gflx(nx,ny), rnflx(nx,ny)
  REAL :: evaprg(nx,ny), evaprtr(nx,ny), evaprr(nx,ny)
  REAL :: relief
  REAL :: radsw(nx,ny)
  
  ! Conductive fluxes (W/m2)
  REAL :: q_cond_PVf, q_cond_PVb
  REAL :: R_PVf, R_PVb
  
  ! PV power (W/m2)
  REAL :: q_dot_PV, eta_PV, T_cell_C
  
  ! Heat capacities (J/m2/K)
  REAL :: C_PVf, C_cell, C_PVb
  
  INTEGER :: i, j, is
  REAL, PARAMETER :: pi = 3.14159265359
  
  ! Local potential temperature copies (to avoid modifying INTENT(IN) arrays)
  REAL :: pt_nat(nx,ny), pt_can(nx,ny)
  REAL :: pt_PVb(nx,ny), pt_PVf(nx,ny), pt_PVa(nx,ny)
  
  !=====================================================================
  ! 1. RADIATION - Temperature-dependent!
  !=====================================================================
  
  ! DEBUG: VÃ©rifier les entrÃ©es
  !PRINT *, '=== DEBUG compute_all_forcings ==='
  !PRINT *, 'T_nat(1,1) = ', T_nat(1,1)
  !PRINT *, 'T_PVf(1,1) = ', T_PVf(1,1)
  !PRINT *, 'T_PVb(1,1) = ', T_PVb(1,1)
  !PRINT *, 'T_can(1,1) = ', T_can(1,1)
  !PRINT *, 'T_PVa(1,1) = ', T_PVa(1,1)
  !PRINT *, 'cdh_nat(1,1) = ', cdh_nat(1,1)
  !PRINT *, 'cdm_nat(1,1) = ', cdm_nat(1,1)
  !PRINT *, 'cdq_nat(1,1) = ', cdq_nat(1,1)
  !PRINT *, 'U_can(1,1) = ', U_can(1,1)
  !PRINT *, 'qvref(1,1) = ', qvref(1,1)
  
  ! View factors (geometry only, computed once)
  CALL compute_view_factor_matrix(L, a, w, h_PV, beta, F_matrix)
  CALL compute_sky_view_factors(F_matrix, F_sky)
  
  !PRINT *, ''
  ! ===== FIN DEBUG =====

  ! Surface optical properties
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
  
  ! Direct SW radiation
  CALL compute_S_ground(swin(1,1), L, a, beta, h_PV, lambda, phi, S_direct(1))
  CALL compute_S_PVf(swin(1,1), L, a, beta, h_PV, lambda, phi, S_direct(2))
  CALL compute_S_PVb(swin(1,1), L, a, beta, h_PV, lambda, phi, S_direct(3))
  S_direct(4) = S_direct(3)
  S_direct(5) = S_direct(2)
  S_direct(6) = S_direct(3)
  S_direct(7) = S_direct(1)
  S_direct(8) = S_direct(1)
  
  ! Use input diffuse radiation
  S_diffuse = S_diffuse_in
  L_down = lwin(1,1)
  
  ! SW radiosity
  CALL solve_SW_radiosities(8, alpha_surf, F_matrix, F_sky, S_direct, &
                            S_diffuse, J_SW, info_rad)
  CALL compute_SW_irradiances(8, F_matrix, F_sky, S_direct, S_diffuse, &
                             J_SW, G_SW)
  CALL compute_SW_net_fluxes(8, alpha_surf, G_SW, q_SW)
  
  ! LW radiosity (TEMPERATURE-DEPENDENT!)
  CALL solve_LW_radiosities(8, epsilon_surf, T_surf, F_matrix, F_sky, &
                           L_down, J_LW, info_rad)
  CALL compute_LW_irradiances(8, F_matrix, F_sky, J_LW, L_down, G_LW)
  CALL compute_LW_net_fluxes(8, G_LW, J_LW, q_LW)
  
  ! Extract fluxes
  q_dot_SW_nat = q_SW(1)
  q_dot_LW_nat = q_LW(1)
  q_dot_SW_PVf = q_SW(2)
  q_dot_LW_PVf = q_LW(2)
  q_dot_SW_PVb = q_SW(3)
  q_dot_LW_PVb = q_LW(3)
  G_SW_PVf = G_SW(2)
  G_SW_nat = G_SW(1)
  
  ! Cell absorbs SW at front
  q_dot_SW_cell = q_dot_SW_PVf
  
  ! Calculate area-weighted net radiation for the whole park
  Rn_ref_out = (L * (q_dot_SW_PVf + q_dot_LW_PVf) + &
                L * (q_dot_SW_PVb + q_dot_LW_PVb) + &
                (a + w) * (q_dot_SW_nat + q_dot_LW_nat)) / (a + w)
  
  !=====================================================================
  ! 2. CONVECTION - Use sfcflx directly
  !=====================================================================
  
  ! (1) Natural surface to Canopy air - already have cdh_nat from pv_land
  ! convert to potential temperatures (LOCAL copies, don't modify inputs)
  DO j = 1, ny-1
    DO i = 1, nx-1
      pt_nat(i,j) = T_nat(i,j) * (p0/pbar(i,j,2))**rddcp
      pt_can(i,j) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt, pprt, qv, pt_can, pbar, qvbar, &
             U_can, rhosfc, pt_nat, qvref, &
             cdm_nat, cdh_nat, cdq_nat, &
             usflx_dummy, vsflx_dummy, ptsflx_nat, qvsflx_nat)
  
  ! (2) PV back to Canopy air
  ! convert to potential temperatures (LOCAL copies)
  DO j = 1, ny-1
    DO i = 1, nx-1
      pt_PVb(i,j) = T_PVb(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt, pprt, qv, pt_can, pbar, qvbar, &
             U_can, rhosfc, pt_PVb, qvref, &
             cdm_PVb, cdh_PVb, cdq_PVb, &
             usflx_dummy, vsflx_dummy, ptsflx_PVb, qvsflx_PVb)
  
  ! (3) PV front to Air above PV
  ! convert to potential temperatures (LOCAL copies)
  DO j = 1, ny-1
    DO i = 1, nx-1
      pt_PVf(i,j) = T_PVf(i,j) * (p0/pbar(i,j,2))**rddcp
      pt_PVa(i,j) = T_PVa(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt, pprt, qv, pt_PVa, pbar, qvbar, &
             U_PVa, rhosfc, pt_PVf, qvref, &
             cdm_PVf, cdh_PVf, cdq_PVf, &
             usflx_dummy, vsflx_dummy, ptsflx_PVf, qvsflx_PVf)
  
  !=====================================================================
  ! 3. SOIL FORCINGS - Use soilebm_frc directly
  !=====================================================================
  
  ! Calculate f34 and qvsata (needed by soilebm_frc)
  CALL getqvs(nx, ny, 1, 1, nx-1, 1, ny-1, 1, 1, psfc, tref, qvsata)
  
  DO j = 1, ny-1
    DO i = 1, nx-1
      IF (tref(i,j) >= 302.15) THEN
        f34(i,j) = 1.0 / (1.0 + EXP(1.18 * (tref(i,j) - 314.00)))
      ELSE
        f34(i,j) = 1.0 / (1.0 + EXP(-0.41 * (tref(i,j) - 282.05)))
      END IF
    END DO
  END DO
  
  ! Setup for soilebm_frc
  radsw = G_SW_nat  ! Total SW irradiance on natural surface (from radiosity)
  ! Use net radiative fluxes from radiosity calculations
  rnflx = q_dot_SW_nat + q_dot_LW_nat  ! Net radiation for natural surface
  relief = 0.0  ! Simplified
  
  ! Call soilebm_frc for each soil type
  ! Note: soilebm_frc expects T_can as real temperature (not potential)
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
  
  ! Extract diagnostics
  !H_nat_out = shflx
  !LE_nat_out = lhflx
  H_nat_out = -cp*ptsflx_nat(1,1)
  LE_nat_out = -lathv*qvsflx_nat(1,1)
  G_nat_out = gflx

  !=====================================================================
  ! 4. PV NODE FORCINGS
  !=====================================================================
  
  DO j = 1, ny-1
    DO i = 1, nx-1
      
      ! Convective fluxes (from sfcflx results)
      ! H > 0 means surface is warmer than air → surface LOSES heat
      H_nat = -cp*ptsflx_nat(i,j)
      H_PVb = -cp*ptsflx_PVb(i,j)
      H_PVf = -cp*ptsflx_PVf(i,j)
      
      !===================================================================
      ! Conduction
      !===================================================================
      R_PVf = delta_glass / k_glass + delta_EVA_top / k_EVA
      R_PVb = delta_EVA_bot / k_EVA + delta_backsheet / k_backsheet
      
      q_cond_PVf = (T_cell(i,j) - T_PVf(i,j)) / R_PVf
      q_cond_PVb = (T_cell(i,j) - T_PVb(i,j)) / R_PVb
      
      !===================================================================
      ! PV electrical power
      !===================================================================
      T_cell_C = T_cell(i,j) - 273.15
      eta_PV = eta_STC * (1.0 - gamma_PV * (T_cell_C - 25.0))
      q_dot_PV = eta_PV * G_SW_PVf
      
      !===================================================================
      ! Heat capacities
      !===================================================================
      C_PVf = rho_cp_PVf * delta_PVf
      C_cell = rho_cp_cell * delta_cell
      C_PVb = rho_cp_PVb * delta_PVb
      
      !===================================================================
      ! Forcing terms: dT/dt = Q_total / C  (K/s)
      !===================================================================
      
      ! Front surface: radiation + conduction - convection loss
      frc_PVf(i,j) = (q_dot_LW_PVf + q_cond_PVf - H_PVf) / C_PVf
      
      ! Cell: SW absorbed - conduction out - PV power
      frc_cell(i,j) = (q_dot_SW_cell - q_cond_PVf - q_cond_PVb - q_dot_PV) / C_cell
      
      ! Back surface: radiation + conduction - convection loss
      frc_PVb(i,j) = (q_dot_LW_PVb + q_cond_PVb - H_PVb) / C_PVb
      
    END DO
  END DO

END SUBROUTINE compute_all_forcings
