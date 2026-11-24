!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE PV_NODES_FORCING                ######
!######                                                      ######
!######   Calculate forcing terms (dT/dt) for PV nodes       ######
!######   Called TWICE per timestep for RK2 integration      ######
!######   Like soilebm_frc in soilebm3d.f90                  ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE pv_nodes_forcing(nx, ny, nz, &
                           L, a, beta, h_PV, h_clear, w, &
                           z0_PV, albedo, &
                           T_nat, T_PVf, T_cell, T_PVb, T_can, T_PVa, &
                           uref, swin, lwin, psfc, rhosfc, &
                           U_can, U_PVa, &
                           lambda, phi, &
                           delta_PVf, delta_cell, delta_PVb, &
                           rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                           delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                           k_glass, k_EVA, k_backsheet, &
                           eta_STC, gamma_PV, &
                           frc_PVf, frc_cell, frc_PVb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate forcing terms (dT/dt in K/s) for the 3 PV thermal nodes.
!  
!  This subroutine computes:
!    - Radiation (SW + LW) - depends on current temperatures
!    - Convection - depends on current temperatures  
!    - Conduction within PV module
!    - PV electrical power conversion
!
!  Analogous to soilebm_frc - will be called twice per timestep:
!    1st call: with T^n to get frc^n
!    2nd call: with T^(n+1/2) to get frc^(n+1/2)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  INCLUDE 'phycst.inc'
  
  ! Grid
  INTEGER, INTENT(IN) :: nx, ny, nz
  
  ! Geometry
  REAL, INTENT(IN) :: L, a, beta, h_PV, h_clear, w
  REAL, INTENT(IN) :: z0_PV, albedo
  
  ! Current temperatures (input state)
  REAL, INTENT(IN) :: T_nat(nx,ny), T_PVf(nx,ny), T_cell(nx,ny), T_PVb(nx,ny)
  REAL, INTENT(IN) :: T_can(nx,ny), T_PVa(nx,ny)
  
  ! Forcing
  REAL, INTENT(IN) :: uref(nx,ny), swin(nx,ny), lwin(nx,ny)
  REAL, INTENT(IN) :: psfc(nx,ny), rhosfc(nx,ny)
  REAL, INTENT(IN) :: U_can(nx,ny), U_PVa(nx,ny)
  
  ! Solar geometry
  REAL, INTENT(IN) :: lambda, phi
  
  ! PV thermal parameters
  REAL, INTENT(IN) :: delta_PVf, delta_cell, delta_PVb
  REAL, INTENT(IN) :: rho_cp_PVf, rho_cp_cell, rho_cp_PVb
  REAL, INTENT(IN) :: delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet
  REAL, INTENT(IN) :: k_glass, k_EVA, k_backsheet
  REAL, INTENT(IN) :: eta_STC, gamma_PV
  
  ! OUTPUT: forcing terms (K/s)
  REAL, INTENT(OUT) :: frc_PVf(nx,ny), frc_cell(nx,ny), frc_PVb(nx,ny)
  
  ! Local - Radiation
  REAL :: F_matrix(8,8), F_sky(8)
  REAL :: S_direct(8), S_diffuse, L_down
  REAL :: alpha_surf(8), epsilon_surf(8), T_surf(8)
  REAL :: J_SW(8), G_SW(8), q_SW(8)
  REAL :: J_LW(8), G_LW(8), q_LW(8)
  INTEGER :: info_rad
  
  ! Radiative fluxes per unit area (W/mÂ²)
  REAL :: q_dot_SW_PVf, q_dot_LW_PVf
  REAL :: q_dot_SW_cell
  REAL :: q_dot_SW_PVb, q_dot_LW_PVb
  REAL :: G_SW_PVf
  
  ! Drag coefficients
  REAL :: cdh_PVb, cdh_PVf
  
  ! Convective fluxes (W/mÂ²)
  REAL :: H_PVf, H_PVb
  
  ! Conductive fluxes (W/mÂ²)
  REAL :: q_cond_PVf, q_cond_PVb
  REAL :: R_PVf, R_PVb
  
  ! PV power (W/mÂ²)
  REAL :: q_dot_PV, eta_PV, T_cell_C
  
  ! Heat capacities (J/mÂ²/K)
  REAL :: C_PVf, C_cell, C_PVb
  
  INTEGER :: i, j
  REAL, PARAMETER :: pi = 3.14159265359
  REAL, PARAMETER :: kappa = 0.387
  REAL, PARAMETER :: Pr = 0.71
  REAL :: l_interface, Cu, Ctheta
  
  !=====================================================================
  ! 1. RADIATION - Temperature-dependent!
  !=====================================================================
  
  ! View factors (geometry only, computed once)
  CALL compute_view_factor_matrix(L, a, w, h_PV, beta, F_matrix)
  CALL compute_sky_view_factors(F_matrix, F_sky)
  
  ! Surface optical properties
  alpha_surf(1) = albedo       ! nat,n
  alpha_surf(2) = 0.10         ! PVf,n
  alpha_surf(3) = 0.15         ! PVb,n
  alpha_surf(4) = alpha_surf(3)
  alpha_surf(6) = alpha_surf(3)
  alpha_surf(5) = alpha_surf(2)
  alpha_surf(7) = alpha_surf(1)
  alpha_surf(8) = alpha_surf(1)

  epsilon_surf = 0.95
  
  ! Current temperatures for radiosity
  T_surf(1) = T_nat(1,1)
  T_surf(2) = T_PVf(1,1)
  T_surf(3) = T_PVb(1,1)
  T_surf(4) = T_surf(3)
  T_surf(6) = T_surf(3)
  T_surf(5) = T_surf(2)
  T_surf(7) = T_surf(1)
  T_surf(8) = T_surf(1)
  
  ! Direct SW radiation
  CALL compute_S_ground(swin(1,1), L, a, beta, h_PV, lambda, phi, S_direct(1))
  CALL compute_S_PVf(swin(1,1), L, a, beta, h_PV, lambda, phi, S_direct(2))
  CALL compute_S_PVb(swin(1,1), L, a, beta, h_PV, lambda, phi, S_direct(3))
  S_direct(4) = S_direct(3)
  S_direct(6) = S_direct(3)
  S_direct(5) = S_direct(2)
  S_direct(7) = S_direct(1)
  S_direct(8) = S_direct(1)
  
  S_diffuse = swin(1,1) * 0.2
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
  
  ! Extract fluxes for PV nodes
  q_dot_SW_PVf = q_SW(2)
  q_dot_LW_PVf = q_LW(2)
  q_dot_SW_PVb = q_SW(3)
  q_dot_LW_PVb = q_LW(3)
  G_SW_PVf = G_SW(2)
  
  ! Cell absorbs SW at front (internal layer, no LW)
  q_dot_SW_cell = q_dot_SW_PVf
  
  !=====================================================================
  ! 2. DRAG COEFFICIENTS (neutral, simple)
  !=====================================================================
  
  ! PVb to canopy air (Eq. 6.7-6.9)
  l_interface = 0.50
  Cu = kappa / LOG(l_interface / z0_PV)
  Ctheta = kappa / (Pr * LOG(l_interface / z0_PV))
  cdh_PVb = Cu * Ctheta
  
  ! PVf to air above PV (Eq. 6.11-6.13)
  l_interface = 0.20
  Cu = kappa / LOG(l_interface / z0_PV)
  Ctheta = kappa / (Pr * LOG(l_interface / z0_PV))
  cdh_PVf = Cu * Ctheta
  
  !=====================================================================
  ! 3. CONVECTION - Temperature-dependent!
  !=====================================================================
  
  DO j = 1, ny
    DO i = 1, nx
      
      ! Eq. 6.10: H_PVb = rho * cp * cdh * U * (T_PVb - T_can)
      H_PVb = rhosfc(i,j) * cp * cdh_PVb * U_can(i,j) * &
              (T_PVb(i,j) - T_can(i,j))
      
      ! Eq. 6.14: H_PVf = rho * cp * cdh * U * (T_PVf - T_PVa)
      H_PVf = rhosfc(i,j) * cp * cdh_PVf * U_PVa(i,j) * &
              (T_PVf(i,j) - T_PVa(i,j))
      
      !===================================================================
      ! 4. CONDUCTION - Temperature-dependent!
      !===================================================================
      
      ! Thermal resistances (Eq. 9.5)
      R_PVf = delta_glass / k_glass + delta_EVA_top / k_EVA
      R_PVb = delta_EVA_bot / k_EVA + delta_backsheet / k_backsheet
      
      ! Conductive fluxes (Eq. 9.4)
      q_cond_PVf = (T_cell(i,j) - T_PVf(i,j)) / R_PVf
      q_cond_PVb = (T_cell(i,j) - T_PVb(i,j)) / R_PVb
      
      !===================================================================
      ! 5. PV ELECTRICAL POWER - Temperature-dependent!
      !===================================================================
      
      ! Eq. 10.2: eta = eta_STC * [1 - gamma*(T_cell - 25Â°C)]
      T_cell_C = T_cell(i,j) - 273.15
      eta_PV = eta_STC * (1.0 - gamma_PV * (T_cell_C - 25.0))
      
      ! Eq. 10.1: q_PV = eta * G_SW_PVf
      q_dot_PV = eta_PV * G_SW_PVf
      
      !===================================================================
      ! 6. HEAT CAPACITIES (constant)
      !===================================================================
      
      C_PVf = rho_cp_PVf * delta_PVf
      C_cell = rho_cp_cell * delta_cell
      C_PVb = rho_cp_PVb * delta_PVb
      
      !===================================================================
      ! 7. FORCING TERMS: dT/dt = Q_total / C  (K/s)
      !===================================================================
      
      ! Eq. 9.1: Front surface
      ! C_PVf * dT_PVf/dt = q_cond + q_conv + q_LW
      frc_PVf(i,j) = (q_cond_PVf + H_PVf + q_dot_LW_PVf) / C_PVf
      
      ! Eq. 9.2: Cell (SW absorbed, conduction out, PV power out)
      ! C_cell * dT_cell/dt = q_SW - q_cond_PVf - q_cond_PVb - q_PV
      frc_cell(i,j) = (q_dot_SW_cell - q_cond_PVf - q_cond_PVb - q_dot_PV) / C_cell
      
      ! Eq. 9.3: Back surface
      ! C_PVb * dT_PVb/dt = q_cond + q_conv + q_LW
      frc_PVb(i,j) = (q_cond_PVb + H_PVb + q_dot_LW_PVb) / C_PVb
      
    END DO
  END DO

END SUBROUTINE pv_nodes_forcing
