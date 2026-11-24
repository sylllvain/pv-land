!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE PV_THERMAL_PARAMS            ######
!######                                                      ######
!######   Set thermal parameters for PV module               ######
!######   Based on Table 3 of PV-LAND PDF (Section 9)        ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE pv_thermal_params(delta_PVf, delta_cell, delta_PVb, &
                             rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                             k_glass, k_EVA, k_backsheet, &
                             delta_glass, delta_EVA_top, &
                             delta_EVA_bot, delta_backsheet, &
                             eta_STC, gamma_PV)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Returns thermal parameters for the PV module 3-node model.
!  Values from Table 3 of the PV-LAND documentation (Khyani et al. 2023).
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!    delta_PVf         Thickness of front surface layer (m)
!    delta_cell        Thickness of cell layer (m)
!    delta_PVb         Thickness of back surface layer (m)
!    rho_cp_PVf        Volumetric heat capacity - front (J/mÃ‚Â³/K)
!    rho_cp_cell       Volumetric heat capacity - cell (J/mÃ‚Â³/K)
!    rho_cp_PVb        Volumetric heat capacity - back (J/mÃ‚Â³/K)
!    k_glass           Thermal conductivity of glass (W/m/K)
!    k_EVA             Thermal conductivity of EVA (W/m/K)
!    k_backsheet       Thermal conductivity of backsheet (W/m/K)
!    delta_glass       Thickness of glass sub-layer (m)
!    delta_EVA_top     Thickness of top EVA sub-layer (m)
!    delta_EVA_bot     Thickness of bottom EVA sub-layer (m)
!    delta_backsheet   Thickness of backsheet sub-layer (m)
!    eta_STC           PV efficiency at STC (-)
!    gamma_PV          Temperature coefficient (1/K)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  REAL, INTENT(OUT) :: delta_PVf, delta_cell, delta_PVb
  REAL, INTENT(OUT) :: rho_cp_PVf, rho_cp_cell, rho_cp_PVb
  REAL, INTENT(OUT) :: k_glass, k_EVA, k_backsheet
  REAL, INTENT(OUT) :: delta_glass, delta_EVA_top
  REAL, INTENT(OUT) :: delta_EVA_bot, delta_backsheet
  REAL, INTENT(OUT) :: eta_STC, gamma_PV
  
  ! Equivalent layer thicknesses (m) - Table 3
  delta_PVf  = 3.2e-3
  delta_cell = 0.5e-3
  delta_PVb  = 0.5e-3
  
  ! Volumetric heat capacities (J/mÃ‚Â³/K) - Table 3
  rho_cp_PVf  = 1.6e6
  rho_cp_cell = 1.6e6
  rho_cp_PVb  = 1.8e6
  
  ! Thermal conductivities (W/m/K)
  k_glass     = 1.8
  k_EVA       = 0.33
  k_backsheet = 0.20
  
  ! Sub-layer thicknesses (m) - for resistance calculation
  delta_glass     = 2.2e-3
  delta_EVA_top   = 1.0e-3
  delta_EVA_bot   = 0.5e-3
  delta_backsheet = 0.5e-3
  
  ! PV electrical parameters
  eta_STC  = 0.18      ! Efficiency at 25°C, 1000 W/m²
  gamma_PV = 0.0045    ! Temperature coefficient (1/K)

END SUBROUTINE pv_thermal_params


!##################################################################
!##################################################################
!######                                                      ######
!######         SUBROUTINE COMPUTE_THERMAL_RESISTANCES       ######
!######                                                      ######
!######   Calculate thermal resistances R_PVf and R_PVb      ######
!######   Equation (9.5) from PV-LAND PDF                  ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_thermal_resistances(L, delta_glass, delta_EVA_top, &
                                      delta_EVA_bot, delta_backsheet, &
                                      k_glass, k_EVA, k_backsheet, &
                                      R_PVf, R_PVb)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute thermal resistances per unit depth for the PV module.
!  Equation (9.5) from PV-LAND documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    L                 Panel width (m)
!    delta_glass       Glass layer thickness (m)
!    delta_EVA_top     Upper EVA layer thickness (m)
!    delta_EVA_bot     Lower EVA layer thickness (m)
!    delta_backsheet   Backsheet layer thickness (m)
!    k_glass           Thermal conductivity of glass (W/m/K)
!    k_EVA             Thermal conductivity of EVA (W/m/K)
!    k_backsheet       Thermal conductivity of backsheet (W/m/K)
!
!  OUTPUT:
!
!    R_PVf             Thermal resistance - front surface (KÃ‚Â·m/W)
!    R_PVb             Thermal resistance - back surface (KÃ‚Â·m/W)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  REAL, INTENT(IN)  :: L
  REAL, INTENT(IN)  :: delta_glass, delta_EVA_top
  REAL, INTENT(IN)  :: delta_EVA_bot, delta_backsheet
  REAL, INTENT(IN)  :: k_glass, k_EVA, k_backsheet
  REAL, INTENT(OUT) :: R_PVf, R_PVb
  
  ! Equation (9.5) - front resistance (glass + upper EVA)
  R_PVf = delta_glass / (k_glass) + &
          delta_EVA_top / (k_EVA)
  
  ! Equation (9.5) - back resistance (lower EVA + backsheet)
  R_PVb = delta_EVA_bot / (k_EVA) + &
          delta_backsheet / (k_backsheet)

END SUBROUTINE compute_thermal_resistances


!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE COMPUTE_PV_POWER               ######
!######                                                      ######
!######   Calculate PV electrical power production           ######
!######   Equations (10.1) and (10.2) from PV-LAND PDF       ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_pv_power(G_SW_PVf, T_cell, eta_STC, gamma_PV, &
                           q_dot_PV)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate electrical power extracted from the PV cell.
!  Equations (10.1) and (10.2) from PV-LAND documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    G_SW_PVf          Shortwave irradiance on front surface (W/mÃ‚Â²)
!    T_cell            Cell temperature (K)
!    eta_STC           PV efficiency at STC (-)
!    gamma_PV          Temperature coefficient (1/K)
!
!  OUTPUT:
!
!    q_dot_PV          Electrical power density (W/mÃ‚Â²)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  REAL, INTENT(IN)  :: G_SW_PVf, T_cell
  REAL, INTENT(IN)  :: eta_STC, gamma_PV
  REAL, INTENT(OUT) :: q_dot_PV
  
  REAL :: T_cell_C, eta_PV
  
  ! Convert to Celsius
  T_cell_C = T_cell - 273.15
  
  ! Equation (10.2) - temperature-dependent efficiency
  eta_PV = eta_STC * (1.0 - gamma_PV * (T_cell_C - 25.0))
  
  ! Equation (10.1) - electrical power
  q_dot_PV = eta_PV * G_SW_PVf

END SUBROUTINE compute_pv_power


!##################################################################
!##################################################################
!######                                                      ######
!######       SUBROUTINE UPDATE_PV_TEMPERATURES              ######
!######                                                      ######
!######   Time integration of PV module temperatures         ######
!######   Equations (9.1), (9.2), (9.3) from PV-LAND         ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE update_pv_temperatures(T_PVf, T_cell, T_PVb, &
                                  q_dot_SW_PVf, q_dot_LW_PVf, &
                                  q_dot_conv_PVf, &
                                  q_dot_LW_PVb, q_dot_conv_PVb, &
                                  G_SW_PVf, L, dt, &
                                  delta_PVf, delta_cell, delta_PVb, &
                                  rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                                  delta_glass, delta_EVA_top, &
                                  delta_EVA_bot, delta_backsheet, &
                                  k_glass, k_EVA, k_backsheet, &
                                  eta_STC, gamma_PV)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Update PV module temperatures using explicit Euler integration.
!  Equations (9.1), (9.2), (9.3) from PV-LAND documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT/OUTPUT:
!
!    T_PVf             Temperature of front surface (K) [INOUT]
!    T_cell            Temperature of cell layer (K) [INOUT]
!    T_PVb             Temperature of back surface (K) [INOUT]
!
!  INPUT:
!
!    q_dot_SW_PVf      Net SW flux at front surface (W/m)
!    q_dot_LW_PVf      Net LW flux at front surface (W/m)
!    q_dot_conv_PVf    Convective flux at front surface (W/m)
!    q_dot_LW_PVb      Net LW flux at back surface (W/m)
!    q_dot_conv_PVb    Convective flux at back surface (W/m)
!    G_SW_PVf          Irradiance on front surface (W/mÃ‚Â²)
!    L                 Panel width (m)
!    dt                Time step (s)
!    delta_*           Layer thicknesses (m)
!    rho_cp_*          Volumetric heat capacities (J/mÃ‚Â³/K)
!    k_*               Thermal conductivities (W/m/K)
!    eta_STC           PV efficiency at STC (-)
!    gamma_PV          Temperature coefficient (1/K)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  REAL, INTENT(INOUT) :: T_PVf, T_cell, T_PVb
  REAL, INTENT(IN)    :: q_dot_SW_PVf, q_dot_LW_PVf, q_dot_conv_PVf
  REAL, INTENT(IN)    :: q_dot_LW_PVb, q_dot_conv_PVb
  REAL, INTENT(IN)    :: G_SW_PVf, L, dt
  REAL, INTENT(IN)    :: delta_PVf, delta_cell, delta_PVb
  REAL, INTENT(IN)    :: rho_cp_PVf, rho_cp_cell, rho_cp_PVb
  REAL, INTENT(IN)    :: delta_glass, delta_EVA_top
  REAL, INTENT(IN)    :: delta_EVA_bot, delta_backsheet
  REAL, INTENT(IN)    :: k_glass, k_EVA, k_backsheet
  REAL, INTENT(IN)    :: eta_STC, gamma_PV
  
  REAL :: C_PVf, C_cell, C_PVb
  REAL :: R_PVf, R_PVb
  REAL :: q_dot_cond_PVf, q_dot_cond_PVb
  REAL :: q_dot_PV
  REAL :: dE_PVf, dE_cell, dE_PVb
  
  ! Heat capacities per unit depth (J/K/m)
  C_PVf  = rho_cp_PVf  * delta_PVf
  C_cell = rho_cp_cell * delta_cell
  C_PVb  = rho_cp_PVb  * delta_PVb
  
  ! Thermal resistances
  CALL compute_thermal_resistances(L, delta_glass, delta_EVA_top, &
                                  delta_EVA_bot, delta_backsheet, &
                                  k_glass, k_EVA, k_backsheet, &
                                  R_PVf, R_PVb)
  
  ! Conductive fluxes (Equation 9.4)
  q_dot_cond_PVf = (T_cell - T_PVf) / R_PVf
  q_dot_cond_PVb = (T_cell - T_PVb) / R_PVb
  
  ! PV electrical power
  CALL compute_pv_power(G_SW_PVf, T_cell, eta_STC, gamma_PV, &
                       q_dot_PV)
  
  ! Equation (9.1) - Front surface energy balance
  dE_PVf = q_dot_cond_PVf + q_dot_conv_PVf + q_dot_LW_PVf
  T_PVf = T_PVf + (dE_PVf / C_PVf) * dt
  
  ! Equation (9.2) - Cell energy balance (SW absorbed here)
  dE_cell = q_dot_SW_PVf - q_dot_cond_PVf - q_dot_cond_PVb - q_dot_PV
  T_cell = T_cell + (dE_cell / C_cell) * dt
  
  ! Equation (9.3) - Back surface energy balance
  dE_PVb = q_dot_cond_PVb + q_dot_conv_PVb + q_dot_LW_PVb
  T_PVb = T_PVb + (dE_PVb / C_PVb) * dt

END SUBROUTINE update_pv_temperatures
