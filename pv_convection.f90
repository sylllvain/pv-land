!##################################################################
!##################################################################
!######                                                      ######
!######      SUBROUTINE SOLVE_AIR_NODE_TEMPERATURES          ######
!######                                                      ######
!######   Solve for T_can and T_PVa using Eq. 2.3            ######
!######   C_can = C_PVa = 0 => q_conv = 0                    ######
!######   Direct algebraic solution with coefficients at t   ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE solve_air_node_temperatures(nx, ny, &
                                       T_nat, T_PVb, T_PVf, T_ref, &
                                       U_can, U_PVa, U_ref, &
                                       cdh_nat, cdh_PVb, cdh_can, &
                                       cdh_PVf, cdh_PVa, &
                                       L, a, w, &
                                       T_can, T_PVa)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Solve for air node temperatures T_can and T_PVa
!  using energy balance with C=0 (Eq. 2.3).
!
!  From the condition q_conv = 0 (Table 2):
!
!  For T_can:
!    H_nat + (L/(w+a))*H_PVb - (a/(w+a))*H_can–ref = 0
!
!  For T_PVa:
!    (L/w)*H_PVf - H_PVa–ref = 0
!
!  Using H = ρ*cp*Cdh*U*(T_surface - T_air), we get:
!
!  T_can = [Cdh_nat*U_can*T_nat + (L/(w+a))*Cdh_PVb*U_can*T_PVb 
!           + (a/(w+a))*Cdh_can*U_ref*T_ref]
!        / [Cdh_nat*U_can + (L/(w+a))*Cdh_PVb*U_can 
!           + (a/(w+a))*Cdh_can*U_ref]
!
!  T_PVa = [(L/w)*Cdh_PVf*U_PVa*T_PVf + Cdh_PVa*U_ref*T_ref]
!        / [(L/w)*Cdh_PVf*U_PVa + Cdh_PVa*U_ref]
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  :: nx, ny
  REAL, INTENT(IN)     :: T_nat(nx,ny), T_PVb(nx,ny)
  REAL, INTENT(IN)     :: T_PVf(nx,ny), T_ref(nx,ny)
  REAL, INTENT(IN)     :: U_can(nx,ny), U_PVa(nx,ny), U_ref(nx,ny)
  REAL, INTENT(IN)     :: cdh_nat(nx,ny), cdh_PVb(nx,ny)
  REAL, INTENT(IN)     :: cdh_can(nx,ny), cdh_PVf(nx,ny)
  REAL, INTENT(IN)     :: cdh_PVa(nx,ny)
  REAL, INTENT(IN)     :: L, a, w
  REAL, INTENT(OUT)    :: T_can(nx,ny), T_PVa(nx,ny)
  
  REAL :: numerator_can, denominator_can
  REAL :: numerator_PVa, denominator_PVa
  REAL :: ratio_L_wa, ratio_a_wa, ratio_L_w
  INTEGER :: i, j
  
  ! Pre-compute geometric ratios
  ratio_L_wa = L / (w + a)
  ratio_a_wa = a / (w + a)
  ratio_L_w = L / w
  
  DO j = 1, ny-1
    DO i = 1, nx-1
    
      !-----------------------------------------------------------------
      ! Solve for T_can
      !-----------------------------------------------------------------
      numerator_can = cdh_nat(i,j) * U_can(i,j) * T_nat(i,j) + &
                     ratio_L_wa * cdh_PVb(i,j) * U_can(i,j) * T_PVb(i,j) + &
                     ratio_a_wa * cdh_can(i,j) * U_ref(i,j) * T_ref(i,j)
      
      denominator_can = cdh_nat(i,j) * U_can(i,j) + &
                       ratio_L_wa * cdh_PVb(i,j) * U_can(i,j) + &
                       ratio_a_wa * cdh_can(i,j) * U_ref(i,j)
      
      T_can(i,j) = numerator_can / denominator_can
      
      !-----------------------------------------------------------------
      ! Solve for T_PVa
      !-----------------------------------------------------------------
      numerator_PVa = ratio_L_w * cdh_PVf(i,j) * U_PVa(i,j) * T_PVf(i,j) + &
                     cdh_PVa(i,j) * U_ref(i,j) * T_ref(i,j)
      
      denominator_PVa = ratio_L_w * cdh_PVf(i,j) * U_PVa(i,j) + &
                       cdh_PVa(i,j) * U_ref(i,j)
      
      T_PVa(i,j) = numerator_PVa / denominator_PVa
      
    END DO
  END DO

END SUBROUTINE solve_air_node_temperatures
