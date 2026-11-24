SUBROUTINE solve_LW_radiosities(N, epsilon, T, F, F_sky, L_down, &
                                 J_LW, info)
!-----------------------------------------------------------------------
! Solve the linear system for longwave radiosities (Eq. 5.4 from PDF):
!   [I - (I - diag(ÃŽÂµ))Ãƒâ€”F] Ãƒâ€” J_LW = [ÃŽÂµ_i ÃÆ’ T_i^4] + [I - diag(ÃŽÂµ)]Ãƒâ€”F_sÃƒâ€”LÃ¢â€ â€œ
!
! Uses LAPACK routine DGESV (Double precision General matrix SolVe)
!
! Input:
!   N         : Number of surfaces
!   epsilon(N): Emissivity of each surface (0 to 1)
!   T(N)      : Temperature of each surface (K)
!   F(N,N)    : View factor matrix F(i,j) = view factor from i to j
!   F_sky(N)  : Sky view factor of each surface
!   L_down    : Incoming longwave radiation from atmosphere (W/mÃ‚Â²)
!
! Output:
!   J_LW(N)   : Longwave radiosity of each surface (W/mÃ‚Â²)
!   info      : 0 if successful, /=0 if error
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: epsilon(N), T(N), F(N,N), F_sky(N)
  REAL, INTENT(IN) :: L_down
  REAL, INTENT(OUT) :: J_LW(N)
  INTEGER, INTENT(OUT) :: info
  
  ! Local variables
  REAL :: A(N,N), b(N)
  INTEGER :: i, j
  INTEGER :: ipiv(N)
  REAL, PARAMETER :: sigma = 5.67e-8  ! Stefan-Boltzmann constant (W/mÃ‚Â²/KÃ¢ÂÂ´)
  
  ! Initialize
  info = 0
  A = 0.0
  b = 0.0
  
  ! Build matrix A = I - (I - diag(ÃŽÂµ))Ãƒâ€”F
  ! Start with identity
  DO i = 1, N
    A(i,i) = 1.0
  END DO
  
  ! Subtract (I - diag(ÃŽÂµ))Ãƒâ€”F = F - diag(ÃŽÂµ)Ãƒâ€”F
  DO i = 1, N
    DO j = 1, N
      A(i,j) = A(i,j) - (1.0 - epsilon(i)) * F(i,j)
    END DO
  END DO
  
  ! Build right-hand side b = [ÃŽÂµ_i ÃÆ’ T_i^4] + [I - diag(ÃŽÂµ)]Ãƒâ€”F_sÃƒâ€”LÃ¢â€ â€œ
  DO i = 1, N
    ! Emission term
    b(i) = epsilon(i) * sigma * T(i)**4
    
    ! Reflected atmospheric radiation term
    b(i) = b(i) + (1.0 - epsilon(i)) * F_sky(i) * L_down
  END DO
  
  ! Solve linear system A Ãƒâ€” J_LW = b using LAPACK DGESV
  !CALL DGESV(N, 1, A, N, ipiv, b, N, info)
  CALL SGESV(N, 1, A, N, ipiv, b, N, info) !modif sylvain
  
  IF (info /= 0) THEN
    PRINT *, 'ERROR in solve_LW_radiosities: DGESV failed with info=', info
    IF (info < 0) THEN
      PRINT *, '  Argument ', -info, ' had an illegal value'
    ELSE
      PRINT *, '  Matrix is singular (U(', info, ',', info, ') = 0)'
    END IF
    J_LW = 0.0
    RETURN
  END IF
  
  ! Solution is now in b
  J_LW = b

END SUBROUTINE solve_LW_radiosities


SUBROUTINE solve_LW_radiosities_lu(N, epsilon, T, LU, ipiv, F_sky, L_down, &
                                   J_LW, info)
!-----------------------------------------------------------------------
! Same as solve_LW_radiosities, but with the matrix I - (I-diag(eps)) F
! already factorised by factor_radiosity_matrix (pvland_radiosity_SW.f90).
! [SYLVAIN]
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: N
  REAL,    INTENT(IN)  :: epsilon(N), T(N), LU(N,N), F_sky(N)
  INTEGER, INTENT(IN)  :: ipiv(N)
  REAL,    INTENT(IN)  :: L_down
  REAL,    INTENT(OUT) :: J_LW(N)
  INTEGER, INTENT(OUT) :: info
  INTEGER :: i
  REAL, PARAMETER :: sigma = 5.67e-8  ! Stefan-Boltzmann constant (W/m2/K4)

  ! Right-hand side b = [eps_i sigma T_i^4] + [I - diag(eps)] x F_s x L_down
  DO i = 1, N
    J_LW(i) = epsilon(i) * sigma * T(i)**4
    J_LW(i) = J_LW(i) + (1.0 - epsilon(i)) * F_sky(i) * L_down
  END DO

  CALL SGETRS('N', N, 1, LU, N, ipiv, J_LW, N, info)

  IF (info /= 0) THEN
    PRINT *, 'ERROR in solve_LW_radiosities_lu: SGETRS failed with info=', info
    J_LW = 0.0
  END IF

END SUBROUTINE solve_LW_radiosities_lu


SUBROUTINE compute_LW_irradiances(N, F, F_sky, J_LW, L_down, G_LW)
!-----------------------------------------------------------------------
! Compute total longwave irradiance on each surface (Eq. 5.2 from PDF):
!   G_LW(i) = sum_j F(i,j)*J_LW(j) + F_sky(i)*L_down
!
! Input:
!   N         : Number of surfaces
!   F(N,N)    : View factor matrix
!   F_sky(N)  : Sky view factors
!   J_LW(N)   : Longwave radiosities (W/mÃ‚Â²)
!   L_down    : Incoming LW radiation from atmosphere (W/mÃ‚Â²)
!
! Output:
!   G_LW(N)   : Total LW irradiance on each surface (W/mÃ‚Â²)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: F(N,N), F_sky(N)
  REAL, INTENT(IN) :: J_LW(N), L_down
  REAL, INTENT(OUT) :: G_LW(N)
  
  ! Local variables
  INTEGER :: i, j
  REAL :: sum_radiation
  
  DO i = 1, N
    ! Radiation from other surfaces
    sum_radiation = 0.0
    DO j = 1, N
      sum_radiation = sum_radiation + F(i,j) * J_LW(j)
    END DO
    
    ! Add atmospheric radiation
    G_LW(i) = sum_radiation + F_sky(i) * L_down
  END DO
  
END SUBROUTINE compute_LW_irradiances





SUBROUTINE compute_LW_net_fluxes(N, G_LW, J_LW, q_LW)
!-----------------------------------------------------------------------
! Compute net longwave radiative flux for each surface (Eq. 5.1 from PDF):
!   q_LW(i) = G_LW(i) - J_LW(i)
!
! Input:
!   N        : Number of surfaces
!   G_LW(N)  : Total LW irradiance (W/mÃ‚Â²)
!   J_LW(N)  : LW radiosity (W/mÃ‚Â²)
!
! Output:
!   q_LW(N)  : Net LW flux (W/mÃ‚Â²) [positive = gaining energy]
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: G_LW(N), J_LW(N)
  REAL, INTENT(OUT) :: q_LW(N)
  
  ! Local variables
  INTEGER :: i
  
  DO i = 1, N
    q_LW(i) = G_LW(i) - J_LW(i)
  END DO
  
END SUBROUTINE compute_LW_net_fluxes




SUBROUTINE compute_LW_net_fluxes_single_reflection(N, epsilon, T, F, F_sky, &
                                               L_down, q_LW)
!-----------------------------------------------------------------------
! Compute net longwave radiative flux with at most one reflection
!
! Implements:
! q_LW(i) = ÃŽÂµ_i [ sum_j ÃŽÂµ_j F_ij ÃÆ’ T_j^4
!               + sum_j (1-ÃŽÂµ_j) F_ij sum_k ÃŽÂµ_k F_jk ÃÆ’ T_k^4
!               + F_sky(i) L_down
!               + sum_j (1-ÃŽÂµ_j) F_ij F_sky(j) L_down
!               - ÃÆ’ T_i^4 ]
!
! Output q_LW > 0 means net gain of LW energy by surface i
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: epsilon(N), T(N)
  REAL, INTENT(IN) :: F(N,N), F_sky(N)
  REAL, INTENT(IN) :: L_down
  REAL, INTENT(OUT) :: q_LW(N)

  INTEGER :: i, j, k
  REAL :: term1, term2, term3, term4
  REAL :: sigma
  REAL :: E(N)

  sigma = 5.67e-8

  ! Precompute ÃŽÂµ ÃÆ’ T^4 for all surfaces
  DO i = 1, N
    E(i) = epsilon(i) * sigma * T(i)**4
  END DO

  DO i = 1, N

    ! (1) Direct emission from other surfaces
    term1 = 0.0
    DO j = 1, N
      term1 = term1 + F(i,j) * E(j)
    END DO

    ! (2) One reflection of surface emission
    term2 = 0.0
    DO j = 1, N
      DO k = 1, N
        term2 = term2 + (1.0 - epsilon(j)) * F(i,j) * &
                            F(j,k) * E(k)
      END DO
    END DO

    ! (3) Direct atmospheric radiation
    term3 = F_sky(i) * L_down

    ! (4) One reflection of atmospheric radiation
    term4 = 0.0
    DO j = 1, N
      term4 = term4 + (1.0 - epsilon(j)) * F(i,j) * &
                          F_sky(j) * L_down
    END DO

    ! Net LW flux
    q_LW(i) = epsilon(i) * ( term1 + term2 + term3 + term4 &
                             - sigma * T(i)**4 )

  END DO

END SUBROUTINE compute_LW_net_fluxes_single_reflection




SUBROUTINE validate_LW_conservation(N, surface_areas, q_LW, J_LW, &
                                     F_sky, total_net, &
                                     total_emitted_to_sky)
!-----------------------------------------------------------------------
! Check energy balance for LW radiation
! 
! Computes:
!   - Total net LW flux (absorbed - emitted)
!   - Total emitted to sky (energy leaving the system)
!
! Input:
!   N              : Number of surfaces
!   surface_areas(N): Area of each surface (m)
!   q_LW(N)        : Net LW flux for each surface (W/mÃ‚Â²)
!   J_LW(N)        : LW radiosity leaving each surface (W/mÃ‚Â²)
!   F_sky(N)       : Sky view factors
!
! Output:
!   total_net             : Total net LW flux (W/m)
!   total_emitted_to_sky  : Total emitted to sky (W/m)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: surface_areas(N), q_LW(N), J_LW(N), F_sky(N)
  REAL, INTENT(OUT) :: total_net, total_emitted_to_sky
  
  ! Local variables
  INTEGER :: i
  
  ! Total net LW flux (only for pattern n: surfaces 1,2,3)
  total_net = 0.0
  DO i = 1, 3
    total_net = total_net + surface_areas(i) * q_LW(i)
  END DO
  
  ! Total emitted to sky (energy leaving the system)
  total_emitted_to_sky = 0.0
  DO i = 1, 3
    total_emitted_to_sky = total_emitted_to_sky + &
                           surface_areas(i) * J_LW(i) * F_sky(i)
  END DO
  
END SUBROUTINE validate_LW_conservation


SUBROUTINE print_LW_results(N, surface_names, T, G_LW, J_LW, q_LW, q_LW_single_refl)
!-----------------------------------------------------------------------
! Print formatted results for LW radiation
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  CHARACTER(LEN=20), INTENT(IN) :: surface_names(N)
  REAL, INTENT(IN) :: T(N), G_LW(N), J_LW(N), q_LW(N), q_LW_single_refl(N)
  
  ! Local variables
  INTEGER :: i
  
  PRINT *, ''
  PRINT *, '=============================================='
  PRINT *, '     LONGWAVE RADIATION RESULTS'
  PRINT *, '=============================================='
  PRINT *, ''
  PRINT '(A20,A15,3A17, A25)', 'Surface', 'T (K)', 'G_LW (W/mÃ‚Â²)', &
                      'J_LW (W/mÃ‚Â²)', 'q_LW (W/mÃ‚Â²)', 'q_LW_single_refl (W/mÃ‚Â²)'
  PRINT *, '------------------------------------------------------------------------'
  
  DO i = 1, N
    PRINT '(A20,F15.2,4F15.2)', TRIM(surface_names(i)), T(i), &
                                 G_LW(i), J_LW(i), q_LW(i), q_LW_single_refl(i)
  END DO
  
  PRINT *, '------------------------------------------------------------------------'
  PRINT *, ''
  
END SUBROUTINE print_LW_results


SUBROUTINE check_LW_physical_validity(N, surface_names, epsilon, T, &
                                       G_LW, J_LW, q_LW, all_valid)
!-----------------------------------------------------------------------
! Check physical validity of LW radiation results
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  CHARACTER(LEN=20), INTENT(IN) :: surface_names(N)
  REAL, INTENT(IN) :: epsilon(N), T(N), G_LW(N), J_LW(N), q_LW(N)
  LOGICAL, INTENT(OUT) :: all_valid
  
  ! Local variables
  INTEGER :: i
  LOGICAL :: valid
  REAL, PARAMETER :: sigma = 5.67e-8
  REAL, PARAMETER :: tol = 1.0
  REAL :: exitance
  
  all_valid = .TRUE.
  
  PRINT *, ''
  PRINT *, '=============================================='
  PRINT *, '     PHYSICAL VALIDITY CHECKS (LW)'
  PRINT *, '=============================================='
  PRINT *, ''
  
  ! Check 1: All radiosities should be positive
  PRINT *, '1. Checking radiosities (J_LW >= 0)...'
  DO i = 1, N
    valid = (J_LW(i) >= -tol)
    IF (.NOT. valid) THEN
      PRINT '(A,A,A,F10.3)', '   FAILED for ', TRIM(surface_names(i)), &
                              ': J_LW = ', J_LW(i)
      all_valid = .FALSE.
    END IF
  END DO
  IF (all_valid) PRINT *, '   PASSED'
  
  
  ! Check 2: Energy balance q_LW = G_LW - J_LW
  PRINT *, '4. Checking q_LW = G_LW - J_LW...'
  DO i = 1, N
    valid = (ABS(q_LW(i) - (G_LW(i) - J_LW(i))) < tol)
    IF (.NOT. valid) THEN
      PRINT '(A,A)', '   FAILED for ', TRIM(surface_names(i))
      all_valid = .FALSE.
    END IF
  END DO
  IF (all_valid) PRINT *, '   PASSED'
  
  PRINT *, ''
  IF (all_valid) THEN
    PRINT *, 'ALL CHECKS PASSED Ã¢Å“â€œ'
  ELSE
    PRINT *, 'SOME CHECKS FAILED Ã¢Å“â€”'
  END IF
  PRINT *, '=============================================='
  PRINT *, ''
  
END SUBROUTINE check_LW_physical_validity
