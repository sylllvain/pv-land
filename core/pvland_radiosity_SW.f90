!=======================================================================
! SW_RADIATION_SUBROUTINES.F90
!
! Subroutines for shortwave radiation calculation using radiosity method
! Based on PV-LAND model (Section 4)
! SINGLE PRECISION VERSION
!
! Author: Academic research code
! Date: 2026
!=======================================================================

SUBROUTINE solve_SW_radiosities(N, alpha, F, F_sky, S_direct, &
                                 S_diffuse, J_SW, info)
!-----------------------------------------------------------------------
! Solve the linear system for shortwave radiosities.
! 
! From PDF Eq. 4.1: J_SW,i = ÃŽÂ±_i Ãƒâ€” G_SW,i
! From PDF Eq. 4.2: G_SW,i = SÃ¢â€¡â€œ_i + F_isÃƒâ€”SÃ¢â€ â€œ + ÃŽÂ£_j F_iÃ¢â€ â€™j Ãƒâ€” J_SW,j
!
! Substituting (4.1) into (4.2) and rearranging:
!   [I - diag(ÃŽÂ±)Ãƒâ€”F] Ãƒâ€” J_SW = diag(ÃŽÂ±) Ãƒâ€” [SÃ¢â€¡â€œ + F_sÃƒâ€”SÃ¢â€ â€œ]
!
! Uses LAPACK routine DGESV (Double precision General matrix SolVe)
! DGESV solves A*X = B using LU decomposition with partial pivoting
!
! Input:
!   N         : Number of surfaces
!   alpha(N)  : Albedo of each surface (0 to 1, 0=black, 1=mirror)
!   F(N,N)    : View factor matrix F(i,j) = view factor from i to j
!   F_sky(N)  : Sky view factor of each surface
!   S_direct(N): Direct SW irradiance on each surface (W/mÃ‚Â²)
!   S_diffuse : Diffuse SW irradiance from sky (W/mÃ‚Â²)
!
! Output:
!   J_SW(N)   : Shortwave radiosity of each surface (W/mÃ‚Â²)
!   info      : 0 if successful, /=0 if error
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: alpha(N), F(N,N), F_sky(N)
  REAL, INTENT(IN) :: S_direct(N), S_diffuse
  REAL, INTENT(OUT) :: J_SW(N)
  INTEGER, INTENT(OUT) :: info
  
  ! Local variables
  REAL :: A(N,N), b(N)
  INTEGER :: i, j
  INTEGER :: ipiv(N)  ! Pivot indices for DGESV
  
  ! Initialize
  info = 0
  A = 0.0
  b = 0.0
  
  ! Build matrix A = I - diag(alpha) Ãƒâ€” F
  ! Start with identity matrix
  DO i = 1, N
    A(i,i) = 1.0
  END DO
  
  ! Subtract diag(alpha) Ãƒâ€” F
  DO i = 1, N
    DO j = 1, N
      A(i,j) = A(i,j) - alpha(i) * F(i,j)
    END DO
  END DO
  
  ! Build right-hand side b = diag(alpha) Ãƒâ€” [S_direct + F_sky Ãƒâ€” S_diffuse]
  DO i = 1, N
    b(i) = alpha(i) * (S_direct(i) + F_sky(i) * S_diffuse)
  END DO
  
  ! Solve linear system A Ãƒâ€” J_SW = b using LAPACK DGESV
  ! DGESV arguments:
  !   N      : Order of matrix A
  !   1      : Number of right-hand sides (NRHS)
  !   A      : Matrix (will be overwritten with LU factorization)
  !   N      : Leading dimension of A (LDA)
  !   ipiv   : Pivot indices (output)
  !   b      : Right-hand side on input, solution on output
  !   N      : Leading dimension of b (LDB)
  !   info   : 0 on success, <0 if argument error, >0 if singular
  !CALL DGESV(N, 1, A, N, ipiv, b, N, info)
  CALL SGESV(N, 1, A, N, ipiv, b, N, info) !modif sylvain
  
  IF (info /= 0) THEN
    PRINT *, 'ERROR in solve_SW_radiosities: DGESV failed with info=', info
    IF (info < 0) THEN
      PRINT *, '  Argument ', -info, ' had an illegal value'
    ELSE
      PRINT *, '  Matrix is singular (U(', info, ',', info, ') = 0)'
    END IF
    J_SW = 0.0
    RETURN
  END IF
  
  ! Solution is now in b
  J_SW = b
  
END SUBROUTINE solve_SW_radiosities


SUBROUTINE compute_SW_irradiances(N, F, F_sky, S_direct, S_diffuse, &
                                   J_SW, G_SW)
!-----------------------------------------------------------------------
! Compute total shortwave irradiance on each surface (Eq. 4.2):
!   G_SW(i) = S_direct(i) + F_sky(i)*S_diffuse + sum_j F(i,j)*J_SW(j)
!
! Input:
!   N         : Number of surfaces
!   F(N,N)    : View factor matrix
!   F_sky(N)  : Sky view factors
!   S_direct(N): Direct SW irradiance (W/mÃ‚Â²)
!   S_diffuse : Diffuse SW from sky (W/mÃ‚Â²)
!   J_SW(N)   : Shortwave radiosities (W/mÃ‚Â²)
!
! Output:
!   G_SW(N)   : Total SW irradiance on each surface (W/mÃ‚Â²)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: F(N,N), F_sky(N)
  REAL, INTENT(IN) :: S_direct(N), S_diffuse, J_SW(N)
  REAL, INTENT(OUT) :: G_SW(N)
  
  ! Local variables
  INTEGER :: i, j
  REAL :: sum_reflected
  
  DO i = 1, N
    ! Direct + diffuse from sky
    G_SW(i) = S_direct(i) + F_sky(i) * S_diffuse
    
    ! Add reflected radiation from other surfaces
    sum_reflected = 0.0
    DO j = 1, N
      sum_reflected = sum_reflected + F(i,j) * J_SW(j)
    END DO
    
    G_SW(i) = G_SW(i) + sum_reflected
  END DO
  
END SUBROUTINE compute_SW_irradiances


SUBROUTINE compute_SW_net_fluxes(N, alpha, G_SW, q_SW)
!-----------------------------------------------------------------------
! Compute net shortwave radiative flux absorbed by each surface (Eq. 4.4):
!   q_SW(i) = (1 - alpha(i)) * G_SW(i)
!
! Input:
!   N        : Number of surfaces
!   alpha(N) : Albedo of each surface
!   G_SW(N)  : Total SW irradiance (W/mÃ‚Â²)
!
! Output:
!   q_SW(N)  : Net SW flux absorbed (W/mÃ‚Â²)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: alpha(N), G_SW(N)
  REAL, INTENT(OUT) :: q_SW(N)
  
  ! Local variables
  INTEGER :: i
  
  DO i = 1, N
    q_SW(i) = (1.0 - alpha(i)) * G_SW(i)
  END DO
  
END SUBROUTINE compute_SW_net_fluxes



SUBROUTINE validate_SW_conservation(N, surface_areas, q_SW, J_SW, &
                                     F_sky, total_absorbed, &
                                     total_reflected_to_sky)
!-----------------------------------------------------------------------
! Check energy conservation for SW radiation
! 
! Energy balance: absorbed + reflected to sky = incident from sky + direct
! Here we compute:
!   - Total absorbed by all surfaces
!   - Total reflected back to sky (lost from system)
!
! Input:
!   N              : Number of surfaces
!   surface_areas(N): Area of each surface (m)
!   q_SW(N)        : Net SW flux absorbed by each surface (W/mÃ‚Â²)
!   J_SW(N)        : SW radiosity leaving each surface (W/mÃ‚Â²)
!   F_sky(N)       : Sky view factors
!
! Output:
!   total_absorbed        : Total energy absorbed by all surfaces (W/m)
!   total_reflected_to_sky: Total energy reflected to sky (W/m)
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  REAL, INTENT(IN) :: surface_areas(N), q_SW(N), J_SW(N), F_sky(N)
  REAL, INTENT(OUT) :: total_absorbed, total_reflected_to_sky
  
  ! Local variables
  INTEGER :: i
  
  ! Total absorbed energy (only for pattern n: surfaces 1,2,3)
  total_absorbed = 0.0
  DO i = 1, 3
    total_absorbed = total_absorbed + surface_areas(i) * q_SW(i)
  END DO
  
  ! Total reflected to sky (energy leaving the system)
  total_reflected_to_sky = 0.0
  DO i = 1, 3
    total_reflected_to_sky = total_reflected_to_sky + &
                             surface_areas(i) * J_SW(i) * F_sky(i)
  END DO
  
END SUBROUTINE validate_SW_conservation


SUBROUTINE print_SW_results(N, surface_names, G_SW, J_SW, q_SW)
!-----------------------------------------------------------------------
! Print formatted results for SW radiation
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  CHARACTER(LEN=20), INTENT(IN) :: surface_names(N)
  REAL, INTENT(IN) :: G_SW(N), J_SW(N), q_SW(N)
  
  ! Local variables
  INTEGER :: i
  
  PRINT *, ''
  PRINT *, '=============================================='
  PRINT *, '     SHORTWAVE RADIATION RESULTS'
  PRINT *, '=============================================='
  PRINT *, ''
  PRINT '(A20,3A15)', 'Surface', 'G_SW (W/mÃ‚Â²)', 'J_SW (W/mÃ‚Â²)', 'q_SW (W/mÃ‚Â²)'
  PRINT *, '--------------------------------------------------------------'
  
  DO i = 1, N
    PRINT '(A20,3F15.2)', TRIM(surface_names(i)), G_SW(i), J_SW(i), q_SW(i)
  END DO
  
  PRINT *, '--------------------------------------------------------------'
  PRINT *, ''
  
END SUBROUTINE print_SW_results


SUBROUTINE check_SW_physical_validity(N, surface_names, alpha, &
                                       G_SW, J_SW, q_SW, &
                                       S_total, all_valid)
!-----------------------------------------------------------------------
! Check physical validity of SW radiation results
!-----------------------------------------------------------------------
  IMPLICIT NONE
  
  ! Arguments
  INTEGER, INTENT(IN) :: N
  CHARACTER(LEN=20), INTENT(IN) :: surface_names(N)
  REAL, INTENT(IN) :: alpha(N), G_SW(N), J_SW(N), q_SW(N)
  REAL, INTENT(IN) :: S_total
  LOGICAL, INTENT(OUT) :: all_valid
  
  ! Local variables
  INTEGER :: i
  LOGICAL :: valid
  REAL, PARAMETER :: tol = 1.0e-6
  
  all_valid = .TRUE.
  
  PRINT *, ''
  PRINT *, '=============================================='
  PRINT *, '     PHYSICAL VALIDITY CHECKS'
  PRINT *, '=============================================='
  PRINT *, ''
  
  ! Check 1: All radiosities should be non-negative
  PRINT *, '1. Checking radiosities (J_SW >= 0)...'
  DO i = 1, N
    valid = (J_SW(i) >= -tol)
    IF (.NOT. valid) THEN
      PRINT '(A,A,A,F10.3)', '   FAILED for ', TRIM(surface_names(i)), &
                              ': J_SW = ', J_SW(i)
      all_valid = .FALSE.
    END IF
  END DO
  IF (all_valid) PRINT *, '   PASSED'
  
  ! Check 2: Irradiances should not exceed reasonable bounds
  PRINT *, '2. Checking irradiances (G_SW <= 3*S_total)...'
  DO i = 1, N
    valid = (G_SW(i) <= 3.0 * S_total + tol)
    IF (.NOT. valid) THEN
      PRINT '(A,A,A,F10.2,A,F10.2)', '   WARNING for ', &
                              TRIM(surface_names(i)), &
                              ': G_SW = ', G_SW(i), ' > 3*S_total = ', &
                              3.0*S_total
      ! Not a hard failure, just warning
    END IF
  END DO
  PRINT *, '   PASSED (warnings above if any)'
  
  ! Check 3: Absorbed flux should match (1-alpha)*G_SW
  PRINT *, '3. Checking absorbed flux consistency...'
  DO i = 1, N
    valid = (ABS(q_SW(i) - (1.0 - alpha(i)) * G_SW(i)) < tol * S_total)
    IF (.NOT. valid) THEN
      PRINT '(A,A)', '   FAILED for ', TRIM(surface_names(i))
      all_valid = .FALSE.
    END IF
  END DO
  IF (all_valid) PRINT *, '   PASSED'
  
  ! Check 4: Radiosity should equal alpha * G_SW (from Eq. 4.1)
  PRINT *, '4. Checking radiosity definition (J_SW = alpha*G_SW)...'
  DO i = 1, N
    valid = (ABS(J_SW(i) - alpha(i) * G_SW(i)) < tol * S_total)
    IF (.NOT. valid) THEN
      PRINT '(A,A,A,F10.2,A,F10.2)', '   FAILED for ', &
                              TRIM(surface_names(i)), &
                              ': J_SW = ', J_SW(i), &
                              ' vs alpha*G_SW = ', alpha(i)*G_SW(i)
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
  
END SUBROUTINE check_SW_physical_validity
