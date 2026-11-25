!##################################################################
!##################################################################
!######                                                      ######
!######        SUBROUTINE compute_view_factor_matrix         ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_view_factor_matrix(L, a, w, h_PV, beta, F_matrix)

!=========================================================================
! SUBROUTINE: compute_view_factor_matrix
!
! PURPOSE:
!   Computes the full 8x8 view factor matrix for radiative heat transfer
!   between solar panel surfaces and ground surfaces in patterns n-1, n, n+1
!
! INPUTS:
!   L       (REAL) - Panel length [m]
!   a       (REAL) - Spacing between panels [m]
!   w       (REAL) - Horizontal projection of the panel, L cos(beta) [m]
!                       (the ground surface of a pattern is a + w wide)
!   h_PV    (REAL) - Height of the panel top edge [m]
!   beta    (REAL) - Tilt angle [radians]
!
! OUTPUTS:
!   F_matrix (REAL(8,8)) - View factor matrix where:
!     Row/Col 1: nat,n     (ground surface of pattern n)
!     Row/Col 2: PVf,n     (front of panel n)
!     Row/Col 3: PVb,n     (back of panel n)
!     Row/Col 4: PVb,n-1   (back of panel n-1)
!     Row/Col 5: PVf,n-1   (front of panel n-1)
!     Row/Col 6: PVb,n+1   (back of panel n+1)
!     Row/Col 7: nat,n-1   (ground surface of pattern n-1)
!     Row/Col 8: nat,n+1   (ground surface of pattern n+1)
!
!=========================================================================

    IMPLICIT NONE
    REAL, INTENT(IN) :: L           ! Panel length
    REAL, INTENT(IN) :: a           ! Spacing between panels
    REAL, INTENT(IN) :: w           ! Horizontal projection of the panel
    REAL, INTENT(IN) :: h_PV        ! Height of the panel top edge
    REAL, INTENT(IN) :: beta        ! Tilt angle (radians)
    REAL, INTENT(OUT) :: F_matrix(8,8)  ! View factor matrix
    
    REAL :: h_clear
    REAL :: F_PVb_to_PVf, F_nat_to_PVb_n, F_nat_to_PVf_nm1
    REAL :: F_nat_to_PVb_nm1, F_nat_to_PVb_np1
    
    ! Initialize matrix to zero
    F_matrix = 0.0
    
    ! Compute clearance height
    h_clear = h_PV - L * SIN(beta)

    ! =====================================================================
    ! View factors between panels
    ! =====================================================================
    F_PVb_to_PVf = (SQRT(a**2 + (L*SIN(beta))**2) + &
                    SQRT((2*w + a)**2 + (L*SIN(beta))**2) - &
                    2*(w + a)) / (2*L)
    
    ! Symmetry: F(PVb,n -> PVf,n-1) = F(PVf,n-1 -> PVb,n) = F(PVf,n -> PVb,n+1) = F(PVb,n+1 -> PVf,n)
    F_matrix(3, 5) = F_PVb_to_PVf  ! PVb,n -> PVf,n-1
    F_matrix(5, 3) = F_PVb_to_PVf  ! PVf,n-1 -> PVb,n
    F_matrix(2, 6) = F_PVb_to_PVf  ! PVf,n -> PVb,n+1
    F_matrix(6, 2) = F_PVb_to_PVf  ! PVb,n+1 -> PVf,n
    
    ! =====================================================================
    ! Ground of pattern n to back of panel n
    ! =====================================================================
    CALL compute_F_nat_to_PVb_n(L, a, w, h_clear, beta, F_nat_to_PVb_n)
    F_matrix(1, 3) = F_nat_to_PVb_n  ! nat,n -> PVb,n
    F_matrix(3, 1) = (a + w) / L * F_nat_to_PVb_n  ! PVb,n -> nat,n (reciprocity)
    
    ! =====================================================================
    ! Ground of pattern n to front of panel n-1
    ! =====================================================================
    CALL compute_F_nat_to_PVf_nm1(L, a, w, h_clear, beta, F_nat_to_PVf_nm1)
    F_matrix(1, 5) = F_nat_to_PVf_nm1  ! nat,n -> PVf,n-1
    F_matrix(8, 2) = F_nat_to_PVf_nm1  ! nat,n+1 -> PVf,n (periodicity)
    F_matrix(5, 1) = (a + w) / L * F_nat_to_PVf_nm1  ! PVf,n-1 -> nat,n (reciprocity)
    F_matrix(2, 8) = (a + w) / L * F_nat_to_PVf_nm1  ! PVf,n -> nat,n+1 (reciprocity)
    
    ! =====================================================================
    ! Ground of pattern n to back of panel n-1
    ! =====================================================================
    CALL compute_F_nat_to_PVb_nm1(L, a, w, h_clear, beta, F_nat_to_PVb_nm1)
    F_matrix(1, 4) = F_nat_to_PVb_nm1  ! nat,n -> PVb,n-1
    F_matrix(8, 3) = F_nat_to_PVb_nm1  ! nat,n+1 -> PVb,n (periodicity)
    F_matrix(4, 1) = (a + w) / L * F_nat_to_PVb_nm1  ! PVb,n-1 -> nat,n (reciprocity)
    F_matrix(3, 8) = (a + w) / L * F_nat_to_PVb_nm1  ! PVb,n -> nat,n+1 (reciprocity)
    
    ! =====================================================================
    ! Ground of pattern n to back of panel n+1
    ! =====================================================================
    CALL compute_F_nat_to_PVb_np1(L, a, w, h_PV, h_clear, beta, F_nat_to_PVb_np1)
    F_matrix(1, 6) = F_nat_to_PVb_np1  ! nat,n -> PVb,n+1
    F_matrix(7, 3) = F_nat_to_PVb_np1  ! nat,n-1 -> PVb,n (periodicity)
    F_matrix(6, 1) = (a + w) / L * F_nat_to_PVb_np1  ! PVb,n+1 -> nat,n (reciprocity)
    F_matrix(3, 7) = (a + w) / L * F_nat_to_PVb_np1  ! PVb,n -> nat,n-1 (reciprocity)
    
    ! The rows of the peripheral surfaces 4-8 are incomplete (only the factors
    ! towards pattern n are filled): only the surfaces of pattern n (1-3) have
    ! a complete row, see compute_sky_view_factors.

END SUBROUTINE compute_view_factor_matrix



!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE compute_F_nat_to_PVb_n          ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_F_nat_to_PVb_n(L, a, w, h_clear, beta, F)

!=========================================================================
! SUBROUTINE: compute_F_nat_to_PVb_n
!
! PURPOSE:
!   Computes the view factor from ground surface of pattern n to the
!   back surface of panel n
!
! INPUTS:
!   L        (REAL) - Panel length [m]
!   a        (REAL) - Spacing between panels [m]
!   w        (REAL) - Horizontal projection of the panel, L cos(beta) [m]
!   h_clear  (REAL) - Clearance height [m]
!   beta     (REAL) - Tilt angle [radians]
!
! OUTPUTS:
!   F        (REAL) - View factor from nat,n to PVb,n
!
!=========================================================================

    IMPLICIT NONE
    REAL, INTENT(IN) :: L, a, w, h_clear, beta
    REAL, INTENT(OUT) :: F
    REAL :: term1, term2, term3, term4
    REAL, PARAMETER :: PI = 3.14159265
    
    term1 = SQRT( L*L + 2.0*L*h_clear*SIN(beta) + h_clear*h_clear )
    term2 = SQRT((a+w)*(a+w) + h_clear*h_clear)
    term3 = SQRT((a+w)*(a+w) + h_clear*h_clear + L*L + 2*h_clear*L*SIN(beta) &
        		- 2*(a+w)*L*COS(beta))
    term4 = ABS(h_clear)

    F = ( term1 + term2 - term3 - term4 ) / ( 2.0*(a + w) )

	F = MAX(0.0, F)

END SUBROUTINE compute_F_nat_to_PVb_n




!##################################################################
!##################################################################
!######                                                      ######
!######         SUBROUTINE compute_F_nat_to_PVf_nm1          ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_F_nat_to_PVf_nm1(L, a, w, h_clear, beta, F)

!=========================================================================
! SUBROUTINE: compute_F_nat_to_PVf_nm1
!
! PURPOSE:
!   Computes the view factor from ground surface of pattern n to the
!   front surface of panel n-1
!
! INPUTS:
!   L        (REAL) - Panel length [m]
!   a        (REAL) - Spacing between panels [m]
!   w        (REAL) - Horizontal projection of the panel, L cos(beta) [m]
!   h_clear  (REAL) - Clearance height [m]
!   beta     (REAL) - Tilt angle [radians]
!
! OUTPUTS:
!   F        (REAL) - View factor from nat,n to PVf,n-1
!
!=========================================================================

    IMPLICIT NONE
    REAL, INTENT(IN) :: L, a, w, h_clear, beta
    REAL, INTENT(OUT) :: F
    REAL :: term1, term2, term3, delta_x0
    REAL, PARAMETER :: PI = 3.14159265
    REAL, PARAMETER :: eps = 0.01

    ! Handle special case: beta = 0 (horizontal panels)
    IF (ABS(beta) < eps) THEN
        F = 0.0
        RETURN
    END IF
    
    ! Compute delta_x0 locally
    delta_x0 = h_clear / TAN(beta)
    
    ! When delta_x0 >= a+w, the view factor vanishes
    IF (ABS(delta_x0) > a+w-eps) THEN
        F = 0.0
    ELSE
		term1 = L
		term2 = SQRT( (a + w)*(a + w) + h_clear*h_clear )
		term3 = SQRT( (a + w)*(a + w) + L*L + h_clear*h_clear &
		              + 2.0*L*h_clear*SIN(beta) &
		              + 2.0*(a + w)*L*COS(beta))
		F = (term1 + term2 - term3) / (2.0*(a+w))
	END IF

	F = MAX(0.0, F)

END SUBROUTINE compute_F_nat_to_PVf_nm1




!##################################################################
!##################################################################
!######                                                      ######
!######         SUBROUTINE compute_F_nat_to_PVb_nm1          ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_F_nat_to_PVb_nm1(L, a, w, h_clear, beta, F)

!=========================================================================
! SUBROUTINE: compute_F_nat_to_PVb_nm1
!
! PURPOSE:
!   Computes the view factor from ground surface of pattern n to the
!   back surface of panel n-1
!
! INPUTS:
!   L        (REAL) - Panel length [m]
!   a        (REAL) - Spacing between panels [m]
!   w        (REAL) - Horizontal projection of the panel, L cos(beta) [m]
!   h_clear  (REAL) - Clearance height [m]
!   beta     (REAL) - Tilt angle [radians]
!
! OUTPUTS:
!   F        (REAL) - View factor from nat,n to PVb,n-1
!
!=========================================================================

    IMPLICIT NONE
    REAL, INTENT(IN) :: L, a, w, h_clear, beta
    REAL, INTENT(OUT) :: F
    REAL :: term1, term2, term3, term4, delta_x0
    REAL, PARAMETER :: PI = 3.14159265
	REAL, PARAMETER :: eps = 0.01

    ! Handle special case: beta = 0 (horizontal panels)
    IF (ABS(beta) < eps) THEN
        term1 = SQRT(L*L+(a+w)*(a+w) &
					+ h_clear*h_clear + 2.0*L*(a+w))
		term2 = h_clear
		term3 = SQRT(L*L+ h_clear*h_clear)
		term4 = SQRT(h_clear*h_clear + (a+w)*(a+w))
		F = (term1 + term2 - term3 - term4) / (2.0*(a+w))
		F = MAX(0.0, F)
        RETURN
    END IF
    
    ! Compute delta_x0 locally
    delta_x0 = h_clear / TAN(beta)	

	! When delta_x0 >= a+w, other formula
    IF (ABS(delta_x0) < a+w-eps) THEN
        term1 = SQRT( L*L + 2.0*L*h_clear*SIN(beta) + h_clear*h_clear )
    	F = ( L + h_clear - term1 ) / ( 2.0*(a + w) )
    ELSE
		term1 = SQRT(L*L+(a+w)*(a+w)+2.0*L*h_clear*SIN(beta) &
					+ h_clear*h_clear + 2.0*L*(a+w)*COS(beta))
		term2 = h_clear
		term3 = SQRT(L*L+2.0*L*h_clear*SIN(beta)+ h_clear*h_clear)
		term4 = SQRT(h_clear*h_clear + (a+w)*(a+w))
		F = (term1 + term2 - term3 - term4) / (2.0*(a+w))
    END IF

	F = MAX(0.0, F)

END SUBROUTINE compute_F_nat_to_PVb_nm1




!##################################################################
!##################################################################
!######                                                      ######
!######         SUBROUTINE compute_F_nat_to_PVb_np1          ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_F_nat_to_PVb_np1(L, a, w, h_PV, h_clear, beta, F)

!=========================================================================
! SUBROUTINE: compute_F_nat_to_PVb_np1
!
! PURPOSE:
!   Computes the view factor from ground surface of pattern n to the
!   back surface of panel n+1
!
! INPUTS:
!   L        (REAL) - Panel length [m]
!   a        (REAL) - Spacing between panels [m]
!   w        (REAL) - Horizontal projection of the panel, L cos(beta) [m]
!   h_PV     (REAL) - Height of the panel top edge [m]
!   h_clear  (REAL) - Clearance height [m]
!   beta     (REAL) - Tilt angle [radians]
!
! OUTPUTS:
!   F        (REAL) - View factor from nat,n to PVb,n+1
!
!=========================================================================

    IMPLICIT NONE
    REAL, INTENT(IN) :: L, a, w, h_PV, h_clear, beta
    REAL, INTENT(OUT) :: F
    REAL :: term1, term2, term3, term4
    REAL :: h_sin, L_h_sin, x1, x2, delta_x0, delta_x1
    REAL, PARAMETER :: eps = 0.01
    
    IF (ABS(beta) < eps) THEN
        delta_x1 = HUGE(1.0)
    ELSE
        delta_x1 = a * h_clear / (h_PV - h_clear)
    END IF
    term1 = SQRT(a*a + h_PV*h_PV)
    term2 = SQRT(4.0*(a+w)**2 + h_clear*h_clear)
    term4 = SQRT((a+w)**2 + h_clear*h_clear)
    IF (delta_x1 < a + w) THEN
        term3 = term4 + SQRT(a*a + (h_PV - h_clear)**2)
    ELSE
        term3 = SQRT((2.0*a + w)**2 + h_PV*h_PV)
    END IF
    F = MAX(0.0, (term1 + term2 - term3 - term4) / (2.0*(a + w)))

END SUBROUTINE compute_F_nat_to_PVb_np1




!##################################################################
!##################################################################
!######                                                      ######
!######         SUBROUTINE compute_sky_view_factors          ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_sky_view_factors(F_matrix, F_sky)

!=========================================================================
! SUBROUTINE: compute_sky_view_factors
!
! PURPOSE:
!   Computes the sky-view factors for all 8 surfaces. The sky-view factor
!   is the fraction of the hemisphere not obstructed by other surfaces,
!   calculated as: F_sky(i) = 1 - sum_j F(i->j). Here, only the sky-view 
!   factors of surfaces of pattern n are correct.
!
! INPUTS:
!   F_matrix (REAL(8,8)) - Complete view factor matrix between surfaces
!
! OUTPUTS:
!   F_sky    (REAL)   - Sky-view factors for each surface
!
!=========================================================================

    IMPLICIT NONE
    REAL, INTENT(IN) :: F_matrix(8,8)
    REAL, INTENT(OUT) :: F_sky(8)
    INTEGER :: i
    
    DO i = 1, 8
        F_sky(i) = 1.0 - SUM(F_matrix(i,:))
    END DO
    
END SUBROUTINE compute_sky_view_factors
