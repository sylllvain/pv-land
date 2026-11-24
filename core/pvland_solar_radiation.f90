!===============================================================================
! FILE: solar_radiation.f90
!
! Collection of subroutines to compute visibility factors and direct solar
! radiation received by ground and PV panels in a periodic PV array.
!===============================================================================


!===============================================================================
! SUBROUTINE: compute_f_ground
!
! Purpose: Compute the visibility factor for the ground surface (f_nat),
!          using the projection and union of intervals method.
!
! Inputs:
!   L      - Panel length [m]
!   a      - Inter-row spacing (gap between panels) [m]
!   beta   - Panel tilt angle from horizontal [rad]
!   h_PV   - Panel top height (altitude of upper edge) [m]
!   lambda - Solar inclination angle from vertical [rad]
!   phi    - Relative solar azimuth angle [rad]
!
! Outputs:
!   f_nat  - Ground visibility factor (illuminated fraction) [-]
!===============================================================================
SUBROUTINE compute_f_ground(L, a, beta, h_PV, lambda, phi, f_nat)

  IMPLICIT NONE

  ! Mathematical constant
  REAL, PARAMETER :: PI = 3.141592653589793
  
  !-------------------------
  ! Input variables
  !-------------------------
  REAL, INTENT(IN) :: L        ! Panel length [m]
  REAL, INTENT(IN) :: a        ! Inter-row spacing [m]
  REAL, INTENT(IN) :: beta     ! Panel tilt angle from horizontal [rad]
  REAL, INTENT(IN) :: h_PV     ! Panel top height [m]
  REAL, INTENT(IN) :: lambda   ! Solar inclination from vertical [rad]
  REAL, INTENT(IN) :: phi      ! Relative solar azimuth angle [rad]

  !-------------------------
  ! Output variable
  !-------------------------
  REAL, INTENT(OUT) :: f_nat   ! Ground visibility factor [-]

  !-------------------------
  ! Local variables
  !-------------------------
  REAL :: w            ! Panel horizontal projection (width) [m]
  REAL :: H            ! Panel vertical extent (height) [m]
  REAL :: h_clear      ! Clearance between ground and panel bottom [m]
  REAL :: period       ! Horizontal period of the array (w + a) [m]
  
  REAL :: tan_omega    ! Tangent of the profile angle
  REAL :: x_bot        ! Shadow of the bottom edge (n=0) [m]
  REAL :: x_top        ! Shadow of the top edge (n=0) [m]
  REAL :: s_min, s_max ! Bounds of the shadow for n=0 [m]
  REAL :: s_min_n, s_max_n ! Bounds for row n [m]
  
  REAL :: x_sample     ! Current ground point being tested [m]
  INTEGER :: n         ! Row index for loop
  INTEGER :: i         ! Sample index
  INTEGER, PARAMETER :: num_samples = 10000 ! Precision of the integration
  INTEGER :: shaded_count
  LOGICAL :: is_shaded
  
  REAL, PARAMETER :: lambda_cutoff = 85.0 * (PI / 180.0)

  !-------------------------
  ! Panel geometry
  !-------------------------
  w = L * COS(beta)
  H = L * SIN(beta)
  h_clear = h_PV - H
  period = w + a

  !-------------------------
  ! Physical Cut-off
  !-------------------------
  IF (ABS(lambda) >= lambda_cutoff) THEN
     f_nat = 0.0
     RETURN
  END IF

  !-------------------------
  ! Shadow Projections (for n=0)
  !-------------------------
  ! Profile angle tangent: tan(lambda) * sin(phi)
  tan_omega = TAN(lambda) * SIN(phi)
  
  ! Coordinates of shadows for row n=0
  ! Panel tilted to the left: pivot at 0, top at -w
  x_bot = -h_clear * tan_omega
  x_top = -w - h_PV * tan_omega
  
  s_min = MIN(x_bot, x_top)
  s_max = MAX(x_bot, x_top)

  !-------------------------
  ! Visibility factor (exact union of intervals)
  !-------------------------
  ! Each row casts on the ground a shadow [s_min, s_max] of length
  ! s_max - s_min, repeated with the array period. The shaded fraction
  ! of one ground period is therefore exactly MIN(1, (s_max-s_min)/period),
  ! whatever the shadow offset. This replaces the former sampling of
  ! 10000 points x 7 rows (same result to 1e-4, ~86% of the PV-LAND
  ! CPU time), which also missed shadows shifted by more than 3 periods
  ! when the sun is very low.  [SYLVAIN]
  f_nat = 1.0 - MIN(1.0, (s_max - s_min) / period)

  !-------------------------
  ! Former sampling method, kept for reference
  !-------------------------
  ! shaded_count = 0
  ! DO i = 0, num_samples - 1
  !    x_sample = (REAL(i) / REAL(num_samples - 1)) * period
  !    is_shaded = .FALSE.
  !    DO n = -3, 3
  !       s_min_n = s_min + REAL(n) * period
  !       s_max_n = s_max + REAL(n) * period
  !       IF (x_sample >= s_min_n .AND. x_sample <= s_max_n) THEN
  !          is_shaded = .TRUE.
  !          EXIT
  !       END IF
  !    END DO
  !    IF (is_shaded) shaded_count = shaded_count + 1
  ! END DO
  ! f_nat = REAL(num_samples - shaded_count) / REAL(num_samples)

END SUBROUTINE compute_f_ground


!===============================================================================
! SUBROUTINE: compute_f_PVb
!
! Purpose: Compute the visibility factor for the back side of a PV panel
!===============================================================================
SUBROUTINE compute_f_PVb(L, a, beta, h_PV, lambda, phi, f_PVb)

  IMPLICIT NONE

  ! Numerical threshold to avoid division by zero
  REAL, PARAMETER :: EPS_TOL = 1.0e-10

  !-------------------------
  ! Input variables
  !-------------------------
  REAL, INTENT(IN) :: L        ! Panel length [m]
  REAL, INTENT(IN) :: a        ! Inter-row spacing [m]
  REAL, INTENT(IN) :: beta     ! Panel tilt angle [rad]
  REAL, INTENT(IN) :: h_PV     ! Panel top height [m]
  REAL, INTENT(IN) :: lambda   ! Solar inclination from vertical [rad]
  REAL, INTENT(IN) :: phi      ! Relative solar azimuth angle [rad]

  !-------------------------
  ! Output variable
  !-------------------------
  REAL, INTENT(OUT) :: f_PVb   ! Back-side visibility factor [-]

  !-------------------------
  ! Local variables
  !-------------------------
  REAL :: w           ! Panel horizontal projection (width) [m]
  REAL :: H           ! Panel vertical extent (height) [m]
  REAL :: x_int_b     ! Intersection abscissa of shading ray [m]
  REAL :: numerator   ! Numerator of intersection expression
  REAL :: denominator ! Denominator of intersection expression

  !-------------------------
  ! Panel geometry
  !-------------------------
  w = L * COS(beta)
  H = L * SIN(beta)

  !-------------------------
  ! Degenerate case: rays nearly vertical or in cross-section
  !-------------------------
  IF (ABS(TAN(lambda) * SIN(phi)) < EPS_TOL) THEN
    f_PVb = 0.0

  ! check if the panel is vertical (threshold 1.5 near pi/2)
  ELSE IF (beta >= 1.5) THEN
	! in this case, the same visibility factor is set for both faces
	! and the dot product makes one of them to vanish
	f_PVb = ABS(a / (L * TAN(lambda)))
	f_PVb = MIN(1.0, MAX(0.0, f_PVb))
  ELSE
	!-------------------------
	! Intersection calculation
	!-------------------------
	numerator   = H * (w + a) / w
	denominator = 1.0 / (TAN(lambda) * SIN(phi)) + H / w

	x_int_b = numerator / denominator

	!-------------------------
	! Visibility factor
	!-------------------------
	f_PVb = (x_int_b - (w + a)) / w
	f_PVb = MIN(1.0, MAX(0.0, f_PVb)) 
  END IF



END SUBROUTINE compute_f_PVb


!===============================================================================
! SUBROUTINE: compute_f_PVf
!
! Purpose: Compute the visibility factor for the front side of a PV panel
!===============================================================================
SUBROUTINE compute_f_PVf(L, a, beta, h_PV, lambda, phi, f_PVf)

  IMPLICIT NONE

  ! Numerical threshold to avoid division by zero
  REAL, PARAMETER :: EPS_TOL = 1.0e-10

  !-------------------------
  ! Input variables
  !-------------------------
  REAL, INTENT(IN) :: L, a, beta, h_PV, lambda, phi

  !-------------------------
  ! Output variable
  !-------------------------
  REAL, INTENT(OUT) :: f_PVf   ! Front-side visibility factor [-]

  !-------------------------
  ! Local variables
  !-------------------------
  REAL :: w           ! Panel horizontal projection (width) [m]
  REAL :: H           ! Panel vertical extent (height) [m]
  REAL :: x_int_f     ! Intersection abscissa of shading ray [m]
  REAL :: numerator   ! Numerator of intersection expression
  REAL :: denominator ! Denominator of intersection expression

  !-------------------------
  ! Panel geometry
  !-------------------------
  w = L * COS(beta)
  H = L * SIN(beta)

  !-------------------------
  ! Degenerate case
  !-------------------------
  IF (ABS(TAN(lambda) * SIN(phi)) < EPS_TOL) THEN
    f_PVf = 1.0

  ! check if the panel is vertical (threshold 1.5 near pi/2)
  ELSE IF (beta >= 1.5) THEN
	! in this case, the same visibility factor is set for both faces
	! and the dot product makes one of them to vanish
	f_PVf = ABS(a / (L * TAN(lambda)))
	f_PVf = MIN(1.0, MAX(0.0, f_PVf))
  ELSE
	!-------------------------
	! Intersection calculation
	!-------------------------
	numerator = 2.0 * (w + a) / (TAN(lambda) * SIN(phi)) + H * (w + a) / w
	denominator = 1.0 / (TAN(lambda) * SIN(phi)) + H / w

	x_int_f = numerator / denominator
	! PRINT *, ' x_int_f    = ', x_int_f

	!-------------------------
	! Visibility factor
	!-------------------------
	f_PVf = (x_int_f - (w + a)) / w
	f_PVf = MIN(1.0, MAX(0.0, f_PVf))
  END IF

END SUBROUTINE compute_f_PVf


!===============================================================================
! SUBROUTINE: compute_S_ground
!===============================================================================
SUBROUTINE compute_S_ground(S_down, L, a, beta, h_PV, lambda, phi, S_nat)

  IMPLICIT NONE

  ! Input variables
  REAL, INTENT(IN) :: S_down   ! Direct solar radiation on horizontal surface [W/mÂ²]
  REAL, INTENT(IN) :: L, a, beta, h_PV, lambda, phi

  ! Output variable
  REAL, INTENT(OUT) :: S_nat   ! Direct solar radiation on ground [W/mÂ²]

  ! Local variable
  REAL :: f_nat                ! Ground visibility factor [-]

  CALL compute_f_ground(L, a, beta, h_PV, lambda, phi, f_nat)
  !PRINT *, ' f_nat    = ', f_nat

  S_nat = f_nat * S_down

END SUBROUTINE compute_S_ground


!===============================================================================
! SUBROUTINE: compute_S_PVb
!===============================================================================
SUBROUTINE compute_S_PVb(S_down, L, a, beta, h_PV, lambda, phi, S_PVb)

  IMPLICIT NONE

  ! Input variables
  REAL, INTENT(IN) :: S_down, L, a, beta, h_PV, lambda, phi

  ! Output variable
  REAL, INTENT(OUT) :: S_PVb   ! Direct solar radiation on PV back side [W/mÂ²]

  ! Local variables
  REAL :: f_PVb           ! Back-side visibility factor [-]
  REAL :: n_dot_s         ! Dot product between back normal and sun direction [-]
  REAL :: beam_irradiance ! Beam irradiance normal to sun rays [W/mÂ²]

  CALL compute_f_PVb(L, a, beta, h_PV, lambda, phi, f_PVb)
  !PRINT *, ' f_PVb    = ', f_PVb

  n_dot_s = -SIN(beta) * SIN(lambda) * SIN(phi) - COS(beta) * COS(lambda)

  beam_irradiance = S_down / COS(lambda)

  S_PVb = MAX(0.0, beam_irradiance * f_PVb * n_dot_s)

END SUBROUTINE compute_S_PVb


!===============================================================================
! SUBROUTINE: compute_S_PVf
!===============================================================================
SUBROUTINE compute_S_PVf(S_down, L, a, beta, h_PV, lambda, phi, S_PVf)

  IMPLICIT NONE

  ! Input variables
  REAL, INTENT(IN) :: S_down, L, a, beta, h_PV, lambda, phi

  ! Output variable
  REAL, INTENT(OUT) :: S_PVf   ! Direct solar radiation on PV front side [W/mÂ²]

  ! Local variables
  REAL :: f_PVf           ! Front-side visibility factor [-]
  REAL :: n_dot_s         ! Dot product between front normal and sun direction [-]
  REAL :: beam_irradiance ! Beam irradiance normal to sun rays [W/mÂ²]

  CALL compute_f_PVf(L, a, beta, h_PV, lambda, phi, f_PVf)
  !PRINT *, ' f_PVf    = ', f_PVf

  n_dot_s = SIN(beta) * SIN(lambda) * SIN(phi) + COS(beta) * COS(lambda)

  beam_irradiance = S_down / COS(lambda)

  S_PVf = MAX(0.0, beam_irradiance * f_PVf * n_dot_s)

END SUBROUTINE compute_S_PVf

