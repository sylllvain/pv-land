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
  REAL :: lo, hi       ! Shadow of row n clipped to [0, period] [m]
  REAL :: cur_end      ! Right end of the part of [0, period] already merged [m]
  REAL :: shaded_len   ! Measure of S_total within [0, period] [m]
  INTEGER :: n            ! Row index for loop
  INTEGER, PARAMETER :: n_rows = 3 ! Neighbouring rows on each side (N_rows)

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
  ! Exact measure of S_total within [0, period] (union of intervals)
  !-------------------------
  ! S_n = [s_min, s_max] + n*period. The S_n have the same width and are
  ! shifted by one period, so they are already sorted by increasing n:
  ! the union is merged in a single pass, keeping the right end (cur_end)
  ! of what has been counted so far to avoid counting overlaps twice.
  shaded_len = 0.0
  cur_end    = 0.0

  DO n = -n_rows, n_rows
     lo = MAX(s_min + REAL(n) * period, 0.0)
     hi = MIN(s_max + REAL(n) * period, period)
     lo = MAX(lo, cur_end)
     IF (hi > lo) THEN
        shaded_len = shaded_len + (hi - lo)
        cur_end    = hi
     END IF
  END DO

  !-------------------------
  ! Visibility factor
  !-------------------------
  ! f_nat = 1 - |S_total ∩ [0, period]| / period
  f_nat = 1.0 - shaded_len / period

END SUBROUTINE compute_f_ground


!===============================================================================
! SUBROUTINE: compute_f_PVb
!
! Purpose: Compute the visibility factor (sunlit fraction) of the back side
!          of a PV panel, accounting for the shadow of the neighbouring row.
!          Same inputs as compute_f_ground.
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

  ! check if the panel is vertical (threshold 1.5d0 near pi/2)
  ELSE IF (beta >= 1.5) THEN
	! in this case, the same visibility factor is set for both faces
	! and the dot product makes one of them to vanish
	f_PVb = ABS(a / (L * TAN(lambda) * SIN(phi)))
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
! Purpose: Compute the visibility factor (sunlit fraction) of the front side
!          of a PV panel, accounting for the shadow of the neighbouring row.
!          Same inputs as compute_f_ground.
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

  ! check if the panel is vertical (threshold 1.5d0 near pi/2)
  ELSE IF (beta >= 1.5) THEN
	! in this case, the same visibility factor is set for both faces
	! and the dot product makes one of them to vanish
	f_PVf = ABS(a / (L * TAN(lambda) * SIN(phi)))
	f_PVf = MIN(1.0, MAX(0.0, f_PVf))
  ELSE
	!-------------------------
	! Intersection calculation
	!-------------------------
	numerator = 2.0 * (w + a) / (TAN(lambda) * SIN(phi)) + H * (w + a) / w
	denominator = 1.0 / (TAN(lambda) * SIN(phi)) + H / w

	x_int_f = numerator / denominator

	!-------------------------
	! Visibility factor
	!-------------------------
	f_PVf = (x_int_f - (w + a)) / w
	f_PVf = MIN(1.0, MAX(0.0, f_PVf))
  END IF

END SUBROUTINE compute_f_PVf


!===============================================================================
! SUBROUTINE: compute_S_ground
!
! Purpose: Direct solar radiation received by the ground, f_nat * S_down
!===============================================================================
SUBROUTINE compute_S_ground(S_down, L, a, beta, h_PV, lambda, phi, S_nat)

  IMPLICIT NONE

  ! Input variables
  REAL, INTENT(IN) :: S_down   ! Direct solar radiation on horizontal surface [W/m2]
  REAL, INTENT(IN) :: L, a, beta, h_PV, lambda, phi

  ! Output variable
  REAL, INTENT(OUT) :: S_nat   ! Direct solar radiation on ground [W/m2]

  ! Local variable
  REAL :: f_nat                ! Ground visibility factor [-]

  CALL compute_f_ground(L, a, beta, h_PV, lambda, phi, f_nat)

  S_nat = f_nat * S_down

END SUBROUTINE compute_S_ground


!===============================================================================
! SUBROUTINE: compute_S_PVb
!
! Purpose: Direct solar radiation received by the PV back side:
!          beam irradiance (S_down / cos(lambda)) * f_PVb * cos(incidence)
!===============================================================================
SUBROUTINE compute_S_PVb(S_down, L, a, beta, h_PV, lambda, phi, S_PVb)

  IMPLICIT NONE

  ! Input variables
  REAL, INTENT(IN) :: S_down, L, a, beta, h_PV, lambda, phi

  ! Output variable
  REAL, INTENT(OUT) :: S_PVb   ! Direct solar radiation on PV back side [W/m2]

  ! Local variables
  REAL :: f_PVb           ! Back-side visibility factor [-]
  REAL :: n_dot_s         ! Dot product between back normal and sun direction [-]
  REAL :: beam_irradiance ! Beam irradiance normal to sun rays [W/m2]

  CALL compute_f_PVb(L, a, beta, h_PV, lambda, phi, f_PVb)

  n_dot_s = -SIN(beta) * SIN(lambda) * SIN(phi) - COS(beta) * COS(lambda)

  beam_irradiance = S_down / COS(lambda)

  S_PVb = MAX(0.0, beam_irradiance * f_PVb * n_dot_s)

END SUBROUTINE compute_S_PVb


!===============================================================================
! SUBROUTINE: compute_S_PVf
!
! Purpose: Direct solar radiation received by the PV front side:
!          beam irradiance (S_down / cos(lambda)) * f_PVf * cos(incidence)
!===============================================================================
SUBROUTINE compute_S_PVf(S_down, L, a, beta, h_PV, lambda, phi, S_PVf)

  IMPLICIT NONE

  ! Input variables
  REAL, INTENT(IN) :: S_down, L, a, beta, h_PV, lambda, phi

  ! Output variable
  REAL, INTENT(OUT) :: S_PVf   ! Direct solar radiation on PV front side [W/m2]

  ! Local variables
  REAL :: f_PVf           ! Front-side visibility factor [-]
  REAL :: n_dot_s         ! Dot product between front normal and sun direction [-]
  REAL :: beam_irradiance ! Beam irradiance normal to sun rays [W/m2]

  CALL compute_f_PVf(L, a, beta, h_PV, lambda, phi, f_PVf)

  n_dot_s = SIN(beta) * SIN(lambda) * SIN(phi) + COS(beta) * COS(lambda)

  beam_irradiance = S_down / COS(lambda)

  S_PVf = MAX(0.0, beam_irradiance * f_PVf * n_dot_s)

END SUBROUTINE compute_S_PVf

