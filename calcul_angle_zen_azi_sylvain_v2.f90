module solar_park_angles
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

  !----------------------------------------------------------
  ! Angle zenithal geographique theta (rad) + declinaison + angle horaire
  !----------------------------------------------------------
  subroutine solar_zenith(doy, hour, lat_rad, theta, cos_theta, decl, ha)
    integer,  intent(in)  :: doy
    real(dp), intent(in)  :: hour, lat_rad
    real(dp), intent(out) :: theta, cos_theta, decl, ha

    decl = (23.45_dp * PI / 180.0_dp) * &
           sin(2.0_dp * PI * real(284 + doy, dp) / 365.0_dp)

    ha = (hour - 12.0_dp) * 15.0_dp * PI / 180.0_dp

    cos_theta = sin(lat_rad) * sin(decl) + cos(lat_rad) * cos(decl) * cos(ha)

    if (cos_theta < 0.0_dp) cos_theta = 0.0_dp

    theta = acos(cos_theta)

  end subroutine solar_zenith

  !----------------------------------------------------------
  ! Azimut solaire geographique gamma (deg), Nord = 0, sens horaire
  !----------------------------------------------------------
  function solar_azimuth(decl, ha, lat_rad) result(az_deg)
    real(dp), intent(in) :: decl, ha, lat_rad
    real(dp) :: az_deg
    real(dp) :: sin_az, cos_az, az

    sin_az = sin(ha)
    cos_az = cos(ha) * sin(lat_rad) - tan(decl) * cos(lat_rad)

    az = atan2(sin_az, cos_az)

    az_deg = mod(az * 180.0_dp / PI + 180.0_dp, 360.0_dp)
    if (az_deg < 0.0_dp) az_deg = az_deg + 360.0_dp

  end function solar_azimuth

  !----------------------------------------------------------
  ! Projection dans le repere du parc : lambda (signe) et phi (deg)
  !----------------------------------------------------------
  subroutine compute_lambda_phi(zenith_deg, azimuth_deg, az_park_deg, &
                                 lambda_deg, phi_deg)
    real(dp), intent(in)  :: zenith_deg, azimuth_deg, az_park_deg
    real(dp), intent(out) :: lambda_deg, phi_deg

    real(dp) :: theta, gamma, alpha
    real(dp) :: sN, sE, sZ
    real(dp) :: s_x, s_y, s_z
    real(dp) :: h, sign_x, ratio

    theta = zenith_deg  * PI / 180.0_dp
    gamma = azimuth_deg * PI / 180.0_dp
    alpha = az_park_deg * PI / 180.0_dp

    ! Etape 5 : vecteur solaire dans le repere geographique
    sN = sin(theta) * cos(gamma)
    sE = sin(theta) * sin(gamma)
    sZ = cos(theta)

    ! Etape 6 : rotation vers le repere du parc
    s_x = cos(alpha) * sN + sin(alpha) * sE
    s_y = -sin(alpha) * sN + cos(alpha) * sE
    s_z = sZ

    ! Etape 7 : lambda = angle zenithal signe (repere parc)
    sign_x = 1.0_dp
    if (s_x < 0.0_dp) sign_x = -1.0_dp
    if (s_x == 0.0_dp) sign_x = 0.0_dp

    lambda_deg = acos(max(-1.0_dp, min(1.0_dp, s_z))) * sign_x * 180.0_dp / PI

    ! Etape 8 : phi = azimut relatif a la direction des rangees
    h = sqrt(s_x**2 + s_y**2)
    if (h > 0.0_dp) then
      ratio = s_y / h
    else
      ratio = s_y / 1.0_dp
    end if
    ratio = max(-1.0_dp, min(1.0_dp, ratio))

    phi_deg = acos(ratio) * 180.0_dp / PI

  end subroutine compute_lambda_phi

end module solar_park_angles


!==================================================================
program angles_solaires_parc
  use solar_park_angles
  implicit none

  ! -----------------------------
  ! Constantes site / parc
  ! -----------------------------
  real(dp), parameter :: LAT_DEG  = 44.9_dp   ! latitude du site [deg]
  real(dp), parameter :: LAT      = LAT_DEG * PI / 180.0_dp
  real(dp), parameter :: AZ_PARK  = 177.0_dp  ! azimut meteo avant du parc [deg]

  integer,  parameter :: DOY_INPUT = 185        ! ex: 185 = 04 juillet (solstice d'ete)
  real(dp), parameter :: DT_MIN    = 30.0_dp

  integer  :: i, h, m, unit_out, n_steps
  real(dp) :: hour, theta, cos_theta, decl, ha
  real(dp) :: zenith_deg, azimuth_deg, lambda_deg, phi_deg

  ! -----------------------------
  ! Ouverture du fichier de sortie
  ! -----------------------------
  open(newunit=unit_out, file='angles_parc_jour.csv', status='replace', action='write')
  write(unit_out, '(A)') 'heure,zenith_deg,azimuth_deg,lambda_deg,phi_deg'

  n_steps = int(24.0_dp * 60.0_dp / dt_min)  ! nombre de pas sur 24h

  ! -----------------------------
  ! Boucle temporelle sur la journee
  ! -----------------------------
  do i = 0, N_STEPS - 1

    hour = real(i, dp) * DT_MIN / 60.0_dp ! heure decimale (0 a 23.5)

    ! Repere geographique
    call solar_zenith(DOY_INPUT, hour, LAT, theta, cos_theta, decl, ha)
    zenith_deg  = theta * 180.0_dp / PI
    azimuth_deg = solar_azimuth(decl, ha, LAT)

    ! Repere du parc
    call compute_lambda_phi(zenith_deg, azimuth_deg, AZ_PARK, lambda_deg, phi_deg)

    h = int(hour)
    m = int((hour - real(h, dp)) * 60.0_dp)

    write(unit_out, '(I2.2,":",I2.2,4(",",F9.4))') &
         h, m, zenith_deg, azimuth_deg, lambda_deg, phi_deg

  end do

  close(unit_out)

  print *, 'Calcul termine. Resultats ecrits dans angles_parc_jour.csv'
  print *, 'Jour de l annee (doy) = ', DOY_INPUT
  print *, 'Nombre de pas de temps = ', n_steps

end program angles_solaires_parc
