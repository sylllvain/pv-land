module solar_position
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: PI = 3.14159265358979323846_dp
  real(dp), parameter :: LAT = 44.9_dp * PI / 180.0_dp  ! radians (Salaunes)
  real(dp), parameter :: I_SC = 1367.0_dp                ! W/m^2

contains

  subroutine solar_zenith(doy, hour_in, theta, cos_theta)
    integer, intent(in)  :: doy
    real(dp), intent(in) :: hour_in
    real(dp), intent(out):: theta, cos_theta
    real(dp) :: hour, decl, ha

    hour = hour_in - 0.25_dp

    decl = 23.45_dp * PI / 180.0_dp * &
           sin(2.0_dp * PI * real(284 + doy, dp) / 365.0_dp)

    ha = (hour - 12.0_dp) * 15.0_dp * PI / 180.0_dp

    cos_theta = sin(LAT) * sin(decl) + cos(LAT) * cos(decl) * cos(ha)

    if (cos_theta < 0.0_dp) cos_theta = 0.0_dp

    theta = acos(cos_theta)

  end subroutine solar_zenith

  function solar_azimuth(doy, hour_in, lat_rad) result(az_deg)
    integer, intent(in)  :: doy
    real(dp), intent(in) :: hour_in, lat_rad
    real(dp) :: az_deg
    real(dp) :: hour, decl, ha, sin_az, cos_az, az

    hour = hour_in - 0.25_dp

    decl = 23.45_dp * PI / 180.0_dp * &
           sin(2.0_dp * PI * real(284 + doy, dp) / 365.0_dp)

    ha = (hour - 12.0_dp) * 15.0_dp * PI / 180.0_dp

    sin_az = sin(ha)
    cos_az = cos(ha) * sin(lat_rad) - tan(decl) * cos(lat_rad)

    az = atan2(sin_az, cos_az)

    az_deg = mod(az * 180.0_dp / PI + 180.0_dp, 360.0_dp)
    if (az_deg < 0.0_dp) az_deg = az_deg + 360.0_dp

  end function solar_azimuth

end module solar_position


!==================================================================
program solar_day_variation
  use solar_position
  implicit none

  integer  :: doy_input          ! jour de l'annee choisi
  integer  :: n_steps, i
  real(dp) :: dt_min, hour_frac
  real(dp) :: theta, cos_theta, az_deg
  integer  :: h, m
  integer  :: unit_out

  ! -----------------------------
  ! Parametres utilisateur
  ! -----------------------------
  doy_input = 185        ! ex: 185 = 04 juillet (solstice d'ete)
  dt_min    = 30.0_dp    ! pas de temps en minutes

  n_steps = int(24.0_dp * 60.0_dp / dt_min)  ! nombre de pas sur 24h

  ! -----------------------------
  ! Ouverture du fichier de sortie
  ! -----------------------------
  open(newunit=unit_out, file='solar_angles_day.txt', status='replace', action='write')
  write(unit_out, '(A)') 'heure  zenith_deg  azimuth_deg  cos_theta'

  ! -----------------------------
  ! Boucle temporelle sur la journee
  ! -----------------------------
  do i = 0, n_steps - 1

    hour_frac = real(i, dp) * dt_min / 60.0_dp   ! heure decimale (0 a 23.5)

    call solar_zenith(doy_input, hour_frac, theta, cos_theta)
    az_deg = solar_azimuth(doy_input, hour_frac, LAT)

    h = int(hour_frac)
    m = int((hour_frac - real(h, dp)) * 60.0_dp)

    write(unit_out, '(I2.2,":",I2.2,3X,F8.3,3X,F8.3,3X,F8.4)') &
         h, m, theta * 180.0_dp / PI, az_deg, cos_theta

  end do

  close(unit_out)

  print *, 'Calcul termine. Resultats ecrits dans solar_angles_day.txt'
  print *, 'Jour de l annee (doy) = ', doy_input
  print *, 'Nombre de pas de temps = ', n_steps

end program solar_day_variation
