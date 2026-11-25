!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE pvland_out                        #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE module_pvland_out

!------------------------------------------------------------------------
!
! AUTHOR: ???
!
! MODIFICATIONS:
!
!
!------------------------------------------------------------------------

!   USE module_pvland_out

!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
  REAL, ALLOCATABLE :: T_can(:,:), T_PVa(:,:)
  REAL, ALLOCATABLE :: T_PVf(:,:), T_cell(:,:), T_PVb(:,:), ptsfc_eff(:,:), qvsfc_eff(:,:)
  REAL, ALLOCATABLE :: U_can(:,:), U_PVa(:,:)
  ! PV-LAND state of version 2, saved in the ARPS restart files
  REAL, ALLOCATABLE :: q_can(:,:)     ! Specific humidity of the canopy air (kg/kg); < 0: not set yet
  REAL, ALLOCATABLE :: L_obk(:,:)     ! Obukhov length (m), 0 = neutral

  REAL, ALLOCATABLE   :: H_ref(:,:),Rn_ref(:,:)
  REAL, ALLOCATABLE   :: H_nat(:,:), LE_nat(:,:), G_nat(:,:), Rn_nat(:,:)
  REAL, ALLOCATABLE   :: SW_out_sky(:,:), LW_out_sky(:,:)
  
  REAL, ALLOCATABLE   :: H_PVa_ref(:,:), H_can_ref(:,:)
  REAL, ALLOCATABLE   :: H_PVb(:,:), H_PVf(:,:)
  REAL, ALLOCATABLE   :: G_SW_nat(:,:), J_SW_nat(:,:)
  REAL, ALLOCATABLE   :: G_LW_nat(:,:), J_LW_nat(:,:)

  REAL, ALLOCATABLE   :: qLW_PVf(:,:), qLW_PVb(:,:)
  REAL, ALLOCATABLE   :: qcond_PVf(:,:), qcond_PVb(:,:),&
	qSW_PVf(:,:), qSW_PVb(:,:), qPV(:,:)
  REAL, ALLOCATABLE   :: JLW_PVf(:,:), JLW_PVb(:,:)


END MODULE module_pvland_out

!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE pvland_in                         #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE module_pvland_in

!------------------------------------------------------------------------
!
! AUTHOR: ???
!
! MODIFICATIONS:
!
!
!------------------------------------------------------------------------

!   USE module_pvland_in

!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
  REAL, ALLOCATABLE, DIMENSION (:,:) :: pv_park
  REAL, ALLOCATABLE, DIMENSION (:,:) :: a, beta, L, width, z_PVa,&
  h_PV, h_clear, lambda_f_park, l_nat_can, theta_0,&
  albedo_nat, albedo_PVf, albedo_PVb, &
                             emiss_nat, emiss_PVf, emiss_PVb, z_can

  REAL, ALLOCATABLE, DIMENSION (:,:)   :: delta_PVf, delta_cell, delta_PVb, &
                        rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                        k_glass, k_EVA, k_backsheet, &
                        delta_glass, delta_EVA_top, &
                        delta_EVA_bot, delta_backsheet, &
                        eta_STC, gamma_PV

  REAL, ALLOCATABLE   :: z0_park(:,:), zd_park(:,:)
  REAL, ALLOCATABLE   :: AZ_PARK(:,:)

END MODULE module_pvland_in

!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE rad_split                         #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

module rad_split
   implicit none
   private
   public :: split_sw

   real, parameter :: PI = 3.14159265358979323846

contains

   integer function day_of_year(year, month, day)
      implicit none
      integer, intent(in) :: year, month, day
      integer :: i
      integer, dimension(12) :: mdays
      logical :: leap

      mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      leap = (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0))
      if (leap) mdays(2) = 29

      day_of_year = day
      do i = 1, month-1
         day_of_year = day_of_year + mdays(i)
      end do

   end function day_of_year


   subroutine split_sw(SWin, year, month, day, hour, minute, &
                       lat_deg, lon_deg, tz_offset, &
                       SW_diffus, SW_direct, I_SC)

      implicit none

      real, intent(in) :: SWin
      integer, intent(in) :: year, month, day
      integer, intent(in) :: hour, minute
      real, intent(in) :: lat_deg, lon_deg
      integer, intent(in) :: tz_offset
      real, intent(out) :: SW_diffus, SW_direct
      real, intent(in), optional :: I_SC

      real :: solar_const
      integer :: doy

      real :: lat
      real :: B, EoT, LSTM, TC
      real :: hour_legal, hour_solar
      real :: decl, ha
      real :: cos_theta
      real :: E0, G0
      real :: kt
      real :: fd

      if (present(I_SC)) then
         solar_const = I_SC
      else
         solar_const = 1367.
      end if

      if (SWin <= 0.) then
         SW_diffus = 0.
         SW_direct = 0.
         return
      end if

      lat = lat_deg * PI / 180.
      doy = day_of_year(year, month, day)

      ! True solar time
      B = 2. * PI * (doy - 81.) / 364.

      EoT = 9.87 * sin(2.*B) - &
            7.53 * cos(B) - &
            1.5  * sin(B)

      LSTM = 15. * tz_offset
      TC = 4. * (lon_deg - LSTM) + EoT

      hour_legal = real(hour,8) + real(minute,8)/60.
      hour_solar = hour_legal + TC/60.

      ! Solar geometry
      decl = (23.45 * PI/180.) * &
             sin(2.*PI*(284.+doy)/365.)

      ha = ((hour_solar-12.)*15.) * PI/180.

      cos_theta = sin(lat)*sin(decl) + &
                  cos(lat)*cos(decl)*cos(ha)

      ! Sun very low: all diffuse
      if (cos_theta < 0.15) then
         SW_diffus = SWin
         SW_direct = 0.
         return
      end if

      ! Extraterrestrial radiation
      E0 = 1. + 0.033*cos(2.*PI*doy/365.)
      G0 = solar_const * E0 * cos_theta

      ! Clearness index
      kt = SWin / G0

      kt = max(0., min(kt,1.))

      ! Erbs diffuse fraction
      if (kt <= 0.22) then

         fd = 1. - 0.09*kt

      else if (kt <= 0.80) then

         fd = 0.9511                    &
            - 0.1604*kt                 &
            + 4.388*kt**2               &
            -16.638*kt**3               &
            +12.336*kt**4

      else

         fd = 0.165

      end if

      fd = max(0., min(fd,1.))

      SW_diffus = fd * SWin
      SW_direct = SWin - SW_diffus

   end subroutine split_sw

end module rad_split

!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE rad_pvland                        #########
!#########                                                      #########
!#########                      Adapted by                      #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

module rad_pvland
   implicit none

  REAL, ALLOCATABLE :: szenith_deg(:,:), sazimuth_deg(:,:)

end module rad_pvland

!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE pvland_rad_cache                  #########
!#########                                                      #########
!########################################################################
!########################################################################
!
! Quantities of the 8-surface radiosity model that do not depend on the
! surface temperatures, computed once (INIT_PVLAND, or every time step
! when vf_update_opt = 1, see module pvland_nml) or once per change of the solar forcing
! (compute_all_forcings), instead of at every call.
! [SYLVAIN]

module pvland_rad_cache
   implicit none

  INTEGER, PARAMETER :: nsurf_pv = 8      ! number of radiating surfaces

  ! vf_update_opt (when to compute the view factors): see module pvland_nml

  ! Geometry only (computed in INIT_PVLAND)
  REAL,    ALLOCATABLE :: F_pv(:,:,:,:)     ! view factors   (8,8,nx,ny)
  REAL,    ALLOCATABLE :: Fsky_pv(:,:,:)    ! sky view factors (8,nx,ny)
  ! LU factorisations of the radiosity matrices (computed in INIT_PVLAND)
  !   SW: I - diag(albedo) F        LW: I - diag(1-emissivity) F
  REAL,    ALLOCATABLE :: LU_SW_pv(:,:,:,:), LU_LW_pv(:,:,:,:)   ! (8,8,nx,ny)
  INTEGER, ALLOCATABLE :: ipiv_SW_pv(:,:,:), ipiv_LW_pv(:,:,:)   ! (8,nx,ny)

  ! Shortwave solution, recomputed only when its inputs change
  LOGICAL, ALLOCATABLE :: sw_valid(:,:)                   ! (nx,ny)
  REAL,    ALLOCATABLE :: sw_key(:,:,:)                   ! (4,nx,ny): swin, S_diffuse, lambda, phi
  REAL,    ALLOCATABLE :: J_SW_pv(:,:,:), G_SW_pv(:,:,:), q_SW_pv(:,:,:)   ! (8,nx,ny)

end module pvland_rad_cache

!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE pvland_nml                        #########
!#########                                                      #########
!########################################################################
!########################################################################
!
! Input parameters of PV-LAND, read from the namelist &pvland:
!   ARPS:    in arps.input (initpara, only when sfcphy = 5), after &soil_ebm
!   offline: in the file pvland.nml (pvland_offline_main)
! The values below are the defaults, used for any parameter missing from
! the namelist. INIT_PVLAND copies them into the (nx,ny) arrays of
! module_pvland_in.  [SYLVAIN]

MODULE pvland_nml

  IMPLICIT NONE

  ! ---- Options ----
  INTEGER :: vf_update_opt = 0      ! View factors, sky view factors, radiosity LU factors and tilt-derived geometry: 0 = computed once at initialisation (fixed panels), 1 = recomputed at every PV-LAND time step (moving panels)

  ! ---- Location of the PV park (offline: the single column is in the park or not) ----
  INTEGER :: pv_park_opt = 0        ! 0 = PV park over the whole domain, 1 = rectangle pv_park_i1:pv_park_i2 x pv_park_j1:pv_park_j2, 2 = no PV park (natural ground only, with the PV-LAND soil physics)
  INTEGER :: pv_park_i1 = 1         ! First global x index (scalar point) of the PV park (pv_park_opt = 1)
  INTEGER :: pv_park_i2 = 99999     ! Last global x index of the PV park, clipped to the domain (pv_park_opt = 1)
  INTEGER :: pv_park_j1 = 1         ! First global y index (scalar point) of the PV park (pv_park_opt = 1)
  INTEGER :: pv_park_j2 = 99999     ! Last global y index of the PV park, clipped to the domain (pv_park_opt = 1)

  ! ---- Park geometry ----
  REAL :: pv_L = 3.23               ! Panel length along the slope (m)
  REAL :: pv_a = 3.4                ! Inter-row spacing (m)
  REAL :: pv_beta_deg = 24.0        ! Panel tilt angle from the horizontal (deg)
  REAL :: pv_h_PV = 1.9             ! Height of the panel top edge above the ground (m)
  REAL :: pv_theta_0_deg = -3.0     ! Park orientation angle used by the wind profile (deg)
  REAL :: pv_az_park = 177.0        ! Meteorological azimuth faced by the panel fronts (deg)

  ! ---- Natural vegetation of the park (wind profile, version 2) ----
  REAL :: pv_zd_nat = 0.285         ! Displacement height of the natural vegetation alone (m); its roughness z0 is roufns
  REAL :: pv_ztz0_nat = 0.1         ! Ratio zt/z0 of thermal to momentum roughness, natural ground -> canopy air (-)
  REAL :: pv_ztz0_air = 1.0         ! Ratio zt/z0, canopy air and air above the panels -> reference level (-)

  ! ---- Radiative properties ----
  REAL :: pv_albedo_nat = 0.22      ! Albedo of the natural ground surface (-)
  REAL :: pv_albedo_PVf = 0.11      ! Albedo of the panel front (-)
  REAL :: pv_albedo_PVb = 0.11      ! Albedo of the panel back (-)
  REAL :: pv_emiss_nat = 0.98       ! Emissivity of the natural ground surface (-)
  REAL :: pv_emiss_PVf = 0.81       ! Emissivity of the panel front (-)
  REAL :: pv_emiss_PVb = 0.845      ! Emissivity of the panel back (-)

  ! ---- Convection on the panel faces ----
  REAL :: pv_C_for_f = 6.80         ! Forced convection, front face: h = C_for_f * U_PVa**p_for_f (W/m2/K)
  REAL :: pv_C_for_b = 6.80         ! Forced convection, back face: h = C_for_b * U_can**p_for_b (W/m2/K)
  REAL :: pv_p_for_f = 0.80         ! Wind exponent of the forced convection, front face (-)
  REAL :: pv_p_for_b = 0.80         ! Wind exponent of the forced convection, back face (-)
  REAL :: pv_fac_h_free = 2.50      ! Scaling of the free convection of Kaplani & Kaplanis (2014) (-)

  ! ---- Obukhov length (computed by PV-LAND from its own fluxes) ----
  REAL :: pv_tau_L = 300.0          ! Relaxation time of 1/L towards the instantaneous value (s), 0 = no relaxation
  REAL :: pv_L_floor = 0.5          ! Minimum |L| (m)

  ! ---- PV module, 3-node thermal model ----
  REAL :: pv_delta_PVf = 3.8e-3     ! Equivalent thickness of the front node (m)
  REAL :: pv_delta_cell = 0.225e-3  ! Equivalent thickness of the cell node (m)
  REAL :: pv_delta_PVb = 1.1e-3     ! Equivalent thickness of the back node (m)
  REAL :: pv_rho_cp_PVf = 1.580e6   ! Volumetric heat capacity of the front node (J/m3/K)
  REAL :: pv_rho_cp_cell = 1.577e6  ! Volumetric heat capacity of the cell node (J/m3/K)
  REAL :: pv_rho_cp_PVb = 1.776e6   ! Volumetric heat capacity of the back node (J/m3/K)
  REAL :: pv_k_glass = 1.8          ! Thermal conductivity of the glass (W/m/K)
  REAL :: pv_k_EVA = 0.35           ! Thermal conductivity of the EVA (W/m/K)
  REAL :: pv_k_backsheet = 0.20     ! Thermal conductivity of the backsheet (W/m/K)
  REAL :: pv_delta_glass = 3.2e-3   ! Thickness of the glass sub-layer (m)
  REAL :: pv_delta_EVA_top = 0.6e-3 ! Thickness of the top EVA sub-layer (m)
  REAL :: pv_delta_EVA_bot = 0.6e-3 ! Thickness of the bottom EVA sub-layer (m)
  REAL :: pv_delta_backsheet = 0.5e-3 ! Thickness of the backsheet sub-layer (m)
  REAL :: pv_eta_STC = 0.18         ! PV efficiency at standard test conditions (25 C, 1000 W/m2) (-)
  REAL :: pv_gamma_PV = 0.0045      ! Temperature coefficient of the PV efficiency (1/K)

  ! ---- Soil of the PV-LAND points (soilebm_frc_pv/evapflx_pv; a negative value keeps the ARPS table, soilcst.inc) ----
  REAL :: pv_cgv = 3.0e-5           ! Thermal coefficient of the vegetation cgv (K m2/J); ARPS: 1.5e-5
  REAL :: pv_wsat = 0.315           ! Saturated volumetric moisture content (m3/m3), all soil types
  REAL :: pv_wfc = 0.135            ! Field capacity (m3/m3), all soil types
  REAL :: pv_wwlt = 0.05            ! Wilting point (m3/m3), all soil types
  REAL :: pv_rsmin = 70.0           ! Minimum stomatal resistance (s/m), all vegetation types
  REAL :: pv_wref = 0.219           ! Moisture threshold of the transpiration stress and of the bare-soil evaporation (m3/m3)
  INTEGER :: pv_drain_opt = 1       ! Drainage of the deep reservoir: 0 = off, 1 = on (c3 * max(w2 - wfc, 0) / tau)
  REAL :: pv_c3_drain = 0.1         ! Drainage coefficient c3 (-)

  NAMELIST /pvland/ vf_update_opt,                                       &
           pv_park_opt, pv_park_i1, pv_park_i2, pv_park_j1, pv_park_j2,  &
           pv_L, pv_a, pv_beta_deg, pv_h_PV, pv_theta_0_deg, pv_az_park, &
           pv_zd_nat, pv_ztz0_nat, pv_ztz0_air,                          &
           pv_albedo_nat, pv_albedo_PVf, pv_albedo_PVb,                  &
           pv_emiss_nat, pv_emiss_PVf, pv_emiss_PVb,                     &
           pv_C_for_f, pv_C_for_b, pv_p_for_f, pv_p_for_b, pv_fac_h_free, &
           pv_tau_L, pv_L_floor,                                         &
           pv_delta_PVf, pv_delta_cell, pv_delta_PVb,                    &
           pv_rho_cp_PVf, pv_rho_cp_cell, pv_rho_cp_PVb,                 &
           pv_k_glass, pv_k_EVA, pv_k_backsheet,                         &
           pv_delta_glass, pv_delta_EVA_top, pv_delta_EVA_bot,           &
           pv_delta_backsheet, pv_eta_STC, pv_gamma_PV,                  &
           pv_cgv, pv_wsat, pv_wfc, pv_wwlt, pv_rsmin, pv_wref,          &
           pv_drain_opt, pv_c3_drain

END MODULE pvland_nml
