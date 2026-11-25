!##################################################################
!##################################################################
!######                                                      ######
!######              PROGRAM PVLAND_OFFLINE_MAIN             ######
!######                                                      ######
!######   Offline driver of the PV-LAND surface model:       ######
!######   - read the measured forcing (CSV)                  ######
!######   - build a single-column ARPS-like atmospheric state ######
!######   - INIT_PVLAND, then one call of the ARPS PV-LAND    ######
!######     surface step (pvland_sfcphy) per dtsfc step      ######
!######   - write the CSV outputs                            ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  The PV-LAND physics (core/, including the surface step pvland_sfcphy
!  of pvland_driver.f90) and the ARPS routines (arps_lib/) are copies of
!  the ARPS sources (see scripts/sync_from_arps.py). Only this file and
!  offline/pvland_park_map.f90 are specific to the offline model.
!  pvland_sfcphy is called with meas_opt = 1: diffuse SW, sun angles in
!  the park frame and Obukhov length are taken from the measurements.
!
!  Usage:  ./pv_land.exe [forcing.csv [pvland.nml]]
!          default forcing file: recorded_data_adjusted_daily_summer_2025.csv
!          default PV-LAND parameter file: pvland.nml (namelist &pvland,
!          the same as in ARPS arps.input; pv_park_opt = 2: no PV panels)
!
!  Grid: nx = ny = 2 (ARPS loops run from 1 to nx-1, so a single
!  column (1,1) is computed), nz = 3 (level 2 = reference level).
!
!-----------------------------------------------------------------------

PROGRAM pvland_offline_main

  USE module_pvland_in,  ONLY : h_PV, albedo_nat, emiss_nat
  USE module_pvland_out
  USE pvland_nml,        ONLY : pvland      ! namelist &pvland (PV park parameters)
  USE pvland_sfcphy_work, ONLY : uref_diag => uref,                    &
                                 cdh_can_diag => cdh_can,              &
                                 cdh_PVa_diag => cdh_PVa

  IMPLICIT NONE

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'

  ! ====== Grid parameters ======
  INTEGER, PARAMETER :: nx = 2, ny = 2, nz = 3
  INTEGER, PARAMETER :: nzsoil = 2
  INTEGER, PARAMETER :: nstyps = 1
  INTEGER, PARAMETER :: ntot_max = 100000  ! Upper bound for array allocation
  INTEGER :: ntot                           ! Actual number of forcing rows

  REAL, PARAMETER :: pi = 3.14159265359
  REAL, PARAMETER :: dt_forcing = 1800.     ! Forcing interval (s)

  ! ====== Site parameters (natural surface, soil, initial state) ======
  ! The PV park geometry and optical properties are set in INIT_PVLAND
  ! (core/pvland.f90), as in ARPS.
  REAL, PARAMETER :: z0_nat     = 0.05      ! Natural surface roughness (m) -> roufns
  REAL, PARAMETER :: veg_site   = 0.95      ! Vegetation fraction
  INTEGER, PARAMETER :: vegtyp_site  = 3
  REAL, PARAMETER :: lai_site   = 2.
  INTEGER, PARAMETER :: soiltyp_site = 6
  REAL, PARAMETER :: tdeep_ini  = 20.95 + 273.15   ! Initial deep soil temperature (K)
  REAL, PARAMETER :: wetsfc_ini = 0.11      ! Initial surface soil moisture (m3/m3)
  REAL, PARAMETER :: wetdp_ini  = 0.10      ! Initial deep soil moisture (m3/m3)

  ! ====== File I/O ======
  INTEGER :: temp, ios, t, it, k
  REAL :: frac
  CHARACTER(LEN=500) :: line
  CHARACTER(LEN=256) :: forcing_file    ! forcing CSV (1st command-line argument)
  INTEGER :: nstep_frc                  ! number of PV-LAND steps (dtsfc) per forcing interval
  CHARACTER(LEN=256) :: nml_file        ! PV-LAND namelist file (2nd command-line argument)

  ! ====== CSV input arrays ======
  REAL, ALLOCATABLE :: uref_input(:), udir_input(:)
  REAL, ALLOCATABLE :: tref_input(:), qvref_input(:), rain_input(:)
  REAL, ALLOCATABLE :: swout_input(:), swin_input(:)
  REAL, ALLOCATABLE :: lwout_input(:), lwin_input(:)
  REAL, ALLOCATABLE :: pref_input(:), rhoref_input(:)
  REAL, ALLOCATABLE :: phi_input(:), lambda_input(:), S_diffuse_input(:)
  REAL, ALLOCATABLE :: H_meas(:), LE_meas(:), Rn_meas(:), G_meas(:)
  REAL, ALLOCATABLE :: Lobk_input(:)    ! Obukhov length L [m]

  ! ====== ARPS-like atmospheric state (single column) ======
  REAL :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
  REAL :: wdspsfc(nx,ny)                 ! surface wind speed (used by pvland_sfcphy outside the park only)
  REAL :: ptprt(nx,ny,nz), pprt(nx,ny,nz), qv(nx,ny,nz)
  REAL :: ptbar(nx,ny,nz), pbar(nx,ny,nz), qvbar(nx,ny,nz)
  REAL :: zp(nx,ny,nz)
  REAL :: prcrate(nx,ny)
  REAL :: radsw(nx,ny), rnflx(nx,ny), radlwin(nx,ny)
  REAL :: psfc(nx,ny), rhosfc(nx,ny)
  REAL :: uref, udir, tref, qvref        ! interpolated forcing (scalars)

  ! ====== Surface properties and soil state (ARPS arrays) ======
  INTEGER :: soiltyp(nx,ny,nstyps), vegtyp(nx,ny)
  REAL :: stypfrct(nx,ny,nstyps)
  REAL :: lai(nx,ny), roufns(nx,ny), veg(nx,ny)
  REAL :: tsoil(nx,ny,nzsoil,0:nstyps), qsoil(nx,ny,nzsoil,0:nstyps)
  REAL :: wetcanp(nx,ny,0:nstyps), qvsfc(nx,ny,0:nstyps)
  REAL :: snowdpth(nx,ny)

  ! ====== Surface fluxes and drag coefficients (ARPS arrays) ======
  REAL :: cdm(nx,ny), cdh(nx,ny), cdq(nx,ny)
  REAL :: usflx(nx,ny), vsflx(nx,ny), ptsflx(nx,ny), qvsflx(nx,ny)
  REAL :: tem1(nx,ny)
  REAL :: pt1(nx,ny), tair(nx,ny), qvair(nx,ny)   ! set by pvland_sfcphy

  ! ====== Measured forcing replacing ARPS diagnostics ======
  REAL :: S_diffuse_in(nx,ny), lambda_in(nx,ny), phi_in(nx,ny), L_obk_in(nx,ny)

  !=====================================================================
  ! INITIALIZATION
  !=====================================================================

  WRITE(*,*) '========================================='
  WRITE(*,*) 'PV-LAND offline model (ARPS PV-LAND physics)'
  WRITE(*,*) '========================================='
  WRITE(*,*) ''

  ! ARPS control parameters used by the PV-LAND code
  dtsfc       = 2.
  nstep_frc   = NINT(dt_forcing/dtsfc)   ! PV-LAND steps per forcing interval
  nsfcst      = 1                         ! soil sub-steps per pvland_sfcphy call (ARPS: dtbig/dtsfc)
  moist       = 1
  initopt     = 1
  soilinitopt = 0
  nstyp       = nstyps
  nstep       = 1

  ! Surface properties
  veg(:,:)        = veg_site
  vegtyp(:,:)     = vegtyp_site
  lai(:,:)        = lai_site
  soiltyp(:,:,:)  = soiltyp_site
  stypfrct(:,:,:) = 1.0
  roufns(:,:)     = z0_nat
  snowdpth(:,:)   = 0.

  !=====================================================================
  ! READ CSV INPUT
  !=====================================================================

  forcing_file = 'recorded_data_adjusted_daily_summer_2025.csv'
  IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, forcing_file)

  WRITE(*,*) 'Reading forcing data from CSV: ', TRIM(forcing_file)

  ALLOCATE(uref_input(ntot_max), udir_input(ntot_max), tref_input(ntot_max), &
           qvref_input(ntot_max), rain_input(ntot_max), swout_input(ntot_max), &
           swin_input(ntot_max), lwout_input(ntot_max), lwin_input(ntot_max), &
           pref_input(ntot_max), rhoref_input(ntot_max), phi_input(ntot_max), &
           lambda_input(ntot_max), S_diffuse_input(ntot_max), H_meas(ntot_max), &
           LE_meas(ntot_max), Rn_meas(ntot_max), G_meas(ntot_max), Lobk_input(ntot_max))

  OPEN(UNIT=10, FILE=TRIM(forcing_file), STATUS='old')
  READ(10,*)  ! Skip header

  temp = 1
  DO
    READ(10,'(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    READ(line,*,IOSTAT=ios) t, uref_input(temp), udir_input(temp), &
                            tref_input(temp), qvref_input(temp), rain_input(temp), &
                            swout_input(temp), swin_input(temp), &
                            lwout_input(temp), lwin_input(temp), &
                            pref_input(temp), rhoref_input(temp), &
                            H_meas(temp), LE_meas(temp), G_meas(temp), Rn_meas(temp), &
                            S_diffuse_input(temp), lambda_input(temp), phi_input(temp), &
                            Lobk_input(temp)
    tref_input(temp) = tref_input(temp) + 273.15
    IF (ios /= 0) EXIT
    temp = temp + 1
    IF (temp > ntot_max) EXIT
  END DO
  ntot = temp - 1
  CLOSE(10)

  WRITE(*,*) 'CSV read: ', ntot, ' timesteps'
  WRITE(*,*) ''

  !=====================================================================
  ! INITIAL STATE (as ARPS: soil state, then INIT_PVLAND)
  !=====================================================================

  CALL set_atmosphere(1, 1, 0.)

  tsoil(:,:,1,:) = tref_input(1)       ! T_nat
  tsoil(:,:,2,:) = tdeep_ini           ! T_soil
  qsoil(:,:,1,:) = wetsfc_ini
  qsoil(:,:,2,:) = wetdp_ini
  wetcanp(:,:,:) = 0.
  qvsfc(:,:,:)   = qvref_input(1)

  ! Arrays allocated by arps.f90 in ARPS
  ALLOCATE(T_PVf(nx,ny), T_cell(nx,ny), T_PVb(nx,ny), ptsfc_eff(nx,ny), qvsfc_eff(nx,ny))

  ! PV-LAND parameters: namelist &pvland (same as in ARPS arps.input)
  nml_file = 'pvland.nml'
  IF (COMMAND_ARGUMENT_COUNT() >= 2) CALL GET_COMMAND_ARGUMENT(2, nml_file)
  OPEN(UNIT=11, FILE=TRIM(nml_file), STATUS='old', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'ERROR: PV-LAND parameter file not found: ', TRIM(nml_file)
    STOP 1
  END IF
  READ(11, NML=pvland, IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'ERROR reading namelist &pvland in ', TRIM(nml_file), ' (iostat =', ios, ')'
    STOP 1
  END IF
  CLOSE(11)
  WRITE(*,*) 'PV-LAND parameters read from ', TRIM(nml_file), ':'
  WRITE(*,NML=pvland)

  CALL init_pvland(nx,ny,nz,ptbar,ptprt,pbar,pprt,qvsfc(:,:,1),tsoil(:,:,1,1),tem1)

  ! INIT_PVLAND fills the computed column (1,1) only: copy it to the
  ! (unused) border points so that no array holds undefined values.
  T_can(:,:)     = T_can(1,1)
  T_PVa(:,:)     = T_PVa(1,1)
  T_PVf(:,:)     = T_PVf(1,1)
  T_cell(:,:)    = T_cell(1,1)
  T_PVb(:,:)     = T_PVb(1,1)
  ptsfc_eff(:,:) = ptsfc_eff(1,1)
  qvsfc_eff(:,:) = qvsfc_eff(1,1)

  ! Vertical grid: reference level (scalar level 2) at 3.3*h_PV
  zp(:,:,1) = -2.0 * 3.3 * h_PV(1,1)
  zp(:,:,2) =  0.0
  zp(:,:,3) =  2.0 * 3.3 * h_PV(1,1)

  !=====================================================================
  ! OPEN OUTPUT CSVs
  !=====================================================================

  OPEN(UNIT=20, FILE='pv_land_output.csv', STATUS='replace')
  WRITE(20,'(A)') 'timestep,H_meas,H_model,LE_meas,LE_model,' // &
                  'Rn_meas,Rn_model,G_meas,G_model,' // &
                  'T_nat,T_soil,T_PVf,T_cell,T_PVb,T_can,T_PVa,T_ref,' // &
                  'qsoil_surf,qsoil_deep,' // &
                  'SW_out_sky,LW_out_sky,H_nat,H_PVa_ref,H_can_ref,' // &
                  'G_SW_nat,J_SW_nat,G_LW_nat,J_LW_nat'

  OPEN(UNIT=30, FILE='pv_diag_energy_balance.csv', STATUS='replace')
  WRITE(30,'(A)') 'timestep,SWin,' // &
                  'T_PVf,T_PVb,T_cell,T_can,T_PVa,T_ref,' // &
                  'qLW_PVf,qLW_PVb,' // &
                  'H_PVf,H_PVb,' // &
                  'qcond_PVf,qcond_PVb,' // &
                  'qSW_PVf,qSW_PVb,' // &
                  'qPV,' // &
                  'cdh_can,cdh_PVa,' // &
                  'U_can,U_ref,' // &
                  'JLW_PVf,JLW_PVb'

  !=====================================================================
  ! MAIN TIME LOOP: forcing intervals (30 min), dtsfc sub-steps
  !=====================================================================

  WRITE(*,*) 'Starting time loop...'

  DO temp = 1, ntot-1

    IF (MOD(temp,48) == 0) WRITE(*,'(A,I6,A,I6)') '  Timestep ', temp, ' / ', ntot-1

    ! Obukhov length: one value per forcing interval (NOT interpolated, since
    ! L passes through +-infinity at neutral). L>0 stable, L<0 unstable.
    L_obk_in(:,:) = Lobk_input(temp)
    prcrate(:,:)  = rain_input(temp)

    DO it = 1, nstep_frc
      frac = REAL(it) / REAL(nstep_frc)
      CALL set_atmosphere(temp, temp+1, frac)
      nstep = nstep + 1

      ! Outside the PV park (pv_park_opt = 2): inputs that ARPS provides to
      ! the natural-ground scheme (sfcphy = 4)
      !   - surface wind speed, as in sfcphysics
      !   - net radiation of the ground, as the ARPS radiation scheme, with
      !     the albedo and emissivity of the natural ground (&pvland)
      wdspsfc(:,:) = MAX(vsfcmin, SQRT(u(:,:,2)**2 + v(:,:,2)**2 + w(:,:,2)**2))
      rnflx(:,:)   = (1.0 - albedo_nat(:,:)) * radsw(:,:)                     &
                   + emiss_nat(:,:) * (radlwin(:,:) - sbcst*tsoil(:,:,1,0)**4)
      CALL pvland_sfcphy(nx,ny,nz,nzsoil,nstyps,                       &
           u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,zp,prcrate,            &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           radsw,rnflx,radlwin,psfc,rhosfc,wdspsfc,0.0,0.0,             &
           cdm,cdh,cdq,usflx,vsflx,ptsflx,qvsflx,                       &
           pt1,tair,qvair,                                              &
           1,S_diffuse_in,lambda_in,phi_in,L_obk_in)
    END DO

    !-------------------------------------------------------------------
    ! Write outputs at the end of the forcing interval
    !-------------------------------------------------------------------
    WRITE(20, '(I0,27(",",F0.4))') temp+1, &
        H_meas(temp+1), H_ref(1,1), &
        LE_meas(temp+1), LE_nat(1,1), &
        Rn_meas(temp+1), Rn_ref(1,1), &
        G_meas(temp+1), G_nat(1,1), &
        tsoil(1,1,1,1) - 273.15, &
        tsoil(1,1,2,1) - 273.15, &
        T_PVf(1,1) - 273.15, &
        T_cell(1,1) - 273.15, &
        T_PVb(1,1) - 273.15, &
        T_can(1,1) - 273.15, &
        T_PVa(1,1) - 273.15, &
        tref_input(temp+1) - 273.15, &
        qsoil(1,1,1,1), &
        qsoil(1,1,2,1), &
        SW_out_sky(1,1), &
        LW_out_sky(1,1), &
        H_nat(1,1), &
        H_PVa_ref(1,1), &
        H_can_ref(1,1), &
        G_SW_nat(1,1), &
        J_SW_nat(1,1), &
        G_LW_nat(1,1), &
        J_LW_nat(1,1)

    WRITE(30, '(I0,22(",",ES13.5))') temp+1, &
        swin_input(temp+1), &
        T_PVf(1,1) - 273.15, &
        T_PVb(1,1) - 273.15, &
        T_cell(1,1) - 273.15, &
        T_can(1,1) - 273.15, &
        T_PVa(1,1) - 273.15, &
        tref_input(temp+1) - 273.15, &
        qLW_PVf(1,1), &
        qLW_PVb(1,1), &
        H_PVf(1,1), &
        H_PVb(1,1), &
        qcond_PVf(1,1), &
        qcond_PVb(1,1), &
        qSW_PVf(1,1), &
        qSW_PVb(1,1), &
        qPV(1,1), &
        cdh_can_diag(1,1), &
        cdh_PVa_diag(1,1), &
        U_can(1,1), &
        uref_diag(1,1), &
        JLW_PVf(1,1), &
        JLW_PVb(1,1)

  END DO

  CLOSE(20)
  CLOSE(30)

  WRITE(*,*) ''
  WRITE(*,*) '========================================='
  WRITE(*,*) 'Simulation complete!'
  WRITE(*,*) 'Output: pv_land_output.csv'
  WRITE(*,*) 'Output: pv_diag_energy_balance.csv'
  WRITE(*,*) '========================================='

CONTAINS

  !---------------------------------------------------------------------
  ! Build the ARPS-like single-column state from the forcing, linearly
  ! interpolated between rows k1 and k2 (frac = 0 -> k1, 1 -> k2).
  ! Wind: u = U sin(dir), v = U cos(dir), so that DRAG_PVLAND recovers
  ! dir = ATAN2(u,v) (same convention as the former offline model).
  !---------------------------------------------------------------------
  SUBROUTINE set_atmosphere(k1, k2, frac)
    INTEGER, INTENT(IN) :: k1, k2
    REAL,    INTENT(IN) :: frac

    uref  = uref_input(k1)  + (uref_input(k2)  - uref_input(k1))  * frac
    udir  = udir_input(k1)  + (udir_input(k2)  - udir_input(k1))  * frac
    tref  = tref_input(k1)  + (tref_input(k2)  - tref_input(k1))  * frac
    qvref = qvref_input(k1) + (qvref_input(k2) - qvref_input(k1)) * frac

    psfc(:,:)    = pref_input(k1)   + (pref_input(k2)   - pref_input(k1))   * frac
    rhosfc(:,:)  = rhoref_input(k1) + (rhoref_input(k2) - rhoref_input(k1)) * frac
    radsw(:,:)   = swin_input(k1)   + (swin_input(k2)   - swin_input(k1))   * frac
    radlwin(:,:) = lwin_input(k1)   + (lwin_input(k2)   - lwin_input(k1))   * frac
    S_diffuse_in(:,:) = S_diffuse_input(k1) + (S_diffuse_input(k2) - S_diffuse_input(k1)) * frac
    lambda_in(:,:)    = lambda_input(k1)    + (lambda_input(k2)    - lambda_input(k1))    * frac
    phi_in(:,:)       = phi_input(k1)       + (phi_input(k2)       - phi_input(k1))       * frac

    pbar(:,:,:)  = psfc(1,1)
    pprt(:,:,:)  = 0.
    ptbar(:,:,:) = tref * (p0/psfc(1,1))**rddcp
    ptprt(:,:,:) = 0.
    qvbar(:,:,:) = qvref
    qv(:,:,:)    = qvref
    u(:,:,:)     = uref * SIN(udir * pi / 180.0)
    v(:,:,:)     = uref * COS(udir * pi / 180.0)
    w(:,:,:)     = 0.
  END SUBROUTINE set_atmosphere

END PROGRAM pvland_offline_main
