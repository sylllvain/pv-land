!##################################################################
!##################################################################
!######                                                      ######
!######              PROGRAM PV_LAND_MAIN                    ######
!######                                                      ######
!######   Main driver program:                               ######
!######   - Read CSV forcing data                            ######
!######   - Loop over timesteps                              ######
!######   - Call pv_land subroutine                          ######
!######   - Write CSV output                                 ######
!######                                                      ######
!######   COMBINED version: writes both main output CSV      ######
!######   (with nat diagnostics) AND PV diag CSV             ######
!######                                                      ######
!##################################################################
!##################################################################

PROGRAM pv_land_main

  IMPLICIT NONE

  ! ====== Grid parameters ======
  INTEGER, PARAMETER :: nx = 2, ny = 2, nz = 3
  INTEGER, PARAMETER :: nzsoil = 2
  INTEGER, PARAMETER :: nstyps = 1
  INTEGER, PARAMETER :: ntot_max = 100000  ! Upper bound for array allocation
  INTEGER :: ntot                           ! Actual number of timesteps (read from CSV)
  
  ! ====== PV park geometry ======
  REAL(8), PARAMETER :: L = 3.23d0           ! Panel width (m)
  REAL(8), PARAMETER :: a = 3.4d0            ! Inter-row spacing (m)
  REAL(8), PARAMETER :: beta_deg = 24.0d0    ! Tilt angle (degrees)
  REAL(8), PARAMETER :: h_PV = 1.9d0         ! Panel height (m)
  REAL(8), PARAMETER :: theta_0_deg = -3.0d0 ! Park orientation angle (degrees)
  REAL(8)            :: h_clear              ! Clearance height (m)
  REAL(8)            :: beta_rad             ! Tilt angle (rad)
  REAL(8)            :: theta_0              ! Park orientation angle (rad)

  ! ====== Surface parameters ======
  REAL(8), PARAMETER :: z0_nat = 0.05d0        ! Natural surface roughness (m)
  REAL(8), PARAMETER :: z0_PV_back = 0.0002d0   ! PV panel back roughness (m)
  REAL(8), PARAMETER :: z0_PV_front = 0.0002d0  ! PV panel front roughness (m)
  
  REAL, PARAMETER :: pi = 3.14159265359

  ! Optical properties
  REAL(8), PARAMETER :: albedo_nat = 0.22d0  ! Natural surface albedo
  REAL(8), PARAMETER :: albedo_PVf = 0.11d0  ! PV front albedo
  REAL(8), PARAMETER :: albedo_PVb = 0.11d0  ! PV back albedo
  REAL(8), PARAMETER :: emiss_nat = 0.98d0   ! Natural surface emissivity
  REAL(8), PARAMETER :: emiss_PVf = 0.79d0   ! PV front emissivity
  REAL(8), PARAMETER :: emiss_PVb = 0.80d0   ! PV back emissivity
  
  ! ====== File I/O ======
  INTEGER :: temp, ios, t
  CHARACTER(LEN=500) :: line
  
  ! ====== CSV input arrays ======
  REAL(8) :: uref_input(ntot_max), udir_input(ntot_max)
  REAL(8) :: tref_input(ntot_max), qvref_input(ntot_max)
  REAL(8) :: rain_input(ntot_max)
  REAL(8) :: swout_input(ntot_max), swin_input(ntot_max)
  REAL(8) :: lwout_input(ntot_max), lwin_input(ntot_max)
  REAL(8) :: pref_input(ntot_max), rhoref_input(ntot_max)
  REAL(8) :: phi_input(ntot_max), lambda_input(ntot_max), S_diffuse_input(ntot_max)
  REAL(8) :: H_meas(ntot_max), LE_meas(ntot_max), Rn_meas(ntot_max), G_meas(ntot_max)
  REAL(8) :: Lobk_input(ntot_max)          ! Obukhov length L [m] (last CSV column)
  
  ! ====== Interpolation arrays ======
  REAL(8) :: urefi(nx,ny,2), udiri(nx,ny,2)
  REAL(8) :: trefi(nx,ny,2), qvrefi(nx,ny,2)
  REAL(8) :: prcrate(nx,ny)
  REAL(8) :: L_obk_in(nx,ny)               ! Obukhov length passed to pv_land [m]
  REAL(8) :: swouti(nx,ny,2), swini(nx,ny,2)
  REAL(8) :: lwouti(nx,ny,2), lwini(nx,ny,2)
  REAL(8) :: prefi(nx,ny,2), rhorefi(nx,ny,2)
  REAL(8) :: phii(nx,ny,2), lambdai(nx,ny,2), S_diffusei(nx,ny,2)
  
  ! ====== State variables ======
  REAL(8) :: T_nat(nx,ny), T_soil(nx,ny)           ! Temperatures (K)
  REAL(8) :: T_can(nx,ny), T_PVa(nx,ny)
  REAL(8) :: T_PVf(nx,ny), T_cell(nx,ny), T_PVb(nx,ny)
  REAL(8) :: qsoil_surf(nx,ny), qsoil_deep(nx,ny)  ! Soil moisture (m3/m3)
  REAL(8) :: wetcanp(nx,ny,0:nstyps)                ! Canopy water amount
  REAL(8) :: qvsfc(nx,ny,0:nstyps)                  ! Effective specific humidity at surface
  
  ! ====== Surface properties ======
  INTEGER :: soiltyp(nx,ny,nstyps), vegtyp(nx,ny)
  REAL(8) :: stypfrct(nx,ny,nstyps)
  REAL(8) :: veg(nx,ny), lai(nx,ny)
  
  ! ====== Output variables ======
  REAL(8) :: H_ref_out(nx,ny), LE_ref_out(nx,ny)
  REAL(8) :: Rn_ref_out(nx,ny), G_out(nx,ny)
  REAL(8) :: SW_out_sky(nx,ny), LW_out_sky(nx,ny)
  REAL(8) :: H_nat_out(nx,ny)
  
  ! Nat surface diagnostics (from _nat version)
  REAL(8) :: H_PVa_ref_out(nx,ny), H_can_ref_out(nx,ny)
  REAL(8) :: G_SW_nat_out(nx,ny), J_SW_nat_out(nx,ny)
  REAL(8) :: G_LW_nat_out(nx,ny), J_LW_nat_out(nx,ny)
  
  ! ====== PV energy balance diagnostics (from _pv version) ======
  REAL(8) :: diag_qLW_PVf(nx,ny), diag_qLW_PVb(nx,ny)
  REAL(8) :: diag_H_PVf(nx,ny), diag_H_PVb(nx,ny)
  REAL(8) :: diag_qcond_PVf(nx,ny), diag_qcond_PVb(nx,ny)
  REAL(8) :: diag_qSW_PVf(nx,ny), diag_qSW_PVb(nx,ny)
  REAL(8) :: diag_qPV(nx,ny)
  REAL(8) :: diag_cdh_can(nx,ny), diag_cdh_PVa(nx,ny)
  REAL(8) :: diag_T_can(nx,ny), diag_T_PVa(nx,ny)
  REAL(8) :: diag_U_can(nx,ny), diag_U_ref(nx,ny)
  REAL(8) :: diag_JLW_PVf(nx,ny), diag_JLW_PVb(nx,ny)
  
  INCLUDE 'timelvls.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'
  
  !=====================================================================
  ! INITIALIZATION
  !=====================================================================
  
  WRITE(*,*) '========================================='
  WRITE(*,*) 'PV-LAND Model (COMBINED version)'
  WRITE(*,*) '========================================='
  WRITE(*,*) ''
  

  beta_rad = beta_deg * pi / 180.0d0
  h_clear = h_PV - L * SIN(beta_rad)
  theta_0 = theta_0_deg * pi / 180.0d0

  ! Surface properties
  veg(:,:) = 0.95d0
  vegtyp(:,:) = 3
  lai(:,:) = 2d0
  soiltyp(:,:,:) = 6
  stypfrct(:,:,:) = 1.0d0
  
  ! Time stepping
  dtsfc = 2d0
  moist = 1
  
  !=====================================================================
  ! READ CSV INPUT
  !=====================================================================
  
  WRITE(*,*) 'Reading forcing data from CSV...'
  
  OPEN(UNIT=10, FILE='recorded_data_adjusted_daily_summer_2025.csv', STATUS='old')
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
    
	
    tref_input(temp) = tref_input(temp) + 273.15d0
    IF (ios /= 0) EXIT
    temp = temp + 1
    IF (temp > ntot_max) EXIT
  END DO
  ntot = temp - 1
  CLOSE(10)

  WRITE(*,*) 'CSV read: ', ntot, ' timesteps'
  WRITE(*,*) ''
  
  !=====================================================================
  ! INITIALIZE STATE VARIABLES
  !=====================================================================
  
  T_nat(:,:)  = tref_input(1)
  T_soil(:,:) = 20.95+273.15
  T_can(:,:)  = tref_input(1)
  T_PVa(:,:)  = tref_input(1)
  T_PVf(:,:)  = tref_input(1)
  T_cell(:,:) = tref_input(1)
  T_PVb(:,:)  = tref_input(1)
  
  qsoil_surf(:,:) = 0.11
  qsoil_deep(:,:) = 0.10
  wetcanp(:,:,:) = 0
  qvsfc(:,:,:) = qvref_input(1)
  
  !=====================================================================
  ! OPEN OUTPUT CSVs
  !=====================================================================
  
  ! Main output CSV (union of _nat and _pv columns)
  OPEN(UNIT=20, FILE='pv_land_output.csv', STATUS='replace')
  WRITE(20,'(A)') 'timestep,H_meas,H_model,LE_meas,LE_model,' // &
                  'Rn_meas,Rn_model,G_meas,G_model,' // &
                  'T_nat,T_soil,T_PVf,T_cell,T_PVb,T_can,T_PVa,T_ref,' // &
                  'qsoil_surf,qsoil_deep,' // &
                  'SW_out_sky,LW_out_sky,H_nat,H_PVa_ref,H_can_ref,' // &
                  'G_SW_nat,J_SW_nat,G_LW_nat,J_LW_nat'
  
  ! Diagnostic CSV for PV energy balance analysis (from _pv version)
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
  ! MAIN TIME LOOP (LARGE TIMESTEPS - 30 min)
  !=====================================================================
  
  WRITE(*,*) 'Starting time loop...'
  WRITE(*,*) ''
  
  DO temp = 1, ntot-1

	print *, "Start large loop, temp = ", temp
    IF (MOD(temp,6) == 0) THEN
      WRITE(*,'(A,I3,A,I3)') '  Timestep ', temp, ' / ', ntot-1
    END IF
    
    ! Setup interpolation arrays (two time levels)
    urefi(:,:,1) = uref_input(temp)
    urefi(:,:,2) = uref_input(temp+1)
    udiri(:,:,1) = udir_input(temp)
    udiri(:,:,2) = udir_input(temp+1)
    trefi(:,:,1) = tref_input(temp)
    trefi(:,:,2) = tref_input(temp+1)
    qvrefi(:,:,1) = qvref_input(temp)
    qvrefi(:,:,2) = qvref_input(temp+1)
    prcrate(:,:) = rain_input(temp)
    ! Obukhov length: single value per large timestep (NOT interpolated, since L
    ! passes through +-infinity at neutral and linear interpolation would be
    ! spurious). Sign convention: L>0 stable, L<0 unstable (matches profil.f90).
    L_obk_in(:,:) = Lobk_input(temp)
    swouti(:,:,1) = swout_input(temp)
    swouti(:,:,2) = swout_input(temp+1)
    swini(:,:,1) = swin_input(temp)
    swini(:,:,2) = swin_input(temp+1)
    lwouti(:,:,1) = lwout_input(temp)
    lwouti(:,:,2) = lwout_input(temp+1)
    lwini(:,:,1) = lwin_input(temp)
    lwini(:,:,2) = lwin_input(temp+1)
    prefi(:,:,1) = pref_input(temp)
    prefi(:,:,2) = pref_input(temp+1)
    rhorefi(:,:,1) = rhoref_input(temp)
    rhorefi(:,:,2) = rhoref_input(temp+1)
    phii(:,:,1) = phi_input(temp)
    phii(:,:,2) = phi_input(temp+1)
    lambdai(:,:,1) = lambda_input(temp)
    lambdai(:,:,2) = lambda_input(temp+1)
    S_diffusei(:,:,1) = S_diffuse_input(temp)
    S_diffusei(:,:,2) = S_diffuse_input(temp+1)

    !-------------------------------------------------------------------
    ! Call PV-LAND subroutine (combined version)
    !-------------------------------------------------------------------
    CALL pv_land(nx, ny, nz, nzsoil, nstyps, &
                L, a, beta_deg, h_PV, h_clear, theta_0, &
                z0_nat, z0_PV_back, z0_PV_front, &
                albedo_nat, albedo_PVf, albedo_PVb, &
                emiss_nat, emiss_PVf, emiss_PVb, &
                urefi, udiri, trefi, qvrefi, prcrate, L_obk_in, &
                swouti, swini, lwouti, lwini, prefi, rhorefi, &
                phii, lambdai, S_diffusei, &
                T_nat, T_soil, T_can, T_PVa, &
                T_PVf, T_cell, T_PVb, &
                qsoil_surf, qsoil_deep, wetcanp, qvsfc, &
                soiltyp, vegtyp, lai, veg, stypfrct, &
                H_ref_out, LE_ref_out, Rn_ref_out, G_out, &
                SW_out_sky, LW_out_sky, H_nat_out, &
                H_PVa_ref_out, H_can_ref_out, &
                G_SW_nat_out, J_SW_nat_out, &
                G_LW_nat_out, J_LW_nat_out, &
                diag_qLW_PVf, diag_qLW_PVb, &
                diag_H_PVf, diag_H_PVb, &
                diag_qcond_PVf, diag_qcond_PVb, &
                diag_qSW_PVf, diag_qSW_PVb, &
                diag_qPV, &
                diag_cdh_can, diag_cdh_PVa, &
                diag_T_can, diag_T_PVa, &
                diag_U_can, diag_U_ref, &
                diag_JLW_PVf, diag_JLW_PVb)
    
	
    !-------------------------------------------------------------------
    ! Write main output
    !-------------------------------------------------------------------
    WRITE(20, '(I0,27(",",F0.4))') temp+1, &
		H_meas(temp+1), H_ref_out(1,1), &
		LE_meas(temp+1), LE_ref_out(1,1), &
		Rn_meas(temp+1), Rn_ref_out(1,1), &
		G_meas(temp+1), G_out(1,1), &
		T_nat(1,1) - 273.15d0, &
		T_soil(1,1) - 273.15d0, &
		T_PVf(1,1) - 273.15d0, &
		T_cell(1,1) - 273.15d0, &
		T_PVb(1,1) - 273.15d0, &
		T_can(1,1) - 273.15d0, &
		T_PVa(1,1) - 273.15d0, &
		tref_input(temp+1) - 273.15d0, &
		qsoil_surf(1,1), &
		qsoil_deep(1,1), &
		SW_out_sky(1,1), &
		LW_out_sky(1,1), &
		H_nat_out(1,1), &
		H_PVa_ref_out(1,1), &
		H_can_ref_out(1,1), &
		G_SW_nat_out(1,1), &
		J_SW_nat_out(1,1), &
		G_LW_nat_out(1,1), &
		J_LW_nat_out(1,1)

    ! Write PV energy balance diagnostics
    WRITE(30, '(I0,22(",",ES13.5))') temp+1, &
        swin_input(temp+1), &
        T_PVf(1,1) - 273.15d0, &
        T_PVb(1,1) - 273.15d0, &
        T_cell(1,1) - 273.15d0, &
        diag_T_can(1,1) - 273.15d0, &
        diag_T_PVa(1,1) - 273.15d0, &
        tref_input(temp+1) - 273.15d0, &
        diag_qLW_PVf(1,1), &
        diag_qLW_PVb(1,1), &
        diag_H_PVf(1,1), &
        diag_H_PVb(1,1), &
        diag_qcond_PVf(1,1), &
        diag_qcond_PVb(1,1), &
        diag_qSW_PVf(1,1), &
        diag_qSW_PVb(1,1), &
        diag_qPV(1,1), &
        diag_cdh_can(1,1), &
        diag_cdh_PVa(1,1), &
        diag_U_can(1,1), &
        diag_U_ref(1,1), &
        diag_JLW_PVf(1,1), &
        diag_JLW_PVb(1,1)


  END DO
  
  CLOSE(20)
  CLOSE(30)
  
  !=====================================================================
  ! FINALIZE
  !=====================================================================
  
  WRITE(*,*) ''
  WRITE(*,*) '========================================='
  WRITE(*,*) 'Simulation complete!'
  WRITE(*,*) 'Output: pv_land_output.csv'
  WRITE(*,*) 'Output: pv_diag_energy_balance.csv'
  WRITE(*,*) '========================================='

END PROGRAM pv_land_main
