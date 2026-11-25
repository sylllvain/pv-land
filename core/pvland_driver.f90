!########################################################################
!########################################################################
!#########                                                      #########
!#########             MODULE pvland_sfcphy_work                #########
!#########                                                      #########
!########################################################################
!########################################################################
!
! Work arrays of SUBROUTINE pvland_sfcphy, allocated at the first call
! and kept between calls. They are in a module (instead of local SAVE
! arrays) so that drivers can read diagnostics such as uref, cdh_can
! or cdh_PVa (offline model output).  [SYLVAIN]

MODULE pvland_sfcphy_work

  IMPLICIT NONE

  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdm_u,cdm_v
  REAL, ALLOCATABLE, DIMENSION (:,:) :: z0_PV

  REAL, ALLOCATABLE, DIMENSION (:,:) :: uref,T_nat,T_soil,qsoil_surf,&
	qsoil_deep,S_diffuse,S_direct,wrmax
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_nat,cdq_nat,cdm_nat
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_PVa,cdq_PVa,cdm_PVa
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_PVb,cdq_PVb,cdm_PVb
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_PVf,cdq_PVf,cdm_PVf
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_can,cdq_can,cdm_can

  REAL, ALLOCATABLE, DIMENSION (:,:) :: rhs_nat,rhs_soil,rhs_PVf,rhs_cell,&
	rhs_PVb,rhs_wsfc,rhs_wdp,rhs_wcnp
  REAL, ALLOCATABLE, DIMENSION (:,:) :: frc_nat,frc_soil,frc_PVf,frc_cell,&
	frc_PVb,frc_wsfc,frc_wdp,frc_wcnp
  REAL, ALLOCATABLE, DIMENSION (:,:) :: T_nat_half,T_soil_half,T_PVf_half,&
	T_cell_half,T_PVb_half,qsoil_surf_half,qsoil_deep_half
  REAL, ALLOCATABLE, DIMENSION (:,:,:) :: wetcanp_half

  REAL, ALLOCATABLE, DIMENSION (:,:) :: SW_out_sky_loc, LW_out_sky_loc
  REAL, ALLOCATABLE, DIMENSION (:,:) :: dum_qLW_PVf, dum_qLW_PVb, &
                             dum_H_PVf, dum_H_PVb, &
                             dum_qcond_PVf, dum_qcond_PVb, &
                             dum_qSW_PVf, dum_qSW_PVb, &
                             dum_qPV, &
                             dum_cdh_can, dum_cdh_PVa, &
                             dum_T_can, dum_T_PVa, &
                             dum_JLW_PVf, dum_JLW_PVb, &
                             dum_G_SW_nat, dum_J_SW_nat, &
                             dum_G_LW_nat, dum_J_LW_nat
  REAL, ALLOCATABLE, DIMENSION (:,:) :: lambda_deg,phi_deg
  ! outside the park: surface state used by dragcoef and fluxes at t^n
  ! (as sfcphy = 4)
  REAL, ALLOCATABLE, DIMENSION (:,:) :: ptsfc_drag,qvsfc_drag,ptsflx_n,qvsflx_n
  ! version 2: friction velocity at z_ref, canopy-air humidity budget
  REAL, ALLOCATABLE, DIMENSION (:,:) :: u_star_pv, q_can_loc, gcan_loc, &
                             vpd_f3_loc, vpd_ref_loc

END MODULE pvland_sfcphy_work

!########################################################################
!########################################################################
!#########                                                      #########
!#########             SUBROUTINE pvland_sfcphy                 #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE pvland_sfcphy(nx,ny,nz,nzsoil,nstyps,                         &
           u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,zp,prcrate,            &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           radsw,rnflx,radlwin,psfc,rhosfc,wdspsfc,gumove,gvmove,       &
           cdm,cdh,cdq,usflx,vsflx,ptsflx,qvsflx,                       &
           pt1,tair,qvair,                                              &
           meas_opt,S_diffuse_in,lambda_in,phi_in,diag_out)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  One PV-LAND surface time step (dtsfc): exchange coefficients, surface
!  fluxes and RK2 (predictor-corrector, as in soilebm) integration of
!  the soil, PV panel and air node temperatures and of soil moisture.
!
!  Called by ARPS (sfcphysics, sfcphy == 5) and by the PV-LAND offline
!  model: both use exactly this routine.
!
!  meas_opt = 0 (ARPS): the diffuse solar radiation (split_sw) and the
!                 sun angles in the park frame (compute_lambda_phi) are
!                 computed; S_diffuse_in, lambda_in, phi_in are not used.
!  meas_opt = 1 (offline): they are taken from S_diffuse_in (W/m2),
!                 lambda_in, phi_in (rad, park frame).
!  The Obukhov length L_obk is computed by PV-LAND from its own fluxes
!  (obukhov_update, version 2), in ARPS and offline.
!
!  Outside the PV park (pv_park = 0) the surface layer is computed as
!  with sfcphy = 4 (dragcoef + sfcflx): reference wind wdspsfc, surface
!  air tair, qvair averaged over levels 1-2, drag coefficients from tsoil
!  and qvsfc, atmospheric fluxes at the beginning of the step (t^n),
!  fluxes weighted by stypfrct over the soil types, qvsfc updated; the
!  soil itself uses the PV-LAND soil model (soilebm_frc_pv, version 2:
!  drainage, pv_wref, namelist soil parameters), with nsfcst sub-steps.
!
!  AUTHOR: moved from sfcphy3d.f90 (sfcphy == 5 block)  [SYLVAIN]
!
!-----------------------------------------------------------------------
!
  USE module_pvland_in
  USE module_pvland_out
  USE rad_pvland
  USE rad_split
  USE pvland_sfcphy_work
  USE pvland_nml, ONLY : vf_update_opt, pv_zd_nat, pv_tau_L, pv_L_floor, &
                         pv_wsat

  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nzsoil,nstyps

  REAL :: u     (nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
  REAL :: ptprt (nx,ny,nz), pprt(nx,ny,nz), qv(nx,ny,nz)
  REAL :: ptbar (nx,ny,nz), pbar(nx,ny,nz), qvbar(nx,ny,nz)
  REAL :: zp    (nx,ny,nz)  ! Height of the w-points (m)
  REAL :: prcrate(nx,ny)    ! Precipitation rate (kg/(m**2*s))

  INTEGER :: soiltyp (nx,ny,nstyps), vegtyp(nx,ny)
  REAL    :: stypfrct(nx,ny,nstyps)
  REAL    :: lai(nx,ny), roufns(nx,ny), veg(nx,ny)

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps)  ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)         ! Canopy water amount
  REAL :: snowdpth(nx,ny)                 ! Snow depth
  REAL :: qvsfc  (nx,ny,0:nstyps)         ! Effective S.H. at sfc.

  REAL :: radsw  (nx,ny)    ! Incoming solar radiation
  REAL :: rnflx  (nx,ny)    ! Net radiation (used outside the park only)
  REAL :: radlwin(nx,ny)    ! Incoming longwave radiation
  REAL :: psfc   (nx,ny)    ! Surface pressure (Pascal)
  REAL :: rhosfc (nx,ny)    ! Air density at the surface
  REAL :: wdspsfc(nx,ny)    ! Surface wind speed of sfcphysics (used outside the park)
  REAL :: gumove, gvmove    ! Domain translation speed (x, y)

  REAL :: cdm(nx,ny), cdh(nx,ny), cdq(nx,ny)   ! Drag coefficients (output)
  REAL :: usflx(nx,ny), vsflx(nx,ny)           ! Momentum fluxes (in: previous step, out)
  REAL :: ptsflx(nx,ny), qvsflx(nx,ny)         ! Heat/moisture fluxes (in: previous step, out)

  REAL :: pt1    (nx,ny)    ! Potential temperature at 1st model level AGL (output)
  REAL :: tair   (nx,ny)    ! Temperature of surface air (output)
  REAL :: qvair  (nx,ny)    ! Specific humidity of surface air (output)

  INTEGER :: meas_opt           ! 0: computed (ARPS), 1: measured forcing (offline)
  REAL :: S_diffuse_in(nx,ny)   ! measured diffuse SW (W/m2)        (meas_opt = 1)
  REAL :: lambda_in(nx,ny)      ! measured lambda, park frame (rad) (meas_opt = 1)
  REAL :: phi_in(nx,ny)         ! measured phi, park frame (rad)    (meas_opt = 1)
  LOGICAL :: diag_out           ! .TRUE.: compute all the diagnostics of the step (output step);
                                ! .FALSE.: end-of-step fluxes only (faster)
!
!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: pi = 3.14159265359

  REAL :: z_ref,usuf,vsuf,Uplus,tema
  INTEGER :: i, j, k, is, it
  REAL :: pts_type, qvs_type   ! heat/moisture fluxes of one soil type, end of step

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'
  INCLUDE 'soilcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ! PV-LAND soil parameter (namelist &pvland), as in soilebm_frc_pv
  IF (pv_wsat > 0.) wsat(:) = pv_wsat

  ! work arrays (module pvland_sfcphy_work) allocated at the first call only
  IF (.NOT. ALLOCATED(cdm_u)) THEN
  ALLOCATE(cdm_u(nx,ny),cdm_v(nx,ny),z0_PV(nx,ny))
  ALLOCATE(uref(nx,ny), T_nat(nx,ny), T_soil(nx,ny), qsoil_surf(nx,ny), &
	qsoil_deep(nx,ny),S_diffuse(nx,ny),S_direct(nx,ny),wrmax(nx,ny))
  ALLOCATE(cdh_nat(nx,ny),cdq_nat(nx,ny),cdm_nat(nx,ny))
  ALLOCATE(cdh_PVa(nx,ny),cdq_PVa(nx,ny),cdm_PVa(nx,ny))
  ALLOCATE(cdh_PVb(nx,ny),cdq_PVb(nx,ny),cdm_PVb(nx,ny))
  ALLOCATE(cdh_PVf(nx,ny),cdq_PVf(nx,ny),cdm_PVf(nx,ny))
  ALLOCATE(cdh_can(nx,ny),cdq_can(nx,ny),cdm_can(nx,ny))
  
  ALLOCATE(rhs_nat(nx,ny),rhs_soil(nx,ny),rhs_PVf(nx,ny),rhs_cell(nx,ny),&
	rhs_PVb(nx,ny),rhs_wsfc(nx,ny),rhs_wdp(nx,ny),rhs_wcnp(nx,ny))
  ALLOCATE(frc_nat(nx,ny),frc_soil(nx,ny),frc_PVf(nx,ny),frc_cell(nx,ny),&
	frc_PVb(nx,ny),frc_wsfc(nx,ny),frc_wdp(nx,ny),frc_wcnp(nx,ny))
  ALLOCATE(T_nat_half(nx,ny),T_soil_half(nx,ny),T_PVf_half(nx,ny),&
	T_cell_half(nx,ny),T_PVb_half(nx,ny),qsoil_surf_half(nx,ny),qsoil_deep_half(nx,ny),wetcanp_half(nx,ny,0:nstyps))
	
  ALLOCATE(SW_out_sky_loc(nx,ny), LW_out_sky_loc(nx,ny))	
  ALLOCATE(dum_qLW_PVf(nx,ny), dum_qLW_PVb(nx,ny), &
                             dum_H_PVf(nx,ny), dum_H_PVb(nx,ny), &
                             dum_qcond_PVf(nx,ny), dum_qcond_PVb(nx,ny), &
                             dum_qSW_PVf(nx,ny), dum_qSW_PVb(nx,ny), &
                             dum_qPV(nx,ny), &
                             dum_cdh_can(nx,ny), dum_cdh_PVa(nx,ny), &
                             dum_T_can(nx,ny), dum_T_PVa(nx,ny), &
                             dum_JLW_PVf(nx,ny), dum_JLW_PVb(nx,ny), &
                             dum_G_SW_nat(nx,ny), dum_J_SW_nat(nx,ny), &
                             dum_G_LW_nat(nx,ny), dum_J_LW_nat(nx,ny))
  ALLOCATE(lambda_deg(nx,ny),phi_deg(nx,ny))
  ALLOCATE(ptsfc_drag(nx,ny),qvsfc_drag(nx,ny),ptsflx_n(nx,ny),qvsflx_n(nx,ny))
  ALLOCATE(u_star_pv(nx,ny), q_can_loc(nx,ny), gcan_loc(nx,ny), &
           vpd_f3_loc(nx,ny), vpd_ref_loc(nx,ny))
  ! canopy air not set by INIT_PVLAND (initopt = 1) or by a restart file:
  ! humidity of the first model level
  WHERE (q_can(:,:) < 0.0) q_can(:,:) = qv(:,:,2)
  END IF


    ! Moving panels (vf_update_opt = 1): geometry derived from the tilt
    ! beta (view factors, radiosity LU factors, h_clear, z_can, ...)
    ! recomputed at every time step. A tracking law updating beta(i,j)
    ! would be applied just before this call.  [SYLVAIN]
    IF (vf_update_opt == 1) CALL pvland_geometry(nx,ny)

    DO j=1,ny-1
      DO i=1,nx-1
        pt1(i,j) = (ptbar(i,j,2)+ptprt(i,j,2))
        tair(i,j) = (ptbar(i,j,2)+ptprt(i,j,2))               &
                        * (psfc(i,j)/p0)**rddcp   ! Temperature
        qvair(i,j) = qv(i,j,2)
        
      END DO
    END DO

    ! outside the park, surface air as sfcphy = 4 (mean of levels 1-2)
    DO j=1,ny-1
      DO i=1,nx-1
        IF (pv_park(i,j).ne.1.) THEN
        tair(i,j) = 0.5 * ( (ptbar(i,j,1)+ptprt(i,j,1))                 &
                          + (ptbar(i,j,2)+ptprt(i,j,2)) )               &
                        * (psfc(i,j)/p0)**rddcp   ! Temperature
        qvair(i,j) = 0.5 * ( qv(i,j,1) + qv(i,j,2) )
        ENDIF
      END DO
    END DO

    tsoil(:,:,:,0) = 0.0
    qsoil(:,:,:,0) = 0.0
    wetcanp(:,:,0) = 0.0
    qvsfc(:,:,0)   = 0.0


    
    DO j=1,ny-1
      DO i=1,nx-1
      
      IF (pv_park(i,j).eq.1.) THEN
      
      IF (meas_opt == 0) THEN
      ! split solar radiation between direct and diffuse
      CALL split_sw(radsw(i,j), year, month, day, hour, minute, &
                       latitud, longitud, 1, &
                       S_diffuse(i,j), S_direct(i,j), 1367.)
                       
      ! Transform zenithal and azimutal angles accoridng to PV park orientation
      call compute_lambda_phi(szenith_deg(i,j), sazimuth_deg(i,j)+180., AZ_PARK(i,j), lambda_deg(i,j), phi_deg(i,j))
      ELSE
      ! measured forcing (offline model)
      S_diffuse(i,j)  = S_diffuse_in(i,j)
      S_direct(i,j)   = radsw(i,j) - S_diffuse(i,j)
      lambda_deg(i,j) = lambda_in(i,j) * 180. / pi
      phi_deg(i,j)    = phi_in(i,j)    * 180. / pi
      ENDIF
      

      ENDIF !pv_park
      

      END DO
    END DO
    
    !ATTENTION TEST
    !radsw=645.78
    !radlwin=407.69
    !S_diffuse=347.032867
    !S_direct=radsw-S_diffuse
    !u=-3.1634
    !v=-1.3095
    !w=0.
    !tair=296.844
    !qvair=0.013668
    !pt1=295.71524
    !psfc=101341.992
    !pbar=101341.992
    !T_can=297.7133
    !T_PVa=297.8247
    !T_PVb=310.9537
    !T_PVf=310.1482
    !T_cell=311.5155
    !tsoil(:,:,1,:)=298.2200
    !tsoil(:,:,2,:)=294.3197
    !qsoil(:,:,1,:)=0.103495
    !qsoil(:,:,2,:)=0.099695
    !phi_deg=2.34709883*180./pi
    !lambda_deg=0.49792251*180./pi
    !rhosfc=1.179566
    !L_obk=-29.1319
    
    ! fluxes accumulated over the soil types (weights stypfrct)
    usflx (:,:) = 0.0
    vsflx (:,:) = 0.0
    ptsflx(:,:) = 0.0
    qvsflx(:,:) = 0.0

    DO is=nstyp,1,-1

      DO i=1,nx-1
        DO j=1,ny-1
			T_nat(i,j) = tsoil(i,j,1,is)
			T_soil(i,j) = tsoil(i,j,2,is)
			qsoil_surf(i,j) = qsoil(i,j,1,is)
			qsoil_deep(i,j) = qsoil(i,j,2,is)			
        END DO
      END DO
    !-------------------------------------------------------------------
    ! 2. COMPUTE LOCAL EXCHANGE COEFFICIENTS cdm, cdh and cdq
    !-------------------------------------------------------------------

    CALL DRAG_PVLAND(nx,ny,nz,zp,u,v,w,pbar,pt1,ptsfc_eff,L_obk,&
        T_nat,rhosfc,roufns,soiltyp(1,1,is),uref,qvair,&
		cdh_nat,cdq_nat,cdm_nat,&
		cdh_PVa,cdq_PVa,cdm_PVa,&
		cdh_PVb,cdq_PVb,cdm_PVb,&
		cdh_PVf,cdq_PVf,cdm_PVf,&
		cdh_can,cdq_can,cdm_can,&
		z0_PV,u_star_pv)

    ! outside the park, as sfcphy = 4: wind wdspsfc, z0 = roufns (z0_PV),
    ! surface state of the soil type is
      DO j=1,ny-1
        DO i=1,nx-1
        IF (pv_park(i,j).eq.1.) THEN
          ptsfc_drag(i,j) = ptsfc_eff(i,j)
          qvsfc_drag(i,j) = qvsfc_eff(i,j)
        ELSE
          uref(i,j) = wdspsfc(i,j)
          ptsfc_drag(i,j) = tsoil(i,j,1,is)*(p0/psfc(i,j))**rddcp
          qvsfc_drag(i,j) = qvsfc(i,j,is)
        ENDIF
        END DO
      END DO

        CALL dragcoef(nx,ny,nz,zp,                                      &
             ptsfc_drag,qvsfc_drag,soiltyp(1,1,is),z0_PV,uref,pt1,  &
             cdm,cdh,cdq)
    
    !-------------------------------------------------------------------
    ! 3. COMPUTE USFLX, VSFLX, PTSFLX, QVSFLX
    !-------------------------------------------------------------------

      DO j=1,ny-1
        DO i=1,nx-1

      IF (pv_park(i,j).eq.1.) THEN
      
			z_ref = zp(i,j,2) + 0.5*(zp(i,j,3)-zp(i,j,2))

			! ATTENTION, il faut prendre en compte la thermique
			CALL compute_U_plus(z_ref, pi/2., a(i,j), beta(i,j)*180./pi, L(i,j), &
			   h_clear(i,j), theta_0(i,j), roufns(i,j), pv_zd_nat, L_obk(i,j), Uplus)
			cdm_u(i,j) = (1. / Uplus)**2
			
			CALL compute_U_plus(z_ref, 0., a(i,j), beta(i,j)*180./pi, L(i,j), &
			   h_clear(i,j), theta_0(i,j), roufns(i,j), pv_zd_nat, L_obk(i,j), Uplus)
			cdm_v(i,j) = (1. / Uplus)**2
			
	  ELSE
			cdm_u(i,j) = cdm(i,j)
			cdm_v(i,j) = cdm(i,j)
	  ENDIF
	  
        END DO
      END DO

	  DO j=1,ny-1
      DO i=2,nx-1
			usuf = 0.5*(u(i,j,2)+u(i,j,1))+gumove
			
			usflx(i,j)  = usflx(i,j)  + 0.125*(rhosfc(i,j)+rhosfc(i-1,j))* &
				 (cdm_u(i,j)+cdm_u(i-1,j))*(uref(i,j)+uref(i-1,j))*                         &
                 sign(1.0, usuf)*MAX(abs(usuf),vsfcmin)*stypfrct(i,j,is)
      END DO
      END DO
	  DO j=2,ny-1
      DO i=1,nx-1
			vsuf = 0.5*(v(i,j,2)+v(i,j,1))+gvmove
			
			vsflx(i,j)  = vsflx(i,j)  + 0.125*(rhosfc(i,j)+rhosfc(i,j-1))* &
				 (cdm_v(i,j)+cdm_v(i,j-1))*(uref(i,j)+uref(i,j-1))*                         &
                 sign(1.0, vsuf)*MAX(abs(vsuf),vsfcmin)*stypfrct(i,j,is)
      END DO
      END DO

      DO j=1,ny-1
        DO i=1,nx-1

        ! heat and moisture fluxes at t^n (sfcflx), used outside the park
        tema = rhosfc(i,j)*MAX(uref(i,j),vsfcmin)
        ptsflx_n(i,j)=cdh(i,j)*tema*(pt1(i,j)-ptsfc_drag(i,j))
        qvsflx_n(i,j)=cdq(i,j)*tema*(qv(i,j,2)-qvsfc_drag(i,j))
 	  
        END DO
      END DO

    !-------------------------------------------------------------------
    ! nsfcst sub-steps of dtsfc per call (dtbig = nsfcst*dtsfc), with the
    ! exchange coefficients of the call, as soilebm
    !-------------------------------------------------------------------

    DO it = 1, nsfcst

    !-------------------------------------------------------------------
    ! 4. SOLVE AIR NODE TEMPERATURES with current T_nat, T_PVb, T_PVf
    !-------------------------------------------------------------------

    CALL solve_air_node_temperatures(nx, ny, pv_park, &
                                     T_nat, T_PVb, T_PVf, tair, &
                                     U_can, U_PVa, uref, &
                                     cdh_nat, cdh_PVb, cdh_can, &
                                     cdh_PVf, cdh_PVa, &
                                     L, a, width, &
                                     T_can, T_PVa)
    
    !-------------------------------------------------------------------
    ! 5. UNIFIED RK2 INTEGRATION FOR ALL NODES
    !-------------------------------------------------------------------
    
    ! STEP 1: Calculate forcing at t^n
    CALL compute_all_forcings(nx, ny, nz, nzsoil, pv_park, rnflx, &
                             L, a, beta, h_PV, h_clear, width, &
                             roufns, albedo_nat, albedo_PVf, albedo_PVb, &
                             emiss_nat, emiss_PVf, emiss_PVb, &
                             soiltyp(:,:,is), vegtyp, lai, veg, &
                             T_nat, T_soil, T_PVf, T_cell, T_PVb, &
                             T_can, T_PVa, &
                             qsoil_surf, qsoil_deep, wetcanp(:,:,is), snowdpth, qvsfc(:,:,is), &
                             uref, tair, qvair, radsw, radlwin, psfc, rhosfc, prcrate, &
                             U_can, U_PVa, uref, &
                             lambda_deg*pi/180., phi_deg*pi/180., S_diffuse, &
                             delta_PVf, delta_cell, delta_PVb, &
                             rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                             delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                             k_glass, k_EVA, k_backsheet, &
                             eta_STC, gamma_PV, &
                             zp, u, v, w, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                             cdh, cdh_nat, cdh_PVb, cdh_PVf, &
                             cdm, cdm_nat, cdm_PVb, cdm_PVf, &
                             cdq, cdq_nat, cdq_PVb, cdq_PVf, &
                             cdq_can, q_can, &
                             frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                             frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                             H_nat, LE_nat, G_nat, Rn_nat, Rn_ref, &
                             SW_out_sky_loc, LW_out_sky_loc, &
                             dum_qLW_PVf, dum_qLW_PVb, &
                             dum_H_PVf, dum_H_PVb, &
                             dum_qcond_PVf, dum_qcond_PVb, &
                             dum_qSW_PVf, dum_qSW_PVb, &
                             dum_qPV, &
                             dum_cdh_can, dum_cdh_PVa, &
                             dum_T_can, dum_T_PVa, &
                             dum_JLW_PVf, dum_JLW_PVb, &
                             dum_G_SW_nat, dum_J_SW_nat, &
                             dum_G_LW_nat, dum_J_LW_nat, &
                             q_can_loc, gcan_loc, vpd_f3_loc, vpd_ref_loc, .TRUE.)
    
    ! Predictor step: T^(n+1/2) = T^n + (dt/2) * frc^n
    ! Same predictor-corrector as soilebm (soilebm3d.f90)  [SYLVAIN]
    DO j = 1, ny-1
      DO i = 1, nx-1
      IF (pv_park(i,j).eq.1.) THEN
        T_PVf_half(i,j)  = T_PVf(i,j)  + 0.5*dtsfc * frc_PVf(i,j)
        T_cell_half(i,j) = T_cell(i,j) + 0.5*dtsfc * frc_cell(i,j)
        T_PVb_half(i,j)  = T_PVb(i,j)  + 0.5*dtsfc * frc_PVb(i,j)
       ENDIF
        T_nat_half(i,j)  = T_nat(i,j)  + 0.5*dtsfc * frc_nat(i,j)
        T_soil_half(i,j) = T_soil(i,j) + 0.5*dtsfc * frc_soil(i,j)
        ! Moisture half-step (clamped to [0, wsat] like ARPS soilebm)
        qsoil_surf_half(i,j) = qsoil_surf(i,j) + 0.5*dtsfc * frc_wsfc(i,j)
        qsoil_surf_half(i,j) = MAX(qsoil_surf_half(i,j), 0.0)
        qsoil_surf_half(i,j) = MIN(qsoil_surf_half(i,j), wsat(soiltyp(i,j,is)))
        qsoil_deep_half(i,j) = qsoil_deep(i,j) + 0.5*dtsfc * frc_wdp(i,j)
        qsoil_deep_half(i,j) = MAX(qsoil_deep_half(i,j), 0.0)
        qsoil_deep_half(i,j) = MIN(qsoil_deep_half(i,j), wsat(soiltyp(i,j,is)))
      END DO
    END DO
    wetcanp_half = wetcanp
    DO j = 1, ny-1
      DO i = 1, nx-1
        wetcanp_half(i,j,is) = wetcanp(i,j,is) + 0.5*dtsfc * frc_wcnp(i,j)
        wetcanp_half(i,j,is) = MAX(wetcanp_half(i,j,is), 0.0)
        wetcanp_half(i,j,is) = MIN(wetcanp_half(i,j,is), wrmax(i,j))
      END DO
    END DO
    
    ! Save forcing at t^n
    rhs_nat  = frc_nat
    rhs_soil = frc_soil
    rhs_PVf  = frc_PVf
    rhs_cell = frc_cell
    rhs_PVb  = frc_PVb
    rhs_wsfc = frc_wsfc
    rhs_wdp  = frc_wdp
    rhs_wcnp = frc_wcnp
    
    ! STEP 2: Update air node temps for half-step state, then calc forcing
    CALL solve_air_node_temperatures(nx, ny, pv_park, &
                                     T_nat_half, T_PVb_half, T_PVf_half, tair, &
                                     U_can, U_PVa, uref, &
                                     cdh_nat, cdh_PVb, cdh_can, &
                                     cdh_PVf, cdh_PVa, &
                                     L, a, width, &
                                     T_can, T_PVa)
    
    CALL compute_all_forcings(nx, ny, nz, nzsoil, pv_park, rnflx, &
                             L, a, beta, h_PV, h_clear, width, &
                             roufns, albedo_nat, albedo_PVf, albedo_PVb, &
                             emiss_nat, emiss_PVf, emiss_PVb, &
                             soiltyp(:,:,is), vegtyp, lai, veg, &
                             T_nat_half, T_soil_half, T_PVf_half, T_cell_half, T_PVb_half, &
                             T_can, T_PVa, &
                             qsoil_surf_half, qsoil_deep_half, wetcanp_half(:,:,is), snowdpth, qvsfc(:,:,is), &
                             uref, tair, qvair, radsw, radlwin, psfc, rhosfc, prcrate, &
                             U_can, U_PVa, uref, &
                             lambda_deg*pi/180., phi_deg*pi/180., S_diffuse, &
                             delta_PVf, delta_cell, delta_PVb, &
                             rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                             delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                             k_glass, k_EVA, k_backsheet, &
                             eta_STC, gamma_PV, &
                             zp, u, v, w, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                             cdh, cdh_nat, cdh_PVb, cdh_PVf, &
                             cdm, cdm_nat, cdm_PVb, cdm_PVf, &
                             cdq, cdq_nat, cdq_PVb, cdq_PVf, &
                             cdq_can, q_can, &
                             frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                             frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                             H_nat, LE_nat, G_nat, Rn_nat, Rn_ref, &
                             SW_out_sky_loc, LW_out_sky_loc, &
                             dum_qLW_PVf, dum_qLW_PVb, &
                             dum_H_PVf, dum_H_PVb, &
                             dum_qcond_PVf, dum_qcond_PVb, &
                             dum_qSW_PVf, dum_qSW_PVb, &
                             dum_qPV, &
                             dum_cdh_can, dum_cdh_PVa, &
                             dum_T_can, dum_T_PVa, &
                             dum_JLW_PVf, dum_JLW_PVb, &
                             dum_G_SW_nat, dum_J_SW_nat, &
                             dum_G_LW_nat, dum_J_LW_nat, &
                             q_can_loc, gcan_loc, vpd_f3_loc, vpd_ref_loc, .TRUE.)
    
    ! Canopy-air humidity of this small timestep (F3 of the next call) and
    ! Obukhov length from the fluxes just computed (lagged by one dtsfc,
    ! explicit u* <-> L feedback), as in PV_LAND_final:
    !   H_ref  = area-weighted fluxes PVa -> ref and canopy air -> ref
    !   LE_ref = latent heat flux of the natural surface
    ! with T_can, T_PVa of the half-step air-node solve.  [SYLVAIN]
    DO j = 1, ny-1
      DO i = 1, nx-1
      IF (pv_park(i,j).eq.1.) THEN
        q_can(i,j) = q_can_loc(i,j)
        tema = rhosfc(i,j) * cp * uref(i,j) *                            &
               ( cdh_PVa(i,j) * (T_PVa(i,j) - tair(i,j)) * width(i,j) / (width(i,j) + a(i,j)) &
               + cdh_can(i,j) * (T_can(i,j) - tair(i,j)) * a(i,j) / (width(i,j) + a(i,j)) )
        z_ref = zp(i,j,2) + 0.5*(zp(i,j,3)-zp(i,j,2))
        CALL obukhov_update(u_star_pv(i,j), tema, LE_nat(i,j),          &
                            tair(i,j), qvair(i,j), psfc(i,j), rhosfc(i,j), &
                            cp, lathv, p0, rddcp,                        &
                            MAX(z_ref - zd_park(i,j), 0.1), dtsfc, pv_tau_L, &
                            pv_L_floor, L_obk(i,j))
      ENDIF
      END DO
    END DO

    ! Average the two forcings
    DO j = 1, ny-1
      DO i = 1, nx-1
      IF (pv_park(i,j).eq.1.) THEN
        rhs_PVf(i,j)  = 0.5 * (rhs_PVf(i,j)  + frc_PVf(i,j))
        rhs_cell(i,j) = 0.5 * (rhs_cell(i,j) + frc_cell(i,j))
        rhs_PVb(i,j)  = 0.5 * (rhs_PVb(i,j)  + frc_PVb(i,j))
      ENDIF
        rhs_nat(i,j)  = 0.5 * (rhs_nat(i,j)  + frc_nat(i,j))
        rhs_soil(i,j) = 0.5 * (rhs_soil(i,j) + frc_soil(i,j))
        rhs_wsfc(i,j) = 0.5 * (rhs_wsfc(i,j) + frc_wsfc(i,j))
        rhs_wdp(i,j)  = 0.5 * (rhs_wdp(i,j)  + frc_wdp(i,j))
        rhs_wcnp(i,j) = 0.5 * (rhs_wcnp(i,j) + frc_wcnp(i,j))
      END DO
    END DO
    
    ! Corrector step: T^(n+1) = T^n + dt * avg(frc)
    DO j = 1, ny-1
      DO i = 1, nx-1
      IF (pv_park(i,j).eq.1.) THEN
        T_PVf(i,j)  = T_PVf(i,j)  + dtsfc * rhs_PVf(i,j)
        T_cell(i,j) = T_cell(i,j) + dtsfc * rhs_cell(i,j)
        T_PVb(i,j)  = T_PVb(i,j)  + dtsfc * rhs_PVb(i,j)
      ENDIF
        T_nat(i,j)  = T_nat(i,j)  + dtsfc * rhs_nat(i,j)
        T_soil(i,j) = T_soil(i,j) + dtsfc * rhs_soil(i,j)
        ! Moisture corrector (clamped to [0, wsat] like ARPS soilebm)
        qsoil_surf(i,j) = qsoil_surf(i,j) + dtsfc * rhs_wsfc(i,j)
        qsoil_surf(i,j) = MAX(qsoil_surf(i,j), 0.0)
        qsoil_surf(i,j) = MIN(qsoil_surf(i,j), wsat(soiltyp(i,j,is)))
        qsoil_deep(i,j) = qsoil_deep(i,j) + dtsfc * rhs_wdp(i,j)
        qsoil_deep(i,j) = MAX(qsoil_deep(i,j), 0.0)
        qsoil_deep(i,j) = MIN(qsoil_deep(i,j), wsat(soiltyp(i,j,is)))
        wetcanp(i,j,is) = wetcanp(i,j,is) + dtsfc * rhs_wcnp(i,j)
        wetcanp(i,j,is) = MAX(wetcanp(i,j,is), 0.0)
        wetcanp(i,j,is) = MIN(wetcanp(i,j,is), wrmax(i,j))
      END DO
    END DO

    END DO  ! it = 1, nsfcst

  ! Solve air node temps with final T_nat, T_PVb, T_PVf
  CALL solve_air_node_temperatures(nx, ny, pv_park, &
                                   T_nat, T_PVb, T_PVf, tair, &
                                   U_can, U_PVa, uref, &
                                   cdh_nat, cdh_PVb, cdh_can, &
                                   cdh_PVf, cdh_PVa, &
                                   L, a, width, &
                                   T_can, T_PVa)
  
  ! Recompute all forcings at final state (we only use the flux outputs)
  CALL compute_all_forcings(nx, ny, nz, nzsoil, pv_park, rnflx, &
                           L, a, beta, h_PV, h_clear, width, &
                           roufns, albedo_nat, albedo_PVf, albedo_PVb, &
                           emiss_nat, emiss_PVf, emiss_PVb, &
                           soiltyp(:,:,is), vegtyp, lai, veg, &
                           T_nat, T_soil, T_PVf, T_cell, T_PVb, &
                           T_can, T_PVa, &
                           qsoil_surf, qsoil_deep, wetcanp(:,:,is), snowdpth, qvsfc(:,:,is), &
                           uref, tair, qvair, radsw, radlwin, psfc, rhosfc, prcrate, &
                           U_can, U_PVa, uref, &
                           lambda_deg*pi/180., phi_deg*pi/180., S_diffuse, &
                           delta_PVf, delta_cell, delta_PVb, &
                           rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                           delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                           k_glass, k_EVA, k_backsheet, &
                           eta_STC, gamma_PV, &
                           zp, u, v, w, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                           cdh, cdh_nat, cdh_PVb, cdh_PVf, &
                           cdm, cdm_nat, cdm_PVb, cdm_PVf, &
                           cdq, cdq_nat, cdq_PVb, cdq_PVf, &
                             cdq_can, q_can, &
                           frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                           frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                           H_nat, LE_nat, G_nat, Rn_nat, Rn_ref, &
                           SW_out_sky, LW_out_sky, &
                             dum_qLW_PVf, dum_qLW_PVb, &
                             dum_H_PVf, dum_H_PVb, &
                             dum_qcond_PVf, dum_qcond_PVb, &
                             dum_qSW_PVf, dum_qSW_PVb, &
                             dum_qPV, &
                             dum_cdh_can, dum_cdh_PVa, &
                             dum_T_can, dum_T_PVa, &
                             dum_JLW_PVf, dum_JLW_PVb, &
                             dum_G_SW_nat, dum_J_SW_nat, &
                             dum_G_LW_nat, dum_J_LW_nat, &
                             q_can_loc, gcan_loc, vpd_f3_loc, vpd_ref_loc, diag_out)
  ! frc_* are discarded; H_nat, LE_nat, G_nat, Rn_ref_out, SW/LW_out_sky matter.
  ! With diag_out = .FALSE. only the soil part is computed (LE_nat, H_nat,
  ! G_nat for the fluxes); the radiation and PV diagnostics keep the
  ! values of the last output step.  [SYLVAIN]
  IF (diag_out) THEN
  
  qLW_PVf = dum_qLW_PVf
  qLW_PVb = dum_qLW_PVb
  H_PVf = dum_H_PVf
  H_PVb = dum_H_PVb
  qcond_PVf = dum_qcond_PVf
  qcond_PVb = dum_qcond_PVb
  qSW_PVf = dum_qSW_PVf
  qSW_PVb = dum_qSW_PVb
  qPV = dum_qPV
  !T_can = dum_T_can !modif sylvain
  !T_PVa = dum_T_PVa !modif sylvain
  JLW_PVf = dum_JLW_PVf
  JLW_PVb = dum_JLW_PVb
  G_SW_nat = dum_G_SW_nat
  J_SW_nat = dum_J_SW_nat
  G_LW_nat = dum_G_LW_nat
  J_LW_nat = dum_J_LW_nat
  ENDIF  ! diag_out
  
  
  !=====================================================================
  ! COMPUTE OUTPUT FLUXES (from final-state diagnostics)
  !=====================================================================
  
  ! H_ref = weighted combination of H_PVa_ref and H_can_ref
  DO j = 1, ny-1
    DO i = 1, nx-1
    IF (PV_PARK(i,j).eq.1.) THEN
      H_PVa_ref(i,j) = rhosfc(i,j) * cp * cdh_PVa(i,j) * uref(i,j) * &
                  (T_PVa(i,j) - tair(i,j))
      H_can_ref(i,j) = rhosfc(i,j) * cp * cdh_can(i,j) * uref(i,j) * &
                  (T_can(i,j) - tair(i,j))
      
      H_ref(i,j) = H_PVa_ref(i,j) * width(i,j) / (width(i,j) + a(i,j)) + &
					H_can_ref(i,j) * a(i,j) / (width(i,j) + a(i,j))
	ELSE
      H_ref(i,j) = H_nat(i,j)
	ENDIF
    END DO
  END DO
  
      ! Fluxes to the atmosphere, accumulated over the soil types:
      !   PV park : from the end-of-step state (H_ref, LE_nat)
      !   outside : at t^n (ptsflx_n, qvsflx_n), as sfcflx with sfcphy = 4
      ! ptsfc_eff, qvsfc_eff: effective surface values giving the
      ! end-of-step fluxes; qvsfc: effective surface humidity (as soilebm)
      DO j=1,ny-1
        DO i=1,nx-1

        pts_type = - (H_ref(i,j)/cp)*(p0/psfc(i,j))**rddcp
        qvs_type = - (LE_nat(i,j)/lathv)
        tema = rhosfc(i,j)*MAX(uref(i,j),vsfcmin)
        ptsfc_eff(i,j) = pt1(i,j)-pts_type/(cdh(i,j)*tema)
        qvsfc_eff(i,j) = qv(i,j,2)-qvs_type/(cdq(i,j)*tema)

        IF (pv_park(i,j).eq.1.) THEN
          ptsflx(i,j) = ptsflx(i,j) + pts_type*stypfrct(i,j,is)
          qvsflx(i,j) = qvsflx(i,j) + qvs_type*stypfrct(i,j,is)
          qvsfc(i,j,is) = qvsfc_eff(i,j)
        ELSE
          ptsflx(i,j) = ptsflx(i,j) + ptsflx_n(i,j)*stypfrct(i,j,is)
          qvsflx(i,j) = qvsflx(i,j) + qvsflx_n(i,j)*stypfrct(i,j,is)
          qvsfc(i,j,is) = qvair(i,j)                                    &
                        + LE_nat(i,j)/(lathv*rhosfc(i,j)*cdq(i,j)*uref(i,j))
        ENDIF
 	  
        END DO
      END DO
  
        DO i=1,nx-1
        DO j=1,ny-1
    tsoil(i,j,1,is) = T_nat(i,j)
    tsoil(i,j,2,is) = T_soil(i,j)
    qsoil(i,j,1,is) = qsoil_surf(i,j)
    qsoil(i,j,2,is) = qsoil_deep(i,j)
        END DO
      END DO
 
!-----------------------------------------------------------------------
!  If we want to use the sensible and latent heat fluxes calculated
!  inside the soil model (soilebm) in the atmospheric model, doing the
!  folloing do loop. Otherwise, the flux calculated in sfcflx will be
!  used, which is based on qvsfc one time step earlier. Otherwise,
!  they should be the same.
!-----------------------------------------------------------------------
!
!     DO i=1,nx-1
!       DO j=1,ny-1
!         ptsflx(i,j) = ptsflx(i,j) + shflx(i,j)*stypfrct(i,j,is)
!         qvsflx(i,j) = qvsflx(i,j) + lhflx(i,j)*stypfrct(i,j,is)
!       END DO
!     END DO

      DO k=1,nzsoil
        DO j=1,ny-1
          DO i=1,nx-1
            tsoil  (i,j,k,0) = tsoil  (i,j,k,0)                         &
                    + tsoil  (i,j,k,is)*stypfrct(i,j,is)
            qsoil  (i,j,k,0) = qsoil  (i,j,k,0)                         &
                    + qsoil  (i,j,k,is)*stypfrct(i,j,is)
          END DO
        END DO
      END DO

      DO j=1,ny-1
        DO i=1,nx-1
          wetcanp(i,j,0) = wetcanp(i,j,0)                               &
                         + wetcanp(i,j,is)*stypfrct(i,j,is)
          qvsfc  (i,j,0) = qvsfc  (i,j,0)                               &
                         + qvsfc  (i,j,is)*stypfrct(i,j,is)
        END DO
      END DO

    END DO  ! loop over is

  RETURN
END SUBROUTINE pvland_sfcphy
