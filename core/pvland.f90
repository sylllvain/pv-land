!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INIT_PVLAND                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE INIT_PVLAND(nx,ny,nz,ptbar,ptprt,pbar,pprt,qvsfc,tsoil,&
			tem1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  
!
!-----------------------------------------------------------------------
!
!  AUTHOR: 
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  USE module_pvland_in
  USE module_pvland_out
  USE rad_pvland
  USE pvland_rad_cache
  USE pvland_nml
  
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  REAL :: ptbar(nx,ny,nz),ptprt(nx,ny,nz),pbar(nx,ny,nz),pprt(nx,ny,nz),&
	tsoil(nx,ny),qvsfc(nx,ny),tem1(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: beta_deg,theta_0_deg,psfc
  REAL :: t_air                ! temperature of the first model level (K)
  REAL, PARAMETER :: pi = 3.14159265359
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ALLOCATE(szenith_deg(nx,ny), sazimuth_deg(nx,ny))

  ALLOCATE(pv_park(nx,ny))
  
  ALLOCATE(a(nx,ny), beta(nx,ny), L(nx,ny), width(nx,ny), z_PVa(nx,ny),&
  h_PV(nx,ny), h_clear(nx,ny), lambda_f_park(nx,ny), l_nat_can(nx,ny),&
  theta_0(nx,ny), albedo_nat(nx,ny),&
  albedo_PVf(nx,ny), albedo_PVb(nx,ny), &
  emiss_nat(nx,ny), emiss_PVf(nx,ny), emiss_PVb(nx,ny), z_can(nx,ny))

  ALLOCATE(delta_PVf(nx,ny), delta_cell(nx,ny), delta_PVb(nx,ny), &
   rho_cp_PVf(nx,ny), rho_cp_cell(nx,ny), rho_cp_PVb(nx,ny), &
   k_glass(nx,ny), k_EVA(nx,ny), k_backsheet(nx,ny), &
   delta_glass(nx,ny), delta_EVA_top(nx,ny), &
   delta_EVA_bot(nx,ny), delta_backsheet(nx,ny), &
   eta_STC(nx,ny), gamma_PV(nx,ny))

  ! T_can, T_PVa: allocated in arps.f90 (restart state)
  ALLOCATE(U_can(nx,ny), U_PVa(nx,ny))

  ALLOCATE(H_ref(nx,ny), Rn_ref(nx,ny))
  ALLOCATE(H_nat(nx,ny), LE_nat(nx,ny), G_nat(nx,ny), Rn_nat(nx,ny))
  ALLOCATE(SW_out_sky(nx,ny), LW_out_sky(nx,ny))
  
  ALLOCATE(H_PVa_ref(nx,ny), H_can_ref(nx,ny))
  ALLOCATE(H_PVb(nx,ny), H_PVf(nx,ny))
  ALLOCATE(G_SW_nat(nx,ny), J_SW_nat(nx,ny))
  ALLOCATE(G_LW_nat(nx,ny), J_LW_nat(nx,ny))

  ALLOCATE(qLW_PVf(nx,ny), qLW_PVb(nx,ny))
  ALLOCATE(qcond_PVf(nx,ny), qcond_PVb(nx,ny),&
	qSW_PVf(nx,ny), qSW_PVb(nx,ny), qPV(nx,ny))
  ALLOCATE(JLW_PVf(nx,ny), JLW_PVb(nx,ny))

  ALLOCATE(z0_park(nx,ny), zd_park(nx,ny))
  ALLOCATE(AZ_PARK(nx,ny))

  ! ====== PV park location ======
  ! Host-model specific (pvland_park_map.f90): ARPS grid + MPI in ARPS,
  ! whole domain in the offline model.  [SYLVAIN]
  CALL pvland_park_map(nx,ny,pv_park,tem1)
  
  ! ====== PV park parameters: namelist &pvland (module pvland_nml) ======
  ! Read by initpara (ARPS, arps.input) or by the offline main program
  ! (pvland.nml) before this call.  [SYLVAIN]

  ! Park geometry
  L = pv_L                 ! Panel length along the slope (m)
  a = pv_a                 ! Inter-row spacing (m)
  beta_deg = pv_beta_deg   ! Tilt angle (degrees)
  h_PV = pv_h_PV           ! Height of the panel top edge (m)
  theta_0_deg = pv_theta_0_deg ! Park orientation angle (degrees)
  AZ_PARK  = pv_az_park    ! Meteorological azimuth faced by the panel fronts (deg)

  
  ! Optical properties
  albedo_nat = pv_albedo_nat  ! Natural surface albedo
  albedo_PVf = pv_albedo_PVf  ! PV front albedo
  albedo_PVb = pv_albedo_PVb  ! PV back albedo
  emiss_nat = pv_emiss_nat    ! Natural surface emissivity
  emiss_PVf = pv_emiss_PVf    ! PV front emissivity
  emiss_PVb = pv_emiss_PVb    ! PV back emissivity
  
  ! Geometry: tilt and orientation. The quantities derived from the tilt
  ! (h_clear, width, lambda_f_park, z_can, z_PVa, view factors) are
  ! computed in pvland_geometry, called below.
  ! written as deg*(pi/180.) to give the same rounding as the former
  ! constant values folded by the compiler (24*pi/180 differs by 1 ulp)
  beta = beta_deg * (pi / 180.0)
  theta_0 = theta_0_deg * (pi / 180.0)

  ! (l_nat_can = z_can is set in pvland_geometry)
  
  ! PV module thermal and electrical parameters
  delta_PVf       = pv_delta_PVf
  delta_cell      = pv_delta_cell
  delta_PVb       = pv_delta_PVb
  rho_cp_PVf      = pv_rho_cp_PVf
  rho_cp_cell     = pv_rho_cp_cell
  rho_cp_PVb      = pv_rho_cp_PVb
  k_glass         = pv_k_glass
  k_EVA           = pv_k_EVA
  k_backsheet     = pv_k_backsheet
  delta_glass     = pv_delta_glass
  delta_EVA_top   = pv_delta_EVA_top
  delta_EVA_bot   = pv_delta_EVA_bot
  delta_backsheet = pv_delta_backsheet
  eta_STC         = pv_eta_STC
  gamma_PV        = pv_gamma_PV

  ! ====== Radiosity model: temperature-independent quantities ======
  ! View factors and LU factorisations of the SW and LW radiosity
  ! matrices depend on geometry and optical properties only.
  !   vf_update_opt = 0 : computed once here (fixed panels)
  !   vf_update_opt = 1 : recomputed at every PV-LAND time step
  !                       (pvland_sfcphy), for moving panels, e.g. sun
  !                       tracking: beta may then change in time.  [SYLVAIN]

  ALLOCATE(F_pv(nsurf_pv,nsurf_pv,nx,ny), Fsky_pv(nsurf_pv,nx,ny))
  ALLOCATE(LU_SW_pv(nsurf_pv,nsurf_pv,nx,ny), LU_LW_pv(nsurf_pv,nsurf_pv,nx,ny))
  ALLOCATE(ipiv_SW_pv(nsurf_pv,nx,ny), ipiv_LW_pv(nsurf_pv,nx,ny))
  ALLOCATE(sw_valid(nx,ny), sw_key(4,nx,ny))
  ALLOCATE(J_SW_pv(nsurf_pv,nx,ny), G_SW_pv(nsurf_pv,nx,ny), q_SW_pv(nsurf_pv,nx,ny))

  CALL pvland_geometry(nx,ny)
                        
    ! Initial PV-LAND state (initopt = 1). With a restart (initopt /= 1)
    ! these variables, allocated in arps.f90, are read from the restart
    ! file (rst3d.f90).
    !   - air nodes (T_can, T_PVa) and PV nodes (T_PVf, T_cell, T_PVb):
    !     temperature of the first model level
    !   - effective surface values ptsfc_eff, qvsfc_eff of the soil
    !   - neutral Obukhov length; canopy-air humidity set at the first
    !     call of pvland_sfcphy (q_can < 0)
    IF (initopt.eq.1) THEN
    L_obk(:,:) = 0.0
    q_can(:,:) = -1.0
    DO j=1,ny-1
      DO i=1,nx-1
        psfc = 0.5 * ( pprt(i,j,1) + pbar(i,j,1)                         &
                     + pprt(i,j,2) + pbar(i,j,2) )
        t_air = (ptbar(i,j,2)+ptprt(i,j,2)) * (psfc/p0)**rddcp
        T_can(i,j)  = t_air
        T_PVa(i,j)  = t_air
        T_PVf(i,j)  = t_air
        T_cell(i,j) = t_air
        T_PVb(i,j)  = t_air
        ptsfc_eff(i,j) = tsoil(i,j)*(p0/psfc)**rddcp
        qvsfc_eff(i,j) = qvsfc(i,j)
      END DO
    END DO
    ENDIF
    
  RETURN
END SUBROUTINE init_pvland

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PVLAND_GEOMETRY            ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE pvland_geometry(nx,ny)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Quantities of the PV park that depend on the panel tilt beta:
!    - h_clear, width, lambda_f_park, z_can, l_nat_can, z_PVa = (h_PV+h_clear)/2
!      (h_PV, the height of the panel top edge, is kept fixed);
!    - view factors F_pv, sky view factors Fsky_pv and LU factorisations
!      of the SW and LW radiosity matrices (module pvland_rad_cache).
!  The cached SW solution is invalidated.
!
!  Called once by INIT_PVLAND, and at every PV-LAND time step by
!  pvland_sfcphy when vf_update_opt = 1 (moving panels).
!
!  Surfaces: 1=nat_n, 2=PVf_n, 3=PVb_n, 4=PVb_{n-1}, 5=PVf_{n-1},
!            6=PVb_{n+1}, 7=nat_{n-1}, 8=nat_{n+1}
!
!  AUTHOR: moved from INIT_PVLAND  [SYLVAIN]
!
!-----------------------------------------------------------------------
!
  USE module_pvland_in
  USE pvland_rad_cache

  IMPLICIT NONE

  INTEGER :: nx,ny

  INTEGER :: i, j
  REAL :: coef(8)              ! albedo or 1-emissivity of the 8 surfaces
  INTEGER :: info_rad

  h_clear = h_PV - L * SIN(beta)
  width = L * COS(beta)
  
  ! Park-scale roughness: frontal area density (geometry-only)
  DO j=1,ny
  DO i=1,nx
  CALL compute_lambda_f(a(i,j), L(i,j), beta(i,j), lambda_f_park(i,j))
  ENDDO
  ENDDO
  
  ! Heights for wind profile
  z_can = (h_PV + h_clear) / 4.0

  ! Eq. 6.3-6.5: natural ground to canopy air
  l_nat_can  = z_can

  ! air above the panels: mid-height of the panels (version 2)
  z_PVa = (h_PV + h_clear) / 2.0

  ! View factors and LU factorisations of the radiosity matrices
  DO j=1,ny-1
  DO i=1,nx-1
    CALL compute_view_factor_matrix(L(i,j), a(i,j), width(i,j), h_PV(i,j), beta(i,j), &
                                    F_pv(:,:,i,j))
    CALL compute_sky_view_factors(F_pv(:,:,i,j), Fsky_pv(:,i,j))

    coef = (/ albedo_nat(i,j), albedo_PVf(i,j), albedo_PVb(i,j), albedo_PVb(i,j), &
              albedo_PVf(i,j), albedo_PVb(i,j), albedo_nat(i,j), albedo_nat(i,j) /)
    CALL factor_radiosity_matrix(nsurf_pv, coef, F_pv(:,:,i,j), &
                                 LU_SW_pv(:,:,i,j), ipiv_SW_pv(:,i,j), info_rad)

    coef = 1.0 - (/ emiss_nat(i,j), emiss_PVf(i,j), emiss_PVb(i,j), emiss_PVb(i,j), &
                    emiss_PVf(i,j), emiss_PVb(i,j), emiss_nat(i,j), emiss_nat(i,j) /)
    CALL factor_radiosity_matrix(nsurf_pv, coef, F_pv(:,:,i,j), &
                                 LU_LW_pv(:,:,i,j), ipiv_LW_pv(:,i,j), info_rad)
  ENDDO
  ENDDO

  ! the SW solution depends on the geometry: recompute it at the next call
  sw_valid(:,:) = .FALSE.

  RETURN
END SUBROUTINE pvland_geometry

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DRAG_PVLAND                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE DRAG_PVLAND(nx,ny,nz,zp,u,v,w,pbar,pt1,ptsfc,L_obk,&
		T_nat,rhosfc,z0_nat,soiltyp,uref,qvref,&
		cdh_nat,cdq_nat,cdm_nat,&
		cdh_PVa,cdq_PVa,cdm_PVa,&
		cdh_PVb,cdq_PVb,cdm_PVb,&
		cdh_PVf,cdq_PVf,cdm_PVf,&
		cdh_can,cdq_can,cdm_can,&
		z0_park_arr,u_star)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Wind profile in and above the PV park and exchange coefficients of
!  the PV-LAND interfaces (version 2, from PV_LAND_final):
!    - park z0 and d from the wind direction, the frontal density and the
!      roughness z0_nat / displacement zd_nat of the natural vegetation;
!    - U_can, U_PVa from the stability-dependent profile (u* at z_ref),
!      floored at vsfcmin;
!    - natural ground -> canopy air, canopy air -> reference, air above
!      the panels -> reference: dragcoef_pv (thermal roughness zt/z0 =
!      pv_ztz0_nat for the ground, pv_ztz0_air for the air nodes);
!    - panel faces: mixed convection, forced h = C_for*U**p_for and
!      free convection (Kaplani & Kaplanis 2014) scaled by pv_fac_h_free.
!  Outside the PV park only uref (and z0_park_arr = z0_nat) are used.
!
!-----------------------------------------------------------------------
!
  USE module_pvland_in
  USE module_pvland_out, L_obk_state => L_obk   ! L_obk is the dummy argument here
  USE pvland_nml, ONLY : pv_zd_nat, pv_ztz0_nat, pv_ztz0_air,          &
                         pv_C_for_f, pv_C_for_b, pv_p_for_f, pv_p_for_b, &
                         pv_fac_h_free
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz       ! The number grid points in 3 directions
  REAL, INTENT(IN) :: zp     (nx,ny,nz) ! The height of the terrain.

  REAL, INTENT(IN) :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL, INTENT(IN) :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL, INTENT(IN) :: w      (nx,ny,nz) ! Total w-velocity (m/s)
  REAL, INTENT(IN) :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL, INTENT(IN) :: pt1  (nx,ny) ! Potential temperature (K) at reference level
  REAL, INTENT(IN) :: ptsfc  (nx,ny) ! Potential surface temperature (K)
  REAL, INTENT(IN) :: qvref      (nx,ny) ! ref qvair (kg/kg)
  REAL, INTENT(IN) :: T_nat      (nx,ny) ! soil surface temperature (K)
  REAL, INTENT(IN) :: rhosfc (nx,ny)    ! Base-state air density at the surface
  REAL, INTENT(IN) :: z0_nat(nx,ny)     ! Roughness length of the natural surface (m)
  INTEGER, INTENT(IN) :: soiltyp      (nx,ny) ! 
  REAL, INTENT(IN) :: L_obk(nx,ny)      ! Obukhov length (m), 0 = neutral
  
  REAL, INTENT(OUT) :: uref      (nx,ny) ! ref velocity (m/s)
  REAL, INTENT(OUT) :: cdh_nat(nx,ny),cdq_nat(nx,ny),cdm_nat(nx,ny),&
		cdh_PVa(nx,ny),cdq_PVa(nx,ny),cdm_PVa(nx,ny),&
		cdh_PVb(nx,ny),cdq_PVb(nx,ny),cdm_PVb(nx,ny),&
		cdh_PVf(nx,ny),cdq_PVf(nx,ny),cdm_PVf(nx,ny),&
		cdh_can(nx,ny),cdq_can(nx,ny),cdm_can(nx,ny),&
		z0_park_arr(nx,ny)
  REAL, INTENT(OUT) :: u_star(nx,ny)    ! friction velocity at z_ref (m/s), PV park only
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL, PARAMETER :: pi = 3.14159265359
  REAL :: tema,temb,z_ref
  REAL :: Uplus
  REAL, PARAMETER :: Pr = 0.71
  
  ! (SAVE: allocated at the first call only)  [SYLVAIN]
  REAL, ALLOCATABLE, SAVE :: udir(:,:),zp_dc(:,:,:),                     &
	wdspsfc(:,:), zp_dc_nat(:,:,:), pt_nat(:,:),                          &
	pt_can(:,:), pt_PVa(:,:), z0_sfc_arr(:,:)
	
  ! Dummy soil type for dragcoef_pv
  INTEGER, ALLOCATABLE, SAVE :: soiltyp_dummy(:,:)

  !===================================================================
  ! PV panel convection (both faces):
  !   forced: h_forced = C_for * U**p_for (pv_C_for_*, pv_p_for_*)
  !   free  : Kaplani & Kaplanis (2014), scaled by pv_fac_h_free
  !   mixing: Kaplani Eq.(20), m = 3
  ! Air properties at T_air + 0.75*(T_surface - T_air).
  !===================================================================
  REAL :: Tprop      ! [K]      property temperature, T_air + 0.75*(T_surface - T_air)
  REAL :: kair       ! [W/m/K]  air thermal conductivity at Tprop
  REAL :: muair      ! [Pa.s]   air dynamic viscosity at Tprop (Sutherland)
  REAL :: nuair      ! [m2/s]   air kinematic viscosity, mu/rho
  REAL :: alair      ! [m2/s]   air thermal diffusivity, k/(rho*cp)
  REAL :: dTs        ! [K]      signed surface-to-air temperature difference
  REAL :: vface      ! [m/s]    wind speed at the face (U_can or U_PVa)
  REAL :: h_cv       ! [W/m2/K] convective heat transfer coefficient
  REAL :: cc_pr      ! [-]      Prandtl function of Kaplani Eq.(8), precomputed once
  REAL :: Tf_na      ! [K]      Kaplani expansion temperature, Ta + 0.25*(Ts-Ta)
  REAL :: Ra_na      ! [-]      Rayleigh number on the physical length L (Kaplani Eq.(6))
  REAL :: theta_na   ! [rad]    plate angle from the vertical (Kaplani convention)
  REAL :: Grc_na     ! [-]      critical Grashof number (Kaplani, laminar-turbulent transition)
  REAL :: Rac_na     ! [-]      critical Rayleigh number, Grc*Pr
  REAL :: Nu_nat     ! [-]      natural-convection Nusselt (Kaplani, on L)
  REAL :: h_free     ! [W/m2/K] free (natural) convection coefficient
  REAL :: h_forced   ! [W/m2/K] forced-convection coefficient
  LOGICAL :: unstab  ! [-]      .TRUE. if the buoyancy arrangement is unstable

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'sfcphycst.inc'
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! Prandtl function of the Kaplani Eq.(8) branch
  cc_pr = ( 1.0 + (0.492/Pr)**(9.0/16.0) )**(8.0/27.0)

  IF (.NOT. ALLOCATED(soiltyp_dummy)) THEN
  ALLOCATE(soiltyp_dummy(nx,ny))
  soiltyp_dummy = 11

    ALLOCATE (udir(nx,ny),zp_dc(nx,ny,3),                               &
		wdspsfc(nx,ny),zp_dc_nat(nx,ny,nz),                             &
		pt_nat(nx,ny),pt_can(nx,ny),pt_PVa(nx,ny),z0_sfc_arr(nx,ny))
  END IF
    
    DO j=1,ny-1
      DO i=1,nx-1
    !-------------------------------------------------------------------
    ! 1b. COMPUTE PARK-SCALE z0 AND d FROM WIND DIRECTION
    !-------------------------------------------------------------------
    IF (PV_PARK(i,j).eq.1.) THEN
      udir(i,j)=ATAN2((u(i,j,2)+u(i+1,j,2)),(v(i,j,2)+v(i,j+1,2)))
      if (udir(i,j) < 0.0) udir(i,j) = udir(i,j) + 2*pi
    CALL z0_theta(udir(i,j), lambda_f_park(i,j), h_PV(i,j), beta(i,j)*180./pi, theta_0(i,j), &
                  z0_nat(i,j), pv_zd_nat, z0_park(i,j))
    CALL d_theta(udir(i,j), lambda_f_park(i,j), h_PV(i,j), theta_0(i,j), pv_zd_nat, zd_park(i,j))
    z0_park_arr(i,j) = z0_park(i,j)
    ELSE
    z0_park_arr(i,j) = z0_nat(i,j)
    udir(i,j) = 0.            ! not used outside the park (defined value)
    ENDIF

    !-------------------------------------------------------------------
    ! 2. WIND PROFILE
    !-------------------------------------------------------------------
      tema = 0.5*(u(i+1,j,2)+u(i,j,2))
      temb = 0.5*(v(i,j+1,2)+v(i,j,2))
      uref(i,j)=MAX(vsfcmin,                                         &
                   SQRT(tema*tema+temb*temb+w(i,j,2)*w(i,j,2)))
	  enddo !i
	  enddo !j
		       
    DO j=1,ny-1
      DO i=1,nx-1
    IF (PV_PARK(i,j).eq.1.) THEN
    z_ref = zp(i,j,2) + 0.5*(zp(i,j,3)-zp(i,j,2))
    
    CALL compute_U_plus(z_ref, udir(i,j), a(i,j), beta(i,j)*180./pi, &
		L(i,j), h_clear(i,j), theta_0(i,j), z0_nat(i,j), pv_zd_nat, L_obk(i,j), Uplus)
    u_star(i,j) = uref(i,j) / Uplus
    
    CALL compute_U_plus(z_can(i,j), udir(i,j), a(i,j), beta(i,j)*180./pi, &
		L(i,j), h_clear(i,j), theta_0(i,j), z0_nat(i,j), pv_zd_nat, L_obk(i,j), Uplus)
    U_can(i,j) = MAX(u_star(i,j) * Uplus, vsfcmin)   ! same wind floor as in sfcflx
    
    CALL compute_U_plus(z_PVa(i,j), udir(i,j), a(i,j), beta(i,j)*180./pi, &
		L(i,j), h_clear(i,j), theta_0(i,j), z0_nat(i,j), pv_zd_nat, L_obk(i,j), Uplus)
    U_PVa(i,j) = MAX(u_star(i,j) * Uplus, vsfcmin)
    ELSE
    u_star(i,j) = 0.
    U_can(i,j)  = uref(i,j)
    U_PVa(i,j)  = uref(i,j)
    ENDIF

      END DO
    END DO

    !-------------------------------------------------------------------
    ! 3. DRAG COEFFICIENTS - computed with current temperatures
    !-------------------------------------------------------------------    

    zp_dc_nat(:,:,1) = zp(:,:,1)
    zp_dc_nat(:,:,2) = 0.0
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt_nat(i,j) = T_nat(i,j) * (p0/pbar(i,j,2))**rddcp
        pt_can(i,j) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
      END DO
    END DO

    ! --- natural ground -> canopy air: interface l_nat_can, wind U_can ---
    z0_sfc_arr(:,:)  = z0_nat
    zp_dc_nat(:,:,3) = 2.0*l_nat_can        ! => z1 = l_nat_can
    CALL dragcoef_pv(nx, ny, nz, zp_dc_nat, pt_nat, qvref, soiltyp_dummy, &
                  z0_sfc_arr, pv_ztz0_nat, U_can, pt_can, cdm_nat, cdh_nat, cdq_nat)

    ! --- cdh_PVb: mixed convection on the BACK face (faces the ground) ---
    cdm_PVb(:,:) = 0.0   ! momentum not used for the PV faces
    DO j = 1, ny-1
      DO i = 1, nx-1
        ! signed surface-to-air temperature difference  [K]
        dTs    = T_PVb(i,j) - T_can(i,j)
        ! property temperature, 3/4 of the way from the air to the face
        Tprop  = T_can(i,j) + 0.75*dTs
        ! air thermal conductivity k(Tprop)  [W/m/K]
        kair   = 1.5207e-11*Tprop**3 - 4.8574e-8*Tprop**2 + 1.0184e-4*Tprop - 3.9333e-4
        ! air dynamic viscosity mu(Tprop), Sutherland's law  [Pa.s]
        muair  = 1.458e-6 * Tprop**1.5 / (Tprop + 110.4)
        ! air kinematic viscosity nu = mu/rho  [m2/s]
        nuair  = muair / rhosfc(i,j)
        ! air thermal diffusivity alpha = k/(rho*cp)  [m2/s]
        alair  = kair / (rhosfc(i,j) * cp)
        ! buoyancy arrangement: back face points DOWN, so a COLD face (night) is unstable
        unstab = (dTs < 0.0)
        ! forced convection, canopy wind
        vface    = U_can(i,j)
        h_forced = pv_C_for_b * vface**pv_p_for_b
        ! ---- Natural convection: Kaplani & Kaplanis (2014), Sec. 2.1 ----
        Tf_na    = T_can(i,j) + 0.25*dTs
        Ra_na    = 9.81 / Tf_na * MAX(ABS(dTs), 0.05) * L(i,j)**3 / (nuair*alair)
        theta_na = pi/2.0 - beta(i,j)         ! angle from the VERTICAL (Kaplani convention)
        IF (unstab) THEN
          ! up-facing hot-plate arrangement (day front / night back): Kaplani Eq.(7)
          IF (beta(i,j) > pi/6.0) THEN        ! theta < 60 deg
            Grc_na = 1.327e10 * EXP(-3.708*theta_na)
            Rac_na = Grc_na * Pr
            Nu_nat = 0.13*MAX(Ra_na**(1.0/3.0) - Rac_na**(1.0/3.0), 0.0)  &
                   + 0.56*( Rac_na*MAX(SIN(beta(i,j)),0.0) )**0.25
          ELSE                           ! theta >= 60 deg (near horizontal, turbulent plume)
            Nu_nat = 0.13 * Ra_na**(1.0/3.0)
          END IF
        ELSE
          ! down-facing hot-plate arrangement (day back / night front): Kaplani Eq.(8)-(9)
          IF (beta(i,j) >= pi/6.0) THEN       ! beta >= 30 deg : Kaplani Eq.(8), all Ra
            Nu_nat = ( 0.825 + 0.387*( Ra_na*MAX(SIN(beta(i,j)),0.0) )**(1.0/6.0) / cc_pr )**2
          ELSE IF (beta(i,j) > 2.0*pi/180.0) THEN  ! 2 < beta < 30 deg : Kaplani Eq.(9), inclined
            Nu_nat = 0.56 * ( Ra_na*MAX(SIN(beta(i,j)),0.0) )**0.25
          ELSE                           ! beta <= 2 deg : Kaplani Eq.(9), horizontal plate
            Nu_nat = 0.58 * Ra_na**0.2
          END IF
        END IF
        Nu_nat = MAX(Nu_nat, 0.0)
        h_free   = pv_fac_h_free * Nu_nat * kair / L(i,j)
        h_cv     = ( h_forced**3 + h_free**3 )**(1.0/3.0)
        ! convert h to the model's exchange coefficient: H = rho*cp*cdh*U*dT
        cdh_PVb(i,j) = h_cv / ( rhosfc(i,j) * cp * U_can(i,j))
      END DO
    END DO
    ! moisture coefficient follows heat (dry surface: no latent flux anyway)
    cdq_PVb(:,:) = cdh_PVb(:,:)
    
    ! --- cdh_PVf : mixed convection on the FRONT face (faces the sky) ---
    DO j = 1, ny-1
      DO i = 1, nx-1
        pt_PVa(i,j) = T_PVa(i,j) * (p0/pbar(i,j,2))**rddcp
      END DO
    END DO
    cdm_PVf(:,:) = 0.0   ! momentum not used for the PV faces
    DO j = 1, ny-1
      DO i = 1, nx-1
        dTs    = T_PVf(i,j) - T_PVa(i,j)
        Tprop  = T_PVa(i,j) + 0.75*dTs
        kair   = 1.5207e-11*Tprop**3 - 4.8574e-8*Tprop**2 + 1.0184e-4*Tprop - 3.9333e-4
        muair  = 1.458e-6 * Tprop**1.5 / (Tprop + 110.4)
        nuair  = muair / rhosfc(i,j)
        alair  = kair / (rhosfc(i,j) * cp)
        ! buoyancy arrangement: front face points UP, so a HOT face (day) is unstable
        unstab = (dTs > 0.0)
        ! forced convection, wind above the panels
        vface    = U_PVa(i,j)
        h_forced = pv_C_for_f * vface**pv_p_for_f
        Tf_na    = T_PVa(i,j) + 0.25*dTs
        Ra_na    = 9.81 / Tf_na * MAX(ABS(dTs), 0.05) * L(i,j)**3 / (nuair*alair)
        theta_na = pi/2.0 - beta(i,j)         ! angle from the VERTICAL (Kaplani convention)
        IF (unstab) THEN
          IF (beta(i,j) > pi/6.0) THEN        ! theta < 60 deg
            Grc_na = 1.327e10 * EXP(-3.708*theta_na)
            Rac_na = Grc_na * Pr
            Nu_nat = 0.13*MAX(Ra_na**(1.0/3.0) - Rac_na**(1.0/3.0), 0.0)  &
                   + 0.56*( Rac_na*MAX(SIN(beta(i,j)),0.0) )**0.25
          ELSE                           ! theta >= 60 deg (near horizontal, turbulent plume)
            Nu_nat = 0.13 * Ra_na**(1.0/3.0)
          END IF
        ELSE
          IF (beta(i,j) >= pi/6.0) THEN       ! beta >= 30 deg : Kaplani Eq.(8), all Ra
            Nu_nat = ( 0.825 + 0.387*( Ra_na*MAX(SIN(beta(i,j)),0.0) )**(1.0/6.0) / cc_pr )**2
          ELSE IF (beta(i,j) > 2.0*pi/180.0) THEN  ! 2 < beta < 30 deg : Kaplani Eq.(9), inclined
            Nu_nat = 0.56 * ( Ra_na*MAX(SIN(beta(i,j)),0.0) )**0.25
          ELSE                           ! beta <= 2 deg : Kaplani Eq.(9), horizontal plate
            Nu_nat = 0.58 * Ra_na**0.2
          END IF
        END IF
        Nu_nat = MAX(Nu_nat, 0.0)
        h_free   = pv_fac_h_free * Nu_nat * kair / L(i,j)
        h_cv     = ( h_forced**3 + h_free**3 )**(1.0/3.0)
        cdh_PVf(i,j) = h_cv / ( rhosfc(i,j) * cp * U_PVa(i,j) )
      END DO
    END DO
    cdq_PVf(:,:) = cdh_PVf(:,:)

    ! --- canopy air -> reference, with stability ---
    ! zp_dc(3) is set so that dragcoef_pv's internal z1 = z_ref - max(z_can, zd_park).
    zp_dc(:,:,1) = zp(:,:,1)
    zp_dc(:,:,2) = 0.0

    DO j = 1, ny-1
      DO i = 1, nx-1
        z_ref = zp(i,j,2) + 0.5*(zp(i,j,3)-zp(i,j,2))
        zp_dc(i,j,3) = 2.0 * (z_ref - max(z_can(i,j), zd_park(i,j)))
        pt_can(i,j) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
        wdspsfc(i,j) = uref(i,j)
      END DO
    END DO
    
    CALL dragcoef_pv(nx, ny, nz, zp_dc, &
                 pt_can, qvref, soiltyp_dummy, z0_park_arr, pv_ztz0_air, &
                 wdspsfc, pt1, &
                 cdm_can, cdh_can, cdq_can)
    
    ! --- air above the panels -> reference ---
    DO j = 1, ny-1
      DO i = 1, nx-1
        z_ref = zp(i,j,2) + 0.5*(zp(i,j,3)-zp(i,j,2))
        zp_dc(i,j,3) = 2.0 * (z_ref - max(z_PVa(i,j), zd_park(i,j)))
      END DO
    END DO
    
    CALL dragcoef_pv(nx, ny, nz, zp_dc, &
                 pt_PVa, qvref, soiltyp_dummy, z0_park_arr, pv_ztz0_air, &
                 wdspsfc, pt1, &
                 cdm_PVa, cdh_PVa, cdq_PVa)
    
  RETURN
END SUBROUTINE drag_pvland


!##################################################################
!##################################################################
!######                                                      ######
!######      SUBROUTINE SOLVE_AIR_NODE_TEMPERATURES          ######
!######                                                      ######
!######   Solve for T_can and T_PVa using Eq. 2.3            ######
!######   C_can = C_PVa = 0 => q_conv = 0                    ######
!######   Direct algebraic solution with coefficients at t   ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE solve_air_node_temperatures(nx, ny, PV_PARK, &
                                       T_nat, T_PVb, T_PVf, T_ref, &
                                       U_can, U_PVa, U_ref, &
                                       cdh_nat, cdh_PVb, cdh_can, &
                                       cdh_PVf, cdh_PVa, &
                                       L, a, w, &
                                       T_can, T_PVa)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Solve for air node temperatures T_can and T_PVa
!  using energy balance with C=0 (Eq. 2.3).
!
!  From the condition q_conv = 0 (Table 2):
!
!  For T_can:
!    H_nat + (L/(w+a))*H_PVb - (a/(w+a))*H_can–ref = 0
!
!  For T_PVa:
!    (L/w)*H_PVf - H_PVa–ref = 0
!
!  Using H = ρ*cp*Cdh*U*(T_surface - T_air), we get:
!
!  T_can = [Cdh_nat*U_can*T_nat + (L/(w+a))*Cdh_PVb*U_can*T_PVb 
!           + (a/(w+a))*Cdh_can*U_ref*T_ref]
!        / [Cdh_nat*U_can + (L/(w+a))*Cdh_PVb*U_can 
!           + (a/(w+a))*Cdh_can*U_ref]
!
!  T_PVa = [(L/w)*Cdh_PVf*U_PVa*T_PVf + Cdh_PVa*U_ref*T_ref]
!        / [(L/w)*Cdh_PVf*U_PVa + Cdh_PVa*U_ref]
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  :: nx, ny
  REAL, INTENT(IN)     :: PV_PARK(nx,ny)
  REAL, INTENT(IN)     :: T_nat(nx,ny), T_PVb(nx,ny)
  REAL, INTENT(IN)     :: T_PVf(nx,ny), T_ref(nx,ny)
  REAL, INTENT(IN)     :: U_can(nx,ny), U_PVa(nx,ny), U_ref(nx,ny)
  REAL, INTENT(IN)     :: cdh_nat(nx,ny), cdh_PVb(nx,ny)
  REAL, INTENT(IN)     :: cdh_can(nx,ny), cdh_PVf(nx,ny)
  REAL, INTENT(IN)     :: cdh_PVa(nx,ny)
  REAL, INTENT(IN)     :: L(nx,ny), a(nx,ny), w(nx,ny)
  REAL, INTENT(OUT)    :: T_can(nx,ny), T_PVa(nx,ny)
  
  REAL :: numerator_can, denominator_can
  REAL :: numerator_PVa, denominator_PVa
  REAL :: ratio_L_wa, ratio_a_wa, ratio_L_w
  INTEGER :: i, j
  
  DO j = 1, ny-1
    DO i = 1, nx-1
    
    IF (PV_PARK(i,j).eq.1.) THEN
    
  ! Pre-compute geometric ratios
  ratio_L_wa = L(i,j) / (w (i,j)+ a(i,j))
  ratio_a_wa = a(i,j) / (w(i,j) + a(i,j))
  ratio_L_w = L(i,j) / w(i,j)
  
      !-----------------------------------------------------------------
      ! Solve for T_can
      !-----------------------------------------------------------------
      numerator_can = cdh_nat(i,j) * U_can(i,j) * T_nat(i,j) + &
                     ratio_L_wa * cdh_PVb(i,j) * U_can(i,j) * T_PVb(i,j) + &
                     ratio_a_wa * cdh_can(i,j) * U_ref(i,j) * T_ref(i,j)
      
      denominator_can = cdh_nat(i,j) * U_can(i,j) + &
                       ratio_L_wa * cdh_PVb(i,j) * U_can(i,j) + &
                       ratio_a_wa * cdh_can(i,j) * U_ref(i,j)
      
      T_can(i,j) = numerator_can / denominator_can
      
      !-----------------------------------------------------------------
      ! Solve for T_PVa
      !-----------------------------------------------------------------
      numerator_PVa = ratio_L_w * cdh_PVf(i,j) * U_PVa(i,j) * T_PVf(i,j) + &
                     cdh_PVa(i,j) * U_ref(i,j) * T_ref(i,j)
      
      denominator_PVa = ratio_L_w * cdh_PVf(i,j) * U_PVa(i,j) + &
                       cdh_PVa(i,j) * U_ref(i,j)
      
      T_PVa(i,j) = numerator_PVa / denominator_PVa
      
      ENDIF
      
    END DO
  END DO

END SUBROUTINE solve_air_node_temperatures

!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE COMPUTE_ALL_FORCINGS            ######
!######                                                      ######
!######   Calculate forcing terms (dT/dt) for ALL nodes      ######
!######   Called TWICE per timestep for RK2 integration       ######
!######   Uses soilebm_frc and sfcflx                        ######
!######                                                      ######
!######   COMBINED version: outputs both PV diagnostics      ######
!######   (diag_qLW_*, diag_H_*, etc.) AND nat surface       ######
!######   irradiances/radiosities (G_SW_nat, J_SW_nat, etc.) ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE compute_all_forcings(nx, ny, nz, nzsoil, pv_park, rnflx_arps, &
                               L, a, beta, h_PV, h_clear, w, &
                               z0_nat, &
                               albedo_nat, albedo_PVf, albedo_PVb, &
                               emiss_nat, emiss_PVf, emiss_PVb, &
                               soiltyp, vegtyp, lai, veg, &
                               T_nat, T_soil, T_PVf, T_cell, T_PVb, &
                               T_can, T_PVa, &
                               wetsfc, wetdp, wetcanp, snowdpth, qvsfc, &
                               uref, tref, qvref, swin, lwin, psfc, rhosfc, precip, &
                               U_can, U_PVa, U_ref, &
                               lambda, phi, S_diffuse_in, &
                               delta_PVf, delta_cell, delta_PVb, &
                               rho_cp_PVf, rho_cp_cell, rho_cp_PVb, &
                               delta_glass, delta_EVA_top, delta_EVA_bot, delta_backsheet, &
                               k_glass, k_EVA, k_backsheet, &
                               eta_STC, gamma_PV, &
                               zp, u_wind, v_wind, w_wind, ptprt, pprt, qv, ptbar, pbar, qvbar, &
                               cdh, cdh_nat, cdh_PVb, cdh_PVf, &
                               cdm, cdm_nat, cdm_PVb, cdm_PVf, &
                               cdq, cdq_nat, cdq_PVb, cdq_PVf, &
                               cdq_can, q_can_prev, &
                               frc_nat, frc_soil, frc_PVf, frc_cell, frc_PVb, &
                               frc_wsfc, frc_wdp, frc_wcnp, wrmax, &
                               H_nat_out, LE_nat_out, G_nat_out, Rn_nat_out, Rn_ref_out, &
                               SW_out_sky_out, LW_out_sky_out, &
                               diag_qLW_PVf, diag_qLW_PVb, &
                               diag_H_PVf, diag_H_PVb, &
                               diag_qcond_PVf, diag_qcond_PVb, &
                               diag_qSW_PVf, diag_qSW_PVb, &
                               diag_qPV, &
                               diag_cdh_can_out, diag_cdh_PVa_out, &
                               diag_T_can_out, diag_T_PVa_out, &
                               diag_JLW_PVf_out, diag_JLW_PVb_out, &
                               G_SW_nat_out, J_SW_nat_out, &
                               G_LW_nat_out, J_LW_nat_out, &
                               q_can_out, gcan_out, vpd_f3_out, vpd_ref_out, full)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate forcing terms (dT/dt in K/s) for all thermal nodes
!  (natural surface, soil, PV front/cell/back).
!
!  COMBINED version outputting:
!    - PV energy balance diagnostics (from _pv version)
!    - Natural surface irradiances/radiosities (from _nat version)
!
!-----------------------------------------------------------------------
!
  USE pvland_rad_cache   ! view factors, LU factors, SW solution  [SYLVAIN]

  IMPLICIT NONE
  
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
  
  REAL, INTENT(IN) :: pv_park(nx,ny)
  REAL, INTENT(IN) :: rnflx_arps(nx,ny)
  
  ! Grid
  INTEGER, INTENT(IN) :: nx, ny, nz, nzsoil
  
  ! Geometry
  REAL, INTENT(IN) :: L(nx,ny), a(nx,ny), beta(nx,ny), h_PV(nx,ny), h_clear(nx,ny), w(nx,ny)
  REAL, INTENT(IN) :: z0_nat(nx,ny)
  
  ! Optical properties
  REAL, INTENT(IN) :: albedo_nat(nx,ny), albedo_PVf(nx,ny), albedo_PVb(nx,ny)
  REAL, INTENT(IN) :: emiss_nat(nx,ny), emiss_PVf(nx,ny), emiss_PVb(nx,ny)
  
  ! Soil/vegetation properties
  INTEGER, INTENT(IN) :: soiltyp(nx,ny), vegtyp(nx,ny)
  REAL, INTENT(IN) :: lai(nx,ny), veg(nx,ny)
  
  ! Current temperatures (input state)
  REAL, INTENT(IN) :: T_nat(nx,ny), T_soil(nx,ny)
  REAL, INTENT(IN) :: T_PVf(nx,ny), T_cell(nx,ny), T_PVb(nx,ny)
  REAL, INTENT(IN) :: T_can(nx,ny), T_PVa(nx,ny)
  
  ! Soil moisture
  REAL, INTENT(IN) :: wetsfc(nx,ny), wetdp(nx,ny)
  REAL, INTENT(IN) :: wetcanp(nx,ny), qvsfc(nx,ny)
  REAL, INTENT(IN) :: snowdpth(nx,ny)
  
  ! Atmospheric forcing
  REAL, INTENT(IN) :: uref(nx,ny), tref(nx,ny), qvref(nx,ny)
  REAL, INTENT(IN) :: swin(nx,ny), lwin(nx,ny)
  REAL, INTENT(IN) :: psfc(nx,ny), rhosfc(nx,ny), precip(nx,ny)
  REAL, INTENT(IN) :: U_can(nx,ny), U_PVa(nx,ny), U_ref(nx,ny)
  
  ! Solar geometry
  REAL, INTENT(IN) :: lambda(nx,ny), phi(nx,ny), S_diffuse_in(nx,ny)
  
  ! PV thermal parameters
  REAL, INTENT(IN) :: delta_PVf(nx,ny), delta_cell(nx,ny), delta_PVb(nx,ny)
  REAL, INTENT(IN) :: rho_cp_PVf(nx,ny), rho_cp_cell(nx,ny), rho_cp_PVb(nx,ny)
  REAL, INTENT(IN) :: delta_glass(nx,ny), delta_EVA_top(nx,ny), delta_EVA_bot(nx,ny), delta_backsheet(nx,ny)
  REAL, INTENT(IN) :: k_glass(nx,ny), k_EVA(nx,ny), k_backsheet(nx,ny)
  REAL, INTENT(IN) :: eta_STC(nx,ny), gamma_PV(nx,ny)
  
  ! Atmospheric state (for sfcflx)
  REAL, INTENT(IN) :: zp(nx,ny,nz)
  REAL, INTENT(IN) :: u_wind(nx,ny,nz), v_wind(nx,ny,nz), w_wind(nx,ny,nz)
  REAL, INTENT(IN) :: ptprt(nx,ny,nz), pprt(nx,ny,nz), qv(nx,ny,nz)
  REAL, INTENT(IN) :: ptbar(nx,ny,nz), pbar(nx,ny,nz), qvbar(nx,ny,nz)
  
  ! Drag coefficients (computed in pv_land, once per small timestep)
  REAL, INTENT(IN) :: cdh(nx,ny)
  REAL, INTENT(IN) :: cdh_nat(nx,ny)
  REAL, INTENT(IN) :: cdh_PVb(nx,ny)
  REAL, INTENT(IN) :: cdh_PVf(nx,ny)
  REAL, INTENT(IN) :: cdm(nx,ny)
  REAL, INTENT(IN) :: cdm_nat(nx,ny)
  REAL, INTENT(IN) :: cdm_PVb(nx,ny)
  REAL, INTENT(IN) :: cdm_PVf(nx,ny)
  REAL, INTENT(IN) :: cdq(nx,ny)
  REAL, INTENT(IN) :: cdq_nat(nx,ny)
  REAL, INTENT(IN) :: cdq_PVb(nx,ny)
  REAL, INTENT(IN) :: cdq_PVf(nx,ny)
  REAL, INTENT(IN) :: cdq_can(nx,ny)      ! Moisture drag coeff. canopy air -> reference
  REAL, INTENT(IN) :: q_can_prev(nx,ny)   ! Canopy-air humidity of the previous small timestep (kg/kg), for F3
  
  ! OUTPUT: temperature forcing terms (K/s)
  REAL, INTENT(OUT) :: frc_nat(nx,ny), frc_soil(nx,ny)
  REAL, INTENT(OUT) :: frc_PVf(nx,ny), frc_cell(nx,ny), frc_PVb(nx,ny)
  
  ! OUTPUT: soil moisture forcing terms
  REAL, INTENT(OUT) :: frc_wsfc(nx,ny), frc_wdp(nx,ny), frc_wcnp(nx,ny)
  REAL, INTENT(OUT) :: wrmax(nx,ny)
  
  ! Diagnostic outputs
  REAL, INTENT(OUT) :: H_nat_out(nx,ny), LE_nat_out(nx,ny), G_nat_out(nx,ny)
  REAL, INTENT(OUT) :: Rn_ref_out(nx,ny), Rn_nat_out(nx,ny) !modif sylvain
  REAL, INTENT(OUT) :: SW_out_sky_out(nx,ny), LW_out_sky_out(nx,ny)
  
  ! PV energy balance diagnostics (W/m2) -- from _pv version
  REAL, INTENT(OUT) :: diag_qLW_PVf(nx,ny), diag_qLW_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_H_PVf(nx,ny), diag_H_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_qcond_PVf(nx,ny), diag_qcond_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_qSW_PVf(nx,ny), diag_qSW_PVb(nx,ny)
  REAL, INTENT(OUT) :: diag_qPV(nx,ny)
  REAL, INTENT(OUT) :: diag_cdh_can_out(nx,ny), diag_cdh_PVa_out(nx,ny)
  REAL, INTENT(OUT) :: diag_T_can_out(nx,ny), diag_T_PVa_out(nx,ny)
  REAL, INTENT(OUT) :: diag_JLW_PVf_out(nx,ny), diag_JLW_PVb_out(nx,ny)
  
  ! Natural surface irradiances/radiosities -- from _nat version
  REAL, INTENT(OUT) :: G_SW_nat_out(nx,ny), J_SW_nat_out(nx,ny)
  REAL, INTENT(OUT) :: G_LW_nat_out(nx,ny), J_LW_nat_out(nx,ny)
  ! Canopy-air humidity diagnostics (PV park; outside: q_can = q_ref, gcan = 0)
  REAL, INTENT(OUT) :: q_can_out(nx,ny)    ! Canopy-air specific humidity (kg/kg)
  REAL, INTENT(OUT) :: gcan_out(nx,ny)     ! Conductance canopy air -> reference [kg m-2 s-1]
  REAL, INTENT(OUT) :: vpd_f3_out(nx,ny)   ! VPD used in F3 [kPa]
  REAL, INTENT(OUT) :: vpd_ref_out(nx,ny)  ! qvsat(T_air) - q_ref [kPa]
  ! .TRUE.: all forcings and diagnostics. .FALSE.: only the soil part
  ! (H_nat, LE_nat, G_nat, q_can, soil-moisture forcings), used for the
  ! end-of-step fluxes when no output is written: the LW radiosities, the
  ! fluxes and the energy balance of the PV faces, the radiation
  ! diagnostics and frc_nat (it needs rnflx) are not computed; the SW part
  ! comes from the cache (unchanged state of the sun).  [SYLVAIN]
  LOGICAL, INTENT(IN) :: full
  
  ! Local - Radiation (8-surface radiosity system)
  REAL :: F_matrix(8,8), F_sky(8)
  REAL :: S_direct(8), L_down
  REAL :: alpha_surf(8), epsilon_surf(8), T_surf(8)
  REAL :: J_SW(8), G_SW(8), q_SW(8)
  REAL :: J_LW(8), G_LW(8), q_LW(8)
  INTEGER :: info_rad
  
  ! Radiative fluxes per unit area (W/m2)
  ! (SAVE: allocated at the first call only)  [SYLVAIN]
  REAL , ALLOCATABLE, DIMENSION (:,:), SAVE :: q_dot_SW_PVf, q_dot_LW_PVf
  REAL, ALLOCATABLE, DIMENSION (:,:), SAVE :: q_dot_SW_cell
  REAL, ALLOCATABLE, DIMENSION (:,:), SAVE :: q_dot_SW_PVb, q_dot_LW_PVb
  REAL, ALLOCATABLE, DIMENSION (:,:), SAVE :: q_dot_SW_nat, q_dot_LW_nat
  REAL, ALLOCATABLE, DIMENSION (:,:), SAVE :: G_SW_PVf, G_SW_nat
  
  ! Convective fluxes (W/m2) - from sfcflx
  REAL :: ptsflx_nat(nx,ny), qvsflx_nat(nx,ny)
  REAL :: ptsflx_PVb(nx,ny), qvsflx_PVb(nx,ny)
  REAL :: ptsflx_PVf(nx,ny), qvsflx_PVf(nx,ny)
  REAL :: usflx_dummy(nx,ny), vsflx_dummy(nx,ny)
  REAL :: H_PVf, H_PVb, H_nat
  
  ! Local 3D arrays for sfcflx
  REAL :: ptbar_loc(nx,ny,nz), ptprt_loc(nx,ny,nz)
  REAL :: ptsfc_loc(nx,ny)
  REAL :: uref_nat(nx,ny),tref_nat(nx,ny),cdhref_nat(nx,ny),cdqref_nat(nx,ny)
  
  ! Soil variables (for soilebm_frc)
  REAL :: frc_tsfc(nx,ny), frc_tdeep(nx,ny)
  REAL :: qvsata(nx,ny), qvsata_pv(nx,ny),qvsat_local(nx,ny), f34(nx,ny)
  REAL :: ct(nx,ny), shflx(nx,ny), lhflx(nx,ny), gflx(nx,ny), rnflx(nx,ny)
  REAL :: evaprg(nx,ny), evaprtr(nx,ny), evaprr(nx,ny)
  REAL :: relief
  REAL :: radsw(nx,ny)
  REAL :: f3
  REAL :: VPD_kPa
  
  ! Conductive fluxes (W/m2)
  REAL :: q_cond_PVf, q_cond_PVb
  REAL :: R_PVf, R_PVb
  
  ! PV electrical power (W/m2)
  REAL :: q_dot_PV, eta_PV, T_cell_C
  
  ! Heat capacities (J/m2/K)
  REAL :: C_PVf, C_cell, C_PVb
  
  INTEGER :: i, j, is
  REAL, PARAMETER :: pi = 3.14159265359
  
  ! Debug call counter
  INTEGER, SAVE :: call_count = 0
  
  !=====================================================================
  ! 1. RADIATION (temperature-dependent)
  !=====================================================================

IF (.NOT. ALLOCATED(q_dot_SW_PVf)) &
ALLOCATE (q_dot_SW_PVf(nx,ny), q_dot_LW_PVf(nx,ny),q_dot_SW_cell(nx,ny),&
 q_dot_SW_PVb(nx,ny), q_dot_LW_PVb(nx,ny),q_dot_SW_nat(nx,ny), &
 q_dot_LW_nat(nx,ny),G_SW_PVf(nx,ny), G_SW_nat(nx,ny))

    DO j = 1, ny-1
      DO i = 1, nx-1
  
      IF (pv_park(i,j).eq.1.) THEN

  ! View factors (geometry only): computed once in INIT_PVLAND  [SYLVAIN]
  !CALL compute_view_factor_matrix(L(i,j), a(i,j), w(i,j), h_PV(i,j), beta(i,j), F_matrix)
  !CALL compute_sky_view_factors(F_matrix, F_sky)
  F_matrix = F_pv(:,:,i,j)
  F_sky    = Fsky_pv(:,i,j)
  
  ! Surface optical properties (8-surface periodic array layout)
  ! Surfaces: 1=nat_n, 2=PVf_n, 3=PVb_n, 4=PVb_{n-1}, 5=PVf_{n-1},
  !           6=PVb_{n+1}, 7=nat_{n-1}, 8=nat_{n+1}
  alpha_surf(1) = albedo_nat(i,j)       ! nat,n
  alpha_surf(2) = albedo_PVf(i,j)       ! PVf,n
  alpha_surf(3) = albedo_PVb(i,j)       ! PVb,n
  alpha_surf(4) = albedo_PVb(i,j)       ! PVb,n-1  
  alpha_surf(5) = albedo_PVf(i,j)       ! PVf,n-1
  alpha_surf(6) = albedo_PVb(i,j)       ! PVb,n+1
  alpha_surf(7) = albedo_nat(i,j)       ! nat,n-1
  alpha_surf(8) = albedo_nat(i,j)       ! nat,n+1  

  epsilon_surf(1) = emiss_nat(i,j)      ! nat,n
  epsilon_surf(2) = emiss_PVf(i,j)      ! PVf,n
  epsilon_surf(3) = emiss_PVb(i,j)      ! PVb,n
  epsilon_surf(4) = emiss_PVb(i,j)      ! PVb,n-1  
  epsilon_surf(5) = emiss_PVf(i,j)      ! PVf,n-1
  epsilon_surf(6) = emiss_PVb(i,j)      ! PVb,n+1
  epsilon_surf(7) = emiss_nat(i,j)      ! nat,n-1
  epsilon_surf(8) = emiss_nat(i,j)      ! nat,n+1  

  T_surf(1) = T_nat(i,j)           ! nat,n
  T_surf(2) = T_PVf(i,j)           ! PVf,n
  T_surf(3) = T_PVb(i,j)           ! PVb,n
  T_surf(4) = T_PVb(i,j)           ! PVb,n-1 
  T_surf(5) = T_PVf(i,j)           ! PVf,n-1
  T_surf(6) = T_PVb(i,j)           ! PVb,n+1
  T_surf(7) = T_nat(i,j)           ! nat,n-1
  T_surf(8) = T_nat(i,j)           ! nat,n+1 
  
  ! SW part: independent of the temperatures, so recomputed only when its
  ! inputs change (compute_all_forcings is called 3 times per time step
  ! with the same solar forcing).  [SYLVAIN]
  IF ( sw_valid(i,j) .AND. sw_key(1,i,j) == swin(i,j)         &
       .AND. sw_key(2,i,j) == S_diffuse_in(i,j)               &
       .AND. sw_key(3,i,j) == lambda(i,j)                     &
       .AND. sw_key(4,i,j) == phi(i,j) ) THEN
  J_SW = J_SW_pv(:,i,j)
  G_SW = G_SW_pv(:,i,j)
  q_SW = q_SW_pv(:,i,j)
  ELSE

  ! Direct SW on each surface (diffuse is handled separately in radiosity)
  CALL compute_S_ground(swin(i,j) - S_diffuse_in(i,j), L(i,j), a(i,j), beta(i,j), h_PV(i,j), lambda(i,j), phi(i,j), S_direct(1))
  CALL compute_S_PVf(swin(i,j) - S_diffuse_in(i,j), L(i,j), a(i,j), beta(i,j), h_PV(i,j), lambda(i,j), phi(i,j), S_direct(2))
  CALL compute_S_PVb(swin(i,j) - S_diffuse_in(i,j), L(i,j), a(i,j), beta(i,j), h_PV(i,j), lambda(i,j), phi(i,j), S_direct(3))
  S_direct(4) = S_direct(3)
  S_direct(5) = S_direct(2)
  S_direct(6) = S_direct(3)
  S_direct(7) = S_direct(1)
  S_direct(8) = S_direct(1)
  
  ! SW radiosity (matrix factorised once in INIT_PVLAND)  [SYLVAIN]
  !CALL solve_SW_radiosities(8, alpha_surf, F_matrix, F_sky, S_direct, &
  !                          S_diffuse_in(i,j), J_SW, info_rad)
  CALL solve_SW_radiosities_lu(8, alpha_surf, LU_SW_pv(:,:,i,j), ipiv_SW_pv(:,i,j), &
                            F_sky, S_direct, S_diffuse_in(i,j), J_SW, info_rad)
  CALL compute_SW_irradiances(8, F_matrix, F_sky, S_direct, S_diffuse_in(i,j), &
                             J_SW, G_SW)
  CALL compute_SW_net_fluxes(8, alpha_surf, G_SW, q_SW)

  J_SW_pv(:,i,j) = J_SW
  G_SW_pv(:,i,j) = G_SW
  q_SW_pv(:,i,j) = q_SW
  sw_key(:,i,j)  = (/ swin(i,j), S_diffuse_in(i,j), lambda(i,j), phi(i,j) /)
  sw_valid(i,j)  = .TRUE.
  ENDIF
  
  ! LW radiosity (temperature-dependent; matrix factorised once in INIT_PVLAND)
  !CALL solve_LW_radiosities(8, epsilon_surf, T_surf, F_matrix, F_sky, &
  !                         Lwin(i,j), J_LW, info_rad)
  IF (full) THEN
  CALL solve_LW_radiosities_lu(8, epsilon_surf, T_surf, LU_LW_pv(:,:,i,j), ipiv_LW_pv(:,i,j), &
                           F_sky, Lwin(i,j), J_LW, info_rad)
  CALL compute_LW_irradiances(8, F_matrix, F_sky, J_LW, Lwin(i,j), G_LW)
  CALL compute_LW_net_fluxes(8, G_LW, J_LW, q_LW)
  
  ! Store LW radiosities for PV diagnostic output (from _pv)
  diag_JLW_PVf_out(i,j) = J_LW(2)
  diag_JLW_PVb_out(i,j) = J_LW(3)
  
  ! Extract per-surface net fluxes
  q_dot_SW_nat(i,j) = q_SW(1)
  q_dot_LW_nat(i,j) = q_LW(1)
  q_dot_SW_PVf(i,j) = q_SW(2)
  q_dot_LW_PVf(i,j) = q_LW(2)
  q_dot_SW_PVb(i,j) = q_SW(3)
  q_dot_LW_PVb(i,j) = q_LW(3)
  G_SW_PVf(i,j) = G_SW(2)
  G_SW_nat(i,j) = G_SW(1)
  
  ! Store nat surface irradiances and radiosities for output (from _nat)
  G_SW_nat_out(i,j) = G_SW(1)   ! SW irradiance received by nat (SW in nat)
  J_SW_nat_out(i,j) = J_SW(1)   ! SW radiosity leaving nat (SW out nat)
  G_LW_nat_out(i,j) = G_LW(1)   ! LW irradiance received by nat (LW in nat)
  J_LW_nat_out(i,j) = J_LW(1)   ! LW radiosity leaving nat (LW out nat)
  Rn_nat_out(i,j) = q_dot_SW_nat(i,j) + q_dot_LW_nat(i,j) !modif sylvain
  
  ! Cell absorbs SW transmitted through front glass
  q_dot_SW_cell(i,j) = q_dot_SW_PVf(i,j)
  
  ! Area-weighted net radiation for the park
  Rn_ref_out(i,j) = (L(i,j) * (q_dot_SW_PVf(i,j) + q_dot_LW_PVf(i,j)) + &
                L(i,j) * (q_dot_SW_PVb(i,j) + q_dot_LW_PVb(i,j)) + &
                (a(i,j) + w(i,j)) * (q_dot_SW_nat(i,j) + q_dot_LW_nat(i,j))) / (a(i,j) + w(i,j))
  
  ! Outgoing fluxes to sky (as seen by a horizontal radiometer above the park)
  SW_out_sky_out(i,j) = (L(i,j) * F_sky(2) * J_SW(2) + &
                    L(i,j) * F_sky(3) * J_SW(3) + &
                    (a(i,j) + w(i,j)) * F_sky(1) * J_SW(1)) / (a(i,j) + w(i,j))
  
  LW_out_sky_out(i,j) = (L(i,j) * F_sky(2) * J_LW(2) + &
                    L(i,j) * F_sky(3) * J_LW(3) + &
                    (a(i,j) + w(i,j)) * F_sky(1) * J_LW(1)) / (a(i,j) + w(i,j))
  
  ! Setup radiation inputs for soilebm_frc
  ENDIF  ! full
  radsw(i,j) = G_SW(1)  ! Total SW irradiance on natural surface (from radiosity)
  rnflx(i,j) = q_dot_SW_nat(i,j) + q_dot_LW_nat(i,j)  ! Net radiation for natural surface

  ELSE
  
  Rn_nat_out(i,j) = rnflx_arps(i,j)
  Rn_ref_out(i,j) = Rn_nat_out(i,j)

  ! Setup radiation inputs for soilebm_frc
  radsw(i,j) = swin(i,j)  ! Total SW irradiance on natural surface (from radiosity)
  rnflx(i,j) = rnflx_arps(i,j)  ! Net radiation for natural surface
  
  ENDIF
  
  ENDDO
  ENDDO

  !=====================================================================
  ! 2. CONVECTION - Using sfcflx with proper 3D arrays
  !=====================================================================
  
  IF (full) THEN   ! fluxes of the PV faces
  ptprt_loc = 0.0
  
  ! (1) Natural surface -> Canopy air
  !DO j = 1, ny-1
  !  DO i = 1, nx-1
  !    ptbar_loc(i,j,:) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
  !    ptsfc_loc(i,j)   = T_nat(i,j) * (p0/pbar(i,j,2))**rddcp
  !  END DO
  !END DO
  
  !CALL sfcflx(nx, ny, nz, &
  !           u_wind, v_wind, w_wind, ptprt_loc, pprt, qv, ptbar_loc, pbar, qvbar, &
  !           U_can, rhosfc, ptsfc_loc, qvref, &
  !           cdm_nat, cdh_nat, cdq_nat, &
  !           usflx_dummy, vsflx_dummy, ptsflx_nat, qvsflx_nat)
  
  ! (2) PV back -> Canopy air
  DO j = 1, ny-1
    DO i = 1, nx-1
      ptbar_loc(i,j,:) = T_can(i,j) * (p0/pbar(i,j,2))**rddcp
      ptsfc_loc(i,j) = T_PVb(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt_loc, pprt, qv, ptbar_loc, pbar, qvbar, &
             U_can, rhosfc, ptsfc_loc, qvref, &
             cdm_PVb, cdh_PVb, cdq_PVb, &
             usflx_dummy, vsflx_dummy, ptsflx_PVb, qvsflx_PVb)
  
  ! (3) PV front -> Air above PV
  DO j = 1, ny-1
    DO i = 1, nx-1
      ptbar_loc(i,j,:) = T_PVa(i,j) * (p0/pbar(i,j,2))**rddcp
      ptsfc_loc(i,j)   = T_PVf(i,j) * (p0/pbar(i,j,2))**rddcp
    END DO
  END DO
  
  CALL sfcflx(nx, ny, nz, &
             u_wind, v_wind, w_wind, ptprt_loc, pprt, qv, ptbar_loc, pbar, qvbar, &
             U_PVa, rhosfc, ptsfc_loc, qvref, &
             cdm_PVf, cdh_PVf, cdq_PVf, &
             usflx_dummy, vsflx_dummy, ptsflx_PVf, qvsflx_PVf)
  ENDIF  ! full
  
  !=====================================================================
  ! 3. SOIL FORCINGS (via soilebm_frc)
  !=====================================================================
  
  ! Humidity of the canopy air (air under the panels, version 2): water
  ! vapour from the natural surface goes through the canopy air, a node
  ! without storage like T_can, then to the reference level through
  !    gcan = (a/(w+a)) rho_ref cdq_can U_ref   [kg m-2 s-1]
  ! q_can is computed in soilebm_frc_pv/evapflx_pv. F3 (Alfieri et al.
  ! 2008) uses the vapour pressure deficit of the canopy air,
  ! qvsat(T_can) - q_can, with q_can of the previous small timestep.
  ! Outside the park: no canopy-air node (gcan = 0, q_can = q_ref) and F3,
  ! F4 from the air at the reference level.  [SYLVAIN]
  CALL getqvs(nx, ny, 1, 1, nx-1, 1, ny-1, 1, 1, psfc, T_can, qvsata_pv)
  CALL getqvs(nx, ny, 1, 1, nx-1, 1, ny-1, 1, 1, psfc, Tref, qvsata)
  
  DO j = 1, ny-1
    DO i = 1, nx-1
    
      IF (pv_park(i,j).eq.1.) THEN

      gcan_out(i,j) = a(i,j) / (w(i,j) + a(i,j)) * rhosfc(i,j) * cdq_can(i,j) * U_ref(i,j)

      ! f3: Alfieri et al. (2008), Eq. 7, g = 0.1914 1/kPa (Noah LSM)
      VPD_kPa = (qvsata_pv(i,j) - q_can_prev(i,j)) * psfc(i,j) / 622.0
      vpd_f3_out(i,j)  = VPD_kPa
      vpd_ref_out(i,j) = (qvsata_pv(i,j) - qvref(i,j)) * psfc(i,j) / 622.0
      f3 = MAX(0.25, MIN(1.0, 1.0 - 0.1914*VPD_kPa))
      
      ! f4: stomatal response to temperature
      IF (T_can(i,j) >= 302.15) THEN
        f34(i,j) = 1.0 / (1.0 + EXP(1.18 * (T_can(i,j) - 314.00)))
      ELSE
        f34(i,j) = 1.0 / (1.0 + EXP(-0.41 * (T_can(i,j) - 282.05)))
      END IF
      
      uref_nat(i,j) = U_can(i,j)
      tref_nat(i,j) = T_can(i,j)
      cdhref_nat(i,j) = cdh_nat(i,j)
      cdqref_nat(i,j) = cdq_nat(i,j)
      
      ELSE

      gcan_out(i,j) = 0.0

      ! f3: Alfieri et al. (2008), air at the reference level
      VPD_kPa = (qvsata(i,j) - qvref(i,j)) * psfc(i,j) / 622.0
      vpd_f3_out(i,j)  = VPD_kPa
      vpd_ref_out(i,j) = VPD_kPa
      f3 = MAX(0.25, MIN(1.0, 1.0 - 0.1914*VPD_kPa))
      
      ! f4: stomatal response to temperature
      IF (Tref(i,j) >= 302.15) THEN
        f34(i,j) = 1.0 / (1.0 + EXP(1.18 * (Tref(i,j) - 314.00)))
      ELSE
        f34(i,j) = 1.0 / (1.0 + EXP(-0.41 * (Tref(i,j) - 282.05)))
      END IF

      uref_nat(i,j) = uref(i,j)
      tref_nat(i,j) = tref(i,j)
      cdhref_nat(i,j) = cdh(i,j)
      cdqref_nat(i,j) = cdq(i,j)

      ENDIF

      f34(i,j) = f34(i,j) * f3
      
    END DO
  END DO
  
  relief = 0.0
  
  CALL soilebm_frc_pv(nx, ny, nz, soiltyp, vegtyp, lai, veg, &
                    T_nat, T_soil, wetsfc, wetdp, wetcanp, snowdpth, &
                    qvsfc, uref_nat, psfc, rhosfc, precip, &
                    tref_nat, qvref, cdhref_nat, cdqref_nat, radsw, rnflx, &
                    shflx, lhflx, gflx, ct, &
                    evaprg, evaprtr, evaprr, qvsat_local, qvsata, f34, &
                    frc_tsfc, frc_tdeep, frc_wsfc, frc_wdp, frc_wcnp, &
                    wrmax, relief, gcan_out, q_can_out)
  
  ! Extract soil forcings
  frc_nat = frc_tsfc
  frc_soil = frc_tdeep
  
  ! Extract output fluxes
  !H_nat_out  = -cp * ptsflx_nat
  H_nat_out  = shflx
  LE_nat_out = lhflx
  G_nat_out = gflx
  
  ! --- H_nat COMPARISON: sfcflx vs soilebm_frc ---
  !IF (call_count <= 15) THEN
  !  WRITE(*,'(A)') '--- H_nat COMPARISON ---'
  !  WRITE(*,'(A,F12.4)') '  sfcflx   H_nat  = ', -cp * ptsflx_nat(1,1)
  !  WRITE(*,'(A,F12.4)') '  soilebm  shflx  = ', shflx(1,1)
  !  WRITE(*,'(A,F12.4)') '  diff             = ', -cp * ptsflx_nat(1,1) - shflx(1,1)
  !  WRITE(*,'(A)') '------------------------'
  !END IF
  
  ! --- LE COMPARISON: sfcflx vs soilebm_frc ---
  !call_count = call_count + 1
  !IF (call_count <= 5) THEN
  !  WRITE(*,'(A)') '--- LE COMPARISON (nat surface) ---'
  !  WRITE(*,'(A,F12.4)') '  sfcflx  LE_nat  = ', -lathv * qvsflx_nat(1,1)
  !  WRITE(*,'(A,F12.4)') '  soilebm lhflx   = ', lhflx(1,1)
  !  WRITE(*,'(A,F12.4)') '  sfcflx  qvsflx  = ', qvsflx_nat(1,1)
  !  WRITE(*,'(A,F12.4)') '  soilebm evaprg  = ', evaprg(1,1)
  !  WRITE(*,'(A,F12.4)') '  soilebm evaprtr = ', evaprtr(1,1)
  !  WRITE(*,'(A,F12.4)') '  soilebm evaprr  = ', evaprr(1,1)
  !  WRITE(*,'(A)') '-----------------------------------'
  !END IF

  !=====================================================================
  ! 4. PV NODE FORCINGS
  !=====================================================================
  
  IF (full) THEN   ! energy balance of the PV nodes
  DO j = 1, ny-1
    DO i = 1, nx-1
      
      IF (pv_park(i,j).eq.1.) THEN

      ! Convective fluxes from sfcflx (positive = sfc loses heat)
      !H_nat = -cp * ptsflx_nat(i,j)
      H_PVb = -cp * ptsflx_PVb(i,j)
      H_PVf = -cp * ptsflx_PVf(i,j)
      
      ! Conduction through PV layers
      R_PVf = delta_glass(i,j) / k_glass(i,j) + delta_EVA_top(i,j) / k_EVA(i,j)
      R_PVb = delta_EVA_bot(i,j) / k_EVA(i,j) + delta_backsheet(i,j) / k_backsheet(i,j)
      
      q_cond_PVf = (T_cell(i,j) - T_PVf(i,j)) / R_PVf
      q_cond_PVb = (T_cell(i,j) - T_PVb(i,j)) / R_PVb
      
      ! PV electrical power output
      T_cell_C = T_cell(i,j) - 273.15
      eta_PV = eta_STC(i,j) * (1.0 - gamma_PV(i,j) * (T_cell_C - 25.0))
      q_dot_PV = eta_PV * G_SW_PVf(i,j)
      
      ! Heat capacities
      C_PVf = rho_cp_PVf(i,j) * delta_PVf(i,j)
      C_cell = rho_cp_cell(i,j) * delta_cell(i,j)
      C_PVb = rho_cp_PVb(i,j) * delta_PVb(i,j)
      
      ! Forcing terms: dT/dt = Q_net / C  (K/s)
      
      ! Front glass: LW radiation + conduction from cell - convective loss
      frc_PVf(i,j) = (q_dot_LW_PVf(i,j) + q_cond_PVf - H_PVf) / C_PVf
      
      ! Cell: SW absorbed - conduction to front and back - electrical power
      frc_cell(i,j) = (q_dot_SW_cell(i,j) - q_cond_PVf - q_cond_PVb - q_dot_PV) / C_cell
      
      !IF (i==1 .AND. j==1 .AND. call_count <= 4) THEN
      !  WRITE(*,'(A)') '--- CELL BALANCE ---'
      !  WRITE(*,'(A,F12.3)') '  q_SW_cell    = ', q_dot_SW_cell
      !  WRITE(*,'(A,F12.3)') '  q_cond_PVf   = ', q_cond_PVf
      !  WRITE(*,'(A,F12.3)') '  q_cond_PVb   = ', q_cond_PVb
      !  WRITE(*,'(A,F12.3)') '  q_dot_PV     = ', q_dot_PV
      !  WRITE(*,'(A,F12.3)') '  C_cell       = ', C_cell
      !  WRITE(*,'(A,ES12.4)') '  frc_cell     = ', frc_cell(1,1)
      !  WRITE(*,'(A,F12.3)') '  T_cell       = ', T_cell(1,1)
      !  WRITE(*,'(A,F12.3)') '  T_PVf        = ', T_PVf(1,1)
      !  WRITE(*,'(A,F12.3)') '  T_PVb        = ', T_PVb(1,1)
      !  WRITE(*,'(A)') '--------------------'
      !END IF
      
      ! Back surface: LW radiation + conduction from cell - convective loss
      frc_PVb(i,j) = (q_dot_LW_PVb(i,j) + q_dot_SW_PVb(i,j) + q_cond_PVb - H_PVb) / C_PVb
      
      ! Store PV energy balance diagnostics (from _pv)
      diag_qLW_PVf(i,j)  = q_dot_LW_PVf(i,j)
      diag_qLW_PVb(i,j)  = q_dot_LW_PVb(i,j)
      diag_H_PVf(i,j)    = H_PVf
      diag_H_PVb(i,j)    = H_PVb
      diag_qcond_PVf(i,j) = q_cond_PVf
      diag_qcond_PVb(i,j) = q_cond_PVb
      diag_qSW_PVf(i,j)  = q_dot_SW_PVf(i,j)
      diag_qSW_PVb(i,j)  = q_dot_SW_PVb(i,j)
      diag_qPV(i,j)      = q_dot_PV
      diag_cdh_can_out(i,j)  = cdh_nat(i,j)
      diag_cdh_PVa_out(i,j)  = cdh_PVf(i,j)
      diag_T_can_out(i,j) = T_can(i,j)
      diag_T_PVa_out(i,j) = T_PVa(i,j)
      
      ENDIF
      
    END DO
  END DO

  ENDIF  ! full

!DEALLOCATE (q_dot_SW_PVf, q_dot_LW_PVf,q_dot_SW_cell,&    ! kept (SAVE)  [SYLVAIN]
! q_dot_SW_PVb, q_dot_LW_PVb,q_dot_SW_nat, &
! q_dot_LW_nat,G_SW_PVf, G_SW_nat)

END SUBROUTINE compute_all_forcings

  !----------------------------------------------------------
  ! Projection dans le repere du parc : lambda (signe) et phi (deg)
  !----------------------------------------------------------
 subroutine compute_lambda_phi(zenith_deg, azimuth_deg, az_park_deg, &
                                 lambda_deg, phi_deg)
    real, intent(in)  :: zenith_deg, azimuth_deg, az_park_deg
    real, intent(out) :: lambda_deg, phi_deg

    real :: theta, gamma, alpha
    real :: sN, sE, sZ
    real :: s_x, s_y, s_z

   real, parameter :: PI = 3.14159265358979323846

    theta = zenith_deg  * PI / 180.0
    gamma = azimuth_deg * PI / 180.0
    alpha = az_park_deg * PI / 180.0

    ! Etape 5 : vecteur solaire dans le repere geographique
    sN = sin(theta) * cos(gamma)
    sE = sin(theta) * sin(gamma)
    sZ = cos(theta)

    ! Etape 6 : rotation vers le repere du parc
    s_x = cos(alpha) * sN + sin(alpha) * sE
    s_y = -sin(alpha) * sN + cos(alpha) * sE
    s_z = sZ

    ! Etape 7 : lambda = angle zenithal du soleil (non signe), version 2
    lambda_deg = acos(max(-1.0, min(1.0, s_z))) * 180.0 / PI
    ! Etape 8 : phi = azimut du soleil dans le repere du parc, compte depuis
    ! la direction des rangees (y) vers la normale avant des panneaux (x) :
    !   sin(lambda)*sin(phi) = s_x,  sin(lambda)*cos(phi) = s_y
    ! (convention de PV_LAND_final : tan(lambda)*sin(phi) = tangente de
    ! l'angle de profil ; le produit n.s = sin(beta) sin(lambda) sin(phi)
    ! + cos(beta) cos(lambda) est inchange)  [SYLVAIN]
    if (s_x == 0.0 .and. s_y == 0.0) then
      phi_deg = 0.0
    else
      phi_deg = atan2(s_x, s_y) * 180.0 / PI
    end if

 end subroutine compute_lambda_phi
