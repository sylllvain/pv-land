!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE PVLAND_SFCPHY_OFFLINE         ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE pvland_sfcphy_offline(nx,ny,nz,nzsoil,nstyps,               &
           u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,zp,prcrate,            &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           radsw,rnflx,radlwin,psfc,rhosfc,gumove,                      &
           cdm,cdh,cdq,usflx,vsflx,ptsflx,qvsflx,                       &
           S_diffuse_in,lambda_in,phi_in,L_obk_in,                      &
           uref_diag,cdh_can_diag,cdh_PVa_diag)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  One PV-LAND surface time step (dtsfc) for the offline model.
!
!  The executable part is a COPY of the PV-LAND block of ARPS
!  (sfcphy3d.f90, SUBROUTINE sfcphysics, "ELSE IF ( sfcphy == 5 )"),
!  generated automatically. Only the lines tagged [OFFLINE] differ:
!  the offline model is forced by MEASURED values instead of
!    - split_sw            -> measured diffuse SW     (S_diffuse_in)
!    - compute_lambda_phi  -> measured lambda, phi    (lambda_in, phi_in, rad)
!    - L_obk from ptsflx   -> measured Obukhov length (L_obk_in)
!  and the original ARPS lines are kept commented out just above them.
!  Three lines tagged [OFFLINE DIAG] save local diagnostics (uref,
!  cdh_can, cdh_PVa) for the offline output files.
!
!  Any change of the PV-LAND physics must be done in ARPS (or ported
!  to ARPS) and this file regenerated, see README / Makefile.
!
!-----------------------------------------------------------------------
!
  USE module_pvland_in
  USE module_pvland_out
  USE rad_pvland
  USE rad_split

  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nzsoil,nstyps

  REAL :: u     (nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
  REAL :: ptprt (nx,ny,nz), pprt(nx,ny,nz), qv(nx,ny,nz)
  REAL :: ptbar (nx,ny,nz), pbar(nx,ny,nz), qvbar(nx,ny,nz)
  REAL :: zp    (nx,ny,nz)
  REAL :: prcrate(nx,ny)

  INTEGER :: soiltyp (nx,ny,nstyps), vegtyp(nx,ny)
  REAL    :: stypfrct(nx,ny,nstyps)
  REAL    :: lai(nx,ny), roufns(nx,ny), veg(nx,ny)

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps)
  REAL :: wetcanp(nx,ny,0:nstyps)
  REAL :: snowdpth(nx,ny)
  REAL :: qvsfc  (nx,ny,0:nstyps)

  REAL :: radsw  (nx,ny)    ! Incoming solar radiation
  REAL :: rnflx  (nx,ny)    ! Net radiation (used outside the park only)
  REAL :: radlwin(nx,ny)    ! Incoming longwave radiation
  REAL :: psfc   (nx,ny)    ! Surface pressure (Pascal)
  REAL :: rhosfc (nx,ny)    ! Air density at the surface
  REAL :: gumove            ! Domain translation speed (0 offline)

  REAL :: cdm(nx,ny), cdh(nx,ny), cdq(nx,ny)
  REAL :: usflx(nx,ny), vsflx(nx,ny), ptsflx(nx,ny), qvsflx(nx,ny)

  REAL :: S_diffuse_in(nx,ny)   ! [OFFLINE] measured diffuse SW (W/m2)
  REAL :: lambda_in(nx,ny)      ! [OFFLINE] measured lambda, park frame (rad)
  REAL :: phi_in(nx,ny)         ! [OFFLINE] measured phi, park frame (rad)
  REAL :: L_obk_in(nx,ny)       ! [OFFLINE] measured Obukhov length (m)

  REAL :: uref_diag(nx,ny), cdh_can_diag(nx,ny), cdh_PVa_diag(nx,ny) ! [OFFLINE] diagnostics
!
!-----------------------------------------------------------------------
!  Local variables of sfcphysics used by the PV-LAND block
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: pi = 3.14159265359

  REAL :: pt1    (nx,ny)    ! Potential temperature at 1st model level AGL
  REAL :: tair   (nx,ny)    ! Temperature of surface air
  REAL :: qvair  (nx,ny)    ! Specific humidity of surface air

  REAL :: z_ref,usuf,vsuf,Uplus
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdm_u,cdm_v
  REAL, ALLOCATABLE, DIMENSION (:,:) :: z0_PV

  REAL, ALLOCATABLE, DIMENSION (:,:) :: uref,T_nat,T_soil,qsoil_surf,&
	qsoil_deep,S_diffuse,S_direct,wrmax
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_nat,cdq_nat,cdm_nat
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_PVa,cdq_PVa,cdm_PVa
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_PVb,cdq_PVb,cdm_PVb
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_PVf,cdq_PVf,cdm_PVf
  REAL, ALLOCATABLE, DIMENSION (:,:) :: cdh_can,cdq_can,cdm_can
  REAL, ALLOCATABLE, DIMENSION (:,:) :: L_obk

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

  INTEGER :: i, j, k, is
  REAL :: tema
  REAL :: nsfcstsave

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
!=== BEGIN copy of ARPS sfcphy3d.f90, sfcphy == 5 block ================
!@@ARPS_SFCPHY5_BLOCK@@
!=== END copy of ARPS sfcphy3d.f90, sfcphy == 5 block ==================

  RETURN
END SUBROUTINE pvland_sfcphy_offline
