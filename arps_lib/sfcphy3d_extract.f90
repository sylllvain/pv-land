!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE sfcflx (do not edit) ====
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SFCFLX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sfcflx(nx,ny,nz,                                             &
           u,v,w,ptprt,pprt,qv,ptbar,pbar,qvbar,                        &
           wdspsfc,rhosfc,ptsfc,qvsfc, cdm,cdh,cdq,                     &
           usflx,vsflx,ptsflx,qvsflx)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the surface momentum, heat and surface moisture fluxes.
!  These fluxes will be passed into turbulent mixing subroutines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/07/1992.
!
!  MODIFICATION HISTORY:
!
!  2/24/1999 (M. Xue)
!  Simplified codes that calculate the sfc fluxes.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        z component of velocity at a given time level (m/s)
!
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    qvbar    Base state specific humidity (kg/kg)
!
!    wdspsfc  Wind speed at the surface
!    rhosfc   Base-state air density at the surface
!    ptsfc    Potential temperature at the ground level (K)
!    qvsfc    Effective specific humidity at sfc.
!    cdm      Drag coefficient for momentum defined as scalar point
!    cdh      Drag coefficient for heat
!    cdq      Drag coefficient for moisture
!
!  OUTPUT:
!
!    usflx    Surface flux of u-momentum
!    vsflx    Surface flux of v-momentum
!    ptsflx   Surface flux of heat (K*kg/(m**2*s))
!    qvsflx   Surface flux of moisture (K*kg/(m**2*s))
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions
  REAL :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL :: w      (nx,ny,nz) ! Total w-velocity (m/s)

  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: qv     (nx,ny,nz) ! Water vapor specific humidity (kg/kg)

  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: qvbar  (nx,ny,nz) ! Base state specific humidity (kg/kg)

  REAL :: ptsfc  (nx,ny)    ! Potential temperature at the ground level (K)
  REAL :: rhosfc (nx,ny)    ! Surface air density (kg/m**3)
  REAL :: wdspsfc(nx,ny)    ! Wind speed at the surface
  REAL :: qvsfc  (nx,ny)    ! Effective specific humidity at surface

  REAL :: cdm    (nx,ny)    ! Drag coefficient for surface momentum flux
  REAL :: cdh    (nx,ny)    ! Drag coefficient for surface heat flux
  REAL :: cdq    (nx,ny)    ! Drag coefficient for surface moisture flux

  REAL :: usflx  (nx,ny)    ! surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx  (nx,ny)    ! surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx (nx,ny)    ! surface flux of heat (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)    ! surface flux of moisture (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: usuf, vsuf,qvsfc0
  REAL :: tema
  REAL :: gumove,gvmove
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!  Include file globcst.inc passes cdmlnd,cdmwtr, cdhlnd,cdmwtr,
!  cdqlnd,cdmwtr and UMOVE and VMOVE into this subroutine.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( grdtrns == 0 ) THEN
    gumove = 0.0
    gvmove = 0.0
  ELSE
    gumove = umove
    gvmove = vmove
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute the surface flux of u-momentum, i.e. the surface
!  drag in x-direction, according to a simple drag law:
!
!  usflx = rhobar *bar( (w'u') ) = rhobar *Cdm*V*u
!
!  where Cdm is the drag coefficient, V the wind speed at the first
!  scalar level above ground, and u the x-component of wind at the
!  same level.
!
!  usflx is defined one-half dz below the u-point, i.e. at ground level.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=2,nx-1
      usuf = 0.5*(u(i,j,2)+u(i,j,1))+gumove
      usflx(i,j)=0.125*(cdm(i,j)+cdm(i-1,j))*                           &
                 (rhosfc(i,j)+rhosfc(i-1,j))*                           &
                 (wdspsfc(i,j)+wdspsfc(i-1,j))*                         &
                 sign(1.0, usuf)*MAX(abs(usuf),vsfcmin)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the surface flux of v-momentum, i.e. the surface
!  drag in y-direction, according to a simple drag law:
!
!  vsflx = - rhobar *bar( (w'v') ) = rhobar *Cdm*V*v
!
!  where Cdm is the drag coefficient, V the wind speed at the first
!  scalar level above ground, and v the y-component of wind at the
!  same level.
!
!  vsflx is defined one-half dz below the v-point, i.e. at ground level.
!
!-----------------------------------------------------------------------
!
  DO j=2,ny-1
    DO i=1,nx-1
      vsuf = gvmove + 0.5*(v(i,j,2)+v(i,j,1))
      vsflx(i,j)=0.125*(cdm(i,j)+cdm(i,j-1))*                           &
                 (rhosfc(i,j)+rhosfc(i,j-1))*                           &
                 (wdspsfc(i,j)+wdspsfc(i,j-1))*                         &
                 sign(1.0,vsuf)*MAX(abs(vsuf),vsfcmin)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Compute the surface heat flux according to a simple drag law:
!
!  ptsflx = - rhobar *bar( (w' pt') ) = rhobar *Cdh*V*(pt-ptsfc0)
!
!  where Cdh is the drag coefficient, V the wind speed and pt the
!  potential temperatureat at the first scalar level above ground.
!
!  ptsfc0 is defined as the potential temperature at ground leval.
!  ptsfc0 = ptbar at surface
!
!  ptsflx is defined one-half dz below the scalar point, i.e. at
!  ground level.
!
!  Similarly, the vertical moisture flux at the surface is defined as
!
!  qvsflx = - rhobar *bar( (w' qv') ) = rhobar *Cdq*V*(qv-qvsfc0)
!
!  where qv is the specific humidity at the first model scalar level.
!  qvsfc0 is defined as the mixing ratio at the ground and taken as
!  the base state moisture, qvsfc0 = qvbar, at ground level.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      IF ( sfcphy == 1 .OR. sfcphy == 3 ) THEN
        qvsfc0=0.5*(qvbar(i,j,1)+qvbar(i,j,2))
      ELSE
        qvsfc0=qvsfc(i,j)
      END IF
      tema = rhosfc(i,j)*MAX(wdspsfc(i,j),vsfcmin)
      ptsflx(i,j)=cdh(i,j)*tema*(ptprt(i,j,2)+ptbar(i,j,2)-ptsfc(i,j))
      qvsflx(i,j)=cdq(i,j)*tema*(qv(i,j,2)-qvsfc0)
    END DO
  END DO

  RETURN
END SUBROUTINE sfcflx

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE dragcoef (do not edit) ====
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE DRAGCOEF                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dragcoef(nx,ny,nz,                                           &
           zp, ptsfc,qvsfc,soiltyp,roufns,wdspsfc,pt1,                  &
           cdm,cdh,cdq)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the surface momentum, heat and surface moisture fluxes
!  using a stability dependent formulation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/07/1992.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  6/4/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  9/14/1992 (M. Xue)
!  Different drag coefficients defined for momentum, temperature and
!  moisture. Definition of ptsfc and qvsfc changed. qvbar added to
!  the argument list.
!
!  9/28/1992 (M. Xue)
!  Sign correction in all the flux formulations.
!
!  8/28/1993 (M. Xue and Hao jin)
!  Modified the flux formulations considering the terrain effect.
!
!  9/10/1993 (V. Wong, X. Song and N. Lin)
!  Stability dependent formulation used for usflx, vsflx, and ptsflx.
!  The different flux formulations come from subroutines USTARC and
!  PTSTARC.
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  12/7/1994 (Y. Liu)
!  Cleaned up an empty DO loop with a label 360.
!
!  01/28/1995 (V. Wong and X. Song)
!  Add the option of stability and roughness dependent surface layer
!  parameterization for both land and sea surfaces.
!
!  03/02/1995 (Y. Liu)
!  Added the minimum surface total wind speed to guarantee the basic
!  heat and moisture fluxes.
!
!  02/07/96 (V.Wong and X.Song)
!  Changed the calculation of roughness zo over the sea according to
!  Clarke(1970).  The new formulation is insensitive to the value of
!  the depth. of the first model layer above the ground.
!  Added kvwtr to denote the Von Karman constant over the sea.
!  Set a lower limiter, zolimit, for zo, and an upper limiter, z1limit,
!  for depth of the surface layer z1.
!
!  03/17/1997 (Yuhe Liu, Vince Wong and Jing Tian)
!  Added an option for constant cdh (cdq) over water.
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!  02/25/1999 (M. Xue)
!  Isolated from original SFCFLXSD routine. The above listing is for
!  SFCFLXSD.
!
!  01/09/2002 (Y. Wang)
!  Changed the DO loop lower bound from i=2, j=2 to i=1, j=1 in
!  the calculation of drag coefficient for momentum.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!
!    ptsfc    Potential temperature at the ground level (K)
!    qvsfc    Effective specific humidity at sfc
!    soiltyp  Soil type at each point
!    roufns   Surface roughness
!
!    wdspsfc  Surface wind speed (at first level about ground) (m/s)
!    pt1      Potential temperature at the first model level AGL
!
!  OUTPUT:
!
!    cdm      Drag coefficient for surface momentum flux (nondimensional)
!    cdh      Drag coefficient for surface heat flux (nondimensional)
!    cdq      Drag coefficient for surface moisture flux (nondimensional)
!
!  Local work arrays:
!
!    c_u      Cu
!    c_pt     Cpt
!    c_uneu   Cu for neutral stratification
!    c_ptneu  Cpt for neutral stratification
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  REAL :: ptsfc (nx,ny)        ! Potential temperature at the ground level (K)
  REAL :: qvsfc (nx,ny)        ! Effective surface specific humidity

  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point
  REAL :: roufns(nx,ny)        ! Surface roughness length

  REAL :: wdspsfc(nx,ny)       ! Wind speed at the surface
  REAL :: pt1   (nx,ny)        ! Potential temperature at the first model level AGL

  REAL :: cdm   (nx,ny)        ! Drag coefficient for surface momentum flux
  REAL :: cdh   (nx,ny)        ! Drag coefficient for surface heat flux
  REAL :: cdq   (nx,ny)        ! Drag coefficient for surface moisture flux

  REAL :: c_u   (nx,ny)        ! Cu
  REAL :: c_pt  (nx,ny)        ! Cpt
  REAL :: c_uneu(nx,ny)        ! Cu for neutral stratification
  REAL :: c_ptneu(nx,ny)       ! Cpt for neutral stratification
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: check
  REAL :: z1,z1drou,z1droup
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!  Include file globcst.inc passes cdq (moisture transfer coeffcient)
!  and UMOVE and VMOVE into this subroutine.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny-1
    DO i=1,nx-1
!
! Calculate C_u and C_pt at the neutral state
!
      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)
      z1drou = z1/roufns(i,j)

      z1droup = (z1+roufns(i,j))/roufns(i,j)

      c_uneu (i,j) = kv/LOG(z1droup)
      c_ptneu(i,j) = c_uneu(i,j) /prantl0l


    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  A general procedure for computing the drag coefficients:
!
!  1.  Obtain the following quantities:
!      a. The first scalar level above ground level (z1)
!      b. Potential temperature at the height of surface roughness
!         (roufns) - given.
!      c. Potential temperature (pt1) at the first scalar level above
!         ground.
!      d. C_u and C_pt at neutral state in the land case from
!             c_uneu =kv/log10(z1/roughness)
!             c_ptneu=c_uneu/prandtl_number
!         In the water case, call subroutines CUNEUWTR and CPTNEUWTR.
!
!  2.  Given the information from (1), compute the Bulk Richardson
!      number, defined as
!
!         BulkRi = (g/ptsfc) * ( (pt1-ptsfc)*(z1-roufns) / V**2 )
!      where V is the wind speed at the surface.
!
!  3.  After calculation of the BulkRi, then:
!      a. Given the BulkRi and the information in (1), compute C_u at
!         both land and water cases (see subroutines CUC and CUCWTR).
!
!      b. From C_u and C_pt, compute draf coefficients cdh and cdq.
!
!-----------------------------------------------------------------------
!
!  Calculate drag coefficient (defined at scalar point) for momentum
!  fluxes.
!
!-----------------------------------------------------------------------
!
  CALL cuc(nx,ny,nz,1,nx-1,1,ny-1,zp,                                   &
           roufns,wdspsfc,ptsfc,pt1,c_uneu, c_u)

  IF (wtrexist == 1) THEN

    CALL cuneuwtr(nx,ny,nz,1,nx-1,1,ny-1,soiltyp,wdspsfc,               &
                  c_uneu)

    CALL cucwtr  (nx,ny,nz,1,nx-1,1,ny-1,zp,soiltyp,wdspsfc,            &
                  ptsfc,pt1,c_uneu, c_u)

  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      check=ptsfc(i,j)-pt1(i,j)
!
!    Put lower and upper limits on C_u
!
      IF(check >= 0.0) THEN  ! Unstable case
        c_u(i,j)=MIN( c_u(i,j), 2.0*c_uneu(i,j))
      ELSE                   ! Stable case
        c_u(i,j)=MAX( c_u(i,j), blimit*c_uneu(i,j))
      END IF
!
!-----------------------------------------------------------------------
!
!  c_dm=c_U**2
!
!  Change Von Karman constant from 0.4 to 0.35 and since the drag
!  coefficient Cdm is proportional to the square of the constant,
!  the final usflx should be multiplied by (0.35/0.40)**2 = 0.765625
!
!-----------------------------------------------------------------------
!
      cdm(i,j) = 0.765625 * c_u(i,j) * c_u(i,j)

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  From C_u and C_pt, compute C_dh = C_u * C_pt and pt flux,
!
!-----------------------------------------------------------------------
!
  CALL cptc(nx,ny,nz,1,nx-1,1,ny-1,zp,roufns,wdspsfc,                   &
            ptsfc,pt1,c_ptneu,c_pt,c_u)

  IF ( wtrexist == 1 .AND. cdhwtropt == 0 ) THEN

    CALL cptneuwtr(nx,ny,nz,1,nx-1,1,ny-1,soiltyp,wdspsfc,c_uneu,       &
                   c_ptneu)

    CALL cptcwtr(nx,ny,nz,1,nx-1,1,ny-1,zp,soiltyp,wdspsfc,             &
                 ptsfc,pt1,c_ptneu, c_pt)

  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      check=ptsfc(i,j)-pt1(i,j)
!
!    Impose lower and upper limits on C_u and C_pt.
!

!      IF (soilmodel_option == 1) THEN

      IF ( soiltyp(i,j) == 13 .AND. cdhwtropt == 1 ) THEN
        cdh(i,j) = cdhwtr
      ELSE
        IF(check >= 0) THEN
          c_pt(i,j) = MIN(c_pt(i,j),3.3333*c_ptneu(i,j) )
        ELSE
          c_pt(i,j) = MAX(c_pt(i,j),blimit*c_ptneu(i,j) )
        END IF

        cdh(i,j) = c_u(i,j)*c_pt(i,j)
      END IF
!
!-----------------------------------------------------------------------
!
!  Change Von Karman constant from 0.4 to 0.35 and since the drag
!  coefficient Cdh is proportional to the square of the constant,
!  the final Cdh should be multiplied by (0.35/0.40)**2 = 0.765625
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! MODIF SYLVAIN FOR PVLAND
!-----------------------------------------------------------------------
      !cdh(i,j) = 0.765625*cdh(i,j)
      !cdq(i,j) = 0.7*cdh(i,j)
	  cdh(i,j) = 0.936*cdh(i,j)
	  cdq(i,j) = cdh(i,j)
!-----------------------------------------------------------------------
! FIN MODIF SYLVAIN FOR PVLAND
!-----------------------------------------------------------------------

!      ELSE IF (soilmodel_option == 2) THEN
!        IF ( soiltyp(i,j) == 13 .AND. cdhwtropt == 1 ) THEN
!          cdh(i,j) = cdhwtr
!        ELSE
!          cdh(i,j) = c_pt(i,j)
!          cdq(i,j) = 0.7*cdh(i,j)
!        END IF
!      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE dragcoef

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE cuc (do not edit) ====
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    SUBROUTINE CUC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cuc(nx,ny,nz,ibgn,iend,jbgn,jend,zp,roufns,wspd,ptsfc,       &
           pt1, c_uneu,c_u)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute C_u (friction velocity / wind speed). The quantity C_U is used by
!  the subroutine SFCFLX to obtain surface fluxes for both unstable
!  and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong, X. Song and N. Lin
!  9/10/1993
!
!  For the unstable case (Bulk Richardson number bulkri < 0 ), the
!  formulation is based on the paper "On the Analytical Solutions of
!  Flux-Profile relationships for the Atmospheric Surface Layer" by
!  D.W. Byun in J. of Applied Meteorlogy, July 1990, pp. 652-657.
!
!  For the stable case, the formulation is based on the
!  paper "A Short History of the Operational PBL - Parameterization
!  at ECMWF" by J.F.Louis, M. Tiedtke and J.F. Geleyn in " Workshop
!  on Planetary Boundary Layer Parameterization", a publication
!  by the European Centre for Medium Range Weather Forecasts,
!  25-27 Nov. 1981.
!
!  MODIFICATION HISTORY:
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  2/27/95 (V. Wong and X. Song)
!
!  02/07/96 (V.Wong and X.Song)
!  Set an upper limiter, z1limit, for depth of the surface layer z1.
!
!  05/01/97 (V. Wong and X. Tan)
!  Changed the computation of stabp, use the formulation similar to
!  one used by Bynn(1990), the momentum and thermal roughness
!  lengths are different.
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    roufns   Surface roughness
!
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!    c_uneu   Friction velocity at neutral state
!
!  OUTPUT:
!
!    ustar    Friction velocity
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate
                               ! defined at w-point of staggered grid.
  REAL :: roufns(nx,ny)        ! Surface roughness length

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground
                               ! level (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the
                               ! 1st scalar
                               ! point above ground level, (K)
  REAL :: c_uneu(nx,ny)      ! Frictional velocity (m/s) at neutral

  REAL :: c_u(nx,ny)        ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: dz0                  ! z1-roufns, where roufns is defined as
                               ! surface roughness length
  REAL :: bulkri               ! Bulk Richardson number

  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)
  REAL :: x1,x0,psim           ! Intermediate parameters needed
  REAL :: z1drou,qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! During  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb

  REAL :: zt                   ! Thermal roughness length
  REAL :: dzt
  REAL :: ztdrou

  REAL :: z1droup, ztdroup

  REAL, PARAMETER :: epsilon = 1.0E-6

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Following Byun (1990).
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)

      zt = ztz0*roufns(i,j)           ! Original Zt formulation

!---------------------------------------------------------------------
!     Test of new thermal roughness equation (Chen/Dudhia 2001) (JAB)
!---------------------------------------------------------------------
! The following was commented out temporarily because variable psim
! is used before it is initilized. It need to be taken care of later.
!
!      dz0 = z1-roufns(i,j)
!      z1droup = (z1+roufns(i,j))/roufns(i,j)
!
!      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
!               (wspd(i,j)*wspd(i,j))
!
!      IF (bulkri <= 0.0) THEN
!        c_u(i,j) =kv/(LOG(z1droup)-psim)
!      ELSE
!        a=kv/LOG(z1droup)
!        b=5.0
!        d=5.0
!        c=SQRT(1.0+d*bulkri)
!        c_u(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)
!      END IF
!
!       zt = (0.4*c_u(i,j)/0.000024) + 100.0
!
!       IF (abs(zt) > epsilon ) THEN
!         zt = 1.0/zt
!       ELSE IF (abs(zt) <= epsilon ) THEN
!         zt = ztz0*roufns(i,j)
!       END IF
!
!----------------------------------------------------------------------

      dzt = z1-zt
      ztdrou = z1/zt
      ztdroup = (z1+zt)/zt

      dz0 = z1-roufns(i,j)
      z1drou = z1/roufns(i,j)
      z1droup = (z1+roufns(i,j))/roufns(i,j)

      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
               (wspd(i,j)*wspd(i,j))

      IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: See equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        bulkri = MAX (bulkri,-10.0)

        sb =bulkri/prantl0l

        qb=oned9*(c1l+c2l*sb*sb)
        pb=oned54*(c3l+c4l*sb*sb)

        qb3pb2=qb**3-pb*pb
        c7 = (z1*dzt*LOG(z1droup)*LOG(z1droup))/(dz0*dz0*LOG(ztdroup))

        IF( qb3pb2 >= 0.0 ) THEN

          sqrtqb = SQRT(qb)
          tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

          thetab=ACOS(tempan)
          stabp =c7*(-2.0*sqrtqb*COS(thetab/3.0)+c5l)

        ELSE

          tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
          stabp =c7*(-(tb+qb/tb)+c5l)

        END IF
!
!-----------------------------------------------------------------------
!
!  According to equation (14) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        c8=gammaml*stabp
        x1=(1. - c8)**0.25
        x0=(1. - c8/z1drou)**0.25

        psim=2.0*LOG((1.0+x1)/(1.0+x0))+LOG((1.+x1*x1)/(1.+x0*x0))-     &
             2.0*ATAN(x1)+2.0*ATAN(x0)

!
!-----------------------------------------------------------------------
!
!  Compute C_u via equation (10) in Byun (1981).
!
!-----------------------------------------------------------------------
!
        c_u(i,j) =kv/(LOG(z1droup)-psim)

      ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case:
!
!-----------------------------------------------------------------------
!
        a=kv/LOG(z1droup)
        b=5.0
        d=5.0
        c=SQRT(1.0+d*bulkri)

        c_u(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)

      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cuc

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE cptc (do not edit) ====
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE CPTC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cptc(nx,ny,nz,ibgn,iend,jbgn,jend,                           &
           zp,roufns,wspd,ptsfc,pt1, c_ptneu,c_pt,c_u)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute c_pt (a nondimensional temperature scale)
!
!  The quantity c_pt is used by the subroutine SFCFLX to obtain
!  surface fluxes for both the unstable and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong, X. Song and N. Lin
!  9/10/1993
!
!  For the unstable case (Bulk Richardson number bulkri < 0 ), the
!  formulation is based on the paper "On the Analytical Solutions of
!  Flux-Profile relationships for the Atmospheric Surface Layer" by
!  D.W. Byun in J. of Applied Meteorlogy, July 1990, pp. 652-657.
!
!  For stable case, the formulation is based on the paper "A Short
!  History of the Operational PBL - Parameterization at ECMWF" by
!  J.F.Louis, M. Tiedtke and J.F. Geleyn in "Workshop on Planetary
!  Boundary Layer Parameterization", a publication by European Centre
!  for Medium Range Weather Forecasts, 25-27 Nov. 1981.
!
!  MODIFICATION HISTORY:
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  2/27/95 (V. Wong and X. Song)
!
!  02/07/96 (V.Wong and X.Song)
!  Set an upper limiter, z1limit, for depth of the surface layer z1.
!
!  05/01/97 (V. Wong and X. Tan)
!  Changed the computation of temperature scale
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    roufns   Surface roughness
!
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!    c_ptneu  Temperature scale (K) at neutral state, defined by
!             surface heat flux / friction velocity
!
!  OUTPUT:
!
!    c_pt     Temperature scale (K), defined by
!             surface heat flux / friction velocity
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins
  INTEGER :: iend              ! i-index where evaluation ends
  INTEGER :: jbgn              ! j-index where evaluation begins
  INTEGER :: jend              ! j-index where evaluation ends

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate
                               ! defined at w-point of staggered grid.
  REAL :: roufns(nx,ny)        ! Surface roughness length

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the
                               ! ground level (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the
                               ! 1st scalar
                               ! point above ground level, (K)
  REAL :: c_ptneu(nx,ny)       ! Normalized temperature scale (K)
                               ! at neutral state

  REAL :: c_pt  (nx,ny)        ! Temperature scale (K)

  REAL :: c_u   (nx,ny)        ! Friction velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                  ! The height of 1st scalar point above
                              ! ground level (m)
  REAL :: dz0                 ! z1-roufns, where roufns is defined as
                              ! surface roughness length
  REAL :: bulkri              ! Bulk Richardson number

  REAL :: stabp               ! Monin-Obukhov STABility Parameter
                              ! (zeta)

  REAL :: y1,y0,psih          ! Intermediate parameters needed
  REAL :: z1drou,qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb  ! During  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb

  REAL :: zt                  ! Thermal roughness length
  REAL :: dzt,ztdrou

  REAL :: z1droup,ztdroup

  REAL, PARAMETER :: epsilon = 1.0E-6

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Following Byun (1990).
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)

!      zt = ztz0*roufns(i,j)

!---------------------------------------------------------------------
!     Test of new thermal roughness equation (Chen/Dudhia 2001) (JAB)
!---------------------------------------------------------------------
       zt = (0.4*c_u(i,j)/0.000024) + 100.0

       IF (abs(zt) > epsilon ) THEN
         zt = 1.0/zt
       ELSE IF (abs(zt) <= epsilon ) THEN
         zt = ztz0*roufns(i,j)
       END IF

!----------------------------------------------------------------------

      dzt = z1-zt
      ztdrou = z1/zt
      ztdroup = (z1+zt)/zt

      dz0 = z1-roufns(i,j)
      z1drou = z1/roufns(i,j)
      z1droup = (z1+roufns(i,j))/roufns(i,j)

      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
               (wspd(i,j)*wspd(i,j))


!-------------------------------------------------------------------
!     Soilmodel_option = 2 (Implicit)
!--------------------------------------------------------------------

!     IF (soilmodel_option == 2) THEN
!
!     IF (bulkri <= 0.0) THEN
!---------------------------------------------------------------------
!     Unstable case
!---------------------------------------------------------------------
!
!        stabp = 1.0 - ( (15.0 * bulkri) /                         &
!          ( 1.0 + 7.5 * (10.0*kv*kv/(LOG(z1drou)*LOG(ztdrou))) *  &
!          sqrt(-bulkri * z1drou) ) )
!
!     ELSE
!---------------------------------------------------------------------
!     Stable case
!---------------------------------------------------------------------
!
!        stabp = EXP(-bulkri)
!
!     END IF
!
!     c_pt(i,j) = kv*kv*stabp / (LOG(z1drou)*LOG(ztdrou))
!
!     END IF      !soilmodel_option = 2
!-------------------------------------------------------------------



!-------------------------------------------------------------------
!      Soilmodel_option = 1  (Force-Restore)
!------------------------------------------------------------------

!      IF (soilmodel_option == 1) THEN


      IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: See equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        bulkri = MAX (bulkri,-10.0)

        sb =bulkri/prantl0l

        qb=oned9*(c1l+c2l*sb*sb)
        pb=oned54*(c3l+c4l*sb*sb)

        qb3pb2=qb**3-pb*pb
        c7 = (z1*dzt*LOG(z1droup)*LOG(z1droup))/(dz0*dz0*LOG(ztdroup))

        IF( qb3pb2 >= 0.0 ) THEN

          sqrtqb = SQRT(qb)
          tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

          thetab=ACOS(tempan)
          stabp =c7*(-2.0*SQRT(qb)*COS(thetab/3.0)+c5l)

        ELSE

          tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
          stabp =c7*(-(tb+qb/tb)+c5l)

        END IF
!
!-----------------------------------------------------------------------
!
!  According to equation (15) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        c8=gammahl*stabp
        y1=SQRT(1.0 - c8)
        y0=SQRT(1.0 - c8/ztdrou)

        psih=2.0*LOG((y1+1.0)/(y0+1.0))
!
!-----------------------------------------------------------------------
!
!  Compute c_pt via equation (11) in Byun (1981).
!
!-----------------------------------------------------------------------
!
        c_pt(i,j)=kv / (prantl0l*(LOG(ztdroup)-psih))

      ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case: See Louis et al (1981).
!
!-----------------------------------------------------------------------
!
        a=kv/LOG(ztdroup)
        b=5.0
        d=5.0
        c=SQRT(1.0+d*bulkri)

        c_pt(i,j) = SQRT(1.0+2.0*b*bulkri/c)
        c_pt(i,j) = a*c_pt(i,j)/(prantl0l*(1.0+3.0*b*bulkri*c))

      END IF
!    END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cptc

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE cucwtr (do not edit) ====
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CUCWTR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cucwtr(nx,ny,nz,ibgn,iend,jbgn,jend,zp,soiltyp,              &
           wspd,ptsfc,pt1,c_uwtrneu,c_uwtr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute C_u (friction velocity / wind speed) at sea case. The quantity
!  C_uwtr is used by the subroutine SFCFLX to obtain surface fluxes for both
!  unstable and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!  MODIFICATION HISTORY:
!
!  2/27/95 (V.W. and X.S.)
!
!  1/12/96 (V.W. and X.S.)
!  Changed the calculation related to zo.
!  Added kvwtr to denote the Von Karman constant over the sea;
!  Set a lower limiter for zo, zolimit, and an upper limiter for z1, z1limit.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!  OUTPUT:
!
!    c_uwtr   Friction velocity
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point.

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground level
                               ! (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the 1st scalar
                               ! point above ground level, (K)
  REAL :: c_uwtrneu(nx,ny)        ! Frictional velocity (m/s)
  REAL :: c_uwtr(nx,ny)        ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: zo                   ! Defined as surface momentum roughness
                               ! length
  REAL :: zt                   ! Defined as surface heat transfer
                               ! roughness length
  REAL :: dzo                  ! z1-zo
  REAL :: dzt                  ! z1-zt
  REAL :: z1dzo                ! z1/zo
  REAL :: z1dzt                ! z1/zt
  REAL :: xcdn                 ! Cdn (sea)
  REAL :: xcdh                 ! Hot-wired value of Cdh (sea)
  REAL :: bulkri               ! Bulk Richardson number

  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)

  REAL :: x1,psim           ! Intermediate parameters needed
  REAL :: qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! during  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xcdh = 1.1E-3

  DO j=jbgn,jend
    DO i=ibgn,iend

      IF (soiltyp (i,j) == 13) THEN

        xcdn = (0.4+0.079*wspd(i,j)) * 1.e-3
        z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
        z1  = MIN(z1,z1limit)
!
! calculate zo and zt
!
        zo =0.032*xcdn*wspd(i,j)*wspd(i,j)/9.8
        zo = MAX(zo,zolimit)
        zt = z1 * EXP( -kvwtr*SQRT(xcdn)/(prantl0w*xcdh) )

        dzo = z1-zo
        z1dzo = z1/zo
        dzt = z1-zt
        z1dzt = z1/zt

        bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dzo*dzo/          &
                 (dzt*wspd(i,j)*wspd(i,j))

        IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: A modified formulation, which is similar to
!  equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
          bulkri = MAX (bulkri,-10.0)

          sb =bulkri/prantl0w

          qb=oned9*(c1w+c2w*sb*sb)
          pb=oned54*(c3w+c4w*sb*sb)

          qb3pb2=qb**3-pb*pb
          c7 = (z1*dzt/(dzo*dzo))*(ALOG(z1dzo)*ALOG(z1dzo)/ALOG(z1dzt))

          IF( qb3pb2 >= 0.0 ) THEN

            sqrtqb = SQRT(qb)
            tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

            thetab=ACOS(tempan)
            stabp =c7*(-2.0*SQRT(qb)*COS(thetab/3.0)+c5w)

          ELSE

            tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
            stabp =c7*(-(tb+qb/tb)+c5w)

          END IF
!
!-----------------------------------------------------------------------
!
!  According to a modified equation, which is similar to equation (14)
!  in Byun (1990).
!
!-----------------------------------------------------------------------
!
          c8=gammamw*stabp
          x1=(1. - c8)**0.25

          psim=2.0*ALOG(0.5*(1.0+x1))+ALOG(0.5*(1.+x1*x1))-             &
               2.0*ATAN(x1)+ASIN(1.)

!
!-----------------------------------------------------------------------
!
!  Compute c_uwtr via equation (10) in Byun (1981).
!
!-----------------------------------------------------------------------
!
          c_uwtr(i,j) =kv/(LOG(z1dzo)-psim)

        ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case: See Louis et al (1981).
!
!-----------------------------------------------------------------------
!
          a=kv/LOG(z1dzo)
          b=5.0
          d=5.0
          c=SQRT(1.0+d*bulkri)

          c_uwtr(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)

        END IF

      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cucwtr

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE cptcwtr (do not edit) ====

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CPTCWTR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cptcwtr(nx,ny,nz,ibgn,iend,jbgn,jend,zp,soiltyp,wspd,        &
           ptsfc,pt1,c_ptwtrneu,c_ptwtr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute c_ptwtr (the product of ustar and ptstar) at sea case.
!  Note:  a temperature scale defined as surface heat flux divided
!  by the friction velocity.
!
!  The quantity c_ptwtr is used by the subroutine SFCFLX to obtain
!  surface fluxes for both the unstable and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!  MODIFICATION HISTORY:
!
!  2/27/95 (V.W. and X.S.)
!
!  1/12/96 (V.W. and X.S.)
!  Changed the calculation related to zo over the sea.
!  Added kvwtr to denote the Von Karman constant over the sea;
!  Set a lower limiter for zo, zolimit, and an upper limiter for z1, z1limit.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    iend     i-index where evaluation ends.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!
!  OUTPUT:
!
!    c_ptwtr  The product of ustar and ptstar. ptstar is temperature
!             scale (K), defined by surface heat flux / friction velocity
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point.

  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: ptsfc(nx,ny)         ! Potential temperature at the ground level
                               ! (K)
  REAL :: pt1   (nx,ny)        ! Potential temperature at the 1st scalar
                               ! point above ground level, (K)

  REAL :: c_ptwtrneu(nx,ny)       ! Product of ustar and ptstar
  REAL :: c_ptwtr(nx,ny)       ! Product of ustar and ptstar
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: zo                   ! Defined as surface momentum roughness
                               ! length
  REAL :: zt                   ! Defined as surface heat transfer
                               ! roughness length
  REAL :: dzo                  ! z1-zo
  REAL :: dzt                  ! z1-zt
  REAL :: z1dzo                ! z1/zo
  REAL :: z1dzt                ! z1/zt
  REAL :: xcdn                 ! Cdn (sea)
  REAL :: xcdh                 ! Hot-wired value of Cdh (sea)

  REAL :: bulkri               ! Bulk Richardson number
  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)

  REAL :: y1,y0,psih           ! Intermediate parameters needed
  REAL :: qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! during  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xcdh = 1.1E-3

  DO j=jbgn,jend
    DO i=ibgn,iend
      IF ( soiltyp(i,j) == 13) THEN

        xcdn = (0.4+0.07*wspd(i,j)) * 1.e-3
        z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
        z1  = MIN(z1,z1limit)
!
! calculate zo and zt
!
        zo =0.032*xcdn*wspd(i,j)*wspd(i,j)/9.8
        zo = MAX(zo,zolimit)
        zt = z1 * EXP( -kv*SQRT(xcdn)/(prantl0w*xcdh) )

        dzo = z1-zo
        z1dzo = z1/zo
        dzt = z1-zt
        z1dzt = z1/zt

        bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dzo*dzo/          &
                 (dzt*wspd(i,j)*wspd(i,j))

        IF (bulkri <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: A modified formulation, which is similar to
!  equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
          bulkri = MAX (bulkri,-10.0)

          sb =bulkri/prantl0w

          qb=oned9*(c1w+c2w*sb*sb)
          pb=oned54*(c3w+c4w*sb*sb)

          qb3pb2=qb**3-pb*pb
          c7 = (z1*dzt/(dzo*dzo))*(ALOG(z1dzo)*ALOG(z1dzo)/ALOG(z1dzt))

          IF( qb3pb2 >= 0.0 ) THEN

            sqrtqb = SQRT(qb)
            tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

            thetab=ACOS(tempan)
            stabp =c7*(-2.0*SQRT(qb)*COS(thetab/3.0)+c5w)

          ELSE

            tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
            stabp =c7*(-(tb+qb/tb)+c5w)

          END IF
!
!-----------------------------------------------------------------------
!
!  According to a modified equation, which is similar to equation (14)
!  in Byun (1990).
!
!-----------------------------------------------------------------------
!
          c8=gammamw*stabp
          y1=SQRT(1.0 - c8)
          y0=SQRT(1.0 - c8/z1dzt)

          psih=2.0*ALOG((y1+1.0)/(y0+1.0))
!
!-----------------------------------------------------------------------
!
!  Compute ptstar via equation (11) in Byun (1981).
!
!-----------------------------------------------------------------------
!
          c_ptwtr(i,j)=kv / (prantl0w*(ALOG(z1dzt)-psih))

        ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case: With the modified formulation in Louis et al (1981).
!
!-----------------------------------------------------------------------
!
          a=kv*kv/(prantl0w*LOG(z1dzo)*LOG(z1dzt))
          b=5.0
          d=5.0
          c=SQRT(1.0+d*bulkri)

          c_ptwtr(i,j) = SQRT(1.0+2.0*b*bulkri/c)
          c_ptwtr(i,j) = a*c_ptwtr(i,j)/(prantl0l*(1.0+3.0*b*bulkri*c))

        END IF
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cptcwtr

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE cuneuwtr (do not edit) ====

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CUNEUWTR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cuneuwtr(nx,ny,nz,ibgn,iend,jbgn,jend,soiltyp,               &
           wspd, c_uneu)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute c_uneu (friction velocity) at the neutral state. The quantity
!  c_uneu is used by the subroutine SFCFLX to obtain surface fluxes.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    soiltyp  Soil type at each point
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!
!
!  OUTPUT:
!
!    c_uneu   Friction velocity at sea case.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  INTEGER :: soiltyp(nx,ny)    ! Soil type at each point.
  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)


  REAL :: c_uneu   (nx,ny)     ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      IF (soiltyp(i,j) == 13) THEN
        c_uneu(i,j) = SQRT ((0.4+0.079*wspd(i,j)) * 1.e-3)
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cuneuwtr

!==== Extracted from ARPS sfcphy3d.f90 : SUBROUTINE cptneuwtr (do not edit) ====
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CPTNEUWTR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cptneuwtr(nx,ny,nz,ibgn,iend,jbgn,jend,                      &
           soiltyp,wspd,c_uneu,c_ptwtrneu)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong and X. Song
!  8/04/1994
!
!  MODIFICATION:
!
!  03/17/1997 (Yuhe Liu and Vince Wong)
!  Fixed a bug. Variable c_ptwtrneu was calculated in its inverse.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins
!    iend     i-index where evaluation ends
!    jbgn     j-index where evaluation begins
!    jend     j-index where evaluation ends
!
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!
!  OUTPUT:
!
!    c_ptwtrneu   Temperature scale (K), defined by
!             surface heat flux / friction velocity
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER :: ibgn              ! i-index where evaluation begins.
  INTEGER :: iend              ! i-index where evaluation ends.
  INTEGER :: jbgn              ! j-index where evaluation begins.
  INTEGER :: jend              ! j-index where evaluation ends.

  INTEGER :: soiltyp (nx,ny)
  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL :: c_uneu (nx,ny)       !
  REAL :: c_ptwtrneu(nx,ny)    ! Temperature scale (K)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: xcdh
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xcdh = 1.14E-3

  DO j=jbgn,jend
    DO i=ibgn,iend

      IF (soiltyp(i,j) == 13) THEN
        c_ptwtrneu(i,j) = xcdh/c_uneu(i,j)
      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cptneuwtr

