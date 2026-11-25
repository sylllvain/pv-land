!########################################################################
!########################################################################
!#########                                                      #########
!#########   PV-LAND exchange coefficients (version 2)          #########
!#########                                                      #########
!########################################################################
!########################################################################
!
! dragcoef_pv, cuc_pv, cptc_pv: copies of dragcoef, cuc, cptc
! (sfcphy3d.f90) from PV_LAND_final, used only for the PV-LAND interfaces
! (natural ground -> canopy air, canopy air and air above the panels ->
! reference), so that the ARPS surface layer (sfcphy = 2, 4 and outside
! the PV park) is unchanged:
!   - thermal roughness zt = ztz0_in*z0 given as argument, also used in
!     the neutral heat coefficient c_ptneu (instead of c_uneu/prantl0l)
!     and in cptc (instead of the zt(u*) formula of ARPS);
!   - psim of cuc uses (z1+z0)/z0 (z1droup) for x0.
! Converted to single precision (ARPS).  [SYLVAIN]

	SUBROUTINE dragcoef_pv(nx,ny,nz,                                           &
		       zp, ptsfc,qvsfc,soiltyp,roufns,ztz0_in,wdspsfc,pt1,          &
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
	!  PV-LAND
	!  Added the input ztz0_in, ratio of thermal to momentum roughness
	!  length, passed to CUC and CPTC (zt = ztz0_in*roufns in both). It
	!  replaces the constant ztz0 of sfcphycst.inc and the Chen & Dudhia
	!  (2001) thermal roughness formerly used in CPTC, so that the caller
	!  sets zt per interface. The neutral c_ptneu, which sets the bounds
	!  imposed on c_pt, uses the same zt.
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
	!    ztz0_in  Ratio of thermal to momentum roughness length, zt/z0
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
		REAL, INTENT(IN) :: ztz0_in  ! Ratio zt/z0 of thermal to momentum roughness length [-]

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
		REAL :: zt                   ! Thermal roughness length, ztz0_in*roufns [m]
		REAL :: ztdroup              ! (z1+zt)/zt [-]
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

				zt      = ztz0_in*roufns(i,j)
				ztdroup = (z1+zt)/zt

				c_uneu (i,j) = kv/LOG(z1droup)
				! Neutral C_pt with the same thermal roughness zt as CPTC, i.e.
				! the neutral limit of the c_pt that CPTC returns, so that the
				! bounds imposed on c_pt below are relative to its own neutral
				! value. Identical to c_uneu/prantl0l when zt = z0 (ztz0_in = 1).
				c_ptneu(i,j) = kv/(prantl0l*LOG(ztdroup))


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
		!             c_uneu =kv/log((z1+z0)/z0)
		!             c_ptneu=kv/(prandtl_number*log((z1+zt)/zt)), zt = ztz0_in*z0
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
		CALL cuc_pv(nx,ny,nz,1,nx-1,1,ny-1,zp,                                   &
			   roufns,ztz0_in,wdspsfc,ptsfc,pt1,c_uneu, c_u)

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
			!  Change Von Karman constant from 0.4d0 to 0.35d0 and since the drag
			!  coefficient Cdm is proportional to the square of the constant,
			!  the final usflx should be multiplied by (0.35d0/0.40d0)**2 = 0.765625d0
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
		CALL cptc_pv(nx,ny,nz,1,nx-1,1,ny-1,zp,roufns,ztz0_in,wdspsfc,           &
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
			!  Change Von Karman constant from 0.4d0 to 0.387d0 and since the drag
			!  coefficient Cdh is proportional to the square of the constant,
			!  the final Cdh should be multiplied by (0.387d0/0.40d0)**2 = 0.936d0
			!
			!-----------------------------------------------------------------------
			!
			  cdh(i,j) = 0.936*cdh(i,j)
			  !cdq(i,j) = 0.7d0*cdh(i,j)
			  cdq(i,j) = cdh(i,j)

			!      ELSE IF (soilmodel_option == 2) THEN
			!        IF ( soiltyp(i,j) == 13 .AND. cdhwtropt == 1 ) THEN
			!          cdh(i,j) = cdhwtr
			!        ELSE
			!          cdh(i,j) = c_pt(i,j)
			!          cdq(i,j) = 0.7d0*cdh(i,j)
			!        END IF
			!      END IF
			END DO
		END DO

		RETURN
	END SUBROUTINE dragcoef_pv




	!
	!##################################################################
	!##################################################################
	!######                                                      ######
	!######                    SUBROUTINE cuc_pv                    ######
	!######                                                      ######
	!######                     Developed by                     ######
	!######     Center for Analysis and Prediction of Storms     ######
	!######                University of Oklahoma                ######
	!######                                                      ######
	!##################################################################
	!##################################################################
	!

	SUBROUTINE cuc_pv(nx,ny,nz,ibgn,iend,jbgn,jend,zp,roufns,ztz0_in,wspd,     &
		       ptsfc,pt1, c_uneu,c_u)
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
	!  PV-LAND
	!  The thermal roughness length used in the stability parameter is
	!  zt = ztz0_in*roufns, with ztz0_in passed by DRAGCOEF (same zt as
	!  in CPTC).
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
	!    ztz0_in  Ratio of thermal to momentum roughness length, zt/z0
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
	  REAL, INTENT(IN) :: ztz0_in  ! Ratio zt/z0 of thermal to momentum roughness length [-]

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

	  REAL, PARAMETER :: eps_tol = 1.0e-6

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

		  zt = ztz0_in*roufns(i,j)        ! zt/z0 set by the caller (DRAGCOEF)

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
	!      IF (bulkri <= 0.0d0) THEN
	!        c_u(i,j) =kv/(LOG(z1droup)-psim)
	!      ELSE
	!        a=kv/LOG(z1droup)
	!        b=5.0d0
	!        d=5.0d0
	!        c=SQRT(1.0d0+d*bulkri)
	!        c_u(i,j) = a/SQRT(1.0d0+2.0d0*b*bulkri/c)
	!      END IF
	!
	!       zt = (0.4d0*c_u(i,j)/0.000024d0) + 100.0d0
	!
	!       IF (abs(zt) > eps_tol ) THEN
	!         zt = 1.0d0/zt
	!       ELSE IF (abs(zt) <= eps_tol ) THEN
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
		    x0=(1. - c8/z1droup)**0.25

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
	END SUBROUTINE cuc_pv




	!##################################################################
	!##################################################################
	!######                                                      ######
	!######                   SUBROUTINE cptc_pv                    ######
	!######                                                      ######
	!######                     Developed by                     ######
	!######     Center for Analysis and Prediction of Storms     ######
	!######                University of Oklahoma                ######
	!######                                                      ######
	!##################################################################
	!##################################################################
	!

	SUBROUTINE cptc_pv(nx,ny,nz,ibgn,iend,jbgn,jend,                           &
		       zp,roufns,ztz0_in,wspd,ptsfc,pt1, c_ptneu,c_pt,c_u)
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
	!  PV-LAND
	!  The thermal roughness length is zt = ztz0_in*roufns, with ztz0_in
	!  passed by DRAGCOEF (same zt as in CUC). The Chen & Dudhia (2001)
	!  thermal roughness is no longer used, see the comment in the code.
	!  The argument c_u is therefore no longer used here.
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
	!    ztz0_in  Ratio of thermal to momentum roughness length, zt/z0
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
	  REAL, INTENT(IN) :: ztz0_in  ! Ratio zt/z0 of thermal to momentum roughness length [-]

	  REAL :: wspd  (nx,ny)        ! Surface wind speed (m/s)

	  REAL :: ptsfc(nx,ny)         ! Potential temperature at the
		                           ! ground level (K)
	  REAL :: pt1   (nx,ny)        ! Potential temperature at the
		                           ! 1st scalar
		                           ! point above ground level, (K)
	  REAL :: c_ptneu(nx,ny)       ! Normalized temperature scale (K)
		                           ! at neutral state

	  REAL :: c_pt  (nx,ny)        ! Temperature scale (K)

	  REAL :: c_u   (nx,ny)        ! Cu = u*/U [-] (no longer used here)
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

	  REAL, PARAMETER :: eps_tol = 1.0e-6

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

		  zt = ztz0_in*roufns(i,j)        ! zt/z0 set by the caller (DRAGCOEF)

	!---------------------------------------------------------------------
	!     Test of new thermal roughness equation (Chen/Dudhia 2001) (JAB)
	!---------------------------------------------------------------------
	! No longer used: the formula expects u* in m/s but was given
	! c_u = u*/U (dimensionless), and it describes the molecular sublayer
	! over a solid surface, so it does not apply to the air-air interfaces
	! (canopy air -> reference, air above panels -> reference).
	!
	!       zt = (0.4d0*c_u(i,j)/0.000024d0) + 100.0d0
	!
	!       IF (abs(zt) > eps_tol ) THEN
	!         zt = 1.0d0/zt
	!       ELSE IF (abs(zt) <= eps_tol ) THEN
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


	!-------------------------------------------------------------------
	!     Soilmodel_option = 2 (Implicit)
	!--------------------------------------------------------------------

	!     IF (soilmodel_option == 2) THEN
	!
	!     IF (bulkri <= 0.0d0) THEN
	!---------------------------------------------------------------------
	!     Unstable case
	!---------------------------------------------------------------------
	!
	!        stabp = 1.0d0 - ( (15.0d0 * bulkri) /                         &
	!          ( 1.0d0 + 7.5d0 * (10.0d0*kv*kv/(LOG(z1drou)*LOG(ztdrou))) *  &
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
	END SUBROUTINE cptc_pv

