!########################################################################
!########################################################################
!#########                                                      #########
!#########   PV-LAND soil model (version 2, PV_LAND_final)      #########
!#########                                                      #########
!########################################################################
!########################################################################
!
! soilebm_frc_pv and evapflx_pv: copies of soilebm_frc and evapflx
! (soilebm3d.f90) with the PV-LAND version 2 physics, used only by
! PV-LAND (sfcphy = 5), so that the ARPS soil model (sfcphy = 3, 4) is
! unchanged:
!   - evaporation towards the canopy-air humidity q_can (conductance gcan
!     canopy air -> reference; gcan = 0 gives q_can = qvair),
!   - water reaching the soil (1-veg)P + runoff of the interception
!     reservoir, drainage of the deep reservoir (pv_drain_opt, pv_c3_drain),
!   - moisture threshold pv_wref instead of wfc, no factor 1.03 in the
!     deep-soil restore term,
!   - soil and vegetation parameters from the namelist &pvland (pv_cgv,
!     pv_wsat, pv_wfc, pv_wwlt, pv_rsmin; negative = ARPS table).
! Converted to single precision (ARPS).  [SYLVAIN]



!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE soilebm_frc_pv               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE soilebm_frc_pv(nx,ny,nz,soiltyp,vegtyp,lai,veg,           &
     tsfc,tdeep,wetsfc,wetdp,wetcanp,snowdpth,qvsfc,windsp,psfc,  &
     rhoa,precip,tair,qvair, cdha,cdqa,radsw, rnflx,shflx,lhflx,  &
     gflx,ct,evaprg,evaprtr,evaprr,qvsat,qvsata,f34,              &
     frc_tsfc,frc_tdeep,frc_wsfc,frc_wdp,frc_wcnp,wrmax,relief,   &
     gcan,qcan)

!------------------------------------------------------------------
!
!  AUTHOR: Diandong Ren (dd_ren@rossby.metr.ou.edu)
!          Ming Xue     (mxue@ou.edu)
!  12/08/2001
!
!  2002/02/15 (Yunheng Wang)
!  Changed wrmax to a 2-D array which was a bug found during mpi testing.
!
!  PV-LAND
!  Water vapour leaves the natural surface into the canopy air (air
!  under the panels), whose humidity qcan is a node without storage,
!  like T_can. gcan is the conductance canopy air -> reference and qvair
!  the reference humidity. qcan and the evaporation terms are computed
!  in evapflx (see there); lhflx is the flux into the canopy air, equal
!  to the flux at the reference level since the node has no storage.
!
!------------------------------------------------------------------
!
!  PURPOSE:
!  To calculate the right hand side forcing terms in the
!  soil-vegetation model.
!------------------------------------------------------------------
  USE pvland_nml   ! PV-LAND soil parameters  [SYLVAIN]


  IMPLICIT NONE

  INTEGER :: nx,ny,nz
!  REAL :: eta_fac,wmax,c1max,sigma_sqr

  REAL :: windsp(nx,ny)         ! Wind speed
  REAL :: psfc  (nx,ny)         ! Surface pressure
  REAL :: rhoa  (nx,ny)         ! Air density near the surface
  REAL :: precip(nx,ny)         ! Precipitation rate at the surface

  INTEGER :: soiltyp(nx,ny)     ! Soil type at each point
  INTEGER :: vegtyp (nx,ny)     ! Vegetation type at each point

  REAL :: lai    (nx,ny)        ! Leaf Area Index
  REAL :: veg    (nx,ny)        ! Vegetation fraction

  REAL :: tsfc   (nx,ny)        ! Temperature at ground surface (K)
  REAL :: tdeep  (nx,ny)        ! Deep soil temperature (K)
  REAL :: wetsfc (nx,ny)        ! Surface soil moisture
  REAL :: wetdp  (nx,ny)        ! Deep soil moisture
  REAL :: wetcanp(nx,ny)        ! Canopy water amount
  REAL :: qvsfc  (nx,ny)        ! Effective humidity at surface

  REAL :: snowdpth(nx,ny)       ! Snow depth (m)

  REAL :: tair   (nx,ny)        ! Surface air temperature (K)
  REAL :: qvair  (nx,ny)        ! Reference specific humidity q_ref (kg/kg)
  REAL :: cdha   (nx,ny)        ! Surface drag coeff. for heat
  REAL :: cdqa   (nx,ny)        ! Surface drag coeff. for moisture

  REAL :: radsw  (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx  (nx,ny)        ! Radiation flux at surface
  REAL :: shflx  (nx,ny)        ! Sensible heat flux at surface
  REAL :: lhflx  (nx,ny)        ! Latent heat flux at surface
  REAL :: gflx   (nx,ny)        ! Ground diffusive heat flux
  REAL :: ct     (nx,ny)        ! Soil thermal coefficient

  REAL :: evaprg (nx,ny)        ! Evaporation
  REAL :: evaprtr(nx,ny)        ! Transpiration from leaves
  REAL :: evaprr (nx,ny)        ! Direct evaporation from leaves
  REAL :: f34    (nx,ny)        ! Resistance factor of F3*F4
  REAL :: qvsata (nx,ny)        ! qvsat(tair) (kg/kg)
  REAL :: qvsat  (nx,ny)        !
  REAL :: gcan   (nx,ny)        ! Conductance canopy air -> reference per unit
                                ! ground area [kg m-2 s-1]
  REAL :: qcan   (nx,ny)        ! Specific humidity of the canopy air (kg/kg), output

  REAL :: frc_tsfc(nx,ny)      ! Right hand side forcing for tsfc eq.
  REAL :: frc_tdeep(nx,ny)       ! Right hand side forcing for tsoil eq.
  REAL :: frc_wsfc(nx,ny)       ! Right hand side forcing for wetsfc eq. 
  REAL :: frc_wdp(nx,ny)        ! Right hand side forcing for wetsp eq. 
  REAL :: frc_wcnp(nx,ny)       ! Right hand side forcing for wetcanp eq.
  REAL :: wrmax(nx, ny)         ! Maximum value for canopy moisture, wetcanp

  REAL :: c1wg        ! Coefficient in the surface moisture Eq. of Wg
  REAL :: c2wg        ! Coefficient in the surface moisture Eq. of Wg

  REAL :: pi          ! Pi

  PARAMETER (pi = 3.141592654)

!  REAL :: tau         ! Seconds of a day = 24. * 3600.
!  PARAMETER (tau = 86400.)

  REAL :: dtsfc2      ! Length of half time step in SFCEBM,
                      ! dtsfc2 = dtsfc/2.
  REAL :: log100      ! Constant: alog(100)
                      ! dependent distance from the earth to the sun
  REAL :: cg          ! Soil thermal coefficient for bare ground
  REAL :: wgeq        ! Surface moisture when gravity balances the capillarity
  REAL :: wr2max      ! Tendency to reach the maximum wrmax
  REAL :: runoff      ! Runoff of the interception reservoir.
  REAL :: pnet        ! Residual of precip. and evap.
  REAL :: vegp        ! Precip. intercepted by vegetation
  REAL :: pg          ! Liquid water reaching the soil, (1-veg)P + runoff

!-----------------------------------------------------------------------
!  Gravitational drainage at the bottom of the d2 layer
!  (Mahfouf and Noilhan 1996, J. Appl. Meteor. 35, Eq. 4):
!     dw2/dt = ... - C3/tau * max(0, w2 - wfc)
!  C3 = 0.206 : Table 1 of MN96, sandy clay loam, d2 = 1 m.
!  use_drainage = .FALSE. gives back the original scheme (no drainage).
!-----------------------------------------------------------------------
  REAL :: drain       ! Drainage tendency of w2 (s-1)

  INTEGER :: i, j

  REAL :: tema
  REAL :: eta, c1max, wmax, sig2 


  REAL :: relief       ! Difference between seasonal average skin and deep soil 
                      ! temperature  (t_skin - t_deeplayer).
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! PV-LAND soil parameters (namelist &pvland); a negative value keeps
  ! the ARPS table of soilcst.inc  [SYLVAIN]
  IF (pv_cgv   > 0.) cgv      = pv_cgv
  IF (pv_wsat  > 0.) wsat(:)  = pv_wsat
  IF (pv_wfc   > 0.) wfc(:)   = pv_wfc
  IF (pv_wwlt  > 0.) wwlt(:)  = pv_wwlt
  IF (pv_rsmin > 0.) rsmin(:) = pv_rsmin


  log100    = ALOG(100.0)
  dtsfc2    = dtsfc/2.0

  DO j = 1, ny-1
    DO i = 1, nx-1
      tema = MAX( wetdp(i,j), wwlt(soiltyp(i,j)) )
  
      IF (snowdpth(i,j) >= snowdepth_crit) THEN !snow cover
        ct(i,j) = cg_snow
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)                     &
                          *snowflxfac/(tau*ct(i,j))
      ELSE
        cg = cgsat(soiltyp(i,j))                                      &
            * ( wsat(soiltyp(i,j))/tema )                             &
            **( bslope(soiltyp(i,j))/log100 )
        ! Test
        !cgv = 2E-05    
            
        ct(i,j) = cg * cgv / ( (1.0-veg(i,j)) * cgv + veg(i,j) * cg )  
        !ct(i,j) = 5*ct(i,j)   
        gflx(i,j) = 2.0*pi*(tsfc(i,j)-tdeep(i,j)-relief)/(tau*ct(i,j))
      END IF
  
      shflx(i,j) = rhoa(i,j) * cp * cdha(i,j) * windsp(i,j)           &
                  * ( tsfc(i,j) - tair(i,j) )  
    END DO
  END DO

  IF ( moist == 0 ) THEN

    DO j = 1, ny-1
      DO i = 1, nx-1
         evaprg(i,j) = 0.0
         evaprtr(i,j) = 0.0
         evaprr(i,j) = 0.0
         lhflx(i,j) = 0.0
         qcan(i,j) = qvair(i,j)
      END DO
    END DO

  ELSE

    CALL getqvs(nx,ny,1, 1,nx-1,1,ny-1,1,1, psfc,tsfc,qvsat)

    CALL evapflx_pv(nx,ny,radsw,f34,cdqa,windsp,psfc,rhoa,qvair,         &
                   soiltyp,vegtyp,lai,veg,tsfc,wetsfc,wetdp,wetcanp,    &
                   snowdpth,evaprg,evaprtr,evaprr,lhflx,qvsat,gcan,qcan)
  END IF

  DO j = 1, ny-1
    DO i = 1, nx-1
      evaprg (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprg(i,j)
      evaprtr(i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprtr(i,j)
      evaprr (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*evaprr(i,j)
      lhflx  (i,j) = rhoa(i,j)*cdqa(i,j)*windsp(i,j)*lhflx(i,j)
      lhflx(i,j) = lhflx(i,j) * lathv       ! Latent heat flux
    END DO
  END DO

  DO j = 1, ny-1
    DO i = 1, nx-1
      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
        frc_tsfc (i,j) = 0.
        frc_tdeep(i,j) = 0.
      ELSE
        frc_tsfc(i,j) = ct(i,j)*(rnflx(i,j)-shflx(i,j)-lhflx(i,j)-gflx(i,j))
        IF ( snowdpth(i,j) >= snowdepth_crit ) THEN
          frc_tdeep(i,j)= (tsfc(i,j)-tdeep(i,j)-relief)*tauinv*snowflxfac ! Snow cover
        ELSE
          frc_tdeep(i,j)= (tsfc(i,j) - tdeep(i,j)-relief) * tauinv      
        END IF
      END IF
    END DO
  END DO

  IF ( moist /= 0 ) THEN   

    DO j = 1, ny-1                                       !HYDROLOGY
      DO i = 1, nx-1
        wrmax(i,j) = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! meter
        IF ( soiltyp(i,j) == 12 .OR.  soiltyp(i,j) == 13 .OR.           &
             snowdpth(i,j) >= snowdepth_crit ) THEN
          frc_wsfc(i,j) = 0.
          frc_wdp(i,j)  = 0.
          frc_wcnp(i,j) = 0.
        ELSE

          tema = MAX( wetsfc(i,j), wwlt(soiltyp(i,j)) )

          c1wg = 0.4* c1sat(soiltyp(i,j))                               &
                 * ( wsat(soiltyp(i,j)) / tema )                        &
                 **( bslope(soiltyp(i,j)) / 2.0 + 1.0)

!--------------------------------------------------------------------------
! Replacement Cl to improve dry soils (NM1996 - A.3)   (JAB)
!--------------------------------------------------------------------------

          IF (wetsfc(i,j) < wwlt(soiltyp(i,j)) ) THEN

            eta = (-0.01815 * tsfc(i,j) + 6.41 ) * wwlt(soiltyp(i,j)) + &
                  ( 0.0065 * tsfc(i,j) - 1.4 )
      
            c1max = (1.19 * wwlt(soiltyp(i,j)) - 5.09)*0.01*tsfc(i,j)   &
                   +( 1.464 * wwlt(soiltyp(i,j)) + 17.86 )
      
            wmax = eta * wwlt(soiltyp(i,j))
      
            sig2 = - ( (2 * alog ( 0.01 / c1max )) / ( wmax*wmax) )
      
            c1wg = c1max * EXP ( -0.5* ( (wetsfc(i,j)-wmax)*(wetsfc(i,j)-wmax)*sig2 ) )
          ENDIF

!----------------------------------------------------------------------------

          c2wg = c2ref(soiltyp(i,j)) * wetdp(i,j)                       &
                / ( wsat(soiltyp(i,j)) - wetdp(i,j) + wetsml )
          wgeq = wetdp(i,j) - wsat(soiltyp(i,j))                        &
                  * awgeq(soiltyp(i,j))                                 &
                  * ( wetdp(i,j) / wsat(soiltyp(i,j)) )                 &
                  ** pwgeq(soiltyp(i,j))                                &
                  * ( 1.0 - ( wetdp(i,j) / wsat(soiltyp(i,j)) )         &
                  ** ( 8 * pwgeq(soiltyp(i,j)) ) )

          wr2max = ( wrmax(i,j) - wetcanp(i,j) ) / dtsfc2
          vegp   = veg(i,j) * precip(i,j)
          pnet   = vegp - evaprr(i,j)
          tema   = pnet - wr2max * rhow
          runoff = MAX( tema, 0.0 )
          vegp   = vegp - runoff
          frc_wcnp(i,j) = ( vegp - evaprr(i,j) ) / rhow

          pg = ( 1.0 - veg(i,j) ) * precip(i,j) + runoff

          frc_wsfc(i,j) = c1wg * ( pg - evaprg(i,j) ) / (rhow * d1) - c2wg * ( wetsfc(i,j) - wgeq ) * tauinv
          frc_wdp(i,j)  = ( pg - evaprg(i,j) - evaprtr(i,j) ) / ( rhow * d2 )

          IF ( pv_drain_opt == 1 ) THEN
            IF (wetdp(i,j) < (wsat(soiltyp(i,j))-0.04)) THEN
            	drain = pv_c3_drain * MAX( wetdp(i,j) - wfc(soiltyp(i,j)), 0.0 ) * tauinv
            ELSE
                drain = 0.0
            END IF
            frc_wdp(i,j) = frc_wdp(i,j) - drain
          END IF
        END IF

      END DO
    END DO           

  END IF

  RETURN
END SUBROUTINE soilebm_frc_pv

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE evapflx_pv                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE evapflx_pv(nx,ny, radsw, f34, cdqa,                             &
           windsp,psfc,rhoa,qvair,                                      &
           soiltyp,vegtyp,lai,veg,                                      &
           tsfc,wetsfc,wetdp,wetcanp,snowdpth,                          &
           evaprg,evaprtr,evaprr,qvflx,qvsat,gcan,qcan)
!
!-----------------------------------------------------------------------
!
!  PV-LAND: humidity of the canopy air
!
!  The evaporation terms below are written against the humidity of the
!  canopy air qcan (air under the panels), not against the reference
!  humidity qvair. The canopy air is a node without storage, like T_can:
!  the vapour leaving the surface,
!     E_nat = g_nat (qvsats - qcan),
!     g_nat = rhoa [ (1-veg) rhgs/Ra + veg (1-delta)/(Ra+Rs) + veg delta/Ra ]
!           = rhoa cdq windsp [ (1-veg) rhgs + veg ((1-delta) Ra/(Ra+Rs) + delta) ],
!  equals the vapour leaving the canopy air towards the reference level,
!     E_nat = gcan (qcan - qvair),
!  hence
!     qcan = (g_nat qvsats + gcan qvair) / (g_nat + gcan).
!  gcan [kg m-2 s-1] is given by the caller (canopy air -> reference,
!  per unit ground area). With gcan <= 0 there is no air node and
!  qcan = qvair (original behaviour).
!
!  Dew: the branch (evaporation, or dew with no stomatal resistance,
!  delta = 1) is chosen with the sign of qvsats - qvair, as before. Since
!  qcan is a weighted mean of qvsats and qvair, qvsats - qcan has the same
!  sign, so g_nat is computed in that branch and no iteration is needed.
!
!  The terms below (Eg', Etr', Er', qvflx') are therefore written with
!  qcan in place of qvair.
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate:
!
!  1. Evaporation from ground surface,
!
!        evaprg = rhoa * cdq * windsp * evaprg'
!
!  where
!
!        evaprg' = (1.-veg) * (rhgs * qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  2. Direct evaporation from the fraction delta of the foliage covered
!     by intercepted water.
!
!        Er = rhoa * cdq * windsp * Er'
!
!        Er' = delta * veg * (qvsats-qvair)
!
!  3. Transpiration of the remaining part (1-delta) of leaves,
!
!        Etr = rhoa * cdq * windsp * Etr'
!
!  where
!
!        Etr' = veg * (1-delta) * Ra/(Ra+Rs) * (qvsats-qvair )
!
!  and Ra is aerodynamic resistance and Rs is the surface resistance
!
!                     1
!        Ra = ----------------
!               cdq * windsp
!
!                   Rsmin
!        Rs = ------------------
!              LAI*F1*F2*F3*F4
!
!
!               f + Rsmin/Rsmax
!        F1 = -------------------
!                   f + 1
!
!                   Rg    2
!        f = 0.55 ----- -----
!                  Rgl   LAI
!
!             -  1,                             Wfc < W2
!             |
!             |    W2 - Wwlt
!        F2 = -  ------------,              Wwlt <= W2 <= Wfc
!             |   Wfc - Wwlt
!             |
!             -  0,                             W2 < Wwlt
!
!
!               1-0.06*(qvsats-qvair),   qvsats-qvair <= 12.5 g/kg
!        F3 = {
!               0.25,                      otherwise
!
!
!        F4 = 1 - 0.0016 * (298-tair)**2
!
!  4. Water vapor flux, qvflx,
!
!        qvflx = rhoa * cdq * windsp * qvflx'
!
!        qvflx' = (Eg' + Etr' + Er') = (qvsfc - qvair)
!
!     where qvsfc is the effective surface specific humidity
!
!        (qvsfc - qvair) = (Eg' + Etr' + Er')
!
!  This subroutine will solve Eg', Etr', Er', and qvflx'
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Vince Wong
!  4/20/94
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    radsw    Incoming solar radiation
!
!    f34      f3*f4;
!             f3 = Fractional conductance of atmospheric vapor pressure
!             f4 = Fractional conductance of air temperature
!
!    cdqa     Array for surface drag coefficient for water vapor
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    veg      Vegetation fraction
!
!    tsfc     Surface soil temperature (K)
!    wetsfc   Surface soil moisture
!    wetdp    Deep soil moisture
!    wetcanp  Vegetation moisture
!
!    psfc     Surface pressure
!    qvair    Specific humidity at the reference level
!    windsp     Wind speed near the surface
!    rhoa     Air density near the surface
!    gcan     Conductance canopy air -> reference [kg m-2 s-1]
!             (<= 0: no air node, qcan = qvair)
!
!  OUTPUT:
!
!    evaprp   Evaporation from groud surface
!    evaprtr  Transpiration of the remaining part (1-delta) of leaves
!    evaprr   Direct evaporation from the fraction delta
!    qvflx    Water vapor flux
!    qcan     Specific humidity of the canopy air (kg/kg)
!
!  WORK ARRAY:
!
!
!-----------------------------------------------------------------------
!
  USE pvland_nml   ! PV-LAND soil parameters  [SYLVAIN]

  IMPLICIT NONE
!
  INTEGER :: nx,ny
!
  REAL :: radsw  (nx,ny)

  REAL :: f34    (nx,ny)

  REAL :: cdqa   (nx,ny)

  INTEGER :: soiltyp(nx,ny)
  INTEGER :: vegtyp (nx,ny)
  REAL :: lai    (nx,ny)
  REAL :: veg    (nx,ny)

  REAL :: tsfc   (nx,ny)
  REAL :: wetsfc (nx,ny)
  REAL :: wetdp  (nx,ny)
  REAL :: wetcanp(nx,ny)
  REAL :: snowdpth(nx,ny)
!
  REAL :: psfc   (nx,ny)
  REAL :: qvair  (nx,ny)
  REAL :: windsp   (nx,ny)
  REAL :: rhoa   (nx,ny)

  REAL :: evaprg (nx,ny)
  REAL :: evaprtr(nx,ny)
  REAL :: evaprr (nx,ny)
  REAL :: qvflx  (nx,ny)
  REAL :: qvsat  (nx,ny)
  REAL :: gcan   (nx,ny)       ! Conductance canopy air -> reference [kg m-2 s-1]
  REAL :: qcan   (nx,ny)       ! Specific humidity of the canopy air (kg/kg)

!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc and phycst.inc
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
!
!-----------------------------------------------------------------------
!
!  Local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: pi
  PARAMETER ( pi = 3.141592654 )
!
  INTEGER :: i, j

  REAL :: wrmax       ! Maximum value for canopy moisture, wetcanp
  REAL :: rstcoef     ! Coefficient of resistance
  REAL :: delta
  REAL :: rhgs        ! R.H. at ground surface
!
  REAL :: tema
  REAL :: temb

  REAL :: pterm
  REAL :: mterm
  REAL :: waf         ! Available soil moisture fraction
  REAL :: bw          ! Half point of the available soil mstr frctn curve
  REAL :: ps          ! Shelter factor
  REAL :: gnat        ! Conductance of the natural surface [kg m-2 s-1]

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! PV-LAND soil parameters (namelist &pvland); a negative value keeps
  ! the ARPS table of soilcst.inc  [SYLVAIN]
  IF (pv_cgv   > 0.) cgv      = pv_cgv
  IF (pv_wsat  > 0.) wsat(:)  = pv_wsat
  IF (pv_wfc   > 0.) wfc(:)   = pv_wfc
  IF (pv_wwlt  > 0.) wwlt(:)  = pv_wwlt
  IF (pv_rsmin > 0.) rsmin(:) = pv_rsmin

  DO j = 1, ny-1
    DO i = 1, nx-1

      IF ( soiltyp(i,j) == 12 .OR. soiltyp(i,j) == 13 ) THEN
!
!-----------------------------------------------------------------------
!
!  Over water and ice, qvsfc should be saturated
!
!-----------------------------------------------------------------------
!
        gnat = rhoa(i,j) * cdqa(i,j) * windsp(i,j)
        IF ( gcan(i,j) > 0.0 ) THEN
          qcan(i,j) = ( gnat*qvsat(i,j) + gcan(i,j)*qvair(i,j) ) / ( gnat + gcan(i,j) )
        ELSE
          qcan(i,j) = qvair(i,j)
        END IF

        qvflx(i,j) = qvsat(i,j) - qcan(i,j)

        evaprg (i,j) = qvflx(i,j)
        evaprtr(i,j) = 0.0
        evaprr (i,j) = 0.0

      ELSE                                  ! over land

        wrmax = 0.2 * 1.e-3 * veg(i,j) * lai(i,j)     ! in meter
!
!-----------------------------------------------------------------------
!
!  In order to calculate the qv flux at the surface, we need to
!  calculate some parameters like the resistance coefficient,
!
!      rstcoef = rsta / (rsts + rsta)
!
!  where rsta is the aerodynamic resistance and rsts is the surface
!  resistances.
!
!      rsta = 1 / ( cdq * Va )
!
!      rsts = rsmin(vtyp) / ( lai(vtyp) * f1 * f2 * f3 * f4 )
!
!  f3 * f4 is time-independent and has been calculated previously
!  and stored in f34(i,j)
!
!-----------------------------------------------------------------------
!
!  Calculate f1.
!
!        f = 0.55*(radsw/rgl(vtyp))*(2./lai(vtyp))       ! NP, Eq. 34
!        f1 = ( rsmin(vtyp)/rsmax + f ) / ( 1. + f )     ! NP, Eq. 34
!
!  Note: the incoming solar radiation radsw is stored in radsw(i,j).
!
!-----------------------------------------------------------------------
!
        IF ( lai(i,j) == 0. ) THEN        ! No vegetation, I/W
          rstcoef = 1.
        ELSE
!          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )               &
!               * ( 2.0 / lai(i,j) )

          temb = 0.55 * ( radsw(i,j) / rgl(vegtyp(i,j)) )               &
               * ( 2.0 )       !(XP2001 - Eq.8)   (JAB) 

          rstcoef = ( rsmin(vegtyp(i,j)) / rsmax + temb )               &
                  / ( 1.0 + temb )
        END IF
!
!-----------------------------------------------------------------------
!
!  Calculate f2 and f1*f2.
!
!-----------------------------------------------------------------------
!
!        pterm = .5 + SIGN( .5, wetdp(i,j) - wfc(soiltyp(i,j)) )
!        mterm = SIGN( .5, wetdp(i,j) - wwlt(soiltyp(i,j)) )              &
!              - SIGN( .5, wetdp(i,j) - wfc(soiltyp(i,j)) )
!
!        rstcoef = rstcoef * ( pterm + mterm                             &
!                * ( wetdp(i,j)-wwlt(soiltyp(i,j)) )                     &
!                / ( wfc(soiltyp(i,j)) - wwlt(soiltyp(i,j)) ) )
        waf = ( wetdp(i,j) - wwlt(soiltyp(i,j)) )                 &
              / ( pv_wref - wwlt(soiltyp(i,j)) )

        ! F2 : soil moisture stress function (XP2001, Eq. 9)
        ! Logistic form: F2 = 1/{1 + exp[-5*(waf - bw)]}
        ! bw must be normalized (= 1/3) so that F2 = 0.5 at
        ! w2 = (wfc-wwlt)/3 + wwlt for ALL soil types.
        ! (The original code used bw in absolute units, which
        !  only works for silty clay where bw_abs ~ 1/3.)
        rstcoef = rstcoef / ( 1.0 + EXP( -5.0 * (waf - 1.0/3.0) ) )
!             (XP2001 - Eq.9, bw normalized)                        (JAB)

        
        
        ! F2: linear (Noilhan and Planton 1989; Jacquemin and Noilhan 1990).
        ! Thresholds: pv_wref = Theta_ref of sand from Chen and Dudhia (2001),
        ! scaled to wsat = 0.31; wwlt = measured late-summer floor.
        !waf = ( wetdp(i,j) - wwlt(soiltyp(i,j)) )                 &
        !      / ( pv_wref - wwlt(soiltyp(i,j)) )
        !waf = MAX( MIN( waf, 1.0 ), 0.0 )
        !rstcoef = rstcoef * waf
!
!-----------------------------------------------------------------------
!
!  Calculate lai*f1*f2*f3*f4 where f3*f4 is stored in f34(i,j).
!
!-----------------------------------------------------------------------
!
!        rstcoef = lai(i,j)*rstcoef*f34(i,j)      ! lai*f1*f2*f3*f4

        ps = 0.3 * lai(i,j) + 0.7                 ! XP2001             (JAB)
        rstcoef = lai(i,j)*rstcoef*f34(i,j)/ps    ! lai*f1*f2*f3*f4/ps (JAB)
!
!-----------------------------------------------------------------------
!
!  Calculate the resistance coefficient, rsta/(rsts+rsta)
!
!        rsts = rsmin(vtyp)/(lai(i,j)*f1*f2*f3*f4) ! Sfc. resistance
!        rsta = 1./(cdh*va)                 ! NP, between Eq. 32 & 33
!        rstcoef = rsta/(rsta+rsts)
!                = 1/(1+rsts/rsta)
!
!-----------------------------------------------------------------------
!
        tema = rsmin(vegtyp(i,j)) * cdqa(i,j) * windsp(i,j)

        IF ( ABS(rstcoef) > 1.0E-30 ) THEN
          rstcoef = 1.0 / (1.0 + tema/rstcoef)
        END IF
!
!-----------------------------------------------------------------------
!
!  1. evaprg'
!
!        evaprg' = (1.-veg) * (rhgs*qvsats - qvair)
                                ! Evaporation from the ground
                                ! NP, Eq. 27
!
!  evaprg will be stored for current and future use to
!  calculate the latent heat flux and soil moisture transports.
!
!-----------------------------------------------------------------------
!
        pterm = .5 + SIGN( .5, wetsfc(i,j)-1.1*pv_wref )

        rhgs = pterm                                                    &
             + (1.-pterm) * ( 0.25 * ( 1.0 - COS( wetsfc(i,j)           &
                              * pi / (1.1*pv_wref))) ** 2 )

!      IF (snowdpth(i,j) .ge. snowdepth_crit) rhgs=1.0 !Snow cover

!-----------------------------------------------------------------------
!
!  Fraction delta of the foliage covered by intercepted water, and
!  branch evaporation / dew. The branch is chosen with the reference
!  humidity qvair (see header): for dew, delta = 1 (no stomatal
!  resistance).
!
!-----------------------------------------------------------------------
!
        IF ( wrmax == 0.0 ) THEN
          delta = 0.0
        ELSE
          delta = ( wetcanp(i,j) / wrmax ) ** 0.66666667
        END IF

        pterm = .5 + SIGN( .5, qvsat(i,j) - qvair(i,j) )
        delta = pterm * delta + ( 1. - pterm )
!
!-----------------------------------------------------------------------
!
!  Conductance of the natural surface g_nat [kg m-2 s-1] in this branch,
!  and humidity of the canopy air (node without storage):
!     qcan = (g_nat qvsats + gcan qvair) / (g_nat + gcan)
!
!-----------------------------------------------------------------------
!
        gnat = rhoa(i,j) * cdqa(i,j) * windsp(i,j)                      &
             * ( ( 1.0 - veg(i,j) ) * rhgs                              &
               + veg(i,j) * ( ( 1.0 - delta ) * rstcoef + delta ) )

        IF ( gcan(i,j) > 0.0 ) THEN
          qcan(i,j) = ( gnat*qvsat(i,j) + gcan(i,j)*qvair(i,j) )        &
                    / ( gnat + gcan(i,j) )
        ELSE
          qcan(i,j) = qvair(i,j)
        END IF

!        evaprg(i,j) = ( 1.0 - veg(i,j) )                                &
!                    * ( rhgs * qvsat(i,j) - qvair(i,j) )

        evaprg(i,j) = ( 1.0 - veg(i,j) )                                &
                    *  rhgs * ( qvsat(i,j) - qcan(i,j) )  !XP2001 (JAB)

!
!-----------------------------------------------------------------------
!
!  2. Transpiration of the remaining part (1-delta) of leaves, Etr',
!
!        Etr' = (1-delta) * veg
!             * Ra/(Ra+Rs) * ( qvsats - qvair )
!
!  3. Direct evaporation from the fraction delta, Er'
!
!        Er' = delta * veg * ( qvsats - qvair )
!
!  Er' and Etr' are stored for future use to calculate the latent heat
!  flux and soil moisture transports.
!
!  (delta and the branch have been set above, before g_nat.)
!
!-----------------------------------------------------------------------
!
        tema = veg(i,j) * ( qvsat(i,j) - qcan(i,j) )

        evaprtr(i,j) = ( 1.0 - delta ) * rstcoef * tema
        evaprr (i,j) = delta * tema
!
!-----------------------------------------------------------------------
!
!  4. Water vapor flux, qvflx',
!
!        qvflx' = evaprg' + evaprtr' + evaprr'
!
!  qvflx will be saved for the future use to calculate the latent
!  heat flux.
!
!-----------------------------------------------------------------------
!
        qvflx(i,j) = evaprg(i,j) + evaprtr(i,j) + evaprr(i,j)
      END IF
                                                   ! NP, expl. Eq. 26-27
    END DO
  END DO

  RETURN
END SUBROUTINE evapflx_pv  

