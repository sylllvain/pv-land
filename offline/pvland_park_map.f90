!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE PVLAND_PARK_MAP              ######
!######                (offline PV-LAND model)               ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pvland_park_map(nx,ny,pv_park,tem1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Location of the PV park, called by INIT_PVLAND (core/pvland.f90).
!  Offline version, same namelist options as ARPS (&pvland, module
!  pvland_nml) applied to the offline grid (computed column (1,1)):
!    pv_park_opt = 0 : the column is in the PV park
!    pv_park_opt = 1 : in the park if (1,1) is inside the rectangle
!                      pv_park_i1:pv_park_i2 x pv_park_j1:pv_park_j2
!    pv_park_opt = 2 : no PV park (natural ground only, as ARPS sfcphy = 4)
!  The ARPS version (src/arps/pvland_park_map.f90) builds the park map
!  on the ARPS grid with MPI.
!
!-----------------------------------------------------------------------
!
  USE pvland_nml, ONLY : pv_park_opt, pv_park_i1, pv_park_i2,         &
                         pv_park_j1, pv_park_j2

  IMPLICIT NONE

  INTEGER :: nx,ny
  REAL :: pv_park(nx,ny)
  REAL :: tem1(nx,ny)          ! unused offline (same interface as ARPS)

  INTEGER :: i, j

  DO j = 1, ny
    DO i = 1, nx
      IF (pv_park_opt == 0) THEN
        pv_park(i,j) = 1.
      ELSE IF (pv_park_opt == 1 .AND. i >= pv_park_i1 .AND. i <= pv_park_i2 &
                                .AND. j >= pv_park_j1 .AND. j <= pv_park_j2) THEN
        pv_park(i,j) = 1.
      ELSE
        pv_park(i,j) = 0.
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE pvland_park_map
