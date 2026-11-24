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
!  Offline version: the single-column domain is entirely covered by
!  the PV park. The ARPS version (src/arps/pvland_park_map.f90) builds
!  the park map on the ARPS grid with MPI.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny
  REAL :: pv_park(nx,ny)
  REAL :: tem1(nx,ny)          ! unused offline (same interface as ARPS)

  pv_park(:,:) = 1.

  RETURN
END SUBROUTINE pvland_park_map
