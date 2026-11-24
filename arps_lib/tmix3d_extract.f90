!==== Extracted from ARPS tmix3d.f90 : SUBROUTINE tridiag2 (do not edit) ====
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TRIDIAG2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tridiag2(nx,ny,nz,ibgn,iend,jbgn,jend,                       &
           kbgn,kend,a,b,c,d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the solution for a tridiagonal equations.
!  Note:
!       a: left of main diagonal, useful range [kbgn+1,kend];
!       b: main diagonal, useful range [kbgn,kend];
!       c: right of main diagonal, useful range [kbgn,kend-1];
!       d: right hand side of equations.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: X. Song, Donghai Wang and M. Xue
!  4/2/96
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where operation begins.
!    iend     i-index where operation ends.
!    jbgn     j-index where operation begins.
!    jend     j-index where operation ends.
!    kbgn     k-index where operation begins.
!    kend     k-index where operation ends.
!
!    a        left of main diagonal, useful range [kbgn+1,kend]
!    b        main diagonal, useful range [kbgn,kend]
!    c        right of main diagonal, useful range [kbgn,kend-1]
!    d        right hand side of equations
!
!  OUTPUT:
!
!    d        The solution
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE             ! Force explicit declarations

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: c(nx,ny,nz)       ! Right of main diagonal
  REAL :: a(nx,ny,nz)       ! Left of main diagonal
  REAL :: d(nx,ny,nz)       ! Right hand side of equations
  REAL :: b(nx,ny,nz)       ! Main diagonal

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
  INTEGER :: i,j,k
  REAL :: r
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn+1,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        r=a(i,j,k)/b(i,j,k-1)
        b(i,j,k)=b(i,j,k)-r*c(i,j,k-1)
        d(i,j,k)=d(i,j,k)-r*d(i,j,k-1)
      END DO
    END DO
  END DO

  DO j=jbgn,jend
    DO i=ibgn,iend
      d(i,j,kend)=d(i,j,kend)/b(i,j,kend)
    END DO
  END DO

  DO k=kend-1,kbgn,-1
    DO j=jbgn,jend
      DO i=ibgn,iend
        d(i,j,k)=(d(i,j,k)-c(i,j,k)*d(i,j,k+1))/b(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE tridiag2

