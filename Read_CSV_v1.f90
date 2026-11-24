PROGRAM read_csv

  IMPLICIT NONE

  INTEGER :: ntot

  ! ====== File I/O ======
  INTEGER :: temp, ios, t
  CHARACTER(LEN=500) :: filename
  
  ! ====== CSV input arrays ======
  REAL, ALLOCATABLE, DIMENSION (:) :: uref_input, udir_input, &
                            tref_input, qvref_input, rain_input, &
                            swout_input, swin_input, &
                            lwout_input, lwin_input, &
                            pref_input, rhoref_input, &
                            H_meas, LE_meas, G_meas, Rn_meas, &
                            S_diffuse_input, lambda_input, phi_input, &
                            Lobk_input

  !=====================================================================
  ! READ CSV INPUT
  !=====================================================================
  
  filename='./Data_mesures/recorded_data.csv'
  
  WRITE(*,*) 'Reading forcing data from CSV...'
  
  OPEN(UNIT=10, FILE=filename, STATUS='old')
  READ(10,*)  ! Skip header
  
  temp = 1
  DO
    READ(10,'(A)', IOSTAT=ios)
    IF (ios /= 0) EXIT
    temp = temp + 1
  END DO
  ntot = temp - 1
  CLOSE(10)
  
  print *,'ntot=',ntot

  ALLOCATE (uref_input(ntot), udir_input(ntot), &
                            tref_input(ntot), qvref_input(ntot), rain_input(ntot), &
                            swout_input(ntot), swin_input(ntot), &
                            lwout_input(ntot), lwin_input(ntot), &
                            pref_input(ntot), rhoref_input(ntot), &
                            H_meas(ntot), LE_meas(ntot), G_meas(ntot), Rn_meas(ntot), &
                            S_diffuse_input(ntot), lambda_input(ntot), phi_input(ntot), &
                            Lobk_input(ntot))
                            
  OPEN(UNIT=10, FILE=filename, STATUS='old')
  READ(10,*)  ! Skip header
  
  DO temp=1,ntot
    READ(10,*) t, uref_input(temp), udir_input(temp), &
                            tref_input(temp), qvref_input(temp), rain_input(temp), &
                            swout_input(temp), swin_input(temp), &
                            lwout_input(temp), lwin_input(temp), &
                            pref_input(temp), rhoref_input(temp), &
                            H_meas(temp), LE_meas(temp), G_meas(temp), Rn_meas(temp), &
                            S_diffuse_input(temp), lambda_input(temp), phi_input(temp), &
                            Lobk_input(temp)
    
	
    tref_input(temp) = tref_input(temp) + 273.15
  END DO
  CLOSE(10)

  !=====================================================================
  ! WRITE TECPLOT FILE
  !=====================================================================

  OPEN(UNIT=20, FILE='pv_land_input.plt', STATUS='replace')
  WRITE(20,*)' TITLE="TITLE"'
  WRITE(20,'(A)') 'VARIABLES="timestep","uref_input", "udir_input", &
                            "tref_input", "qvref_input", "rain_input", &
                            "swout_input", "swin_input", &
                            "lwout_input", "lwin_input", &
                            "pref_input", "rhoref_input", &
                            "H_meas", "LE_meas", "G_meas", "Rn_meas", &
                            "S_diffuse_input", "lambda_input", "phi_input", &
                            "Lobk_input"'
  WRITE (20,*)'ZONE T="CHAMPS" I=', ntot

  DO temp = 1, ntot
      WRITE(20, *) temp, &
		uref_input(temp), udir_input(temp), &
                            tref_input(temp), qvref_input(temp), rain_input(temp), &
                            swout_input(temp), swin_input(temp), &
                            lwout_input(temp), lwin_input(temp), &
                            pref_input(temp), rhoref_input(temp), &
                            H_meas(temp), LE_meas(temp), G_meas(temp), Rn_meas(temp), &
                            S_diffuse_input(temp), lambda_input(temp), phi_input(temp), &
                            Lobk_input(temp)
  ENDDO
  CLOSE(20)
  !=====================================================================
  ! FINALIZE
  !=====================================================================

  DEALLOCATE (uref_input, udir_input, &
                            tref_input, qvref_input, rain_input, &
                            swout_input, swin_input, &
                            lwout_input, lwin_input, &
                            pref_input, rhoref_input, &
                            H_meas, LE_meas, G_meas, Rn_meas, &
                            S_diffuse_input, lambda_input, phi_input, &
                            Lobk_input)
  
  WRITE(*,*) ''
  WRITE(*,*) '========================================='
  WRITE(*,*) 'Simulation complete!'
  WRITE(*,*) '========================================='

END PROGRAM read_csv
