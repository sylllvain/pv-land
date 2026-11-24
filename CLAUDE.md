# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

PV-LAND: a Fortran 90 land-surface model for a ground-mounted photovoltaic park. It runs offline (not coupled to ARPS) and is forced by measured meteorological data. It computes the energy balance of the natural ground surface, the soil, and the PV panels (front glass / cell / back), and compares modelled fluxes (H, LE, Rn, G) with measurements. A large part of the soil/surface physics is reused from the ARPS model (`soilebm3d.f90`, `dragcoef_sfcflx.f90`, `thermolib3d.f90`, `tridiag2.f90`, and the `*.inc` common-block files). The user writes in French.

## Build and run

```bash
make            # builds pv_land.exe with ifx (Intel oneAPI)
make clean      # remove *.o, *.mod, executable
make rebuild
make run        # build then ./pv_land.exe
```

- The compiler is `ifx` with `-g -traceback -check bounds`. gfortran flags are kept commented out in the Makefile; switching to gfortran also requires `LDFLAGS = -llapack -lblas` plus the system LAPACK packages.
- LAPACK (`SGESV`, used in `radiosity_SW.f90` / `radiosity_LW.f90`) comes from Intel MKL (`-qmkl=sequential`). The oneAPI environment (`setvars.sh`) is usually not sourced locally, so the Makefile hardcodes `MKLROOT` and `INTEL_LIB` under `/opt/intel/oneapi/.../latest` and sets an rpath to them. Override these variables if the install lives somewhere else.
- Header dependencies are not tracked: after editing a `.inc` file, run `make rebuild`.
- On the HPC cluster the run is submitted with `sbatch batch_curta.sh` (Slurm; it loads `compiler/2023.0.0` and runs `./pv_land.exe`).
- There is no test suite. Validate a change by running the model and comparing `H_model`/`LE_model`/... with the `*_meas` columns in `pv_land_output.csv`. The radiosity modules also have `validate_*_conservation` and `check_*_physical_validity` routines.

## Inputs / outputs (hardcoded file names in `pv_land_main.f90`)

- Input: `recorded_data_adjusted_daily_summer_2025.csv`. Columns are read by position, in this order: `t, wind_speed, wind_dir, T6 (°C, converted to K), specific_humidity, Rain, SWout, SWin, LWout, LWin, air_pressure, air_density, H, LE, G, RN, SW_diffus, lambda, phi, L (Obukhov length)`. Rows are 30 min apart. `lambda`/`phi` are solar angles in the park frame, precomputed by the standalone `calcul_angle_zen_azi_sylvain_v2.f90` program.
- Outputs: `pv_land_output.csv` (modelled vs measured fluxes, natural-surface diagnostics) and `pv_diag_energy_balance.csv` (PV panel energy-balance terms).
- Park geometry, optical properties and roughness are `PARAMETER`s at the top of `pv_land_main.f90`. Changing them means recompiling.

## Architecture

Call flow per 30-min forcing interval:

1. `pv_land_main` reads the CSV and loops over the forcing rows. For each interval it fills 2-time-level arrays (`*i(nx,ny,2)`), then calls `pv_land`. The grid is a dummy `nx=ny=2`; the physics is effectively single-column.
2. `pv_land` (`pv_land.f90`) sub-steps with `nsfcst = 1800/dtsfc` small steps (`dtsfc` is set in main and passed through ARPS common blocks in `globcst.inc`). At each sub-step it:
   - linearly interpolates the forcing in time;
   - gets wind-direction-dependent park roughness/displacement (`z0_theta`, `d_theta`) and wind speeds at the canopy (`can`), panel-air (`PVa`) and reference heights from the analytic PV-park wind profile in `profil.f90` (`compute_U_plus`, with MOST stability correction from the Obukhov length);
   - computes surface/air exchange coefficients with ARPS `dragcoef` (`dragcoef_sfcflx.f90`) for nat, PV front/back, can→ref and PVa→ref;
   - solves the two air-node temperatures `T_can`/`T_PVa` diagnostically (`pv_convection.f90`);
   - integrates all prognostic nodes (`T_nat`, `T_soil`, soil moisture, `T_PVf`, `T_cell`, `T_PVb`) together with an RK2 predictor/corrector. Both stages call `compute_all_forcings`.
3. `compute_all_forcings` builds the tendencies:
   - shortwave: geometric view factors (`view_factors.f90`), direct-beam fractions on each surface (`solar_radiation.f90`), multiple reflections via an N-surface radiosity linear system (`radiosity_SW.f90`);
   - longwave: same radiosity approach (`radiosity_LW.f90`);
   - ground/soil: ARPS `soilebm_frc` (`soilebm3d.f90`);
   - panels: conduction through the glass/EVA/cell/backsheet stack and electrical power removal (`pv_thermal.f90`).

The radiating surfaces are the ground (`nat`), the PV front (`PVf`) and the PV back (`PVb`), plus the sky. Naming uses the suffixes `_nat`, `_PVf`, `_PVb`, `_cell`, `_can`, `_PVa`, `_ref` consistently across files.

Conventions inherited from ARPS: fixed-form-style upper-case Fortran, implicit interfaces (no modules, so argument lists must match exactly between caller and callee), and global constants/state shared through `INCLUDE`d common blocks (`globcst.inc`, `phycst.inc`, `sfcphycst.inc`, `soilcst.inc`, `grid.inc`, `timelvls.inc`).

## Files not in the build

The Makefile `SOURCES` list is the source of truth. These files are not built: `compute_all_forcings_old.f90`, `pv_nodes_forcing.f90`, `Read_CSV_v1.f90`, and the `calcul_angle_*` standalone programs. `files_new/` holds alternative versions of some sources; `pv_land.f90`, `pv_land_main.f90` and `dragcoef_sfcflx.f90` there differ from the ones at the top level. Edit the top-level files unless asked otherwise.
