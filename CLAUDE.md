# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

PV-LAND offline: a Fortran 90 land-surface model for a ground-mounted photovoltaic park, forced by measured meteorological data. It computes the energy balance of the natural ground surface, the soil and the PV panels (front glass / cell / back), and compares modelled fluxes (H, LE, Rn, G) with measurements. The user writes in French.

**ARPS is the reference.** The same PV-LAND code runs coupled in ARPS (`../arps5.3.4_PVLAND`, `sfcphy = 5`; GitHub `sylllvain/arps_pvland`). This repository is only an offline driver around exact copies of the ARPS sources, so that PV-LAND changes can be moved between both models without manual porting. See `README.md` (French) for the user-facing summary.

## Layout and the "do not edit" rule

- `core/`: PV-LAND physics, verbatim copies of ARPS `src/arps/module_pvland.f90`, `pvland.f90` (`INIT_PVLAND`, `DRAG_PVLAND`, `solve_air_node_temperatures`, `compute_all_forcings`, `compute_lambda_phi`), `pvland_profil.f90`, `pvland_view_factors.f90`, `pvland_solar_radiation.f90`, `pvland_radiosity_SW/LW.f90`, and `pvland_driver.f90` (module `pvland_sfcphy_work` + `SUBROUTINE pvland_sfcphy`: one PV-LAND time step, called by ARPS `sfcphy3d.f90` for `sfcphy == 5` and by the offline main program).
- `arps_lib/`: ARPS routines used by PV-LAND: `soilebm3d.f90`, `thermolib3d.f90` (verbatim), `sfcphy3d_extract.f90` (`sfcflx`, `dragcoef`, `cuc`, `cptc`, `cucwtr`, `cptcwtr`, `cuneuwtr`, `cptneuwtr` extracted from ARPS `sfcphy3d.f90`), `tmix3d_extract.f90` (`tridiag2`), and `include/*.inc`.
- Outside the park (`pv_park = 0`), `pvland_sfcphy` reproduces ARPS `sfcphy = 4` (dragcoef + sfcflx + soilebm): wind `wdspsfc`, air `tair`/`qvair` averaged over levels 1-2, drag coefficients from `tsoil`/`qvsfc`, fluxes to the atmosphere at t^n, `nsfcst` soil sub-steps, fluxes weighted by `stypfrct`, `qvsfc` updated (checked on test_sylvain: equal to sfcphy = 4 to rounding). In the offline model `nsfcst = 1` and the main loop count is `nstep_frc`.
- `pvland_sfcphy(..., meas_opt, S_diffuse_in, lambda_in, phi_in, L_obk_in)`: ARPS passes `meas_opt = 0` (diffuse SW from `split_sw`, lambda/phi from `compute_lambda_phi`, Obukhov length from the previous-step fluxes, bounded to zeta >= -5); the offline model passes `meas_opt = 1` (measured values; the zeta bound still applies). The offline output reads `uref`, `cdh_can`, `cdh_PVa` from module `pvland_sfcphy_work`.
- `offline/pvland_offline_main.f90` and `offline/pvland_park_map.f90`: the only hand-written offline-specific Fortran.
- `legacy/`: old offline sources, not built.

Never edit `core/` or `arps_lib/` by hand. Change ARPS, then `make sync-arps`. `make check-arps` (exit 1 on difference) verifies everything matches ARPS.

Host-specific hook: `INIT_PVLAND` calls `pvland_park_map(nx,ny,pv_park,tem1)`; ARPS provides an MPI/grid version (`src/arps/pvland_park_map.f90`, driven by `pv_park_opt`), the offline one applies the same `pv_park_opt` to the computed column (1,1): 0 = park, 1 = park if (1,1) is in the rectangle, 2 = no park (natural ground only, as ARPS `sfcphy = 4`; the main program then provides the net radiation `rnflx` from `albedo_nat`/`emiss_nat` and the wind `wdspsfc`, as ARPS does; the PV columns of the output files are then meaningless).

## Build and run

```bash
make              # build pv_land.exe (objects and .mod in build/)
make run          # full forcing, ~2 min
make run-short    # 5-day forcing recorded_data_short_5days.csv, ~15 s
make check-arps   # compare with ../arps5.3.4_PVLAND (ARPS_DIR=... to override)
make sync-arps    # copy from ARPS
make rebuild
```

- Compiler `ifx`, `FFLAGS = -O2 -traceback` (debug line `-g -traceback -check bounds` kept commented). REAL is single precision, as in ARPS. gfortran flags are commented out; gfortran needs `MODFLAG = -J$(BUILDDIR)` and `-llapack -lblas`.
- LAPACK (`SGESV`, radiosity solvers) comes from Intel MKL (`-qmkl=sequential`); the Makefile hardcodes `MKLROOT`/`INTEL_LIB` under `/opt/intel/oneapi/.../latest` with an rpath, since `setvars.sh` is often not sourced.
- Compiler warnings about `cdm_PVb`/`cdm_PVf` (INTENT(OUT) not set in `DRAG_PVLAND`) and subnormal constants in `pvland_profil.f90` come from the ARPS code and are known.
- On the cluster: `sbatch batch_curta.sh`.
- No test suite. Reference outputs of the pre-restructuring model are kept locally (git-ignored) in `ref_8d92348/{full,short_5days}{,_O2}/`. Validate a change with `make run-short` and compare `pv_land_output.csv` with the reference, then with the `*_meas` columns.

## Inputs / outputs

- Forcing: first command-line argument, default `recorded_data_adjusted_daily_summer_2025.csv`. Columns read by position: `t, wind_speed, wind_dir, T6 (°C), specific_humidity, Rain, SWout, SWin, LWout, LWin, air_pressure, air_density, H, LE, G, RN, SW_diffus, lambda, phi, L (Obukhov length)`. Rows 30 min apart. `lambda`/`phi` (rad, park frame) are precomputed by the standalone `calcul_angle_zen_azi_sylvain_v2.f90`.
- Outputs: `pv_land_output.csv` and `pv_diag_energy_balance.csv` (same columns as before the restructuring).
- PV park parameters: namelist `&pvland` (module `pvland_nml` in `core/module_pvland.f90`, which holds the defaults and a one-line description of each parameter). Offline: file `pvland.nml` (2nd command-line argument, default `pvland.nml`), read by `pvland_offline_main` before `INIT_PVLAND`. ARPS: block `&pvland` of `arps.input` after `&soil_ebm`, read by `initpara` only when `sfcphy = 5`. `INIT_PVLAND` copies the scalars into the `(nx,ny)` arrays of `module_pvland_in`. `pv_park_opt`/`pv_park_i1..j2` place the park on the ARPS grid; offline, `pv_park_opt = 2` runs the column without PV panels. Everything derived from the tilt `beta` (h_clear, width, lambda_f_park, z_can, z_PVa, view factors, radiosity LU factors) is computed by `pvland_geometry`: once in `INIT_PVLAND` when `vf_update_opt = 0`, or at every `pvland_sfcphy` step when `vf_update_opt = 1` (moving panels; no tracking law for `beta` is implemented yet). Site parameters (z0_nat, soil/vegetation types, initial soil state, `dtsfc = 2 s`): top of `offline/pvland_offline_main.f90`.

## Offline driver

`pvland_offline_main` builds an ARPS-like single column on a `nx = ny = 2`, `nz = 3` grid (ARPS loops run 1..nx-1, so only column (1,1) is computed; scalar level 2 is the reference level at 3.3·h_PV). Per 30-min forcing interval it sub-steps `nsfcst = 1800/dtsfc` times: linear interpolation of the forcing, `u = U sin(dir)`, `v = U cos(dir)` (so `DRAG_PVLAND` recovers `dir = ATAN2(u,v)`), then `pvland_sfcphy` (meas_opt = 1), which is the ARPS sequence: `DRAG_PVLAND` (wind profile, exchange coefficients), `solve_air_node_temperatures`, RK2 predictor (dt/2)/corrector as in `soilebm` with `compute_all_forcings` at both stages, final diagnostics (`H_ref`, ...). State lives in ARPS arrays (`tsoil`, `qsoil`, `wetcanp`, `qvsfc`) and in the modules `module_pvland_in/out` (`T_PVf`, `T_can`, ...).

Behaviour inherited from ARPS: `uref = MAX(vsfcmin = 0.1 m/s, |V|)` in `DRAG_PVLAND`.

Conventions inherited from ARPS: upper-case Fortran, mostly implicit interfaces (argument lists must match exactly), global constants through `INCLUDE`d common blocks.

## Version control

Git repository, branch `main`, private GitHub repo `sylllvain/pv-land` (remote `origin` over SSH; `gh` is not installed). `.gitignore` excludes `build/`, build products, `slurm-*.out`/`*.log`, run outputs (`pv_land_output.csv`, `pv_diag_energy_balance.csv` in any directory, `*.plt`), `*.zip`, `pdf/`, `*.swp`. Input data are tracked.
