# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

PV-LAND offline: a Fortran 90 land-surface model for a ground-mounted photovoltaic park, forced by measured meteorological data. It computes the energy balance of the natural ground surface, the soil and the PV panels (front glass / cell / back), and compares modelled fluxes (H, LE, Rn, G) with measurements. The user writes in French.

**ARPS is the reference.** The same PV-LAND code runs coupled in ARPS (`../arps5.3.4_PVLAND`, `sfcphy = 5`; GitHub `sylllvain/arps_pvland`). This repository is only an offline driver around exact copies of the ARPS sources, so that PV-LAND changes can be moved between both models without manual porting. See `README.md` (French) for the user-facing summary.

## Layout and the "do not edit" rule

- `core/`: PV-LAND physics, verbatim copies of ARPS `src/arps/module_pvland.f90`, `pvland.f90` (`INIT_PVLAND`, `DRAG_PVLAND`, `solve_air_node_temperatures`, `compute_all_forcings`, `compute_lambda_phi`), `pvland_profil.f90`, `pvland_view_factors.f90`, `pvland_solar_radiation.f90`, `pvland_radiosity_SW/LW.f90`, `pvland_pv_thermal.f90`.
- `arps_lib/`: ARPS routines used by PV-LAND: `soilebm3d.f90`, `thermolib3d.f90` (verbatim), `sfcphy3d_extract.f90` (`sfcflx`, `dragcoef`, `cuc`, `cptc`, `cucwtr`, `cptcwtr`, `cuneuwtr`, `cptneuwtr` extracted from ARPS `sfcphy3d.f90`), `tmix3d_extract.f90` (`tridiag2`), and `include/*.inc`.
- `offline/pvland_sfcphy_offline.f90`: one PV-LAND time step. **Generated** by `scripts/sync_from_arps.py` from the `ELSE IF ( sfcphy == 5 )` block of ARPS `sfcphy3d.f90`, inserted into `scripts/pvland_sfcphy_offline.tmpl.f90` (declarations). Only lines tagged `[OFFLINE]` differ from ARPS: measured diffuse SW instead of `split_sw`, measured lambda/phi instead of `compute_lambda_phi`, measured Obukhov length instead of the flux-based one; `[OFFLINE DIAG]` lines save `uref`, `cdh_can`, `cdh_PVa` for the output files.
- `offline/pvland_offline_main.f90` and `offline/pvland_park_map.f90`: the only hand-written offline-specific Fortran.
- `legacy/`: old offline sources, not built.

Never edit `core/`, `arps_lib/` or the generated driver by hand. Change ARPS (or the template / the substitutions in `make_offline_step()` of the sync script), then `make sync-arps`. `make check-arps` (exit 1 on difference) verifies everything matches ARPS.

Host-specific hook: `INIT_PVLAND` calls `pvland_park_map(nx,ny,pv_park,tem1)`; ARPS provides an MPI/grid version (`src/arps/pvland_park_map.f90`), the offline one sets `pv_park = 1`.

## Build and run

```bash
make              # build pv_land.exe (objects and .mod in build/)
make run          # full forcing, ~2 min
make run-short    # 5-day forcing recorded_data_short_5days.csv, ~15 s
make check-arps   # compare with ../arps5.3.4_PVLAND (ARPS_DIR=... to override)
make sync-arps    # copy/generate from ARPS
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
- PV park geometry and optical properties: in `INIT_PVLAND` (`core/pvland.f90`, i.e. in ARPS). Site parameters (z0_nat, soil/vegetation types, initial soil state, `dtsfc = 2 s`): top of `offline/pvland_offline_main.f90`.

## Offline driver

`pvland_offline_main` builds an ARPS-like single column on a `nx = ny = 2`, `nz = 3` grid (ARPS loops run 1..nx-1, so only column (1,1) is computed; scalar level 2 is the reference level at 3.3·h_PV). Per 30-min forcing interval it sub-steps `nsfcst = 1800/dtsfc` times: linear interpolation of the forcing, `u = U sin(dir)`, `v = U cos(dir)` (so `DRAG_PVLAND` recovers `dir = ATAN2(u,v)`), then `pvland_sfcphy_offline`, which is the ARPS sequence: `DRAG_PVLAND` (wind profile, exchange coefficients), `solve_air_node_temperatures`, RK2 predictor (dt/2)/corrector as in `soilebm` with `compute_all_forcings` at both stages, final diagnostics (`H_ref`, ...). State lives in ARPS arrays (`tsoil`, `qsoil`, `wetcanp`, `qvsfc`) and in the modules `module_pvland_in/out` (`T_PVf`, `T_can`, ...).

Behaviour inherited from ARPS: `uref = MAX(vsfcmin = 0.1 m/s, |V|)` in `DRAG_PVLAND`.

Conventions inherited from ARPS: upper-case Fortran, mostly implicit interfaces (argument lists must match exactly), global constants through `INCLUDE`d common blocks.

## Version control

Git repository, branch `main`, private GitHub repo `sylllvain/pv-land` (remote `origin` over SSH; `gh` is not installed). `.gitignore` excludes `build/`, build products, `slurm-*.out`/`*.log`, run outputs (`pv_land_output.csv`, `pv_diag_energy_balance.csv` in any directory, `*.plt`), `*.zip`, `pdf/`, `*.swp`. Input data are tracked.
