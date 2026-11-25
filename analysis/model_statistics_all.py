#!/usr/bin/env python3
"""
Model performance statistics for PV-LAND
=========================================
Computes Bias, RMSE, MAE, nRMSE(%), R² for each variable pair
(measured vs modeled), split into: All / Day / Night

Day/Night criterion: SWin > 10 W/m²

Outlier filtering: rate-of-change filter applied to the MEASURED fluxes
only (eddy-covariance H/LE/Rn/G). The model output is never filtered:
large half-hour-to-half-hour changes in the model are a legitimate
response to the forcing (e.g. broken-cloud shortwave) and excluding them
would bias the statistics in the model's favour. A timestep is flagged
as a measurement spike if the absolute change between consecutive steps
exceeds SPIKE_SIGMA * std(diff) for that measured series. Flagged points
are printed to the console and excluded from statistics (both the
measured and the corresponding modelled value are dropped for that
variable, so the comparison stays paired).

Usage:  python3 model_statistics.py

Input files (same directory):
  - pv_land_output.csv
  - measured_data_summer.csv
  - recorded_data_adjusted_summer_2025.csv

Output:
  - Terminal table
  - model_statistics_output.csv
"""

import numpy as np
import pandas as pd
import sys

# =====================================================================
# PARAMETERS
# =====================================================================
SW_THRESHOLD = 10.0    # W/m², day/night criterion
ABS_MAX      = 9000.0  # filter extreme absolute values (keep as safety net)

# Time window — set to None to use the full data range
T_START = None     # first timestep to include (inclusive), or None for auto
T_END   = None     # last  timestep to include (inclusive), or None for auto

# Wind-direction filter for H and LE measurements ONLY
# Measurements acquired with wind_dir in [WIND_DIR_MIN, WIND_DIR_MAX] (degrees)
# are excluded from H and LE statistics (flow distortion / bad footprint).
# Set FILTER_VARS = set() to disable.
WIND_DIR_MIN = 175.0
WIND_DIR_MAX = 255.0
FILTER_VARS  = {"H", "LE"}   # variable names (as defined in `variables` below)

# Rate-of-change spike filter: flag timestep t if |x[t] - x[t-1]| > SPIKE_SIGMA * std(diff)
SPIKE_SIGMA  = 5.0     # number of std deviations to flag as spike

# Columns to apply the rate-of-change filter on.
# MEASURED fluxes only: the eddy-covariance signal contains isolated,
# non-physical jumps (sensor icing/condensation, intermittent turbulence
# at low wind, etc.) that a rate-of-change test can reasonably catch.
# The model output is deliberately NOT included here (see module
# docstring): a modelled spike is not necessarily an artefact, so
# despiking it before computing statistics would not be defensible.
SPIKE_COLS = ["H_meas", "LE_meas", "G_meas", "Rn_meas"]

# For each measured column that is flagged, the paired model column is
# also excluded from the statistics for that timestep, so that Bias/RMSE/
# R2 are always computed on the same set of paired points. This does NOT
# alter the model series itself, nor is it used to detect model spikes.
PAIRED_MODEL_COL = {
    "H_meas":  "H_model",
    "LE_meas": "LE_model",
    "G_meas":  "G_model",
    "Rn_meas": "Rn_model",
}

# =====================================================================
# READ DATA
# =====================================================================
print("Reading pv_land_output.csv ...")
df_mod = pd.read_csv("pv_land_output_2025.csv")
for col in df_mod.columns:
    df_mod[col] = pd.to_numeric(df_mod[col], errors="coerce")
    df_mod = df_mod.replace(-9999, np.nan)
print(f"  -> {len(df_mod)} timesteps read.")

print("Reading measured_data_summer.csv ...")
df_meas = pd.read_csv("measured_data_2025.csv")
for col in df_meas.columns:
    df_meas[col] = pd.to_numeric(df_meas[col], errors="coerce")
print(f"  -> {len(df_meas)} timesteps read.")

print("Reading recorded_data_summer_2025.csv ...")
df_forc = pd.read_csv("recorded_data_2025.csv")
for col in df_forc.columns:
    df_forc[col] = pd.to_numeric(df_forc[col], errors="coerce")
print(f"  -> {len(df_forc)} timesteps read.")

# =====================================================================
# MERGE ON TIMESTEP
# =====================================================================
df = df_mod.merge(df_forc[["t", "SWin"]], left_on="timestep", right_on="t",
                  how="inner", suffixes=("", "_forc"))

df = df.merge(
    df_meas[["t", "T_PVf", "T_PVb", "T_nat", "wind_dir"]].rename(
        columns={"T_PVf": "T_PVf_meas", "T_PVb": "T_PVb_meas", "T_nat": "T_nat_meas"}
    ),
    left_on="timestep", right_on="t", how="left", suffixes=("", "_meas2")
)

print(f"\nMerged dataset: {len(df)} timesteps")
print(f"  with measured temperatures: {df['T_PVf_meas'].notna().sum()} timesteps")

# =====================================================================
# TIME WINDOW FILTER
# =====================================================================
if T_START is not None or T_END is not None:
    n_before = len(df)
    if T_START is not None:
        df = df[df["timestep"] >= T_START]
    if T_END is not None:
        df = df[df["timestep"] <= T_END]
    df = df.reset_index(drop=True)
    ts_min = int(df["timestep"].min()) if len(df) else None
    ts_max = int(df["timestep"].max()) if len(df) else None
    print(f"\nTime window applied: T_START={T_START}, T_END={T_END}")
    print(f"  -> {len(df)} timesteps kept (from {n_before})"
          f"{f', range [{ts_min}, {ts_max}]' if ts_min is not None else ''}")
    if len(df) == 0:
        print("ERROR: empty dataset after time-window filtering.")
        sys.exit(1)

# =====================================================================
# RATE-OF-CHANGE SPIKE FILTER
# =====================================================================
# Per-column spike masks (measured columns only). Union kept only for
# reporting the total count; the exclusion itself is applied per pair
# below, not with a single global mask.
spike_masks = {}   # col -> boolean Series, True = spike

print()
print("=" * 70)
print("  SPIKE DETECTION on MEASURED fluxes only")
print("  (rate-of-change > {:.0f} × std of diff)".format(SPIKE_SIGMA))
print("=" * 70)

any_spike_found = False

for col in SPIKE_COLS:
    if col not in df.columns:
        continue
    series = df[col].values.astype(float)
    diff   = np.diff(series, prepend=series[0])   # diff[i] = series[i] - series[i-1]
    std_diff = np.nanstd(diff)
    if std_diff < 1e-10:
        spike_masks[col] = pd.Series(False, index=df.index)
        continue
    threshold = SPIKE_SIGMA * std_diff
    is_spike  = np.abs(diff) > threshold
    spike_masks[col] = pd.Series(is_spike, index=df.index)
    n_spikes  = is_spike.sum()
    if n_spikes > 0:
        timesteps_flagged = df.loc[is_spike, "timestep"].values
        values_flagged    = series[is_spike]
        diffs_flagged     = diff[is_spike]
        print(f"\n  Column '{col}': {n_spikes} spike(s) detected (threshold = ±{threshold:.1f} W/m²)")
        for ts, val, dv in zip(timesteps_flagged, values_flagged, diffs_flagged):
            print(f"    timestep={int(ts):5d}  value={val:9.2f}  Δ={dv:+9.2f}")
        any_spike_found = True

if not any_spike_found:
    print("  No spikes detected.")

n_excluded = int(pd.concat(spike_masks.values(), axis=1).any(axis=1).sum()) if spike_masks else 0
print(f"\n  Total timesteps with at least one flagged measured value: {n_excluded} / {len(df)}")
print("=" * 70)

# Apply: a flagged measured value is excluded from the comparison for
# that variable, together with the paired model value at the same
# timestep (so Bias/RMSE/R2 stay computed on matched pairs). Nothing
# else is touched: the model series itself is never modified or
# excluded on its own account, and measured/modelled values for the
# OTHER variables at that same timestep are left untouched.
df_filtered = df.copy()
for meas_col, is_spike in spike_masks.items():
    mod_col = PAIRED_MODEL_COL.get(meas_col)
    df_filtered.loc[is_spike, meas_col] = np.nan
    if mod_col is not None and mod_col in df_filtered.columns:
        df_filtered.loc[is_spike, mod_col] = np.nan

# =====================================================================
# DEFINE VARIABLE PAIRS
# =====================================================================
variables = [
    # (name, unit, obs_column, mod_column)
    ("H",     "W/m2",  "H_meas",     "H_model"),
    ("LE",    "W/m2",  "LE_meas",    "LE_model"),
    ("Rn",    "W/m2",  "Rn_meas",    "Rn_model"),
    ("G",     "W/m2",  "G_meas",     "G_model"),
    ("T_PVf", "degC",  "T_PVf_meas", "T_PVf"),
    ("T_PVb", "degC",  "T_PVb_meas", "T_PVb"),
    ("T_nat", "degC",  "T_nat_meas", "T_nat"),
]

# =====================================================================
# STATISTICS FUNCTION
# =====================================================================
def compute_stats(obs, mod):
    """Compute Bias, RMSE, MAE, nRMSE(%), R² from obs and mod arrays."""
    n = len(obs)
    if n == 0:
        return {"N": 0, "Mean_obs": np.nan, "Mean_mod": np.nan,
                "Bias": np.nan, "RMSE": np.nan, "MAE": np.nan,
                "nRMSE_pct": np.nan, "R2": np.nan}

    mean_obs = np.mean(obs)
    mean_mod = np.mean(mod)
    diff     = mod - obs
    bias     = np.mean(diff)
    rmse     = np.sqrt(np.mean(diff**2))
    mae      = np.mean(np.abs(diff))

    if abs(mean_obs) > 1e-10:
        nrmse = 100.0 * rmse / abs(mean_obs)
    else:
        nrmse = np.nan

    ss_tot = np.sum((obs - mean_obs)**2)
    ss_res = np.sum(diff**2)
    if ss_tot > 1e-10:
        r2 = 1.0 - ss_res / ss_tot
    else:
        r2 = np.nan

    return {"N": n, "Mean_obs": mean_obs, "Mean_mod": mean_mod,
            "Bias": bias, "RMSE": rmse, "MAE": mae,
            "nRMSE_pct": nrmse, "R2": r2}

# =====================================================================
# COMPUTE AND DISPLAY
# =====================================================================
periods = [
    ("All",   df_filtered["SWin"].notna()),
    ("Day",   df_filtered["SWin"] > SW_THRESHOLD),
    ("Night", df_filtered["SWin"] <= SW_THRESHOLD),
]

results = []

# =====================================================================
# WIND-DIRECTION FILTER (pre-compute the boolean mask once)
# =====================================================================
# True where wind direction is inside the excluded sector.
# NaN wind_dir → treated as valid (not excluded).
if "wind_dir" in df_filtered.columns:
    wd = df_filtered["wind_dir"]
    in_excluded_sector = (wd >= WIND_DIR_MIN) & (wd <= WIND_DIR_MAX)
else:
    in_excluded_sector = pd.Series(False, index=df_filtered.index)

n_wind_excluded = int(in_excluded_sector.sum())

hdr = f"{'Variable':>10s}{'Unit':>7s}{'Period':>8s}{'N_val':>8s}" \
      f"{'Mean_obs':>10s}{'Mean_mod':>10s}{'Bias':>10s}{'RMSE':>10s}" \
      f"{'MAE':>10s}{'nRMSE(%)':>11s}{'R2':>10s}"
sep = "-" * len(hdr)

print()
print("=" * len(hdr))
print("                  PV-LAND MODEL PERFORMANCE STATISTICS")
print("=" * len(hdr))
print(f"  Day criterion: SWin > {SW_THRESHOLD:.0f} W/m2")
print(f"  Spike filter:  rate-of-change > {SPIKE_SIGMA:.0f} × std(diff), "
      f"applied to MEASURED fluxes only ({n_excluded} timesteps affected); "
      f"model output is never filtered.")
if T_START is not None or T_END is not None:
    tw_start = T_START if T_START is not None else "auto"
    tw_end   = T_END   if T_END   is not None else "auto"
    print(f"  Time window:   T_START={tw_start}, T_END={tw_end}  "
          f"({len(df)} timesteps)")
if FILTER_VARS:
    print(f"  Wind-dir filter: {','.join(sorted(FILTER_VARS))} excluded when "
          f"wind_dir ∈ [{WIND_DIR_MIN:.0f}°, {WIND_DIR_MAX:.0f}°]  "
          f"({n_wind_excluded} timesteps in sector)")
print()
print(hdr)
print(sep)

for vname, vunit, obs_col, mod_col in variables:
    for pname, pmask in periods:
        mask = pmask.copy()
        if obs_col in df_filtered.columns and mod_col in df_filtered.columns:
            mask = mask & df_filtered[obs_col].notna() & df_filtered[mod_col].notna()
            mask = mask & (df_filtered[obs_col].abs() < ABS_MAX) & (df_filtered[mod_col].abs() < ABS_MAX)
            # Wind-direction filter (only for variables in FILTER_VARS)
            if vname in FILTER_VARS:
                mask = mask & (~in_excluded_sector)
        else:
            mask = pd.Series(False, index=df_filtered.index)

        obs_arr = df_filtered.loc[mask, obs_col].values
        mod_arr = df_filtered.loc[mask, mod_col].values

        s = compute_stats(obs_arr, mod_arr)

        nrmse_str = f"{s['nRMSE_pct']:11.2f}" if not np.isnan(s['nRMSE_pct']) else "        N/A"
        r2_str    = f"{s['R2']:10.4f}" if not np.isnan(s['R2']) else "       N/A"
        print(f"{vname:>10s}{vunit:>7s}{pname:>8s}{s['N']:8d}"
              f"{s['Mean_obs']:10.2f}{s['Mean_mod']:10.2f}{s['Bias']:10.3f}"
              f"{s['RMSE']:10.3f}{s['MAE']:10.3f}{nrmse_str}{r2_str}")

        results.append({
            "Variable": vname, "Unit": vunit, "Period": pname,
            **s
        })

    print(sep)

# =====================================================================
# WRITE CSV
# =====================================================================
df_out = pd.DataFrame(results)
df_out.to_csv("model_statistics_output_2025.csv", index=False, float_format="%.6f")
print()
print("Results written to: model_statistics_output.csv")
print()
