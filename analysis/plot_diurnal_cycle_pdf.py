#!/usr/bin/env python3
"""
plot_diurnal_cycle_pdf.py
=========================
Mean diurnal cycle of the measured and simulated turbulent fluxes (H, LE)
over the 2025 evaluation period -> plot_diurnal_cycle.pdf
(Sect. "Surface energy flux validation", Fig. "fig:diurnal").

Layout: 2 rows x 1 column
    rows    : H, LE
    column  : 2025
(The layout follows the content of FILES: with several years, one column per
year is drawn.)
Gray bands: the two windows discussed in the text (06:00-10:00, 15:00-19:00).

The values quoted in the text still come from the statistics script. This
script draws the figure and prints the window means of (model - measurement),
so that the figure can be checked against the text. The line "UTC+0 center"
must reproduce the windows printed by compute_balance_stats.py.

INPUT FILES (same files as the other figure scripts)
    pv_land_output_2025.csv   model output; columns timestep, H_meas, H_model,
                              LE_meas, LE_model, Rn_meas, G_meas, ...
                              (H_meas and LE_meas are the closure-adjusted fluxes
                              copied by the driver from the forcing file, row for
                              row, so model and measurement share the same timestep)
    measured_data_2025.csv    measured data; only t and wind_dir are used here

FILTERS (same as Sect. "Error metrics" of the article)
    - missing values (-9999, or |x| > 1e4) are set to NaN when the files are read,
      before any mean, difference or standard deviation is computed;
    - H and LE: half-hours with wind from WIND_DIR_MIN to WIND_DIR_MAX removed
      (pine stand); missing wind direction is kept, as in the statistics.
    - rate-of-change filter on the MEASURED fluxes only, each series on its own
      (as in model_statistics_all.py): a half-hour of H (resp. LE) is removed
      when |x[t] - x[t-1]| > SPIKE_SIGMA * std(increments) for H_meas (resp.
      LE_meas). The model value is removed with it, so that means are always
      computed on matched pairs.

TIME CONVENTION
    Records are 30-min averages time-stamped in UTC at the END of the interval.
    Timestep t = 1 is the interval 00:00-00:30 UTC on 2 July: with this mapping
    the clear-day centroid of the incoming shortwave radiation falls within a
    few minutes of solar noon (about 12:07 UTC at Salaunes).
    STAMP = "center" places a record at the center of its interval,
    STAMP = "end" at its end; TZ_OFFSET_H shifts UTC to local time (2 = CEST).

PAIRING OF MODEL AND MEASUREMENT
    The driver writes, at timestep t, the model state reached at the end of the
    interval AND the measured flux of the same record t (read from the same row
    of the forcing file). The forcing (temperature, wind, radiation...) and the
    measured fluxes are both 30-min averages time-stamped at the end of the
    interval, and the model treats each forcing record as the value at its time
    stamp, so model and measurement are compared on the same footing: no shift
    is applied. Averaging the model over the interval would add a spurious
    delay and must not be done.

USAGE
    python3 plot_diurnal_cycle_pdf.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# --- Computer Modern, the font of the thesis. To go back to the default font, delete this block.
import dataclasses as _dc
from matplotlib import font_manager as _fm
for _f in list(_fm.fontManager.ttflist):        # bold text (fontweight="bold") uses cmb10
    if _f.name == "cmb10":
        _fm.fontManager.ttflist.append(_dc.replace(_f, name="cmr10", weight=700))
plt.rcParams.update({
    "font.family": ["cmr10", "DejaVu Serif"],   # DejaVu only for characters that cmr10 lacks (en dash, arrows...)
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,        # proper minus sign on the tick labels
    "axes.unicode_minus": False,
})
# Do not type the characters ° ² ³ µ ± in a label: cmr10 prints them as ligatures (fl, ffl...).
# Use $^\circ$, $^2$... instead.
# --- end of font block

# ──────────────────────────────────────────────
# Input / output
# ──────────────────────────────────────────────
FILES = {
    2025: ("pv_land_output_2025.csv", "measured_data_2025.csv"),
}
OUT_PDF = "plot_diurnal_cycle.pdf"

# ──────────────────────────────────────────────
# Time convention (see docstring)
# ──────────────────────────────────────────────
TZ_OFFSET_H = 0.0          # 0 = UTC (convention of the article), 2 = CEST
STAMP = "center"           # "center" or "end"
WINDOWS = [(6.0, 10.0), (15.0, 19.0)]   # [start, end) in hours

# ──────────────────────────────────────────────
# Filters (must match the statistics of the article)
# ──────────────────────────────────────────────
WIND_DIR_MIN = 175.0
WIND_DIR_MAX = 255.0
WIND_FILTER_VARS = {"H", "LE"}
SPIKE_SIGMA = 5.0
SPIKE_COLS = ["H_meas", "LE_meas", "Rn_meas", "G_meas"]
MISSING_VALUE = -9999.0    # sentinel written for missing data in the CSV files
BAD_VALUE_ABS = 1e4        # any |x| above this is also treated as missing

# ──────────────────────────────────────────────
# Style
# ──────────────────────────────────────────────
FS_AXLABEL = 24
FS_TICK = 21
FS_LEGEND = 22
FS_PANEL = 23           # panel letters a), b)... (normal weight, above the panel, left)
FS_COLTITLE = 25        # column title (year), bold, above the panel letter
SHOW_COLTITLE = True    # False: no year above the first panel (the caption says it)
PANEL_PAD = 8           # points between the top of the axes and the panel letter
# The column title sits above the panel letter: its distance to the top of the
# axes must exceed the height of the letter line (pad + one line of text).
COLTITLE_PAD = PANEL_PAD + 1.25 * FS_PANEL + 4
LW_LINE = 3.5
MARKER_SIZE = 55
COL_MEAS = "tab:blue"
COL_MODEL = "tab:orange"
WINDOW_COLOR = "0.85"
WIDTH_PER_COLUMN = 8.5   # inches per year column (the original two-column figure was 15 wide)
HEIGHT = 9.5

plt.rcParams.update({
    "axes.labelsize": FS_AXLABEL,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "legend.fontsize": FS_LEGEND,
})


# ══════════════════════════════════════════════
# Data
# ══════════════════════════════════════════════
def read_numeric(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce")
    # Missing values (-9999 sentinels, absurd values) become NaN HERE, before any
    # statistic: otherwise they would be averaged, and they would also inflate the
    # standard deviation of the increments and disable the rate-of-change filter.
    # Note that abs(-9999) < 1e4, so the sentinel needs its own test.
    df = df.mask((df == MISSING_VALUE) | (df.abs() > BAD_VALUE_ABS))
    return df


def load_year(model_csv, meas_csv):
    mod = read_numeric(model_csv)
    meas = read_numeric(meas_csv)
    df = mod.merge(meas[["t", "wind_dir"]].rename(columns={"t": "timestep"}),
                   on="timestep", how="left")
    return df.sort_values("timestep").reset_index(drop=True)


def valid_masks(df):
    """Boolean mask of the half-hours kept for H and for LE."""
    # Spike test on each measured series ON ITS OWN, as in model_statistics_all.py:
    # a spike on Rn or G does not remove H or LE at the same half-hour.
    def spike_of(col):
        x = df[col].to_numpy(dtype=float)
        dx = np.diff(x, prepend=x[0])
        s = np.nanstd(dx)
        if not (np.isfinite(s) and s > 1e-10):
            return np.zeros(len(df), dtype=bool)
        with np.errstate(invalid="ignore"):
            return np.abs(dx) > SPIKE_SIGMA * s
    wd = df["wind_dir"]
    in_sector = ((wd >= WIND_DIR_MIN) & (wd <= WIND_DIR_MAX)).to_numpy()
    masks = {}
    for v in ("H", "LE"):
        ok = (df[f"{v}_meas"].notna() & df[f"{v}_model"].notna()).to_numpy() & ~spike_of(f"{v}_meas")
        if v in WIND_FILTER_VARS:
            ok &= ~in_sector
        masks[v] = ok
    return masks


def hour_of_day(timestep, tz=TZ_OFFSET_H, stamp=STAMP):
    shift = 0.5 if stamp == "center" else 0.0
    return np.mod((np.asarray(timestep, dtype=float) - shift) * 0.5 + tz, 24.0)


def diurnal(df, ok, var, tz=TZ_OFFSET_H, stamp=STAMP):
    """Mean measured and modeled value per half-hour slot."""
    h = np.round(hour_of_day(df["timestep"], tz, stamp), 3)
    d = pd.DataFrame({"h": h[ok],
                      "meas": df[f"{var}_meas"].to_numpy()[ok],
                      "model": df[f"{var}_model"].to_numpy()[ok]})
    return d.groupby("h").agg(meas=("meas", "mean"), model=("model", "mean"))


def window_means(df, ok, var, tz, stamp):
    h = hour_of_day(df["timestep"], tz, stamp)
    diff = df[f"{var}_model"].to_numpy(dtype=float) - df[f"{var}_meas"].to_numpy(dtype=float)
    out = []
    for a, b in WINDOWS:
        k = ok & (h >= a) & (h < b) & np.isfinite(diff)
        out.append((np.mean(diff[k]) if k.any() else np.nan, int(k.sum())))
    return out


# ══════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════
def main():
    data = {}
    for year, (f_mod, f_meas) in FILES.items():
        df = load_year(f_mod, f_meas)
        data[year] = (df, valid_masks(df))

    # ---- check against the values quoted in the text -------------------
    print("Mean (model - measurement) in the windows, W m-2  [N half-hours]")
    win_lbl = "   ".join(f"{a:05.2f}-{b:05.2f}" for a, b in WINDOWS)
    for tz, stamp in [(0.0, "center"), (0.0, "end"), (2.0, "center"), (2.0, "end")]:
        tag = f"UTC{tz:+.0f} {stamp:<6s}"
        if tz == TZ_OFFSET_H and stamp == STAMP:
            tag += " (used for the figure)"
        print(f"\n  {tag}     windows: {win_lbl}")
        for year, (df, masks) in data.items():
            for v in ("H", "LE"):
                vals = window_means(df, masks[v], v, tz, stamp)
                txt = "   ".join(f"{m:+7.1f} [{n:4d}]" for m, n in vals)
                print(f"    {year} {v:<2s}  {txt}")

    # ---- figure --------------------------------------------------------
    years = list(data.keys())
    # squeeze=False: axes is ALWAYS a 2-D array (rows, columns), even with one
    # single column (np.atleast_2d would turn the 2 axes into a 1 x 2 array).
    fig, axes = plt.subplots(2, len(years), figsize=(WIDTH_PER_COLUMN * len(years), HEIGHT),
                             sharex=True, sharey="row", squeeze=False,
                             layout="constrained")
    letters = "abcdefgh"
    tz_lbl = "UTC" if TZ_OFFSET_H == 0 else f"UTC{TZ_OFFSET_H:+g}"
    for j, year in enumerate(years):
        df, masks = data[year]
        for i, v in enumerate(("H", "LE")):
            ax = axes[i, j]
            c = diurnal(df, masks[v], v)
            for a, b in WINDOWS:
                ax.axvspan(a, b, color=WINDOW_COLOR, zorder=0, lw=0)
            ax.plot(c.index, c["model"], lw=LW_LINE, color=COL_MODEL, zorder=2)
            ax.scatter(c.index, c["meas"], s=MARKER_SIZE, color=COL_MEAS, zorder=3)
            ax.grid(True, alpha=0.3)
            ax.set_title(f"{letters[i * len(years) + j]})", loc="left",
                         fontsize=FS_PANEL, fontweight="normal", pad=PANEL_PAD)
            if i == 0 and SHOW_COLTITLE:
                # column title: bold, centered on the column, above the letters;
                # it applies to every panel of the column
                ax.annotate(str(year), xy=(0.5, 1.0), xycoords="axes fraction",
                            xytext=(0.0, COLTITLE_PAD), textcoords="offset points",
                            ha="center", va="bottom", fontsize=FS_COLTITLE,
                            fontweight="bold", annotation_clip=False)
        axes[-1, j].set_xlim(0, 24)
        axes[-1, j].set_xticks(range(0, 25, 4))
        axes[-1, j].set_xlabel(f"Hour ({tz_lbl})")

    axes[0, 0].set_ylabel(r"$H$ (W m$^{-2}$)")
    axes[1, 0].set_ylabel(r"$LE$ (W m$^{-2}$)")

    handles = [
        plt.Line2D([0], [0], marker="o", color=COL_MEAS, ls="None", markersize=11, label="Measure"),
        plt.Line2D([0], [0], color=COL_MODEL, lw=LW_LINE, label="Model"),
    ]
    # strip kept free under the axes, so that the legend does not touch the "Hour" label
    fig.get_layout_engine().set(rect=(0.0, 0.085, 1.0, 0.915))
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.savefig(OUT_PDF)
    plt.close(fig)
    print(f"\nSaved: {OUT_PDF}")


if __name__ == "__main__":
    main()
