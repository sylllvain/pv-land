#!/usr/bin/env python3
"""
Simple plotting script — generates plot_fluxes_summer.pdf only.
Four panels: H, LE, Rn, G (model vs meas).
Single shared legend below the last panel.
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime, timedelta

# ──────────────────────────────────────────────
# Style configuration — tweak freely
# ──────────────────────────────────────────────
FS_TITLE    = 18   # figure suptitle
FS_AXLABEL  = 28   # x/y axis labels
FS_TICK     = 22   # tick numbers
FS_LEGEND   = 25   # legend text
LW_LINE     = 4  # data line width
LW_GRID     = 0.8  # grid line width
LABELPAD_X  = 15   # space between x-axis label and tick numbers (points)
MARKER_SIZE = 20   # scatter marker size

# Time window — set to None to use the full data range
#T_START = 337     # first timestep to plot (e.g. 500), or None for auto
#T_END   = 816     # last  timestep to plot (e.g. 1200), or None for auto

T_START = None
T_END = None

# Wind-direction filter for H and LE measurements ONLY
# Measurements acquired with wind_dir in [WIND_DIR_MIN, WIND_DIR_MAX] (degrees)
# are highlighted in FILTER_COLOR instead of the default measurement color.
# Set FILTER_COLS = set() to disable.
WIND_DIR_MIN = 175.0
WIND_DIR_MAX = 255.0
FILTER_COLS  = {"H_meas", "LE_meas"}
FILTER_COLOR = "lightgrey"        # color for points inside the excluded sector
FILTER_LABEL = "measure (removed sector)"

# Timestep-to-date mapping.
# Timestep 1 corresponds to T_ORIGIN + DT_MINUTES (i.e. first recorded instant).
# With T_ORIGIN = 2025-07-01 00:30 and DT_MINUTES = 30, t=1 -> 2025-07-01 01:00.
T_ORIGIN   = datetime(2025, 7, 2, 0, 30)
DT_MINUTES = 30
DATE_FMT   = "%d/%m"   # tick label format, e.g. "08/07"

# Apply globally via rcParams so every axes created afterwards inherits them
plt.rcParams.update({
    "axes.titlesize":   FS_AXLABEL,
    "axes.labelsize":   FS_AXLABEL,
    "xtick.labelsize":  FS_TICK,
    "ytick.labelsize":  FS_TICK,
    "legend.fontsize":  FS_LEGEND,
    "figure.titlesize": FS_TITLE,
})
# ──────────────────────────────────────────────

# ──────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────
df_model = pd.read_csv("pv_land_output_2025.csv")
for col in df_model.columns:
    df_model[col] = pd.to_numeric(df_model[col], errors="coerce")
df_model = df_model.mask(df_model.abs() > 1e4)
t_model = df_model["timestep"]

df_meas = pd.read_csv("measured_data_2025.csv")
for col in df_meas.columns:
    df_meas[col] = pd.to_numeric(df_meas[col], errors="coerce")
t_meas = df_meas["t"]

# Align wind_dir with df_model by timestep so that H_meas / LE_meas points
# can be colored differently depending on wind direction.
df_model = df_model.merge(
    df_meas[["t", "wind_dir"]].rename(columns={"t": "timestep"}),
    on="timestep", how="left",
)

STEPS_PER_DAY = 48
DAYS_PER_PAGE = 10
STEPS_PER_PAGE = STEPS_PER_DAY * DAYS_PER_PAGE

t_min = int(min(t_model.min(), t_meas.min()))
t_max = int(max(t_model.max(), t_meas.max()))
if T_START is not None:
    t_min = max(t_min, T_START)
if T_END is not None:
    t_max = min(t_max, T_END)


# ══════════════════════════════════════════════
# plot_fluxes_summer.pdf
# ══════════════════════════════════════════════
# Per-panel configuration:
#   ylim : (ymin, ymax) or None for auto-scaling
#   pad  : spacing between y-axis label and tick numbers (in points)
fluxes = [
    # (meas_col,  model_col,  name, ylabel,       ylim, pad)
    ("Rn_meas", "Rn_model", "Rn", r"$\mathrm{R_n}$", None, 12),
    ("H_meas",  "H_model",  "H",  "H",  (-80, 450), 12),
    ("LE_meas", "LE_model", "LE", "LE", (-50,350), 12),
    ("G_meas",  "G_model",  "G",  "G",  None, 23),
]

with PdfPages("plot_fluxes_2025.pdf") as pdf:
    current = t_min
    while current <= t_max:
        end = min(current + STEPS_PER_PAGE - 1, t_max)

        mask = (t_model >= current) & (t_model <= end)
        dc = df_model.loc[mask]

        if dc.empty:
            current = end + 1
            continue

        tc = dc["timestep"]

        d1 = (current - t_min) // STEPS_PER_DAY + 1
        d2 = (end - t_min) // STEPS_PER_DAY + 1

        fig, axes = plt.subplots(nrows=4, sharex=True, figsize=(18, 12))
        #fig.suptitle(f"Energy fluxes — timesteps {current}–{end} (days {d1}–{d2})",
        #             fontsize=FS_TITLE, fontweight="bold")

        for ax, (meas, model, name, ylabel, ylim, pad) in zip(axes, fluxes):
            # Split measurements into two groups if a wind-dir filter applies
            if meas in FILTER_COLS and "wind_dir" in dc.columns:
                wd = dc["wind_dir"]
                in_sector = (wd >= WIND_DIR_MIN) & (wd <= WIND_DIR_MAX)
                keep = ~in_sector
                # Valid points (default color)
                ax.scatter(tc[keep], dc.loc[keep, meas],
                           s=MARKER_SIZE, color="tab:blue", alpha=0.6, zorder=3)
                # Points inside the excluded sector (highlight color)
                ax.scatter(tc[in_sector], dc.loc[in_sector, meas],
                           s=MARKER_SIZE, color=FILTER_COLOR, alpha=0.6, zorder=3)
            else:
                ax.scatter(tc, dc[meas], s=MARKER_SIZE, color="tab:blue",
                           alpha=0.6, zorder=3)
            ax.plot(tc, dc[model], lw=LW_LINE, color="tab:orange")
            ax.set_ylabel(ylabel, labelpad=pad)
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.grid(True, alpha=0.3, lw=LW_GRID)

        # Build daily ticks: one per day at 01:00 (timestep 1, 49, 97, ...).
        # Major ticks = day boundaries (no labels, just graduations).
        # Minor ticks = mid-day positions (carry the date labels).
        first_day_tick = current + ((1 - current) % STEPS_PER_DAY)
        tick_positions = list(range(first_day_tick, end + 1, STEPS_PER_DAY))

        # Midpoint between two consecutive day boundaries = middle of the day.
        # Also include the last day even if its closing boundary is beyond `end`.
        boundaries = tick_positions + [tick_positions[-1] + STEPS_PER_DAY]
        mid_positions = [(boundaries[i] + boundaries[i + 1]) / 2
                         for i in range(len(boundaries) - 1)]
        # Label each midpoint with the date of the *starting* boundary.
        mid_labels = [(T_ORIGIN + t * timedelta(minutes=DT_MINUTES)).strftime(DATE_FMT)
                      for t in boundaries[:-1]]

        ax_x = axes[-1]
        ax_x.set_xticks(tick_positions)                # major ticks, no labels
        ax_x.set_xticklabels([""] * len(tick_positions))
        ax_x.set_xticks(mid_positions, minor=True)
        ax_x.set_xticklabels(mid_labels, minor=True)
        ax_x.tick_params(axis="x", which="minor", length=0)  # no mark, label only
        ax_x.set_xlim(current, end)

        # Light vertical day separators on every panel
        for ax in axes:
            for v in tick_positions:
                ax.axvline(v, color="grey", ls=":", lw=0.5, alpha=0.5)

        axes[-1].set_xlabel("Date", labelpad=LABELPAD_X)
 
        # Single legend outside, below the four panels.
        # Since meas/model style is shared across rows, use generic handles.
        legend_handles = [
            plt.Line2D([0], [0], marker="o", color="tab:blue", linestyle="None",
                       markersize=(MARKER_SIZE ** 0.8), label="measure"),
            plt.Line2D([0], [0], marker="o", color=FILTER_COLOR, linestyle="None",
                       markersize=(MARKER_SIZE ** 0.8), label=FILTER_LABEL),
            plt.Line2D([0], [0], color="tab:orange", lw=LW_LINE, label="model"),
        ]
        fig.legend(handles=legend_handles, loc="lower center", ncol=3,
                   fontsize=FS_LEGEND, frameon=True, bbox_to_anchor=(0.5, 0))
        fig.tight_layout(rect=[0, 0.08, 1, 0.99])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        current = end + 1

print("Saved: plot_fluxes_summer.pdf")


# ══════════════════════════════════════════════
# plot_wind_summer.pdf
# ══════════════════════════════════════════════
# Same layout/logic as above, but only measured data (wind_dir, wind_speed).
wind_vars = [
    # (col,         name,    ylabel,           ylim, pad)
    ("wind_speed", "wind_speed", "wind speed", None, 12),
    ("wind_dir",   "wind_dir",   "wind dir",   (0, 360), 12),
]

t_min_w = int(t_meas.min())
t_max_w = int(t_meas.max())
if T_START is not None:
    t_min_w = max(t_min_w, T_START)
if T_END is not None:
    t_max_w = min(t_max_w, T_END)

with PdfPages("plot_wind_summer.pdf") as pdf:
    current = t_min_w
    while current <= t_max_w:
        end = min(current + STEPS_PER_PAGE - 1, t_max_w)

        mask = (t_meas >= current) & (t_meas <= end)
        dc = df_meas.loc[mask]

        if dc.empty:
            current = end + 1
            continue

        tc = dc["t"]

        fig, axes = plt.subplots(nrows=len(wind_vars), sharex=True, figsize=(18, 8))

        for ax, (col, name, ylabel, ylim, pad) in zip(axes, wind_vars):
            ax.scatter(tc, dc[col], s=MARKER_SIZE, color="tab:blue",
                       alpha=0.6, zorder=3)
            ax.set_ylabel(ylabel, labelpad=pad)
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.grid(True, alpha=0.3, lw=LW_GRID)

        # Build daily ticks: one per day at 01:00 (timestep 1, 49, 97, ...).
        first_day_tick = current + ((1 - current) % STEPS_PER_DAY)
        tick_positions = list(range(first_day_tick, end + 1, STEPS_PER_DAY))

        boundaries = tick_positions + [tick_positions[-1] + STEPS_PER_DAY]
        mid_positions = [(boundaries[i] + boundaries[i + 1]) / 2
                         for i in range(len(boundaries) - 1)]
        mid_labels = [(T_ORIGIN + t * timedelta(minutes=DT_MINUTES)).strftime(DATE_FMT)
                      for t in boundaries[:-1]]

        ax_x = axes[-1]
        ax_x.set_xticks(tick_positions)
        ax_x.set_xticklabels([""] * len(tick_positions))
        ax_x.set_xticks(mid_positions, minor=True)
        ax_x.set_xticklabels(mid_labels, minor=True)
        ax_x.tick_params(axis="x", which="minor", length=0)
        ax_x.set_xlim(current, end)

        for ax in axes:
            for v in tick_positions:
                ax.axvline(v, color="grey", ls=":", lw=0.5, alpha=0.5)

        axes[-1].set_xlabel("Date", labelpad=LABELPAD_X)

        legend_handles = [
            plt.Line2D([0], [0], marker="o", color="tab:blue", linestyle="None",
                       markersize=(MARKER_SIZE ** 0.8), label="measure"),
        ]
        fig.legend(handles=legend_handles, loc="lower center", ncol=1,
                   fontsize=FS_LEGEND, frameon=True, bbox_to_anchor=(0.5, 0))
        fig.tight_layout(rect=[0, 0.08, 1, 0.99])
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        current = end + 1

print("Saved: plot_wind_summer.pdf")
