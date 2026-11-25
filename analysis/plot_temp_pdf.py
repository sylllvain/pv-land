#!/usr/bin/env python3
"""
Simple plotting script — generates plot_temp_summer.pdf only.
Three panels: T_PVf, T_PVb, T_nat (model vs meas), labelled a), b), c).
Single shared legend below the last panel.
"""

import pandas as pd
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
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime, timedelta

# ──────────────────────────────────────────────
# Style configuration — tweak freely
# ──────────────────────────────────────────────
FS_TITLE    = 18   # figure suptitle
FS_AXLABEL  = 26   # x/y axis labels
FS_TICK     = 26   # tick numbers
FS_LEGEND   = 24   # legend text
FS_PANEL    = 26   # panel label a), b), c) above each panel
LW_LINE     = 4  # data line width
LW_GRID     = 0.8  # grid line width
LABELPAD_X  = 15   # space between x-axis label and tick numbers (points)
MARKER_SIZE = 25   # scatter marker size

# Vertical position of the panel label, in axes coordinates
# (1.0 = top of the panel frame; increase to push it further up)
PANEL_LABEL_Y = 1.03
# Vertical gap between panels, as a fraction of the panel height.
# Needs to be large enough for the panel labels to fit in.
PANEL_HSPACE  = 0.25 

# Time window — set to None to use the full data range
T_START = 337     # first timestep to plot (e.g. 500), or None for auto
T_END   = 816     # last  timestep to plot (e.g. 1200), or None for auto

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
# plot_temp_summer.pdf
# ══════════════════════════════════════════════
# Per-panel configuration:
#   label        : panel label drawn above the panel, e.g. "a)".
#                  Append text if a title is wanted, e.g. "a) front surface".
#   meas_cols    : list of measurement column names (1 or more)
#   meas_colors  : list of matplotlib colors, same length as meas_cols
#   meas_labels  : list of legend labels, same length as meas_cols
#   model_col    : model column name (from pv_land_output.csv), or None
#   ylabel       : y-axis label
#   ylim         : (ymin, ymax) or None for auto-scaling
#   yticks       : list of y tick positions
#   pad          : y-label padding (points)
temps = [
    {
        "label":       "a)",
        "meas_cols":   ["T_PVf"],
        "meas_colors": ["tab:blue"],
        "meas_labels": ["Measure"],
        "model_col":   "T_PVf",
        "ylabel":      r"$T_\mathrm{PVf}$ ($^\circ$C)",
        "ylim":        None,
        "yticks":      [0, 15, 30, 45, 60],
        "pad":         12,
    },
    {
        "label":       "b)",
        "meas_cols":   ["T_PVb"],
        "meas_colors": ["tab:blue"],
        "meas_labels": ["Measure"],
        "model_col":   "T_PVb",
        "ylabel":      r"$T_\mathrm{PVb}$ ($^\circ$C)",
        "ylim":        None,
        "yticks":      [0, 15, 30, 45, 60],
        "pad":         12,
    },
    {
        "label":       "c)",
        "meas_cols":   ["T_nat1", "T_nat2"],
        "meas_colors": ["goldenrod", "saddlebrown"],      # tweak colors here
        "meas_labels": [r"$T_\mathrm{nat,interrow}$ (measure)", r"$T_\mathrm{nat,under\ PV}$ (measure)"],
        "model_col":   "T_nat",
        "ylabel":      r"$T_\mathrm{nat}$ ($^\circ$C)",
        "ylim":        None,
        "yticks":      [0, 10, 20, 30, 40],
        "pad":         12,
    },
]

with PdfPages("plot_temp_2025.pdf") as pdf:
    current = t_min
    while current <= t_max:
        end = min(current + STEPS_PER_PAGE - 1, t_max)

        mask_mod = (t_model >= current) & (t_model <= end)
        mask_mea = (t_meas  >= current) & (t_meas  <= end)
        dc = df_model.loc[mask_mod]
        dm = df_meas.loc[mask_mea]

        if dc.empty and dm.empty:
            current = end + 1
            continue

        tc = dc["timestep"] if not dc.empty else pd.Series(dtype=float)
        tm = dm["t"]        if not dm.empty else pd.Series(dtype=float)

        d1 = (current - t_min) // STEPS_PER_DAY + 1
        d2 = (end - t_min) // STEPS_PER_DAY + 1

        fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(18, 12))
        #fig.suptitle(f"Temperatures — timesteps {current}–{end} (days {d1}–{d2})",
        #             fontsize=FS_TITLE, fontweight="bold")

        for ax, p in zip(axes, temps):
            # Measurement scatter(s) — one or more series
            if not dm.empty:
                for col, color in zip(p["meas_cols"], p["meas_colors"]):
                    if col in dm.columns:
                        ax.scatter(tm, dm[col], s=MARKER_SIZE, color=color,
                                   alpha=0.6, zorder=3)
            # Model line
            if p["model_col"] is not None and not dc.empty and p["model_col"] in dc.columns:
                ax.plot(tc, dc[p["model_col"]], lw=LW_LINE, color="tab:orange")
            ax.set_ylabel(p["ylabel"], labelpad=p["pad"])
            if p["ylim"] is not None:
                ax.set_ylim(p["ylim"])
            ax.set_yticks(p["yticks"])
            ax.grid(True, alpha=0.3, lw=LW_GRID)
            # Panel label above the frame, left aligned with the y axis
            ax.text(0.0, PANEL_LABEL_Y, p["label"], transform=ax.transAxes,
                    fontsize=FS_PANEL, fontweight="normal",
                    va="bottom", ha="left")

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
        mid_labels = [(T_ORIGIN + t * timedelta(minutes=DT_MINUTES)).strftime(DATE_FMT)
                      for t in boundaries[:-1]]

        ax_x = axes[-1]
        ax_x.set_xticks(tick_positions)
        ax_x.set_xticklabels([""] * len(tick_positions))
        ax_x.set_xticks(mid_positions, minor=True)
        ax_x.set_xticklabels(mid_labels, minor=True)
        ax_x.tick_params(axis="x", which="minor", length=0)
        ax_x.set_xlim(current, end)

        # Light vertical day separators on every panel
        for ax in axes:
            for v in tick_positions:
                ax.axvline(v, color="grey", ls=":", lw=0.5, alpha=0.5)

        axes[-1].set_xlabel("Date", labelpad=LABELPAD_X)

        # Single legend outside, below the three panels.
        # Collect the distinct (label, color) pairs across all meas series,
        # preserving order of first appearance, then add the model.
        legend_entries = []
        seen = set()
        for p in temps:
            for lbl, col in zip(p["meas_labels"], p["meas_colors"]):
                if lbl not in seen:
                    legend_entries.append((lbl, col))
                    seen.add(lbl)
        legend_handles = [
            plt.Line2D([0], [0], marker="o", color=col, linestyle="None",
                       markersize=11, label=lbl)
            for lbl, col in legend_entries
        ]
        legend_handles.append(
            plt.Line2D([0], [0], color="tab:orange", lw=LW_LINE, label="Model")
        )
        fig.legend(handles=legend_handles, loc="lower center",
                   ncol=len(legend_handles),
                   fontsize=FS_LEGEND, frameon=False, bbox_to_anchor=(0.5, 0.004))
        fig.tight_layout(rect=[0, 0.08, 1, 0.99])
        # Applied after tight_layout, which would otherwise reset it:
        # leaves room above each panel for the a)/b)/c) labels.
        fig.subplots_adjust(hspace=PANEL_HSPACE)
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        current = end + 1

print("Saved: plot_temp_2025.pdf")
