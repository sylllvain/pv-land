#!/usr/bin/env python3
"""
Standalone plotting script — generates plot_energy_balance_surfaces_summer.pdf.
Three rows labelled a), b), c): PVf, PVb, nat — each row shows the surface
energy balance with radiative terms grouped as a net radiation Rn, plus the
residual (dashed black line) which corresponds to the energy needed for
temperature variations (C · dT/dt).

Sign convention (positive = into the surface / stored):
  PVf : Rn_PVf + qcond_PVf - H_PVf          = C_PVf · dT/dt   -> residual
        with Rn_PVf = qSW_PVf + qLW_PVf
  PVb : Rn_PVb + qcond_PVb - H_PVb          = C_PVb · dT/dt   -> residual
        with Rn_PVb = qSW_PVb + qLW_PVb
  nat : Rn_nat - H_nat - LE_nat - G_nat     ≈ C_nat · dT/dt   -> residual
        with Rn_nat = (G_SW_nat - J_SW_nat) + (G_LW_nat - J_LW_nat)
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
import numpy as np

# ──────────────────────────────────────────────
# Style configuration — tweak freely
# ──────────────────────────────────────────────
FS_TITLE    = 18   # figure suptitle
FS_AXLABEL  = 26   # x/y axis labels
FS_TICK     = 26   # tick numbers
FS_LEGEND   = 24   # legend text
FS_PANEL    = 26   # panel label a), b), c) above each panel
LW_LINE     = 3.5  # data line width
LW_GRID     = 0.8  # grid line width
LABELPAD_X  = 15   # space between x-axis label and tick numbers (points)
LW_RESIDUAL = 2.5  # residual line width

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

df_diag = pd.read_csv("pv_diag_energy_balance_2025.csv")
for col in df_diag.columns:
    df_diag[col] = pd.to_numeric(df_diag[col], errors="coerce")
df_diag = df_diag.mask(df_diag.abs() > 1e4)

# Merge diag (PV surface fluxes) with model output (nat surface fluxes)
df_bal = pd.merge(
    df_diag[["timestep", "qSW_PVf", "qLW_PVf", "qcond_PVf", "H_PVf",
                         "qSW_PVb", "qLW_PVb", "qcond_PVb", "H_PVb",
                         "qPV"]],
    df_model[["timestep", "G_SW_nat", "J_SW_nat", "G_LW_nat", "J_LW_nat",
                          "H_nat", "LE_model", "G_model"]],
    on="timestep", how="inner",
)

# Net radiation and aliases
df_bal["Rn_PVf"]   = df_bal["qSW_PVf"] + df_bal["qLW_PVf"]
df_bal["Rn_PVb"]   = df_bal["qSW_PVb"] + df_bal["qLW_PVb"]
df_bal["Rn_nat"]   = (df_bal["G_SW_nat"] - df_bal["J_SW_nat"]) + \
                     (df_bal["G_LW_nat"] - df_bal["J_LW_nat"])
df_bal["LE_nat"]   = df_bal["LE_model"]
df_bal["G_nat"]    = df_bal["G_model"]

# Residuals = balance (≈ C · dT/dt, energy needed for T variations).
# Sign convention: positive = into the surface / system.
# PVf+cell system: q_cond_PVf is internal and cancels out.
#   balance: Rn_PVf - q_cond_PVb - H_PVf - q_PV
# PVb (single surface): Rn_PVb + q_cond_PVb - H_PVb
# nat (single surface): Rn_nat - H_nat - LE_nat - G_nat
df_bal["res_PVf_cell"] = df_bal["Rn_PVf"] - df_bal["qcond_PVb"] - df_bal["H_PVf"] - df_bal["qPV"]
df_bal["res_PVb"]      = df_bal["Rn_PVb"] + df_bal["qcond_PVb"] - df_bal["H_PVb"]
df_bal["res_nat"]      = df_bal["Rn_nat"] - df_bal["H_nat"] - df_bal["LE_nat"] - df_bal["G_nat"]

t_bal = df_bal["timestep"]
t_bal_min = int(t_bal.min())
t_bal_max = int(t_bal.max())
if T_START is not None:
    t_bal_min = max(t_bal_min, T_START)
if T_END is not None:
    t_bal_max = min(t_bal_max, T_END)

STEPS_PER_DAY  = 48
DAYS_PER_PAGE  = 10
STEPS_PER_PAGE = STEPS_PER_DAY * DAYS_PER_PAGE


# ══════════════════════════════════════════════
# Per-panel configuration
# ══════════════════════════════════════════════
# Each row: (label, terms, residual_col, ylim, pad, yticks)
#   label        : text drawn above the panel, e.g. "a) PVf + cell"
#   terms        : list of (column, label, linestyle, color)
#   residual_col : name of the residual column to plot as dashed black line
#   ylim         : (ymin, ymax) or None for auto-scaling
#   pad          : spacing between y-axis label and tick numbers (in points)
#   yticks       : list/range of y-tick positions, or None for automatic
surface_rows = [
    ("a) PVf + cell", [
        ("Rn_PVf",    r"$R_{\mathrm{n}}$",     "-",  "tab:red"),
        ("qcond_PVb", r"$\dot{q}_\mathrm{cond}$",  "-.", "tab:green"),
        ("H_PVf",     r"$H$",    "--", "tab:blue"),
        ("qPV",       r"$\dot{q}_\mathrm{PV}$",   ":",  "tab:orange"),
    ], "res_PVf_cell", (-100,1100), 22, np.linspace(0,1000,5)),
    ("b) PVb", [
        ("Rn_PVb",    r"$R_{\mathrm{n}}$",    "-",  "tab:red"),
        ("qcond_PVb", r"$\dot{q}_\mathrm{cond}$", "-.", "tab:green"),
        ("H_PVb",     r"$H$",     "--", "tab:blue"),
    ], "res_PVb", (-140,400), 17, np.linspace(-100,300,5)),
    ("c) nat", [
        ("Rn_nat", r"$R_{\mathrm{n}}$", "-",  "tab:red"),
        ("G_nat",  r"$G$",  "-.", "tab:green"),
        ("H_nat",  r"$H$",  "--", "tab:blue"),
        ("LE_nat", r"$LE$", ":",  "tab:purple"),
    ], "res_nat", (-120,550), 38, np.linspace(0,450,4)),
]


with PdfPages("plot_energy_balance_surfaces_2025.pdf") as pdf:
    current = t_bal_min
    while current <= t_bal_max:
        end = min(current + STEPS_PER_PAGE - 1, t_bal_max)

        mask = (t_bal >= current) & (t_bal <= end)
        db = df_bal.loc[mask]

        if db.empty:
            current = end + 1
            continue

        tb = db["timestep"]
        d1 = (current - t_bal_min) // STEPS_PER_DAY + 1
        d2 = (end - t_bal_min) // STEPS_PER_DAY + 1

        fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(18, 12))
        #fig.suptitle(f"Surface energy balance — timesteps {current}–{end} (days {d1}–{d2})",
        #             fontsize=FS_TITLE, fontweight="bold")

        for ax, (label, terms, res_col, ylim, pad, yticks) in zip(axes, surface_rows):
            for col, term_label, ls, color in terms:
                ax.plot(tb, db[col], label=term_label, lw=LW_LINE, ls=ls, color=color)
            # Residual = energy needed for temperature variations (C·dT/dt)
            ax.plot(tb, db[res_col], label=r"$C \cdot \frac{\partial T}{\partial t}$", lw=LW_RESIDUAL,
                    ls="--", color="black")
            ax.axhline(0.0, color="k", lw=0.6, alpha=0.5)
            ax.set_ylabel(r"Flux (W m$^{-2}$)", labelpad=pad)
            if ylim is not None:
                ax.set_ylim(ylim)
            if yticks is not None:
                ax.set_yticks(yticks)
            # Panel label above the frame, left aligned with the y axis
            ax.text(0.0, PANEL_LABEL_Y, label, transform=ax.transAxes,
                    fontsize=FS_PANEL, fontweight="normal",
                    va="bottom", ha="left")
            ax.grid(True, alpha=0.3, lw=LW_GRID)
            ax.legend(loc="upper right", ncol=len(terms) + 1)

        # Build daily ticks: one per day at 01:00 (timestep 1, 49, 97, ...).
        first_day_tick = current + ((1 - current) % STEPS_PER_DAY)
        tick_positions = list(range(first_day_tick, end + 1, STEPS_PER_DAY))

        # Midpoint between two consecutive day boundaries = middle of the day.
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

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        # Applied after tight_layout, which would otherwise reset it:
        # leaves room above each panel for the a)/b)/c) labels.
        fig.subplots_adjust(hspace=PANEL_HSPACE)
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        current = end + 1

print("Saved: plot_energy_balance_surfaces_2025.pdf")
