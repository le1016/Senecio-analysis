# -*- coding: utf-8 -*-
"""
Whole-chromosome scatter plot of raw BetaScan beta values.

Reads a 4-column, tab-separated betascan output file
    chrom  start  end  beta
and draws every raw point of the chosen chromosome as a dot
(position vs beta), with a horizontal threshold line and
above-threshold points in red.

Usage:
    python plot_chr_scatter.py [data_file] [chrom] [threshold]

Defaults: data_file = Chrall_bt.txt, chrom = Chr07, threshold = 26.44518
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
data_file = sys.argv[1] if len(sys.argv) > 1 else "Chrall_bt.txt"
chrom     = sys.argv[2] if len(sys.argv) > 2 else "Chr07"
threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 26.44518

COL_ABOVE = "#d62728"    # red
COL_BELOW = "#D2D2D2"    # steel blue

# ----------------------------------------------------------------------
# Read only the requested chromosome
# ----------------------------------------------------------------------
prefix = chrom + "\t"
pos, beta = [], []
with open(data_file, "r") as fh:
    for line in fh:
        if not line.startswith(prefix):
            continue
        parts = line.rstrip("\n").split("\t")
        pos.append(int(parts[1]))
        beta.append(float(parts[3]))

if not pos:
    sys.exit(f"No points found for chromosome {chrom} in {data_file}.")

pos  = np.array(pos, dtype=float) / 1e6     # bp -> Mb
beta = np.array(beta, dtype=float)

above = beta > threshold
n_above = int(above.sum())
print(f"{chrom}: {len(beta)} raw points; above threshold {threshold}: {n_above}")
print(f"Beta range: {beta.min():.4f} - {beta.max():.4f}; "
      f"position range: {pos.min():.2f} - {pos.max():.2f} Mb")

# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(16, 4), dpi=300)

ax.scatter(pos[~above], beta[~above], s=3, c=COL_BELOW,
           alpha=0.3, edgecolors="none", label=f"β ≤ {threshold}")
ax.scatter(pos[above],  beta[above],  s=6, c=COL_ABOVE,
           alpha=0.9, edgecolors="none", label=f"β > {threshold}  (n = {n_above})")

ax.axhline(threshold, color="black", linestyle="--", linewidth=1.2,
           label=f"Threshold = {threshold}")

ax.set_xlabel(f"{chrom} position (Mb)")
ax.set_ylabel("β value")
ax.set_title(f"BetaScan β values — {chrom} scatter (all raw points, "
             f"threshold = {threshold})")

ax.set_xlim(pos.min(), pos.max())
ax.margins(y=0.05)
ax.grid(alpha=0.25, linestyle=":")
ax.legend(loc="upper right", frameon=True, fontsize=10)

fig.tight_layout()
thr_tag = str(threshold)
out_png = f"{chrom}_betascan_scatter_thr{thr_tag}.png"
out_pdf = out_png.replace(".png", ".pdf")
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print(f"Saved: {out_png}  (also {out_pdf})")
