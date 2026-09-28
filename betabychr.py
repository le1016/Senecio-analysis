#!/usr/bin/env python3
# ============================================================
# 全基因组 beta 值散点图 —— 染色体纵向堆叠 (每条一行) 版
# 纯 matplotlib (无 numba/datashader 依赖), 适合 500 万+ 位点
# 同时输出 PNG 和 可编辑 SVG
# 依赖: pip install pandas numpy matplotlib
# ============================================================

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                      # 服务器无显示环境必备
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# 关键: 让 SVG 里的文字保持为可编辑文本(而非转成路径),
#       这样在 Illustrator / Inkscape 里能直接改字号、字体、内容
plt.rcParams["svg.fonttype"] = "none"

# ------------------------------------------------------------
# ★ 配置区
# ------------------------------------------------------------
INFILE    = "Chrall_bt.txt"      # Tab 分隔, 无表头, 4 列: CHR START END BETA
THRESHOLD = 33.9388                 # ← 阈值线; 多条写成 [4.0, 5.0]
OUT_PNG   = "beta_by_chr.png"
OUT_SVG   = "beta_by_chr.svg"    # ← 可编辑矢量图

# ------------------------------------------------------------
# 1. 读取数据
# ------------------------------------------------------------
df = pd.read_csv(INFILE, sep="\t", header=None,
                 names=["CHR", "START", "END", "BETA"])
df["POS"] = df["START"]
df["CHRNUM"] = df["CHR"].str.extract(r"(\d+)").astype(int)
print(f"位点总数: {len(df):,}")

thr_list = THRESHOLD if isinstance(THRESHOLD, (list, tuple)) else [THRESHOLD]
top_thr  = max(thr_list)

chroms = sorted(df["CHRNUM"].unique())     # 自动识别有几条染色体 (你的是 19)
n      = len(chroms)
xmax   = df["POS"].max() / 1e6             # 全基因组统一 x 刻度(Mb), 各行对齐
ymin, ymax = df["BETA"].min(), df["BETA"].max()
print(f"染色体数: {n}")

# ------------------------------------------------------------
# 2. 纵向堆叠: n 行 1 列, 共用 x 轴(同一 bp 刻度, 便于跨染色体对齐)
# ------------------------------------------------------------
fig, axes = plt.subplots(n, 1, figsize=(12, 0.55 * n + 1),
                         sharex=True, sharey=True)
if n == 1:
    axes = [axes]

for ax, c in zip(axes, chroms):
    sub  = df[df["CHRNUM"] == c]
    pos  = sub["POS"].to_numpy() / 1e6     # 转 Mb
    beta = sub["BETA"].to_numpy()
    sig  = beta >= top_thr

    # 背景点(蓝) + 超阈值点(红); rasterized=True 让密集点云栅格化, SVG 不会爆炸
    ax.scatter(pos[~sig], beta[~sig], s=4, c="#4C78A8", marker=".",
               linewidths=0, rasterized=True)
    ax.scatter(pos[sig],  beta[sig],  s=6, c="#d62728", marker=".",
               linewidths=0, rasterized=True)

    # ★ 阈值线 (每行都画)
    for t in thr_list:
        ax.axhline(t, color="red", linestyle="--", linewidth=0.6)

    # 左侧标注染色体名, 右侧/上下边框去掉, y 轴只留 2-3 个刻度
    ax.set_ylabel(sub["CHR"].iloc[0], rotation=0, ha="right", va="center",
                  fontsize=9)
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.tick_params(axis="y", labelsize=7)
    ax.set_xlim(0, xmax)
    ax.set_ylim(ymin, ymax)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

axes[-1].set_xlabel("Position (Mb)")
fig.suptitle(f"Genome-wide SNP beta by chromosome ({len(df):,} loci)",
             y=0.995, fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.99])

# ------------------------------------------------------------
# 3. 输出: PNG (预览) + SVG (可编辑矢量)
# ------------------------------------------------------------
fig.savefig(OUT_PNG, dpi=300)
fig.savefig(OUT_SVG, dpi=300)             # ← 可编辑 SVG; 框架/文字矢量, 点云为内嵌栅格
print(f"已保存: {OUT_PNG}  和  {OUT_SVG}")
