#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 SMC++ 结果画平滑的 Ne(t) 曲线（对数-对数坐标）。

两种输入（二选一）:
  1) --csv 3pop.marginal.csv
     推荐：读 smc++ plot --csv 导出的曲线（对常数模型和样条模型都适用）。
  2) model.final.json 文件列表
     仅适用于常数模型（JSON 中含 "epochs"）；样条模型请用 --csv。

用法（CSV 模式）:
  python3 plot_smoothed_ne.py --csv plots/3pop.marginal.csv \
      --gen-time 1 --out ne_smoothed.pdf --labels NJ XS Ssc \
      --glaciation --glaciation-range 19000 26500 --glaciation-label LGM \
      --split-table plots/3pop.split_time.tsv

说明:
  - 平滑仅影响展示，不改变推断结果；浅色线=原始曲线，深色线=log-log 平滑。
  - 需要 matplotlib；scipy 可选（没有则用滑动中值平滑）。
  - 冰期阴影带/分化时间线均为绝对年代（年）。
"""
import json
import os
import sys
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_marginal(path, gen_time, time_factor):
    """读取常数模型 JSON，返回 (t_years, Ne)。样条模型请用 load_csv。"""
    with open(path, "r", encoding="utf-8") as fh:
        m = json.load(fh)["model"]
    if "epochs" not in m:
        raise ValueError(
            "模型 JSON 中没有 'epochs' 键（可能是样条模型），可用键: %s。"
            "请改用 --csv 模式读取 smc++ 导出的 CSV。" % list(m.keys()))
    n0 = float(m["N0"])
    epochs = m["epochs"]  # [[t, ne], ...]：t 以 N0 世代为单位，ne 以 N0 为单位
    t = np.array([0.0] + [float(e[0]) for e in epochs])
    ne = np.array([1.0] + [float(e[1]) for e in epochs])
    return t * n0 * time_factor * gen_time, ne * n0


def load_csv(path):
    """读取 smc++ plot --csv 输出，返回 {标签: (t, Ne)}。
    smc++ 的 CSV 列名是 label,x,y,plot_type,plot_num：
      x = 时间(年/世代，随 -g 而定)，y = Ne，label = 群体/模型名。
    自动识别列名；识别不到时按前 3 列处理。"""
    import pandas as pd
    df = pd.read_csv(path)
    cols = {str(c).strip().lower(): c for c in df.columns}
    def pick(*names):
        for n in names:
            if n in cols:
                return cols[n]
        return None
    tcol = pick("time", "years", "year", "generations", "generation", "x")
    necol = pick("ne", "n", "effective_ne", "effective_population_size", "y")
    mcol = pick("model", "population", "pop", "label", "lineage")
    # 防止把 label 列误当时间列（smc++ 格式时 x 才是时间）
    if mcol is not None and tcol == mcol:
        tcol = pick("x") or df.columns[1]
        necol = pick("y") or df.columns[2]
    if tcol is None:
        tcol = df.columns[0]
    if necol is None:
        necol = df.columns[1]
    print("[csv] time列=%s  Ne列=%s  model列=%s" % (tcol, necol, mcol))
    # 只保留曲线路径行（去掉可能的其它 plot_type）
    if "plot_type" in cols:
        df = df[df[cols["plot_type"]].astype(str) == "path"]
    if mcol is None:
        return {"default": (df[tcol].to_numpy(float), df[necol].to_numpy(float))}
    out = {}
    for lab, grp in df.groupby(mcol):
        g = grp.sort_values(tcol)
        out[str(lab)] = (g[tcol].to_numpy(float), g[necol].to_numpy(float))
    return out


def smooth_line(t, ne, method="auto", window=7, s=None):
    """在 log-log 空间平滑曲线。返回 (t_smooth, ne_smooth)。"""
    mask = t > 0
    lt = np.log10(t[mask])
    ln = np.log10(ne[mask])
    order = np.argsort(lt)
    lt, ln = lt[order], ln[order]

    if method == "median":
        half = window // 2
        ln_s = np.array([np.median(ln[max(0, i - half):i + half + 1])
                         for i in range(len(ln))])
        return 10 ** lt, 10 ** ln_s

    try:
        from scipy.interpolate import UnivariateSpline
        if s is None:
            s = 0.5 * len(lt)
        spl = UnivariateSpline(lt, ln, k=3, s=s)
        tt = np.linspace(lt.min(), lt.max(), 300)
        return 10 ** tt, 10 ** spl(tt)
    except ImportError:
        return smooth_line(t, ne, method="median", window=window)


def main():
    ap = argparse.ArgumentParser(description="SMC++ Ne 曲线平滑绘图")
    ap.add_argument("models", nargs="*", help="model.final.json（常数模型，配合 --labels 用）")
    ap.add_argument("--csv", default=None,
                    help="smc++ plot --csv 导出的 CSV（推荐；样条模型必须用它）")
    ap.add_argument("--gen-time", type=float, default=1.0, help="年/代，默认 1")
    ap.add_argument("--time-factor", type=float, default=2.0,
                    help="JSON 模式的 N0 时间因子，默认 2；与 smc++ PDF 不符时改 1")
    ap.add_argument("--labels", nargs="+", default=None, help="每条曲线的标签")
    ap.add_argument("--split-table", default=None,
                    help="3pop.split_time.tsv：画出分化时间竖虚线")
    ap.add_argument("--glaciation", action="store_true",
                    help="标出末次冰期(LGM)阴影带（绝对年代，单位：年）")
    ap.add_argument("--glaciation-range", nargs=2, type=float,
                    default=[19000.0, 26500.0],
                    help="冰期时间范围(年)，默认 LGM 19–26.5 ka")
    ap.add_argument("--glaciation-label", default="LGM", help="冰期阴影带标签")
    ap.add_argument("--out", default="ne_smoothed.pdf")
    ap.add_argument("--smooth-method", default="auto", choices=["auto", "median"])
    ap.add_argument("--smooth-s", type=float, default=None,
                    help="scipy 样条平滑参数（越大越平滑）")
    ap.add_argument("--no-smooth", action="store_true",
                    help="不画显示层平滑线，只画原始曲线（模型已是样条时推荐）")
    ap.add_argument("--koala", action="store_true",
                    help="考拉风格：细线=每个群体，粗线=均值，阴影=CI/范围带")
    ap.add_argument("--ci-dir", default=None,
                    help="自助法 CI 目录（含 ci.<标签>.tsv，列: time mean lo hi）；"
                         "给出时每个群体画均值粗线+CI 阴影")
    ap.add_argument("--xlim", nargs=2, type=float, default=None)
    ap.add_argument("--ylim", nargs=2, type=float, default=None)
    args = ap.parse_args()

    if not args.csv and not args.models:
        sys.exit("需要 --csv 或至少一个 model.final.json")
    if args.csv:
        curves = load_csv(args.csv)
    else:
        labels = args.labels or [m.split("/")[-2] for m in args.models]
        if len(labels) != len(args.models):
            sys.exit("--labels 数量必须与模型文件数量一致")
        curves = {}
        for fn, lab in zip(args.models, labels):
            curves[lab] = load_marginal(fn, args.gen_time, args.time_factor)

    labels = args.labels or list(curves.keys())
    if len(labels) != len(curves):
        sys.exit("--labels 数量必须与曲线数量一致")

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.Set1(np.linspace(0, 1, len(curves)))

    for (lab, (t, ne)), c in zip(curves.items(), colors):
        if args.koala or args.ci_dir:
            # 考拉风格：细线 = 每个群体（原始曲线）
            if len(t) > 2 and np.all(t[:-1] < t[1:]):
                ax.step(t, ne, where="post", color=c, lw=1.2, alpha=0.65,
                        label="_nolegend_" if args.ci_dir else lab)
            else:
                ax.plot(t, ne, color=c, lw=1.2, alpha=0.65,
                        label="_nolegend_" if args.ci_dir else lab)
        elif args.no_smooth:
            # 只画原始曲线（模型为样条时 CSV 本身就是光滑曲线）
            if len(t) > 2 and np.all(t[:-1] < t[1:]):
                ax.step(t, ne, where="post", color=c, lw=2.2, label=lab)
            else:
                ax.plot(t, ne, color=c, lw=2.2, label=lab)
        else:
            # 原始浅色阶梯 + 深色平滑线
            if len(t) > 2 and np.all(t[:-1] < t[1:]):
                ax.step(t, ne, where="post", color=c, lw=1, alpha=0.35,
                        label="_nolegend_")
            ts, nes = smooth_line(t, ne, method=args.smooth_method, s=args.smooth_s)
            ax.plot(ts, nes, color=c, lw=2.5, label=lab)

    if args.ci_dir:
        # 自助法 CI：每个群体一条均值粗线 + CI 阴影（ci.<标签>.tsv: time mean lo hi）
        for lab, c in zip(curves.keys(), colors):
            p = os.path.join(args.ci_dir, "ci.%s.tsv" % lab)
            if os.path.exists(p):
                d = np.loadtxt(p, skiprows=1)
                if d.ndim == 1:
                    d = d.reshape(1, -1)
                ax.fill_between(d[:, 0], d[:, 2], d[:, 3], color=c, alpha=0.25,
                                linewidth=0)
                ax.plot(d[:, 0], d[:, 1], color=c, lw=2.5, label="%s mean" % lab)
    elif args.koala:
        # 快速版（无自助）：三群体 log-Ne 均值粗线 + min-max 范围带
        t_ok = [t for t, _ in curves.values()]
        lo_t = min(np.min(t[t > 0]) for t in t_ok if np.any(t > 0))
        hi_t = max(np.max(t) for t in t_ok)
        grid = np.logspace(np.log10(lo_t), np.log10(hi_t), 300)
        mats = []
        for t, ne in curves.values():
            ok = t > 0
            mats.append(np.interp(np.log10(grid), np.log10(t[ok]), np.log10(ne[ok])))
        L = np.vstack(mats)
        ax.fill_between(grid, 10 ** L.min(axis=0), 10 ** L.max(axis=0),
                        color="grey", alpha=0.20, linewidth=0)
        ax.plot(grid, 10 ** L.mean(axis=0), color="black", lw=3.0, label="Mean")

    if args.split_table:
        with open(args.split_table) as fh:
            header = fh.readline()
            cols = header.rstrip("\n").split("\t")
            yi = cols.index("split_years")
            p1 = cols.index("population_1")
            p2 = cols.index("population_2")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                y = float(row[yi])
                ax.axvline(y, color="grey", ls="--", lw=0.8, alpha=0.7)
                ax.text(y, ax.get_ylim()[0] * 1.1, f"{row[p1]}-{row[p2]}",
                        rotation=90, fontsize=7, color="grey")

    if args.glaciation:
        gs, ge = args.glaciation_range
        ax.axvspan(gs, ge, color="steelblue", alpha=0.15, zorder=0)
        ax.text(gs, ax.get_ylim()[1] * 0.96, args.glaciation_label,
                fontsize=10, color="steelblue", ha="left", va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Years ago" if args.gen_time != 1 else "Generations ago")
    ax.set_ylabel("Effective population size Ne")
    if args.ci_dir:
        title = "SMC++ marginal Ne histories (bootstrap 95% CI)"
    elif args.koala or args.no_smooth:
        title = "SMC++ marginal Ne histories"
    else:
        title = "SMC++ marginal Ne histories (smoothed line = display smoothing)"
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    if args.xlim:
        ax.set_xlim(*args.xlim)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print("saved:", args.out)


if __name__ == "__main__":
    main()
