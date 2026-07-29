#!/home/bater/miniconda3/envs/pytwm-env/bin/python
"""
plot_PM_self_force.py

Plots PM self-force components and magnitude from PM_self_force_diag.csv.

Usage examples
--------------
# 1D: fix dx=0, dy=0, vary dz
python plot_PM_self_force.py --fix dx=0.0 dy=0.0 --vary dz

# 1D: fix dx=0.2, dz=0.4, vary dy
python plot_PM_self_force.py --fix dx=0.2 dz=0.4 --vary dy

# 3D: fix dx=0, vary dy and dz, plot |F|
python plot_PM_self_force.py --fix dx=0.0 --vary dy dz --plot3d Fmod

# 3D: all four force components with fix dz=0
python plot_PM_self_force.py --fix dz=0.0 --vary dx dy
"""

import argparse
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

AXES  = ["dx", "dy", "dz"]
FORCE = ["Fx", "Fy", "Fz", "Fmod"]
LABELS = {
    "Fx":   r"$F_x$",
    "Fy":   r"$F_y$",
    "Fz":   r"$F_z$",
    "Fmod": r"$|\mathbf{F}|$",
}
COL_IDX = {name: i for i, name in enumerate(AXES + FORCE)}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot PM self-force from PM_self_force_diag.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--csv",    default="PM_self_force_diag.csv")
    p.add_argument("--fix",    nargs="+", default=[], metavar="AXIS=VALUE",
                   help="Fixed axes, e.g. --fix dx=0.0 dy=0.2")
    p.add_argument("--vary",   nargs="+", required=True, metavar="AXIS",
                   help="Varying axis/axes: one name for 1D, two for 3D")
    p.add_argument("--plot3d", default=None, choices=FORCE,
                   help="Single column for 3D surface (omit to plot all four)")
    p.add_argument("--tol",    type=float, default=1e-6)
    p.add_argument("--no-show", action="store_true",
                   help="Save PNG without opening a window")
    return p.parse_args()


# ---------------------------------------------------------------------------
# CSV loader (numpy only)
# ---------------------------------------------------------------------------

def load_csv(path):
    with open(path) as fh:
        header = fh.readline().strip().split(";")
    data = np.loadtxt(path, delimiter=";", skiprows=1)
    return header, data


def col(data, name):
    return data[:, COL_IDX[name]]


# ---------------------------------------------------------------------------
# Filter rows where fixed axes match within tol
# ---------------------------------------------------------------------------

def filter_data(data, fixed, tol):
    mask = np.ones(len(data), dtype=bool)
    for ax, val in fixed.items():
        mask &= np.abs(col(data, ax) - val) < tol
    subset = data[mask]
    if len(subset) == 0:
        avail = {ax: np.unique(col(data, ax)) for ax in AXES}
        sys.exit(
            f"No rows match {fixed} (tol={tol}).\n"
            + "\n".join(f"  {ax}: {v}" for ax, v in avail.items())
        )
    return subset


def parse_fix(fix_args):
    result = {}
    for token in fix_args:
        if "=" not in token:
            sys.exit(f"--fix expects axis=value, got: {token!r}")
        k, v = token.split("=", 1)
        if k not in AXES:
            sys.exit(f"Unknown axis {k!r}. Choose from {AXES}.")
        result[k] = float(v)
    return result


def fix_label(fixed):
    return "  ".join(f"{k}={v:.3g}" for k, v in fixed.items())


def save_and_show(fig, path, args):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Saved: {path}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# 1D plot: four subplots (Fx, Fy, Fz, |F|)
# ---------------------------------------------------------------------------

def plot_1d(data, vary_ax, fixed, args):
    subset = filter_data(data, fixed, args.tol)
    order  = np.argsort(col(subset, vary_ax))
    subset = subset[order]
    x = col(subset, vary_ax)

    fig, axes = plt.subplots(4, 1, figsize=(7, 10), sharex=True)
    fig.suptitle(
        f"PM self-force  —  {args.csv}\n"
        f"Fixed: {fix_label(fixed)}    Varying: {vary_ax}",
        fontsize=10,
    )
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for ax_plot, fc, color in zip(axes, FORCE, colors):
        ax_plot.plot(x, col(subset, fc), "o-", color=color, lw=1.5, ms=5)
        ax_plot.set_ylabel(LABELS[fc], fontsize=11)
        ax_plot.axhline(0, color="gray", lw=0.6, ls="--")
        ax_plot.grid(True, alpha=0.3)
    axes[-1].set_xlabel(rf"$\delta_{vary_ax[-1]}$", fontsize=12)
    plt.tight_layout()

    out = (
        f"PM_self_force_1D_{vary_ax}"
        f"_fix_{'_'.join(f'{k}{v:.2g}' for k, v in fixed.items())}.png"
    )
    save_and_show(fig, out, args)


# ---------------------------------------------------------------------------
# 3D surface plot for one force column
# ---------------------------------------------------------------------------

def plot_3d_one(data, vary_axes, fixed, fc, args):
    subset = filter_data(data, fixed, args.tol)
    ax1, ax2 = vary_axes

    u_vals = np.unique(col(subset, ax1))
    v_vals = np.unique(col(subset, ax2))
    U, V   = np.meshgrid(v_vals, u_vals)   # V=ax1 rows, U=ax2 cols
    Z      = np.full(U.shape, np.nan)

    for row in subset:
        i = np.searchsorted(u_vals, row[COL_IDX[ax1]])
        j = np.searchsorted(v_vals, row[COL_IDX[ax2]])
        Z[i, j] = row[COL_IDX[fc]]

    fig = plt.figure(figsize=(9, 6))
    ax3d = fig.add_subplot(111, projection="3d")
    surf = ax3d.plot_surface(U, V, Z, cmap=cm.viridis, alpha=0.85, edgecolor="none")
    ax3d.scatter(U.ravel(), V.ravel(), Z.ravel(), color="k", s=8, zorder=5)

    ax3d.set_xlabel(rf"$\delta_{ax2[-1]}$", fontsize=11)
    ax3d.set_ylabel(rf"$\delta_{ax1[-1]}$", fontsize=11)
    ax3d.set_zlabel(LABELS[fc], fontsize=11)
    ax3d.set_title(
        f"PM self-force: {LABELS[fc]}\n"
        f"Fixed: {fix_label(fixed)}    Varying: {ax1}, {ax2}",
        fontsize=10,
    )
    fig.colorbar(surf, ax=ax3d, shrink=0.5, pad=0.1, label=LABELS[fc])
    plt.tight_layout()

    out = (
        f"PM_self_force_3D_{fc}_{ax1}_{ax2}"
        f"_fix_{'_'.join(f'{k}{v:.2g}' for k, v in fixed.items())}.png"
    )
    save_and_show(fig, out, args)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    if args.no_show:
        matplotlib.use("Agg")

    header, data = load_csv(args.csv)
    if header != AXES + FORCE:
        sys.exit(f"Unexpected columns: {header}\nExpected: {AXES + FORCE}")

    fixed = parse_fix(args.fix)

    for v in args.vary:
        if v not in AXES:
            sys.exit(f"Unknown vary axis {v!r}. Choose from {AXES}.")
        if v in fixed:
            sys.exit(f"Axis {v!r} is in both --fix and --vary.")

    if len(args.vary) == 1:
        # 1D: four-panel plot
        plot_1d(data, args.vary[0], fixed, args)

    elif len(args.vary) == 2:
        # 3D: one column or all four
        cols_3d = [args.plot3d] if args.plot3d else FORCE
        for fc in cols_3d:
            plot_3d_one(data, args.vary, fixed, fc, args)

    else:
        sys.exit("Provide 1 or 2 --vary axes.")
