#!/usr/bin/env python3
"""Visualize FEM Laplace solver output (output/result.dat)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULT_FILE = ROOT / "output" / "result.dat"
FIGURE_FILE = ROOT / "output" / "result.png"


def load_result(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    x = data[:, 1]
    y = data[:, 2]
    phi = data[:, 3]
    return x, y, phi


def plot_solution(x: np.ndarray, y: np.ndarray, phi: np.ndarray, save_path: Path) -> None:
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    nx = len(x_unique)
    ny = len(y_unique)

    X = x.reshape(ny, nx)
    Y = y.reshape(ny, nx)
    Z = phi.reshape(ny, nx)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    ax0 = axes[0]
    levels = np.linspace(Z.min(), Z.max(), 21)
    contour = ax0.contourf(X, Y, Z, levels=levels, cmap="viridis")
    ax0.contour(X, Y, Z, levels=levels, colors="k", linewidths=0.3, alpha=0.35)
    fig.colorbar(contour, ax=ax0, label=r"$\phi$")
    ax0.set_xlabel("$x$")
    ax0.set_ylabel("$y$")
    ax0.set_title("FEM solution (contour)")
    ax0.set_aspect("equal")

    ax1 = axes[1]
    surface = ax1.pcolormesh(X, Y, Z, shading="gouraud", cmap="viridis")
    fig.colorbar(surface, ax=ax1, label=r"$\phi$")
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$y$")
    ax1.set_title("FEM solution (surface map)")
    ax1.set_aspect("equal")

    fig.suptitle("2D Laplace equation (finite element method)", fontsize=14)
    fig.savefig(save_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not RESULT_FILE.exists():
        raise FileNotFoundError(f"Result file not found: {RESULT_FILE}")

    x, y, phi = load_result(RESULT_FILE)
    FIGURE_FILE.parent.mkdir(parents=True, exist_ok=True)
    plot_solution(x, y, phi, FIGURE_FILE)
    print(f"Loaded {len(phi)} nodes from {RESULT_FILE}")
    print(f"phi range: [{phi.min():.6f}, {phi.max():.6f}]")
    print(f"Figure saved to {FIGURE_FILE}")


if __name__ == "__main__":
    main()
