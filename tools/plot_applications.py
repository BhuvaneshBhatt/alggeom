"""Regenerate the application guide's real contour and cusp illustration."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def build_figure():
    """Show the half-response contour and the effect of changing a real parameter range."""
    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "alggeom-applications"})
    figure, (response, cusp) = plt.subplots(1, 2, figsize=(10, 4.8), layout="constrained")
    blue, orange = "#176a99", "#b55222"

    angles = np.linspace(0, 2 * np.pi, 401)
    response.fill(np.cos(angles), np.sin(angles), color=blue, alpha=0.08)
    response.plot(np.cos(angles), np.sin(angles), color=blue, linewidth=2.5)
    response.scatter([0], [0], color=blue, s=18)
    response.text(0, 0.2, "$R<1/2$", ha="center", color=blue)
    response.text(0.73, -1.24, "$R=1/2$", color=blue)
    response.set(
        title="A two-input response reduces to one radius",
        xlabel="$x$",
        ylabel="$y$",
        xlim=(-1.5, 1.5),
        ylim=(-1.5, 1.5),
        aspect="equal",
    )

    parameter = np.linspace(0, 1.1, 300)
    cusp.plot(
        parameter**2,
        parameter**3,
        color=blue,
        linewidth=2.5,
        label=r"$s\geq 0$: original real trace",
    )
    cusp.plot(
        parameter**2,
        -(parameter**3),
        color=orange,
        linestyle="--",
        linewidth=2.5,
        label="$s<0$: additional real branch",
    )
    cusp.scatter([0, 1], [0, 1], color=blue, s=22, zorder=3)
    cusp.annotate(
        r"$t=\pm 1$",
        xy=(1, 1),
        xytext=(0.55, 1.22),
        arrowprops={"arrowstyle": "->", "color": blue},
        color=blue,
    )
    cusp.set(
        title="The cusp parameter must retain its range",
        xlabel="$u=s^2$",
        ylabel="$v=s^3$",
        xlim=(-0.12, 1.35),
        ylim=(-1.45, 1.45),
        aspect="equal",
    )
    handles, labels = cusp.get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside lower center", ncols=2, fontsize=9, frameon=False)
    for axis in (response, cusp):
        axis.grid(color="#d6dee3", linewidth=0.6, alpha=0.7)
        axis.set_axisbelow(True)
        axis.spines[["top", "right"]].set_visible(False)
    return figure


def main():
    """Write the SVG used in the guide and an optional PNG for visual inspection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", type=Path, help="optional PNG output path")
    options = parser.parse_args()
    destination = Path(__file__).resolve().parents[1] / "docs/images/application_curves.svg"
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure = build_figure()
    figure.savefig(destination, metadata={"Date": None, "Creator": "alggeom"})
    if options.preview:
        figure.savefig(options.preview, dpi=160)
    plt.close(figure)
    print(destination)


if __name__ == "__main__":
    main()
