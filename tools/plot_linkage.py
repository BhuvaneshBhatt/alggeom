"""Draw the independently derived two-link workspace and its two regular elbow positions."""

import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def main():
    """Write a vector illustration of the reachable disk and the fiber over (1,0)."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", type=Path)
    args = parser.parse_args()
    matplotlib.use("Agg")
    theta = np.linspace(0, 2 * np.pi, 401)
    figure, (workspace, fiber) = plt.subplots(1, 2, figsize=(9, 4), layout="constrained")
    workspace.fill(2 * np.cos(theta), 2 * np.sin(theta), color="#d8e8ff")
    workspace.plot(2 * np.cos(theta), 2 * np.sin(theta), color="#2459a6")
    workspace.scatter([0, 1, 3], [0, 0, 0], c=["#2459a6", "#18804b", "#ac3434"])
    workspace.annotate("fully folded", (0, 0), xytext=(-1.7, -0.8), arrowprops={"arrowstyle": "->"})
    workspace.annotate("unreachable", (3, 0), xytext=(1.4, 1.4), arrowprops={"arrowstyle": "->"})
    workspace.set(
        xlim=(-2.4, 3.5),
        ylim=(-2.4, 2.4),
        title=r"Reachable endpoints: $x^2+y^2\leq4$",
        xlabel="x",
        ylabel="y",
    )
    h = np.sqrt(3) / 2
    for sign, color in ((1, "#18804b"), (-1, "#7a4bb3")):
        fiber.plot(
            [0, 0.5, 1],
            [0, sign * h, 0],
            "-o",
            color=color,
            label="upper elbow" if sign == 1 else "lower elbow",
        )
    fiber.annotate(r"$(1,0)$", (1, 0), xytext=(0.95, 0.16))
    fiber.set(
        xlim=(-0.25, 1.4),
        ylim=(-1.2, 1.2),
        title="Two unit-link configurations",
        xlabel="x",
        ylabel="y",
    )
    fiber.legend(loc="upper right", fontsize=8)
    for ax in (workspace, fiber):
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    destination = Path(__file__).resolve().parents[1] / "docs/images/linkage-workspace.svg"
    destination.parent.mkdir(exist_ok=True)
    figure.savefig(destination, metadata={"Date": None})
    if args.preview:
        figure.savefig(args.preview, dpi=120)
    plt.close(figure)
    print(destination)


if __name__ == "__main__":
    main()
