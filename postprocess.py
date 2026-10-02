"""Post-process ProjectProgram output.

Reads the per-temperature output files written by ProjectProgram, computes
steady-state surface coverage and oxidation probability, and makes plots.

Oxidation probability is defined as
    epsilon = (carbon atoms removed per step) / (O atoms hitting the surface per step).

Usage:
    python postprocess.py --flux 50 [--dir .] [--steady-fraction 0.2] [--show-T 1200]

Requires numpy and matplotlib.
"""
import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np


def load_case(directory, T):
    def read(prefix):
        return np.loadtxt(os.path.join(directory, f"{prefix}{T}.txt"))

    return {
        "total": read("surfCovTot"),
        "O": read("surfCovO"),
        "CO": read("surfCovCO"),
        "carbon": read("carbonFlux"),
    }


def find_temperatures(directory):
    temps = []
    for path in glob.glob(os.path.join(directory, "surfCovTot*.txt")):
        match = re.search(r"surfCovTot(\d+)\.txt$", path)
        if match:
            temps.append(int(match.group(1)))
    return sorted(temps)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--flux", type=float, required=True,
                        help="O atoms hitting the surface per time step (particle_flux input)")
    parser.add_argument("--dir", default=".", help="directory with the output files")
    parser.add_argument("--steady-fraction", type=float, default=0.2,
                        help="fraction of the final time steps averaged for steady state")
    parser.add_argument("--dt", type=float, default=1e-6, help="time step size (s)")
    parser.add_argument("--show-T", type=int, default=None,
                        help="temperature to plot as a coverage-vs-time history")
    parser.add_argument("--out", default="figures", help="output folder for plots")
    args = parser.parse_args()

    temps = find_temperatures(args.dir)
    if not temps:
        raise SystemExit(f"No surfCovTot*.txt files found in {args.dir}")
    os.makedirs(args.out, exist_ok=True)

    rows = []
    for T in temps:
        case = load_case(args.dir, T)
        n = len(case["total"])
        start = int(n * (1.0 - args.steady_fraction))
        # carbonFlux[0] is the initial placeholder; skip it for the average
        carbon = case["carbon"][max(start, 1):]
        rows.append((T,
                     case["total"][start:].mean(),
                     case["O"][start:].mean(),
                     case["CO"][start:].mean(),
                     carbon.mean() / args.flux))
    data = np.array(rows)

    header = "T_K  theta_total  theta_O  theta_CO  oxidation_probability"
    np.savetxt(os.path.join(args.out, "steady_state.txt"), data, header=header, fmt="%.6g")
    print(header)
    for row in data:
        print("  ".join(f"{v:.4g}" for v in row))

    # Steady-state coverage and oxidation probability vs temperature
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(data[:, 0], data[:, 1], "s", label="Total")
    ax1.plot(data[:, 0], data[:, 2], "o", label="O")
    ax1.plot(data[:, 0], data[:, 3], "^", label="CO")
    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel(r"Steady-state coverage $\theta$")
    ax1.legend()
    ax2.plot(data[:, 0], data[:, 4], "s")
    ax2.set_xlabel("Temperature (K)")
    ax2.set_ylabel(r"Oxidation probability $\epsilon$")
    fig.suptitle(f"{args.flux:g} O atoms per time step")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "steady_state.png"), dpi=200)

    # Coverage history for one temperature
    T_show = args.show_T if args.show_T in temps else temps[len(temps) // 2]
    case = load_case(args.dir, T_show)
    t_us = np.arange(len(case["total"])) * args.dt * 1e6
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(t_us, case["total"], label="Total")
    ax.plot(t_us, case["O"], label="O")
    ax.plot(t_us, case["CO"], label="CO")
    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel(r"Coverage $\theta$")
    ax.set_title(f"T = {T_show} K")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, f"coverage_history_T{T_show}.png"), dpi=200)

    # Langmuir verification plot, if the test output is present
    for candidate in ("langmuir_verification.txt",
                      os.path.join("build", "langmuir_verification.txt")):
        if os.path.exists(candidate):
            v = np.loadtxt(candidate)
            fig, ax = plt.subplots(figsize=(5, 4))
            ax.plot(v[:, 0], v[:, 1], ".", ms=2, label="Simulation")
            ax.plot(v[:, 0], v[:, 2], "-", label="Analytical")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel(r"Coverage $\theta$")
            ax.set_title(r"Langmuir verification ($r_A = 1.6$, $r_D = 2$)")
            ax.legend()
            fig.tight_layout()
            fig.savefig(os.path.join(args.out, "langmuir_verification.png"), dpi=200)
            break

    print(f"Plots written to {args.out}/")


if __name__ == "__main__":
    main()
