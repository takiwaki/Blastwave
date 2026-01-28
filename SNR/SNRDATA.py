import os
import glob
import sys
import csv
from dataclasses import dataclass
from typing import List, Dict, Any, Tuple

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

# -----------------------------
# Matplotlib style (your settings)
# -----------------------------
fnameforfig = "sans-serif"
fsizeforfig = 14
fsizeforlabel = 16
plt.rcParams["font.family"] = fnameforfig
plt.rcParams["font.size"] = fsizeforfig
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True

cmapudc = cycler(color=["#ff2800", "#0041ff", "#35a16B", "#faf500", "#66ccff",
                        "#ff99a0", "#ff9900", "#9a0079", "#663300"])
plt.rcParams["axes.prop_cycle"] = cmapudc
cmap = ["#ff2800", "#0041ff", "#35a16B", "#faf500", "#66ccff",
        "#ff99a0", "#ff9900", "#9a0079", "#663300"]


@dataclass
class SimCurve:
    t_year: np.ndarray
    r_pc: np.ndarray
    Msw: np.ndarray
    kTshock: np.ndarray
    Vshock: np.ndarray
    Lbol: np.ndarray
    label: str


def main():
    # --- observation points ---
    obs = read_obs_data("SNR.csv")
    obs = calculate_radius_pc(obs)

    sim_name_files = [
        [r"$M_{ejecta}=10\,M_\odot, E_{exp}=10^{51}\,erg$", "../HYD1D/analysis/t-prof.dat"]
    ]

    sim_list = read_simulation_curves(sim_name_files)

    plot_age_radius(obs, sim_list, outputfile="Age-Radius.png")


# -----------------------------
# Observation (NO pandas)
# -----------------------------
def _to_float(s: str) -> float:
    try:
        return float(s)
    except Exception:
        return float("nan")


def read_obs_data(path: str) -> Dict[str, Any]:
    """
    Read SNR.csv without pandas.

    Returns a dict of numpy arrays:
      - "Name" (list[str])
      - "age [kyr]" (np.ndarray)
      - "size x[arcmin]" (np.ndarray)
      - "size y[arcmin]" (np.ndarray)
      - "Distance [kpc]" (np.ndarray)
      - plus computed "radius [pc]" after calculate_radius_pc()
    """
    try:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise RuntimeError("CSV header not found")

            names: List[str] = []
            age_kyr: List[float] = []
            size_x: List[float] = []
            size_y: List[float] = []
            dist_kpc: List[float] = []

            for row in reader:
                names.append((row.get("Name") or "").strip())
                age_kyr.append(_to_float(row.get("age [kyr]", "")))
                size_x.append(_to_float(row.get("size x[arcmin]", "")))
                size_y.append(_to_float(row.get("size y[arcmin]", "")))
                dist_kpc.append(_to_float(row.get("Distance [kpc]", "")))

    except OSError:
        print("cannot open " + path)
        sys.exit(1)

    return {
        "Name": names,
        "age [kyr]": np.asarray(age_kyr, dtype=float),
        "size x[arcmin]": np.asarray(size_x, dtype=float),
        "size y[arcmin]": np.asarray(size_y, dtype=float),
        "Distance [kpc]": np.asarray(dist_kpc, dtype=float),
    }


def calculate_radius_pc(obs: Dict[str, Any]) -> Dict[str, Any]:
    # arcmin -> rad
    arcmin = 0.000290888
    dtheta = (obs["size x[arcmin]"] + obs["size y[arcmin]"]) / 4.0 * arcmin
    obs["radius [pc]"] = dtheta * obs["Distance [kpc]"] * 1000.0
    return obs


# -----------------------------
# Simulation curves (NO pandas)
# -----------------------------
def read_one_sim_file(path: str, modelname: str) -> SimCurve:
    try:
        # Expect 2 columns: t_year, r_pc
        data = np.loadtxt(path, comments="#", dtype=float)
        if data.ndim == 1:
            # single row -> shape (2,) -> (1,2)
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            raise RuntimeError(f"Expected >=2 columns but got {data.shape[1]} columns")
        t = data[:, 0]
        r = data[:, 1]
        M = data[:, 2]
        T = data[:, 3]
        V = data[:, 4]
        L = data[:, 5]
    except Exception as e:
        raise RuntimeError(f"Failed to read simulation file: {path}\n{e}")

    return SimCurve(t_year=t, r_pc=r, Msw=M, kTshock=T, Vshock=V,Lbol=L,label=modelname)


def read_simulation_curves(name_file_list: List[List[str]]) -> List[SimCurve]:
    """
    Read many simulation files. Returns list[SimCurve].
    Each curve has: t_year, r_pc, label
    """
    curves: List[SimCurve] = []
    for name, f in name_file_list:
        try:
            curves.append(read_one_sim_file(f, name))
        except RuntimeError as e:
            print(e)
    return curves


# -----------------------------
# Plot
# -----------------------------
def plot_age_radius(obs: Dict[str, Any], sim_list: List[SimCurve], outputfile: str = "Age-Radius.png"):
    # Observation
    x_obs = obs["age [kyr]"] * 1000.0  # year
    y_obs = obs["radius [pc]"]
    labels = obs["Name"]

    fig = plt.figure(figsize=(6.4, 5.2))
    ax = fig.add_subplot(1, 1, 1)

    # Observed points
    ax.scatter(x_obs, y_obs, marker="o", color=cmap[1], label="Observed SNR")

    # Labels for observed points (optional)
    try:
        from adjustText import adjust_text
        texts = [plt.text(float(x_obs[i]), float(y_obs[i]), labels[i],
                          ha="center", va="center") for i in range(len(labels))]
        adjust_text(texts)
    except Exception:
        # adjustText not installed -> skip
        pass

    # Simulation curves
    for curve in sim_list:
        order = np.argsort(curve.t_year)  # sort by time in case the file is unsorted
        ax.plot(curve.t_year[order], curve.r_pc[order], lw=2, label=curve.label)

    ax.grid(color="lightgray")
    ax.set_xlabel(r"Age [year]", fontsize=fsizeforlabel)
    ax.set_ylabel(r"Radius [pc]", fontsize=fsizeforlabel)
    ax.legend(fontsize=10, frameon=False)
    fig.tight_layout()
    fig.savefig(outputfile)
    print("saved:", outputfile)


if __name__ == "__main__":
    main()
