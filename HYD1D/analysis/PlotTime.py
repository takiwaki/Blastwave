#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path
from typing import Dict

import numpy as np
import matplotlib.pyplot as plt


def read_tprof_dat(path) -> Dict[str, np.ndarray]:
    """
    Read a time-series data file (t-prof.dat) without using pandas.
    """
    path = Path(path)
    print("reading ",path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    with path.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("Input file is empty.")

    header = lines[0].strip()
    if not header.startswith("#"):
        raise ValueError("The first line must start with '#' and contain column names.")

    columns = header.lstrip("#").strip().split()
    ncol = len(columns)

    values = []
    for lineno, line in enumerate(lines[1:], start=2):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) != ncol:
            raise ValueError(
                f"Line {lineno}: expected {ncol} columns, got {len(parts)}"
            )

        values.append([float(x) for x in parts])

    if not values:
        raise ValueError("No numerical data found in the file.")

    arr = np.array(values, dtype=float)

    time_cols = [i for i, c in enumerate(columns) if c.startswith("time")]
    if time_cols:
        it = time_cols[0]
        order = np.argsort(arr[:, it])
        arr = arr[order]

    data = {}
    for i, name in enumerate(columns):
        data[name] = arr[:, i]

    return data

def plot_time_vs_etot(data: Dict[str, np.ndarray], output_png: str = "t-E.png") -> None:
    """
    Plot total energy as a function of time.
    """
    time_col = "time[year]"
    energy_col = "Etot[erg]"

    if time_col not in data:
        raise KeyError(f"Column '{time_col}' not found in data.")
    if energy_col not in data:
        raise KeyError(f"Column '{energy_col}' not found in data.")

    t = data[time_col]
    e = data[energy_col]

    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(t, e, linewidth=2)

    ax.set_xlabel(r"$t\ [{\rm year}]$")
    ax.set_ylabel(r"$E_{\rm tot}\ [{\rm erg}]$")

    fig.tight_layout()
    fig.savefig(output_png, dpi=200)
    print(f"Figure saved to: {output_png}")


def main():
    input_file = Path("./t-prof.dat")

    data = read_tprof_dat(input_file)

    plot_time_vs_etot(data, output_png="./t-E.png")


if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        print(f"[ERROR] {err}", file=sys.stderr)
        sys.exit(1)
