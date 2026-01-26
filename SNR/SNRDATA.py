import os
import glob
import sys
import pandas as pd
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


def main():
    # --- observation points ---
    dfobs = read_obs_data("SNR.csv")
    dfobs = calculate_radius_pc(dfobs)

    sim_name_files = [ [r"$M_{ejecta}=10\,M_\odot, E_{exp}=10^{51}\,erg$","../HYD1D/analysis/t-prof.dat"]]
    
    dfsim_list = read_simulation_curves(sim_name_files)

    plot_age_radius(dfobs, dfsim_list, outputfile="Age-Radius.png")


# -----------------------------
# Observation
# -----------------------------
def read_obs_data(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=",", skiprows=0, header=0)
    except IOError:
        print("cannot open " + path)
        sys.exit(1)


def calculate_radius_pc(df: pd.DataFrame) -> pd.DataFrame:
    # arcmin -> rad
    arcmin = 0.000290888
    dtheta = (df["size x[arcmin]"] + df["size y[arcmin]"]) / 4.0 * arcmin
    df["radius [pc]"] = dtheta * df["Distance [kpc]"] * 1000.0
    return df


# -----------------------------
# Simulation curves
# -----------------------------
def read_one_sim_file(path: str, modelname) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, delim_whitespace=True, comment="#", header=None)
        t = df.iloc[:, 0].to_numpy()
        r = df.iloc[:, 1].to_numpy()

    except Exception as e:
        raise RuntimeError(f"Failed to read simulation file: {path}\n{e}")

    out = pd.DataFrame({"t_year": t, "r_pc": r})
    out["label"] = modelname
    return out


def read_simulation_curves(name_file_list):
    """
    Read many simulation files. Returns list[DataFrame].
    Each DF has: t_year, r_pc, label
    """
    dfs = []
    for name,f in name_file_list:
        try:
            dfs.append(read_one_sim_file(f,name))
        except RuntimeError as e:
            print(e)
    return dfs


# -----------------------------
# Plot
# -----------------------------
def plot_age_radius(dfobs: pd.DataFrame, dfsim_list, outputfile="Age-Radius.png"):
    # Observation
    x_obs = dfobs["age [kyr]"] * 1000.0  # year
    y_obs = dfobs["radius [pc]"]
    text = dfobs["Name"]

    fig = plt.figure(figsize=(6.4, 5.2), layout="tight")
    ax = fig.add_subplot(1, 1, 1)

    # Observed points
    ax.scatter(x_obs, y_obs, marker="o", color=cmap[1], label="Observed SNR")

    # Labels for observed points (optional)
    try:
        from adjustText import adjust_text
        texts = [plt.text(x_obs.iloc[i], y_obs.iloc[i], text.iloc[i],
                          ha="center", va="center") for i in range(len(x_obs))]
        adjust_text(texts)
    except Exception:
        # adjustText not installed -> skip
        pass

    # Simulation curves
    for i, df in enumerate(dfsim_list):
        # sort by time in case the file is unsorted
        df_sorted = df.sort_values("t_year")
        ax.plot(df_sorted["t_year"], df_sorted["r_pc"],
                lw=2, label=f"{df_sorted['label'].iloc[0]}")

    ax.grid(color="lightgray")
    ax.set_xlabel(r"Age [year]", fontsize=fsizeforlabel)
    ax.set_ylabel(r"Radius [pc]", fontsize=fsizeforlabel)
    ax.legend(fontsize=10, frameon=False)

    fig.savefig(outputfile)
    print("saved:", outputfile)

if __name__ == "__main__":
    main()
