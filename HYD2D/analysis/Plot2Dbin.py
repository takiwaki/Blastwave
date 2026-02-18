#!/usr/bin/python
import sys
import glob
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class FluidSnapshot:
    time: float
    rad: np.ndarray
    the: np.ndarray
    rho: np.ndarray
    pre: np.ndarray
    vr: np.ndarray
    vt: np.ndarray

###################
# setup for matplotlib
###################
fsizeforfig = 14
fsizeforlabel = 16
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": fsizeforfig,
    "axes.labelsize": fsizeforlabel,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "lines.linewidth": 2.5,
    "axes.linewidth": 1.5,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
})

def EnsureDir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def GetTimeIndex(path="../bindata", prefix="unf", suffix=".dat"):
    pattern = os.path.join(path, f"{prefix}*.dat")
    files = glob.glob(pattern)

    timesteps = []
    for f in files:
        name = os.path.basename(f)
        m = re.match(rf"{prefix}(\d{{5}}){suffix}", name)
        if m:
            timesteps.append(int(m.group(1)))

    return sorted(timesteps)

def ReadData(step: int, path="../bindata") -> FluidSnapshot:
    metafile = path + "/unf%05d.dat" % (step)
    try:
        with open(metafile, "r") as inputf:
            line = inputf.readline()
            item = line.split()
            t = float(item[1])

            line = inputf.readline()
            item = line.split()
            nrad = int(item[1])
            nradgs = int(item[2])

            line = inputf.readline()
            item = line.split()
            nthe = int(item[1])
            nthegs = int(item[2])
    except IOError:
        print(" cannot open " + metafile)
        sys.exit()

    datafile = path + "/bin%05d.dat" % (step)
    print("Reading", datafile)
    try:
        fp = open(datafile, "rb")
    except IOError:
        print(" cannot open " + datafile)
        sys.exit()

    nr = nrad + 2 * nradgs
    nt = nthe + 2 * nthegs

    pc = 3.085677581e18  # parsec [cm]

    # radial grid
    nvar = 3
    grid = np.fromfile(fp, np.float64, nr * nvar).reshape(nvar, nr)
    # radb = grid[0, nradgs:-(nradgs + 1)] / pc  # unused
    rada = grid[1, nradgs:] / pc
    # dvr = grid[2, nradgs:-(nradgs + 1)]       # unused

    # theta grid
    grid = np.fromfile(fp, np.float64, nt * nvar).reshape(nvar, nt)
    # theb = grid[0, nthegs:-(nthegs + 1)]       # unused
    thea = grid[1, nthegs:]
    # dvth = grid[2, nthegs:-(nthegs + 1)]      # unused

    # hydro
    nvar = 6
    Qhyd = np.fromfile(fp, np.float64, nr * nt * nvar).reshape(nvar, nt, nr)
    rho = Qhyd[0, nthegs:-(nthegs), nradgs:-(nradgs)]
    v1  = Qhyd[1, nthegs:-(nthegs), nradgs:-(nradgs)]
    v2  = Qhyd[2, nthegs:-(nthegs), nradgs:-(nradgs)]
    # v3  = Qhyd[3, nthegs:-(nthegs), nradgs:-(nradgs)]  # unused
    pre = Qhyd[4, nthegs:-(nthegs), nradgs:-(nradgs)]
    # ei  = Qhyd[5, nthegs:-(nthegs), nradgs:-(nradgs)]  # unused

    fp.close()

    rho = rho / 1.66053e-24  # g/cm^3 => 1/cm^3
    rad2D, the2D = np.meshgrid(rada, thea)

    year = 365.0e0 * 24 * 60 * 60  # sec
    tyear = t / year

    return FluidSnapshot(
        time=tyear,
        rad=rad2D,
        the=the2D,
        rho=rho,
        pre=pre,
        vr=v1,
        vt=v2,
    )

def MakeAllOutputs(
    steps,
    bindata_path="../bindata",
    outputdir="figures/",
    moviesdir="movies/",
    fps=15,
    dpi=300,
    save_png=True,
    save_mp4=True,
):
    """
    I/O is performed only once: snapshots are loaded first, then reused.
    Per-frame vmin/vmax (clim) + colorbar update_normal is applied.
    Outputs:
      figures/dentwoXXXXX.png ... and dentwo.mp4/gif (same for pretwo/veltwo)
    """
    EnsureDir(outputdir)

    # ---- 1) Load all snapshots (ONLY I/O pass) ----
    snapshots = []
    for step in steps:
        snapshots.append(ReadData(step, path=bindata_path))

    # Optional: reduce memory (comment in if needed)
    # for s in snapshots:
    #     s.rho = s.rho.astype(np.float32, copy=False)
    #     s.pre = s.pre.astype(np.float32, copy=False)
    #     s.vr  = s.vr.astype(np.float32, copy=False)

    # rmax (can be max over snapshots; no extra I/O)
    rmax_global = max(float(s.rad.max()) for s in snapshots)

    from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

    def build_anim(field_name: str, title: str, out_stem: str):
        s0 = snapshots[0]
        rad = s0.rad
        the = s0.the

        fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
        ax.set_title(title)
        ax.set_thetalim(-np.pi, np.pi)
        ax.set_rlim(0, rmax_global)
        ax.set_theta_zero_location("N")
        ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
        ax.set_yticks([])

        # texts
        time_text = ax.text(0.0, 0.0, "", transform=ax.transAxes, ha="left")
        rad_text  = ax.text(1.0, 0.0, r"$R=$" + f"{rmax_global:.0f}" + " pc",
                            transform=ax.transAxes, ha="right")

        arr0 = getattr(s0, field_name)
        # initial per-frame clim
        vmin0 = float(np.nanmin(arr0))
        vmax0 = float(np.nanmax(arr0))
        if not np.isfinite(vmin0): vmin0 = 0.0
        if not np.isfinite(vmax0): vmax0 = 1.0
        if vmin0 == vmax0: vmax0 = vmin0 + 1e-12

        im1 = ax.pcolormesh(the,  rad, arr0, shading="auto", vmin=vmin0, vmax=vmax0)
        im2 = ax.pcolormesh(-the, rad, arr0, shading="auto", vmin=vmin0, vmax=vmax0)
        cbar = fig.colorbar(im1)

        fig.tight_layout()

        # QuadMesh update: set_array wants flattened values
        def update(i: int):
            s = snapshots[i]
            arr = getattr(s, field_name)

            vmin = float(np.nanmin(arr))
            vmax = float(np.nanmax(arr))
            if not np.isfinite(vmin): vmin = 0.0
            if not np.isfinite(vmax): vmax = 1.0
            if vmin == vmax: vmax = vmin + 1e-12

            im1.set_array(arr.ravel())
            im2.set_array(arr.ravel())
            im1.set_clim(vmin, vmax)
            im2.set_clim(vmin, vmax)
            cbar.update_normal(im1)  # update colorbar to new clim

            time_text.set_text(r"$T=$" + f"{s.time:.0f}" + " year")
            return (im1, im2, time_text, rad_text)

        # ---- PNG output using same update() ----
        if save_png:
            for i, step in enumerate(steps):
                update(i)
                png_path = os.path.join(outputdir, f"{out_stem}{step:05d}.png")
                print("Writing", png_path)
                fig.savefig(png_path, dpi=dpi)

        # ---- Animations ----
        ani = FuncAnimation(fig, update, frames=len(steps), interval=1000/fps, blit=False)

        if save_mp4:
            mp4_path = os.path.join(moviesdir, "ani"+out_stem + ".mp4")
            try:
                writer = FFMpegWriter(fps=fps, metadata={"artist": "matplotlib"})
                print("Writing", mp4_path)
                ani.save(mp4_path, writer=writer, dpi=dpi)
            except Exception as e:
                print("MP4 output failed (ffmpeg not available?):", e)

        plt.close(fig)

    build_anim("rho", r"$\rho [{\rm 1/cm^3}]$", "dentwo")
    build_anim("pre", r"$p [{\rm erg/cm^3}]$", "pretwo")
    build_anim("vr",  r"$v_r [{\rm cm/s}]$",   "veltwo")

def Main():
    steps = GetTimeIndex()
    if len(steps) == 0:
        print("No timesteps found.")
        return

    MakeAllOutputs(
        steps,
        bindata_path="../bindata",
        outputdir="figures/",
        moviesdir="movies/",
        fps=15,
        dpi=300,
        save_png=True,
        save_mp4=True,
    )

if __name__ == "__main__":
    Main()
