#!/usr/bin/python

import sys
import glob
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import math
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
fsizeforfig=14
fsizeforlabel=16
plt.rcParams.update({
    # font
    "font.family": "sans-serif",
    "font.size": fsizeforfig,
    "axes.labelsize": fsizeforlabel,
    # ticks
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,

    # line / axis
    "lines.linewidth": 2.5,
    "axes.linewidth": 1.5,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
})

def Main():
    steps=GetTimeIndex()
    for step in steps:
        fdata = ReadData(step)
        PlotRadTheData(step,fdata)

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

def ReadData(step:int,path="../bindata")-> FluidSnapshot:
    metafile = path + "/unf%05d.dat"%(step)
    input = metafile
    try:
        inputf= open(input, 'r')
        line = inputf.readline() # Read the fisrt line
        item= line.split()
        t = float(item[1])
        line = inputf.readline() # Read the fisrt line
        item= line.split()
        nrad = int(item[1])
        nradgs = int(item[2])
        line = inputf.readline() # Read the fisrt line
        item= line.split()
        nthe = int(item[1])
        nthegs = int(item[2])
        #print("t:",t)
        #print("grid:",nrad+2*nradgs,nthe+2*nthegs)
        inputf.close()
    except IOError:
        print(" cannot open " + input )
        sys.exit()

    datafile = path + "/bin%05d.dat"%(step)
    input = datafile
    print("Reading ",input)
    try:
        fp = open(datafile,"rb")
    except IOError:
        print(" cannot open " + input )
        sys.exit()
    nr = nrad+2*nradgs
    nt = nthe+2*nthegs

    pc  = 3.085677581e18   # parsec [cm]
    nvar = 3
    grid = np.fromfile(fp,np.float64,nr*nvar).reshape(nvar,nr)
    radb = grid[0,nradgs:-(nradgs+1)]/pc
    rada = grid[1,nradgs:]/pc
    dvr  = grid[2,nradgs:-(nradgs+1)]
    
    grid = np.fromfile(fp,np.float64,nt*nvar).reshape(nvar,nt)
    theb = grid[0,nthegs:-(nthegs+1)]
    thea = grid[1,nthegs:]
    dvth = grid[2,nthegs:-(nthegs+1)]

    nvar = 6
    Qhyd = np.fromfile(fp,np.float64,nr*nt*nvar).reshape(nvar,nt,nr)
    rho = Qhyd[0,nthegs:-(nthegs),nradgs:-(nradgs)]
    v1   = Qhyd[1,nthegs:-(nthegs),nradgs:-(nradgs)]
    v2  = Qhyd[2,nthegs:-(nthegs),nradgs:-(nradgs)]
    v3  = Qhyd[3,nthegs:-(nthegs),nradgs:-(nradgs)]
    pre = Qhyd[4,nthegs:-(nthegs),nradgs:-(nradgs)]
    ei  = Qhyd[5,nthegs:-(nthegs),nradgs:-(nradgs)]
    rho = rho/ 1.66053e-24 # g/cm^3 => 1/cm^3
    rad2D, the2D = np.meshgrid(rada,thea)

    year = 365.0e0*24*60*60 # sec
    tyear= t/year
    return FluidSnapshot(
        time = tyear,
        rad  = rad2D,
        the  = the2D, 
        rho  = rho, 
        pre  = pre, 
        vr   = v1, 
        vt   = v2,)

def PlotRadTheData(num,fdata):
  from matplotlib import ticker, cm, colors
  time = fdata.time
  rad = fdata.rad
  the = fdata.the
  rho = fdata.rho
  pre = fdata.pre
  vel = fdata.vr
  
  timetxt=r"$T=$"+f"{time:.0f}"+" year"
  rmax = rad.max()
  radtxt=r"$R=$"+f"{rmax:.0f}"+" pc"
  
  outputdatapath = "figures/"
  outputfile=outputdatapath+ "dentwo%05d.png"%(num)
  fig1,ax = plt.subplots(subplot_kw={'projection': 'polar'})
  ax.set_title(r"$\rho [{\rm 1/cm^3}]$")
  ax.set_thetalim(-np.pi, np.pi)
  ax.set_rlim(0, rmax)
  ax.set_theta_zero_location("N")
  im1=ax.pcolormesh( the,rad,rho)
  im2=ax.pcolormesh(-the,rad,rho)
  ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
  fig1.colorbar(im1)
  ax.set_yticks([])
  ax.text(1.0,0.0,radtxt,transform=ax.transAxes,ha="right")
  ax.text(0.0, 0.0,timetxt, transform=ax.transAxes,ha="left")
  fig1.tight_layout()
  print("Writing "+outputfile)
  fig1.savefig(outputfile)
  plt.close()

  outputfile=outputdatapath+ "pretwo%05d.png"%(num)
  fig1,ax = plt.subplots(subplot_kw={'projection': 'polar'})
  ax.set_title(r"$p [{\rm erg/cm^3}]$")
  ax.set_thetalim(-np.pi, np.pi)
  ax.set_rlim(0, rad.max())
  ax.set_theta_zero_location("N")
  im1=ax.pcolormesh( the,rad,pre)
  im2=ax.pcolormesh(-the,rad,pre)
  ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
  fig1.colorbar(im1)
  ax.set_yticks([])
  ax.text(1.0,0.0,radtxt,transform=ax.transAxes,ha="right")
  ax.text(0.0, 0.0,timetxt, transform=ax.transAxes,ha="left")
  fig1.tight_layout()
  print("Writing "+outputfile)
  fig1.savefig(outputfile)
  plt.close()

  outputfile=outputdatapath+ "veltwo%05d.png"%(num)
  fig1,ax = plt.subplots(subplot_kw={'projection': 'polar'})
  ax.set_title(r"$v_r [{\rm cm/s}]$")
  ax.set_thetalim(-np.pi, np.pi)
  ax.set_rlim(0, rad.max())
  ax.set_theta_zero_location("N")
  im1=ax.pcolormesh( the,rad,vel)
  im2=ax.pcolormesh(-the,rad,vel)
  ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
  fig1.colorbar(im1)
  ax.set_yticks([])
  ax.text(1.0,0.0,radtxt,transform=ax.transAxes,ha="right")
  ax.text(0.0, 0.0,timetxt, transform=ax.transAxes,ha="left")
  fig1.tight_layout()
  print("Writing "+outputfile)
  fig1.savefig(outputfile)
  plt.close()

Main()
