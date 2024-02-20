#!/usr/bin/python

import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dir_path = "./output/"
year = 60*60*24*365
pc   = 3.085677581e18 
def Main():
    global dir_path
    files=GetFileList()

    is_initial = True
    for file in files:
        print(file)
        t,rad,rho = ReadData(file)
        t = t/year
        nrad=len(rad)
        time = np.array([t*np.ones(nrad)]).T # horizontarl => vertical array
        rad  = rad/pc
        if (is_initial):
            timemin = t
            Xtd = np.empty((nrad,0),dtype=float)
            Ytd = np.empty((nrad,0),dtype=float)
            Ztd = np.empty((nrad,0),dtype=float)
            is_initial = False
    # making three 2D array
        Xtd = np.append(Xtd, time, axis=1)
        Ytd = np.append(Ytd, rad,  axis=1)
        Ztd = np.append(Ztd, rho,  axis=1)
    timemax=t
    PlotSpaceTimeData(Xtd,Ytd,Ztd,timemin,timemax)

def GetFileList():
    global dir_path
    filenames= dir_path+"onepro*.dat"
    files = glob.glob(filenames)
    files = sorted(files)
    return files

def ReadData(file):
    input = file
    try:
        results = np.genfromtxt(input,skip_header=2,delimiter='    ') # Read numbers
#        print(results)
        inputf= open(input, 'r')
        header= inputf.readline() # Read the fisrt line
        item= header.split()
        t = float(item[2])
        tyr = t/year
        print("t="+str(tyr)+" year")
        inputf.close()

    except IOError:
        print(" cannot open " + input )
        sys.exit()
    rad, rho, pre, vel =  np.split(results,4,1)
    
    return t, rad, rho

def PlotSpaceTimeData(Xtd,Ytd,Ztd,timemin,timemax):
  from matplotlib import ticker, cm, colors
  #######################################
  # Space Time diagram
  #######################################

  fnameforfig="sans-serif"
  fsizeforfig=14
  fsizeforlabel=16

  plt.rcParams['font.family'] = fnameforfig
  plt.rcParams['font.size'] = fsizeforfig
  plt.rcParams['xtick.direction'] = 'in'
  plt.rcParams['ytick.direction'] = 'in'
  plt.rcParams["xtick.minor.visible"] = True
  plt.rcParams["ytick.minor.visible"] = True
  plt.rcParams['xtick.top'] = True
  plt.rcParams['ytick.right'] = True

  figB = plt.figure(figsize=(6,6))
  ax4 = figB.add_subplot(1,1,1)
  cDen = make_colormap({0.:'#ff0000',1.0:'#00ff00',2:'#0000ff',3:'#000000'})
  tpc=ax4.pcolormesh(Xtd, Ytd,Ztd,cmap=cDen,norm=colors.LogNorm(vmin=Ztd.min(), vmax=Ztd.max()))
  ax4.set_title(r"$\rho\ [{\rm 1/cm^2}]$", fontsize=fsizeforlabel)
  ax4.set_ylabel(r"$r\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax4.set_xlabel(r"$t\ [{\rm year}]$", fontsize=fsizeforlabel)
  figB.colorbar(tpc)
  ax4.set_xlim(timemin,timemax)
  ax4.set_ylim(0,9)
#  anat = np.linspace(0.0, 0.4, 100)
#  anar = anat**(2/5)
#  ax4.plot(anat,anar, label="analytic" ,color='white', linewidth=4)
  figB.tight_layout()

# Save the figure.
  figB.savefig('t-r-rho.png')

def make_colormap(colors):
    from matplotlib.colors import LinearSegmentedColormap, ColorConverter
    from numpy import sort

    z  = np.array(sorted(colors.keys()))
    n  = len(z)
    z1 = min(z)
    zn = max(z)
    x0 = (z - z1) / (zn - z1)

    CC = ColorConverter()
    R = []
    G = []
    B = []
    for i in range(n):
        Ci = colors[z[i]]
        if type(Ci) == str:
            RGB = CC.to_rgb(Ci)
        else:
            RGB = Ci
        R.append(RGB[0])
        G.append(RGB[1])
        B.append(RGB[2])

    cmap_dict = {}
    cmap_dict['red']   = [(x0[i],R[i],R[i]) for i in range(len(R))]
    cmap_dict['green'] = [(x0[i],G[i],G[i]) for i in range(len(G))]
    cmap_dict['blue']  = [(x0[i],B[i],B[i]) for i in range(len(B))]
    mymap = LinearSegmentedColormap('mymap',cmap_dict)
    return mymap

Main()
