#!/usr/bin/python

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

dir_path = "./output/"
year= 60*60*24*365
def Main():
    global dir_path
    files=GetFileList()


    time = np.empty(0,dtype=float)
    R    = np.empty(0,dtype=float)
    M    = np.empty(0,dtype=float)
    T    = np.empty(0,dtype=float)
    L    = np.empty(0,dtype=float)
    Eexp = np.empty(0,dtype=float)
    for file in files:
        print(file)
        timeloc,Rloc,Mloc,Tloc,Lloc,Eexploc = ReadData(file)
#        print(str(timeloc)+str(Eexploc))
        time = np.append(time, timeloc)
        R    = np.append(   R,    Rloc)
        M    = np.append(   M,    Mloc)
        T    = np.append(   T,    Tloc)
        L    = np.append(   L,    Lloc)
        Eexp = np.append(Eexp, Eexploc)
    PlotTimeEexp(time,Eexp)

def GetFileList():
    global dir_path
    filenames= dir_path+"tpr*.dat"
    files = glob.glob(filenames)
    files = sorted(files)
    return files

def ReadData(file):
    input = file
    try:
        inputf= open(input, 'r')
        header= inputf.readline() # Read the fisrt line
        item= header.split()
        t = float(item[0])
        R = float(item[1])
        M = float(item[2])
        T = float(item[3])
        L = float(item[4])
        Eexp =float(item[5])
        print(t)
        inputf.close()

    except IOError:
        print(" cannot open " + input )
        sys.exit()
    
    return t, R,M,T,L,Eexp

def PlotTimeEexp(time,Eexp):
  from matplotlib import ticker, cm, colors
  #######################################
  # Space Time diagram
  #######################################
  outputdatapath="./"
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

  outputfile=outputdatapath+ "t-E.png"
  fig1 = plt.figure(figsize=(6,4.5))
  ax1 = fig1.add_subplot(1,1,1)
  time = time/year
  ax1.plot(time,Eexp,linewidth=2)
  ax1.set_xlabel(r"$t\ [{\rm year}]$", fontsize=fsizeforlabel)
  ax1.set_ylabel(r"$E_{\rm exp}\ [{\rm erg}]$", fontsize=fsizeforlabel)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)

Main()
