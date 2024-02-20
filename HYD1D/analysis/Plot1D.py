#!/usr/bin/python

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

dir_path = "./output/"
outputdatapath="./figures/"

year = 60*60*24*365
pc   = 3.085677581e18

def Main():
    import os
    global dir_path
    files=GetFileList()
    #files = [dir_path + "onepro00050.dat"]

    mkdir(outputdatapath)
    is_initial = True
    for file in files:
        print(file)
        time,rad,rho,pre,vel = ReadData(file)
        
        match = re.search(r'\d+', file)
        fileindex = match.group()

        PlotRadData(fileindex,time,rad,rho,pre,vel)

def mkdir(path):
    import os
    if not os.path.isdir(path):
        os.makedirs(path)

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
    
    return t, rad, rho, pre, vel

def PlotRadData(num,time,rad,rho,pre,vel):
  from matplotlib import ticker, cm, colors
  #######################################
  # Space Time diagram
  #######################################
  global outputdatapath
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
  time = time/year
  rad  = rad/pc
  timetxt=r"$T=$"+str(time)+" [year]"

  outputfile=outputdatapath+ "denone"+ num +".png"
  fig1 = plt.figure(figsize=(6,4.5))
  ax1 = fig1.add_subplot(1,1,1)
  ax1.plot(rad,rho,linewidth=2)
  ax1.set_xlabel(r"$r\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax1.set_ylabel(r"$\rho\ [{\rm 1/cm^3}]$", fontsize=fsizeforlabel)
  ax1.text(0.7, 0.7,timetxt, ha='center', transform=ax1.transAxes)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)

  outputfile=outputdatapath+ "preone"+ num +".png"
  fig1 = plt.figure(figsize=(6,4.5))
  ax1 = fig1.add_subplot(1,1,1)
  ax1.plot(rad,pre,linewidth=2)
  ax1.set_xlabel(r"$r\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax1.set_ylabel(r"$p\ [{\rm erg/cm^3}]$", fontsize=fsizeforlabel)
  ax1.text(0.7, 0.7,timetxt, ha='center', transform=ax1.transAxes)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)

  outputfile=outputdatapath+ "velone"+ num +".png"
  fig1 = plt.figure(figsize=(6,4.5))
  ax1 = fig1.add_subplot(1,1,1)
  ax1.plot(rad,vel,linewidth=2)
  ax1.set_xlabel(r"$r\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax1.set_ylabel(r"$v\ [{\rm cm/s}]$", fontsize=fsizeforlabel)
  ax1.text(0.7, 0.7,timetxt, ha='center', transform=ax1.transAxes)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)

if __name__ == "__main__":
    Main()
