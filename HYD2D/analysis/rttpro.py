#!/usr/bin/python

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

dir_path = "./output/"

def Main():
    global dir_path
    files=GetFileList()
#    files = [dir_path + "rtp00050.dat"]

    is_initial = True
    for file in files:
        print(file)
        time,rad,the,rho,pre,vel = ReadData(file)
        
        match = re.search(r'\d+', file)
        fileindex = match.group()

        PlotRadTheData(fileindex,time,rad,the,rho,pre,vel)

def GetFileList():
    global dir_path
    filenames= dir_path+"rtp*.dat"
    files = glob.glob(filenames)
    files = sorted(files)
    return files

def ReadData(file):
    input = file
    try:
        results = np.genfromtxt(input,skip_header=3,delimiter='    ') # Read numbers
#        print(results)
        inputf= open(input, 'r')
        header= inputf.readline() # Read the fisrt line
        item= header.split()
        t = float(item[2])
        print(t)
        header= inputf.readline() # Read the second line
        item= header.split()
#        print(item)
        nrad = int(item[2])
        nthe = int(item[4])
        print(nrad,nthe)
        inputf.close()

    except IOError:
        print(" cannot open " + input )
        sys.exit()
    rad, the, rho, pre, vel =  np.split(results,5,1)
    rad=rad.reshape(nthe,nrad)
    the=the.reshape(nthe,nrad)
    rho=rho.reshape(nthe,nrad)
    pre=pre.reshape(nthe,nrad)
    vel=vel.reshape(nthe,nrad)

    return t, rad, the, rho, pre, vel

def PlotRadTheData(num,time,rad,the,rho,pre,vel):
  from matplotlib import ticker, cm, colors
  #######################################
  # Space Time diagram
  #######################################
  outputdatapath="./figures/"
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

  timetxt=r"$T=$"+str(time)+" [year]"

  outputfile=outputdatapath+ "dnt"+ num +".png"
  fig1,ax = plt.subplots(subplot_kw={'projection': 'polar'})
  ax.set_title(r"$\rho [{\rm 1/cm^3}]$")
  ax.set_thetalim(-np.pi, np.pi)
  ax.set_rlim(0, rad.max())
  ax.set_theta_zero_location("N")
  im1=ax.pcolormesh( the,rad,rho)
  im2=ax.pcolormesh(-the,rad,rho)
  ax.set_xlabel(r"$x\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax.set_ylabel(r"$z\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
  fig1.colorbar(im1)
  ax.text(0.8, 1.0,timetxt, ha='center', transform=ax.transAxes)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)

  outputfile=outputdatapath+ "prt"+ num +".png"
  fig1,ax = plt.subplots(subplot_kw={'projection': 'polar'})
  ax.set_title(r"$p [{\rm erg/cm^3}]$")
  ax.set_thetalim(-np.pi, np.pi)
  ax.set_rlim(0, rad.max())
  ax.set_theta_zero_location("N")
  im1=ax.pcolormesh( the,rad,pre)
  im2=ax.pcolormesh(-the,rad,pre)
  ax.set_xlabel(r"$x\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax.set_ylabel(r"$z\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
  fig1.colorbar(im1)
  ax.text(0.8, 1.0,timetxt, ha='center', transform=ax.transAxes)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)

  outputfile=outputdatapath+ "vlt"+ num +".png"
  fig1,ax = plt.subplots(subplot_kw={'projection': 'polar'})
  ax.set_title(r"$v_r [{\rm cm/s}]$")
  ax.set_thetalim(-np.pi, np.pi)
  ax.set_rlim(0, rad.max())
  ax.set_theta_zero_location("N")
  im1=ax.pcolormesh( the,rad,vel)
  im2=ax.pcolormesh(-the,rad,vel)
  ax.set_xlabel(r"$x\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax.set_ylabel(r"$z\ [{\rm pc}]$", fontsize=fsizeforlabel)
  ax.tick_params(labelbottom=False, labelleft=True, labelright=False, labeltop=False)
  fig1.colorbar(im1)
  ax.text(0.8, 1.0,timetxt, ha='center', transform=ax.transAxes)
  fig1.tight_layout()
  print("output"+outputfile)
  fig1.savefig(outputfile)


Main()
