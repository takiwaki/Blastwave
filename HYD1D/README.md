# 1D Blastwave simulation

[Go to top](../README.md)  

## How to run and analyse

This is the instruction for spring school of division of science.

### login and go to work directory 
First login the server, more.

    ssh <your account>@more.cfca.nao.ac.jp
    
Then go to work directory.

    mkdir /cfca-work/<your account>
    cd /cfca-work/<your account>

Copy the programs.
    
    cp -r /cfca-work/dos04/Blastwave .
   
### compile 
To run the code, you need to compile `Simulation.f90`.
    
    cd Blastwave/HYD1D
    module load intel/2022
    make Simulation.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    qsub pbs_more.sh
    
The simulation data is saved in `bindata/`.

### analysis
Go to analysis server. Here ?? below is 09-14. To analyze the data, let us make `Analysis.x`.
    
    ssh <your account>@an??.cfca.nao.ac.jp

Then go to work directory.
    
    cd /cfca-work/<your account>/Blastwave/HYD1D/analysis
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    ./CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Analysis.x
    
The output is saved in `output/`.
### 1D plots and animation.
If you need 1D snapshots, use the following command. Using `output/onepro*.dat` (1DProfile), image files are made and save as `figures/*.png` (e.g., `denone00050.png`).
    
    gnuplot Plot1D.plt
    ls figures/
    display figures/denone00050.png

Compare the figure with the following one.

    gnuplot
    gnuplot> plot "output/onepro00050.dat" u 1:2 w l
    
    
All snapshots are made by the following command. 
    
    make 1Dsnaps
   
To make movie from the files. Type as follows.

    make movie
   
The movie files in saved in `movies/`. You can see the movie with the following command.

    ls movies/
    mplayer movies/ani???.mp4
   
### See energy conservation

    make t-E.png
    display t-E.png
    
### Space Time Diagram

    make t-r-rho.png
    display t-r-rho.png
    
Compare the figure with the following one.

    gnuplot
    gnuplot> set view map
    gnuplot> splot "t-r-pro.dat" u 1:2:3 w pm3d
    
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
# How to change parameter
Let us try to change the parameters. Before change it. Change the name of the previous directory.
If you are still in `analysis`, change directory.

    cd ../..
    mv HYD1D HYD1D-model1
    cp -r /cfca-work/dos04/Blastwave/HYD1D .
    cd HYD1D

Now you are in `HYD1D`. You can change the number of numerical grid in `module commons` in `Simulation.f90`.
<pre>
      integer,parameter::izones=200
</pre>
Also you can change the timescale and simulation region.
<pre>
      real(8),parameter:: timemax=2.0d3*year
      real(8),parameter:: x1min=0.0d0,x1max=10.0d0*pc
</pre>
In `subroutine GenerateProblem`, you can find the following part.
You can change the paramters as you like.

<pre>
      dr = 8.0d0*(x1a(is+1)-x1a(is)) ! 8 mesh
      write(6,*) "shell length [pc]",dr/pc

! parameters
      Mejecta = 10.0d0 ! M_sun
      EexpThermal = 0.8 ! 10^51 erg
      EexpKinetic = 0.2 ! 10^51 erg
      RhoMedium   = 1.0 ! 1/cm^3
      TMedium     = 1.0d4 ! 10^4 [K]
! circum stellar  medium
      rho2 = RhoMedium * mu ! Interstellar medium 1 [g/cm^3]
      pre2 = rho2* kbol * Tmedium
      vel2 = 0.0d0

! blastwave
      vol  = (4.0*pi/3.0d0*dr**3)-(4.0*pi/3.0d0*x1min**3) ! cm^3
      rho1 = (Mejecta*Msolar)/vol ! g/cm^3
      eexp = EexpThermal*(1.0d51) ! erg
      pre1 = eexp/vol*(gam-1.0d0) ! erg/cm^3
      vel1 = sqrt(EexpKinetic*1.0d51/vol/rho1) ! cm/s

      write(6,*) "Eex= ",eexp/1.0d51,"[10^51 erg]"
      write(6,*) "rho= ",rho1/mu,"[1/cm^3]"
      write(6,*) "vel= ",vel1   ,"[cm/s]"
      write(6,*) "pre= ",pre1   ,"[erg/cm^3]"
</pre>
