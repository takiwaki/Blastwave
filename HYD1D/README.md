# 1D Blastwave simulation

[Go to top](../README.md)  

## How to run

### compile 
This is the instruction for spring school of division of science. First login the server, more.

    ssh <your account>@more.cfca.nao.ac.jp
    
Then copy the source code.

    cd /cfca-work/<your account>
    mkdir /cfca-work/<your account>
    cp -r /cfca-work/dos04/Blastwave .
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
GO to analysis server. Here ?? below is 09-14. To analyze the data, let us make `Analysis.x`.
    
    ssh <your account>@an??.cfca.nao.ac.jp
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
    
    gnuplot radpro.plt
    ls figures/
    display figures/denone00050.png
    
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
    
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
# How to change paramete
You can change the number of numerical grid in  `module commons`.
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

! circum steller  medium
      rho2 = 1.0d0*mu ! Intersteller medium 1 [1/cm^3]
      pre2 = rho2* kbol *1.0d4 ! 10^4 [K]
      vel2 = 0.0d0

! blast wave
      vol  = (4.0*pi/3.0d0*dr**3)-(4.0*pi/3.0d0*x1min**3)
      frac = 0.8d0
      rho1 = (10.0d0*Msolar)/vol
      eexp = frac*(1.0d51) ! erg
      pre1 = eexp/vol*(gam-1.0d0)  
      vel1 = sqrt((1.0d0-frac)*eexp/vol/rho1)

      write(6,*) "Eex= ",frac   ,"[10^51 erg]"
      write(6,*) "rho= ",rho1/mu,"[1/cm^3]"
      write(6,*) "vel= ",vel1   ,"[cm/s]"
      write(6,*) "pre= ",pre1   ,"[erg/cm^3]"
</pre>
