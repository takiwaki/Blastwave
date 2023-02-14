# 1D Blastwave simulation

[Go to top](../README.md)  

## How to run

### compile 
This is the instruction for spring school of division of science. First login the server, more.

    ssh <your account>@more.cfca.nao.ac.jp
    
Then copy the source code.

    cd /cfca-work/<your account>
    cp -r /cfca-work/dos00/Blastwave .
To run the code, you need to compile 'Simulation.f90'.
    
    cd Blastwave/HYD1D
    module load intel
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
    
    ./Analyis.x
    
The output is saved in `output/`.
### 1D plots and animation.
If you need 1D snapshots, use the following command. Using `output/rpr*.dat` (Radial PRofile), image files are made and save as `figures/*.png` (e.g., `den00050.png`).
    
    gnuplot radpro.plt
    ls figures/
    display den00050.png
    
All snapshots are made by the following command. 
    
    make 1Dsnaps
   
To make movie from the files. Type as follows.

    make movies
   
The movie files in saved in `movies/`. You can see the movie with the following command.

    ls movies/
    mplayor movies/ani???.mp4
   
### See energy conservation

    display t-E.png
    
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
