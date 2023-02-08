# 1D Blastwave simulation

[Go to top](../README.md)  

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90'.
    
    make Simulation.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    qsub pbs_more.sh
    
The simulation data is saved in `bindata/`.

### analysis
GO to analysis server. To analyze the data, let us make `Analysis.x`.
    
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Analyis.x
    
The output is saved in `output/`.
### 2D plots and animation.
If you need 2D snapshots. 
    
    make 2Dsnaps
   
Using `output/vor*.dat`, image files are made and save as `figures/vor*.png`.
To make movie from the files. Type as follows.

    make movie
   
The movie files in saved in `movie/anivor`.

### spectrum
To obtain the spectrum
   
      make spectrum
      
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
