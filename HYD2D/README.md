# 2D Blastwave simulation

[Go to top](../README.md)  

## How to run

### compile 
This is the instruction for spring school of division of science. First login the server, more.

    ssh <your account>@more.cfca.nao.ac.jp
    
Then copy the source code.

    cd /cfca-work/<your account>
    cp -r /cfca-work/dos04/Blastwave .
To run the code, you need to compile `Simulation.f90`.
    
    cd Blastwave/HYD2D
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
    cd /cfca-work/<your account>/Blastwave/HYD2D/analysis
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    ./CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Analyis.x
    
The output is saved in `output/`.
### 2D plots and animation.

If you need 2D snapshots, use the following command. Using `output/twopro*.dat` (2D Profile), image files are made and save as `figures/*.png` (e.g., `dentwo00050.png`).
    
    gnuplot Plot2D.plt
    ls figures/
    display dentwo00050.png
    
All snapshots are made by the following command. 
    
    make 2Dsnaps
   
To make movie from the files. Type as follows.

    make movie
   
The movie files in saved in `movies/`. You can see the movie with the following command.

    ls movies/
    mplayor movies/ani???.mp4
    
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
