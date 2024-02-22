# 2D Blastwave simulation

[Go to top](../README.md)  

## How to run and analyse

This is the instruction for spring school of division of science.

### login and go to work directory 
First login the server, `more.cfca.nao.ac.jp`.

    ssh <your account>@more.cfca.nao.ac.jp
    
Then, go to work directory. Make it if that does not exist.

    mkdir /cfca-work/<your account>
    cd /cfca-work/<your account>

Copy the programs.　If you do not did it before. 
    
    cp -r /cfca-work/dos04/Blastwave .

Keep the original program as it is.
    
    cd Blastwave/
    mv HYD2D HYD2D_original
   
### Making your model 
Start the simulation by copying the original file. You can name the directory as you like. `_model1` is an example.
    
    cp -r HYD2D_original HYDD_model1
    cd HYD2D_model1


### compile 
To run the code, you need to compile `Simulation.f90`.
    
    module load intel/2022
    make Simulation.x
    
Then `Simulation.x` is made in this directory.

### run
Let's run the code.
    
    qsub pbs_more.sh
    
The simulation data is saved in `bindata/`.
    
    ls bindata/
    
### analysis
Open another terminal and go to analysis server, `an??.cfca.nao.ac.jp`. Here ?? is 09-14. To analyze the data, let us make `Analysis.x`.
    
    ssh <your account>@an??.cfca.nao.ac.jp

Then go to the work directory. Change `_model1` to the name you used.
    
    cd /cfca-work/<your account>/Blastwave/HYD2D_model1/analysis
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    ./CountBindata.sh
    cat control.dat
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Analyis.x
    
The output is saved in `output/`.
    
    ls output/
    
### 2D plots and animation.
If you need 2D snapshots, use the following command. Using `output/twopro*.dat` (2D Profile), image files are made and save as `figures/*.png` (e.g., `dentwo00050.png`).
    
    gnuplot Plot2D.plt
    ls figures/
    display figures/dentwo00050.png
    display figures/pretwo00050.png
    display figures/veltwo00050.png

Compare the figure with the following one.

    gnuplot
    gnuplot> set view map
    gnuplot> set title "density"
    gnuplot> splot "output/twopro00050.dat" u 1:2:3 w pm3d
    
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
      
If you want th delete all the analysis, type `make allclean`
