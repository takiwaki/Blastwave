# Supernova lightcurve

[Go to top](../README.md)  

## How to run

### compile 
This is the instruction for spring school of division of science. First login the server, an??. ?? is 09-14.

    ssh <your account>@an??.cfca.nao.ac.jp
    
Then copy the source code.

    cp -r /cfca-work/dos04/lightcurve .

To run the code, you need to compile.
    
    module load intel
    ifort lightcurve.f90
    
Then `a.out`is made in this directory.

### run
Let's run the code.
    
    ./a.out
    
### analysis
See the result.
    
    gnuplot
    plot "lightcurve.dat"

### Annotation
To make a figure.

    gnuplot lc.plt
    display lightcurve.png
    
