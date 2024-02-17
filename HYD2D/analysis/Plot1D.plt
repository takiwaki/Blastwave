#####################
# initialize
######################

unset key

# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" size 640,480 
# crop 

if (exist("ifnum")==0 ) ifnum=100

ifnames = sprintf("output/rpr%05d.dat",ifnum)

command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)
flag=0
flag=system(command)
print  ifnames." found"

command = sprintf(" head -n 1 %s | sed 's/#//' ",ifnames)
time   = system(command)
print "time= ".time


set xlabel "Radius [pc]" offset 0,0


##########################
# Pressure
##########################

ofname = sprintf("figures/pre%05d.png",ifnum)
set output ofname

set label 1 time at screen 0.45, screen 0.845
set ylabel "Pressure [erg/cm^3]" offset 0,0

plot  \
  ifnames u ($1):3 w l lw 6 \


##########################
# density
##########################_

ofname = sprintf("figures/den%05d.png",ifnum)
set output ofname

set label 1 time at screen 0.45, screen 0.845
set ylabel "Density [1/cm^3]" offset 0,0

plot  \
  ifnames u ($1):2 w l lw 6 \

##########################
# velocity
##########################_

ofname = sprintf("figures/vel%05d.png",ifnum)
set output ofname

set label 1 time at screen 0.45, screen 0.845
set ylabel "velocity [km/s]" offset 0,0

plot  \
  ifnames u ($1):4 w l lw 6 \




reset
set term pop
