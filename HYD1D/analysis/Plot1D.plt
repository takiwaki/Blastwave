#####################
# initialize
######################

unset key

# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" 
# crop 

if (exist("ifnum")==0 ) ifnum=50

ifnames = sprintf("output/onepro%05d.dat",ifnum)

command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)
flag=0
flag=system(command)
print  ifnames." found"

command = sprintf("awk 'NR==1 {print $3}' %s",ifnames)
time   = system(command)
time = time + 0.0
year=60*60*24*365
time = time/year
timeunit=" year"
timetxt = sprintf("%g",time).timeunit
print "time= ".timetxt

pc=3.085677581e18 
set xlabel "Radius [pc]" offset 0,0

##########################
# Pressure
##########################

ofname = sprintf("figures/preone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Pressure [erg/cm^3]" offset 0,0

plot  \
  ifnames u ($1/pc):3 w l lw 6 \

print "output ".ofname

##########################
# density
##########################_

ofname = sprintf("figures/denone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "Density [1/cm^3]" offset 0,0

plot  \
  ifnames u ($1/pc):2 w l lw 6 \

print "output ".ofname

##########################
# velocity
##########################

ofname = sprintf("figures/velone%05d.png",ifnum)
set output ofname

set label 1 timetxt at screen 0.45, screen 0.845
set ylabel "velocity [km/s]" offset 0,0

plot  \
  ifnames u ($1/pc):4 w l lw 6 \


print "output ".ofname


reset
set term pop
