
##########################################
# parameters
##########################################

# Range of the plot [pc]
srange=10

# Format of the output
pngflag=1

##########################################
# initialize
##########################################

# color palette
#set palette rgbformulae 21,22,23 # black red yellow white
#set palette functions sqrt(gray), gray**3, sin(gray*2*pi) # black blue red yellow
#set palette rgb 33,13,10 # blue-green-yellow-red
#set palette define (0.0 "black",0.05 "blue", 0.5 "red", 0.75 "yellow", 1.0 "yellow")
#set palette define (-1.0 "blue", 0.0 "white", 1.0 "red")

# Lines
set style line 10 lt  1 lw 4 lc rgb "white" 

####################
# Input control
####################

# ifnum : Input File NUMber
if (exist("ifnum")==0 ) ifnum=50

print ifnum
# Scalar
ifnames = sprintf("output/rtp%05d.dat",ifnum)

# Stop Ploting
command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)

flag=0
flag=system(command)
if(flag ne "0") print ifnames." not found"; quit

# Extract Time
print  ifnames." found"

command = sprintf(" head -n 1 %s | sed 's/#//' ",ifnames)
time   = system(command)
timetxt = "time=".time." yr"
print timetxt

# Showing Time
set label timetxt at screen 0.65, screen 0.85


####################
# Output control
####################

# OUTPUT PNG
if (pngflag==1) set term push
if (pngflag==1) set term pngcairo  enhanced font "Helvetica, 12" size 550,500

####################
# Annotation
####################

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set size ratio -1
set view map
unset key

set size 1.0
unset origin


# vertical and horizontal axis
set origin 0.0,0.0
set xlabel "X [pc]" offset 0,0
set xtics 1
set ylabel "Z [pc]" offset 0,0
set ytics 1
vr=srange/10

# Position of color bar
set colorbox horizontal user origin 0.235, 0.87 size 0.5, 0.04
set cbtics offset 0,3.2

# Range of color bar
# range of the variable
#cmin=0
#cmax=5
#set cbrange [cmin:cmax]

####################
# Plot
####################


# Density
ofname = sprintf("figures/dnt%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname


splot [-srange:srange][-srange:srange] \
  ifnames u ( $1*sin($2)):($1*cos($2)):($1<srange?($3):NaN) w pm3d \
, ifnames u (-$1*sin($2)):($1*cos($2)):($1<srange?($3):NaN) w pm3d \


# Pressure
ofname = sprintf("figures/prt%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname

splot [-srange:srange][-srange:srange] \
  ifnames u ( $1*sin($2)):($1*cos($2)):($1<srange?($4):NaN) w pm3d \
, ifnames u (-$1*sin($2)):($1*cos($2)):($1<srange?($4):NaN) w pm3d \

# Velocity
ofname = sprintf("figures/vlt%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname

splot [-srange:srange][-srange:srange] \
  ifnames u ( $1*sin($2)):($1*cos($2)):($1<srange?($5):NaN) w pm3d \
, ifnames u (-$1*sin($2)):($1*cos($2)):($1<srange?($5):NaN) w pm3d \


unset label

##########################################
# finalize
##########################################

unset output
if(pngflag==1)set term pop

undefine ifnum