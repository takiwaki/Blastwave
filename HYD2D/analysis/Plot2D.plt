
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
ifnames = sprintf("output/twopro%05d.dat",ifnum)

# Stop Ploting
command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)

flag=0
flag=system(command)
if(flag ne "0") print ifnames." not found"; quit

# Extract Time
print  ifnames." found"

command = sprintf("awk 'NR==1 {print $3}' %s",ifnames)
time   = system(command)
time = time + 0.0
year=60*60*24*365
time = time/year
timeunit=" year"
timetxt = sprintf("%g",time).timeunit
print "time=".timetxt

# Showing Time
set label timetxt at screen 0.65, screen 0.88


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
set bmargin 1
set lmargin 1
set rmargin 2

set size ratio -1
set view map
unset key

set size 1.0
unset origin


pc=3.085677581e18 
# vertical and horizontal axis
set origin 0.0,0.0
set xlabel "X [pc]" offset 0,1.0
set xtics 1
set xtics offset 0,0.5
set ylabel "Z [pc]" offset 1.0,0
set ytics 1
set ytics offset 0.5,0
vr=srange/10

# Position of color bar
#set colorbox horizontal user origin 0.235, 0.87 size 0.5, 0.04
#set cbtics offset 0,3.2

set colorbox user origin 0.83, 0.13 size 0.03, 0.7
set cbtics offset -0.5,0

# Range of color bar
# range of the variable
#cmin=0
#cmax=5
#set cbrange [cmin:cmax]

####################
# Plot
####################


# Density
set title "log(n [1/cm^3])"
ofname = sprintf("figures/dentwo%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname
set cbrange[-2:*]
set palette defined (0.00 "#440154", 0.25 "#3b528b", 0.50 "#21908d", 0.75 "#5dc863", 1.00 "#fde725" ) # viridis
splot [-srange:srange][-srange:srange] \
  ifnames u ( $1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?(log($3)):NaN) w pm3d \
, ifnames u (-$1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?(log($3)):NaN) w pm3d \


# Pressure
ofname = sprintf("figures/pretwo%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname
set title "log(p [erg/cm^3])"
set cbrange[-20:*]
set palette defined (0.00 "#000004",0.25 "#3b0f70",0.50 "#8c2981",0.75 "#de4968",1.00 "#fe9f6d" ) # magma

splot [-srange:srange][-srange:srange] \
  ifnames u ( $1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?(log($4)):NaN) w pm3d \
, ifnames u (-$1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?(log($4)):NaN) w pm3d \

# Velocity
set title "v_r [km/s]"
norm=1.0e5
ofname = sprintf("figures/veltwo%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname
set cbrange[*:*]
set palette defined (0.00 "#0d0887", 0.25 "#7e03a8", 0.50 "#cc4778", 0.75 "#f89540", 1.00 "#f0f921" ) #plasma

splot [-srange:srange][-srange:srange] \
  ifnames u ( $1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?($5/norm):NaN) w pm3d \
, ifnames u (-$1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?($5/norm):NaN) w pm3d \

# Edot
set title "log(edot [erg/cm^3/s])"
norm=1.0e-21
ofname = sprintf("figures/edttwo%05d.png",ifnum)
print ofname
if (pngflag==1) set output ofname
set cbrange[*:*]
set palette defined (0.00 "#000004",0.25 "#420a68",0.50 "#932667", 0.75 "#dd513a", 1.00 "#fba40a" ) # inferno

splot [-srange:srange][-srange:srange] \
  ifnames u ( $1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?(log($6)):NaN) w pm3d \
, ifnames u (-$1/pc*sin($2)):($1/pc*cos($2)):($1/pc<srange?(log($6)):NaN) w pm3d \


unset label

##########################################
# finalize
##########################################

unset output
if(pngflag==1)set term pop

undefine ifnum