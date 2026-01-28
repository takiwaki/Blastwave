##########################################
# Makefile for Turbulence Studies
##########################################

##########################################
# Programs
##########################################

exe=Simulation.x

######################	
# complile options
######################
fc=ifort -diag-disable=10448 -extend-source
fopt=-g -traceback -O2
foptomp=-qopenmp
#fopt=-g -traceback -check all -fpe0

##########################################
# Timeseqencial data and figures
##########################################

######################
# directory
######################

dird := bindata

#################
# simulation
#################

${exe}: Simulation.f90
	${fc} ${fopt} ${foptomp} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -rf ${dird} ${exe} *.o *.mod *~
