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
fc=ifort -extend-source
fopt=-g -traceback -O2
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
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dird} ${exe} *.o *.mod *~
