##########################################
# Makefile for Blastwave Simulations
##########################################

##########################################
# Programs
##########################################

exe=Simulation.x

######################	
# complile options
######################
fc=ifort
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
.PHONY: compile
compile: ${exe} ## Compile the prgram

${exe}: Simulation.f90
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
.PHONY: clean
clean:
	rm -f ${exe} *.o *.mod *~

.PHONY: help
help: ## Display this help screen
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'
