##########################################
# Makefile for Blastwave Analysis
##########################################

##########################################
# Programs
##########################################

ana=Analysis.x

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

dird := ../bindata
dira := output
dirf := figures
dirm := movies

countfile=control.dat

######################
# files
######################

filename = ${dird}/unf00001.dat
fileexists = $(shell ls | grep ${filename})
ifeq (${fileexists}, ${filename})
    BIN= $(shell ls ${dird}/unf*.dat)
else
    BIN= ${filename}
endif

# 1D radial profile
ONEPRO    = $(patsubst ${dird}/unf%.dat,${dira}/onepro%.dat,$(BIN))
DENONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/denone%.png,$(ONEPRO))
VELONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/velone%.png,$(ONEPRO))
PREONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/preone%.png,$(ONEPRO))
XCMONEPNG = $(patsubst ${dira}/onepro%.dat,${dirf}/xcmone%.png,$(ONEPRO))

# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

denonemovie=${dirm}/anidenone.mp4
velonemovie=${dirm}/anivelone.mp4
preonemovie=${dirm}/anipreone.mp4

movies = ${denonemovie} ${velonemovie} ${preonemovie}
##########################################
# Time evolution
##########################################

timeevo=t-E.png

timeradrho=t-r-rho.png

##########################################
# Proceadures
##########################################

all: ${ana} ${movies} ${timeevo}

.PHONY: all

##########################################
# Time evolution
##########################################
.PHONY: timeevolution
timeevolution: ${timefile} ## time evolution 

${timeevo}: PlotTime.plt t-prof.dat
	gnuplot PlotTime.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh

##########################################
# Space time diagram
##########################################
${timeradrho}: PlotTimeRad.plt t-r-pro.dat
	gnuplot PlotTimeRad.plt

t-r-pro.dat: MakeTimeRad.sh ${RPR}
	./MakeTimeRad.sh

#################
# ffmpeg
#################
.PHONY: movie
movie: ${movies}

${denonemovie}: MakeMovie.sh ${DENONEPNG}
	./MakeMovie.sh denone

${velonemovie}: MakeMovie.sh ${VELONEPNG}
	./MakeMovie.sh velone

${preonemovie}: MakeMovie.sh ${PREONEPNG}
	./MakeMovie.sh preone

#################
# gnuplot
#################
.PHONY: 1Dsnaps
1Dsnaps: ${DENONEPNG} ${VELONEPNG} ${PREONEPNG} ## Make one-dimentional plot

${DENONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

${VELONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

${PREONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

#################
# analysis
#################
run: ${ONEPRO} ## run the analysis

count-number: ${countfile} ## count bin files to prepare the analysis

${countfile}: CountBindata.sh
	./CountBindata.sh

${ONEPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

${TOT}: ${ana} ${BIN} ${countfile}
	./${ana}

#################
# analysis
#################

build: ${ana} ## Compile the analyis code

${ana}: Analysis.f90
	${fc} ${fopt} $< -o ${ana}

#################
# clean up
#################
.PHONY: clean allclean
clean: ## Delete intermediate files to compile
	rm -f ${ana} *.o *.mod *~

allclean: ## Delete data
	rm -fr ${dirm} ${dirf}  ${dira} ${countfile}

#################
# help
#################
.PHONY: help
help: ## Display this help screen
	@grep -E '^[0-9a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'
