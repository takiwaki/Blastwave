##########################################
# Makefile for Turbulence Studies
##########################################

##########################################
# Programs
##########################################

ana=Analysis.x

######################	
# complile options
######################
fc=ifort
#fopt=-g -traceback -O2
fopt=-g -traceback -check all -fpe0

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

TWOPRO    = $(patsubst ${dird}/unf%.dat,${dira}/twopro%.dat,$(BIN))
DENTWOPNG = $(patsubst ${dira}/twopro%.dat,${dirf}/dentwo%.png,$(TWOPRO))
VELTWOPNG = $(patsubst ${dira}/twopro%.dat,${dirf}/veltwo%.png,$(TWOPRO))
PRETWOPNG = $(patsubst ${dira}/twopro%.dat,${dirf}/pretwo%.png,$(TWOPRO))

# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

denonemovie=${dirm}/anidenone.mp4
velonemovie=${dirm}/anivelone.mp4
preonemovie=${dirm}/anipreone.mp4

dentwomovie=${dirm}/anidentwo.mp4
veltwomovie=${dirm}/aniveltwo.mp4
pretwomovie=${dirm}/anipretwo.mp4

movies = ${denonemovie} ${velonemovie} ${preonemovie} ${dentwomovie} ${veltwomovie} ${pretwomovie} 

##########################################
# Time evolution
##########################################

timefile=t-E.png

##########################################
# Proceadures
##########################################

all: ${exe} ${movies} ${timefile}

.PHONY: all

##########################################
# Time evolution
##########################################
.PHONY: timeevolution
timeevolution: ${timefile} ## time evolution 

${timefile}: PlotTime.plt t-prof.dat
	gnuplot PlotTime.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh

#################
# ffmpeg
#################
.PHONY: movie
movie: ${movies} ## Make movie

${denonemovie}: MakeMovie.sh ${DENONEPNG}
	./MakeMovie.sh denone
${velonemovie}: MakeMovie.sh ${VELONEPNG}
	./MakeMovie.sh velone
${preonemovie}: MakeMovie.sh ${PREONEPNG}
	./MakeMovie.sh preone

${dentwomovie}: MakeMovie.sh ${DENTWOPNG}
	./MakeMovie.sh dentwo
${veltwomovie}: MakeMovie.sh ${VELTWOPNG}
	./MakeMovie.sh veltwo
${pretwomovie}: MakeMovie.sh ${PRETWOPNG}
	./MakeMovie.sh pretwo

#################
# gnuplot
#################

.PHONY: 1Dsnaps 2Dsnaps
1Dsnaps: ${DENONEPNG} ${VELONEPNG} ${PREONEPNG} ## Make one-dimensional plot

${DENONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

${VELONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

${PREONEPNG}: MakeSnap.sh Plot1D.plt ${ONEPRO}
	./MakeSnap.sh onepro Plot1D.plt

2Dsnaps: ${DENTWOPNG} ${VELTWOPNG} ${PRETWOPNG}  ## Make two-dimensional plot

${DENTWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh twopro Plot2D.plt
${VELTWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh twopro Plot2D.plt

${PRETWOPNG}: MakeSnap.sh Plot2D.plt ${TWOPRO}
	./MakeSnap.sh twopro Plot2D.plt

#################
# analysis
#################
run: ${ONEPRO} ${TWOPRO} ## run the analysis

count-number: ${countfile} ## count bin files to prepare the analysis

${countfile}: CountBindata.sh
	./CountBindata.sh

${ONEPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

${TWOPRO}: ${ana} ${BIN} ${countfile}
	./${ana}

${TOT}: ${ana} ${BIN} ${countfile}
	./${ana}

#################
# simulation
#################
build: ${ana} ## Compile the analysis code

${ana}: Analysis.f90
	${fc} ${fopt} $< -o ${ana}


#################
# clean up
#################
.PHONY: clean allclean
clean: ## Delete intermediate files to compile
	rm -f  ${ana} *.o *.mod *~

allclean: ## Delete data
	rm -fr ${dirm} ${dirf} ${dira} ${countfile}

.PHONY: help
help: ## Display this help screen
	@grep -E '^[0-9a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'
