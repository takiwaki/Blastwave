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
fc=ifort -extend-source
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
RPR   = $(patsubst ${dird}/unf%.dat,${dira}/rpr%.dat,$(BIN))
DENPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/den%.png,$(RPR))
VELPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/vel%.png,$(RPR))
PREPNG= $(patsubst ${dira}/rpr%.dat,${dirf}/pre%.png,$(RPR))

# intgrated values 
TOT   = $(patsubst ${dird}/unf%.dat,${dira}/tot%.dat,$(BIN))

##########################################
# Movies
##########################################

denmovie=${dirm}/aniden.mp4
velmovie=${dirm}/anivel.mp4
premovie=${dirm}/anipre.mp4

movies = ${denmovie} ${velmovie} ${premovie}
##########################################
# Time evolution
##########################################

timefile=t-E.png

##########################################
# Proceadures
##########################################

all: ${ana} ${movies}

.PHONY: all clean allclean

##########################################
# Time evolution
##########################################
${timefile}: TimeEvolution.plt t-prof.dat
	gnuplot TimeEvolution.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh

#################
# ffmpeg
#################

movie: ${moviefile}

${denmovie}: MakeMovie.sh ${DENPNG}
	./MakeMovie.sh den

${velmovie}: MakeMovie.sh ${VELPNG}
	./MakeMovie.sh vel

${premovie}: MakeMovie.sh ${PREPNG}
	./MakeMovie.sh pre

#################
# gnuplot
#################

${DENPNG}: MakeSnap.sh radpro.plt ${RPR}
	./MakeSnap.sh rpr radpro.plt

${velPNG}: MakeSnap.sh radpro.plt ${RPR}
	./MakeSnap.sh rpr radpro.plt

${PREPNG}: MakeSnap.sh radpro.plt ${RPR}
	./MakeSnap.sh rpr radpro.plt


#################
# analysis
#################

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

${RPR}: ${ana} ${BIN} ${countfile}
	./${ana}

${TOT}: ${ana} ${BIN} ${countfile}
	./${ana}

#################
# analysis
#################

${ana}: Analysis.f90
	${fc} ${fopt} $< -o ${ana}


#################
# clean up
#################
clean:
	rm -f ${ana} *.o *.mod *~

allclean:
	rm -fr ${dirm} ${dirf} ${dira} ${countfile}
