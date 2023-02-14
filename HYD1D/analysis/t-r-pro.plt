
reset

pngflag=1

if(pngflag==1) outflag=1
if(outflag==1) set terminal push

if(pngflag ==1) set terminal pngcairo enhanced font "Helvetica, 18" crop

input="t-r-pro.dat"

set pm3d map

##########################################
# Space Time Diagram
##########################################
if(pngflag ==1)set output "t-r-rho.png"

set xlabel "Time [year]" offset 0,1
set xtic offset 0,0.5
set xtic 200

set ylabel "Radius [pc]" offset 1.0, 0.0
set ytic offset 0.5,0.0

set palette defined ( 0 "white", 1.0 "white", 5.0 "black")

set cbrange [*:10]
set xrange [*:1000]
set yrange [*:5]
splot input  u 1:2:3 notitle \

##########################################
# Finalize
##########################################

reset
if(outflag==1) set terminal pop