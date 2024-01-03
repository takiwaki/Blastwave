#!/bin/bash

gscript=$1
echo "use "$gscript

fstfile=`ls -1 ./output/rpr*.dat 2>/dev/null | head -1`
echo $fstfile
declare -i fstnum=`echo  ${fstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`
echo $fstnum

lstfile=`ls -1 ./output/rpr*.dat 2>/dev/null | tail -1`
echo $lstfile
declare -i lstnum=`echo  ${lstfile##*/} | tr -cd '0123456789\n' |sed -e 's/^0\+\([0-9]\+\)$/\1/'`
echo $lstnum

outfile="t-r-pro.dat"
echo " " > $outfile

for n in $(seq ${fstnum} ${lstnum}); do
file=`printf "./output/rpr%05d.dat\n" "${n}"`
timev=`awk 'NR==1{print($3)}' ${file}`
echo ${timev}

awk 'NR>2{print('${timev}', $1, $2, $3, $4 )}' ${file} > ./temp.dat
echo "" >> ./temp.dat
cat ./temp.dat >> ./${outfile}

done

rm ./temp.dat
