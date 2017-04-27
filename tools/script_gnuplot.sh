#!/bin/bash
 
for number in `seq -w 0 100`
do
  if [ ! -e "output_sod_$number.txt" ]
  then
  echo "File output_sod_$number.txt does not exists, EXITING"
    exit 1
  fi
  echo "Generating output_sod_$number.txt"
  gnuplot -e "set terminal png size 800,800; set output 'output_$number.png'; set xrange [0:1]; set yrange [0:1.5]; plot 'output_sod_$number.txt'"
done

