#!/bin/bash

# TODO: add header and description
 
for number in `seq -w 0 99999`
do
  if [ ! -e "output_sodtube_$number.txt" ]
  then
  echo "File output_sodtube_$number.txt does not exists, EXITING"
    exit 1
  fi
  echo "Generating output_sodtube_$number.txt"
  gnuplot -e "
    set terminal png size 1024,768;
    set output 'output_sodtube_$number.png'; 
    set multiplot layout 2,2; 
    set xrange [0:$1]; 
    set yrange [0:1.1];
    plot 'output_sodtube_$number.txt' using 1:2 title 'Density'; 
    set xrange [0:$1];  
    set yrange [0:1.1];
    plot 'output_sodtube_$number.txt' using 1:3 title 'Pressure'; 
    set xrange [0:$1];  
    set yrange [1.6:3.3];
    plot 'output_sodtube_$number.txt' using 1:4 title 'InternalEnergy'; 
    set xrange [0:$1];
    set yrange [-0.05:1.1];
    plot 'output_sodtube_$number.txt' using 1:5 title 'Velocity'"
done

