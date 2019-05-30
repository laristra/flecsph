#!/bin/bash
#
# Copyright (c) 2017 Triad National Security, LLC
# All rights reserved
#
#  A utility script to plot a variable using gnuplot
#  Requirements:
#  - must be run after slicer.sh
#  - uses xdp.x.dat and xuv.x.dat
#

# Options for plotting:
xmin=0.0
xmax=1.5
ymin=
ymax=

# parse variable name
if [ $# -lt 1 ]; then
  echo "Usage: ./$0 <var> [-st <step>|-t <time>]"
  echo "where var = {rho|p|u|ux}"
  exit -1
fi
var=$1
case $var in
u|rho|p|P|vx) ;;
*)
  echo "ERROR: bad variable name" >> /dev/stderr
  exit -2
esac
shift

# determine step and time limits
xdpfile="xdp.x.dat"
xuvfile="xuv.x.dat"
if [ ! -f $xdpfile ]; then
  echo "ERROR: file $xdpfile does not exist." >> /dev/stderr
  echo "Make sure to run the slicer.sh script first " >> /dev/stderr
  exit -5
fi

linefirst=`grep time $xdpfile| head -n1`
linelast=`grep time $xdpfile| tail -n1`
st_first=`echo $linefirst | awk -F'[= ]*' '{print $3}'`
st_last=`echo $linelast | awk -F'[= ]*' '{print $3}'`
t_first=`echo $linefirst | awk '{print $4}'` 
t_last=`echo $linelast | awk '{print $4}'` 
echo "Read $xdpfile, found the following time / steps ranges:"
echo "steps:"
printf " - first:   %12d\n" $st_first
printf " - last:    %12d\n" $st_last
echo "time:"
printf " - initial: %12.5e\n" $t_first
printf " - final:   %12.5e\n" $t_last

# parse step and time
st=$st_first
tm=$t_first
if [ $# -ne 0 ]; then
  if [ $# -ne 2 ]; then
    echo "ERROR: wrong number of arguments" >> /dev/stderr
    exit -4
  fi
  set -e
  case $1 in
  -st) st=$2 
    bool=`awk 'BEGIN{print ('$st'<'$st_first') || ('$st'>'$st_last')}'`
    if [ $bool -eq 1 ]; then 
      echo "ERROR: step out of range"
      exit -7
    fi
    tm=`grep "step=$st " $xdpfile| awk '{print $4}'`
  ;;
  -t)  tm=$2
    bool=`awk 'BEGIN{print ('$tm'<'$t_first') || ('$tm'>'$t_last')}'`
    if [ $bool -eq 1 ]; then 
      echo "ERROR: time out of range"
      exit -7
    fi
    st=`grep "time" $xdpfile| awk -F'[= ]*' '($5>'$tm'){print $3-1; exit} ($5=='$tm'){print $3; exit}'`
    ;;
  *)
    echo "ERROR: bad argument" >> /dev/stderr
    exit -3
  esac
fi
echo "Using:"
echo " - step = $st"
echo " - time = $tm"

# generate gnuplot script
st0=`printf "%06d" $st`
varplotfile=${var}.x.${st0}.gnu
cat > $varplotfile << EOF
#!`which gnuplot`

set term postscript eps color enhanced
set output "${varplotfile%.gnu}.eps"

set xlabel "x"
set ylabel "$var"
set size 0.7
set xrange [$xmin:$xmax]
set yrange [$ymin:$ymax]
set key top right samplen 1.2
EOF
case $var in 
rho)
  echo 'p "'$xdpfile'" u 1:2 i '$st' w l lw 2 t "t = '$tm'"' >> $varplotfile
  ;;
p|P)
  echo 'p "'$xdpfile'" u 1:3 i '$st' w l lw 2 t "t = '$tm'"' >> $varplotfile
  ;;
u)
  echo 'p "'$xuvfile'" u 1:2 i '$st' w l lw 2 t "t = '$tm'"' >> $varplotfile
  ;;
vx)
  echo 'p "'$xuvfile'" u 1:3 i '$st' w l lw 2 t "t = '$tm'"' >> $varplotfile
  ;;
esac
chmod +x $varplotfile
./$varplotfile
echo "The plot has been saved to a file ${varplotfile%.gnu}.eps"
echo "The script to generate the plot is ${varplotfile}. Modify it to your liking"

