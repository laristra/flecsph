#!/usr/bin/gnuplot -persist

box_length = 2.3
box_width  = 1.0

airfoil_size = 2.0
airfoil_thickness = 0.05
airfoil_camber = 0.1
airfoil_attack_angle = 10.0
airfoil_anchor_x = -1.6
airfoil_anchor_y = -0.3

min(a,b) = ((a>b)?(b):(a))
max(a,b) = ((a>b)?(a):(b))
wmargin = 0.1*airfoil_size
xmin = airfoil_anchor_x
xmax = xmin + airfoil_size + box_length
ymin = -box_width/2
ymax =  box_width/2
f = (2*wmargin + ymax - ymin)/(2*wmargin + xmax - xmin)

set t png size 1024,768 
set o "wt_airfoil.png"

set size ratio -1
set xlab "x"
set ylab "y"
set samples 3000

a = airfoil_attack_angle*pi/180
xrot(x,y) = x*cos(a) - y*sin(a)
yrot(x,y) = x*sin(a) + y*cos(a)
shape(x)  = (x<0?1/0:1)*(x>airfoil_size?1/0:1)*airfoil_thickness*x*sqrt(airfoil_size**2 - x*x);
camber_line(x) = airfoil_camber*sin(pi*x/2.)
set obj 4 rect from xmax - box_length,-box_width/2 to xmax,box_width/2 fc rgb "#0088ff" fs noborder
#set style arrow 16
set arrow from xmax - box_length/2, 0 to xmax - box_length/2*1.3, 0 ls 1 lw 2 lc rgb "#004400" 

set xrange [xmin - wmargin : xmax + wmargin]
set yrange [ymin - wmargin : ymax + wmargin]
set parametric

upper_surface_x(t) = xrot(t, camber_line(t) + shape(t)) + airfoil_anchor_x
upper_surface_y(t) = yrot(t, camber_line(t) + shape(t)) + airfoil_anchor_y
lower_surface_x(t) = xrot(t, camber_line(t) - shape(t)) + airfoil_anchor_x
lower_surface_y(t) = yrot(t, camber_line(t) - shape(t)) + airfoil_anchor_y

p upper_surface_x(t),upper_surface_y(t) w l lc rgb "#000088" lw 2 t "" \
, lower_surface_x(t),lower_surface_y(t) w l lc rgb "#000088" lw 2 t "" \

