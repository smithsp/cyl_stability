set term postscript eps enhanced color solid
set output "finite_elements/bspline.eps"
set size .50,.250
set key off
#set yrange [-.1:.9]
set xrange [0:.4]
set xtics 0,.1,.4
set ytics 0,.2,.6
set title 'Cubic B-Spline
plot 'finite_element_values.txt' using 1:9  with lines 
set term x11
set output
