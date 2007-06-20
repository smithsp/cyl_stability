set term postscript eps enhanced color solid
set output "finite_elements/bspline_deriv.eps"
set size .50,.250
set key off
set yrange [-8:8]
set xrange [0:.4]
set xtics 0,.1,.4
set ytics -8,4,8
set title 'Cubic B-Spline Derivative'
plot 'finite_element_values.txt' using 1:13 with lines 
set term x11
set output
