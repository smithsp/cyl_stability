set term postscript eps enhanced color
set output "finite_elements/left_splines_deriv.eps"
set key outside below
set xrange [0:.4]
set title 'Derivative left end elements '
plot 'finite_element_values.txt' using 1:15 with lines title 'Left0 Bspline Derivative ', \
'finite_element_values.txt' using 1:11 with lines title 'Left1 Bspline Derivative', \
'finite_element_values.txt' using 1:12 with lines title 'Left2 Bspline Derivative', \
'finite_element_values.txt' using 1:13 with lines title 'Middle Bspline Derivative'
set term x11
set output
