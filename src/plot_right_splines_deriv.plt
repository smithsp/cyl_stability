set term postscript eps enhanced color
set output "finite_elements/right_splines_deriv.eps"
set key outside below
set xrange [.5:.9]
set title 'Right end derivative elements'
plot 'finite_element_values.txt' using 1:20 with lines title 'RightN-2 Bspline Derivative', \
'finite_element_values.txt' using 1:14 with lines title 'RightN-1 Bspline Derivative', \
'finite_element_values.txt' using 1:21 with lines title 'RightN Bspline Derivative', \
'finite_element_values.txt' using 1:22 with lines title 'RightN+1 Bspline Derivative'
set term x11
set output
