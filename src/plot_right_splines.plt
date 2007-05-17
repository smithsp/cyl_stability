set term postscript eps enhanced color
set output "finite_elements/right_splines.eps"
set key outside below
set xrange [.5:.9]
set title 'Right end elements'
plot 'finite_element_values.txt' using 1:17 with lines title 'RightN-2 Bspline', \
'finite_element_values.txt' using 1:10 with lines title 'RightN-1 Bspline', \
'finite_element_values.txt' using 1:18 with lines title 'RightN Bspline', \
'finite_element_values.txt' using 1:19 with lines title 'RightN+1 Bspline'
set term x11
set output
