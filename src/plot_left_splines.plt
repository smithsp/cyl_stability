set term postscript eps enhanced color
set output "finite_elements/left_splines.eps"
set key outside below
#set yrange [-.1:.9]
set xrange [0:.4]
set title 'Left end elements'
plot 'finite_element_values.txt' using 1:16 with lines title 'Left0 Bspline', \
'finite_element_values.txt' using 1:7 with lines title 'Left1 Bspline', \
'finite_element_values.txt' using 1:8 with lines title 'Left2 Bspline', \
'finite_element_values.txt' using 1:9 with lines title 'Left3 Bspline'
set term x11
set output
