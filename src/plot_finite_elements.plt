reset
N=140
set term postscript eps enhanced color
set output "finite_elements/right_splines.eps"
set key outside below
unset key
set xrange [(N-4.)/N:1]
set xtics ("" (N-4.)/N,"" (N-3.)/N,"" (N-2.)/N,"" (N-1.)/N,"" 1.)
set ytics (0)
#set title 'Right end elements'
plot 'finite_element_values.txt' using 1:17 with lines title 'RightN-2 Bspline', \
'finite_element_values.txt' using 1:10 with lines title 'RightN-1 Bspline', \
'finite_element_values.txt' using 1:18 with lines title 'RightN Bspline', \
'finite_element_values.txt' using 1:19 with lines title 'RightN+1 Bspline'
set output
set term x11
reset
N=140

set term postscript eps enhanced color solid
set output "finite_elements/right_splines_deriv.eps"
set key outside below
unset key
set xrange [(N-4.)/N:1]
set xtics ("" (N-4.)/N,"" (N-3.)/N,"" (N-2.)/N,"" (N-1.)/N,"" 1.)
set ytics (0)
#set title 'Right end derivative elements'
plot 'finite_element_values.txt' using 1:20 with lines title 'RightN-2 Bspline Derivative', \
'finite_element_values.txt' using 1:14 with lines title 'RightN-1 Bspline Derivative', \
'finite_element_values.txt' using 1:21 with lines title 'RightN Bspline Derivative', \
'finite_element_values.txt' using 1:22 with lines title 'RightN+1 Bspline Derivative'
set output
set term x11
reset
N=140
set term postscript eps enhanced color solid
set output "finite_elements/left_splines.eps"
set key outside below
unset key
#set yrange [-.1:.9]
set xrange [0:4./N]
set xtics ("" 0, "" 1./N, "" 2./N, "" 3./N, "" 4./N)
set ytics (0)
#set title 'Left end elements'
plot 'finite_element_values.txt' using 1:16 with lines title 'Left0 Bspline', \
'finite_element_values.txt' using 1:7 with lines title 'Left1 Bspline', \
'finite_element_values.txt' using 1:8 with lines title 'Left2 Bspline', \
'finite_element_values.txt' using 1:9 with lines title 'Left3 Bspline'
set output
set term x11
reset
N=140
set term postscript eps enhanced color solid
set output "finite_elements/left_splines_deriv.eps"
set key outside below
unset key
set xrange [0:4./N]
set xtics ("" 0, "" 1./N, "" 2./N, "" 3./N, "" 4./N)
set ytics (0)
#set title 'Derivative left end elements '
plot 'finite_element_values.txt' using 1:15 with lines title 'Left0 Bspline Derivative ', \
'finite_element_values.txt' using 1:11 with lines title 'Left1 Bspline Derivative', \
'finite_element_values.txt' using 1:12 with lines title 'Left2 Bspline Derivative', \
'finite_element_values.txt' using 1:13 with lines title 'Middle Bspline Derivative'
set output
set term x11
reset
N=140
set term postscript eps enhanced color solid 22
set output "finite_elements/bspline.eps"
#set size .50,.250
set key off
#set yrange [-.1:.9]
set xrange [0:.4]
set xtics 0,.1,.4
set ytics 0,.2,.6
unset key
set xrange [0:4./N]
set xtics ("" 0, "" 1./N, "" 2./N, "" 3./N, "" 4./N)
set ytics (0)
#set title 'Cubic B-Spline
plot 'finite_element_values.txt' using 1:9  with lines 
set term x11
set output
reset
N=140
set term postscript eps enhanced color solid 22
set output "finite_elements/bspline_deriv.eps"
#set size .50,.250
set key off
#set yrange [-12:12]
set xrange [0:.4]
set xtics 0,.1,.4
set ytics -12,4,12
unset key
set xrange [0:4./N]
set xtics ("" 0, "" 1./N, "" 2./N, "" 3./N, "" 4./N)
set ytics (0)
#set title 'Cubic B-Spline Derivative'
plot 'finite_element_values.txt' using 1:13 with lines 
set output
set term x11
reset
N=140
set xrange [0:.2]
unset key
set xrange [0:2./N]
set xtics ("" 0, "" 1./N, "" 2./N, "" 3./N, "" 4./N)
set ytics (0)
#set title 'Comparison of derivatives of finite elements'
plot 'finite_element_values.txt' using 1:($31*2400) with lines title 'd/dr(Integral finite element)'
replot 'finite_element_values.txt' using 1:($30) with lines title 'rB_{/Symbol q}/B_z d/dr(tent)'
set term postscript enhanced eps color 22
set output 'finite_elements/left_mod_lin_deriv.eps'
replot
set term x11
set output
reset
exit
quit
N=140
set xrange [0:.2]
set key off
set title 'Integral finite element'
plot 'finite_element_values.txt' using 1:28 with lines 
set term postscript enhanced eps color 22
set output 'finite_elements/left_mod_lin.eps'
replot
set term x11
set output
