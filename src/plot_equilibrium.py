import Gnuplot, Numeric, sys
num = int(sys.argv[1])
print 'num = %d'%num
g = Gnuplot.Gnuplot()
g('set term postscript enhanced')
g('set output "equilib%d.ps"'%num)

g('set multiplot')
g('set size .5, 1/3.')

g('set origin .5, 2./3.')
g('set label "Profiles for equilib = %d" at screen .5,1.02 center font "Arial,24"'%num)
g('set xlabel "r"')
g('set title "pressure(r)"')
g('plot "equilib%d.txt" using 1:2 notitle with lines'%num)

g('set origin .5, 1./3.')
g('set title "B_z(r)"')
g('unset label')
g('plot "equilib%d.txt" using 1:3 notitle with lines'%num)

g('set origin 0, 1/3.')
g('set title "B_{/Symbol q}(r)"')
g('plot "equilib%d.txt" using 1:4 notitle with lines'%num)

g('set origin 0, 2/3.')
g('set title "{/Symbol r}(r)"')
g('plot "equilib%d.txt" using 1:5 notitle with lines'%num)

g('set origin 0., 0.')
g('set title "{/Symbol W}(r)*r"')
g('plot "equilib%d.txt" using 1:6 notitle with lines'%num)

g('set origin .5, 0')
g('set title "V_z(r)"')
g('plot "equilib%d.txt" using 1:7 notitle with lines'%num)

g('unset multiplot')
g('set term x11')
g('set output')
