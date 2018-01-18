#! /usr/pppl/python/2.3.5x/bin/python
import Gnuplot, Numeric, sys, os
num = int(sys.argv[1])
try:
  equilib_path = 'equilibria_vcyl'
  filename = 'plots/%d_equilibrium.eps'%(num)
  infile = '%s/%d.txt'%(equilib_path,num)
  f= file(infile,'r')
  s = f.readline()
  f.close()
except IOError:
  equilib_path = 'equilibria_cyl'
  filename = 'plots/%d_equilibrium.eps'%(num)
  infile = '%s/%d.txt'%(equilib_path,num)
  f= file(infile,'r')
  s = f.readline()
  f.close()
  #print infile + ' does not exist'
  #sys.exit()
g = Gnuplot.Gnuplot()
g('set term postscript enhanced')
g('set output "%s"'%filename)

g('set multiplot')
g('set size .5, 1/3.')

g('set origin .5, 2./3.')
g('set label "Equilibrium profiles for run %d" at screen .5,.98 center font "Arial,16"'%(num))
g('set xlabel "r"')
g('set title "pressure(r)"')
g('plot "%s/%d.txt" using 1:2 notitle with lines'%(equilib_path,num))

g('set origin .5, 1./3.')
g('set title "B_z(r)"')
g('unset label')
g('plot "%s/%d.txt" using 1:3 notitle with lines'%(equilib_path,num))

g('set origin 0, 1/3.')
g('set title "B_{/Symbol q}(r)"')
g('plot "%s/%d.txt" using 1:4 notitle with lines'%(equilib_path,num))

g('set origin 0, 2/3.')
g('set title "q(r)"')
g('plot "%s/%d.txt" using 1:8 notitle with lines'%(equilib_path,num))
#g('set title "{/Symbol r}(r)"')
#g('plot "%s/%d.txt" using 1:5 notitle with lines'%(equilib_path,num))

g('set origin 0., 0.')
g('set title "V_{/Symbol q}(r)"')
g('plot "%s/%d.txt" using 1:6 notitle with lines'%(equilib_path,num))

g('set origin .5, 0')
g('set title "V_z(r)"')
g('plot "%s/%d.txt" using 1:7 notitle with lines'%(equilib_path,num))


g('unset multiplot')
g('set term x11')
g('set output')
os.system('gv %s &'%(filename))
