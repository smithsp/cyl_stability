#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, Gnuplot, sys, sets
from get_EV import *
ref = int(sys.argv[1])
d = get_EV(ref,output_folder='output_vcyl')

len_branch = len(d.evals)/6-1
xs = 1.5
ys = 2.5*math.log(len_branch)
dy = ys/len_branch
dx = xs/3.
inds = list(Numeric.argsort(d.evals.real))
inds.reverse()
x = Numeric.arange(0,d.data['a'],1./(d.data['ngrid']*10))
#print inds
#print Numeric.take(d.evals.real,inds)
fns = ['fast', 'alfven', 'slow']
for i in range(3):
  g = Gnuplot.Gnuplot()
  g('set output "plots/%d_all_%s.eps"'%(ref,fns[i]))
  g('set term postscript enhanced eps solid 8')
  g('set size %g,%g'%(xs,ys))
  g('set multiplot')
  g('set size %g,%g'%(dx,dy))
  if len_branch>50:
    g('unset xtics')
    #g('unset ytics')
    #g('set yrange [0:1]')
    g('set xrange [0:%g]'%d.data['a'])
  
  for j in range(len_branch-1,-1,-1):
    y0 = j*dy
    ind = inds[i*len_branch+j+6]
    g('set ylabel "{/Symbol w}=%g"'%(d.evals[ind].real))
    for k in range(3):
      if j==len_branch-1:
        g('set title "{/Symbol x}_%d"'%(k+1))
      else:
        g('unset title')
      g('set origin %g,%g'%(dx*k,y0))
      if k>0:
        g('unset ylabel')
      evecr = (d.get_evec1('%d'%(k+1),ind,x))[0].real
      if max(abs(evecr))==abs(min(evecr)):
        evecr = evecr*(-1)
      zinds = range(len(x))
      if len_branch>50:
        zinds = Numeric.nonzero(abs(evecr[1:]-evecr[:-1])>(max(evecr)-min(evecr))/len_branch)
        zinds = list(sets.Set(list(zinds)+list(zinds+1)+list(zinds-1))
        zinds.sort()
      g.plot(Gnuplot.Data(Numeric.take(x,zinds),Numeric.take(evecr,zinds),with='lines'))
  g('unset multiplot')    
  g('set term x11')
  g('set output')
sys.exit()
