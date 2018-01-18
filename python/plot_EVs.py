#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, Gnuplot, sys, math, multiplot, os
from get_EV import *
num = (sys.argv[1])
evec_num = int(sys.argv[2])
nc = int(sys.argv[3])
coord = (sys.argv[-nc:])
nc = len(coord)
N_fac = 10
EV = get_EV(num) 
ereal = ((EV.evals[evec_num])).real
eimag = ((EV.evals[evec_num])).imag
minx = 0
maxx = EV.data['a']
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  minx = max(0,float(sys.argv[ind+1]))
  maxx = min(EV.data['a'],float(sys.argv[ind+2]))
x = Numeric.arange(minx,maxx,EV.data['a']/(EV.data['ngrid']*N_fac))
print_keys = EV.idata
tt = '{/Symbol w}=%g+%gi, '%(ereal,eimag)
for i in range(len(print_keys)):
  tt = tt+'%s=%g, '%(print_keys[i],EV.data[print_keys[i]])
title=tt[:-2] 
ylabels = []
d = []
for i in range(nc):
  #ylabels.append('{/Symbol x}_%s'%coord[i])
  (y,yl)  = EV.get_evec1(coord[i],evec_num,x)
  ylabels.append(yl)
  dr =  Gnuplot.Data(x,y.real,title='Real',with='lines')
  di =  Gnuplot.Data(x,y.imag,title='Imag',with='lines')
  inds = Numeric.nonzero((EV.grid>=minx)*(EV.grid<=maxx))
  dg = Gnuplot.Data(Numeric.take(EV.grid,inds),0*inds,with='points 1',title='Grid')
  d.append([dr,di])
  if EV.data['vcyl'] and eimag!=0:
    if EV.data['tw'] <= 0 and abs(eimag)>0:
      #print 'here'
      evec_num = int(round(evec_num+eimag/abs(eimag)))
    else:
      EV.imag = True
      del EV.sevecs  
    (y1,yl)  = EV.get_evec1(coord[i],evec_num,x)
    EV.imag=False
    acoef = EV.get_EVrot(minr=-1e7,maxr=1e7)#EV.data['k']*EV.data['Vz0']
    y2 = (acoef*y+1j*acoef*y1).real
    y3 = (acoef*y+1j*acoef*y1).imag
    print acoef, y[-1],y1[-1],y2[-1],y3[-1]
    norm = EV.get_EVnorm()
    y2 = y2*norm
    y3 = y3*norm
    dr = Gnuplot.Data(x,y2,title='Re',with='lines 3')
    di =  Gnuplot.Data(x,y3,title='Im',with='lines 2')  
    inds = Numeric.nonzero((EV.grid>=minx)*(EV.grid<=maxx))
    dg = Gnuplot.Data(Numeric.take(EV.grid,inds),0*inds,with='points 1',title='Grid')
    d[-1]=[dr,di]

evec_num0 = evec_num
fn = 'plots/%s.evecs_%d.eps'%(num,evec_num0)
if 'xr' in sys.argv:
  fn = 'plots/%s.evecs_%d_zoom.eps'%(num,evec_num0)
print 'output to %s'%fn
d.reverse()
ylabels.reverse()
multiplot.multiplot(d,ylabels,fn,title)
os.system('gv %s &'%fn)
sys.exit()
