#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, Gnuplot, sys, math, multiplot
from get_EV import *
num = (sys.argv[1])
coord = sys.argv[2]
evec_num = sys.argv[3:]
for i in range(len(evec_num)):
  evec_num[i] = int(evec_num[i])
nn = len(evec_num)
N_fac = 10
EV = get_EV(num) 
acoef = EV.get_EVrot(minr=.04,maxr=EV.data['k']*EV.data['Vz0'])
norm = EV.get_EVnorm()
ylabels = []
d = []
x = Numeric.arange(0,EV.data['a'],EV.data['a']/(EV.data['ngrid']*N_fac))
tt = '{/Symbol x}_%s, '%coord
print_keys = EV.idata
for i in range(len(print_keys)):
  tt = tt+'%s=%g, '%(print_keys[i],EV.data[print_keys[i]])
title=tt[:-2] 
for i in range(nn):
  EV.imag=False
  try:
    del EV.sevecs
  except:
    pass
  (y,yl)  = EV.get_evec1(coord,evec_num[i],x)
  ereal = ((EV.evals[evec_num[i]])).real
  eimag = ((EV.evals[evec_num[i]])).imag
  ylabels.append( '{/Symbol w}=%g+%gi'%(ereal,eimag))
  dr =  Gnuplot.Data(x,y.real,title='Real',with='lines')
  di =  Gnuplot.Data(x,y.imag,title='Imag',with='lines')
  dg = Gnuplot.Data(EV.grid,[0]*len(EV.grid),with='points 1',title='Grid')
  d.append([dr,di,dg])
  if EV.data['vcyl'] and eimag!=0:
    if EV.data['tw'] <= 0 and abs(eimag)>0:
      evec_num[i] = int(round(evec_num[i]+eimag/abs(eimag)))
    else:
      EV.imag = True
      del EV.sevecs  
    (y1,yl)  = EV.get_evec1(coord,evec_num[i],x)
    y2 = (acoef*y+1j*acoef*y1).real
    y3 = (acoef*y+1j*acoef*y1).imag
    y2 = y2*norm
    y3 = y3*norm
    dr = Gnuplot.Data(x,y2,title='Re',with='lines 3')
    di =  Gnuplot.Data(x,y3,title='Im',with='lines 2')    
    dg = Gnuplot.Data(EV.grid,[0]*len(EV.grid),with='points 1',title='Grid')
    d[-1]=[dr,di,dg]

fn = 'plots/%s.evec%s_multi.eps'%(num,coord)
print 'Output to %s'%fn
multiplot.multiplot(d,ylabels,fn,title)
