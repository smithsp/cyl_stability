#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, Gnuplot, sys, math, multiplot
from get_EV import *
coord = sys.argv[1]
var_key = sys.argv[2]
ylabels = []
d = []
N_fac = 10
num0 = int(sys.argv[3])
EV = get_EV(num0) 
x = Numeric.arange(0,EV.data['a'],EV.data['a']/(EV.data['ngrid']*N_fac))
tt = '{/Symbol x}_%s'%coord
print_keys = EV.idata
for i in range(len(print_keys)):
  if print_keys[i]!=var_key:
    tt = tt+', %s=%g'%(print_keys[i],EV.data[print_keys[i]])
title=tt
for i in range(3,len(sys.argv),2):
  num = int(sys.argv[i])
  evec_num = int(sys.argv[i+1])
  EV = get_EV(num) 
  (y,yl)  = EV.get_evec1(coord,evec_num,x)
  ereal = ((EV.evals[evec_num])).real
  eimag = ((EV.evals[evec_num])).imag
  ylabels.append( '%s=%g'%(var_key,EV.data[var_key]))
  dr =  Gnuplot.Data(x,y.real,title='Real',with='lines')
  di =  Gnuplot.Data(x,y.imag,title='Imag',with='lines')
  dg = Gnuplot.Data(EV.grid,[0]*len(EV.grid),with='points 1',title='Grid')
  d.append([dr,di,dg])
  if EV.data['vcyl'] and eimag!=0:
    if EV.data['tw'] <= 0 and abs(eimag)>0:
      evec_num = int(round(evec_num+eimag/abs(eimag)))
    else:
      EV.imag = True
      del EV.sevecs  
    (y1,yl)  = EV.get_evec1(coord,evec_num,x)
    acoef = 1.0
    ind = -1
    if y[ind] != 0:
      acoef = 1.-1j*y1[ind]/y[ind]  
    y4 = (y+1j*y1).real
    y5 = (y+1j*y1).imag
    y = y4
    y1 = y5
    y2 = (acoef*y+1j*acoef*y1).real
    y3 = (acoef*y+1j*acoef*y1).imag
    dr = Gnuplot.Data(x,y2,title='Re',with='lines')
    di =  Gnuplot.Data(x,y3,title='Im',with='lines')  
    dg = Gnuplot.Data(EV.grid,[0]*len(EV.grid),with='points 1',title='Grid')
    d[-1]=[dr,di,dg]

fn = 'plots/%g_var_%s_evec%s.eps'%(num0,var_key,coord)
print 'output plotted to %s'%fn
multiplot.multiplot(d,ylabels,fn,title)
