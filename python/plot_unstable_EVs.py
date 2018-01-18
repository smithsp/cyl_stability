#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, Gnuplot, sys, math, multiplot, os
from get_EV import *
coord = sys.argv[1]
var_key = sys.argv[2]
unst_ind = int(sys.argv[3])
ylabels = []
d = []

num0 = int(sys.argv[4])
EV = get_EV(num0) 
minx = 0
maxx = EV.data['a']
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  minx = max(0,float(sys.argv[ind+1]))
  maxx = min(EV.data['a'],float(sys.argv[ind+2]))
tt = '{/Symbol x}_%s, %dst unstable'%(coord,unst_ind)
print_keys = EV.idata
for i in range(len(print_keys)):
  if print_keys[i]!=var_key:
    tt = tt+', %s=%g'%(print_keys[i],EV.data[print_keys[i]])
title=tt
amp = []
h_min = []
zero = []
var1 = []
rs = []
evr = []
evi =  []
for i in range(4,len(sys.argv)):
  try:
    num = int(sys.argv[i])
  except:
    break
  EV = get_EV(num) 
  N_fac = 5*max(EV.data['alpha'],1)
  if 'xr_rs' in sys.argv:
    rs1 = EV.data['rs']
    minx = max(0,rs1*.99)
    maxx = min(EV.data['a'],rs1*1.01)
  x = Numeric.arange(minx,maxx,EV.data['a']/(EV.data['ngrid']*N_fac))
  var1.append(EV.data[var_key])
  try:
    h_min.append(min(EV.grid[1:]-EV.grid[0:-1]))
  except AttributeError:
    print num
    continue
  #print num, h_min[-1]**-1
  print num, EV.data[var_key],
  rs.append(EV.data['rs'])
  unst = EV.get_unstable_evals()
  minr = -1e8
  maxr = 1e8
  if 'rwm' in sys.argv:
    minr = 0
    maxr = .1
  rwm_inds = Numeric.nonzero((unst.real>minr)*(unst.real<maxr))
  evec_num = EV.unst_inds[rwm_inds[unst_ind]]
  #print evec_num
  #print EV.evals[evec_num]
  (y,yl)  = EV.get_evec1(coord,evec_num,x)
  ereal = ((EV.evals[evec_num])).real
  eimag = ((EV.evals[evec_num])).imag
  evr.append(ereal)
  evi.append(eimag)
  ylabels.append( '%s=%g'%(var_key,EV.data[var_key]))
  if var_key=='ngrid':
    ylabels[-1] = 'N_{eff}=%g'%(int(EV.data['a']/min(EV.grid[1:]-EV.grid[:-1]))+1)
  dr =  Gnuplot.Data(x,y.real,title='Real',with='lines 3')
  di =  Gnuplot.Data(x,y.imag,title='Imag',with='lines 2')
  inds = Numeric.nonzero((EV.grid>=minx)*(EV.grid<=maxx))
  dg = Gnuplot.Data(Numeric.take(EV.grid,inds),0*inds,with='points 1',title='Grid')
  d.append([dg,dr,di])
  if EV.data['vcyl'] and eimag!=0:
    if EV.data['tw'] <= 0 and abs(eimag)>0:
      evec_num = int(round(evec_num+eimag/abs(eimag)))
    else:
      EV.imag = True
      del EV.sevecs  
    (y1,yl)  = EV.get_evec1(coord,evec_num,x)
    acoef = EV.get_EVrot(minr=minr,maxr=maxr)
    y2 = (acoef*y+1j*acoef*y1).real
    y3 = (acoef*y+1j*acoef*y1).imag
    #print acoef, y[-1],y1[-1],y2[-1],y3[-1]
    norm = EV.get_EVnorm()
    y2 = y2*norm
    y3 = y3*norm
    amp.append(1./((max(y2)+abs(min(y2)))/2.))
    try:
      x_zero = x[Numeric.nonzero(y2[1:]*y2[:-1]<0)[0]]
      ind = Numeric.argmin(Numeric.absolute(EV.grid-x_zero))
      zero.append(EV.grid[ind])
    except:
      zero.append(0)
    dr = Gnuplot.Data(x,y2,title='Re',with='lines 3')
    di =  Gnuplot.Data(x,y3,title='Im',with='lines 2')  
    inds = Numeric.nonzero((EV.grid>=minx)*(EV.grid<=maxx))
    dg = Gnuplot.Data(Numeric.take(EV.grid,inds),0*inds,with='points 1',title='Grid')
    d[-1]=[dg,dr,di]
    if 'nolegend' in sys.argv:
      for i in range(len(d[-1])):
        d[-1][i].set_option(title=None)
    EV.imag=False

fn = 'plots/%g_var_%s_evec%s_%dunstable.eps'%(num0,var_key,coord,unst_ind)
if 'xr' in sys.argv:
  fn = 'plots/%g_var_%s_evec%s_%dunstable_zoom.eps'%(num0,var_key,coord,unst_ind)
print 'output plotted to %s'%fn
if '1P' in sys.argv:#One plot overlayed
  g = Gnuplot.Gnuplot()
  g.itemlist = []
  for i in range(len(d)):
    g.itemlist.append(d[i][1])
  for i in range(len(g.itemlist)):
    g.itemlist[i].set_option(with='lines',title=ylabels[i])
#  g.refresh()
  g('set key right bottom Left reverse')
  g.xlabel('r')
  g.ylabel('Arb. units')
  g.hardcopy(fn,eps=True,fontsize=18,solid=False)
elif 'animate' in sys.argv:
  g = Gnuplot.Gnuplot()
  items = []
  inds = list(Numeric.argsort(h_min))
  inds.reverse()
  for i in range(len(d)):
    items.append(d[inds[i]][1])
    items[i].set_option(with='lines',title=ylabels[inds[i]])
    if i > 0:
      items[i].set_option(with='lines ls 3')
    g.itemlist = [items[i]]
    if i>0:
      g.itemlist.insert(0,items[0])
    g('set key right bottom Left reverse')
    g('set yrange [-50:50]')
    g.hardcopy(fn[:-4]+'animate-%d'%i+'.eps',solid=True,size=(4,3),linewidth=3,dashlength=50,eps=True,fontsize=24)
else:
  multiplot.multiplot(d,ylabels,fn)
#os.system('gv %s &'%fn)
if coord=='3':
  g = Gnuplot.Gnuplot()
  g._clear_queue()
  if var_key=='alpha' or var_key=='ngrid':
    inds = Numeric.argsort(h_min)
    h_min = list(Numeric.take(h_min,inds))
    amp = list(Numeric.take(amp,inds))
    evi = list(Numeric.take(evi,inds))
    g._add_to_queue([Gnuplot.Data(h_min,amp)])
    g._add_to_queue([Gnuplot.Data(h_min,evi,axes='x1y2')])
    m = (amp[1]-amp[0])/(h_min[1]-h_min[0])
    h_min.insert(0,0)
    amp.insert(0,0)
    evi.insert(0,0)
    h_min = Numeric.array(h_min)
    g._add_to_queue([Gnuplot.Data(h_min,m*h_min+amp[1]-m*h_min[1],with='lines 1')])
    print amp[1]-m*h_min[1]
    m = (evi[2]-evi[1])/(h_min[2]-h_min[1])
    g._add_to_queue([Gnuplot.Data(h_min,m*h_min+evi[1]-m*h_min[1],with='lines 1',axes='x1y2')])
    g.xlabel('h_{min}')
    g.ylabel('1/Amplitude (+)')
    g('set y2label "{/Symbol G (\264)}"')
    g('set y2tics auto')
    g('set y2range [-2e-6:1.4e-5]')
    g('set yrange [-.02:.14]')
    g('set ytics nomirror')
    g('set xrange [0:.0011]')
    g('set arrow 1 from .0007,.04 to .0009,.04 fill size .00005,30,90')
    g('set arrow 2 from .0006,.09 to .0004,.09 fill size .00005,30,90')
    #g.refresh()
    fn = 'plots/%g_var_%s_amp_vs_h.eps'%(num0,var_key)
    g.hardcopy(fn,eps=True,solid=False,color=False)
    print 'output plotted to %s'%fn
    g._clear_queue()
    
    
  elif var_key=='b':
    inds = Numeric.argsort(var1)
    var1 = Numeric.take(var1,inds)
    zero = Numeric.take(zero,inds)
    rs = Numeric.take(rs,inds)
    if 'fout' in sys.argv:
      f = open('b_zero_rs.out','w')
      for i in range(len(var1)):
        f.write('%g\t%g\t%g\n'%(var1[i],zero[i],rs[i]))
      f.close()
    #g.plot(Gnuplot.Data(var1,rs-zero))
    #sys.exit()
    bzero=1.6105
    ref_ind = int(len(var1)/2)
    fit = (var1-bzero)**.5*zero[ref_ind]/(var1[ref_ind]-bzero)**.5
    g._add_to_queue([Gnuplot.Data(var1,zero),Gnuplot.Data(var1,fit,with='lines',title='(b-%g)^{.5}*%g'%(bzero,zero[ref_ind]/(var1[ref_ind]-bzero)**.5))])
    g.refresh()
    fn = 'plots/%g_var_%s_vs_zero3.eps'%(num0,var_key)
    g.hardcopy(fn,eps=True)
    print 'output plotted to %s'%fn
    for i in range(len(var1)):
      print var1[i]
      print zero[i]
  

