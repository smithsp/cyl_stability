#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys, multiplot
ref = sys.argv[1]
var_key = sys.argv[2]
scale_fac = 10**4
try:
  scale_fac = float(sys.argv[3])
except:
  print 'scale_fac=',scale_fac
sfr = scale_fac
sfi = scale_fac
if 'sfi' in sys.argv:
  ind = sys.argv.index('sfi')
  sfi = float(sys.argv[ind+1])
if 'sfr' in sys.argv:
  ind = sys.argv.index('sfr')
  sfr = float(sys.argv[ind+1])
linearx = False
lineary = False
if 'linearx' in sys.argv:
  linearx = True
if 'lineary' in sys.argv:
  lineary = True
d1 = get_EV(ref)
dk = dict()
print_keys = d1.idata 
dk_keys = ['t_assemb','num','t_solve',var_key,'epsilo','NN','vcyl','rs','nu']
if var_key=='ngrid' or var_key=='alpha':
  dk_keys+=['nu']
for key in d1.data.keys():
  if key not in dk_keys:
    dk[key] = d1.data[key]
matching_files = []
var_list= []
p0_list = []
eps_list = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
kV_fac = 1.0  
g.xlabel('[Re({/Symbol w})-kV_{z0}]/{/Symbol w}_A')
if 'kVf0' in sys.argv:
  kV_fac = 0.0
  g.xlabel('Re({/Symbol w})/{/Symbol w}_A')
g.ylabel('{/Symbol G}/{/Symbol w}_A')
#g('set size ratio -1')
ttl = var_key
if ttl=='tw':
  ttl = '{/Symbol t}_w'
#print 'Free numbers: ',
#for i in range(1,271,1)+range(696,710):
maxr = []
maxi = []
minnum = 1
maxnum = 2500
if 'minn' in sys.argv:
  minnum = int(sys.argv[sys.argv.index('minn')+1])
if 'maxn' in sys.argv:
  maxnum = int(sys.argv[sys.argv.index('maxn')+1])
dn = 1
if 'dn' in sys.argv:
  ind = sys.argv.index('dn')
  dn = int(sys.argv[ind+1])
if 'SN' in sys.argv:#Specify the shot numbers
  nums = []
  ind = sys.argv.index('SN')
  for num in sys.argv[ind+1:]:
    try:
      nums.append(int(num))
    except:
      break
else:
  nums = range(minnum,maxnum,dn)
for i in nums:
  try:
    #print i,
    d2 = get_EV(i)
    data = d2.data
  except IOError:
    #print '%d, '%i,
    continue
  match = True
  for key in dk.keys(): 
    if var_key=='alpha':
      if key=='rs':
        continue
    if data[key] != dk[key]:
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:
    matching_files.append(data['num'])
    evals = d2.get_numeric_evals() #get_unstable_evals()
    wA = d2.get_wA()
    evalsr = evals.real  #Only for initialization to the right size
    evalsi = evals.imag/wA  #Only for initialization to the right size
    for j in range(len(evals)):
      if data['equilibr']<11 or data['equilibr']==13:
        if linearx:
          evalsr[j] = (evals[j].real-kV_fac*data['k']*data['Vz0'])/wA
        else:
          evalsr[j] = cmath.asinh(sfr*(evals[j].real-kV_fac*data['k']*data['Vz0'])/wA).real
      elif data['equilibr']==12:
        kV = data['k']*data['Vp0']/data['Bt0']*data['Bz0']*data['a']+data['m']*data['Vp0']*data['a']
        kV = kV*kV_fac
        if linearx:
          evalsr[j] = (evals[j].real-kV)/wA
        else:
          evalsr[j] = cmath.asinh(sfr*(evals[j].real-kV)/wA).real
      if not lineary:
        evalsi[j] = cmath.asinh(sfi*(evals[j].imag)).real
    ind = Numeric.argmax(evalsi)
    maxr.append(evalsr[ind])
    maxi.append(evalsi[ind])
    ttl1 = '%s=%2.4g'%(ttl,data[var_key])
    if var_key=='tw':
      if data['tw']==-1:
        ttl1='Ideal Wall'
      elif data['tw']==0:
        ttl1='No Wall'
      else:
        ttl1 = '%s=10^{%g}'%(ttl,(cmath.log10(data[var_key]).real))
    elif var_key=='alpha':
      try:
        shft_grid=d2.grid[1:]-d2.grid[:-1]
        ttl1 = 'N_l: %d; N_h: %d'%(int(data['a']/max(shft_grid)),int(data['a']/min(shft_grid)))
      except:
        pass
    try:
      if 'minmax' in sys.argv:
        inds = d2.get_unst_inds()
        etemp = Numeric.array([d2.evals[inds[0]],d2.evals[inds[-1]],d2.evals[inds[-2]]])
        evalsr = etemp.real
        evalsi = etemp.imag
      g._add_to_queue([Gnuplot.Data(evalsr,evalsi,with='points ',title=ttl1)])#)])#
      #if i == int(ref):
        #g.itemlist[-1].set_option(with='points ps 3')
    except (ValueError):
      continue
    var_list.append(data[var_key])

for i in range(len(matching_files)):
  print matching_files[i],':',var_list[i],';  ',
print ''
inds = Numeric.argsort(var_list)
g.itemlist = list(Numeric.take(g.itemlist,inds))
if min(var_list)==-1:
  g.itemlist.append(g.itemlist[0])
  g.itemlist = g.itemlist[1:]

#inds = Numeric.argsort(maxr)
maxr = Numeric.take(maxr,inds)
maxi = Numeric.take(maxi,inds)

if 'max_line' in sys.argv:
  g._add_to_queue([Gnuplot.Data(maxr,maxi,with='lines')])
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1.data[print_keys[i]])
if 'notitle' not in sys.argv:
  g.title(tt[:-2])

set_ytics = 'set ytics ('
set_xtics = 'set xtics ('
for i in range(2):
  for j in range(-int(cmath.log10(sfr).real)-2,13,1):
    scale = (-1)**i*10**j
    if j%2 != 1:
      set_xtics  = set_xtics + '"" %g 1, '%(cmath.asinh(scale*sfr).real)
    else:
      set_xtics  = set_xtics + '"%g" %g, '%(scale,cmath.asinh(scale*sfr).real)
  for j in range(-int(cmath.log10(sfi).real)+0,13,1):
    scale = (-1)**i*10**j
    if j%2 != 0:
      set_ytics  = set_ytics + '"" %g 1, '%(cmath.asinh(scale*sfi).real)
    else:
      set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*sfi).real)
set_ytics = set_ytics + '"0" 0 )'
set_xtics = set_xtics + '"0" 0 )'
if not lineary:
  g(set_ytics)
if not linearx:
  g(set_xtics)
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  yr1 = float(sys.argv[ind+1])
  yr2 = float(sys.argv[ind+2])
  if lineary:
    g('set yrange [%g:%g]'%(yr1,yr2))
  else:
    g('set yrange [%g:%g]'%(cmath.asinh(sfi*yr1).real,cmath.asinh(sfi*yr2).real))
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  xr1 = float(sys.argv[ind+1])
  xr2 = float(sys.argv[ind+2])
  if linearx:
    g('set xrange [%g:%g]'%(xr1,xr2))
  else:
    g('set xrange [%g:%g]'%(cmath.asinh(sfr*xr1).real,cmath.asinh(sfr*xr2).real))
if 'animate' in sys.argv:
  #print 'here'
  items = g.itemlist[:]
  g('set key top left reverse Left')
  if var_key=='b':
    items.reverse()
  for i in range(len(items)):
    item = items[i]
    g.itemlist = [item]
    g('set rmargin 10')
    if var_key=='Vz0':
      g.itemlist.append(Gnuplot.Data([0,0],[cmath.asinh(sfi*-9e-3).real,cmath.asinh(sfi*9e-3).real],with='lines ls 0',title=None))
    if i > 0:
      g.itemlist.insert(0,items[0])
      g.itemlist[1].set_option(with='points 7')
    g.hardcopy('plots/%d_var_%s_complex_animate-%d.eps'%(min(matching_files),var_key,i),eps=True)
  g.itemlist = items
if 'multi' in sys.argv:
  fn = "plots/%d_var_%s_complex_separate.eps"%(min(matching_files),var_key)
  items = g.itemlist
  itemlist = []
  for i in range(len(g.itemlist)):
    itemlist.append([items[i],Gnuplot.Data([0,0],[cmath.asinh(sfi*-9e-3).real,cmath.asinh(sfi*9e-3).real],with='lines ls 0',title=None)])
  
  #multiplot.multiplot(itemlist,ylabels=['{/Symbol G}']*len(itemlist),fn=fn,xsize=3,ysize=5)
  xs = 3
  ys = 5
  yh = 1./len(itemlist)
  g('unset xlabel')
  g('unset key')
  g('set xtics ("" -.2, "" -.1, "" 0, "" 0.1, "" 0.2, "" 0.3, "" 0.4)')
  g('set bmargin 0')
  g('set rmargin 0')
  g('set lmargin 11')
  g('set tmargin 0') 
  g('set term postscript eps enhanced lw 3 size 3,5 18')
  g('set output "%s"'%(fn))
  g('set multiplot layout %d,1'%(len(itemlist)+1))
  alphabet='abcdefghijklmnop'
  for i in range(len(itemlist)):
    g.itemlist = (itemlist[i])
    g('set label 1 "%s" at -0.2,%g'%(alphabet[i],cmath.asinh(sfi*1e-3).real))
    g.refresh()
  #g('show margin')
  g('unset label 1')
  #g('set size 1,0.001')
  #g('set origin 0,.125')
  g('set bmargin 5')
  g('set tmargin 0')
  g('set xlabel "Re({/Symbol w})"')
  g('set xtics auto')
  g('unset ylabel')
  g('unset ytics')
  g.itemlist = [g.itemlist[1]]
  g.refresh()
  g('unset multiplot')
  g('set output')
  g('set term x11')
  sys.exit()

g('set key center outside  Left reverse') 
if 'IK' in sys.argv:
  if 'kVf0' in sys.argv:
    g('set key top inside right Left reverse')
  else:
    g('set key top inside left Left reverse') 
for i in range(1,13):
  lbl = 'lb%d'%i
  if lbl in sys.argv:
    ind = sys.argv.index(lbl)
    x1 = float(sys.argv[ind+2])
    y1 = float(sys.argv[ind+3])
    if not linearx:
      x1 = cmath.asinh(sfr*x1).real
    if not lineary:
      y1 = cmath.asinh(sfi*y1).real
    g('set label "%s" at %g,%g center'%(sys.argv[ind+1],x1,y1))
for i in range(1,13):
  lbl = 'rlb%d'%i
  if lbl in sys.argv:
    ind = sys.argv.index(lbl)
    x1 = float(sys.argv[ind+2])
    y1 = float(sys.argv[ind+3])
    if not linearx:
      x1 = cmath.asinh(sfr*x1).real
    if not lineary:
      y1 = cmath.asinh(sfi*y1).real
    g('set label "%s" at %g,%g center rotate by 90'%(sys.argv[ind+1],x1,y1))
if ('xr' in sys.argv) or ('yr' in sys.argv):
  fn = 'plots/%d_var_%s_shift_complex_zoom.eps'%(min(matching_files),var_key)
else:
  fn = 'plots/%d_var_%s_shift_complex.eps'%(min(matching_files),var_key)
print 'Plotted to %s'%(fn)
#g.refresh()
#g('unset key')
g.hardcopy(fn,fontsize=18,eps=True,size=(5,3.5))#,eps=True,size=(5,3.5))
  
