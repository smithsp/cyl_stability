#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys
ref = sys.argv[1]
var_key = sys.argv[2]
scale_fac = 10**10
try:
  scale_fac = float(sys.argv[3])
except:
  print 'scale_fac=',scale_fac
d1 = get_EV(ref)
dk = dict()
print_keys = d1.idata #['Vz0','P1','alpha','ngrid','tw','b']
dk_keys = ['t_assemb','num','t_solve',var_key,'epsilo','NN','vcyl']
for key in d1.data.keys():
  if key not in dk_keys:
    dk[key] = d1.data[key]
matching_files = []
var_list= []
p0_list = []
eps_list = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g.xlabel('Re({/Symbol w})')
g.ylabel('{/Symbol G}')
g('set size ratio -1')
g('set key center outside  Left reverse')
ttl = var_key
if ttl=='tw':
  ttl = '{/Symbol t}_w'
#print 'Free numbers: ',
for i in range(101,103)+range(265,272,1)+range(696,700):
#minnum = 1
#maxnum = 1400
#if 'minn' in sys.argv:
#  minnum = int(sys.argv[sys.argv.index('minn')+1])
#if 'maxn' in sys.argv:
#  maxnum = int(sys.argv[sys.argv.index('maxn')+1])
#for i in range(minnum,maxnum,1):
  try:
    d2 = get_EV(i)
    data = d2.data
  except IOError:
    #print '%d, '%i,
    continue
  match = True
  for key in dk.keys(): 
    if key == 'm':
      if abs(data[key]) == abs(dk[key]):
        continue   
    if data[key] != dk[key]:
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    matching_files.append(data['num'])
    evals = d2.get_numeric_evals() #get_unstable_evals()
    evalsr = evals.real  #Only for initialization to the right size
    evalsi = evals.imag  #Only for initialization to the right size
    for j in range(len(evals)):
      evalsr[j] = cmath.asinh(scale_fac*(evals[j].real-0*d2.data['k']*d2.data['Vz0'])).real
      evalsi[j] = cmath.asinh(scale_fac*(evals[j].imag)).real
    ttl1 = '%s=%2.4g'%(ttl,data[var_key])
    if var_key=='tw':
      if data['tw']==-1:
        ttl1='Ideal Wall'
      elif data['tw']==0:
        ttl1='No Wall'
      else:
        ttl1 = '%s=10^{%g}'%(ttl,(cmath.log10(data[var_key]).real))
    try:
      g._add_to_queue([Gnuplot.Data(evalsr,evalsi,with='points ps 2',title=ttl1)])
      if i == int(ref):
        g.itemlist[-1].set_option(with='points ps 3')
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
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1.data[print_keys[i]])
g.title(tt[:-2])
set_ytics = 'set ytics ('
set_xtics = 'set xtics ('
for i in range(2):
  for j in range(-int(cmath.log10(scale_fac).real)+0,13,1):
    scale = (-1)**i*10**j
    if j%2 != 0:
      set_ytics  = set_ytics + '"" %g 1, '%(cmath.asinh(scale*scale_fac).real)
      set_xtics  = set_xtics + '"" %g 1, '%(cmath.asinh(scale*scale_fac).real)
    else:
      set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)
      set_xtics  = set_xtics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)
set_ytics = set_ytics + '"0" 0 )'
set_xtics = set_xtics + '"0" 0 )'
g(set_ytics)
g(set_xtics)
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  g('set yrange [%g:%g]'%(cmath.asinh(scale_fac*float(sys.argv[ind+1])).real,cmath.asinh(scale_fac*float(sys.argv[ind+2])).real))
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g('set xrange [%g:%g]'%(cmath.asinh(scale_fac*float(sys.argv[ind+1])).real,cmath.asinh(scale_fac*float(sys.argv[ind+2])).real))
for i in range(1,13):
  lbl = 'lb%d'%i
  if lbl in sys.argv:
    ind = sys.argv.index(lbl)
    x1 = cmath.asinh(scale_fac*float(sys.argv[ind+2])).real
    y1 = cmath.asinh(scale_fac*float(sys.argv[ind+3])).real
    g('set label "%s" at %g,%g center'%(sys.argv[ind+1],x1,y1))
for i in range(1,13):
  lbl = 'rlb%d'%i
  if lbl in sys.argv:
    ind = sys.argv.index(lbl)
    x1 = cmath.asinh(scale_fac*float(sys.argv[ind+2])).real
    y1 = cmath.asinh(scale_fac*float(sys.argv[ind+3])).real
    g('set label "%s" at %g,%g center rotate by 90'%(sys.argv[ind+1],x1,y1))
fn = 'plots/%d_var_%s_complex.eps'%(min(matching_files),var_key)
print 'Plotted to %s'%(fn)
#g.refresh()
#g('unset key')
g.hardcopy(fn,fontsize=24,eps=True,size=(1.5,1.5))
  
