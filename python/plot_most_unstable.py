#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys, sets
ref = int(sys.argv[1])
var_key = sys.argv[2]
var_key2 = sys.argv[3]
nq = False
Vr = False
if var_key == 'nq':
  nq = True
  if 'z' in sys.argv:
    var_key = 'Bz0'
  else:
    var_key = 'Bt0'
elif var_key == 'Vr':
  Vr = True
  ind = sys.argv.index('r0')
  r0 = float(sys.argv[ind+1])
  if 'z' in sys.argv:
    var_key  = 'Vz0'
  else:
    var_key = 'epsVz'
d1 = get_EV(ref)
print_keys = d1.idata 
if 'pk' in sys.argv:
  print_keys = []
  ind = sys.argv.index('pk')
  for pk in sys.argv[ind+1:]:
    if pk in d1.data.keys():
      print_keys.append(pk)
    else:
      print pk, 'not in', d1.data.keys()
      break
d1 = d1.data
dk = dict()
dk_keys = ['t_assemb','num','t_solve',var_key,var_key2,'epsilo','NN','vcyl']
for key in d1.keys():
  if key not in dk_keys:
    dk[key] = d1[key]
matching_files = []
var_list= []
var2_list = []

g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
g1.xlabel(var_key)
if nq:
  g1.xlabel('nq')
elif Vr:
  g1.xlabel('k V_z(r_0=%g)'%(r0))
if var_key == 'Lend0':
  g1('set xrange [-1:2]')
if 'ly' in sys.argv:
  g1('set log y')
negty = ''
nfacy = 1 
if 'negy' in sys.argv:
  negty = '-'
  nfacy = -1
g1.ylabel('%s{/Symbol G}'%negty)
if 'tw-1' in sys.argv:
  g1.ylabel('%s{/Symbol G}{/Symbol t}_w'%negty)
scale_fac = 10**5
minevalsr = 0
maxevalsr = 0
minn = 1
maxn = 2100
if 'minn' in sys.argv:
  ind = sys.argv.index('minn')
  minn = int(sys.argv[ind+1])
if 'maxn' in sys.argv:
  ind = sys.argv.index('maxn')
  maxn = int(sys.argv[ind+1])
for i in range(minn,maxn):
  try:
    d2 = get_EV(i)
    data = d2.data
  except:
    continue
  match = True
  for key in dk.keys():   
    if data[key] != dk[key]:
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    matching_files.append(data['num'])
    #nq = data['nz']*data['Bz0']/(20*cmath.pi**2*data['Bt0'])
    if nq:
      var_list.append((data['k']*data['Bz0']/data['Bt0']*data['a']))
    elif Vr:
      var_list.append(data['k']*d2.Vz([r0])[0])
    else:
      var_list.append(data[var_key])
    var2_list.append(data[var_key2])
uniq_var2 = list(sets.Set(var2_list))
uniq_var2.sort()
k = 1
if 'var2r' in sys.argv:
  ind = sys.argv.index('var2r')
  var2_list2 = list(sys.argv[ind+1:])
  for i in range(len(var2_list2)):
    try:
      var2_list2[i] = float(var2_list2[i])
    except:
      break
  var2_list2.sort()
  uniq_var2 = var2_list2
if var_key2 == 'tw':
  try:
    ind = uniq_var2.index(-1)
    uniq_var2.append(uniq_var2[ind])
    uniq_var2.pop(ind)
  except:
    pass
ttl = var_key2
if var_key2=='tw':
  ttl = '{/Symbol t}_w'
#print 'uniq_var2 = ',uniq_var2
for v2 in uniq_var2:
  maxunstable = []
  maxvarlist = []
  for i in range(len(matching_files)):
    if var2_list[i]==v2:
      d2 = get_EV(matching_files[i])
      nfacy = 1
      if 'tw-1' in sys.argv and var_key=='tw':
        nfacy = d2.data['tw']
      evals = d2.evals
      if 'wr' in sys.argv: #This is for the range in Re(w)
        ind = sys.argv.index('wr')
        wr1 = float(sys.argv[ind+1])
        wr2 = float(sys.argv[ind+2])
        inds = Numeric.nonzero((evals.real>wr1)*(evals.real<wr2))
        evals = Numeric.take(evals,inds)
      maxunstable.append(max(evals.imag)*nfacy)
      maxvarlist.append(var_list[i])
  if len(maxvarlist)==0:
    print 'Warning!!! The value of %s=%g is not valid'%(var_key2,v2)
    continue
  inds = Numeric.argsort(maxvarlist)
  maxvarlist = Numeric.take(maxvarlist,inds)
  maxunstable = Numeric.take(maxunstable,inds)
  ttl1 = '%s=%g'%(ttl,v2)
  if var_key2=='tw':
    if v2==-1:
      ttl1 = "%s={/Symbol \245}"%(ttl)
    elif v2==0:
      ttl1 = "%s=0"%(ttl)
    else:
      ttl1 = '%s=10^%d'%(ttl,int((cmath.log10(v2)).real))
  elif var_key2=='Vz0':
      if v2==0:
        ttl1 = 'V_{z0}=0'
      else:
        wA = abs(d2.get_wA())
        ttl1 = 'V_{z0}=%2.2g v_{A}'%abs(v2/(wA/data['k']))
      g1('set key box spacing 1.5')
  g1._add_to_queue([Gnuplot.Data(maxvarlist,maxunstable,with='linespoints ps 2',title=ttl1)])
for i in range(len(matching_files)):
  print matching_files[i],':',var_list[i],';  ',
print ''

set_ytics = 'set ytics ('
for i in range(2):
  for j in range(-5,13,1):
    scale = (-1)**i*10**j
    if j%2 != 0:
      set_ytics  = set_ytics + '"" %g 1, '%(cmath.asinh(scale*scale_fac).real)
    else:
      set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)
set_ytics = set_ytics + '"0" 0 )'
if not nq:
  maxfac = 1+0.1*max(var_list)/abs(max(var_list))
  if abs(min(var_list)/max(var_list))<1/(6.0e1):
    g1('set log x')
  else:
    if min(var_list)!=0:
      minfac = 1-0.1*min(var_list)/abs(min(var_list))
      g1('set xrange [%g:%g]'%(min(var_list)*minfac,max(var_list)*maxfac))
    else:
      g1('set xrange [:%g]'%(max(var_list)*maxfac))
g1('set key Left reverse bottom')
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1[print_keys[i]])
g1.title(tt[:-2])
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g1('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  g1('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g1('set ytics auto')
for i in range(1,13):
  lbl = 'lb%d'%i
  if lbl in sys.argv:
    ind = sys.argv.index(lbl)
    x1 = float(sys.argv[ind+2])
    y1 = float(sys.argv[ind+3])
    g1('set label "%s" at %g,%g left'%(sys.argv[ind+1],x1,y1))
#g1('set size .5,1')
#g1.refresh()
if var_key=='tw':
  g1.xlabel('{/Symbol t}_w')
  if 'tw-1' not in sys.argv:
    g1._add_to_queue([Gnuplot.Data(Numeric.sort(Numeric.array(var_list)),2*abs(dk['m'])/Numeric.sort(Numeric.array(var_list)),with='lines')])
if nq:
  fn = 'plots/%d_var_%s.eps'%(min(matching_files),'nq')
elif Vr:
  fn = 'plots/%d_var_Vr_%s_maxg.eps'%(min(matching_files),var_key2)
else:
  fn = 'plots/%d_var_%s_%s_maxg.eps'%(min(matching_files),var_key,var_key2)
fn = '%s_imag.eps'%fn[:-4]
print 'Plotted to %s'%(fn)
if 'kp' in sys.argv:
  ind = sys.argv.index('kp')
  g1('set key %s,%s'%(sys.argv[ind+1],sys.argv[ind+2]))
minv = min(var_list)*.99
maxv = max(var_list)*1.01
x0 = Numeric.array([minv,(2*minv+maxv)/3.,(minv+2*maxv)/3.,maxv])
#x0 = Numeric.array(var_list)
y0 = x0*0
g1._add_to_queue([Gnuplot.Data(x0,y0,with='lines 7')])
if 'notitle' in sys.argv:
  g1('unset title')
g1.hardcopy(fn,eps=True,color=True,size=(7,5),fontsize=24)
