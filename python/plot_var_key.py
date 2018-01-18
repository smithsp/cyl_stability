#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys
#print sys.argv
#sys.exit()
ref = int(sys.argv[1])
var_key = sys.argv[2]
nq = False
if var_key == 'nq':
  nq = True
  if 'z' in sys.argv:
    var_key = 'Bz0'
  else:
    var_key = 'Bt0'
d1 = get_EV(ref)
print_keys = d1.idata 
d1 = d1.data
dk = dict()
dk_keys = ['t_assemb','num','t_solve',var_key,'epsilo','NN','vcyl', 'rs']#,'k','Bt0','Bz0','p0']
for key in d1.keys():
  if key not in dk_keys:
    dk[key] = d1[key]
matching_files = []
var_list= []
p0_list = []
eps_list = []

g = Gnuplot.Gnuplot(persist=1)
g1 = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g1._clear_queue()
g.xlabel(var_key)
g1.xlabel(var_key)
if nq:
  g.xlabel('nq')
  g1.xlabel('nq')
if var_key=='ngrid':
  g1.xlabel('N')
if var_key == 'Lend0':
  g('set xrange [-1:2]')
  g1('set xrange [-1:2]')
if var_key=='rs':
  g.xlabel('r_s')
  g1.xlabel('r_s')
g.ylabel('{/Symbol w}^2')
if 'ly' in sys.argv:
  g1('set log y')
negty = ''
nfacy = 1 
if 'negy' in sys.argv:
  negty = '-'
  nfacy = -1
g1.ylabel('%s{/Symbol G}'%negty)
scale_fac = 10**0.5
if 'sf' in sys.argv:
  ind = sys.argv.index('sf')
  scale_fac = float(sys.argv[ind+1])
#print 'scale_fac=%g'%(scale_fac)
minevalsr = 0
maxevalsr = 0
minnum = 1
maxnum = 2400
if 'minn' in sys.argv:
  minnum = int(sys.argv[sys.argv.index('minn')+1])
if 'maxn' in sys.argv:
  maxnum = int(sys.argv[sys.argv.index('maxn')+1])
dn = 1
if 'dn' in sys.argv:
  ind = sys.argv.index('dn')
  dn = sys.argv[ind+1]
nums = range(minnum,maxnum,dn)
if 'SN' in sys.argv:
  ind = sys.argv.index('SN')
  nums = []
  for num in sys.argv[ind+1:]:
    try:
      nums.append(int(num))
    except:
      break
maxg = []
for i in nums:
  try:
    d2 = get_EV(i)
    data = d2.data
    if var_key=='alpha' and data['rs']>data['a']:
      data['rs']=dk['rs']
      data['alpha']=0
  except:
    continue
  match = True
  for key in dk.keys(): 
    #if key == 'm':
      #if abs(data[key]) == abs(dk[key]):
        #continue   
    if key=='alpha' and data['rs']>data['a']:
      continue
    if data[key] != dk[key]:
      if 'verb' in sys.argv:
        print 'run %d does not match on key %s'%(i,key)
      if key=='P0' and data['equilibr']==12:
        #print data['Vp0'],d1['Vp0']
        if abs((data['P0']-.5*data['a']**2*data['Vp0']**2)-(d1['P0']-.5*d1['a']**2*d1['Vp0']**2))<1e-5:
          continue
        #else:
          #print data['P0']-.5*data['a']**2*data['Vp0']**2,d1['P0']-.5*d1['a']**2*d1['Vp0']**2
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    matching_files.append(data['num'])
    #nq = data['nz']*data['Bz0']/(20*cmath.pi**2*data['Bt0'])
    if nq:
      var_list.append((data['k']*data['Bz0']/data['Bt0']*data['a']))
    else:
      var_list.append(data[var_key])
    p0_list.append(data['P0'])
    eps_list.append(data['eps'])
    evals = d2.get_numeric_evals()
    maxg.append(max(evals.imag))
    evalsr = evals.real  #Only for initialization to the right size
    for j in range(len(evals)):
      if 'lin' in sys.argv:
        if 'err' in sys.argv:
          evalsr[j] = Numeric.sum(d2.BCerror)
        else:
          evalsr[j] = (evals[j]**2).real
      else:
        evalsr[j] = cmath.asinh(scale_fac*(evals[j]**2).real).real
    if minevalsr > min(evalsr):
      minevalsr = min(evalsr)
    if maxevalsr < max(evalsr):
      maxevalsr = max(evalsr)
    g._add_to_queue([Gnuplot.Data([var_list[-1]]*len(evalsr),evalsr,with='points pt 1 ps .6')])
    g1._add_to_queue([Gnuplot.Data([var_list[-1]]*len(evalsr),nfacy*evals.imag,with='points pt 1')])
for i in range(len(matching_files)):
  print matching_files[i],':',var_list[i],';  ',
print ''
inds = Numeric.argsort(var_list)
g.itemlist = list(Numeric.take(g.itemlist,inds))
g1.itemlist = list(Numeric.take(g1.itemlist,inds))
var_list = list(Numeric.take(var_list,inds))
maxg = Numeric.take(maxg,inds)
#print var_list,maxg
#y = 
if var_key=='ngrid':
  g1._add_to_queue([Gnuplot.Data(var_list,1./Numeric.array(var_list)*var_list[-1]*maxg[-1],with='lines 3')])
#for i in range(len(matching_files)):
  #print abs(var_list[i])

set_ytics = 'set ytics ('
for i in range(2):
  for j in range(-int(cmath.log10(scale_fac).real),13,1):
    scale = (-1)**i*10**j
    if j%2 != 0:
      set_ytics  = set_ytics + '"" %g 1, '%(cmath.asinh(scale*scale_fac).real)
    else:
      set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)
set_ytics = set_ytics + '"0" 0 )'
if 'lin' not in sys.argv:
  g(set_ytics)
if not nq:
  if abs(max(var_list))!=0:
    maxfac = 1+0.1*max(var_list)/abs(max(var_list))
    if abs(min(var_list)/max(var_list))<1/(6.0e1) and min(var_list)!=0:
  #g('set xrange [%g:%g]'%(min(var_list)/1.1,max(var_list)*1.1))
      g('set log x')
      g1('set log x')
    else:
      if min(var_list)!=0:
        minfac = 1-0.1*min(var_list)/abs(min(var_list))
        g('set xrange [%g:%g]'%(min(var_list)*minfac,max(var_list)*maxfac))
        g1('set xrange [%g:%g]'%(min(var_list)*minfac,max(var_list)*maxfac))
      else:
        g('set xrange [:%g]'%(max(var_list)*maxfac))
        g1('set xrange [:%g]'%(max(var_list)*maxfac))

if 'lx' in sys.argv:
  g1('set xrange [:]')
  g1('set log x')
#g1('set yrange [.001:]')
g('set yrange [%g:%g]'%(minevalsr-2,maxevalsr+2))
#if min(var_list)!=0:

tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1[print_keys[i]])
if 'notitle' not in sys.argv:
  g1.title(tt[:-2])
  g.title(tt[:-2])  
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g1('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  g1('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g1('set ytics auto')
  g('set yrange [%s:%s]'%(cmath.asinh(float(sys.argv[ind+1])*scale_fac).real,cmath.asinh(float(sys.argv[ind+2])*scale_fac).real))
  if 'lin' in sys.argv:
    g('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
    g('set ytics auto')
#g.refresh()
if var_key=='tw':
  g1._add_to_queue([Gnuplot.Data(Numeric.sort(Numeric.array(var_list)),2*abs(dk['m'])/Numeric.sort(Numeric.array(var_list)),with='lines')])
if 'imag' in sys.argv:
  g1.refresh()
if nq:
  fn = 'plots/%d_var_%s.eps'%(min(matching_files),'nq')
else:
  fn = 'plots/%d_var_%s.eps'%(min(matching_files),var_key)
#g.refresh()
if 'mod' in sys.argv:
  for i in range(10):
    comm = raw_input('Gnuplot command? ')
    if len(comm)==0:
      break
    g(comm)
    g1(comm)
g.hardcopy(fn,fontsize=18,size=[5,5],eps=True)
fn = '%s_imag.eps'%fn[:-4]
#g1.refresh()
g1.hardcopy(fn,eps=True,fontsize=18,solid=False,color=False)
print 'Plotted to %s'%(fn)
