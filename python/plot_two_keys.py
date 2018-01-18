#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys, sets
ref = int(sys.argv[1])
var_key = sys.argv[2]
var_key2 = sys.argv[3]
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
dk_keys = ['t_assemb','num','t_solve',var_key,var_key2,'epsilo','NN','vcyl']
for key in d1.keys():
  if key not in dk_keys:
    dk[key] = d1[key]
matching_files = []
var_list= []
var2_list = []

g = Gnuplot.Gnuplot(persist=1)
g1 = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g1._clear_queue()
g.xlabel(var_key)
g1.xlabel(var_key)
if nq:
  g.xlabel('nq')
  g1.xlabel('nq')
if var_key == 'Lend0':
  g('set xrange [-1:2]')
  g1('set xrange [-1:2]')
g.ylabel('{/Symbol w}^2')
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
      if key=='P0' and data['equilibr']==12:
        #print data['Vp0'],d1['Vp0']
        if abs((data['P0']-.5*data['a']**2*data['Vp0']**2)-(d1['P0']-.5*d1['a']**2*d1['Vp0']))<1e-5:
          continue
        else:
          print i,data['P0']-.5*data['a']**2*data['Vp0']**2,d1['P0']-.5*d1['a']**2*d1['Vp0']
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    matching_files.append(data['num'])
    if nq:
      var_list.append((data['k']*data['Bz0']/data['Bt0']*data['a']))
    else:
      var_list.append(data[var_key])
    var2_list.append(data[var_key2])
uniq_var2 = list(sets.Set(var2_list))
uniq_var2.sort()
k = 1
#print '\n',d1['tw'],data['tw'],'\n'
if 'var2r' in sys.argv:
  ind = sys.argv.index('var2r')
  var2_list2 = []
  for i in range(len(sys.argv[ind+1:])):
    try:
      var2_list2.append(float(sys.argv[ind+1+i]))
    except:
      break
  var2_list2.sort()
  if 'rev' in sys.argv:
    var2_list2.reverse()
  uniq_var2 = var2_list2
if var_key2 == 'tw':
  try:
    ind = uniq_var2.index(-1)
    uniq_var2.append(uniq_var2[ind])
    uniq_var2.pop(ind)
  except:
    pass
print 'Length of uniq_var2=%d'%len(uniq_var2)
for v2 in uniq_var2:
  ttl = True #The title is needed
  for i in range(len(matching_files)):
    #print matching_files[i],
    if var2_list[i]==v2:
      d2 = get_EV(matching_files[i])
      if 'tw-1' in sys.argv and d2.data['tw']!=0:
        nfacy = d2.data['tw']
      evals = d2.evals
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
      if ttl:
        g._add_to_queue([Gnuplot.Data([var_list[i]]*len(evalsr),evalsr,with='points %i '%k,title='%s=%g'%(var_key2,v2))])
        g1._add_to_queue([Gnuplot.Data([var_list[i]]*len(evalsr),nfacy*evals.imag,with='points  %i'%k,title='%s=%g'%(var_key2,v2))])
        ttl = False
      else:
        g._add_to_queue([Gnuplot.Data([var_list[i]]*len(evalsr),evalsr,with='points %i '%k)])
        g1._add_to_queue([Gnuplot.Data([var_list[i]]*len(evalsr),nfacy*evals.imag,with='points %i'%k)])
  k = k+1
#for i in range(len(matching_files)):
  #print matching_files[i],':',var_list[i],';  ',
print matching_files

set_ytics = 'set ytics ('
for i in range(2):
  for j in range(-5,13,1):
    scale = (-1)**i*10**j
    if j%2 != 0:
      set_ytics  = set_ytics + '"" %g 1, '%(cmath.asinh(scale*scale_fac).real)
    else:
      set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)
set_ytics = set_ytics + '"0" 0 )'
if 'lin' not in sys.argv:
  g(set_ytics)
if not nq:
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
#g1('set yrange [.001:]')
g('set yrange [%g:%g]'%(minevalsr-2,maxevalsr+2))
#if min(var_list)!=0:
g('set key outside')
g1('set key Left reverse outside')
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1[print_keys[i]])
g.title(tt[:-2])
g1.title(tt[:-2])
if 'notitle' in sys.argv:
  g('unset title')
  g1('unset title')
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g1('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  if 'lin' in sys.argv:
    g('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
    g('set ytics auto')
  g1('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g1('set ytics auto')
g1.refresh()
if var_key=='tw':
  g1.xlabel('{/Symbol t}_w')
  if 'tw-1' not in sys.argv:
    g1._add_to_queue([Gnuplot.Data(Numeric.sort(Numeric.array(var_list)),2*abs(dk['m'])/Numeric.sort(Numeric.array(var_list)),with='lines')])
if 'imag' in sys.argv:
  g1.refresh()
if nq:
  fn = 'plots/%d_var_%s.eps'%(min(matching_files),'nq')
else:
  fn = 'plots/%d_var_%s_%s.eps'%(min(matching_files),var_key,var_key2)
print 'Plotted to %s'%(fn)
#g.refresh()
g.hardcopy(fn,fontsize=14)
g1.hardcopy('%s_imag.eps'%fn[:-4],eps=True,fontsize=12)#,fontsize=24)#,eps=True)
