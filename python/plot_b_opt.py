#! /usr/pppl/python/2.3.5x/bin/python
import Gnuplot,sys,Numeric,sets,math
from get_EV import *
ref = int(sys.argv[1])
var_key2 = 'ngrid'#'ngrid'#sys.argv[2]
var_key = 'nu'#sys.argv[3]
d1 = get_EV(ref)
maxr = d1.data['k']*d1.data['Vz0']
minr = .04
maxr = 1e10
minr = -1e10
print 'Warning: The filter on the real part of omega is %g<Re(omega)<%g'%(minr,maxr)
hexp = 0
if 'hexp' in sys.argv:
  ind = sys.argv.index('hexp')
  hexp = float(sys.argv[ind+1])
print_keys = d1.idata 
d1 = d1.data
dk = dict()
nq=False
dk_keys = ['t_assemb','num','t_solve',var_key,var_key2,'epsilo','NN','vcyl','rs','nu']
if 'rs' in dk_keys or 'alpha' in dk_keys:
  print 'Shots of incongruent grids are included.'
for key in d1.keys():
  if key not in dk_keys:
    dk[key] = d1[key]
matching_files = []
var_list= []
var2_list = []
max_imag = []
g = Gnuplot.Gnuplot(persist=0)
g._clear_queue()
g.xlabel(var_key)
if var_key=='nu':
  g.xlabel('nu*h**%g'%hexp)
g.ylabel('{/Symbol G}')
minnum = 1
maxnum = 1750
if 'minn' in sys.argv:
  minnum = int(sys.argv[sys.argv.index('minn')+1])
if 'maxn' in sys.argv:
  maxnum = int(sys.argv[sys.argv.index('maxn')+1])
for i in range(minnum,maxnum,1):
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
    if key=='alpha' and data['rs']>data['a']:
      continue
    if data[key] != dk[key]:
      if key=='P0' and data['equilibr']==12:
        if abs((data['P0']-.5*data['a']**2*data['Vp0']**2)-(d1['P0']-.5*d1['a']**2*d1['Vp0']**2))<1e-5:
          continue
      match = False
      break
  if match:    
    matching_files.append(data['num'])
    try:
      h_min = min(d2.grid[1:]-d2.grid[0:-1])
    except:
      h_min = 1./(d2.data['ngrid']-1)
    if nq:
      var_list.append((data['k']*data['Bz0']/data['Bt0']*data['a']))
    elif var_key=='nu':
      var_list.append(abs(data[var_key])*h_min**hexp)
    else:
      var_list.append(data[var_key])
    if var_key2=='ngrid':# and var_key=='nu':
      var2_list.append(round(1./h_min)+1)
    else:
      var2_list.append(data[var_key2])
    evals = d2.get_unstable_evals()
    inds = list(Numeric.nonzero((evals.real<maxr)*(evals.real>minr)))
    max_imag.append(evals[inds[0]].imag)
print matching_files
print var2_list
uniq_var2 = list(sets.Set(var2_list))
uniq_var2.sort()
k = 1
if 'var2r' in sys.argv:
  ind = sys.argv.index('var2r')
  var2_list2 = list(sys.argv[ind+1:])
  for i in range(len(var2_list2)):
    var2_list2[i] = float(var2_list2[i])
  var2_list2.sort()
  uniq_var2 = var2_list2
b_opt = []
opt_b_local = []
v_opt = []

for v2 in uniq_var2:
  max_eval = []
  var_list1 = []
  for i in range(len(matching_files)):
    if var2_list[i]==v2:
      max_eval.append(max_imag[i])
      var_list1.append(var_list[i])
  ind = Numeric.argsort(max_eval)
  inds = list(Numeric.argsort(var_list1))
  inds.reverse()
  for i in range(len(inds)-1):
    if var_key=='nu':
      if max_eval[inds[i]]*2>max_eval[inds[i+1]]:
        continue
      else:
        h_min = 1./v2
        b_opt.append((var_list1[inds[i]]*1./h_min**hexp))
        v_opt.append((h_min))
        break
    else:
      if max_eval[inds[i]]>max_eval[inds[i+1]]:
        continue
      else:
        b_opt.append(var_list1[inds[i]])
        v_opt.append((v2))
        break
  x = var_list1
  x = Numeric.take(var_list1,inds)
  y = Numeric.take(max_eval,inds)
  ux = list(sets.Set(x))
  x_min = []
  y_min = []
  for ux1 in ux:
    x_min.append(ux1)
    inds = Numeric.nonzero(x==ux1)
    y_min.append(min(Numeric.take(y,inds)))
  x = x_min
  y = Numeric.array(y_min)
  if 'dg' in sys.argv: #For difference between growth rate and a reference
    y = abs(y-min(y))/min(y)*100
    g.ylabel('|{/Symbol G}-{/Symbol G}_o|/{/Symbol G}_o (%)')
  g._add_to_queue([Gnuplot.Data(x,y,with='points %i'%k,title='%s=%g'%(var_key2,v2))])
  k = k+1
g('set key Left reverse outside')
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1[print_keys[i]])
g.title(tt[:-2])
if 'ly' in sys.argv:
  g('set log y')
if 'lx' in sys.argv:
  g('set log x')
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  g('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
#g.refresh()
if nq:
  fn = 'plots/%d_var_%s.eps'%(min(matching_files),'nq')
elif var_key=='nu':
  fn = 'plots/%d_var_%sh%g_%s_most_unstable.eps'%(min(matching_files),var_key,hexp,var_key2)
else:
  fn = 'plots/%d_var_%s_%s_most_unstable.eps'%(min(matching_files),var_key,var_key2)
print 'Plotted to %s'%(fn)
#g.refresh()
g.hardcopy(fn,fontsize=12)
g('reset')
if 'lx1' in sys.argv:
  g('set log x')
if 'ly1' in sys.argv:
  g('set log y')
if 'xr1' in sys.argv:
  ind = sys.argv.index('xr1')
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'yr1' in sys.argv:
  ind = sys.argv.index('yr1')
  g('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
#sys.exit()
g.xlabel('h')
g.ylabel('{/Symbol n}_{crit}')
g._clear_queue()
g._add_to_queue([Gnuplot.Data(v_opt,b_opt,with='points')])
g._add_to_queue([Gnuplot.Data(v_opt,2.*Numeric.array(v_opt)**(-2.7)*math.exp(-23.95),with='lines',title='h^{-2.7}e^{-24}')])
g('set key reverse center')
fn = 'plots/%d_var_%s_opt_%s.eps'%(min(matching_files),var_key2,var_key)
print 'Plotted to %s'%(fn)
g.hardcopy(fn,fontsize=12)
#print matching_files
