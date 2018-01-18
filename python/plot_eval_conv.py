#! /usr/pppl/python/2.3.5x/bin/python
import get_EV, Gnuplot, sys, Numeric, sets
ref = sys.argv[1]
unst_ind = int(sys.argv[2])
nexp = 1
if 'nexp' in sys.argv:  
  ind = sys.argv.index('nexp')
  nexp = float(sys.argv[ind+1])
var_key = 'ngrid'
d1 = get_EV.get_EV(ref)
dk = dict()
print_keys = d1.idata 
d1 = d1.data
dk_keys = ['t_assemb','num','t_solve',var_key,'epsilo','NN','vcyl','alpha','rs','nu']
for key in d1.keys():
  if key not in dk_keys:
    dk[key] = d1[key]
matching_files = []
var_list= []
ev_list = []
ev1_list = []
g = Gnuplot.Gnuplot(persist=0)
g1 = Gnuplot.Gnuplot(persist=0)
g._clear_queue()
g1._clear_queue()
g.xlabel('h_{min}^{%g}'%nexp)
g1.xlabel('h_{min}^{%g}'%nexp)
g.ylabel('{/Symbol G}')#%unst_ind)
g1.ylabel('Re({/Symbol w}) of %dst unstable RWM'%unst_ind)
minnum = 2290
maxnum = 2400
if 'minn' in sys.argv:
  minnum = int(sys.argv[sys.argv.index('minn')+1])
if 'maxn' in sys.argv:
  maxnum = int(sys.argv[sys.argv.index('maxn')+1])
for i in range(minnum,maxnum,1):
  try:
    d2 = get_EV.get_EV(i)
    data = d2.data
  except IOError:
    continue
  match = True
  for key in dk.keys(): 
    if data[key] != dk[key]:
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:
    matching_files.append(data['num'])
    evals = d2.get_unstable_evals()
    inds = Numeric.nonzero((evals.real<1e10)*(evals.real>-1e10))
    neff = data['ngrid']-1.
    if data['alpha']!=0:
      try:
        shift_grid = d2.grid[1:]-d2.grid[:-1]
        neff = data['a']/min(shift_grid)
      except:
        print 'for %d, grid was not available'%(i)
        neff = data['ngrid']-1.
    var_list.append(1./neff**nexp)
    ev_list.append(evals[inds[unst_ind]].imag)
    ev1_list.append(evals[inds[unst_ind]].real)
#print matching_files
#print var_list#,ev_list,ev1_list
ux = list(sets.Set(var_list))
x2 = []
y2 = []
z2 = []
uux = []
min_files = []
for uh1 in ux:
  if len(Numeric.nonzero(Numeric.absolute(Numeric.array(uux)-uh1)<1e-8))==0:
    uux.append(uh1)
for x1 in uux:
  x2.append(x1)
  inds = Numeric.nonzero(Numeric.absolute(Numeric.array(var_list)-x1)<1e-8)
  z2.append(min(Numeric.take(ev_list,inds)))
  inds = Numeric.nonzero(Numeric.array(ev_list)==z2[-1])
  min_files.append(Numeric.take(matching_files,inds)[0])
  y2a = Numeric.take(ev1_list,inds)
  y2.append(y2a[0])
var_list = x2
ev_list = z2
ev1_list = y2
#print 1./Numeric.array(var_list)
inds = Numeric.argsort(var_list)
var_list = list(Numeric.take(var_list,inds))
ev_list = list(Numeric.take(ev_list,inds))
ev1_list = list(Numeric.take(ev1_list,inds))
min_files = list(Numeric.take(min_files,inds))
for i in range(len(x2)):
  print var_list[i],':',min_files[i],' ',
print ''
#print var_list[0:2],ev_list[0:2]
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1[print_keys[i]])
#g.title(tt[:-2])
#g1.title(tt[:-2])
g('set xrange [0:%g]'%(max(var_list)*1.1))
g1('set xrange [0:%g]'%(max(var_list)*1.1))
if 'xr' in sys.argv:
  ind = int(sys.argv.index('xr'))
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g1('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'ly' in sys.argv:
  g('set log y')
  g1('set log y')
elif 'xr' not in sys.argv:
  g('set yrange [0:%g]'%(max(ev_list)*1.1))
  g1('set yrange [%g:%g]'%(min(ev1_list)*(1-.1*min(ev1_list)/abs(min(ev1_list))),max(ev1_list)*1.1))
g.itemlist.append(Gnuplot.Data(var_list,ev_list,with='points'))
g1.itemlist.append(Gnuplot.Data(var_list,ev1_list,with='points'))
m=(ev_list[1]-ev_list[0])/(var_list[1]-var_list[0])
m1=(ev1_list[1]-ev1_list[0])/(var_list[1]-var_list[0])
var_list.insert(0,0)
ev_list.insert(0,0)
ev1_list.insert(0,0)
#print var_list
imag_fit = m*Numeric.array(var_list)+ev_list[1]-m*var_list[1]
print imag_fit[0]
real_fit = m1*Numeric.array(var_list)+ev1_list[1]-m1*var_list[1]
g.itemlist.append(Gnuplot.Data(var_list,imag_fit,with='lines'))
g1.itemlist.append(Gnuplot.Data(var_list,real_fit,with='lines',title='linear extrapolation'))
fn = 'plots/%d_im_conv_N%g_%d'%(min(matching_files),nexp,unst_ind)
fn1 = 'plots/%d_re_conv_N%g_%d'%(min(matching_files),nexp,unst_ind)
if 'ly' in sys.argv:
  fn = fn+'_log.eps'
  fn1 = fn1+'_log.eps'
else:
  fn = fn+'.eps'
  fn1 = fn1+'.eps'
print 'output to %s and %s'%(fn,fn1)
g('set key bottom right Left reverse')
g1('set key center right')
g.hardcopy(fn,eps=True)
g1.hardcopy(fn1)
