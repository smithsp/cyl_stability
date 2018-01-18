#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys
ref = sys.argv[1]
var_key = sys.argv[2]
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
if 'lx' in sys.argv:
  g('set log x')
if 'ly' in sys.argv:
  g('set log y')
negty = ''
nfacy = 1 
if 'negy' in sys.argv:
  negty = '-'
  nfacy = -1
negtx = ''
nfacx = 1 
if 'negx' in sys.argv:
  negtx = '-'
  nfacx = -1
g._clear_queue()
g.xlabel('%sRe({/Symbol w})'%(negtx))
g.ylabel('%sIm({/Symbol w})'%(negty))
for i in range(1,1000):
  try:
    d2 = get_EV(i)
    data = d2.data
  except IOError:
    #print 'IOError for i=%d'%i
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
    p0_list.append(data['P0'])
    eps_list.append(data['eps'])
    evals = d2.get_numeric_evals() #get_unstable_evals()
    evalsr = evals.real  #Only for initialization to the right size
    #w0 = Numeric.sum(evalsr)/len(evalsr)
    try:
      g._add_to_queue([Gnuplot.Data(nfacx*evals.real,nfacy*evals.imag,title='%s=%2.2e'%(var_key,(data[var_key])))])
      if i == int(ref):
        g.itemlist[-1].set_option(with='points ps 3')
    except (ValueError):
      continue
    var_list.append(data[var_key])
inds = Numeric.argsort(var_list)
print matching_files
print var_list
g.itemlist = list(Numeric.take(g.itemlist,inds))
tt = ''
for i in range(len(print_keys)):
  if print_keys[i] != var_key:
    tt = tt+'%s=%g, '%(print_keys[i],d1.data[print_keys[i]])
g.title(tt[:-2])
g('set key outside  under')
#g('set xrange [-.5:.5]')
if var_key!='Vz0': 
  pass 
  #g('set xrange [%g:%g]'%(-.5,d1.get_w0()+5.11e-2))
  #g('set yrange [1e-6:1e-1]')
  #g('set size square')
print d1.get_w0()
#if len(matching_files)>1:
fn = 'plots/%d_var_%s_unst.eps'%(min(matching_files),var_key)
print 'Plotted to %s'%(fn)
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  g('set xrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
if 'yr' in sys.argv:
  ind = sys.argv.index('yr')
  g('set yrange [%s:%s]'%(sys.argv[ind+1],sys.argv[ind+2]))
  g('set ytics auto')
g.refresh()
g.hardcopy(fn)
  
