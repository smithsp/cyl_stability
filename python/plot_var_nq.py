#! /usr/pppl/python/2.3.5x/bin/python
from get_EV import *
import Gnuplot, Numeric, cmath, sys
equilib = '3'
dk = {'a': 1.0,  'ngrid': 20, 'rho0': 1.0, 'tw':10,'b': 2.0,'Bt0': 1.0,  'k': 0.20000000000000001,  'm': 2, 'equilibr': int(equilib),  'fe_type':'spline'}
ref = sys.argv[1]
dk = get_EV(ref).data
equilib = dk['equilibr']
matching_files = []
nq_list= []
p0_list = []
eps_list = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue
g.xlabel('nq')
g.ylabel('{/Symbol w}^2')
print dk['equilibr']
scale_fac = 10**7
for i in range(1,1000):
  try:
    d2 = get_EV(i)
    data = d2.data
  except:
    continue
  match = True
  for key in dk.keys(): 
    if key == 'm':
      if abs(data[key]) == abs(dk[key]):
        continue   
    if data[key] != dk[key]:
      match = False
      #print data[key]
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    try:
      nq = (data['k']*data['Bz0']/data['Bt0']*data['a'])*data['m']/abs(data['m'])
      print data['Bt0']
    except ZeroDivisionError:
      continue
    matching_files.append(data['num'])
    #nq = data['nz']*data['Bz0']/(20*cmath.pi**2*data['Bt0'])
    nq_list.append(nq)
    p0_list.append(data['P0'])
    eps_list.append(data['eps'])
    evals = d2.evals
    evalsr = evals.real  #Only for initialization to the right size
    for j in range(len(evals)):
      evalsr[j] = cmath.asinh(scale_fac*(evals[j]**2).real).real
    g._add_to_queue([Gnuplot.Data([nq]*len(evalsr),evalsr,with='points 1')])
#if len(matching_files) != 0:
print matching_files, nq_list
#print nq_list, p0_list, eps_list, matching_files
#g('set size 2.0,1.0')
#g('set xrange [-2.1:-1.9]')
g._add_to_queue([Gnuplot.Data([min(nq_list),max(nq_list)],[cmath.asinh(0).real]*2,with='lines')])
g.title('br = %g, tw = %g'%(dk['b'],dk['tw']))
set_ytics = 'set ytics ('
for i in range(2):
  for j in range(-6,7,2):
    scale = (-1)**i*10**j
    set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)

set_ytics = set_ytics + '"0" 0 )'
#print set_ytics
#sys.exit()
g('set xrange [-2.1:-1.9]')
#g('set yrange [%g:%g]'%(cmath.asinh(-1.e-3*scale_fac).real,cmath.asinh(1.e5*scale_fac).real))
g(set_ytics)
g.refresh()
#else:
  #print 'There were no matching files.  Recheck dk.'
fn = 'plots/%d_var_nq.eps'%(min(matching_files))
print fn
g.hardcopy(fn)
#import os
#os.system('gv %s &'%fn)
