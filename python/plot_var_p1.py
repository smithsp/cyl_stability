from get_EV import *
import sys, Numeric, Gnuplot, cmath
equilib = '3'
dk = {'a': 1.0,'lambd': 3.1760000000000002,  'ngrid': 20, 'rho0': 1.0,   'k': 1.20000000000000001,  'm': 1, 'equilibr': 4, 'nz': 1, 'fe_type':'spline','Vz0':0.0,'P0':0.050000000000000003}
matching_files = []
p1_list= []
p0_list = []
eps_list = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue
g.xlabel('p1')
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
      print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    matching_files.append(data['num'])
    p1_list.append(data['P1'])
    p0_list.append(data['P0'])
    eps_list.append(data['eps'])
    evals = d2.evals
    evalsr = evals.real  #Only for initialization to the right size
    for j in range(len(evals)):
      evalsr[j] = cmath.asinh(scale_fac*(evals[j]**2).real).real
    g._add_to_queue([Gnuplot.Data([data['P1']]*len(evalsr),evalsr)])
#if len(matching_files) != 0:

print p1_list, p0_list, eps_list, matching_files
#g('set size 2.0,1.0')
#g('set xrange [-2.1:-1.9]')
g._add_to_queue([Gnuplot.Data([min(p1_list),max(p1_list)],[cmath.asinh(0).real]*2,with='lines')])
g.title('Vz0=%g'%(dk['Vz0']))
set_ytics = 'set ytics ('
for i in range(2):
  for j in range(-6,7,2):
    scale = (-1)**i*10**j
    set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)

set_ytics = set_ytics + '"0" 0 )'
#print set_ytics
#sys.exit()
g(set_ytics)
g.refresh()
#else:
  #print 'There were no matching files.  Recheck dk.'
#g.hardcopy('equilib%s_varnq_m=%d.eps'%(equilib,dk['m']))
