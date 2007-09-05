from get_EV import *
import Gnuplot, Numeric, cmath
equilib = '3'
dk = {'a': 1.0, 'epsVp': 0.0, 'ngrid': 20, 'rho0': 1.0, 'NN': 126, 'Bz0': 1.0, 'k': 0.20000000000000001, 'epsVz': 0.0, 'm': 2, 'equilibr': int(equilib), 'nz': 1, 'Vz0': 0.0, 'Vp0': 0.0}
matching_files = []
nq_list= []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue
print dk['equilibr']
scale_fac = 10**5
for i in range(1,1000):
  try:
    data = get_data(equilib,'%d'%i)
  except:
    break
  match = True
  for key in dk.keys():
    if data[key] != dk[key]:
      match = False
      break
  if match:
    matching_files.append(i)
    nq = (data['k']*data['Bz0']/data['Bt0']/data['a'])
    #nq = data['nz']*data['Bz0']/(20*cmath.pi**2*data['Bt0'])
    nq_list.append(nq)
    evals = get_evals(equilib,'%d'%i,data)
    evalsr = evals.real  #Only for initialization to the right size
    for j in range(len(evals)):
      evalsr[j] = cmath.asinh(scale_fac*(evals[j]**2).real).real
    g._add_to_queue([Gnuplot.Data([nq]*len(evalsr),evalsr)])
#if len(matching_files) != 0:
g._add_to_queue([Gnuplot.Data([min(nq_list),max(nq_list)],[cmath.asinh(0).real]*2,with='lines')])
g.title('Bz0 = %g, rho0 = %g, nz = %d, ar = %g'%(dk['Bz0'],dk['rho0'],dk['nz'],dk['a']))
set_ytics = 'set ytics ('
for i in range(2):
  for j in range(-4,6,2):
    scale = (-1)**i*10**j
    set_ytics  = set_ytics + '"%g" %g, '%(scale,cmath.asinh(scale*scale_fac).real)

set_ytics = set_ytics + '"0" 0 )'
#print set_ytics
#sys.exit()
g(set_ytics)
g.refresh()
#else:
  #print 'There were no matching files.  Recheck dk.'
g.hardcopy('equilib%s_varnq_m=%d.eps'%(equilib,dk['m']))
