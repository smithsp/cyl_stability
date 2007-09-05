from get_EV import *
import sys
gamma_mt_1 = [1.841183781,3.054236928 ]
gamma_alpha = [1.841183781,5.331442774,8.536316366,11.7060049,14.86358863,18.01552786,21.16436986,24.31132686,27.45705057,30.60192297,33.7461829,3.68899874E+01,40.03344405,43.17662897,46.31959756,49.46239114,52.60504111,55.74757179,58.8900023]
gamma_ind = int(sys.argv[1])
dk = {'a': 1.0, 'rho0': 20.0, 'Bt0': 0.0, 'Bz0': 0.02,  'k': 0.5, 'epsVz': 0.0, 'm': 1, 'equilibr': 1,  'Vz0': 0.0,'s2': 0.083333333333333301}
if dk['m'] != 1:
  print 'm must equal 1'
  sys.exit()
matching_files = []
N = []
error = []
A = dk['a']**2
s2 = dk['s2']
omega1 = dk['Bz0']**2/(2*dk['rho0']*s2/(5./3.)*dk['Bz0']**2*A)*(A*dk['k']**2+gamma_alpha[gamma_ind]**2)*(1+s2)*(1-(1-4*s2*A*dk['k']**2/((1+s2)**2*(A*dk['k']**2+gamma_alpha[gamma_ind]**2)))**(1./2.))
slow_inf = s2*dk['Bz0']**2*dk['k']**2/(dk['rho0']*(1+s2)*s2/(5./3.)*dk['Bz0']**2)
for i in range(1,100):
  try:
    data = get_data('1','%d'%i)
  except:
    break
  match = True
  for key in dk.keys():    
    if key == 's2':
      continue
    if data[key] != dk[key]:        
      match = False
      break
  if match:
    matching_files.append(i)
    N.append(data['ngrid'])
    evals = get_evals('1','%d'%i,data)
    evalsr = Numeric.absolute(evals.real)
    error.append(min(Numeric.absolute(evalsr-omega1**(0.5))/omega1**.5))
N = Numeric.array(N)
inds = Numeric.argsort(N)
N = Numeric.take(N,inds)
error = Numeric.take(error,inds)
g = Gnuplot.Gnuplot(persist=1)
N_exp = 1
#g('set xrange [*:%g]'%((1./(gamma_ind+3))**N_exp))
g('set logscale y')
g.xlabel('N')
g.ylabel('Error')
g.title('Error in slow wave frequency # %d'%gamma_ind)
d1 = [Gnuplot.Data((N)**N_exp,error,with='linespoints ')]
import math
avg_log = []
N_fac = 3
for i in range(gamma_ind+3,len(error)):
  avg_log.append(math.log(error[i])/math.log(N_fac*N[i]))
#print avg_log
avg_log =  Numeric.sum(avg_log)/len(avg_log)
d1.append(Gnuplot.Data(N,(N_fac*N)**(avg_log),with='lines',title='(%g*N)^{%g}'%(N_fac,avg_log)))
g._clear_queue
g._add_to_queue(d1)
g.refresh()
g.hardcopy('slow_wave_converge_%d.eps'%gamma_ind)

