from get_EV import *
import sys, math
gamma_mt_1 = [1.841183781,3.054236928 ]
gamma_alpha = [1.841183781,5.331442774,8.536316366,11.7060049,14.86358863,18.01552786,21.16436986,24.31132686,27.45705057,30.60192297,33.7461829,3.68899874E+01,40.03344405,43.17662897,46.31959756,49.46239114,52.60504111,55.74757179,58.8900023]
gamma_ind = int(sys.argv[1])
dk = {'a': 1.0, 'rho0': 20.0, 'Bt0': 0.0, 'Bz0': 0.02,  'k': 0.5, 'm': 1, 'equilibr': 1,  's2': 0.083333333333333301}
if dk['m'] != 1:
  print 'm must equal 1'
  sys.exit()
matching_files = []
N = []
N_lin = []
error = []
err_lin = []
A = dk['a']**2
s2 = dk['s2']
omega1 = dk['Bz0']**2/(2*dk['rho0']*s2/(5./3.)*dk['Bz0']**2*A)*(A*dk['k']**2+gamma_alpha[gamma_ind]**2)*(1+s2)*(1-(1-4*s2*A*dk['k']**2/((1+s2)**2*(A*dk['k']**2+gamma_alpha[gamma_ind]**2)))**(1./2.))
slow_inf = s2*dk['Bz0']**2*dk['k']**2/(dk['rho0']*(1+s2)*s2/(5./3.)*dk['Bz0']**2)
for i in range(1,1000):
  try:
    d = get_EV(i,output_folder='output_cyl')
    data = d.data
  except:
    continue
  match = True
  for key in dk.keys():    
    if key == 's2':
      continue
    if data[key] != dk[key]:        
      match = False
      break
  if match:
    matching_files.append(i)
    evals = d.evals
    evalsr = Numeric.absolute(evals.real)
    if data['fe_type'] == 'spline':
      N.append(data['ngrid']+1)
      error.append(min(Numeric.absolute(evalsr-omega1**(1))/omega1**1))
    else:
      N_lin.append(data['ngrid']-1)
      err_lin.append(min(Numeric.absolute(evalsr-omega1**(1))/omega1**1))
N = Numeric.array(N)
inds = Numeric.argsort(N)
N = Numeric.take(N,inds)
error = Numeric.take(error,inds)
N_exp = 1
d1 = [Gnuplot.Data((N)**N_exp,error,with='linespoints ',title='Cubic B-spline')]
avg_log = []
N_fac = float(sys.argv[2])
N_ind = Numeric.nonzero(N>gamma_ind+5)
for i in N_ind:#range(max(gamma_ind+5-min(N),0),len(error)):
  avg_log.append(math.log(error[i])/math.log(N_fac*N[i]))
#print avg_log
avg_log =  Numeric.sum(avg_log)/len(avg_log)
N = Numeric.take(N,N_ind)
d1.append(Gnuplot.Data(N,(N_fac*N)**(avg_log),with='lines',title='(%g*N)^{%g}'%(N_fac,avg_log)))

N = N_lin
inds = Numeric.argsort(N)
N = Numeric.take(N,inds)
error = Numeric.take(err_lin,inds)
N_exp = 1
d1.append(Gnuplot.Data((N)**N_exp,error,with='linespoints ',title='Linear'))
avg_log = []
N_fac = float(sys.argv[3])
N_ind = Numeric.nonzero(N>gamma_ind+5)
for i in N_ind:#range(max(gamma_ind+5-min(N),0),len(error)):
  avg_log.append(math.log(error[i])/math.log(N_fac*N[i]))
#print avg_log
avg_log =  Numeric.sum(avg_log)/len(avg_log)
N = Numeric.take(N,N_ind)
d1.append(Gnuplot.Data(N,(N_fac*N)**(avg_log),with='lines',title='(%g*N)^{%g}'%(N_fac,avg_log)))


g = Gnuplot.Gnuplot(persist=1)
#g('set xrange [*:%g]'%((1./(gamma_ind+3))**N_exp))
g('set logscale y')
g('set key bottom left')
g('set key spacing 1.2')
g.xlabel('N')
g.ylabel('Error')
g.title('Error in slow wave frequency # %d'%gamma_ind)
g._clear_queue
g._add_to_queue(d1)
g.refresh()
g.hardcopy('plots/slow_wave_converge_%d.eps'%gamma_ind)

