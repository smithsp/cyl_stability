#! /usr/pppl/python/2.3.5x/bin/python
Bz = 2e-2
rho = 20
R = 1
k = .5
from pygsl import sf
from get_EV import *
import Numeric
import Gnuplot
import sys
m=1
gamma = ([1.841183781, 5.331442774, 8.536316366, 11.7060049, 14.86358863, 18.01552786, 21.16436986, 24.31132686, 27.45705057, 30.60192297, 33.7461829, 3.68899874E+01, 40.03344405, 43.17662897, 46.31959756, 49.46239114, 52.60504111, 55.74757179, 58.8900023])
N = 10
d = get_EV(int(sys.argv[1]))
data = d.data
nphi = data['ngrid']
xgrid = d.grid
x = Numeric.arange(0,data['a'],data['a']/100.)
x[0] = x[1]/10000.0
y = Numeric.zeros(len(x),typecode=Numeric.Float)
ygrid = Numeric.zeros(len(xgrid),typecode=Numeric.Float)
dk = data
A = dk['a']**2
s2 = dk['s2']
evals = d.evals.real**2
if data['fe_type'] == 'linconst':
  xgrid = xgrid+xgrid[1]/2.
for gamma_ind in range(6):#min(nphi,len(gamma)-1)):#nphi,2):
  for i in range(len(x)):
    y[i] = sf.bessel_J1(x[i]*gamma[gamma_ind])[0]
  for i in range(len(xgrid)):
    ygrid[i] = sf.bessel_J1(xgrid[i]*gamma[gamma_ind])[0]
  omega1 = dk['Bz0']**2/(2*dk['rho0']*s2/(5./3.)*dk['Bz0']**2*A)*(A*dk['k']**2+gamma[gamma_ind]**2)*(1+s2)*(1-(1-4*s2*A*dk['k']**2/((1+s2)**2*(A*dk['k']**2+gamma[gamma_ind]**2)))**(1./2.))
  min_ind = Numeric.argsort(Numeric.absolute(evals-omega1))[0]
  print 'min_ind=',min_ind
  error = evals[min_ind]-omega1
  #EV =get_evec_3(data,3,min_ind)
  EV = d.get_evec1('3',min_ind,x).real
  fac =((EV[-2]+EV[-2])/2./ygrid[-2])
  fac = EV[-1]/y[-1]
  if data['fe_type'] == 'linconst':
    fac = EV[0]/(ygrid[0]+ygrid[1])*4
  g = Gnuplot.Gnuplot(persist=1)
  #print 'x=',x,'y=',y,'EV=',EV,'fac=',fac
  g.plot(Gnuplot.Data(x,y,with='lines',title='Exact'),Gnuplot.Data(x,EV/fac,title='Computed'),title='{/Symbol w}^2=%g; Analytic=%g'%(evals[min_ind],omega1))
  #if gamma_ind < 10:
    #g.hardcopy('plots/%s.evecz_slow_%d.eps'%(sys.argv[1],min_ind),eps=True)

sys.exit()

f = file('lin_slow_EVs.txt')
ngrid = int(f.readline())
#print ngrid
nphi = int(f.readline())
xgrid = Numeric.zeros(ngrid,typecode=Numeric.Float)
for i in range(ngrid):
  xgrid[i] = float(f.readline())
EVs = Numeric.zeros([nphi,ngrid],typecode=Numeric.Float)
#print Numeric.shape(EVs)
for j in range(nphi):
  b = f.readline()
  for i in range(ngrid):
    EVs[j,i] = float(f.readline())
x = Numeric.arange(0,R,R/100.)
y = Numeric.zeros(len(x),typecode=Numeric.Float)
ygrid = Numeric.zeros(len(xgrid),typecode=Numeric.Float)
for j in range(min(nphi,len(gamma)-1,1)): #nphi
  for i in range(len(x)):
    y[i] = sf.bessel_J1(x[i]*gamma[j])[0]
  for i in range(len(xgrid)):
    ygrid[i] = sf.bessel_J1(xgrid[i]*gamma[j])[0]
  fac = EVs[j,-1]/ygrid[-1]
  g = Gnuplot.Gnuplot(persist=1)
  g.plot(Gnuplot.Data(x,y,with='lines',title='Exact'),Gnuplot.Data(xgrid,EVs[j,:]/fac,with='linespoints',title='Computed'),title='Linear Elements')
  #g.hardcopy('linear/linear_const_ev_%d.eps'%j)
f.close()

f = file('hermite_slow_EVs.txt')
ngrid = int(f.readline())
#print ngrid
nphi = int(f.readline())
xgrid = Numeric.zeros(ngrid,typecode=Numeric.Float)
for i in range(ngrid):
  xgrid[i] = float(f.readline())
EVs = Numeric.zeros([nphi,ngrid],typecode=Numeric.Float)
#print Numeric.shape(EVs)
for j in range(nphi):
  b = f.readline()
  for i in range(ngrid):
    EVs[j,i] = float(f.readline())
x = Numeric.arange(0,R,R/100.)
y = Numeric.zeros(len(x),typecode=Numeric.Float)
ygrid = Numeric.zeros(len(xgrid),typecode=Numeric.Float)
for j in range(min(nphi,len(gamma)-1,3)):#nphi,2):
  for i in range(len(x)):
    y[i] = sf.bessel_J1(x[i]*gamma[j])[0]
  for i in range(len(xgrid)):
    ygrid[i] = sf.bessel_J1(xgrid[i]*gamma[j])[0]
  fac =(EVs[j,-1]/ygrid[-1])
  g = Gnuplot.Gnuplot(persist=1)
  #g.plot(Gnuplot.Data(x,y,with='lines',title='Exact'),Gnuplot.Data(xgrid,EVs[j,:]/fac,title='Computed'),title='Hermite Elements')
f.close()
