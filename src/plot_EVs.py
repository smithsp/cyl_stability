Bz = 2e-2
rho = 20
R = 1
k = .5
from pygsl import sf
import Numeric
import Gnuplot
m=1
gamma = ([1.841183781, 5.331442774, 8.536316366, 11.7060049, 14.86358863, 18.01552786, 21.16436986, 24.31132686, 27.45705057, 30.60192297, 33.7461829, 3.68899874E+01, 40.03344405, 43.17662897, 46.31959756, 49.46239114, 52.60504111, 55.74757179, 58.8900023])
N = 10

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
  #g.hardcopy('linear_const_ev_%d.eps'%j)
f.close()

f = file('spline_slow_EVs.txt')
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
for j in range(min(nphi,len(gamma)-1)):#nphi,2):
  for i in range(len(x)):
    y[i] = sf.bessel_J1(x[i]*gamma[j])[0]
  for i in range(len(xgrid)):
    ygrid[i] = sf.bessel_J1(xgrid[i]*gamma[j])[0]
  fac =(EVs[j,-1]/ygrid[-1])
  g = Gnuplot.Gnuplot(persist=1)
  g.plot(Gnuplot.Data(x,y,with='lines',title='Exact'),Gnuplot.Data(xgrid,EVs[j,:]/fac,title='Computed'),title='Spline Elements')
  #g.hardcopy('spline_ev_%d.eps'%j)
f.close()
import sys
sys.exit()
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
