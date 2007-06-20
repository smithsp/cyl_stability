import Gnuplot, Numeric, math
from pygsl import sf
f = file('spline_var_k.txt')
m_Omega = []
Bt0 = []
kz = []
eps = 1.e-15
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
#g('set logscale y')
g.xlabel('k')
g.ylabel('Re({/Symbol w})')
g('set yrange [-100:100]')
g.title('Illustration of Doppler shift with Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
#g1('set logscale y')
g1.xlabel('k')
g1.ylabel('Im({/Symbol w})')
g1.title('Im({/Symbol w}) for fixed body rotation with Appert coords in Guazzotto approach')
for k in range(20):
  try:
    Bt0.append(float(f.readline()))
    rho0 = float(f.readline())
    kz.append(float(f.readline()))
    m = int(f.readline())
    ar = float(f.readline())
    m_Omega.append(m*math.sqrt(2/rho0)/ar*Bt0[k])
    #if eps[k] == 0:
      #eps[k] = 9.e-7
    npts = int(f.readline())
    evalsr = Numeric.zeros(npts,typecode=Numeric.Float)
    evalsi = Numeric.zeros(npts,typecode=Numeric.Float)
    for i in range(npts):
      evalsr[i] = float(f.readline())
    for i in range(npts):
      evalsi[i] = float(f.readline())
    f.readline()
    dr = Gnuplot.Data([kz[k]]*npts,evalsr,with='points 1')
    di = Gnuplot.Data([kz[k]]*npts,evalsi,with='points 1')
    g._add_to_queue([dr])
    g1._add_to_queue([di])
  except:
    break
    #print eps,max_slow_inf
print 'kz=', kz
print 'Bt0=', Bt0
#g('set size 1.,.5')
g.replot(Gnuplot.Data(Bt0,m_Omega,title='Doppler shift (m*{/Symbol W})',with='lines 2'))
g.hardcopy('bspline/spline_var_k_real.eps',size=(1,.50))
#g1('set size 1.,.5')
g1.replot(Gnuplot.Data(Bt0,m_Omega,title='Doppler shift (m*{/Symbol W})',with='lines 2'))
g1.hardcopy('bspline/spline_var_k_imag.eps',size=(1,.50))
#R = sf.bessel_I1(kz*ar)[0]/(kz*ar*(sf.bessel_I1(kz*ar+eps)[0]-sf.bessel_I1(kz*ar)[0])/eps)
#G = m_Omega*(1-R)
#F = (m_Omega/m)**2*(m**2*(1-3*R)+2/R)
#print 'R = ',R
#print 'Bessel_I1(%g*%g) = '%(kz,ar),sf.bessel_I1(kz*ar)[0]
#print 'Bessel_I1_prime(%g*%g) = '%(kz,ar), (sf.bessel_I1(kz*ar+eps)[0]-sf.bessel_I1(kz*ar)[0])/eps
