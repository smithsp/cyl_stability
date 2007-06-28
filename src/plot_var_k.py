import Gnuplot, Numeric, math, sys
from pygsl import sf
try:
  N = (sys.argv[1])
  delta_k = (sys.argv[2])
  max_k = (sys.argv[3])
  log_k = (sys.argv[4])
except:
  print 'You must provide 4 arguments: N, delta_k, max_k, log_k'
  sys.exit()
#print 'sys.argv=',sys.argv
f = file('spline%s_var_k_%s_%s_%s.txt'%(N,delta_k,max_k,log_k))
m_Omega = []
Bt0 = []
kz = []
eps = 1.e-15
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g.xlabel('k')
g.ylabel('Re({/Symbol w})')
g('set yrange [-100:100]')
g.title('Illustration of Doppler shift with Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
g1.xlabel('k')
g1.ylabel('Im({/Symbol w})')
g1.title('Im({/Symbol w}) for fixed body rotation with Appert coords in Guazzotto approach')
g2 = Gnuplot.Gnuplot(persist=1)
g2.xlabel('Re({/Symbol w})')
g2.ylabel('Im({/Symbol w})')
g2._clear_queue()
g2('set key outside')
minval = 10.
maxval = -10.
for k in range(100):
  try:
    Bt0.append(float(f.readline()))
    rho0 = float(f.readline())
    kz.append(float(f.readline()))
    m = int(f.readline())
    ar = float(f.readline())
  except:
    break
  m_Omega.append(m*math.sqrt(2/rho0)/ar*Bt0[k])
  npts = int(f.readline())
  evalsr = Numeric.zeros(npts,typecode=Numeric.Float)
  evalsi = Numeric.zeros(npts,typecode=Numeric.Float)
  for i in range(npts):
    evalsr[i] = float(f.readline())
  for i in range(npts):
    evalsi[i] = float(f.readline())
  f.readline() 
  inds = Numeric.nonzero(evalsi>1)
  if len(inds)>0:
    maxval = max(max(Numeric.take(evalsr,inds)),maxval)
    minval = min(min(Numeric.take(evalsr,inds)),minval) 
  dr = Gnuplot.Data([kz[k]]*npts,evalsr,with='points 1')
  di = Gnuplot.Data([kz[k]]*npts,evalsi,with='points 1')
  dri = Gnuplot.Data(evalsr,evalsi+k,with='points',title='k=%g'%(kz[k]))
  g._add_to_queue([dr])
  g1._add_to_queue([di])
  g2._add_to_queue([dri])

g.replot(Gnuplot.Data(Bt0,m_Omega,title='Doppler shift (m*{/Symbol W})',with='lines 2'))
g.hardcopy('bspline/spline_var_k_real.eps',size=(1,.50))
g1.refresh()
g1.hardcopy('bspline/spline_var_k_imag.eps',size=(1,.50))
g2('set xrange [%g:%g]'%(minval/10.,maxval*10.))
g2('set logscale x')
g2.title('{/Symbol w} in complex plane with Bt0=%g'%(Bt0[0]))
g2.refresh()
g2.hardcopy('bspline/spline_var_k_complex.eps')
