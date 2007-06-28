import Gnuplot, Numeric, math, sys
from pygsl import sf
try:
  N = (sys.argv[1])
  delta_Bt0 = (sys.argv[2])
  max_Bt0 = (sys.argv[3])
  log_Bt0 = (sys.argv[4])
except:
  print 'You must provide 4 arguments: N, delta_Bt0, max_Bt0, log_Bt0'
  sys.exit()
f = file('spline%s_var_Bt0_%s_%s_%s.txt'%(N,delta_Bt0,max_Bt0,log_Bt0))
m_Omega = []
Bt0 = []
kz = []
eps = 1.e-15
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g.xlabel('Bt0')
g.ylabel('Re({/Symbol w})')
g('set yrange [-100:100]')
g.title('Re({/Symbol w}) for fixed body rotation with Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
g1.xlabel('Bt0')
g1.ylabel('Im({/Symbol w})')
g1.title('Im({/Symbol w}) for fixed body rotation with Appert coords in Guazzotto approach')
g2 = Gnuplot.Gnuplot(persist=1)
g2._clear_queue()
g2.xlabel('Re({/Symbol w})')
g2.ylabel('Im({/Symbol w})')
g2('set key outside')
minval = 10.
maxval = -10.
for k in range(100):
  try:
    Bt0.append(float(f.readline()))
    rho0 = float(f.readline())
    kz = float(f.readline())
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
  dr = Gnuplot.Data([Bt0[k]]*npts,evalsr,with='points 1')
  di = Gnuplot.Data([Bt0[k]]*npts,evalsi,with='points 1')
  dri = Gnuplot.Data(evalsr,evalsi+k,with='points',title='Bt0 = %g'%(Bt0[k]))
  g._add_to_queue([dr])
  g1._add_to_queue([di])
  g2._add_to_queue([dri])

g.replot(Gnuplot.Data(Bt0,m_Omega,title='Doppler shift (m*{/Symbol W})',with='lines 2'))
g.hardcopy('bspline/spline_var_Bt0_real.eps',size=(1,.50))
g1.refresh()
g1.hardcopy('bspline/spline_var_Bt0_imag.eps',size=(1,.50))
g2('set xrange [%g:%g]'%(minval/10.,maxval*10.))
g2('set logscale x')
g2.title('{/Symbol w} in complex plane with k=%g'%(kz))
g2.refresh()
g2.hardcopy('bspline/spline_var_Bt0_complex.eps')
