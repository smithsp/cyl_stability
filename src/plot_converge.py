import Gnuplot, Numeric, math, sys
from pygsl import sf
try:
  min_N = (sys.argv[1])
  max_N = (sys.argv[2])
  delta_N = (sys.argv[3])
except:
  print 'You must provide 3 arguments: min_N, max_N, delta_N'
  sys.exit()
#print 'sys.argv=',sys.argv
f = file('spline_converge_%s_%s_%s.txt'%(min_N,max_N, delta_N))
m_Omega = []
Bt0 = []
kz = []
N = []
eps = 1.e-15
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g.xlabel('Number of grid points')
g.ylabel('Re({/Symbol w})')
#g('set yrange [-100:100]')
g.title('Convergence of Re({/Symbol w}) for Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
g1.xlabel('Number of grid points')
g1.ylabel('Im({/Symbol w})')
g1.title('Convergence of Im({/Symbol w}) for Appert coords in Guazzotto approach')
g2 = Gnuplot.Gnuplot(persist=1)
g2.xlabel('Re({/Symbol w})')
g2.ylabel('Im({/Symbol w})')
g2._clear_queue()
g2('set key outside')
g3 = Gnuplot.Gnuplot(persist=1)
g3.title('Slow branch convergence of Re({/Symbol w}) for Appert coords in Guazzotto approach')
g3.xlabel('Number of grid points')
g3.ylabel('Re({/Symbol w})')
minval = 10.
maxval = -10.
for k in range(100):
  try:
    N= int(f.readline())
    alfven1 = float(f.readline())
    alfven2 = float(f.readline())
    max_slow = float(f.readline())
    slow_cont1 = float(f.readline())
    slow_cont2 = float(f.readline())
    npts = int(f.readline())
  except:
    break
  evalsr = Numeric.zeros(npts,typecode=Numeric.Float)
  evalsi = Numeric.zeros(npts,typecode=Numeric.Float)
  for i in range(npts):
    evalsr[i] = float(f.readline())
  for i in range(npts):
    evalsi[i] = float(f.readline())
  f.readline() 
  inds = Numeric.nonzero(evalsi>0)
  if len(inds)>0:
    maxval = max(max(Numeric.take(evalsr,inds)),maxval)
    minval = min(min(Numeric.take(evalsr,inds)),minval) 
  dr = Gnuplot.Data([N]*npts,evalsr,with='points 1')
  di = Gnuplot.Data([N]*npts,evalsi,with='points 1')
  dri = Gnuplot.Data(evalsr,evalsi+k,with='points',title='N=%i'%(N))
  ds = Gnuplot.Data([N]*npts,(evalsr**2-slow_cont1**2)/slow_cont1**2,with='points 1')
  g._add_to_queue([dr])
  g1._add_to_queue([di])
  g2._add_to_queue([dri])
  g3._add_to_queue([ds])
#g.refresh()
g.hardcopy('bspline/spline_converge_real.eps',size=(1,1))
g('set yrange [%g:*]'%(alfven2))
g.hardcopy('bspline/spline_converge_real_fast.eps',size=(1,1))
N_range = [float(min_N),float(max_N)]
g._add_to_queue([Gnuplot.Data(N_range,[alfven1,alfven1],with='lines 2',title='Alfven1')])
g._add_to_queue([Gnuplot.Data(N_range,[alfven2,alfven2],with='lines 3',title='Alfven2')])
g('set yrange [%g:%g]'%(alfven1/1.001,alfven2*1.001))
g('set key outside')
g('unset logscale y')
g.hardcopy('bspline/spline_converge_real_alfven.eps',size=(1,1))
print slow_cont1,max_slow,slow_cont2
g3('set yrange [*:%g]'%((max_slow-slow_cont1)/slow_cont1*2.))
g3('set logscale y')
g3.hardcopy('bspline/spline_converge_real_slow.eps',size=(1,1))
g1.hardcopy('bspline/spline_converge_imag.eps',size=(1,1))
if maxval>0:
  g2('set xrange [%g:%g]'%(maxval/100.,maxval*10.))
#else:
  #g2('set xrange [%g:%g]'%(.0001,1000.))
g2('set logscale x')
g2.title('{/Symbol w} in complex plane')
#g2.refresh()
g2.hardcopy('bspline/spline_converge_complex.eps')
