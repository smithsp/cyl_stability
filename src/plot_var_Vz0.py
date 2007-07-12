import Gnuplot, Numeric, sys
try:
  N = sys.argv[1]
  delta_k = sys.argv[2]
  max_k = sys.argv[3]
  log_k = sys.argv[4]
except:
  print 'You must provide 4 arguments: N, delta_Vz0, max_Vz0, log_Vz0'
  sys.exit()
f = file('spline%s_var_Vz0_%s_%s_%s.txt'%(N,delta_k,max_k,log_k))
kV = []
Vz0 = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
#g('set logscale y')
g.xlabel('Vz0')
g.ylabel('Re({/Symbol w})')
g.title('Re({/Symbol w}): Illustration of Doppler shift with Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
#g('set logscale y')
g.xlabel('Vz0')
g1.ylabel('Im({/Symbol w})')
g1.title('Im({/Symbol w}): Appert coords in Guazzotto approach')
g2 = Gnuplot.Gnuplot(persist=1)
g2.xlabel('Re({/Symbol w})')
g2.ylabel('Im({/Symbol w})')
g2._clear_queue()
g2('set key outside')
for k in range(20):
  try:
    Vz0.append(float(f.readline()))
    kV.append(float(f.readline())*Vz0[k])
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
  dr = Gnuplot.Data([Vz0[k]]*npts,evalsr,with='points 1')
  #print Numeric.sum(evalsr == 0)
  di = Gnuplot.Data([Vz0[k]]*npts,evalsi,with='points 1')
  dri = Gnuplot.Data(evalsr,evalsi+0*k*.1,with='points',title='Vz0=%g'%(Vz0[k]))
  g._add_to_queue([dr])
  g1._add_to_queue([di])
  g2._add_to_queue([dri])
  #print eps,max_slow_inf
g._add_to_queue([Gnuplot.Data(Vz0,kV,title='Doppler shift (k*Vz)',with='lines 2')])
#g('set logscale y')
#g('set yrange [1e-14:1000]')
g.hardcopy('bspline/spline_var_Vz0_real.eps')
g1('set logscale y')
#g1.refresh()
g1.hardcopy('bspline/spline_var_Vz0_imag.eps')
g2('set logscale x')
g2.title('{/Symbol w} in complex plane')
g2.hardcopy('bspline/spline_var_Vz0_complex.eps')
