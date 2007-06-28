import Gnuplot, Numeric, sys
try:
  N = (sys.argv[1])
  delta_k = (sys.argv[2])
  max_k = (sys.argv[3])
except:
  print 'You must provide 3 arguments: N, delta_k, max_k'
  sys.exit()
f = file('spline%s_var_Vz0_%s_%s.txt'%(N,delta_k,max_k))
kV = []
Vz0 = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
#g('set logscale y')
g.xlabel('Vz0')
g.ylabel('({/Symbol w})')
g.title('Re({/Symbol w}): Illustration of Doppler shift with Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
#g('set logscale y')
g.xlabel('Vz0')
g1.ylabel('({/Symbol w})')
g1.title('Im({/Symbol w}): Appert coords in Guazzotto approach')
for k in range(20):
  try:
    Vz0.append(float(f.readline()))
    kV.append(float(f.readline())*Vz0[k])
    npts = int(f.readline())
    evalsr = Numeric.zeros(npts,typecode=Numeric.Float)
    evalsi = Numeric.zeros(npts,typecode=Numeric.Float)
    for i in range(npts):
      evalsr[i] = float(f.readline())
    for i in range(npts):
      evalsi[i] = float(f.readline())
    f.readline()
    dr = Gnuplot.Data([Vz0[k]]*npts,evalsr,with='points 1')
    di = Gnuplot.Data([Vz0[k]]*npts,evalsi,with='points 1')
    g._add_to_queue([dr])
    g1._add_to_queue([di])
  except:
    break
    #print eps,max_slow_inf
g.replot(Gnuplot.Data(Vz0,kV,title='Doppler shift (k*Vz)',with='lines 2'))
g.hardcopy('bspline/spline_var_Vz0_real.eps')
g1.refresh()
g1.hardcopy('bspline/spline_var_Vz0_imag.eps')
