import Gnuplot, Numeric, math
f = file('spline_var_Bt0.txt')
m_Omega = []
Bt0 = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
#g('set logscale y')
g.xlabel('Vz0')
g.ylabel('Re({/Symbol w})')
g.title('Illustration of Doppler shift with Appert coords in Guazzotto approach')
g1 = Gnuplot.Gnuplot(persist=1)
g1._clear_queue()
#g1('set logscale y')
g1.xlabel('Vz0')
g1.ylabel('Im({/Symbol w})')
g1.title('Illustration of Doppler shift with Appert coords in Guazzotto approach')
for k in range(20):
  try:
    Bt0.append(float(f.readline()))
    rho0 = float(f.readline())
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
    dr = Gnuplot.Data([Vz0[k]]*npts,evalsr,with='points 1')
    di = Gnuplot.Data([Vz0[k]]*npts,evalsi,with='points 1')
    g._add_to_queue([dr])
    g1._add_to_queue([di])
  except:
    break
    #print eps,max_slow_inf
g.replot(Gnuplot.Data(Bt0,m_Omega,title='Doppler shift (m*{/Symbol W})',with='lines 2'))
g.hardcopy('bspline/spline_var_Vz0_real.eps')
g1.replot(Gnuplot.Data(Bt0,m_Omega,title='Doppler shift (m*{/Symbol W})',with='lines 2'))
g1.hardcopy('bspline/spline_var_Vz0_imag.eps')
