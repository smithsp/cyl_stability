import Gnuplot, Numeric
f = file('spline_var_Vz0.txt')
kV = []
Vz0 = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
#g('set logscale y')
g.xlabel('Vz0')
g.ylabel('({/Symbol w})')
g.title('Illustration of Doppler shift with Appert coords in Guazzotto approach')
for k in range(20):
  try:
    Vz0.append(float(f.readline()))
    kV.append(float(f.readline())*Vz0[k])
    #if eps[k] == 0:
      #eps[k] = 9.e-7
    npts = int(f.readline())
    evals = Numeric.zeros(npts,typecode=Numeric.Float)
    for i in range(npts):
      evals[i] = float(f.readline())
    f.readline()
    d1 = Gnuplot.Data([Vz0[k]]*npts,evals,with='points 1')
    g._add_to_queue([d1])
  except:
    break
    #print eps,max_slow_inf
g.replot(Gnuplot.Data(Vz0,kV,title='Doppler shift (k*Vz)',with='lines 2'))
g.hardcopy('bspline/spline_var_Vz0.eps')
