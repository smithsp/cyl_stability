import Gnuplot, Numeric
f = file('spline_slow_vareps.txt')
max_slow_inf = []
eps = []
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
g('set logscale xy')
g.xlabel('epsilon')
g.ylabel('({/Symbol w}^2-{/Symbol w}_{/Symbol \245}^2)/{/Symbol w}_{/Symbol \245}^2')
g.title('Illustration of slow discrete->continuum transition')
for k in range(20):
  try:
    eps.append(float(f.readline()))
    max_slow_inf.append(float(f.readline()))
    if eps[k] == 0:
      eps[k] = 9.e-7
    npts = int(f.readline())
    evals = Numeric.zeros(npts,typecode=Numeric.Float)
    for i in range(npts):
      evals[i] = float(f.readline())
    f.readline()
    d1 = Gnuplot.Data([eps[k]]*npts,evals,with='points 1')
    g._add_to_queue([d1])
  except:
    break
    #print eps,max_slow_inf
g.replot(Gnuplot.Data(eps[1:],max_slow_inf[1:],title='max of continuum',with='lines 2'))
g.hardcopy('bspline/spline_slow_vareps.eps')
