#! /usr/pppl/python/2.3.5x/bin/python
import get_EV, Gnuplot, Numeric
n = []
ta = []
ts = []
v = []
for i in range(1,2300):
  try:
    d = get_EV.get_EV(i)
  except:
    continue
  ta.append(d.data['t_assemb'])
  ts.append(d.data['t_solve'])
  n.append(d.data['NN'])
  v.append(d.data['Vz0'])
g = Gnuplot.Gnuplot()
g.title('size of matrix (NNxNN) vs time(s) to solve')
g.ylabel('Time (s)')
g.xlabel('NN')
g('set log xy')
g('set key bottom right')
inds = Numeric.nonzero(Numeric.array(v)==0)
g.plot(Gnuplot.Data(Numeric.take(n,inds),Numeric.take(ts,inds),title='Vz0=0'))
inds = Numeric.nonzero(Numeric.array(v)!=0)
g._add_to_queue([Gnuplot.Data(Numeric.take(n,inds),Numeric.take(ts,inds),title='Vz0<>0')])
for i in range(1,5):
  g.replot(Gnuplot.Data(n,Numeric.array(n,typecode=Numeric.Float)**i*ts[-3]/n[-3]**i,title='NN^%d'%i,with='lines'))
g.hardcopy('plots/NN_tsolve.eps',eps=True)
g.title('size of matrix (NNxNN) vs time (s) to assemble')
g.ylabel('Time (s)')
g.xlabel('NN')
g('set log xy')
g('set key bottom right')
g.plot(Gnuplot.Data(n,ta,title='Data'))
for i in range(1,5):
  g.replot(Gnuplot.Data(n,Numeric.array(n,typecode=Numeric.Float)**(i/2.)*ta[-3]/n[-3]**(i/2.),title='NN^{%g}'%(i/2.),with='lines'))
g.hardcopy('plots/NN_tassemb.eps',eps=True)
print ts[-1], ta[-1]
