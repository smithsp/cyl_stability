#! /usr/pppl/python/2.3.5x/bin/python
import get_EV, Gnuplot, sys
num = sys.argv[1]
d = get_EV.get_EV(num)
print_keys = d.idata 
print_keys.insert(0,'num')
#sigma = -1j*(d.evals-d.data['m'])
sigma = d.get_numeric_evals() # This is only to calculate the necessary indices (d.numinds)
g = Gnuplot.Gnuplot(persist=1)
#g('set xrange [%g:%g]'%(d.data['m']-1,d.data['m']+1))
#g('set xrange [-3:7]')
#g('set yrange [-.11:.11]')
#g('set size .5,1')
g('set log y')
g.xlabel('-Re({/Symbol w})')
g.ylabel('-Im({/Symbol w})')
tt = ''
for i in range(len(print_keys)):
  tt = tt+'%s=%g, '%(print_keys[i],d.data[print_keys[i]])
g.title(tt[:-2])
dd = []
for i in d.numinds:
  dd.append(Gnuplot.Data(-d.evals[i].real,-d.evals[i].imag,with='points pt 1 ps %g'%(d.BCerror[i]**0)))
g.itemlist = dd
g.refresh()
g.hardcopy('plots/%d_evals_BCerr.eps'%d.data['num'],color=False)#,size=[.5,1])
dd = []
for i in d.numinds:
  dd.append(Gnuplot.Data(-d.evals[i].real,-d.evals[i].imag,with='points pt 1 ps %g'%(d.error[i]**0)))
g.itemlist = dd
g.refresh()
g.hardcopy('plots/%d_evals.eps'%d.data['num'],color=False)#,size=[.5,1])
