#! /usr/pppl/python/2.3.5x/bin/python
import get_EV, Gnuplot, sys
num = sys.argv[1]
d = get_EV.get_EV(num)
sigma = -1j*(d.evals-d.data['m'])
g = Gnuplot.Gnuplot(persist=1)
#g('set xrange [%g:%g]'%(d.data['m']-1,d.data['m']+1))
#g('set xrange [-.1:.1]')
g('set yrange [-1:1]')
g('set size .5,1')
g.xlabel('Re({/Symbol s})')
g.ylabel('Im({/Symbol s})')
g.title('Run %d'%d.data['num'])
g.plot(Gnuplot.Data(sigma.real,sigma.imag,with='points 2'))
g.hardcopy('plots/%d_evals.eps'%d.data['num'],size=[.5,1])
