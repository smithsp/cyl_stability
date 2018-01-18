#! /usr/pppl/python/2.3.5x/bin/python
import get_EV,Gnuplot,sys,Numeric
g = Gnuplot.Gnuplot()
g('set xlabel "b"')
g('set ylabel "r_s"')
g.itemlist = [Gnuplot.File("b_zero_rs.out1", using=" 1:2",title=None)]
#g.plot(
#g('plot ')
fn = 'plots/1068_rs_b.eps'
g.hardcopy(fn,eps=True,color=False,solid=False)
