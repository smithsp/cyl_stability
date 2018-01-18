#! /usr/pppl/python/2.3.5x/bin/python
from finite_elements import *
import Gnuplot
import Numeric
import multiplot

N = 20
grid = Numeric.zeros([N],typecode=Numeric.Float)
for i in range(len(grid)):
  grid[i] = i/float(N-1)
x = Numeric.arange(0,1,(1/2.)**9)
d = []
ylabels = []

#This is to plot the integral finite elements
phi = []
p0=0.
Bt0=1.
Bz0=15.
a=1.
phi.append(linear_fa(grid[0],p3=grid[1],p0=p0,Bt0=Bt0,Bz0=Bz0,a=a))
for i in range(1,len(grid)-1):
  phi.append(linear_fa(grid[i],p2=grid[i-1],p3=grid[i+1],p0=p0,Bt0=Bt0,Bz0=Bz0,a=a))
phi.append(linear_fa(grid[-1],p2=grid[-2],p0=p0,Bt0=Bt0,Bz0=Bz0,a=a))

y = Numeric.zeros([len(x),len(phi)],typecode=Numeric.Float)
for j in range(len(phi)):
  for i in range(len(x)):
    y[i,j] = phi[j].val(x[i])
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
for j in range(len(phi)):
  g._add_to_queue([Gnuplot.Data(x,y[:,j],with='lines')])
d.append(g.itemlist)
g.title('Integral finite elements')
ylabels.append('Integral finite elements')
g.refresh()
g.hardcopy('plots/integral_finite_elements.eps',eps=True,color=True)

#This is to plot the integral finite elements' derivative
for j in range(len(phi)):
  for i in range(len(x)):
    y[i,j] = phi[j].val_prime(x[i])
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
for j in range(len(phi)):
  g._add_to_queue([Gnuplot.Data(x,y[:,j],with='lines')])
d.append(g.itemlist)
g.title('Integral finite elements derivative')
ylabels.append('Integral finite elements derivative')
g.refresh()
g.hardcopy('plots/integral_finite_elements_deriv.eps',eps=True,color=True)

#This is to plot the linear elements
phi = []
phi.append(linear(grid[0],p3=grid[1]))
for i in range(1,len(grid)-1):
  phi.append(linear(grid[i],p2=grid[i-1],p3=grid[i+1],))
phi.append(linear(grid[-1],p2=grid[-2]))

for j in range(len(phi)):
  for i in range(len(x)):
    y[i,j] = phi[j].val(x[i])
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
for j in range(len(phi)):
  g._add_to_queue([Gnuplot.Data(x,y[:,j],with='lines')])
d.append(g.itemlist)
g.title('Linear finite elements')
ylabels.append('Linear finite elements')
g.refresh()
g.hardcopy('plots/linear_finite_elements.eps',eps=True,color=True)

#This is to plot the linear elements' derivative
for j in range(len(phi)):
  for i in range(len(x)):
    y[i,j] = phi[j].val_prime(x[i])
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
for j in range(len(phi)):
  g._add_to_queue([Gnuplot.Data(x,y[:,j],with='lines')])
d.append(g.itemlist)
g.title('Linear finite elements derivative')
ylabels.append('Linear finite elements derivative')
g.refresh()
g.hardcopy('plots/linear_finite_elements_deriv.eps',eps=True,color=True)

#This is to plot the constant elements
phi = []
phi.append(const(grid[0],p3=grid[1]))
for i in range(1,len(grid)-1):
  phi.append(const(grid[i],p3=grid[i+1],))

y = Numeric.zeros([len(x),len(phi)],typecode=Numeric.Float)
for j in range(len(phi)):
  for i in range(len(x)):
    y[i,j] = phi[j].val(x[i])
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
for j in range(len(phi)):
  g._add_to_queue([Gnuplot.Data(x,y[:,j],with='lines')])
d.append(g.itemlist)
g.title('Constant finite elements')
ylabels.append('Constant finite elements')
g('set yrange [0:1.1]')
g.refresh()
g.hardcopy('plots/constant_finite_elements.eps',eps=True,color=True)

#This is to plot the constant elements' derivative (which should be 0)
for j in range(len(phi)):
  for i in range(len(x)):
    y[i,j] = phi[j].val_prime(x[i])
g = Gnuplot.Gnuplot(persist=1)
g._clear_queue()
for j in range(len(phi)):
  g._add_to_queue([Gnuplot.Data(x,y[:,j],with='lines')])
d.append(g.itemlist)
g.title('Constant finite elements derivative')
ylabels.append('Constant finite elements derivative')
g('set yrange [0:1.1]')
g.refresh()
g.hardcopy('plots/constant_finite_elements_deriv.eps',eps=True,color=True)

multiplot.multiplot(d,ylabels,title='Finite elements',fn='plots/finite_elements.eps')
