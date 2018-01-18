#! /usr/pppl/python/2.3.5x/bin/python
import get_EV,Numeric,Gnuplot,multiplot,sys,os,math

num = (sys.argv[1])
evec_nums = []
for i in range(2,len(sys.argv)):
  try:
    evec_nums.append(int(sys.argv[i]))
  except:
    break
N_fac = 10
EV = get_EV.get_EV(num) 
minx = 0
maxx = EV.data['a']
if 'xr' in sys.argv:
  ind = sys.argv.index('xr')
  minx = max(0,float(sys.argv[ind+1]))
  maxx = min(EV.data['a'],float(sys.argv[ind+2]))
x = Numeric.arange(minx,maxx,EV.data['a']/(EV.data['ngrid']*N_fac))
coords = ['1','2','3','r']
nc = len(coords)
ne = len(evec_nums)
xs = 1.
dx = xs/nc
ys = 1.
dy = ys/ne
g = Gnuplot.Gnuplot()
fn = 'plots/%s.evecs'%num
inds = Numeric.argsort(evec_nums)
for evec_num in Numeric.take(evec_nums,inds):
  fn = fn+'_%d'%(evec_num)
fn = fn+'.eps'
print 'Output to %s'%fn
g('set output "%s"'%fn)
g('set term postscript eps enhanced dl 2 ')
g('set multiplot')
g('set size %g,%g'%(dx,dy))
g('set xtics ( 0, "" 0.5,"1" 1)')
g('set xrange [0:1]')
for j in range(ne):
  evec_num = evec_nums[j]
  for i in range(nc):
    coord = coords[i]
    if j==ne-1:
      g('set title "{/Symbol x}_%s"'%coord)
    else:
      g('unset title')
    if i == 0:
      er=EV.evals[evec_num].real
      ei=EV.evals[evec_num].imag
      if ei>0:
        yl = '{/Symbol w}=%4.2g+%4.2g i'%(er,ei)
      elif ei<0:
        yl = '{/Symbol w}=%4.2g%4.2g i'%(er,ei)
      elif ei==0:
        yl = '{/Symbol w}=%4.2g+0i'%(er)
      g.ylabel(yl)
    else:
      g('unset ylabel')
    g('set origin %g,%g'%(i*dx,j*dy))
    (y,yl) = EV.get_evec_complex(coord,evec_num,x)
    dr = Gnuplot.Data(x,y.real,with='lines 1')
    di = Gnuplot.Data(x,y.imag,with='lines 3')
    g.plot(dr,di)
g('unset multiplot')
g('set term x11')
g('set output')
