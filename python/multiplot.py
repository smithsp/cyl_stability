#! /usr/pppl/python/2.3.5x/bin/python
import Gnuplot
def multiplot(d=[Gnuplot.Data([0,1],[0,1],with='lines')],ylabels=['x1'],fn="test.eps",title="notitle",fontsize=18,ysize=2.3,xsize=1.3):
  if len(d)>len(ylabels):
    raise ValueError('Length of ylabels is shorter than length of d')
  g = Gnuplot.Gnuplot()
  g('set term postscript enhanced eps color solid lw 2  %d'%(fontsize))
  g('set output "%s"'%fn)
  num = len(d)
  titled = True
  if title=="notitle":
    titled = False
  yh=(ysize-titled*.1)/float(num)
  xh=xsize
  #g('set key outside')
  g('set size %g,%g'%(xsize,ysize))
  g('set multiplot')
  g('set size %g,%g'%(xh,yh))
  for i in range(num):
    if i==num-1 and titled:
      g('set label "%s" at screen %g,%g center font "Helvetica,14"'%(title,xsize/2.,ysize-.08))
    g('set origin 0,%g'%(yh*i))
    g('set ylabel "%s"'%ylabels[i])  
    g._clear_queue()
    g._add_to_queue(d[i])
    g.refresh()
  g('unset multiplot')
  g('set output')
  g('set term x11')

if __name__=='__main__':
  d = []
  d.append([Gnuplot.Data([0,1],[0,1],with='lines'),Gnuplot.Data([1,0],[0,1],with='lines')])
  d.append([Gnuplot.Data([1,0],[0,1],with='lines')])
  multiplot(d,ylabels=['x1','x2','x3'],fn="test2.eps")
