#! /usr/pppl/python/2.3.5x/bin/python
import sys, Gnuplot, Numeric, sets
from get_EV import *
mn = int(sys.argv[1])
mx = int(sys.argv[2])
data_list = dict()
d1 = get_EV(mn)
data = d1.data
for key in data.keys():
  data_list[key] = [data[key]]
for i in range(mn+1,mx):
  try:
    try:
      d1 = get_EV(i)
      data = d1.data
    except:
      d1 = get_EV(i,output_folder='output_cyl')
      data = d1.data
  except:
    continue
  for key in data_list.keys():
    data_list[key].append(data[key])
for key in data_list.keys():
  if len(sets.Set(data_list[key])) ==1:
    continue
  if 't_' in key or 'type' in key or 'kB' in key or 'num' in key or 'NN' in key:
    continue
  #print key
  data_list[key] = Numeric.array(data_list[key])
  g = Gnuplot.Gnuplot(persist=1)
  #g('set log y')
  g.xlabel('Index number')
  g.ylabel(key)
  try:
    if max(data_list[key])/min(data_list[key]) > 100:
      g('set log y')
  except ZeroDivisionError:
    if max(data_list[key]) > 1e3:
      g('set log y')
  g.plot(Gnuplot.Data(data_list['num'],(data_list[key]),with='points'))
  print key,':'
  print data_list[key]
sys.exit()
if len(sets.Set(data_list['ngrid']))==1:
  sys.exit()  
g = Gnuplot.Gnuplot(persist=1)
g.xlabel('N')
g.ylabel('time(s)')
d1 = []
for key in ['t_solve','t_assemb']:
  d1.append(Gnuplot.Data(data_list['ngrid'],data_list[key],with='points',title=key))
#print d1
g.plot(*d1)
