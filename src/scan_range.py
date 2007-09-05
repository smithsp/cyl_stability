import sys, Gnuplot, Numeric, sets
from get_EV import *
equilib = sys.argv[1]
data_list = dict()
data = get_data(equilib,'%d'%(1))
for key in data.keys():
  data_list[key] = [data[key]]
for i in range(2,1000):
  try:
    data = get_data(equilib,'%d'%i)
  except:
    break
  for key in data_list.keys():
    data_list[key].append(data[key])
for key in data_list.keys():
  if len(sets.Set(data_list[key])) ==1:
    continue
  if 't_' in key:
    continue
  data_list[key] = Numeric.array(data_list[key])
  g = Gnuplot.Gnuplot(persist=1)
  g('set log y')
  g.xlabel('Index number')
  g.ylabel(key)
  g.plot(Gnuplot.Data(range(1,len(data_list[key])+1),Numeric.absolute(data_list[key]),with='points'))

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
