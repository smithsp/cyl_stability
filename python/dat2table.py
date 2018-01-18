#! /usr/pppl/python/2.3.5x/bin/python
import get_EV,os
fout = file('table_data.xls','w')
fout.write(os.getcwd()+'\n')
ddict = dict()
for j in range(1,4500):
  try:
    d1 = get_EV.get_EV(j)
    try:
      ddict[d1.data['equilibr']].append(d1)
    except:
      ddict[d1.data['equilibr']]=[d1]
  except:
    continue
dkeys=ddict.keys()
dkeys.sort()
for i in dkeys:
  fout.write('Equilibrium %d \n'%(i))
  header = True
  for d1 in ddict[i]:
    if header:
      header = False
      fout.write('num \t ')
      for key in d1.idata:
        fout.write(key+'\t ')
      fout.write('\n')
    fout.write('%d'%d1.data['num'])
    for key in d1.idata:
      fout.write('\t %s'%(str(d1.data[key])))
    fout.write('\n')
  fout.write('\n')
fout.close()
    
