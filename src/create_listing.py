from get_EV import *
for i in range(1,6):  
  f = file('equilib%d.list'%(i),'w')
  for j in range(1,1000):
    try:
      data = get_data('%d'%i,'%d'%j)
    except:
      break
    f.write(`data`)
    f.write('\n')
    
