import get_EV, Gnuplot
ref = sys.argv[1]
d1 = get_EV(ref)
var_key='ngrid'
dk = dict()
print_keys = d1.idata #['Vz0','P1','alpha','ngrid','tw','b']
dk_keys = ['t_assemb','num','t_solve',var_key,'epsilo','NN','vcyl']
for key in d1.data.keys():
  if key not in dk_keys:
    dk[key] = d1.data[key]
matching_files = []
var_list= []
minnum = 1
maxnum = 1400
if 'minn' in sys.argv:
  minnum = int(sys.argv[sys.argv.index('minn')+1])
if 'maxn' in sys.argv:
  maxnum = int(sys.argv[sys.argv.index('maxn')+1])
for i in range(minnum,maxnum,1):
  try:
    #print i,
    d2 = get_EV(i)
    data = d2.data
  except IOError:
    #print '%d, '%i,
    continue
  match = True
  for key in dk.keys(): 
    #if key == 'm':
      #if abs(data[key]) == abs(dk[key]):
        #continue   
    if data[key] != dk[key]:
      match = False
      #print 'dk does not equal data for run %d on key: %s'%(i,key)
      break
  if match:    
    matching_files.append(data['num'])
