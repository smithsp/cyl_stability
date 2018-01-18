#! /usr/pppl/python/2.3.5x/bin/python
import sys, os, Numeric
flog = file('multi_input.log','a')
flog.write('---------------------------------------\n')
ref = int(raw_input('What is the reference run number? '))
flog.write('%d\n'%ref)
input_folder = 'input'#raw_input('In which folder should the input files be created? ')
ref_fn = '%s/%d.in'%(input_folder,ref)
if not os.path.exists(ref_fn):
  print """
  ERROR: The reference file did not exist.
  Please enter a reference run for which the input file already exists.
  """
  sys.exit()
start = int(raw_input('What is the first run number? '))
flog.write('%d\n'%start)
if ref==start:
  start = start+1
if os.path.exists('%s/%d.in'%(input_folder,start)):
  print 'The starting run number is already in use.'
  ans = raw_input('Do you want to overwrite it? [y/n] ')
  if ans not in ['Y','y']:
    sys.exit()
var = []
for i in range(0,10):
  var.append(raw_input('Which variable do you want to change? '))
  flog.write('%s\n'%var[i])
  if len(var[i])==0:
    var.pop(len(var)-1)
    break
print ''
ref_file = file(ref_fn,'r')
ref_lines = ref_file.readlines()
ref_dict = dict()
ref_vals = ((ref_lines[0].split(' ')[2]).split('/')[0]).split(',')
for i in range(len(ref_vals)):
  ref_dict[ref_vals[i].split('=')[0]]=ref_vals[i].split('=')[1]
if 'auto' in sys.argv:
  totn = 1
  val = []
  for j in range(len(var)):
    val.append([])
  for j in range(len(var)):
    nv = int(raw_input('How many values of %s do you want? '%(var[j])))
    flog.write('%d\n'%nv)
    totn = totn*nv
    for i in range(nv):
      val[j].append(raw_input('Enter value %d of %s. '%(i+1,var[j])))
      flog.write('%s\n'%val[j][i])
  print val
  flog.close()
  #sys.exit()
  fin = start+totn-1
  print fin
  
  indlen = []
  for j in range(len(var)):
    indlen.append(len(val[j]))
  print indlen
  inds = list(Numeric.argsort(indlen))
  print inds
  var_sort = []
  val_sort = []
  indlen = list(Numeric.take(indlen,inds))
  for ind in inds:
    var_sort.append(var[ind])
    val_sort.append(val[ind])
    indlen
  val = val_sort
  var = var_sort
  print var
  print val
  indarray = list(Numeric.indices(indlen))
  for j in range(len(indarray)):
    indarray[j] = Numeric.ravel(indarray[j])
  for i in range(start,fin+1):
    #print i
    for j in range(len(var)):
      #print (indarray[j])
      ind = indarray[j][i-start]
      ref_dict[var[j].upper()]=val[j][ind]
    f = file('%s/%d.in'%(input_folder,i),'w')
    f.write('&cyl_params ')
    key_str=''
    for key in ref_dict.keys():
      key_str=key_str+'%s=%s,'%(key,ref_dict[key])
    key_str=key_str[:-1]+' /\n'
    f.write(key_str)
    f.close()
  f = file('control_params.in','w')
  f.write('&control_params ref=%d, start=%d, fin=%d, verbose=.false./ \n'%(start,start,fin))
  f.close()
else:
  fin = int(raw_input('What is the last run number? '))
  #flog = file('ref%d_%d_%d.log'%(ref,start,fin),'w')
  #flog.write(var)
  #ref_dict['VZ0']='0.8'
  #ref_dict['N']='20'
  #ref_dict['ALPHA']='10.0'
  #ref_dict['RS']='0.0'
  for i in range(start,fin+1):
    val = []
    for j in range(len(var)):
      val.append(raw_input('Enter value for %s for run %d.  '%(var[j],i)))
    for j in range(len(var)):
      ref_dict[var[j].upper()]=val[j]
    f = file('%s/%d.in'%(input_folder,i),'w')
    f.write('&cyl_params ')
    key_str=''
    for key in ref_dict.keys():
      key_str=key_str+'%s=%s,'%(key,ref_dict[key])
    key_str=key_str[:-1]+' /'
    f.write(key_str)
    f.close()
    #flog.write(' %d: %s '%(i,val))
  f = file('control_params.in','w')
  f.write('&control_params ref=%d, start=%d, fin=%d, verbose=.false./ \n'%(start,start,fin))
  f.close()
print ''
#print sys.stdin
