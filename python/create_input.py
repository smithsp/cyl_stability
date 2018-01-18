#! /usr/pppl/python/2.3.5x/bin/python
import sys, os

ref = int(raw_input('What is the reference run number? '))
start = int(raw_input('What is the first run number? '))
if ref==start:
  start = start+1
fin = int(raw_input('What is the last run number? '))
var = raw_input('Which variable do you want to change? ')
input_folder = 'input'#raw_input('In which folder should the input files be created? ')
print ''
flog = file('ref%d_%d_%d.log'%(ref,start,fin),'w')
flog.write(var)
for i in range(start,fin+1):
  val = (raw_input('Enter value for %s for run %d.  '%(var,i)))
  f = file('%s/%d.in'%(input_folder,i),'w')
  f.write('&cyl_params\n')
  f.write('%s=%s /'%(var,val))
  f.close()
  flog.write(' %d: %s '%(i,val))
f = file('control_params.in','w')
f.write('&control_params ref=%d, start=%d, fin=%d, verbose=.false./ \n'%(ref,start,fin))
f.close()
print ''
#print sys.stdin
