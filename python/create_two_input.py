#! /usr/pppl/python/2.3.5x/bin/python
import sys, os

ref = int(raw_input('What is the reference run number? '))
start = int(raw_input('What is the first run number? '))
if ref==start:
  start = start+1
var1 = raw_input('What is the 1st variable you want to vary? ')
var2 = raw_input('What is the 2nd variable you want to vary? ')
ref_var1 = raw_input('What is the value of %s for run %d? '%(var1,ref))
ref_var2 = raw_input('What is the value of %s for run %d? '%(var2,ref))
nvar1 = int(raw_input('How many values of %s do you want to enter? '%var1))
nvar2 = int(raw_input('How many values of %s do you want to enter? '%var2))
input_folder = 'input'
print ''
ntot = nvar1*nvar2
fin = start+ntot
flog = file('ref%d_%d_%d.log'%(ref,start,fin),'w')
flog.write(var1)
flog.write(var2)
var1val = []
var2val = []
for i in range(0,nvar1):
  var1val.append((raw_input('val %s %i? '%(var1,i+1))))
for i in range(0,nvar2):
  var2val.append((raw_input('val %s %i? '%(var2,i+1))))
for i in range(0,ntot):
  j = i/nvar2
  k = i%nvar2
  #print j,k
  if var1val[j]==ref_var1 and var2val[k]==ref_var2:
    continue
  f = file('%s/%d.in'%(input_folder,i+start),'w')
  f.write('&cyl_params\n')
  f.write('%s=%s,%s=%s /'%(var1,var1val[j],var2,var2val[k]))
  f.close()
  flog.write(' %d: %s=%s,%s=%s '%(i,var1,var1val[j],var2,var2val[k]))
f = file('control_params.in','w')
f.write('&control_params ref=%d, start=%d, fin=%d, verbose=.false./ \n'%(ref,start,fin-1))
f.close()
print ''
