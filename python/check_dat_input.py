#! /usr/pppl/python/2.3.5x/bin/python
import get_EV,nl2dict,os,sys
nlf = 'input/'
datf = 'output_vcyl/'

dkeys = {'ALPHA':'alpha','AR':'a','BCROW':'BCrow','BR':'b',\
  'BT0':'Bt0','BZ0':'Bz0','EPS':'eps','EPSILO':'epsilo',\
  'EPSKVA':'epskVa','EPSVP':'epsVp','EPSVZ':'epsVz',\
  'EQUILIB':'equilibr','KAPPA':'kappa','KZ':'k',\
  'LAMBD':'lambd','LEND0':'Lend0','MT':'m','N':'ngrid',\
  'NU':'nu','P0':'P0','P1':'P1','RHO0':'rho0','RS':'rs',\
  'S2':'s2','TW':'tw','VP0':'Vp0','VZ0':'Vz0'}
stray_keys = []
for i in range(100,5000):
  nlfn = '%s%d.in'%(nlf,i)
  if os.path.exists(nlfn):
    dnl = nl2dict.nl2dict(nlfn)
  else:
    continue
  datfn = '%s%d.dat'%(datf,i)
  if os.path.exists(datfn):
    data = get_EV.get_EV(i).data
  else:
    continue
  for key in dnl.keys():
    if key not in dkeys.keys():
      stray_keys.append(key)
      #print '%s not in dkeys.keys()'%(str(key))
      continue
    try:
      nls = float(dnl[key])
      ds = (data[dkeys[key]])
    except:
      continue
    if nls!=ds:
      if key not in ['S2','KZ','EPS']:
        print 'Run %d does not match on key %s: %s vs %s'%(i,str(key),nls,ds)
  #sys.exit()
  
import sets
print 'The stray keys are ', sets.Set(stray_keys)
