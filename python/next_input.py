#! /usr/pppl/python/2.3.5x/bin/python
import sys,os
def next(ll=1,ul=2000):
  ld = os.listdir('./input')
  free = []
  for i in range(ll,ul):
    if '%d.in'%i not in ld:
      free.append(i)
  return free
  #raise 'There are no free input file numbers between %d and %d'%(ll,ul)
if __name__=='__main__':
  print 'here'
  print next(ll=int(sys.argv[1]),ul=int(sys.argv[2]))
