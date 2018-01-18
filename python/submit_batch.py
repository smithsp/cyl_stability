#! /usr/pppl/python/2.3.5x/bin/python
import sys, os
f = file('control_params.in')
lines = f.readlines()
f.close()
llist = lines[0].split()
#print '#PBS -N vcyl%s_%s'%(llist[2][6:-1],llist[3][4:-1])

f = file('batch_job')
bj = f.readlines()
f.close()
bj[0] = '#PBS -N vcyl_%s_%s\n'%(llist[2][6:-1],llist[3][4:-1])
start = int(llist[2][6:-1])
fin = int(llist[3][4:-1])
numprocs=fin-start+1
if 'serial' in sys.argv:
  numprocs = 1
if numprocs>20:
  numperproc = numprocs/20
  numprocs = numprocs/numperproc
if numprocs>1:
  bj[3] = '#PBS -l nodes=%d:ppn=1\n'%(numprocs)
  #bj[5] = '#PBS -q kestrel\n'
bj[30] = 'mpirun -np ${NPROCS} ./vcyl_parallel.exe\n'
f = file('batch_job1','w')
f.writelines(bj)
f.close()
os.system('qsub batch_job1')
