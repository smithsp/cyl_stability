from get_EV import *
import sys
try:
  equilib = sys.argv[1]
  num = sys.argv[2]
  coord = sys.argv[3]
except:
  print 'You must give three command line parameters: equilib, id_num, coordinate_num'
  sys.exit()
data = get_data(equilib,num)
evals = get_evals(equilib,num,data)
inds = Numeric.nonzero(Numeric.absolute(evals.imag)>0)
if len(inds)==0:
  print 'There were no unstable modes.'
  sys.exit()
grid = get_grid(equilib,num,data)
j = 1
for i in Numeric.take(inds,Numeric.argsort(Numeric.take(evals.real,inds))):
  evec = get_evec(equilib,num,data,coord,i)
  evec_num = i
  g = Gnuplot.Gnuplot(persist=1)
  if j%2==1:
    g.ylabel('Real Part of Eigenfunction')
  else:
    g.ylabel('Imaginary Part of Eigenfunction')
  g.title('Run %s: Coordinate %s eigenfunction for {/Symbol w}=%g+%gi'%(num,coord,evals[evec_num].real,evals[evec_num].imag))
  g.plot(Gnuplot.Data(grid,evec,))
  g.hardcopy('plots/equilib%s/%s.evec%s_%d.eps'%(equilib,num,coord,evec_num))
  j = j+1
