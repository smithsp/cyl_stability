import Numeric, Gnuplot, sys
from get_EV import *
equilib = sys.argv[1]
num = sys.argv[2]
coord = sys.argv[3]
evec_num = int(sys.argv[4])
data = get_data(equilib,num) 
evals = get_evals(equilib,num,data)
#for i in range(len(evals
#print Numeric.take(evals.imag,(evals.imag<0))
grid = get_grid(equilib,num,data)
evec = get_evec(equilib,num,data,coord,evec_num)
g = Gnuplot.Gnuplot(persist=1)
g.title('Run %s: Coordinate %s eigenfunction for {/Symbol w}=%g+%gi'%(num,coord,evals[evec_num].real,evals[evec_num].imag))
g.plot(Gnuplot.Data(grid,evec,))
g.hardcopy('plots/equilib%s/%s.evec%s_%d.eps'%(equilib,num,coord,evec_num))
