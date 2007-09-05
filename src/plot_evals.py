from get_EV import *
import Gnuplot, Numeric, sys
try:
  equilib = sys.argv[1]
  num = sys.argv[2]
except:
  print 'You must enter two command line parameters: equilib and num'
  sys.exit()
data = get_data(equilib,num)
evals = get_evals(equilib,num,data)
g = Gnuplot.Gnuplot(persist=1)
g.plot(Gnuplot.Data(evals.real,evals.imag))
