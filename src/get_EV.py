import Numeric, Gnuplot, sys
def get_data(equilib,num):
  f = file('output/equilib%s/%s.dat'%(equilib,num))
  data = dict()
  data['ngrid'] = int(f.readline())
  data['NN'] = int(f.readline())
  data['m'] = int(f.readline())
  data['equilibr'] = int(f.readline())
  if int(equilib) <> data['equilibr']:
    print 'Error in file system.  path equilib=%d, .dat file equilib=%d'%(equilib,data['equilibr'])
  data['t_assemb'] = float(f.readline())
  data['t_solve'] = float(f.readline())
  data['k'] = float(f.readline())
  data['a'] = float(f.readline())
  data['Vz0'] = float(f.readline())
  data['epsVz'] = float(f.readline())
  data['rho0'] = float(f.readline())
  if equilib=='1' or equilib=='2' or equilib=='3' or equilib=='5':
    data['Bz0'] = float(f.readline())
    data['Bt0'] = float(f.readline())
    if equilib=='1' or equilib=='2':
      data['s2'] = float(f.readline())
      if equilib=='2':
        data['eps_rho'] = float(f.readline())
    elif equilib=='3' or equilib=='5':
      data['Vp0'] = float(f.readline())
      data['epsVp'] = float(f.readline())
      data['nz'] = int(f.readline())
  elif equilib=='4':
    data['P0'] = float(f.readline())
    data['P1'] = float(f.readline())
    data['lambd'] = float(f.readline())
  f.close()
  return data
def get_evals(equilib,num,data):
  if data['equilibr'] != int(equilib):
    print "Error: data['equilibr'] <> equilib in get_evals"
    raise ValueError
  fr = file('output/equilib%s/%s.evalsr'%(equilib,num))
  fi = file('output/equilib%s/%s.evalsi'%(equilib,num))
  evals = Numeric.zeros([data['NN']],typecode=Numeric.Complex)
  for i in range(data['NN']):
    evals[i] = float(fr.readline()) + 1j*float(fi.readline())
  fr.close()
  fi.close()
  return evals
def get_grid(equilib,num,data):
  if data['equilibr'] != int(equilib):
    print "Error: data['equilibr'] <> equilib in get_grid"
    raise ValueError
  grid = Numeric.zeros([data['ngrid']],typecode=Numeric.Float)
  fgrid = file('output/equilib%s/%s.grid'%(equilib,num))
  sgrid = (fgrid.readline())
  fgrid.close()
  sgrid = sgrid.strip()
  sgrid = sgrid.split(',')
  for i in range(data['ngrid']):
    grid[i] = float(sgrid[i])
  return grid
def get_evec(equilib,num,data,coord,evec_num):
  if data['equilibr'] != int(equilib):
    print "Error: data['equilibr'] <> equilib in get_evals"
    raise ValueError
  fevecs = file('output/equilib%s/%s.evecs%s'%(equilib,num,coord))
  sevecs = fevecs.readlines()
  sevec = sevecs[evec_num]
  sevec = sevec.strip()
  sevec = sevec.split(',')
  evec = Numeric.zeros([data['ngrid']],typecode=Numeric.Float)
  for i in range(data['ngrid']):
    evec[i] = float(sevec[i])
  return evec
if __name__ == '__main__':
  try:
    equilib = sys.argv[1]
    num = sys.argv[2]
    coord = sys.argv[3]
    evec_num = int(sys.argv[4])
  except:
    print 'You must give four command line parameters: equilib, id_num, coordinate_num, evec_num'
    sys.exit()
  data = get_data(equilib,num)
  evals = get_evals(equilib,num,data)
  #print Numeric.take(evals,Numeric.nonzero(Numeric.absolute(evals.imag)>0))
  grid = get_grid(equilib,num,data)
  evec = get_evec(equilib,num,data,coord,evec_num)
  g = Gnuplot.Gnuplot(persist=1)
  g.title('Run %s: Coordinate %s eigenfunction for {/Symbol w}=%g+%gi'%(num,coord,evals[evec_num].real,evals[evec_num].imag))
  g.plot(Gnuplot.Data(grid,evec,))
  filename = 'plots/equilib%s/%s.evec%s_%d.eps'%(equilib,num,coord,evec_num)
  g.hardcopy(filename)
  print 'Output written to %s'%(filename)
