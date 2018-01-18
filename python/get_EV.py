#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, Gnuplot, sys, finite_elements, math, sets
class get_EV:
  
  def __init__(self, num, output_folder='output_vcyl'):
    self.num = num
    self.output_folder = output_folder
    self.data = self.get_data(num)
    self.idata = self.get_idata()
    self.evals = self.get_evals()
    try:
      self.BCerror = self.get_BCerror()
      self.error = self.get_error()
    except:
      pass
#      print 'No error data for run %s'%(num)
    self.imag = False
    try:
      self.grid = self.get_grid()
    except (IOError, ValueError):
      pass
  def get_data(self,num):
    f = file(self.output_folder + '/%s.dat'%(num))
    data = dict()
    data['fe_type'] = f.readline()[0:-1]
    data['ngrid'] = int(f.readline())
    data['NN'] = int(f.readline())
    data['m'] = int(f.readline())
    data['equilibr'] = int(f.readline())
    data['num'] = int(f.readline())
    data['BCrow'] = int(f.readline())
    data['epsilo'] = float(f.readline())
    data['t_assemb'] = float(f.readline())
    data['t_solve'] = float(f.readline())
    data['k'] = float(f.readline())
    data['a'] = float(f.readline())
    data['b'] = float(f.readline())
    data['tw'] = float(f.readline())
    data['rho0'] = float(f.readline())
    data['Bz0'] = float(f.readline())
    data['Bt0'] = float(f.readline())
    data['s2'] = float(f.readline())
    data['eps'] = float(f.readline())
    data['P0'] = float(f.readline())
    data['P1'] = float(f.readline())
    data['lambd'] = float(f.readline())
    data['Vz0'] = float(f.readline())
    data['epsVz'] = float(f.readline())
    data['Vp0'] = float(f.readline())
    data['epsVp'] = float(f.readline())
    data['rs'] = float(f.readline())
    data['alpha'] = float(f.readline())
    data['epskVa'] = float(f.readline())
    data['nu'] = float(f.readline())
    data['Lend0'] = f.readline()[0:-1]
    if 'T' in data['Lend0']:
      data['Lend0']=1
    else:
      data['Lend0']=0
    data['vcyl'] = f.readline()[0:-1]
    if 'T' in data['vcyl']:
      data['vcyl']=1
    else:
      data['vcyl']=0
    try:
      data['kappa'] = float(f.readline())
    except:
      data['kappa'] = 0.
    if data['equilibr'] == 6:
      data['Bt0'] = (data['P0']*data['eps'])**0.5
    f.close()
    return data
  def get_idata(self):
    basedata = ['ngrid','m','k','b','tw','alpha','rs','nu', 'kappa']
    idata = basedata[0:5]
    if   self.data['equilibr'] == 1:
      idata = idata + ['rho0','s2','Bz0']
    elif self.data['equilibr'] == 2:
      idata = idata + ['rho0','eps','s2','Bz0']
    elif self.data['equilibr'] == 3:
      idata = idata + ['rho0','Bt0','Bz0','Vz0']
    elif self.data['equilibr'] == 10:
      idata = idata + ['rho0','P0','eps','Bz0','Bt0','Vz0','epsVz']
    elif self.data['equilibr']==12:
      idata = idata + ['rho0','P0','eps','Bz0','Bt0','Vp0']
    elif self.data['equilibr']==13:
      idata = idata + ['Bz0','Bt0','Vz0','epsVz','eps']
    idata = idata + basedata[5:]
    return idata
  def get_evals(self):
    data = self.data
    fr = file(self.output_folder + '/%d.evalsr'%(data['num']))
    fi = file(self.output_folder + '/%d.evalsi'%(data['num']))
    evals = Numeric.zeros([data['NN']],typecode=Numeric.Complex)
    for i in range(data['NN']):
      evals[i] = float(fr.readline()) + 1j*float(fi.readline())
    fr.close()
    fi.close()
    return evals
  def get_error(self):
    f = file(self.output_folder + '/%d.err'%(self.data['num']))
    error = Numeric.zeros([self.data['NN']],typecode=Numeric.Float)
    for i in range(len(error)):
      temp = f.readline()
      if '*' in temp:
        error[i] = 5
      else:
        error[i] = float(temp)
    f.close()
    return error
  def get_BCerror(self):
    f = file(self.output_folder + '/%d.BCerr'%(self.data['num']))
    error = Numeric.zeros([self.data['NN']],typecode=Numeric.Float)
    for i in range(len(error)):
      temp = f.readline()
      if '*' in temp:
        error[i] = 5
      else:
        error[i] = float(temp)
    f.close()
    return error
  def get_grid(self):
    data = self.data
    grid = Numeric.zeros([data['ngrid']],typecode=Numeric.Float)
    fgrid = file(self.output_folder + '/%d.grid'%(data['num']))
    sgrid = (fgrid.readline())
    fgrid.close()
    sgrid = sgrid.strip()
    sgrid = sgrid.split(',')
    for i in range(data['ngrid']):
      grid[i] = float(sgrid[i])
    return grid
  def get_phi(self,deriv=False,num=1):
    if self.data['vcyl']:
      nphi = self.data['NN']/6
    else:
      nphi = self.data['NN']/3
    phi = []
    grid = self.grid
    N = self.data['ngrid']
    if self.data['fe_type']=='spline':
      phi0 = finite_elements.bspline(grid[0],p3=grid[1],deriv=deriv)
      phi1 = finite_elements.bspline(grid[0],p3=grid[1],p4=grid[2])
      phi2 = finite_elements.bspline(grid[0],p3=grid[1],p4=grid[2])
      phi3 = finite_elements.bspline(grid[1],p2=grid[0],p3=grid[2],p4=grid[3])
      phi1 = phi1*(phi0.C[1]/phi1.C[1])-phi0
      phi2 = phi2*(phi0.C[2]/phi2.C[2])-phi0
      phi3 = phi3*(phi0.C[2]/phi3.B[2])
      phi3.B[:] = phi3.B-phi0.C
      if num==1:
        phi.append(phi0)
        phi.append(finite_elements.bspline(grid[0],p3=grid[1],p4=grid[2],deriv=deriv))
        phi.append(finite_elements.bspline(grid[1],p2=grid[0],p3=grid[2],p4=grid[3],deriv=deriv))
      elif num==2:
        phi2.deriv=True
        phi3.deriv=True
        phi.append(phi0)
        phi.append(finite_elements.bspline(grid[0],p3=grid[1],p4=grid[2],deriv=deriv))
        #phi.append(phi2)
        phi.append(finite_elements.bspline(grid[1],p2=grid[0],p3=grid[2],p4=grid[3],deriv=deriv))
        #phi.append(phi3)
      elif num==3:
        phi1.deriv=True
        #print phi0.dx
        #phi0.C = Numeric.array([0.,0.,-3.*phi0.dx[2],2.])
        i=2
        phi2shift = finite_elements.bspline(grid[1],p2=grid[0],p3=grid[2],p4=grid[3],deriv=deriv)
        #print phi2shift.A,phi2shift.B,phi2shift.C,phi2shift.D
        phi2shift.shift()
        #print phi2shift.A
        phi0 = finite_elements.bspline(grid[i],p1=grid[i-2],p2=grid[i-1],p3=grid[i+1],p4=grid[i+2],deriv=deriv)+phi2shift
        phi.append(phi0)
        phi.append(phi1)
        #phi.append(finite_elements.bspline(grid[0],p3=grid[1],p4=grid[2],deriv=deriv))
        phi.append(finite_elements.bspline(grid[1],p2=grid[0],p3=grid[2],p4=grid[3],deriv=deriv))
      for i in range(2,N-2):
        phi.append(finite_elements.bspline(grid[i],p1=grid[i-2],p2=grid[i-1],p3=grid[i+1],p4=grid[i+2],deriv=deriv))
      phi.append(finite_elements.bspline(grid[N-2],p1=grid[N-4],p2=grid[N-3],p3=grid[N-1],deriv=deriv))
      phi.append(finite_elements.bspline(grid[N-1],p1=grid[N-3],p2=grid[N-2],deriv=deriv))
      if nphi > self.data['ngrid']:
        if self.data['a'] < self.data['b']:
          phi.append(finite_elements.bspline(grid[N-1],p2=grid[N-2],deriv=deriv))
      return phi
    else:
      if num==3:
        phi.append(finite_elements.const(grid[0],p3=grid[1]))
        for i in range(N-1):
          phi.append(finite_elements.const(grid[i],p3=grid[i+1]))
        return phi
      elif num==1:
        p0=self.data['P0']
        Bt0=self.data['Bt0']
        Bz0=self.data['Bz0']
        a=self.data['a']
        phi.append(finite_elements.linear_fa(grid[0],p3=grid[1],p0=p0,Bt0=Bt0,Bz0=Bz0,a=a))
        for i in range(1,N-1):
          phi.append(finite_elements.linear_fa(grid[i],p2=grid[i-1],p3=grid[i+1],p0=p0,Bt0=Bt0,Bz0=Bz0,a=a))
        phi.append(finite_elements.linear_fa(grid[N-1],p2=grid[N-2],p0=p0,Bt0=Bt0,Bz0=Bz0,a=a))
        return phi
      else:
        phi.append(finite_elements.linear(grid[0],p3=grid[1]))
        for i in range(1,N-1):
          phi.append(finite_elements.linear(grid[i],p2=grid[i-1],p3=grid[i+1]))
        phi.append(finite_elements.linear(grid[N-1],p2=grid[N-2]))
        return phi
  def get_sevecs(self,coord):
    try:
      sevecs = self.sevecs
    except:
      sevecs = []
      for i in range(1,4):
        fn = self.output_folder + '/%d.evecs%d'%(self.data['num'],i)
        if self.imag:
          fn = self.output_folder + '/%d.evecs_imag%d'%(self.data['num'],i)
        fevecs = file(fn)
        sevecs.append(fevecs.readlines())
      self.sevecs = sevecs
    return sevecs[coord]
  def get_evec(self,coord,evec_num):
    sevecs = self.get_sevecs(coord)
    sevec = sevecs[evec_num]
    sevec = sevec.strip()
    sevec = sevec.split(',')
    evec = Numeric.zeros([len(sevec)-1],typecode=Numeric.Float)
    for i in range(len(sevec)-1):
      evec[i] = float(sevec[i])
    return evec
  def get_evec_val(self,coord,evec_num,x,prime=False):
    coord = coord-1  #-1 is necessary for phi and sevecs lists
    evec = self.get_evec(coord,evec_num)
    try:
      phi = self.phi
    except:
      self.phi = [self.get_phi(deriv=False,num=1)]
      self.phi.append(self.get_phi(deriv=False,num=2))
      self.phi.append(self.get_phi(deriv=True,num=3))
    if len(self.phi[coord])!=len(evec):
      print 'You have the wrong length of phi for coord ',coord
    if self.data['vcyl']:
      coord = divmod(coord,3)[1] #This is to get the right finite elements of the phi list
    y = Numeric.zeros([len(x)],typecode=Numeric.Float)
    if prime:
      for i in range(len(x)):
        y[i] = Numeric.sum([evec[j]*self.phi[coord][j].val_prime(x[i]) for j in range(len(evec))])
    else:
      for i in range(len(x)):
        y[i] = Numeric.sum([evec[j]*self.phi[coord][j].val(x[i]) for j in range(len(evec))])
    return y
  def get_evec_r_val(self,evec_num,x):
    y = (self.get_evec_val(1,evec_num,x)+self.data['m']*self.get_evec_val(2,evec_num,x))
    for i in range(len(x)):
      if x[i]!=0:
        y[i] = y[i]/x[i]
      else:
        y[i] = 0
    return y
  def get_EVrot(self,minr=None,maxr=None):
    try:
      return self.EVrot
    except:
      #print minr, maxr
      if minr==None:
        minr = -1e7
      if maxr==None:
        maxr = 1e7
      try:
        return self.EVrot
      except:
        evals = self.get_unstable_evals()
        inds = self.get_unst_inds()
        rwm_inds = Numeric.nonzero((evals.real>minr)*(evals.real<maxr))
        try:
          del self.sevecs
        except:
          pass
        self.imag=False
        #print self.evals[inds[rwm_inds[0]]]
        #print inds[rwm_inds[0]]
        y = self.get_evec_r_val(inds[rwm_inds[0]],[self.data['a']])
        self.imag = True
        del self.sevecs  
        try:
          y1 = self.get_evec_r_val(inds[rwm_inds[0]],[self.data['a']])
          y2 = (y+1j*y1).real
          y3 = (y+1j*y1).imag
          #print y[-1],y1[-1],y2[-1],y3[-1]
          y = y2
          y1 = y3
        except:
          try:
            y = y.real
            y1 = y.imag
          except:
            y = y.real
            y1 = 0
        self.imag = False
        self.EVrot = 1.-1j*y1/y
        self.EVnorm = 1./(self.EVrot*(y+y1*1j)).real
        return self.EVrot
  def get_EVnorm(self,min_r=-1e7,max_r=1e7):
    try:
      return self.EVnorm
    except:
      temp = self.get_EVrot(minr=min_r,maxr=max_r)
      return self.EVnorm
  def get_evec_theta_val(self,evec_num,x):
    return (self.get_evec_val(2,evec_num,x,prime=True)+self.Bt(x)*self.get_evec_val(3,evec_num,x))
  def get_evec_z_val(self,evec_num,x):
    return (-self.Bt(x)/self.Bz(x)*self.get_evec_val(2,evec_num,x,prime=True)\
      +self.Bz(x)*self.get_evec_val(3,evec_num,x))
  def get_evec_ur_val(self,evec_num,x):
    omega = self.evals[evec_num]
    data = self.data
    y = (omega+1j*self.kV(x))*self.get_evec_r_val(evec_num,x)\
      -1j*self.Vp(x)*self.get_evec_theta_val(evec_num,x)
    return y
  def get_evec_utheta_val(self,evec_num,x):
    omega = self.evals[evec_num]
    data = self.data
    return (omega+1j*self.kV(x))*self.get_evec_theta_val(evec_num,x)\
      +1j*self.Vp(x)*self.get_evec_r_val(evec_num,x)
  def get_evec_uz_val(self,evec_num,x):
    omega = self.evals[evec_num]
    data = self.data
    return (omega+1j*self.kV(x))*self.get_evec_z_val(evec_num,x)
  def get_evec_divxi_val(self,evec_num,x):
    return self.get_evec_val(1,evec_num,x,prime=True)\
      +self.data['k']*self.Bt(x)*x/self.Bz(x)*self.get_evec_val(2,evec_num,x,prime=True)
  def get_evec1(self,coord,evec_num,x):
    y = Numeric.zeros([len(x)],typecode=Numeric.Complex)
    yl = '{/Symbol x}_%s'%coord
    if coord in ['1','2','3']:
      y = self.get_evec_val(int(coord),evec_num,x)
    elif coord=='r':
      y = self.get_evec_r_val(evec_num,x)
    elif coord=='t':
      y = self.get_evec_theta_val(evec_num,x)
    elif coord=='z':
      y = self.get_evec_z_val(evec_num,x)
    elif coord=='4':
      y = self.get_evec_ur_val(evec_num,x)
    elif coord=='5':
      y = self.get_evec_utheta_val(evec_num,x)
    elif coord=='6':
      y = self.get_evec_uz_val(evec_num,x)
    elif coord=='d':
      y = self.get_evec_divxi_val(evec_num,x)
      yl = 'r Div({/Symbol x})'
    else:
      print 'The second argument (coord) must be 1, 2, 3, 4, 5, 6, d, r, t, or z'
      sys.exit() 
    return (Numeric.array(y,typecode=Numeric.Complex),yl)

  def get_evec_complex(self,coord,evec_num,x,minr=None,maxr=None):
    self.imag = False
    try:
      del self.sevecs
    except:
      pass
    (y,yl)  = self.get_evec1(coord,evec_num,x)
    eimag = ((self.evals[evec_num])).imag
    if self.data['vcyl'] and eimag!=0:
      if self.data['tw'] <= 0  and abs(eimag)>0:
        evec_num = int(round(evec_num+eimag/abs(eimag)))
      else:
        self.imag = True
        del self.sevecs
      (y1,yl) = self.get_evec1(coord,evec_num,x)
      if minr == None:
        minr = -1e7
      if maxr == None:
        maxr = 1e7
      acoef = self.get_EVrot(minr=minr,maxr=maxr)
      y2 = (acoef*y+1j*acoef*y1).real
      y3 = (acoef*y+1j*acoef*y1).imag
      norm = self.get_EVnorm()
      y2 = y2*norm
      y3 = y3*norm
      return ((y2+1j*y3),yl)
    else:
      return (y,yl)
  def print_evals(self):
    for i in range(len(self.evals)):
      print i,self.evals[i]
  def print_unstable_evals(self):
    for i in range(len(self.evals)):
      if (self.evals[i].imag) > 0:
        print i,self.evals[i]
  def get_unstable_evals(self):
    unst = []
    inds = []
    for i in range(len(self.evals)):      
      if abs(self.evals[i].imag) > 0:
        inds.append(i)
        unst.append(self.evals[i])  
    sort_inds = list(Numeric.argsort(Numeric.array(unst).imag))
    sort_inds.reverse()
    self.unst_inds = Numeric.take(inds,sort_inds)
    self.unst = Numeric.take(unst,sort_inds)
    return Numeric.array(self.unst)
  def get_unst_inds(self):
    try:
      return self.unst_inds
    except:
      temp = self.get_unstable_evals()
      return self.unst_inds
  def get_numeric_evals(self):
    inf = float('inf')
    nan = inf/inf
    numevals = []
    numinds = []
    for i in range(len(self.evals)):
      if abs(self.evals[i])==inf or abs(self.evals[i])==1e20:
        continue
      else:
        numevals.append(self.evals[i])
        numinds.append(i)
    self.numinds = numinds
    return Numeric.array(numevals)
  def plot_evals(self):
    g = Gnuplot.Gnuplot(persist=1)
    g.title('Spectrum for run %d'%(self.data['num']))
    g.xlabel('Re({/Symbol w})')
    g.ylabel('Im({/Symbol w})')
    unst = self.get_numeric_evals()
    g('set xrange [%g:%g]'%(min(unst.real)*.95,max(unst.real)*1.05))
    g.plot(Gnuplot.Data(self.evals.real,self.evals.imag))
    return g
  def get_w0(self):
    self.w0 = abs(self.kV([self.data['rs']]))
    return self.w0[0]
  def Bt(self,r):
    equ = self.data['equilibr']
    if (equ==3 or equ==10 or equ==12):
      return self.data['Bt0']*r/self.data['a']
    elif equ==13:
      Bt0 = self.data['Bt0']
      Bz0 = self.data['Bz0']
      a = self.data['a']
      y = Numeric.zeros([len(r)],typecode=Numeric.Float)
      for i in range(len(y)):
        if r[i]==0:
          y[i] = 0
        else:
          y[i] = Bt0*(a-(a**2-r[i]**2)**2*math.sqrt(1-r[i]**2/a**2)/a**3)/r[i]
      return y
    else:
      return 0.
  def Bz(self,r):
    equ = self.data['equilibr']
    if (equ==3):
      return self.data['Bz0']
    elif (equ==10 or equ==12):
      Bz0 = self.data['Bz0']
      a = self.data['a']
      p0 = self.data['P0']
      Bt0 = self.data['Bt0']
      y = Numeric.zeros([len(r)],typecode=Numeric.Float)
      for i in range(len(y)):
        y[i] = math.sqrt(a**2*(Bz0**2-2*p0+2*Bt0**2)+2*r[i]**2*(p0-Bt0**2))/a
      return y
    elif equ==13:
      Bt0 = self.data['Bt0']
      Bz0 = self.data['Bz0']
      eps = self.data['eps']
      a = self.data['a']
      y = Numeric.zeros([len(r)],typecode=Numeric.Float)
      for i in range(len(y)):
        y[i] = (1./6.)*math.sqrt(-15.*Bt0**2*eps/a**8\
                         *(-25.*a**8+48.*r[i]**2*a**6-36.*r[i]**4*a**4+16.*r[i]**6*a**2-3.*r[i]**8\
                            +8.*a**8*((1.-r[i]**2/a**2))**(3/2.)+24*a**8*math.sqrt((1.-r[i]**2/a**2))\
                            -12*a**8*(math.log(2*a**2+2*math.sqrt((1.-r[i]**2/a**2))*a**2-r[i]**2)-math.log(a**2)))+36*Bz0**2)
      return y
    else:
      return 0.
  def p(self,r):
    equ = self.data['equilibr']
    data = self.data
    if (equ==3):
      return self.data['Bt0']**2*(1-r**2/self.data['a']**2)
    elif (equ==10):
      return self.data['P0']*(1-r**2/self.data['a']**2)
    elif (equ==12):
      return (2*data['P0']-data['Vp0']**2*data['a']**2)/2.*(1-r**2/data['a']**2)
    elif (equ==13):
      Bt0 = self.data['Bt0']
      Bz0 = self.data['Bz0']
      eps = self.data['eps']
      a = self.data['a']
      y = Numeric.zeros([len(r)],typecode=Numeric.Float)
      for i in range(len(y)):
        y[i] = -0.5e1/0.24e2*(1-eps)*Bt0**2/a**8\
              *(-25*a**8+48*r[i]**2*a**6-36*r[i]**4*a**4+16*r[i]**6*a**2-3*r[i]**8\
                +32*a**8*math.sqrt(1-r[i]**2/a**2)-8*a**6*math.sqrt(1-r[i]**2/a**2)*r[i]**2\
                -12*a**8*math.log(2+2*math.sqrt(1/a**2*(a**2-r[i]**2))-r[i]**2/a**2))
      return y
    else:
      return 0.
  def get_q0(self):
    equ = self.data['equilibr']
    data = self.data
    if (equ==3):
      return data['Bz0']*data['k']/data['Bt0']
    elif (equ==10):
      return data['k']*(data['Bz0']**2+2*(data['Bt0']**2-data['P0']))**0.5/data['Bt0']
    elif (equ==13):
      return 2./5./6.*data['k']*data['a']/data['Bt0']*math.sqrt(180*data['Bt0']**2*(1-data['eps'])+36*data['Bz0']**2)
  def Vp(self,r):
    data = self.data
    y = Numeric.zeros([len(r)],typecode=Numeric.Float)
    equ = data['equilibr']
    if equ==11:
      for i in range(len(y)):
        y[i] = data['Vp0']+data['epsVp']*r[i]+data['eps']*r[i]**2
    elif equ==12:
      return data['Vp0']
    return y
  def Vz(self,r):
    data = self.data
    equ = data['equilibr']
    y = Numeric.zeros([len(r)],typecode=Numeric.Float)
    p0=data['P0']
    Bt0=data['Bt0']
    if (equ<11 or equ==13):
      for i in range(len(r)):
        y[i] = data['Vz0']*(1-data['epsVz']*r[i]**2/data['a']**2)
    elif equ==12:
      for i in range(len(r)):
        y[i] = data['Vp0']/Bt0*(data['a']**2*(data['Bz0']**2-2.*p0+2.*Bt0**2)+\
                                  2.*r[i]**2*(p0-Bt0**2))**.5
    else:
      y = 0.
    return y
  def kV(self,r):
    data = self.data
    return data['m']*1j*self.Vp(r)+data['k']*1j*self.Vz(r)
  def get_wA(self):
    data = self.data
    temp = data['k']*self.Bz([data['a']])+data['m']*self.Bt(data['a'])/data['a']
    return abs(temp[0])
def get_all_list(minn=1,maxn=2200):
  d = []
  for i in range(minn,maxn):
    try:
      d.append(get_EV(i))
    except:
      continue
  return d
def get_matching(ref,d=None,nomatch_keys=None):
  if d==None:
    d = get_all_list()
  d1=get_EV(ref)
  dk = dict()
  if nomatch_keys==None:
    nomatch_keys = get_numeric_nomatch_keys()
  for key in d1.data.keys():
    if key not in nomatch_keys:
      dk[key]=d1.data[key]
  matching_files = []
  for d2 in d:
    match = True
    for key in dk.keys():
      if key=='epsVz':
        if d2.data['Vz0']==0:
          continue
      if d2.data[key]!=dk[key]:
        match=False
        break
    if match:
      matching_files.append(d2)
  return matching_files
def get_numeric_nomatch_keys():
  return ['t_assemb','num','t_solve','epsilo','NN','vcyl','alpha','rs','nu']
def get_h0_extrap(ref,d,hexp=1):
  nomatch=get_numeric_nomatch_keys()
  nomatch.append('ngrid')
  d1 = get_matching(ref,d,nomatch)
  h_min = []
  max_eval=[]
  for d2 in d1:
    try:
      h_min.append(min(d2.grid[1:]-d2.grid[0:-1])**hexp)
    except:
      h_min.append(1./(d2.data['ngrid']-1)**hexp)
    max_eval.append(max(d2.evals.imag))
  uh = list(sets.Set(h_min))  
  min_uh = min(uh)
  uuh = []
  for uh1 in uh:
    if len(Numeric.nonzero(Numeric.absolute(Numeric.array(uuh)-uh1)<1e-8))==0:
      uuh.append(uh1)
  min_eval=[]
  #print Numeric.sort(uuh)
  for hm in uuh:
    inds = Numeric.nonzero(Numeric.absolute(Numeric.array(h_min)-hm)<1e-8)
    min_eval.append(min(Numeric.take(max_eval,inds)))
  inds = Numeric.argsort(uuh)
  h_min = list(Numeric.take(uuh,inds))
  min_eval  = list(Numeric.take(min_eval,inds))
  if len(min_eval)>1:
    m=(min_eval[1]-min_eval[0])/(h_min[1]-h_min[0])
  else:
    m=0
  return min_eval[0]-m*h_min[0]
def get_Vmin_gamma(d):
  b_vals = []
  Vz0 = []
  for d1 in d:
    b_vals.append(d1.data['b'])
    Vz0.append(d1.data['Vz0'])
  ub = list(sets.Set(b_vals))
  minV=[]
  V0 = []
  for b in ub:
    inds = Numeric.nonzero(Numeric.array(b_vals)==b)
    Vb = Numeric.take(Vz0,inds)
    minV.append(min(Numeric.absolute(Vb)))
    indV0=Numeric.nonzero(Numeric.absolute(Vb)==minV[-1])
    V0.append(max(d[inds[indV0[0]]].evals.imag))
  return (ub,minV,V0)
def plot_Re_w_minG(ref,var_key,var_key2,d=None,r0=.84,minn=1,maxn=3300,g=None,neg=None,tt=['a'],rind=0,uv2=None,notitle=False):
  [g,matching,d]=plot_Re_w_varG(ref,var_key,var_key2,d=d,r0=r0,minn=minn,maxn=maxn,g=g,neg=neg,tt=tt,rind=rind,uv2=uv2,notitle=notitle,varG=0)
  g.ylabel('{/Symbol w}_{{/Symbol G}_{min}}')
  nums = [matching[i].num for i in range(len(matching))]
  print nums
  inds = Numeric.argsort(nums)
  for i in  range(len(nums)):
    if matching[inds[i]].data['Vz0']!=0:
      uniq_num = nums[inds[i]]
      break
  if type(r0)==type([]):
    fnb = 'plots/%d_var_%s_%s_Re_minG_r0_multi'%(uniq_num,var_key,var_key2)
  else:
    fnb = 'plots/%d_var_%s_%s_Re_minG_r0_%g'%(uniq_num,var_key,var_key2,r0)
  print fnb
  g('set key out')
  g.hardcopy(fnb+'.eps',eps=True)
  g.hardcopy(fnb+'.ps',size=(7,5),fontsize=24,eps=False)
  g.refresh()
  return [g,matching,d]
def plot_Re_w_maxG(ref,var_key,var_key2,d=None,r0=.84,minn=1,maxn=3300,g=None,neg=None,tt=['a'],rind=0,uv2=None,notitle=False):
  [g,matching,d]=plot_Re_w_varG(ref,var_key,var_key2,d=d,r0=r0,minn=minn,maxn=maxn,g=g,neg=neg,tt=tt,rind=rind,uv2=uv2,notitle=notitle,varG=-1)
  g.ylabel('{/Symbol w}_{{/Symbol G}_{max}}')
  nums = [matching[i].num for i in range(len(matching))]
  print nums
  inds = Numeric.argsort(nums)
  for i in  range(len(nums)):
    if matching[inds[i]].data['Vz0']!=0:
      uniq_num = nums[inds[i]]
      break
  if type(r0)==type([]):
    fnb = 'plots/%d_var_%s_%s_Re_maxG_r0_multi'%(uniq_num,var_key,var_key2)
  else:
    fnb = 'plots/%d_var_%s_%s_Re_maxG_r0_%g'%(uniq_num,var_key,var_key2,r0)
  print fnb
  g('set key out')
  g.hardcopy(fnb+'.eps',eps=True)
  g.hardcopy(fnb+'.ps',size=(7,5),fontsize=24,eps=False)
  g.refresh()
  return [g,matching,d]
def plot_Re_w_varG(ref,var_key,var_key2,d=None,r0=.84,minn=1,maxn=3300,g=None,neg=None,tt=['a'],rind=0,uv2=None,notitle=False,varG=-1):
  if d==None:
    d = get_all_list(minn=minn,maxn=maxn)
  if type(ref)!=type(1):
    print 'ref must be of type int'
  if type(var_key)!=type(''):
    print 'var_key must be of type string'
  if type(var_key2)!=type(''):
    print 'var_key2 must be of type string'
  if type(rind)!=type(1):
    print 'rind must be of type int'
  matching = get_matching(ref,d=d,nomatch_keys=get_numeric_nomatch_keys()+[var_key,var_key2])
  if g==None:
    g = Gnuplot.Gnuplot(persist=1)
    g('set key left top Left reverse')
    g.xlabel(var_key)
    g.ylabel('{/Symbol w}_{min{/Symbol G}}')
  xfac = 1
  if var_key=='Vz0':
    wA = abs(matching[0].get_wA())
    xfac = 1./(wA/matching[0].data['k'])
    g.xlabel('V_z(r=0)/v_A')   
  g.itemlist = []
  for i in range(len(matching)):
    if matching[i].data['num']==ref:
      ref_ind = i
      break
  ttl = ''
  for tt1 in tt:
    ttl=ttl+'%s=%g, '%(tt1,matching[ref_ind].data[tt1])
  ttl = ttl[:-2]
  if notitle==False:
    g.title(ttl)#'V_z(r=a)=%g'%matching[0].Vz([matching[0].data['a']]))
  else:
    g('unset title')
  var_list = []
  var2_list = []
  min_eval = []
  nums = [] 
  matching2 = []
  for i in range(len(matching)):
    evals = matching[i].evals
    if neg!=None:
      inds = Numeric.nonzero(evals.real<1e-5)
      evals = Numeric.take(evals,inds)
    ind = Numeric.argsort(evals.imag)
    if evals[ind[varG]].imag==0:
      continue
    min_eval.append(evals[ind[varG]].real)
    var_list.append(matching[i].data[var_key]) 
    var2_list.append(matching[i].data[var_key2])
    nums.append(matching[i].num)
    matching2.append(matching[i])
  matching = matching2
  if len(var2_list)==0:
    return
  if uv2==None:
    uv2 = Numeric.sort(list(sets.Set(var2_list)))
  else:
    uv2 = Numeric.sort(uv2)
  i = 1
  for u2 in uv2:
    inds = Numeric.nonzero(Numeric.array(var2_list)==u2)
    if len(inds)<=1:
      i = i+1
      continue
    v1 = Numeric.take(var_list,inds)*xfac
    ev1 = Numeric.take(min_eval,inds)
    d1 = matching[inds[rind]]
    #print u2
    #print v1,ev1
    g.itemlist.append(Gnuplot.Data(v1,ev1,title='%s=%g'%(var_key2,d1.data[var_key2]),with='points %d'%i))
    if type(r0)==type([]):
      for k in range(len(r0)):
        r01 = r0[k]
        fit = d1.data['k']*Numeric.array([matching[j].Vz([r01])[0] for j in inds])+ev1[rind]-d1.data['k']*d1.Vz([r01])
        g.itemlist.append(Gnuplot.Data(v1,fit,title='k V_z(r=%g)'%r01,with='lines %d'%(k+2)))
      break
    else:
      fit = d1.data['k']*Numeric.array([matching[j].Vz([r0])[0] for j in inds])+ev1[rind]-d1.data['k']*d1.Vz([r0])
      g.itemlist.append(Gnuplot.Data(v1,fit,title='k V_z(r=%g)'%r0,with='lines %d'%i))
    i = i+1
  return [g,matching,d]
  if type(r0)==type([]):
    fnb = 'plots/%d_var_%s_%s_Re_minG_r0_multi'%(min(nums),var_key,var_key2)
  else:
    fnb = 'plots/%d_var_%s_%s_Re_minG_r0_%g'%(min(nums),var_key,var_key2,r0)
  print fnb
  g('set key out')
def plot_stability_boundary(ref,var_key,var_key2,d=None,minn=1,maxn=4000,\
g=None,tt=['a'],notitle=False,r0=.84,wb=None,ref_g=None):
  if d==None:
    d = get_all_list(minn=minn,maxn=maxn)
  if type(ref)!=type(1):
    print 'ref must be of type int'
  if type(var_key)!=type(''):
    print 'var_key must be of type string'
  if type(var_key2)!=type(''):
    print 'var_key2 must be of type string'
  matching = get_matching(ref,d=d,nomatch_keys=get_numeric_nomatch_keys()+[var_key,var_key2])
  if g==None:
    g = Gnuplot.Gnuplot(persist=1)
    g('set key left top Left reverse')
    g.xlabel(var_key)
    g.ylabel(var_key2)
  var_list = []
  var2_list = []
  stable = []
  if ref_g==None:
    ref_g = max(get_EV(ref).evals.imag)
  for mat in matching:
    var_list.append(mat.data[var_key])
    var2_list.append(mat.data[var_key2])
    stable.append(max(mat.evals.imag)<=ref_g)
  stable = Numeric.array(stable)
  inds1 = Numeric.nonzero(stable)
  inds2 = Numeric.nonzero(1-stable)
  g.itemlist=[]
  if len(inds1)==0:
    print 'There were no matching runs that were more stable than the reference shot.'
  else:
    x1=Numeric.take(var_list,inds1)
    y1=Numeric.take(var2_list,inds1)
    g.itemlist.append(Gnuplot.Data(x1,y1,with='points 1 7',title='{/Symbol G}_{max}<=%g'%ref_g))
  if len(inds2)==0:
    print 'All matching runs were more stable than the reference shot.'
  else:
    x2=Numeric.take(var_list,inds2)
    y2=Numeric.take(var2_list,inds2)  
    g.itemlist.append(Gnuplot.Data(x2,y2,with='points 2 1',title='{/Symbol G}_{max}> %g'%ref_g))
  for i in range(len(matching)):
    if matching[i].data['num']==ref:
      ref_ind = i
      break
  minx = min(var_list)
  maxx = max(var_list)
  g('set xrange [%g:%g]'%(minx,maxx))
  if var_key=='epsVz' and var_key2=='Vz0':
    xl = Numeric.arange(minx-abs(minx)*2,maxx+abs(maxx)*2,.01)
    g.itemlist.append(Gnuplot.Data(xl,2./abs(xl),with='lines'))
    if wb==None:
      indmin=Numeric.argmin(matching[ref_ind].evals.imag)
      wb=matching[ref_ind].evals[indmin].real
    print wb
    data = matching[ref_ind].data
    g.itemlist.append(Gnuplot.Data(xl,abs(wb)/abs(data['k']*(1-xl*r0**2)),with='lines'))
  g.refresh()
  ttl = ''
  for tt1 in tt:
    ttl=ttl+'%s=%g, '%(tt1,matching[ref_ind].data[tt1])
  ttl = ttl[:-2]
  if notitle==False:
    g.title(ttl)#'V_z(r=a)=%g'%matching[0].Vz([matching[0].data['a']]))
  else:
    g('unset title')
  fn = 'plots/%d_%s_%s_stability_boundary.eps'%(min([mat.data['num'] for mat in matching]),var_key,var_key2)
  print 'Output to %s'%fn
  g.hardcopy(fn,eps=True,solid=False)
  return [g,matching,d]
def find_match(mat_dict,d=None,minn=1,maxn=4500):
  if d==None:
    for i in range(minn,maxn):
      try:
        ref = get_EV(i)
      except:
        continue
      match=True
      for key in mat_dict.keys():
        if ref.data[key]==mat_dict[key]:
          continue
        else:
          match=False
          break
      if match:
        return i
  else:
    for d1 in d:
      match=True
      for key in mat_dict.keys():
        if d1.data[key]==mat_dict[key]:
          continue
        else:
          match=False
          break
      if match:
        return d1.data['num']
def add_to_all_list(d,minn,maxn):
  for i in range(minn,maxn):
    d.append(get_EV(i))
  return d
if __name__ == '__main__':
  d = get_all_list()
  nomatch_keys = get_numeric_nomatch_keys()
  nomatch_keys.append('Vz0')
  matching_files = get_matching(1102,d,nomatch_keys)
  print get_h0_extrap(1102,d)
  (b,Vz0,ev)= get_Vmin_gamma(d)
  g = Gnuplot.Gnuplot()
  inds = Numeric.nonzero(Numeric.array(Vz0)==0)
  g('set log y')
  g.plot(Gnuplot.Data(Numeric.take(b,inds),Numeric.take(ev,inds)))
