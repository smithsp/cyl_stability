import math
class equ:
  def __init__(self,data):
    self.data = data
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
  def k_from_q(self,q=2.55):
    data = self.data
    return q/(self.Bz([data['a']])[0]/self.Bt(data['a']))
  def qa(self):
    data = self.data
    return data['k']*self.Bz([data['a']])[0]/self.Bt(data['a'])
if __name__=='__main__':
  import get_EV, Numeric, Gnuplot, sys
  ref = int(sys.argv[1])
  d = equ(get_EV.get_EV(ref).data)
  dBt0 = d.data['Bt0']
  dBz0 = d.data['Bz0']
  d.data['P0'] = .3
  dp0 = d.data['P0']
  i=1
  for Bz0 in Numeric.arange(dBz0,dBz0+dBz0*1.001,dBz0/5.):
    d.data['Bz0'] = Bz0
    for Bt0 in Numeric.arange(dBt0-dBt0*.4,dBt0+dBt0/5.,dBt0/10.):
      d.data['Bt0'] = Bt0
      if Bt0**2-d.data['P0']<0 or Bz0**2-2*d.data['P0']+2*Bt0**2<0:
        continue
      for b in Numeric.arange(1.04,1.18,.03):
        #print i
        i = i+1
        print d.data['P0']
        print d.data['Bz0']
        print d.data['Bt0']
        print d.k_from_q(2.55)
        print b
  raise '%d'%(i-1)
