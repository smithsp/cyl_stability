#! /usr/pppl/python/2.3.5x/bin/python
import Numeric, sys, math

class bspline:

  def __init__(self,xj,p1=None,p2=None,p3=None,p4=None,Lend0=False,deriv=False):
    
    self.deriv = deriv
    self.xj = xj
    self.dx = Numeric.zeros([4],typecode=Numeric.Float)
    self.extent = Numeric.zeros([3])
    if p2!=None:
      self.dx[1] = xj-p2
      self.extent[1] = 1
    if p1!=None:
      if p2!=None:
        self.extent[1]=2
        self.dx[0]=p2-p1
    if p3!=None:
      self.extent[2]=1
      self.dx[2] = p3-xj
    if p4!=None:
      if p3!=None:
        self.extent[2]=2
        self.dx[3] = p4-p3
    if min(self.dx)<0:
      raise ValueError
    dx = self.dx
    A = Numeric.zeros([4],typecode=Numeric.Float)
    B = Numeric.zeros([4],typecode=Numeric.Float)
    C = Numeric.zeros([4],typecode=Numeric.Float)
    D = Numeric.zeros([4],typecode=Numeric.Float)
    if self.extent[1]==2 and self.extent[2]==2:
      num1 = 0.2e1*(dx[2]+dx[1])*(dx[2]+dx[3]+dx[1])/dx[0]/(dx[0]+dx[1])/(dx[0]*pow(dx[2],0.2e1)+dx[0]*dx[2]*dx[3]+0.2e1*dx[0]*dx[1]*dx[2]+dx[0]*dx[1]*dx[3]+pow(dx[1],0.2e1)*dx[3]+0.2e1*pow(dx[1],0.2e1)*dx[2]+0.2e1*dx[1]*pow(dx[2],0.2e1)+0.2e1*dx[1]*dx[2]*dx[3])
      A[0] = -pow(xj-dx[0]-dx[1],0.3e1)/0.3e1*num1
      A[1] = pow(xj-dx[0]-dx[1],0.2e1)*num1
      A[2] = (-xj+dx[0]+dx[1])*num1
      A[3] = (0.1e1/0.3e1)*num1
      
      num2 = 0.2e1/(dx[0]+dx[1])/dx[1]/(dx[0]*pow(dx[2],0.2e1)+dx[0]*dx[2]*dx[3]+0.2e1*dx[0]*dx[1]*dx[2]+dx[0]*dx[1]*dx[3]+pow(dx[1],0.2e1)*dx[3]+0.2e1*pow(dx[1],0.2e1)*dx[2]+0.2e1*dx[1]*pow(dx[2],0.2e1)+0.2e1*dx[1]*dx[2]*dx[3])
      B[0] = (dx[3]/0.3e1+xj+0.2e1/0.3e1*dx[2])*pow(dx[1],0.4e1)+(0.2e1/0.3e1*pow(dx[2],0.2e1)+(0.2e1/0.3e1*dx[3]+0.4e1/0.3e1*dx[0])*dx[2]+0.2e1/0.3e1*dx[3]*dx[0]-0.2e1*xj*xj+0.2e1*xj*dx[0])*pow(dx[1],0.3e1)+((-xj+dx[0])*pow(dx[2],0.2e1)+(0.2e1/0.3e1*pow(dx[0],0.2e1)-xj*dx[3]+dx[3]*dx[0]-0.2e1*xj*xj)*dx[2]+pow(dx[0],0.2e1)*dx[3]/0.3e1+pow(xj,0.3e1)+xj*pow(dx[0],0.2e1)-xj*xj*dx[3]-0.3e1*xj*xj*dx[0])*pow(dx[1],0.2e1)+((-xj*dx[0]+pow(dx[0],0.2e1)/0.3e1)*pow(dx[2],0.2e1)+(pow(dx[0],0.2e1)*dx[3]/0.3e1-0.2e1*xj*xj*dx[0]+0.4e1/0.3e1*pow(xj,0.3e1)-xj*dx[0]*dx[3])*dx[2]-xj*xj*dx[0]*dx[3]+0.2e1/0.3e1*pow(xj,0.3e1)*dx[3]-xj*xj*pow(dx[0],0.2e1)+pow(xj,0.3e1)*dx[0])*dx[1]+pow(xj,0.3e1)*pow(dx[2],0.2e1)/0.3e1+(pow(xj,0.3e1)*dx[3]/0.3e1+0.2e1/0.3e1*pow(xj,0.3e1)*dx[0])*dx[2]+pow(xj,0.3e1)*pow(dx[0],0.2e1)/0.3e1+pow(xj,0.3e1)*dx[0]*dx[3]/0.3e1
      B[1] = -pow(dx[1],0.4e1)-0.2e1*pow(dx[1],0.3e1)*dx[0]+0.4e1*xj*pow(dx[1],0.3e1)-pow(dx[0],0.2e1)*pow(dx[1],0.2e1)-0.3e1*xj*xj*pow(dx[1],0.2e1)+0.2e1*xj*pow(dx[1],0.2e1)*dx[3]+0.4e1*xj*dx[2]*pow(dx[1],0.2e1)+0.6e1*xj*dx[0]*pow(dx[1],0.2e1)+pow(dx[1],0.2e1)*dx[2]*dx[3]+pow(dx[1],0.2e1)*pow(dx[2],0.2e1)+dx[1]*dx[0]*pow(dx[2],0.2e1)+dx[1]*dx[0]*dx[2]*dx[3]-0.3e1*xj*xj*dx[0]*dx[1]+0.2e1*dx[1]*xj*pow(dx[0],0.2e1)-0.4e1*xj*xj*dx[2]*dx[1]+0.4e1*xj*dx[0]*dx[2]*dx[1]+0.2e1*xj*dx[0]*dx[1]*dx[3]-0.2e1*xj*xj*dx[1]*dx[3]-xj*xj*pow(dx[2],0.2e1)-0.2e1*xj*xj*dx[0]*dx[2]-xj*xj*dx[2]*dx[3]-xj*xj*dx[0]*dx[3]-xj*xj*pow(dx[0],0.2e1)
      B[2] = -0.2e1*pow(dx[1],0.3e1)-0.3e1*dx[0]*pow(dx[1],0.2e1)+0.3e1*pow(dx[1],0.2e1)*xj-0.2e1*pow(dx[1],0.2e1)*dx[2]-pow(dx[1],0.2e1)*dx[3]+0.3e1*xj*dx[0]*dx[1]+0.4e1*dx[1]*xj*dx[2]+0.2e1*dx[1]*xj*dx[3]-0.2e1*dx[0]*dx[1]*dx[2]-dx[0]*dx[1]*dx[3]-dx[1]*pow(dx[0],0.2e1)+0.2e1*dx[0]*xj*dx[2]+xj*pow(dx[0],0.2e1)+xj*dx[2]*dx[3]+xj*dx[0]*dx[3]+xj*pow(dx[2],0.2e1)
      B[3] = -pow(dx[1],0.2e1)-0.2e1/0.3e1*dx[1]*dx[3]-0.4e1/0.3e1*dx[1]*dx[2]-dx[0]*dx[1]-pow(dx[0],0.2e1)/0.3e1-pow(dx[2],0.2e1)/0.3e1-0.2e1/0.3e1*dx[0]*dx[2]-dx[2]*dx[3]/0.3e1-dx[3]*dx[0]/0.3e1
      B = B*num2
      
      num3 = 0.2e1/(dx[2]+dx[3])/dx[2]/(dx[0]*pow(dx[2],0.2e1)+dx[0]*dx[2]*dx[3]+0.2e1*dx[0]*dx[1]*dx[2]+dx[0]*dx[1]*dx[3]+pow(dx[1],0.2e1)*dx[3]+0.2e1*pow(dx[1],0.2e1)*dx[2]+0.2e1*dx[1]*pow(dx[2],0.2e1)+0.2e1*dx[1]*dx[2]*dx[3])
      C[0] = -xj*xj*dx[0]*dx[2]*dx[3]+xj*pow(dx[1],0.2e1)*dx[2]*dx[3]-0.2e1/0.3e1*pow(xj,0.3e1)*dx[1]*dx[3]-0.3e1*xj*xj*pow(dx[2],0.2e1)*dx[3]+dx[1]*xj*dx[0]*pow(dx[2],0.2e1)-pow(xj,0.3e1)*pow(dx[2],0.2e1)-0.2e1*xj*xj*pow(dx[2],0.2e1)*dx[1]+0.2e1/0.3e1*pow(dx[2],0.3e1)*dx[0]*dx[3]+pow(dx[3],0.2e1)*pow(dx[2],0.2e1)*dx[0]/0.3e1+pow(dx[1],0.2e1)*pow(dx[3],0.2e1)*dx[2]/ 0.3e1 + 0.2e1 / 0.3e1 * dx[0] * dx[1] * pow(dx[2], 0.3e1) + 0.4e1 / 0.3e1 * dx[1] * pow(dx[2], 0.3e1) * dx[3] + 0.2e1 / 0.3e1 * dx[1] * pow(dx[3], 0.2e1) * pow(dx[2], 0.2e1) + 0.2e1 / 0.3e1 * dx[1] * pow(dx[2], 0.4e1) + dx[0] * dx[1] * pow(dx[2], 0.2e1) * dx[3] - 0.2e1 * xj * pow(dx[2], 0.3e1) * dx[3] - xj * xj * pow(dx[3], 0.2e1) * dx[2] - xj * pow(dx[2], 0.2e1) * pow(dx[3], 0.2e1) - pow(xj, 0.3e1) * dx[0] * dx[1] / 0.3e1 + pow(dx[1], 0.2e1) * pow(dx[2], 0.2e1) * dx[3] + dx[0] * dx[1] * pow(dx[3], 0.2e1) * dx[2] / 0.3e1 - 0.2e1 / 0.3e1 * pow(xj, 0.3e1) * dx[0] * dx[2] - 0.2e1 * xj * xj * dx[2] * dx[1] * dx[3] + dx[1] * xj * dx[0] * dx[2] * dx[3] + 0.2e1 / 0.3e1 * pow(dx[1], 0.2e1) * pow(dx[2], 0.3e1) - xj * xj * dx[0] * pow(dx[2], 0.2e1) + pow(dx[1], 0.2e1) * xj * pow(dx[2], 0.2e1) - pow(xj, 0.3e1) * dx[2] * dx[3] - xj * pow(dx[2], 0.4e1) - pow(xj, 0.3e1) * pow(dx[3], 0.2e1) / 0.3e1 + pow(dx[2], 0.4e1) * dx[0] / 0.3e1 - pow(xj, 0.3e1) * dx[0] * dx[3] / 0.3e1 - 0.4e1 / 0.3e1 * pow(xj, 0.3e1) * dx[1] * dx[2] - pow(xj, 0.3e1) * pow(dx[1], 0.2e1) / 0.3e1 - 0.2e1 * xj * xj * pow(dx[2], 0.3e1)
      C[1] = -dx[1]*dx[0]*dx[2]*dx[3]+0.2e1*xj*dx[0]*pow(dx[2],0.2e1)+0.2e1*xj*xj*dx[0]*dx[2]+xj*xj*dx[0]*dx[3]+0.2e1*xj*dx[0]*dx[2]*dx[3]-dx[1]*dx[0]*pow(dx[2],0.2e1)+xj*xj*dx[0]*dx[1]+0.3e1*xj*xj*pow(dx[2],0.2e1)+0.2e1*xj*xj*dx[1]*dx[3]+0.3e1*xj*xj*dx[2]*dx[3]+xj*xj*pow(dx[1],0.2e1)+0.4e1*xj*xj*dx[2]*dx[1]+0.4e1*xj*pow(dx[2],0.3e1)+0.2e1*xj*dx[2] * pow(dx[3], 0.2e1) + 0.6e1 * xj * pow(dx[2], 0.2e1) * dx[3] + xj * xj * pow(dx[3], 0.2e1) + 0.2e1 * pow(dx[2], 0.3e1) * dx[3] + pow(dx[3], 0.2e1) * pow(dx[2], 0.2e1) - pow(dx[1], 0.2e1) * dx[2] * dx[3] + pow(dx[2], 0.4e1) + 0.4e1 * xj * dx[2] * dx[1] * dx[3] + 0.4e1 * xj * pow(dx[2], 0.2e1) * dx[1] - pow(dx[1], 0.2e1) * pow(dx[2], 0.2e1)
      C[2] = -0.2e1*dx[0]*xj*dx[2]-xj*dx[0]*dx[1]-xj*dx[0]*dx[3]-dx[0]*pow(dx[2],0.2e1)-dx[0]*dx[2]*dx[3]-0.3e1*pow(dx[2],0.2e1)*dx[3]-dx[2]*pow(dx[3],0.2e1)-0.2e1*dx[1]*xj*dx[3]-0.2e1*pow(dx[2],0.3e1)-0.2e1*dx[1]*pow(dx[2],0.2e1)-xj*pow(dx[3],0.2e1)-0.2e1*dx[1]*dx[2]*dx[3]-0.4e1*dx[1]*xj*dx[2]-pow(dx[1],0.2e1)*xj-0.3e1*xj*dx[2]*dx[3]-0.3e1*xj*pow(dx[2],0.2e1)
      C[3] = pow(dx[1],0.2e1)/0.3e1+0.2e1/0.3e1*dx[1]*dx[3]+dx[0]*dx[1]/0.3e1+0.4e1/0.3e1*dx[1]*dx[2]+dx[3]*dx[0]/0.3e1+0.2e1/0.3e1*dx[0]*dx[2]+pow(dx[3],0.2e1)/0.3e1+pow(dx[2],0.2e1)+dx[2]*dx[3]
      C = C*num3
      
      num4 = 0.2e1*(dx[2]+dx[1])*(dx[2]+dx[0]+dx[1])/dx[3]/(dx[2]+dx[3])/(dx[0]*pow(dx[2],0.2e1)+dx[0]*dx[2]*dx[3]+0.2e1*dx[0]*dx[1]*dx[2]+dx[0]*dx[1]*dx[3]+pow(dx[1],0.2e1)*dx[3]+0.2e1*pow(dx[1],0.2e1)*dx[2]+0.2e1*dx[1]*pow(dx[2],0.2e1)+0.2e1*dx[1]*dx[2]*dx[3])
      D[0] = pow(xj+dx[2]+dx[3],0.3e1)/0.3e1
      D[1] = -pow(xj+dx[2]+dx[3],0.2e1)
      D[2] = xj+dx[2]+dx[3]
      D[3] = -0.1e1/0.3e1
      D = D*num4
    elif self.extent[1]==1 and self.extent[2]==2:
      num2 = 0.2e1 * pow(dx[1], -0.3e1) / (0.2e1 * dx[1] * dx[2] + dx[1] * dx[3] + 0.2e1 * pow(dx[2], 0.2e1) + 0.2e1 * dx[2] * dx[3])
      B[0] = pow(xj - dx[1], 0.2e1) * (0.3e1 * xj * pow(dx[1], 0.2e1) + 0.4e1 * dx[1] * xj * dx[2] + 0.2e1 * dx[1] * xj * dx[3] + xj * dx[2] * dx[3] + xj * pow(dx[2], 0.2e1) + 0.2e1 * dx[2] * pow(dx[1], 0.2e1) + pow(dx[1], 0.2e1) * dx[3] + 0.2e1 * dx[1] * dx[2] * dx[3] + 0.2e1 * dx[1] * pow(dx[2], 0.2e1)) / 0.3e1
      B[1] = -(xj - dx[1]) * (0.3e1 * xj * pow(dx[1], 0.2e1) + 0.4e1 * dx[1] * xj * dx[2] + 0.2e1 * dx[1] * xj * dx[3] + xj * dx[2] * dx[3] + xj * pow(dx[2], 0.2e1) - pow(dx[1], 0.3e1) + dx[1] * dx[2] * dx[3] + dx[1] * pow(dx[2], 0.2e1))
      B[2] = -0.2e1 * pow(dx[1], 0.3e1) - 0.2e1 * dx[2] * pow(dx[1], 0.2e1) - pow(dx[1], 0.2e1) * dx[3] + 0.3e1 * xj * pow(dx[1], 0.2e1) + 0.2e1 * dx[1] * xj * dx[3] + 0.4e1 * dx[1] * xj * dx[2] + xj * pow(dx[2], 0.2e1) + xj * dx[2] * dx[3]
      B[3] = -pow(dx[1], 0.2e1) - 0.4e1 / 0.3e1 * dx[1] * dx[2] - 0.2e1 / 0.3e1 * dx[1] * dx[3] - dx[2] * dx[3] / 0.3e1 - pow(dx[2], 0.2e1) / 0.3e1
      B = B*num2
      
      num3 = 2 / dx[2] / (dx[2] + dx[3]) / (2 * dx[1] * dx[2] + dx[1] * dx[3] + 2 *  pow(dx[2], 2) + 2 * dx[2] * dx[3]) / dx[1]
      C[0] = pow(dx[1], 0.2e1) * pow(dx[2], 0.2e1) * dx[3] + xj * pow(dx[1], 0.2e1) * pow(dx[2], 0.2e1) + pow(dx[1], 0.2e1) * pow(dx[3], 0.2e1) * dx[2] / 0.3e1 + xj * pow(dx[1], 0.2e1) * dx[2] * dx[3] - pow(xj, 0.3e1) * pow(dx[1], 0.2e1) / 0.3e1 + 0.2e1 / 0.3e1 * pow(dx[1], 0.2e1) * pow(dx[2], 0.3e1) - 0.2e1 / 0.3e1 * pow(xj, 0.3e1) * dx[1] * dx[3] - 0.2e1 * xj * xj * dx[1] * dx[2] * dx[3] - 0.2e1 * xj * xj * dx[1] * pow(dx[2], 0.2e1) + 0.2e1 / 0.3e1 * pow(dx[3], 0.2e1) * dx[1] * pow(dx[2], 0.2e1) + 0.4e1 / 0.3e1 * pow(dx[2], 0.3e1) * dx[3] * dx[1] - 0.4e1 / 0.3e1 * pow(xj, 0.3e1) * dx[1] * dx[2] + 0.2e1 / 0.3e1 * pow(dx[2], 0.4e1) * dx[1] - pow(xj, 0.3e1) * pow(dx[3], 0.2e1) / 0.3e1 - 0.3e1 * xj * xj * pow(dx[2], 0.2e1) * dx[3] - xj * xj * pow(dx[3], 0.2e1) * dx[2] - 0.2e1 * xj * xj * pow(dx[2], 0.3e1) - pow(xj, 0.3e1) * pow(dx[2], 0.2e1) - pow(xj, 0.3e1) * dx[2] * dx[3] - xj * pow(dx[2], 0.4e1) - 0.2e1 * xj * pow(dx[2], 0.3e1) * dx[3] - xj * pow(dx[2], 0.2e1) * pow(dx[3], 0.2e1)
      C[1] = xj * xj * pow(dx[1], 0.2e1) - pow(dx[1], 0.2e1) * pow(dx[2], 0.2e1) - pow(dx[1], 0.2e1) * dx[2] * dx[3] + 0.4e1 * xj * dx[1] * pow(dx[2], 0.2e1) + 0.4e1 * dx[1] * xj * xj * dx[2] + 0.2e1 * dx[1] * xj * xj * dx[3] + 0.4e1 * xj * dx[1] * dx[2] * dx[3] + 0.3e1 * xj * xj * pow(dx[2], 0.2e1) + xj * xj * pow(dx[3], 0.2e1) + 0.2e1 * xj * dx[2] * pow(dx[3], 0.2e1) + 0.3e1 * xj * xj * dx[2] * dx[3] + 0.4e1 * xj * pow(dx[2], 0.3e1) + 0.6e1 * xj * pow(dx[2], 0.2e1) * dx[3] + pow(dx[2], 0.2e1) * pow(dx[3], 0.2e1) + pow(dx[2], 0.4e1) + 0.2e1 * pow(dx[2], 0.3e1) * dx[3]
      C[2] = -xj * pow(dx[1], 0.2e1) - 0.4e1 * dx[1] * xj * dx[2] - 0.2e1 * dx[1] * dx[2] * dx[3] - 0.2e1 * dx[1] * pow(dx[2], 0.2e1) - 0.2e1 * dx[1] * xj * dx[3] - 0.3e1 * xj * pow(dx[2], 0.2e1) - 0.3e1 * pow(dx[2], 0.2e1) * dx[3] - 0.2e1 * pow(dx[2], 0.3e1) - 0.3e1 * xj * dx[2] * dx[3] - xj * pow(dx[3], 0.2e1) - dx[2] * pow(dx[3], 0.2e1)
      C[3] = pow(dx[1], 0.2e1) / 0.3e1 + 0.4e1 / 0.3e1 * dx[1] * dx[2] + 0.2e1 / 0.3e1 * dx[1] * dx[3] + pow(dx[2], 0.2e1) + dx[2] * dx[3] + pow(dx[3], 0.2e1) / 0.3e1
      C = C*num3
      
      num4 = 0.2e1 * pow(dx[2] + dx[1], 0.2e1) / dx[3] / dx[1] / (dx[2] + dx[3]) / (0.2e1 * dx[1] * dx[2] + dx[1] * dx[3] + 0.2e1 * pow(dx[2], 0.2e1) + 0.2e1 * dx[2] * dx[3])
      D[0] = pow(xj + dx[2] + dx[3], 0.3e1) / 0.3e1
      D[1] = -pow(xj + dx[2] + dx[3], 0.2e1)
      D[2] = xj + dx[2] + dx[3]
      D[3] = -0.1e1 / 0.3e1
      D = D*num4
      
    elif self.extent[1]==0 and self.extent[2]==2:
      num3 = 0.2e1*pow(dx[2],-0.3e1)*pow(dx[3],-0.2e1)
      C[0] = -(0.3e1*pow(dx[2],0.4e1)+0.6e1*pow(dx[2],0.3e1)*dx[3]+0.6e1*xj*pow(dx[2],0.3e1)+0.9e1*xj*pow(dx[2],0.2e1)*dx[3]+0.3e1*pow(dx[2],0.2e1)*pow(dx[3],0.2e1)+0.3e1*xj*xj*pow(dx[2],0.2e1)+0.3e1*xj*dx[2]*pow(dx[3],0.2e1)+0.3e1*xj*xj*dx[2]*dx[3]+xj*xj*pow(dx[3],0.2e1))*xj/0.3e1
      C[1] = pow(dx[2],0.4e1)+0.2e1*pow(dx[2],0.3e1)*dx[3]+0.4e1*xj*pow(dx[2],0.3e1)+0.6e1*xj*pow(dx[2],0.2e1)*dx[3]+pow(dx[2],0.2e1)*pow(dx[3],0.2e1)+0.3e1*xj*xj*pow(dx[2],0.2e1)+0.2e1*xj*dx[2]*pow(dx[3],0.2e1)+0.3e1*xj*xj*dx[2]*dx[3]+xj*xj*pow(dx[3],0.2e1)
      C[2] = -0.2e1*pow(dx[2],0.3e1)-0.3e1*xj*pow(dx[2],0.2e1)-0.3e1*dx[3]*pow(dx[2],0.2e1)-0.3e1*xj*dx[2]*dx[3]-pow(dx[3],0.2e1)*dx[2]-xj*pow(dx[3],0.2e1)
      C[3] = pow(dx[2],0.2e1)+dx[2]*dx[3]+pow(dx[3],0.2e1)/0.3e1
      C = C*num3
      
      num4 = 0.2e1*pow(dx[3],-0.3e1)
      D[0] = pow(xj+dx[2]+dx[3],0.3e1)/0.3e1
      D[1] = -pow(xj+dx[2]+dx[3],0.2e1)
      D[2] = xj+dx[2]+dx[3]
      D[3] = -0.1e1/0.3e1
      D = D*num4
    elif self.extent[1]==0 and self.extent[2]==1:

      num3 = 0.2e1 * pow(dx[2], -0.3e1)
      C[0] = pow(xj + dx[2], 0.3e1) / 0.3e1
      C[1] = -pow(xj + dx[2], 0.2e1)
      C[2] = xj + dx[2]
      C[3] = -0.1e1 / 0.3e1
      
      C = C*num3
    elif self.extent[1]==2 and self.extent[2]==1:
      num1 = 0.2e1*pow(dx[2]+dx[1],0.2e1)/dx[0]/(dx[1]+dx[0])/(dx[0]*dx[2]+0.2e1*dx[1]*dx[0]+0.2e1*pow(dx[1],0.2e1)+0.2e1*dx[1]*dx[2])/dx[2]
      A[0] = -pow(xj-dx[1]-dx[0],0.3e1)/0.3e1
      A[1] = pow(xj-dx[1]-dx[0],0.2e1)
      A[2] = dx[1]+dx[0]-xj
      A[3] = 0.1e1/0.3e1
      A = A*num1
      
      num2 = 2/(dx[1]+dx[0])/dx[2]/dx[1]/(dx[0]*dx[2]+2*dx[1]*dx[0]+2*pow(dx[1],2)+2*dx[1]*dx[2])
      B[0] = 0.2e1/0.3e1*pow(dx[1],0.4e1)*dx[2]+pow(dx[1],0.4e1)*xj-0.2e1*xj*xj*pow(dx[1],0.3e1)+0.2e1*xj*pow(dx[1],0.3e1)*dx[0]+0.4e1/0.3e1*pow(dx[1],0.3e1)*dx[0]*dx[2]+0.2e1/0.3e1*pow(dx[2],0.2e1)*pow(dx[1],0.3e1)+pow(dx[2],0.2e1)*dx[0]*pow(dx[1],0.2e1)-0.2e1*pow(dx[1],0.2e1)*xj*xj*dx[2]+xj*pow(dx[0],0.2e1)*pow(dx[1],0.2e1)-0.3e1*xj*xj*dx[0]*pow(dx[1],0.2e1)+0.2e1/0.3e1*pow(dx[1],0.2e1)*pow(dx[0],0.2e1)*dx[2]+pow(xj,0.3e1)*pow(dx[1],0.2e1)-pow(dx[1],0.2e1)*xj*pow(dx[2],0.2e1)+0.4e1/0.3e1*dx[1]*dx[2]*pow(xj,0.3e1)+pow(dx[2],0.2e1)*pow(dx[0],0.2e1)*dx[1]/0.3e1-0.2e1*dx[1]*dx[2]*xj*xj*dx[0]-dx[1]*pow(dx[2],0.2e1)*xj*dx[0]-xj*xj*dx[1]*pow(dx[0],0.2e1)+dx[1]*pow(xj,0.3e1)*dx[0]+0.2e1/0.3e1*dx[2]*pow(xj,0.3e1)*dx[0]+pow(dx[2],0.2e1)*pow(xj,0.3e1)/0.3e1+pow(xj,0.3e1)*pow(dx[0],0.2e1)/0.3e1
      B[1] = -pow(dx[1],0.4e1)+0.4e1*xj*pow(dx[1],0.3e1)-0.2e1*pow(dx[1],0.3e1)*dx[0]-0.3e1*xj*xj*pow(dx[1],0.2e1)+0.4e1*pow(dx[1],0.2e1)*xj*dx[2]+0.6e1*xj*dx[0]*pow(dx[1],0.2e1)-pow(dx[0],0.2e1)*pow(dx[1],0.2e1)+pow(dx[1],0.2e1)*pow(dx[2],0.2e1)+0.4e1*dx[1]*dx[2]*xj*dx[0]-0.4e1*dx[1]*dx[2]*xj*xj-0.3e1*dx[1]*xj*xj*dx[0]+0.2e1*dx[1]*xj*pow(dx[0],0.2e1)+dx[1]*pow(dx[2],0.2e1)*dx[0]-xj*xj*pow(dx[0],0.2e1)-0.2e1*dx[2]*xj*xj*dx[0]-pow(dx[2],0.2e1)*xj*xj
      B[2] = -0.2e1*pow(dx[1],0.3e1)-0.2e1*dx[2]*pow(dx[1],0.2e1)+0.3e1*xj*pow(dx[1],0.2e1)-0.3e1*pow(dx[1],0.2e1)*dx[0]+0.4e1*dx[2]*xj*dx[1]-0.2e1*dx[2]*dx[1]*dx[0]-pow(dx[0],0.2e1)*dx[1]+0.3e1*xj*dx[0]*dx[1]+pow(dx[2],0.2e1)*xj+0.2e1*dx[2]*xj*dx[0]+xj*pow(dx[0],0.2e1)
      B[3] = -pow(dx[1],0.2e1)-0.4e1/0.3e1*dx[1]*dx[2]-dx[1]*dx[0]-pow(dx[0],0.2e1)/0.3e1-0.2e1/0.3e1*dx[0]*dx[2]-pow(dx[2],0.2e1)/0.3e1
      B = B*num2
      
      num3 = 0.2e1*pow(dx[2],-0.3e1)/(dx[0]*dx[2]+0.2e1*dx[1]*dx[0]+0.2e1*pow(dx[1],0.2e1)+0.2e1*dx[1]*dx[2])
      C[0] = -pow(xj+dx[2],0.2e1)*(xj*dx[0]*dx[1]+0.4e1*dx[2]*xj*dx[1]+0.3e1*pow(dx[2],0.2e1)*xj+0.2e1*dx[2]*xj*dx[0]+xj*pow(dx[1],0.2e1)-dx[0]*pow(dx[2],0.2e1)-0.2e1*dx[1]*pow(dx[2],0.2e1)-0.2e1*dx[2]*dx[1]*dx[0]-0.2e1*dx[2]*pow(dx[1],0.2e1))/0.3e1
      C[1] = (xj+dx[2])*(xj*pow(dx[1],0.2e1)+0.4e1*dx[2]*xj*dx[1]+xj*dx[0]*dx[1]+0.2e1*dx[2]*xj*dx[0]+0.3e1*pow(dx[2],0.2e1)*xj+pow(dx[2],0.3e1)-dx[2]*pow(dx[1],0.2e1)-dx[2]*dx[1]*dx[0])
      C[2] = -xj*pow(dx[1],0.2e1)-xj*dx[0]*dx[1]-0.2e1*dx[1]*pow(dx[2],0.2e1)-0.4e1*dx[2]*xj*dx[1]-0.2e1*pow(dx[2],0.3e1)-0.2e1*dx[2]*xj*dx[0]-dx[0]*pow(dx[2],0.2e1)-0.3e1*pow(dx[2],0.2e1)*xj
      C[3] = pow(dx[1],0.2e1)/0.3e1+0.4e1/0.3e1*dx[1]*dx[2]+dx[1]*dx[0]/0.3e1+0.2e1/0.3e1*dx[0]*dx[2]+pow(dx[2],0.2e1)
      C = C*num3
    
    elif self.extent[1]==2 and self.extent[2]==0:
      num1 = 0.2e1 * pow(dx[0], -0.3e1)
      A[0] = -pow(xj - dx[0] - dx[1], 0.3e1) / 0.3e1
      A[1] = pow(xj - dx[0] - dx[1], 0.2e1)
      A[2] = -xj + dx[0] + dx[1]
      A[3] = 0.1e1 / 0.3e1
      A = A*num1
      
      num2 = 0.2e1*pow(dx[1],-0.3e1)*pow(dx[0],-0.2e1)
      B[0] = (3*pow(dx[1],2)*pow(dx[0],2)-3*pow(dx[0],2)*xj*dx[1]+pow(dx[0],2)*xj*xj+6*pow(dx[1],3)*dx[0]-9*dx[0]*xj*pow(dx[1],2)+3*dx[0]*xj*xj*dx[1]-6*xj*pow(dx[1],3)+3*pow(dx[1],4)+3*xj*xj*pow(dx[1],2))*xj/3
      B[1] = -pow(dx[1],2)*pow(dx[0],2)+2*pow(dx[0],2)*xj*dx[1]-pow(dx[0],2)*xj*xj-3*dx[0]*xj*xj*dx[1]+6*dx[0]*xj*pow(dx[1],2)-2*pow(dx[1],3)*dx[0]+4*xj*pow(dx[1],3)-pow(dx[1],4)-3*xj*xj*pow(dx[1],2)
      B[2] = xj*pow(dx[0],2)-pow(dx[0],2)*dx[1]-3*dx[0]*pow(dx[1],2)+3*xj*dx[0]*dx[1]-2*pow(dx[1],3)+3*xj*pow(dx[1],2)
      B[3] = -pow(dx[0],2)/3-dx[0]*dx[1]-pow(dx[1],2)
      B = B*num2
    
    elif self.extent[1]==1 and self.extent[2]==0:
      num2 = 2*pow(dx[1],-0.3e1)
      B[0] = pow((-xj+dx[1]),0.3e1)*num2/0.6e1
      B[1] = pow((-xj+dx[1]),0.2e1)*num2/0.2e1
      B[2] =     (-xj+dx[1])       *num2/0.2e1
      B[3] =                        num2/0.6e1
    else:
      pass
    self.A = A[:]
    self.B = B[:]
    self.C = C[:]
    self.D = D[:]
  def val(self,x):
    if self.deriv:
      return self.val_prime(x)
    else:
      coeff = self.coeff(x)
      return Numeric.sum([coeff[i]*x**i for i in range(4)])
  def val_prime(self,x):
    coeff = self.coeff(x)
    return Numeric.sum([i*coeff[i]*x**(i-1) for i in range(1,4)])
  def coeff(self,x):
    coeff = Numeric.zeros([4],typecode=Numeric.Float)
    dx = self.dx
    xj = self.xj
    if x>xj-dx[0]-dx[1] and x<=xj-dx[1]:
      coeff = self.A
    elif x>xj-dx[1] and x<xj:
      coeff = self.B
    elif x==xj:
      if self.extent[1]>0:
        coeff = self.B
      else:
        coeff = self.C
    elif x>=xj and x<=xj+dx[2]:
      coeff = self.C
    elif x>xj+dx[2] and x<xj+dx[3]+dx[2]:
      coeff = self.D
    else:
      pass  
    return coeff
  def shift(self):
    self.xj = self.xj+self.dx[2]
    xj = self.xj
    self.A[:] = self.B[:]
    self.B[:] = self.C[:]
    self.C[:] = self.D[:]
    self.D[:] = 0
    self.dx[0:2] = self.dx[1:3]
    self.dx[3] = 0.
    #self.extent[0:2] = self.extent[0:2]+Numeric.array( [0,1,-1])
    print 'calling shift'
  def __add__(self,spline2):
    self.A+=spline2.A
    self.B+=spline2.B
    self.C+=spline2.C
    self.D+=spline2.D
    return self
  def __sub__(self,spline2):
    self.A = self.A-spline2.A
    self.B = self.B-spline2.B
    self.C = self.C-spline2.C
    self.D = self.D-spline2.D
    return self
  def __div__(self,num):
    self.A = self.A/num
    self.B = self.B/num
    self.C = self.C/num
    self.D = self.D/num
    return self
  def __mul__(self,num):
    self.A=self.A[:]*num
    self.B=self.B[:]*num
    self.C=self.C[:]*num
    self.D=self.D[:]*num
    return self

class linear:
  
  def __init__(self,xj,p2=None,p3=None):
    self.xj = xj
    self.dx = Numeric.zeros([2],typecode=Numeric.Float)
    self.extent = Numeric.zeros([2])
    if p2!=None:
      self.dx[0] = xj-p2
      self.extent[0] = 1
    if p3!=None:
      self.dx[1] = p3-xj
      self.extent[1] = 1
    dx = self.dx
  def val(self,x):
    if x>self.xj-self.dx[0] and x<=self.xj:
      return 1.+(x-self.xj)/self.dx[0]
    elif x<self.xj+self.dx[1] and x>self.xj:
      return 1-(x-self.xj)/self.dx[1]
    else:
      return 0
  def val_prime(self,x):
    if x>self.xj-self.dx[0] and x<=self.xj:
      return 1./self.dx[0]
    elif x<self.xj+self.dx[1] and x>self.xj:
      return -1./self.dx[1]
    else:
      return 0

class linear_fa:
  
  def __init__(self,xj,p2=None,p3=None,p0=0.,Bt0=1.,Bz0=1.,a=1.,equilib=10):
    self.xj = xj
    self.dx = Numeric.zeros([2],typecode=Numeric.Float)
    self.extent = Numeric.zeros([2])
    self.int_fac = Numeric.zeros([2],typecode=Numeric.Float)
    self.p0=p0
    self.Bt0=Bt0
    self.Bz0=Bz0
    self.a=a
    self.equilib=equilib
    if p2!=None:
      self.dx[0] = xj-p2
      self.extent[0] = 1
      self.int_fac[0] = self.val_func(xj)-self.val_func(p2)
    if p3!=None:
      self.dx[1] = p3-xj
      self.extent[1] = 1
      self.int_fac[1] = -self.val_func(xj)+self.val_func(p3)
  def val(self,x):
    if x>self.xj-self.dx[0] and x<=self.xj:
      return (self.val_func(x)-self.val_func(self.xj-self.dx[0]))/self.int_fac[0]
    elif x<self.xj+self.dx[1] and x>self.xj:
      return (self.val_func(self.xj+self.dx[1])-self.val_func(x))/self.int_fac[1]
    else:
      return 0
  def val_prime(self,x):
    if x>self.xj-self.dx[0] and x<=self.xj:
      return (self.val_prime_func(x)/self.int_fac[0])
    elif x<self.xj+self.dx[1] and x>self.xj:
      return (-self.val_prime_func(x)/self.int_fac[1])
    else:
      return 0
  def val_func(self,r):
    Bt0 = self.Bt0
    Bz0 = self.Bz0
    p0 = self.p0
    a = self.a
    if self.equilib==10:
      return (Bt0*math.atan(math.sqrt(2*Bt0**2-2*p0)*r\
      /math.sqrt(a**2*(Bz0**2-2*p0+2*Bt0**2)+r**2*2*(p0-Bt0**2)))\
      *a**2*(Bz0**2-2*p0+2*Bt0**2)/4./(Bt0**2-p0)/math.sqrt(2*Bt0**2-2*p0)\
      +Bt0*r*math.sqrt(2*(p0-Bt0**2)*r**2+a**2*(Bz0**2-2*p0+2*Bt0**2))/4./(p0-Bt0**2))
    elif self.equilib==3:
      return (1/3.*Bt0*r**3/math.sqrt(a**2*Bz0**2))
    elif self.equilib==13:
      return -Bt0/48./a**4/Bz0\
                            *(15*a**6*math.asin(r/a)-48*a**5*r\
                              +math.sqrt(a**2-r**2)*r\
                                *(33*a**4-26*a**2*r**2+8*r**4))
  def val_prime_func(self,r):
    Bt0 = self.Bt0
    Bz0 = self.Bz0
    p0 = self.p0
    a = self.a
    if self.equilib==10:
      return (Bt0*r**2/math.sqrt(a**2*(Bz0**2-2*p0+2*Bt0**2)+2*r**2*(p0-Bt0**2)))
    elif self.equilib==3:
      return (Bt0*r**2/math.sqrt(a**2*Bz0**2))
    elif self.equilib==13:
      return Bt0*(a-(a**2-r**2)**2*math.sqrt(1-r**2/a**2)/a**3)/Bz0

class const:

  def __init__(self,xj,p3):
    self.xj = xj
    self.dx = p3-xj
    self.extent = [0,1]
    
  def val(self,x):
    if x>=self.xj and x<self.xj+self.dx:
      return 1.
    else:
      return 0.  
  def val_prime(self,x):
    return 0.  
if __name__=='__main__':
  N = 20
  grid = Numeric.zeros([N],typecode=Numeric.Float)
  for i in range(len(grid)):
    grid[i] = i/float(N-1)
  #print grid
  phi = []
  phi0 = bspline(grid[0],p3=grid[1])
  phi1 = bspline(grid[0],p3=grid[1],p4=grid[2])
  phi2 = bspline(grid[0],p3=grid[1],p4=grid[2])
  phi3 = bspline(grid[1],p2=grid[0],p3=grid[2],p4=grid[3])
  phi1 = phi1*(phi0.C[1]/phi1.C[1])-phi0
  phi2 = phi2*(phi0.C[2]/phi2.C[2])-phi0
  phi3 = phi3*(phi0.C[2]/phi3.B[2])
  phi3.B[:] = phi3.B-phi0.C
  phi.append(phi1)
  phi.append(phi2)
  phi.append(phi3)
  for i in range(2,N-2):
    phi.append(bspline(grid[i],p1=grid[i-2],p2=grid[i-1],p3=grid[i+1],p4=grid[i+2]))
  phi.append(bspline(grid[N-2],p1=grid[N-4],p2=grid[N-3],p3=grid[N-1]))
  phi.append(bspline(grid[N-1],p1=grid[N-3],p2=grid[N-2]))
  for i in range(len(phi)):
    phi[i].deriv=True
  #print phi[0].A, phi[0].extent
  xi = Numeric.zeros([len(phi)],typecode=Numeric.Float)
  for i in range(len(phi)):
    xi[i] = i/float(N-1)
  xi[0] = 0
  xi[1]=xi[1]/2.
  xi[-1]=xi[-1]/2.
  import Gnuplot
  x = Numeric.arange(0,1,(1/2.)**9)
  y = Numeric.zeros([len(x),len(phi)],typecode=Numeric.Float)
  for j in range(len(phi)):
    for i in range(len(x)):
      y[i,j] = phi[j].val(x[i])
  g = Gnuplot.Gnuplot(persist=1)
  g._clear_queue()
  for j in range(3):#len(phi)):
    g._add_to_queue([Gnuplot.Data(x,y[:,j])])
  g('set xrange [0:.1]')
  g.refresh()
  g1=Gnuplot.Gnuplot(persist=1)
  g1._clear_queue()
  z = Numeric.zeros([len(x)],typecode=Numeric.Float)
  for i in range(len(x)):
    z[i] = Numeric.sum([xi[j]*phi[j].val(x[i]) for j in range(len(phi))])
  #g1.plot(Gnuplot.Data(x,z))
  
