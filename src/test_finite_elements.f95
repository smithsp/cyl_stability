PROGRAM test_finite_elements
  USE local
  USE finite_elements_module
  USE cyl_funcs_module
  USE cyl_matrix_module
  USE sort_module
  IMPLICIT NONE
  INTEGER, PARAMETER :: N=20
  TYPE(linear), DIMENSION(N) :: phi,phi_mod
  TYPE(constant), DIMENSION(N-1) :: psi, chi
  TYPE(bspline), DIMENSION(0:N+1) :: xi, xi_deriv
  TYPE(bspline) :: xi1, xi2, xi3, xi0
  REAL(r8), DIMENSION(N) :: grid
  REAL(r8) :: temp, dx1, dx2, Bz1, ptemp(1), xj, dx(5), x(300*N)
  INTEGER :: i, nx
  LOGICAL :: Lend0, Rend0
  INTERFACE
    FUNCTION f(x)
      USE local
      REAL(r8), DIMENSION(:) :: x
      REAL(r8), DIMENSION(size(x)) :: f
    END FUNCTION f
  END INTERFACE
  nx = size(x)
  ar = 1.0
  kz = .1
  gamma = 5./3.
  mt = 2
  eps = 0.
  Bz0 = 15.
!  equilib = 4
  P0  = .0
  P1 = .50
  lambd = 2.
  rho0 = 1.
  Bt0 = 1.
  Bz0 = 15
  equilib=10
  grid = (/ (i*ar/real(N-1), i=0,N-1) /)
  grid = new_grid(grid)
  grid(1) = grid(2)/1000.0
  alpha = 1.0
  rs = .77
  !kB = .true.
  CALL calc_rs(grid)
  !grid = new_grid(grid)
  WRITE (*,*) 'rs = ', rs
  WRITE (*,*) sort(grid)
  !STOP
  
  x = (/( i*(maxval(grid)-minval(grid))/real(nx-1) + minval(grid), i=0,nx-1)/)
  Lend0 = .true.
  Rend0 = .true.
  CALL init(phi(2:N-1),grid(2:N-1),p2=grid(1:N-2),p3=grid(3:N))
  CALL init(phi(1),grid(1),p3=grid(2))
  CALL init(phi(N),grid(N),p2=grid(N-1))
  CALL init(xi(1),grid(1),p3=grid(2), p4=grid(3),LendZero=Lend0)
  CALL init(xi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=Lend0)
  CALL init(xi(0),grid(1),p3=grid(2),LendZero=.true.)
  !CALL init(xi(1),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.)
  !CALL init(xi(2),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.)  
  CALL init(xi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
  CALL init(xi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2), p3=grid(N),RendZero=Rend0)
  CALL init(xi(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=Rend0)
  CALL init(xi(N+1),grid(N),p2=grid(N-1),RendZero=Rend0)
  xi_deriv = xi
  xi_deriv%deriv = .true.
  phi_mod = phi
  phi_mod%mod_lin = .true.
  CALL linear_int_fac(phi_mod)
  xi1 = xi(1)*xi(0)%C(1)/xi(1)%C(1)-xi(0)
  xi1%deriv = .true.
  xi0 = shift_spline(xi_deriv(2))
  xi0 = xi_deriv(3)+xi0
  !xi0%C=(/0.,0.,-3*xi0%dx(3),2./)
  xi2 = xi(1)*xi(0)%C(2)/xi(1)%C(2)-xi(0)
  xi2%deriv = .true.
  xi3 = shift_spline(xi_deriv(0))
  !write (0,*) xi3
  xi3 = xi_deriv(2)*(2*xi3%B(2)-6*xi3%B(3)*xi3%xj)/(2*xi_deriv(2)%B(2)-6*xi_deriv(2)%B(3)*xi_deriv(2)%xj)-xi3
  xi3%deriv = .true.
  !write (0,*) xi1
  !write (0,*) xi2
  !write (0,*) xi3
  !WRITE (0,'(4g)') xi(0)%C
  !WRITE (0,'(4g)') xi(1)%A,xi(1)%B,xi(1)%C,xi(1)%D
  !WRITE (0,'(4g)') xi1%A,xi1%B,xi1%C,xi1%D
  !WRITE (0,'(4g)') xi2%A,xi2%B,xi2%C,xi2%D
  !WRITE (0,'(4g)') xi3%A,xi3%B,xi3%C,xi3%D
  WRITE (*,*) val(xi(N-1),grid(N-4:N))
  CALL init(psi(1:N-1),grid(1:N-1),p3=grid(2:N))
  WRITE (*,*) 'grid = ', grid
  WRITE (*,*) phi%dx(2)
  !WRITE (*,*) linear_val(phi,grid-.2)  
  WRITE (*,*) constant_val(psi,grid(2:N)-.5)
  WRITE (*,'(a)') 'bspline value at left endpoint = '
  WRITE (*,'(g)') val(xi(2:4),grid(1))
  WRITE (*,'(a)') 'bspline value at right endpoint = '
  WRITE (*,'(g)')  val(xi(N-5:N+1),grid(N))
  WRITE (*,'(a)') 'bspline_deriv value at right endpoint = '
  WRITE (*,'(g)')  val(xi_deriv(N-5:N+1),grid(N))
  
  WRITE (*,*) 'xi(0)%error = ', xi(0)%error, 'xi(0)%xj = ', xi(0)%xj
  WRITE (*,'(a,4g)')   'xi(0)%dx = ',xi(0)%dx,'xi(0)%A = ', xi(0)%A, 'xi(0)%B = ', xi(0)%B,&
  &'xi(0)%C = ', xi(0)%C, 'xi(0)%D = ', xi(0)%D
  
  WRITE (*,*) 'xi(1)%error = ', xi(1)%error, 'xi(1)%xj = ', xi(1)%xj
  WRITE (*,'(a,4g)')   'xi(1)%dx = ',xi(1)%dx,'xi(1)%A = ', xi(1)%A, 'xi(1)%B = ', xi(1)%B,&
  &'xi(1)%C = ', xi(1)%C, 'xi(1)%D = ', xi(1)%D
  
  WRITE (*,*) 'xi(2)%error = ', xi(2)%error, 'xi(2)%xj = ', xi(2)%xj
  WRITE (*,'(a,4g)')   'xi(2)%dx = ',xi(2)%dx,'xi(2)%A = ', xi(2)%A, 'xi(2)%B = ', xi(2)%B,&
  &'xi(2)%C = ', xi(2)%C, 'xi(2)%D = ', xi(2)%D
  
  WRITE (*,*) 'xi_deriv(2)%error = ', xi_deriv(2)%error, 'xi_deriv(2)%xj = ', xi_deriv(2)%xj
  WRITE (*,'(a,4g)')   'xi_deriv(2)%dx = ',xi_deriv(2)%dx,'xi_deriv(2)%A = ', xi_deriv(2)%A, 'xi_deriv(2)%B = ', xi_deriv(2)%B,&
  &'xi_deriv(2)%C = ', xi_deriv(2)%C, 'xi_deriv(2)%D = ', xi_deriv(2)%D
  WRITE (*,*) 'val(xi_deriv(2),1e-5) = ', val(xi_deriv(2),1e-5)
  
  WRITE (*,*) 'xi(N/2)%error = ', xi(N/2)%error, 'xi(N/2)%xj = ', xi(N/2)%xj
  WRITE (*,'(a,4g)')   'xi(N/2)%dx = ',xi(N/2)%dx,'xi(N/2)%A = ', xi(N/2)%A, 'xi(N/2)%B = ', xi(N/2)%B,&
  &'xi(N/2)%C = ', xi(N/2)%C, 'xi(N/2)%D = ', xi(N/2)%D
  
  WRITE (*,*) 'xi(N-1)%error = ',xi(N-1)%error, 'xi(N-1)%xj = ', xi(N-1)%xj
  WRITE (*,'(a,4g)')   'xi(N-1)%dx = ',xi(N-1)%dx,'xi(N-1)%A = ', xi(N-1)%A, 'xi(N-1)%B = ', xi(N-1)%B,&
  &'xi(N-1)%C = ', xi(N-1)%C, 'xi(N-1)%D = ', xi(N-1)%D
  
  WRITE (*,*) 'xi(N)%error = ',xi(N)%error, 'xi(N)%xj = ', xi(N)%xj
  WRITE (*,'(a,4g)')   'xi(N)%dx = ',xi(N)%dx,'xi(N)%A = ', xi(N)%A, 'xi(N)%B = ', xi(N)%B,&
  &'xi(N)%C = ', xi(N)%C, 'xi(N)%D = ', xi(N)%D
  
  WRITE (*,*) 'xi(N+1)%error = ',xi(N+1)%error, 'xi(N+1)%xj = ', xi(N+1)%xj
  WRITE (*,'(a,4g)')   'xi(N+1)%dx = ',xi(N+1)%dx,'xi(N+1)%A = ', xi(N+1)%A, 'xi(N+1)%B = ', xi(N+1)%B,&
  &'xi(N+1)%C = ', xi(N+1)%C, 'xi(N+1)%D = ', xi(N+1)%D
  WRITE(*,*) 'val(xi_deriv(N),ar)=',val(xi_deriv(N),ar), 'val(xi(N+1),ar)=', val(xi(N+1),ar), 'val(xi_deriv(N+1),ar)=',val(xi_deriv(N+1),ar)
  
  OPEN (1, status='replace',file='finite_element_values.txt')
  WRITE (1,'(31a20)') &
  & '#x value', 'Left constant', 'Right constant',&
  & 'Left Linear', 'Middle Linear', 'Right Linear',&
  & 'Left1 Bspline', 'Left2 Bspline', 'Left3 Bspline','RightN-1 Bspline',&
  & 'Left1 SplDeriv', 'Left2 SplDeriv', 'Left3 SplDeriv', 'RightN-1 SplDeriv',&
  & 'Left0 SplDeriv', 'Left0 Bspline', &
  & 'RightN-2 Bspline', 'RightN Bspline', 'RightN+1 Bspline', &
  & 'RightN-2 SplDeriv', 'RightN SplDeriv', 'RightN+1 SplDeriv',&
  & 'Left1 SplDeriv val0', 'Left1 SplDeriv drv0', &
  & 'Left2 SplDeriv drv0', 'Left0 SplDeriv drv0', &
  & 'Left0 mod_lin', 'Left1 mod_lin', &
  & 'Left0 mod_lin_prime', 'Left1 mod_lin_prime', 'Left1 LinDerivBt/Bz'
  
  DO i=1,nx
    WRITE (1,'(31g20.7)') &
    & x(i), val(psi(1),x(i)), val(psi(size(psi)),x(i)), &
    & val(phi(1),(/x(i)/)), val(phi(2),(/x(i)/)), val(phi(size(phi)),(/x(i)/)), &
    & val(xi(1),x(i)), val(xi(2),x(i)), val(xi(3),x(i)), val(xi(N-1), x(i)), &
    & val(xi_deriv(1),x(i)), val(xi_deriv(2),x(i)),  val(xi_deriv(3),x(i)), val(xi_deriv(N-1),x(i)), &
    & val(xi_deriv(0),x(i)), val(xi(0),x(i)), &
    & val(xi(N-2),x(i)), val(xi(N),x(i)), val(xi(N+1),x(i)), &
    & val(xi_deriv(N-2),x(i)), val(xi_deriv(N),x(i)), val(xi_deriv(N+1),x(i)), &
    & val(xi1,x(i)), val(xi2,x(i)), val(xi3,x(i)), val(xi0,x(i)), &
    & val(phi_mod(1),(/x(i)/)), val(phi_mod(2),(/x(i)/)),&
    & val_prime(phi_mod(1),(/x(i)/)), val_prime(phi_mod(2),(/x(i)/)), val_prime(phi(2),(/x(i)/))*mod_lin_func((/x(i)/))
    
    
  ENDDO
  CLOSE(1)
  epsilo = 1.e-14
  DO i=1,size(phi_mod)
    WRITE (*,*) phi_mod(i)%int_fac(1),'=',(int_mod_lin_func(phi_mod(i)%xj)-int_mod_lin_func(phi_mod(i)%xj-phi_mod(i)%dx(1)))
  ENDDO
  
  WRITE (0,*) (val_prime(phi_mod(N),(/ar/))),(val_prime(phi_mod(N-1),(/0.,ar/))),ar-phi_mod(N-1)%xj,phi_mod(N-1)%dx(2),ar-phi_mod(N)%xj
  
  verbose = .true.
  WRITE (*,*) 'val(phi_mod(2),(/0,.01,.02,.03/)=',val(phi_mod(2),(/0.,.01,.02,.03/)), 'sum(val(phi_mod(2),(/0,.01,.02,.03/)))=',sum(val(phi_mod(2),(/0.,.01,.02,.03/)))
  WRITE (*,*) 'int_func(phi_mod(2),phi_mod(2),x2)=',int_func(phi_mod(2),phi_mod(2),x2)
  WRITE (*,*) 'int_func(phi(2),phi(2),x2)=',int_func(phi(2),phi(2),x2)
  STOP
  WRITE (*,*) int_func(psi(1),psi(1),f)
  WRITE (*,*) int_func(psi(1),psi(2),f)
  WRITE (*,*) int_func(psi(1),psi(1),W11A)
  WRITE (*,*) int_func(psi(1),phi(2),f)
  WRITE (*,*) int_func(psi(1),phi(1),f)
  WRITE (*,*) int_func(psi(1),phi(3),f)
  WRITE (*,*) int_func(phi(1),phi(1),f)
  WRITE (*,*) int_func(phi(1),phi(2),f)
  WRITE (*,*) int_func(phi(2),phi(2),f)
  WRITE (*,*) int_func(psi(1),phi(2),f,deriv2=.true.)
  WRITE (*,*) int_func(psi(2),phi(2),f,deriv2=.true.)
  WRITE (*,*) int_func(phi(2),phi(2),f,deriv1=.true.)
  WRITE (*,*) int_func(phi(2),phi(2),f,deriv2=.true.)
  WRITE (*,*) int_func(phi(2),phi(2),f,deriv1=.true.,deriv2=.true.)
  WRITE (*,*) int_func(phi(1),phi(2),f,deriv1=.true.,deriv2=.true.)
  WRITE (*,*) int_func(phi(1),phi(2),f,deriv2=.true.)
  WRITE (*,*) int_func(phi(2),phi(2),W11B,deriv1=.true.)
  
  WRITE (*,*) '------------------------------------------------------------------'
  WRITE (*,*) 'Assuming f(r)=1'
  dx(2:4) = grid(2:4)-grid(1:3)
  xj = grid(2)
  WRITE (*,'(a,g,a,g)') 'int_func(xi(2),xi(2),f):  ',int_func(xi(2),xi(2),f), ' = ', &
  &-dble(1/dx(2)**2*(-20*dx(3)**3*dx(2)*dx(4)**5-704*xj*dx(3)**6*dx(2)*dx(4)-210*dx(3)**3*dx(2)**2*dx(4)**4-&
  &56*dx(3)**2*dx(2)**4*dx(4)**3-140*dx(3)**2*dx(2)**3*dx(4)**4+28*dx(2)**7*dx(4)**2-264*dx(2)**6*dx(3)*dx(4)*xj-&
  &1904*xj*dx(3)**5*dx(4)*dx(2)**2-1048*xj*dx(3)**5*dx(4)**2*dx(2)-2424*xj*dx(3)**4*dx(4)**2*dx(2)**2-&
  &680*xj*dx(3)**4*dx(4)**3*dx(2)-1300*xj*dx(3)**3*dx(4)**3*dx(2)**2-160*xj*dx(3)**3*dx(4)**4*dx(2)-&
  &2300*xj*dx(2)**4*dx(3)**3*dx(4)-2912*xj*dx(2)**3*dx(3)**3*dx(4)**2-2760*xj*dx(2)**3*dx(3)**4*dx(4)-&
  &1088*xj*dx(2)**5*dx(3)**2*dx(4)-1912*xj*dx(2)**4*dx(3)**2*dx(4)**2-1224*xj*dx(2)**3*dx(3)**2*dx(4)**3-&
  &240*xj*dx(2)**2*dx(3)**2*dx(4)**4-648*xj*dx(2)**5*dx(3)*dx(4)**2-568*xj*dx(2)**4*dx(3)*dx(4)**3-&
  &160*xj*dx(2)**3*dx(3)*dx(4)**4-104*dx(2)**5*xj*dx(4)**3-176*dx(2)**6*dx(3)**2*xj-24*dx(2)**7*dx(3)*xj+&
  &268*dx(2)**6*dx(3)**2*dx(4)+84*dx(2)**7*dx(3)*dx(4)-88*dx(2)**6*dx(4)**2*xj-24*dx(2)**7*dx(4)*xj+&
  &158*dx(2)**6*dx(4)**2*dx(3)-28*dx(3)**3*dx(2)**4*dx(4)**2-420*dx(3)**3*dx(2)**3*dx(4)**3+9*dx(2)**8*dx(3)-&
  &35*dx(4)**4*dx(3)**5-82*dx(4)**3*dx(3)**6-5*dx(4)**5*dx(3)**4-45*dx(3)**8*dx(4)-88*dx(3)**7*dx(4)**2+&
  &134*dx(2)**6*dx(3)**3+56*dx(2)**7*dx(3)**2-5*dx(2)**4*dx(4)**5-24*xj*dx(3)**8-136*dx(3)**6*dx(2)**3-&
  &134*dx(3)**7*dx(2)**2-56*dx(3)**8*dx(2)+136*dx(2)**5*dx(3)**4+9*dx(2)**8*dx(4)+24*dx(2)**6*dx(4)**3-9*dx(3)**9-&
  &108*xj*dx(3)**7*dx(4)-184*xj*dx(3)**6*dx(4)**2-476*dx(3)**5*dx(2)**3*dx(4)-826*dx(3)**5*dx(2)**2*dx(4)**2-&
  &536*dx(3)**6*dx(2)**2*dx(4)-252*dx(3)**7*dx(2)*dx(4)-440*dx(3)**6*dx(2)*dx(4)**2-364*dx(3)**5*dx(2)*dx(4)**3-&
  &644*dx(3)**4*dx(2)**3*dx(4)**2-602*dx(3)**4*dx(2)**2*dx(4)**3-140*dx(3)**4*dx(2)*dx(4)**4+&
  &340*dx(2)**5*dx(3)**3*dx(4)+272*dx(2)**5*dx(3)**2*dx(4)**2-920*xj*dx(2)**4*dx(3)**4-920*xj*dx(2)**3*dx(3)**5-&
  &544*xj*dx(2)**5*dx(3)**3-40*xj*dx(2)**4*dx(4)**4-30*dx(3)**2*dx(2)**2*dx(4)**5+68*dx(3)*dx(2)**5*dx(4)**3-&
  &35*dx(3)*dx(2)**4*dx(4)**4-20*dx(3)*dx(2)**3*dx(4)**5-140*xj*dx(3)**5*dx(4)**3-40*xj*dx(3)**4*dx(4)**4-&
  &544*xj*dx(3)**6*dx(2)**2-176*xj*dx(3)**7*dx(2))/(2*dx(3)**2+2*dx(3)*dx(4)+2*dx(3)*dx(2)+dx(2)*dx(4))/&
  &(2*dx(3)**3+2*dx(2)*dx(3)**2+4*dx(4)*dx(3)**2+3*dx(3)*dx(2)*dx(4)+2*dx(4)**2*dx(3)+dx(4)**2*dx(2)))/0.630D3

  DO i=3,N-2
    WRITE (*,*) 'i = ', i
    dx(1:4) = grid(i-1:i+2)-grid(i-2:i+1)
    xj = grid(i)
    WRITE (*,'(a,g,a,g)') 'int_func(xi(i),xi(i),f):  ', int_func(xi(i),xi(i),f), ' = ', &
    &dble((dx(2)+dx(3))*(dx(3)+dx(2)+dx(1))*(dx(3)+dx(4)+dx(2))*(dx(3)+dx(4)+dx(1)+dx(2))*(((5*dx(2)+5*dx(1))*dx(3)**2+&
    &(15*dx(1)*dx(2)+10*dx(2)**2+5*dx(1)**2)*dx(3)+5*dx(2)**3+10*dx(1)*dx(2)**2+5*dx(2)*dx(1)**2)*dx(4)**3+&
    &((25*dx(2)+25*dx(1))*dx(3)**3+(40*dx(2)**2+(40*xj+60*dx(1))*dx(2)+40*xj*dx(1)+20*dx(1)**2)*dx(3)**2+&
    &(5*dx(2)**3+(80*xj+10*dx(1))*dx(2)**2+120*dx(2)*dx(1)*xj+40*xj*dx(1)**2-5*dx(1)**3)*dx(3)-10*dx(2)**4+(-25*dx(1)+40*xj)*dx(2)**3+&
    &(-20*dx(1)**2+80*xj*dx(1))*dx(2)**2+(40*xj*dx(1)**2-5*dx(1)**3)*dx(2))*dx(4)**2+&
    &((27*dx(2)+27*dx(1))*dx(3)**4+(50*dx(2)**2+(60*xj+75*dx(1))*dx(2)+60*xj*dx(1)+25*dx(1)**2)*dx(3)**3+&
    &(160*dx(2)**2*xj+(240*xj*dx(1)-10*dx(1)**2)*dx(2)+80*xj*dx(1)**2-10*dx(1)**3)*dx(3)**2+&
    &(-30*dx(2)**4+(120*xj-75*dx(1))*dx(2)**3+(-60*dx(1)**2+240*xj*dx(1))*dx(2)**2+(-15*dx(1)**3+120*xj*dx(1)**2)*dx(2))*dx(3)-&
    &9*dx(2)**5+(-27*dx(1)+24*xj)*dx(2)**4+(60*xj*dx(1)-25*dx(1)**2)*dx(2)**3+(40*xj*dx(1)**2-5*dx(1)**3)*dx(2)**2)*dx(4)+&
    &(9*dx(1)+9*dx(2))*dx(3)**5+(20*dx(2)**2+(30*dx(1)+24*xj)*dx(2)+10*dx(1)**2+24*xj*dx(1))*dx(3)**4+&
    &(80*dx(2)**2*xj+(120*xj*dx(1)-5*dx(1)**2)*dx(2)+40*xj*dx(1)**2-5*dx(1)**3)*dx(3)**3+&
    &(-20*dx(2)**4+(80*xj-50*dx(1))*dx(2)**3+(-40*dx(1)**2+160*xj*dx(1))*dx(2)**2+(80*xj*dx(1)**2-10*dx(1)**3)*dx(2))*dx(3)**2+&
    &(-9*dx(2)**5+(-27*dx(1)+24*xj)*dx(2)**4+(60*xj*dx(1)-25*dx(1)**2)*dx(2)**3+(40*xj*dx(1)**2-5*dx(1)**3)*dx(2)**2)*dx(3))/&
    &(dx(3)+dx(4))/(((dx(1)+2*dx(2))*dx(3)+dx(1)*dx(2)+dx(2)**2)*dx(4)+(dx(1)+2*dx(2))*dx(3)**2+(2*dx(1)*dx(2)+2*dx(2)**2)*dx(3))**2/&
    &(dx(1)+dx(2)))/0.630D3
    
    WRITE  (*,'(a, g,a,g)') 'int_func(xi(i),xi(i),f,deriv1=.true.):  ', int_func(xi(i),xi(i),f,deriv1=.true.), ' = ',&
    &-dble((dx(2)+dx(3))*(dx(3)+dx(2)+dx(1))*(dx(3)+dx(4)+dx(2))*(dx(3)+dx(4)+dx(1)+dx(2))*(((10*dx(2)+10*dx(1))*dx(3)**2+&
    &(30*dx(1)*dx(2)+20*dx(2)**2+10*dx(1)**2)*dx(3)+10*dx(2)**3+10*dx(2)*dx(1)**2+20*dx(1)*dx(2)**2)*dx(4)**2+&
    &((15*dx(2)+15*dx(1))*dx(3)**3+(20*dx(1)**2+60*dx(1)*dx(2)+40*dx(2)**2)*dx(3)**2+&
    &(30*dx(2)**3+30*dx(2)*dx(1)**2+60*dx(1)*dx(2)**2)*dx(3)+&
    &10*dx(2)**2*dx(1)**2+15*dx(2)**3*dx(1)+6*dx(2)**4)*dx(4)+(6*dx(2)+6*dx(1))*dx(3)**4+(30*dx(1)*dx(2)+20*dx(2)**2+10*dx(1)**2)*dx(3)**3+&
    &(20*dx(2)**3+40*dx(1)*dx(2)**2+20*dx(2)*dx(1)**2)*dx(3)**2+(10*dx(2)**2*dx(1)**2+15*dx(2)**3*dx(1)+6*dx(2)**4)*dx(3))/&
    &(dx(3)+dx(4))/(((dx(1)+2*dx(2))*dx(3)+dx(1)*dx(2)+dx(2)**2)*dx(4)+(dx(1)+2*dx(2))*dx(3)**2+&
    &(2*dx(1)*dx(2)+2*dx(2)**2)*dx(3))**2/(dx(1)+dx(2)))/0.315D3
    
    WRITE  (*,'(a,g,a,g)') 'int_func(xi(i),xi(i),f,deriv1=.true.,deriv2=.true.):  ',int_func(xi(i),xi(i),f,deriv1=.true.,deriv2=.true.), ' = ',&
    &0.2D1/0.15D2*dble(dx(2)+dx(3))*dble(dx(3)+dx(2)+dx(1))*dble(dx(3)+dx(4)+dx(2))*dble(dx(3)+dx(4)+dx(1)+dx(2))*&
    &dble(((2*dx(2)+2*dx(1))*dx(3)**2+(6*xj*dx(2)+6*xj*dx(1))*dx(3)-dx(2)**3+(4*xj-2*dx(1))*dx(2)**2+6*dx(2)*dx(1)*xj)*dx(4)+&
    &(dx(1)+dx(2))*dx(3)**3+(4*xj*dx(1)+4*xj*dx(2))*dx(3)**2+(-dx(2)**3+(4*xj-2*dx(1))*dx(2)**2+6*dx(2)*dx(1)*xj)*dx(3))/&
    &dble(dx(3)+dx(4))/dble((((dx(1)+2*dx(2))*dx(3)+dx(1)*dx(2)+dx(2)**2)*dx(4)+(dx(1)+2*dx(2))*dx(3)**2+&
    &(2*dx(1)*dx(2)+2*dx(2)**2)*dx(3))**2)/dble(dx(1)+dx(2))
  ENDDO
  DO i=3,N-3
    WRITE (*,*) 'i = ', i
    dx = grid(i-1:i+3)-grid(i-2:i+2)
    xj = grid(i)
    WRITE (*,'(a, g,a,g)') 'int_func(xi(i),xi(i+1),f,deriv2=.true.):  ', int_func(xi(i),xi(i+1),f,deriv2=.true.), ' = ',&
    &dble((dx(3)+dx(4)+dx(2))*(((-12*dx(5)-12*dx(4))*dx(3)**2+(-12*dx(5)**2-24*dx(4)**2-36*dx(4)*dx(5))*dx(3)-&
    &12*dx(4)**3-24*dx(4)**2*dx(5)-12*dx(4)*dx(5)**2)*dx(2)**5+((-20*dx(5)-20*dx(4))*dx(3)**3+&
    &(-50*dx(4)**2+(-70*dx(5)+42*xj-33*dx(1))*dx(4)-33*dx(1)*dx(5)-20*dx(5)**2+42*dx(5)*xj)*dx(3)**2+&
    &(-40*dx(4)**3+(84*xj-66*dx(1)-70*dx(5))*dx(4)**2+(126*dx(5)*xj-99*dx(1)*dx(5)-30*dx(5)**2)*dx(4)-&
    &33*dx(1)*dx(5)**2+42*dx(5)**2*xj)*dx(3)-10*dx(4)**4+(-33*dx(1)+42*xj-20*dx(5))*dx(4)**3+&
    &(-10*dx(5)**2-66*dx(1)*dx(5)+84*dx(5)*xj)*dx(4)**2+(-33*dx(1)*dx(5)**2+42*dx(5)**2*xj)*dx(4))*dx(2)**4+&
    &((2*dx(5)+2*dx(4))*dx(3)**4+(15*dx(4)**2+(35*dx(5)+105*xj-40*dx(1))*dx(4)+105*dx(5)*xj-40*dx(1)*dx(5)+20*dx(5)**2)*dx(3)**3+&
    &(10*dx(4)**3+(280*xj-100*dx(1)+50*dx(5))*dx(4)**2+(40*dx(5)**2+105*xj*dx(1)-140*dx(1)*dx(5)-25*dx(1)**2+420*dx(5)*xj)*dx(4)+&
    &140*dx(5)**2*xj+105*dx(1)*dx(5)*xj-40*dx(1)*dx(5)**2-25*dx(1)**2*dx(5))*dx(3)**2+(2*dx(4)**4+&
    &(210*xj+25*dx(5)-80*dx(1))*dx(4)**3+(-50*dx(1)**2-140*dx(1)*dx(5)+210*xj*dx(1)+420*dx(5)*xj+30*dx(5)**2)*dx(4)**2+&
    &(210*dx(5)**2*xj-75*dx(1)**2*dx(5)+315*dx(1)*dx(5)*xj-60*dx(1)*dx(5)**2)*dx(4)-25*dx(1)**2*dx(5)**2+&
    &105*dx(1)*dx(5)**2*xj)*dx(3)+2*dx(4)**5+(-20*dx(1)+42*xj+9*dx(5))*dx(4)**4+&
    &(10*dx(5)**2+105*xj*dx(1)-40*dx(1)*dx(5)+105*dx(5)*xj-25*dx(1)**2)*dx(4)**3+&
    &(70*dx(5)**2*xj-50*dx(1)**2*dx(5)+210*dx(1)*dx(5)*xj-20*dx(1)*dx(5)**2)*dx(4)**2+&
    &(105*dx(1)*dx(5)**2*xj-25*dx(1)**2*dx(5)**2)*dx(4))*dx(2)**3+((16*dx(5)+16*dx(4))*dx(3)**5+(71*dx(4)**2+&
    &(84*xj+117*dx(5)+14*dx(1))*dx(4)+14*dx(1)*dx(5)+46*dx(5)**2+84*dx(5)*xj)*dx(3)**4+&
    &(85*dx(4)**3+(60*dx(1)+280*xj+200*dx(5))*dx(4)**2+(210*xj*dx(1)+420*dx(5)*xj+110*dx(1)*dx(5)-20*dx(1)**2+115*dx(5)**2)*dx(4)+&
    &50*dx(1)*dx(5)**2-20*dx(1)**2*dx(5)+210*dx(1)*dx(5)*xj+140*dx(5)**2*xj)*dx(3)**3+&
    &(34*dx(4)**4+(110*dx(5)+50*dx(1)+280*xj)*dx(4)**3+(150*dx(1)*dx(5)+560*dx(5)*xj+90*dx(5)**2-50*dx(1)**2+560*xj*dx(1))*dx(4)**2+&
    &(-70*dx(1)**2*dx(5)+840*dx(1)*dx(5)*xj+70*xj*dx(1)**2+100*dx(1)*dx(5)**2+280*dx(5)**2*xj)*dx(4)+70*dx(1)**2*dx(5)*xj+&
    &280*dx(1)*dx(5)**2*xj-20*dx(1)**2*dx(5)**2)*dx(3)**2+&
    &(4*dx(4)**5+(18*dx(5)+14*dx(1)+84*xj)*dx(4)**4+(-40*dx(1)**2+210*dx(5)*xj+70*dx(1)*dx(5)+20*dx(5)**2+420*xj*dx(1))*dx(4)**3+&
    &(140*xj*dx(1)**2+840*dx(1)*dx(5)*xj+70*dx(1)*dx(5)**2+140*dx(5)**2*xj-70*dx(1)**2*dx(5))*dx(4)**2+&
    &(210*dx(1)**2*dx(5)*xj+420*dx(1)*dx(5)**2*xj-30*dx(1)**2*dx(5)**2)*dx(4)+70*dx(1)**2*dx(5)**2*xj)*dx(3)+&
    &4*dx(4)**5*dx(1)+(84*xj*dx(1)-10*dx(1)**2+18*dx(1)*dx(5))*dx(4)**4+&
    &(-20*dx(1)**2*dx(5)+70*xj*dx(1)**2+210*dx(1)*dx(5)*xj+20*dx(1)*dx(5)**2)*dx(4)**3+&
    &(140*dx(1)*dx(5)**2*xj+140*dx(1)**2*dx(5)*xj-10*dx(1)**2*dx(5)**2)*dx(4)**2+70*dx(1)**2*dx(5)**2*xj*dx(4))*dx(2)**2+&
    &((6*dx(5)+6*dx(4))*dx(3)**6+(30*dx(4)**2+(27*dx(1)+21*xj+48*dx(5))*dx(4)+21*dx(5)*xj+18*dx(5)**2+27*dx(1)*dx(5))*dx(3)**5+&
    &(44*dx(4)**3+(84*xj+98*dx(5)+117*dx(1))*dx(4)**2+(54*dx(5)**2+126*dx(5)*xj+12*dx(1)**2+189*dx(1)*dx(5)+126*xj*dx(1))*dx(4)+&
    &126*dx(1)*dx(5)*xj+42*dx(5)**2*xj+72*dx(1)*dx(5)**2+12*dx(1)**2*dx(5))*dx(3)**4+&
    &(22*dx(4)**4+(140*dx(1)+105*xj+65*dx(5))*dx(4)**3+(50*dx(5)**2+210*dx(5)*xj+420*xj*dx(1)+320*dx(1)*dx(5)+45*dx(1)**2)*dx(4)**2+&
    &(630*dx(1)*dx(5)*xj+105*xj*dx(1)**2+105*dx(5)**2*xj+75*dx(1)**2*dx(5)+180*dx(1)*dx(5)**2)*dx(4)+105*dx(1)**2*dx(5)*xj+&
    &30*dx(1)**2*dx(5)**2+210*dx(1)*dx(5)**2*xj)*dx(3)**3+(2*dx(4)**5+(56*dx(1)+9*dx(5)+42*xj)*dx(4)**4+&
    &(105*dx(5)*xj+40*dx(1)**2+10*dx(5)**2+175*dx(1)*dx(5)+420*xj*dx(1))*dx(4)**3+&
    &(70*dx(5)**2*xj+280*xj*dx(1)**2+840*dx(1)*dx(5)*xj+100*dx(1)**2*dx(5)+140*dx(1)*dx(5)**2)*dx(4)**2+&
    &(60*dx(1)**2*dx(5)**2+420*dx(1)**2*dx(5)*xj+420*dx(1)*dx(5)**2*xj)*dx(4)+140*dx(1)**2*dx(5)**2*xj)*dx(3)**2+&
    &(6*dx(4)**5*dx(1)+(12*dx(1)**2+27*dx(1)*dx(5)+126*xj*dx(1))*dx(4)**4+(45*dx(1)**2*dx(5)+315*dx(1)*dx(5)*xj+&
    &210*xj*dx(1)**2+30*dx(1)*dx(5)**2)*dx(4)**3+(210*dx(1)*dx(5)**2*xj+420*dx(1)**2*dx(5)*xj+&
    &40*dx(1)**2*dx(5)**2)*dx(4)**2+210*dx(1)**2*dx(5)**2*xj*dx(4))*dx(3)+2*dx(4)**5*dx(1)**2+&
    &(42*xj*dx(1)**2+9*dx(1)**2*dx(5))*dx(4)**4+(105*dx(1)**2*dx(5)*xj+10*dx(1)**2*dx(5)**2)*dx(4)**3+&
    &70*dx(4)**2*dx(1)**2*dx(5)**2*xj)*dx(2)+(6*dx(1)*dx(5)+6*dx(1)*dx(4))*dx(3)**6+&
    &(30*dx(1)*dx(4)**2+(11*dx(1)**2+48*dx(1)*dx(5)+21*xj*dx(1))*dx(4)+18*dx(1)*dx(5)**2+21*dx(1)*dx(5)*xj+11*dx(1)**2*dx(5))*dx(3)**5+&
    &(44*dx(1)*dx(4)**3+(84*xj*dx(1)+98*dx(1)*dx(5)+46*dx(1)**2)*dx(4)**2+&
    &(126*dx(1)*dx(5)*xj+54*dx(1)*dx(5)**2+42*xj*dx(1)**2+72*dx(1)**2*dx(5))*dx(4)+42*dx(1)**2*dx(5)*xj+&
    &26*dx(1)**2*dx(5)**2+42*dx(1)*dx(5)**2*xj)*dx(3)**4+&
    &(22*dx(1)*dx(4)**4+(105*xj*dx(1)+55*dx(1)**2+65*dx(1)*dx(5))*dx(4)**3+(120*dx(1)**2*dx(5)+210*dx(1)*dx(5)*xj+50*dx(1)*dx(5)**2+&
    &140*xj*dx(1)**2)*dx(4)**2+(65*dx(1)**2*dx(5)**2+105*dx(1)*dx(5)**2*xj+210*dx(1)**2*dx(5)*xj)*dx(4)+70*dx(1)**2*dx(5)**2*xj)*dx(3)**3+&
    &(2*dx(4)**5*dx(1)+(22*dx(1)**2+9*dx(1)*dx(5)+42*xj*dx(1))*dx(4)**4+(140*xj*dx(1)**2+105*dx(1)*dx(5)*xj+10*dx(1)*dx(5)**2+&
    &65*dx(1)**2*dx(5))*dx(4)**3+(50*dx(1)**2*dx(5)**2+280*dx(1)**2*dx(5)*xj+70*dx(1)*dx(5)**2*xj)*dx(4)**2+&
    &140*dx(1)**2*dx(5)**2*xj*dx(4))*dx(3)**2+(2*dx(4)**5*dx(1)**2+(42*xj*dx(1)**2+9*dx(1)**2*dx(5))*dx(4)**4+&
    &(105*dx(1)**2*dx(5)*xj+10*dx(1)**2*dx(5)**2)*dx(4)**3+70*dx(4)**2*dx(1)**2*dx(5)**2*xj)*dx(3))/(dx(4)+dx(5))/&
    &(((2*dx(4)+dx(5))*dx(3)+dx(4)**2+dx(4)*dx(5))*dx(2)+(2*dx(4)+dx(5))*dx(3)**2+(2*dx(4)**2+2*dx(4)*dx(5))*dx(3))/&
    &((2*dx(3)+dx(4))*dx(2)**2+(2*dx(3)**2+(2*dx(1)+2*dx(4))*dx(3)+dx(1)*dx(4))*dx(2)+dx(1)*dx(3)**2+dx(3)*dx(1)*dx(4))/&
    &(dx(1)+dx(2)))/0.315D3
    
  ENDDO
  WRITE (*,*) '------------------------------------------------------------------'
  write (*,*) 'Assuming p(r)=constant, Bz(r)=constant, Bt(r)=0'
  DO i=1,N-1
    write (*,'(g,a,g)') int_func(psi(i),psi(i),W22A), ' = ', 1./4.*kz**2 * Bz((/1.D0/))**2 / mt**2 *((psi(i)%xj+psi(i)%dx)**4-psi(i)%xj**4)+&
    & 1./2.*(Bz((/1.D0/))**2+gamma*p((/1.D0/)))*((psi(i)%xj+psi(i)%dx)**2-psi(i)%xj**2)
    
    xj = phi(i)%xj; dx1 = phi(i)%dx(1); dx2 = phi(i)%dx(2); Bz1 = minval(Bz((/1.D0/))); ptemp = p((/1.D0/))
    
    write (*,'(g,a,g)') int_func(phi(i),phi(i),W11A,deriv1=.true.,deriv2=.true.)+&
    &2*int_func(phi(i),phi(i),W11B,deriv1=.true.)+int_func(phi(i),phi(i),W11C), ' = ', &
    & 0.1D1/dx1/dx2*(0.12D2*dx2*mt**2*Bz1**2*xj+0.12D2*dx1*mt**2*Bz1**2*xj+0.4D1*kz**2*Bz1**2*dx2*mt**2*xj*dx1**2-&
    &kz**2*Bz1**2*dx2*mt**2*dx1**3+0.4D1*kz**2*Bz1**2*dx2*xj*dx1**2-kz**2*Bz1**2*dx2*dx1**3+&
    &0.12D2*dx2*mt**2*gamma*ptemp*xj+0.4D1*kz**2*Bz1**2*dx1*mt**2*xj*dx2**2+kz**2*Bz1**2*dx1*mt**2*dx2**3+&
    &0.4D1*kz**2*Bz1**2*dx1*xj*dx2**2+kz**2*Bz1**2*dx1*dx2**3+0.12D2*dx1*mt**2*gamma*ptemp*xj)/mt**2/0.12D2
    write (*,'(g,a,g)') int_func(phi(i),phi(i+1),W11A,deriv1=.true.,deriv2=.true.)+int_func(phi(i),phi(i+1),W11B,deriv1=.true.)+int_func(phi(i),phi(i+1),W11B,deriv2=.true.)+&
    & int_func(phi(i),phi(i+1),W11C),  ' = ', &
    & -0.1D1/dx2*(-kz**2*Bz1**2*dx2**3+0.6D1*mt**2*Bz1**2*dx2-0.2D1*kz**2*Bz1**2*xj*dx2**2-0.2D1*kz**2*Bz1**2*mt**2*xj*dx2**2+0.12D2*mt**2*Bz1**2*xj+&
    & 0.12D2*mt**2*gamma*ptemp*xj+0.6D1*mt**2*gamma*ptemp*dx2-kz**2*Bz1**2*mt**2*dx2**3)/mt**2/0.12D2
    write (*,'(g,a,g)') int_func(psi(i),phi(i),W12A,deriv2=.true.) + int_func(psi(i),phi(i),W12B), ' = ', &
    & -(0.6D1*kz**2*Bz1**2*xj**2*dx2+0.4D1*kz**2*Bz1**2*xj*dx2**2+kz**2*Bz1**2*dx2**3+0.12D2*mt**2*Bz1**2*xj+0.6D1*mt**2*Bz1**2*dx2+&
    & 0.12D2*mt**2*gamma*ptemp*xj+0.6D1*mt**2*gamma*ptemp*dx2)/mt**2/0.12D2

  ENDDO
END PROGRAM test_finite_elements
FUNCTION f(r)
  USE local
  REAL(r8) , DIMENSION(:) :: r
  REAL(r8), DIMENSION(size(r)) :: f
  f = 1
END FUNCTION f
