MODULE finite_elements_module
  USE local
  IMPLICIT NONE
  REAL(r8) :: epsilo
  REAL(r8) :: it, it1, it2
  
  TYPE linear
    REAL(r8), DIMENSION(2) :: dx
    REAL(r8) :: xj
    INTEGER, DIMENSION(2) :: extent
  END TYPE linear
  
  TYPE constant
    REAL(r8) :: dx, xj
    INTEGER, DIMENSION(2) :: extent
  END TYPE constant
  
  TYPE bspline !Cubic b-spline
    REAL(r8) :: xj
    REAL(r8), DIMENSION(0:3) :: A, B, C, D
    REAL(r8), DIMENSION(4) :: dx
    INTEGER, DIMENSION(2) :: extent
    LOGICAL :: error, deriv
  END TYPE bspline
  
  INTERFACE int_func
    MODULE PROCEDURE int_constant_constant_func, int_constant_linear_func, int_linear_linear_func, int_spline_spline_func
  END INTERFACE int_func
  INTERFACE val
    MODULE PROCEDURE linear_val, constant_val, bspline_val
  END INTERFACE val
  INTERFACE val_prime
    MODULE PROCEDURE linear_val_prime, bspline_val_prime
  END INTERFACE val_prime
  INTERFACE init
    MODULE PROCEDURE linear_init, constant_init, bspline_init
  END INTERFACE init
  
CONTAINS
  
  ELEMENTAL SUBROUTINE linear_init(this,xj,p1,p2,p3,p4,RendZero,deriv)
  ! This function returns a linear element.  p2(optional), p3(optional) : p2<xj<p3
    IMPLICIT NONE
    TYPE(linear), INTENT(OUT) :: this
    REAL(r8), INTENT(IN) :: xj
    REAL(r8), OPTIONAL, INTENT(IN) :: p1, p2, p3, p4
    LOGICAL, INTENT(IN), OPTIONAL :: RendZero, deriv
    this%xj = xj
    this%extent(:) = (/0, 0/)
    this%dx = 0
    IF (present(p2)) THEN
      this%extent(1) = 1
      this%dx(1) = xj-p2
    ENDIF
    IF (present(p3)) THEN
      this%extent(2) = 1
      this%dx(2) = p3-xj
    ENDIF
  END SUBROUTINE linear_init
  
  ELEMENTAL FUNCTION linear_val(this,x)
    IMPLICIT NONE
    TYPE(linear), INTENT(IN) :: this
    REAL(r8) :: linear_val
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: x1
    !
    IF ( this%extent(2).eq.0 ) THEN
      linear_val = 0.
      RETURN
    ENDIF
    x1 = x-this%xj
    IF ( (x1>-this%dx(1)).and.(x1<0) ) THEN
      linear_val = 1.+x1/this%dx(1)
      RETURN
    ELSEIF( (x1>0).and.(x1<this%dx(2)) ) THEN
      linear_val = 1.-x1/this%dx(2)
      RETURN
    ELSEIF ( (x1==0) ) THEN
      linear_val = 1.
      RETURN
    ELSE      
      linear_val = 0.
      RETURN
    ENDIF
    RETURN
  END FUNCTION linear_val
  
  ELEMENTAL FUNCTION linear_val_prime(this,x)
    IMPLICIT NONE
    TYPE(linear), INTENT(IN) :: this
    REAL(r8) :: linear_val_prime
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: x1
    !
    IF ( this%extent(2).eq.0 ) THEN
      linear_val_prime = 0.
      RETURN
    ENDIF
    x1 = x-this%xj
    IF ( (x1>=-this%dx(1)).and.(x1<0) ) THEN
      linear_val_prime = 1./this%dx(1)
      RETURN
    ELSEIF( (x1>0).and.(x1<=this%dx(2)) ) THEN
      linear_val_prime = -1./this%dx(2)
      RETURN
    ELSEIF ( (x1==0) ) THEN
      linear_val_prime = 0 !or should this be 1/2*(1/d1-1/d2)
      RETURN
    ELSE      
      linear_val_prime = 0.
      RETURN
    ENDIF
    RETURN
  END FUNCTION linear_val_prime
  
  ELEMENTAL SUBROUTINE constant_init(this,xj,p1,p2,p3,p4,RendZero,deriv)
    IMPLICIT NONE
    TYPE(constant), INTENT(OUT) :: this
    REAL(r8), INTENT(IN) :: xj
    REAL(r8), INTENT(IN), OPTIONAL :: p1,p2,p3,p4
    LOGICAL, INTENT(IN), OPTIONAL :: RendZero, deriv
    !
    this%xj = xj
    this%dx = p3-xj
    this%extent = (/0,1/)
  END SUBROUTINE constant_init
  
  ELEMENTAL FUNCTION constant_val(this,x) 
    IMPLICIT NONE
    TYPE(constant), INTENT(IN) :: this
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: constant_val
    REAL(r8) :: x1
    !
    x1 = x-this%xj
    IF ( (x1>=0).and.(x1<=this%dx) ) THEN
      constant_val = 1.
    ELSE
      constant_val = 0.
    ENDIF
  END FUNCTION constant_val

  ELEMENTAL SUBROUTINE bspline_init(this,xj1,p1,p2,p3,p4,RendZero,deriv,LendZero)
    IMPLICIT NONE
    TYPE(bspline), INTENT(OUT) :: this
    REAL(r8), INTENT(IN) :: xj1
    REAL(r8), INTENT(IN), OPTIONAL :: p1, p2, p3, p4
    REAL(r8), DIMENSION(4) :: dx, A(0:3), B(0:3), C(0:3), D(0:3)
    REAL(r8) :: denom, numA, numD, denom1, num1, num2, num3, num4, diff2, xj
    LOGICAL, INTENT(IN), OPTIONAL :: RendZero, deriv, LendZero
    LOGICAL :: Rend0, Lend0
    Rend0 = .false.
    Lend0 = .false.
    this%deriv = .false.
    IF(present(RendZero)) Rend0 = RendZero
    IF(present(LendZero)) Lend0 = LendZero
    IF(present(deriv)) this%deriv = deriv
    this%xj = xj1
    xj = xj1
    this%dx = 0.
    this%extent = 0
    this%error = .false.
    IF(present(p2)) THEN
      this%dx(2) = xj-p2
      this%extent(1) = 1
    ENDIF
    IF(present(p1)) THEN
      IF(present(p2)) THEN
        this%extent(1) = 2
        this%dx(1) = p2-p1
      ENDIF
    ENDIF
    IF(present(p3)) THEN
      this%extent(2) = 1
      this%dx(3) = p3-xj
    ENDIF
    IF(present(p4)) THEN
      IF(present(p3)) THEN
        this%extent(2) = 2
        this%dx(4) = p4-p3
      ENDIF
    ENDIF
    IF(minval(this%dx)<0) THEN
      this%error = .true.
      !WRITE(*,*) 'ERROR: The arguments to this must follow the ordering p2<xj<p3<p4 and (p1<p2 or p4<p1).  Exiting!'
      !STOP
    ENDIF
    dx = this%dx
    xj = 0
    A = 0.; B = 0.; C = 0.; D = 0.
    IF((.not.Lend0).and.(.not.Rend0)) THEN
      IF(this%extent(1).eq.1) THEN !This is a left end element that is not zero at the boundary
        dx(1) = dx(2)
      ENDIF
      IF(this%extent(1).eq.0) THEN !This is a left end element that is not zero at the boundary
        dx(1) = dx(3)
        dx(2) = dx(3)
      ENDIF
      IF(this%extent(2).eq.1) THEN !This is a right end element that is not zero at the boundary
        dx(4) = dx(3)
      ENDIF
      IF(this%extent(2).eq.0) THEN !This is a right end element that is not zero at the boundary
        dx(3) = dx(2)
        dx(4) = dx(2)
        IF(this%extent(1).eq.1) THEN
          xj = dx(3)
        ENDIF
      ENDIF
      IF(minval(dx)<=0) THEN
        this%error = .true.
        !WRITE(*,*) 'ERROR: The arguments to this must follow the ordering p1<p2<xj<p3<p4.  Exiting!'
        !STOP
      ENDIF
      !These are the coefficients for x^n, n=0..3
      num1 = 2*(dx(3)+dx(4)+dx(2))*(dx(3)+dx(2))/dx(1)/(dx(1)+dx(2))/&
      & (2*dx(2)*dx(3)**2+dx(3)**2*dx(1)+2*dx(2)*dx(3)*dx(4)+2*dx(3)*dx(1)*dx(2)+&
      & 2*dx(3)*dx(2)**2+dx(3)*dx(1)*dx(4)+dx(2)**2*dx(4)+dx(1)*dx(2)*dx(4))
      A(0) = -(xj - dx(1) - dx(2)) ** 3 / 0.3D1*num1
      A(1) = (xj - dx(1) - dx(2)) ** 2*num1
      A(2) = (-xj + dx(1) + dx(2))*num1
      A(3) = 0.1D1 / 0.3D1*num1
      
      num2=2/(dx(1)+dx(2))/dx(2)/(2*dx(2)*dx(3)**2+&
      & dx(3)**2*dx(1)+2*dx(2)*dx(3)*dx(4)+2*dx(3)*dx(1)*dx(2)+&
      & 2*dx(3)*dx(2)**2+dx(3)*dx(1)*dx(4)+dx(2)**2*dx(4)+&
      & dx(1)*dx(2)*dx(4))
      B(0)=num2*((xj+0.2D1/0.3D1*dx(3)+dx(4)/0.3D1)*dx(2)**4+(0.2D1/0.3D1*dx(3)**2+&
      & (0.4D1/0.3D1*dx(1)+0.2D1/0.3D1*dx(4))*dx(3)+0.2D1/0.3D1*dx(4)*dx(1)-&
      & 0.2D1*xj**2+0.2D1*xj*dx(1))*dx(2)**3+((-xj+dx(1))*dx(3)**2+(dx(4)*dx(1)-&
      & xj*dx(4)+0.2D1/0.3D1*dx(1)**2-0.2D1*xj**2)*dx(3)-&
      & 0.3D1*xj**2*dx(1)+xj**3-xj**2*dx(4)+xj*dx(1)**2+dx(4)*dx(1)**2/0.3D1)*dx(2)**2+&
      & ((-xj*dx(1)+dx(1)**2/0.3D1)*dx(3)**2+&
      & (0.4D1/0.3D1*xj**3-xj*dx(4)*dx(1)+dx(4)*dx(1)**2/0.3D1-0.2D1*xj**2*dx(1))*dx(3)+&
      & xj**3*dx(1)-xj**2*dx(1)**2+0.2D1/0.3D1*xj**3*dx(4)-xj**2*dx(4)*dx(1))*dx(2)+&
      & xj**3*dx(3)**2/0.3D1+(0.2D1/0.3D1*xj**3*dx(1)+xj**3*dx(4)/0.3D1)*dx(3)+&
      & xj**3*dx(1)*dx(4)/0.3D1+xj**3*dx(1)**2/0.3D1)
      B(1)=num2*(-xj**2*dx(1)**2-dx(1)**2*dx(2)**2+0.2D1*xj*dx(2)*dx(1)**2+0.6D1*xj*dx(2)**2*dx(1)-&
      & 0.2D1*dx(2)**3*dx(1)+dx(2)*dx(3)**2*dx(1)+0.2D1*xj*dx(2)*dx(4)*dx(1)+&
      & 0.4D1*xj*dx(2)*dx(3)*dx(1)+dx(2)*dx(3)*dx(4)*dx(1)-0.3D1*xj**2*dx(2)*dx(1)-&
      & 0.2D1*xj**2*dx(3)*dx(1)-xj**2*dx(4)*dx(1)+0.2D1*xj*dx(2)**2*dx(4)-&
      & xj**2*dx(3)**2+dx(2)**2*dx(3)**2-0.4D1*xj**2*dx(2)*dx(3)-0.3D1*xj**2*dx(2)**2+&
      & 0.4D1*xj*dx(2)**3-0.2D1*xj**2*dx(2)*dx(4)-xj**2*dx(3)*dx(4)+&
      & dx(2)**2*dx(3)*dx(4)-dx(2)**4+0.4D1*xj*dx(2)**2*dx(3))
      B(2)=num2*(-0.2D1*dx(2)**3-dx(2)**2*dx(4)-0.3D1*dx(1)*dx(2)**2-0.2D1*dx(3)*dx(2)**2+&
      & 0.3D1*xj*dx(2)**2+0.3D1*xj*dx(2)*dx(1)-0.2D1*dx(3)*dx(1)*dx(2)+&
      & 0.2D1*xj*dx(2)*dx(4)-dx(2)*dx(1)**2+0.4D1*xj*dx(2)*dx(3)-dx(1)*dx(2)*dx(4)+&
      & xj*dx(3)*dx(4)+xj*dx(1)**2+xj*dx(3)**2+0.2D1*xj*dx(3)*dx(1)+xj*dx(4)*dx(1))
      B(3)=num2*(-dx(2)**2-dx(1)*dx(2)-0.2D1/0.3D1*dx(2)*dx(4)-0.4D1/0.3D1*dx(3)*dx(2)-&
      & dx(3)*dx(4)/0.3D1-dx(1)**2/0.3D1-dx(3)**2/0.3D1-&
      & 0.2D1/0.3D1*dx(1)*dx(3)-dx(4)*dx(1)/0.3D1)

      num3=2/(dx(3)+dx(4))/dx(3)/(2*dx(2)*dx(3)**2+&
      & dx(3)**2*dx(1)+2*dx(2)*dx(3)*dx(4)+2*dx(3)*dx(1)*&
      & dx(2)+2*dx(3)*dx(2)**2+dx(3)*dx(1)*dx(4)+&
      & dx(2)**2*dx(4)+dx(1)*dx(2)*dx(4))
      C(0)=num3*(-0.2D1/0.3D1*xj**3*dx(1)*dx(3)-0.4D1/0.3D1*xj**3*dx(2)*dx(3)+&
      & xj*dx(2)**2*dx(3)**2-xj**3*dx(1)*dx(4)/0.3D1-0.2D1/0.3D1*xj**3*dx(2)*dx(4)-&
      & xj**3*dx(3)*dx(4)+0.4D1/0.3D1*dx(2)*dx(3)**3*dx(4)+dx(2)**2*dx(3)**2*dx(4)+&
      & xj*dx(1)*dx(2)*dx(3)**2-0.2D1*xj**2*dx(2)*dx(3)**2-&
      & xj**2*dx(1)*dx(3)**2+xj*dx(1)*dx(2)*dx(3)*dx(4)-xj**3*dx(1)*dx(2)/0.3D1-&
      & 0.2D1*xj*dx(3)**3*dx(4)+0.2D1/0.3D1*dx(1)*dx(3)**3*dx(4)-&
      & 0.3D1*dx(3)**2*xj**2*dx(4)+0.2D1/0.3D1*dx(2)*dx(3)**2*dx(4)**2-&
      & xj*dx(3)**2*dx(4)**2+dx(1)*dx(3)**2*dx(4)**2/0.3D1-&
      & dx(3)*xj**2*dx(4)**2+0.2D1/0.3D1*dx(1)*dx(2)*dx(3)**3+&
      & dx(2)**2*dx(3)*dx(4)**2/0.3D1+xj*dx(2)**2*dx(3)*dx(4)-xj**3*dx(3)**2-&
      & xj**3*dx(2)**2/0.3D1+dx(1)*dx(3)**4/0.3D1-xj**3*dx(4)**2/0.3D1-&
      & 0.2D1*xj**2*dx(3)**3-xj*dx(3)**4+0.2D1/0.3D1*dx(2)*dx(3)**4+&
      & 0.2D1/0.3D1*dx(2)**2*dx(3)**3+dx(1)*dx(2)*dx(3)*dx(4)**2/0.3D1-&
      & dx(3)*xj**2*dx(1)*dx(4)+dx(1)*dx(2)*dx(3)**2*dx(4)-0.2D1*xj**2*dx(2)*dx(4)*dx(3))
      C(1)=num3*(0.2D1*xj*dx(3)**2*dx(1)+xj**2*dx(2)*dx(1)-dx(2)*dx(3)**2*dx(1)-&
      & dx(2)*dx(3)*dx(4)*dx(1)+0.2D1*xj**2*dx(3)*dx(1)+&
      & 0.2D1*xj*dx(3)*dx(4)*dx(1)+xj**2*dx(4)*dx(1)+&
      & 0.4D1*xj*dx(2)*dx(4)*dx(3)+0.3D1*xj**2*dx(3)**2+0.4D1*xj*dx(2)*dx(3)**2-&
      & dx(2)**2*dx(3)**2+0.4D1*xj**2*dx(2)*dx(3)+xj**2*dx(2)**2+dx(3)**4+&
      & 0.2D1*xj**2*dx(2)*dx(4)+0.3D1*xj**2*dx(3)*dx(4)-dx(2)**2*dx(3)*dx(4)+&
      & 0.6D1*xj*dx(4)*dx(3)**2+xj**2*dx(4)**2+0.4D1*xj*dx(3)**3+&
      & 0.2D1*xj*dx(3)*dx(4)**2+dx(3)**2*dx(4)**2+0.2D1*dx(3)**3*dx(4))
      C(2)=num3*(-xj*dx(2)**2-0.2D1*dx(2)*dx(3)*dx(4)-0.4D1*xj*dx(2)*dx(3)-&
      & 0.2D1*xj*dx(2)*dx(4)-xj*dx(2)*dx(1)-0.2D1*dx(2)*dx(3)**2-&
      & dx(3)**2*dx(1)-0.3D1*xj*dx(3)*dx(4)-0.2D1*dx(3)**3-&
      & 0.3D1*dx(3)**2*dx(4)-0.3D1*xj*dx(3)**2-dx(3)*dx(1)*dx(4)-&
      & xj*dx(4)**2-0.2D1*xj*dx(3)*dx(1)-dx(3)*dx(4)**2-xj*dx(4)*dx(1))
      C(3)=num3*(dx(1)*dx(2)/0.3D1+0.2D1/0.3D1*dx(1)*dx(3)+&
      & dx(4)*dx(1)/0.3D1+0.4D1/0.3D1*dx(3)*dx(2)+dx(3)*dx(4)+&
      & 0.2D1/0.3D1*dx(2)*dx(4)+dx(2)**2/0.3D1+dx(3)**2+dx(4)**2/0.3D1)
      
      num4=2*(dx(3)+dx(2))*(dx(3)+dx(1)+dx(2))/dx(4)/(dx(3)+dx(4))/&
      & (2*dx(2)*dx(3)**2+dx(3)**2*dx(1)+2*dx(2)*dx(3)*dx(4)+2*dx(3)*dx(1)*dx(2)+&
      & 2*dx(3)*dx(2)**2+dx(3)*dx(1)*dx(4)+dx(2)**2*dx(4)+dx(1)*dx(2)*dx(4))
      D(0) = (xj + dx(3) + dx(4)) ** 3 / 0.3D1*num4
      D(1) = -(xj + dx(3) + dx(4)) ** 2*num4
      D(2) = (xj + dx(3) + dx(4))*num4
      D(3) = -0.1D1 / 0.3D1*num4
      IF((this%extent(1).eq.1).and.(this%extent(2).eq.0)) THEN
        B = A
        A = 0
        !this%xj = this%xj+dx(3)
      ENDIF

    ELSEIF((this%extent(1).eq.1).and.(this%extent(2).eq.2).and.(Lend0)) THEN
      IF(minval(dx(2:4))<=0) THEN
        this%error = .true.
        !WRITE(*,*) 'ERROR: The arguments to this must follow the ordering p1<p2<xj<p3<p4.  Exiting!'
        !STOP
      ENDIF
      ! These are the coefficients for x^n, n=0..3
      A = 0.     
      ! These are for val(spline,xj-dx[2])=0  diff(spline,x,x)@x=xj-dx[2].ne.0
      num2=2/(2*dx(3)*dx(2)+dx(2)*dx(4)+2*dx(3)**2+2*dx(3)*dx(4))/dx(2)**3
      B(0)=dble((xj-dx(2))**2*(3*xj*dx(2)**2+4*xj*dx(2)*dx(3)+&
      & 2*xj*dx(2)*dx(4)+xj*dx(3)*dx(4)+xj*dx(3)**2+2*dx(3)*dx(2)**2+&
      & dx(2)**2*dx(4)+2*dx(2)*dx(3)*dx(4)+2*dx(2)*dx(3)**2)*num2)/0.3D1
      B(1)=-(xj-dx(2))*(3*xj*dx(2)**2+4*xj*dx(2)*dx(3)+2*xj*dx(2)*dx(4)+&
      & xj*dx(3)*dx(4)+xj*dx(3)**2-&
      & dx(2)**3+dx(2)*dx(3)*dx(4)+dx(2)*dx(3)**2)*num2
      B(2)=(-2*dx(2)**3-dx(2)**2*dx(4)+3*xj*dx(2)**2-2*dx(3)*dx(2)**2+2*xj*dx(2)*dx(4)+&
      & 4*xj*dx(2)*dx(3)+xj*dx(3)*dx(4)+xj*dx(3)**2)*num2
      B(3)=(-dble(dx(2)**2)-0.2D1/0.3D1*dble(dx(2))*dble(dx(4))-&
      & 0.4D1/0.3D1*dble(dx(3))*dble(dx(2))-&
      & dble(dx(3)**2)/0.3D1-dble(dx(3)*dx(4))/0.3D1)*dble(num2)
      
      num3 = 2/(dx(3)+dx(4))/dx(3)/(2*dx(3)*dx(2)+dx(2)*dx(4)+2*dx(3)**2+2*dx(3)*dx(4))/dx(2)
      C(0) = (dx(2)**2*dx(3)**2*dx(4)+xj*dx(2)**2*dx(3)*dx(4)+dx(2)**2*dx(3)*dx(4)**2/0.3D1+&
      & xj*dx(2)**2*dx(3)**2-xj**3*dx(2)**2/0.3D1+&
      & 0.2D1/0.3D1*dx(2)**2*dx(3)**3-0.4D1/0.3D1*xj**3*dx(2)*dx(3)-&
      & 0.2D1*xj**2*dx(2)*dx(3)**2+0.2D1/0.3D1*dx(2)*dx(3)**2*dx(4)**2+&
      & 0.2D1/0.3D1*dx(2)*dx(3)**4-0.2D1/0.3D1*xj**3*dx(2)*dx(4)+&
      & 0.4D1/0.3D1*dx(2)*dx(3)**3*dx(4)-0.2D1*xj**2*dx(2)*dx(4)*dx(3)-&
      & 0.2D1*xj**2*dx(3)**3-0.3D1*dx(3)**2*xj**2*dx(4)-dx(3)*xj**2*dx(4)**2-&
      & xj**3*dx(3)**2-xj**3*dx(4)**2/0.3D1-xj**3*dx(3)*dx(4)-&
      & xj*dx(3)**4-0.2D1*xj*dx(3)**3*dx(4)-xj*dx(3)**2*dx(4)**2)*num3
      C(1) = (xj**2*dx(2)**2-dx(2)**2*dx(3)*dx(4)-dx(2)**2*dx(3)**2+&
      & 0.4D1*xj**2*dx(2)*dx(3)+0.2D1*xj**2*dx(2)*dx(4)+&
      & 0.4D1*xj*dx(2)*dx(3)**2+0.4D1*xj*dx(2)*dx(4)*dx(3)+&
      & 0.3D1*xj**2*dx(3)**2+xj**2*dx(4)**2+0.4D1*xj*dx(3)**3+&
      & 0.6D1*xj*dx(4)*dx(3)**2+0.3D1*xj**2*dx(3)*dx(4)+&
      & 0.2D1*xj*dx(3)*dx(4)**2+0.2D1*dx(3)**3*dx(4)+&
      & dx(3)**4+dx(3)**2*dx(4)**2)*num3
      C(2) = (-xj*dx(2)**2-0.2D1*xj*dx(2)*dx(4)-0.2D1*dx(2)*dx(3)*dx(4)-&
      & 0.2D1*dx(2)*dx(3)**2-0.4D1*xj*dx(2)*dx(3)-&
      & xj*dx(4)**2-0.2D1*dx(3)**3-dx(3)*dx(4)**2-0.3D1*xj*dx(3)*dx(4)-&
      & 0.3D1*xj*dx(3)**2-0.3D1*dx(3)**2*dx(4))*num3
      C(3) = (dx(2)**2/0.3D1+0.2D1/0.3D1*dx(2)*dx(4)+0.4D1/0.3D1*dx(3)*dx(2)+&
      & dx(4)**2/0.3D1+dx(3)*dx(4)+dx(3)**2)*num3
      
      num4 = 2*(dx(3)+dx(2))**2/dx(4)/(dx(3)+dx(4))/dx(2)/&
      & (2*dx(3)*dx(2)+dx(2)*dx(4)+2*dx(3)**2+2*dx(3)*dx(4))
      D(0) = (xj+dx(3)+dx(4))**3*num4/0.3D1
      D(1) = -(xj+dx(3)+dx(4))**2*num4
      D(2) = (xj+dx(3)+dx(4))*num4
      D(3) = -num4/0.3D1
      
    ELSEIF((this%extent(1).eq.2).and.(this%extent(2).eq.1).and.Rend0) THEN
      IF(minval(dx(1:3))<=0) THEN
        this%error = .true.
      ENDIF
      ! These are the coefficients for x^n, n=0..3      
           
      !These are for val(spline,x=xj+dx[3])=0, diff(spline,x)@x=xj+dx[3].eq.0, diff(spline,x,x)@x=xj+dx[3].ne.0
      num1 = 2*(dx(3)+dx(2))**2/dx(1)/(dx(2)+dx(1))/dx(3)/&
      &(dx(1)*dx(3)+2*dx(2)*dx(1)+2*dx(2)**2+2*dx(2)*dx(3))
      A(0) = -(xj - dx(2) - dx(1)) ** 3 * num1 / 0.3D1
      A(1) = (xj - dx(2) - dx(1)) ** 2 * num1
      A(2) = (dx(2) + dx(1) - xj) * num1
      A(3) = num1 / 0.3D1
      
      num2=2/(dx(2)+dx(1))/dx(3)/dx(2)/(dx(1)*dx(3)+&
      & 2*dx(2)*dx(1)+2*dx(2)**2+2*dx(2)*dx(3))
      B(0)=(dx(2)**4*xj+0.2D1/0.3D1*dx(2)**4*dx(3)+&
      & 0.4D1/0.3D1*dx(1)*dx(3)*dx(2)**3+0.2D1*dx(2)**3*xj*dx(1)+&
      & 0.2D1/0.3D1*dx(3)**2*dx(2)**3-0.2D1*xj**2*dx(2)**3+&
      & xj*dx(1)**2*dx(2)**2+0.2D1/0.3D1*dx(1)**2*dx(3)*dx(2)**2-&
      & 0.3D1*xj**2*dx(1)*dx(2)**2+xj**3*dx(2)**2-&
      & 0.2D1*dx(2)**2*xj**2*dx(3)+dx(3)**2*dx(2)**2*dx(1)-&
      & dx(2)**2*xj*dx(3)**2+dx(2)*xj**3*dx(1)+&
      & 0.4D1/0.3D1*dx(2)*dx(3)*xj**3-dx(2)*xj**2*dx(1)**2-&
      & 0.2D1*dx(2)*dx(3)*xj**2*dx(1)+dx(3)**2*dx(2)*dx(1)**2/0.3D1-&
      & dx(2)*dx(3)**2*xj*dx(1)+0.2D1/0.3D1*dx(3)*xj**3*dx(1)+&
      & xj**3*dx(1)**2/0.3D1+dx(3)**2*xj**3/0.3D1)*num2
      B(1)=(-dx(2)**4+0.4D1*xj*dx(2)**3-&
      & 0.2D1*dx(2)**3*dx(1)+0.6D1*xj*dx(1)*dx(2)**2-&
      & dx(1)**2*dx(2)**2+0.4D1*dx(2)**2*xj*dx(3)-&
      & 0.3D1*xj**2*dx(2)**2+dx(2)**2*dx(3)**2+&
      & 0.4D1*dx(2)*dx(3)*xj*dx(1)-0.4D1*dx(2)*dx(3)*xj**2+&
      & 0.2D1*dx(2)*xj*dx(1)**2-0.3D1*dx(2)*xj**2*dx(1)+&
      & dx(2)*dx(3)**2*dx(1)-0.2D1*dx(3)*xj**2*dx(1)-&
      & xj**2*dx(1)**2-dx(3)**2*xj**2)*num2
      B(2)=(-0.2D1*dx(2)**3-0.2D1*dx(3)*dx(2)**2+&
      & 0.3D1*xj*dx(2)**2-0.3D1*dx(2)**2*dx(1)-&
      & dx(1)**2*dx(2)+0.3D1*xj*dx(1)*dx(2)-&
      & 0.2D1*dx(3)*dx(2)*dx(1)+0.4D1*dx(3)*xj*dx(2)+&
      & dx(3)**2*xj+0.2D1*dx(3)*xj*dx(1)+xj*dx(1)**2)*num2
      B(3)=(-dx(2)**2-0.4D1/0.3D1*dx(2)*dx(3)-dx(2)*dx(1)-&
      & 0.2D1/0.3D1*dx(1)*dx(3)-dx(1)**2/0.3D1-dx(3)**2/0.3D1)*num2

      num3 = 2/dx(3)**3/(dx(1)*dx(3)+2*dx(2)*dx(1)+2*dx(2)**2+2*dx(2)*dx(3))
      C(0) = -dble((xj+dx(3))**2*(xj*dx(1)*dx(2)+4*dx(3)*xj*dx(2)+&
      & 2*dx(3)*xj*dx(1)+3*dx(3)**2*xj+xj*dx(2)**2-&
      & dx(1)*dx(3)**2-2*dx(2)*dx(3)**2-2*dx(3)*dx(2)*dx(1)-&
      & 2*dx(3)*dx(2)**2)*num3)/0.3D1
      C(1) = (xj+dx(3))*(xj*dx(2)**2+4*dx(3)*xj*dx(2)+&
      & 3*dx(3)**2*xj+xj*dx(1)*dx(2)+2*dx(3)*xj*dx(1)+dx(3)**3-&
      & dx(3)*dx(2)**2-dx(3)*dx(2)*dx(1))*num3
      C(2) = (-xj*dx(2)**2-xj*dx(1)*dx(2)-2*dx(2)*dx(3)**2-&
      & 4*dx(3)*xj*dx(2)-3*dx(3)**2*xj-2*dx(3)**3-&
      & dx(1)*dx(3)**2-2*dx(3)*xj*dx(1))*num3
      C(3) = (dble(dx(2)**2)/0.3D1+0.4D1/0.3D1*dble(dx(2))*dble(dx(3))+&
      & dble(dx(2)*dx(1))/0.3D1+0.2D1/0.3D1*dble(dx(1))*dble(dx(3))+&
      & dble(dx(3)**2))*dble(num3)

      D = 0.   
     
    ELSEIF((this%extent(1).eq.0).and.(this%extent(2).eq.2).and.(Lend0)) THEN 
      ! This is a left hand element
      ! These are the coefficients for x^n, n=0..3  
      A = 0.
      B = 0.
      
      ! These are for val(spline,xj)=0, diff(spline,x)@x=xj.ne.0, diff(spline,x,x)@x=xj.ne.0
      num3 = 2/dx(3)**3/dx(4)**2
      C(0) = -dble((3*dx(3)**4+6*xj*dx(3)**3+6*dx(3)**3*dx(4)+&
      & 9*xj*dx(4)*dx(3)**2+3*xj**2*dx(3)**2+3*dx(4)**2*dx(3)**2+&
      & 3*xj**2*dx(4)*dx(3)+3*xj*dx(4)**2*dx(3)+xj**2*dx(4)**2)*xj*num3)/0.3D1
      C(1) = (dx(3)**4+4*xj*dx(3)**3+2*dx(3)**3*dx(4)+&
      & 6*xj*dx(4)*dx(3)**2+3*xj**2*dx(3)**2+dx(4)**2*dx(3)**2+&
      & 3*xj**2*dx(4)*dx(3)+2*xj*dx(4)**2*dx(3)+xj**2*dx(4)**2)*num3
      C(2) = (-2*dx(3)**3-3*dx(4)*dx(3)**2-3*xj*dx(3)**2-&
      & 3*xj*dx(4)*dx(3)-dx(4)**2*dx(3)-xj*dx(4)**2)*num3
      C(3) = (dble(dx(3)**2)+dble(dx(3)*dx(4))&
      & +dble(dx(4)**2)/0.3D1)*dble(num3)
      
      num4 = 2/dx(4)**3
      D(0) = (xj+dx(3)+dx(4))**3*num4/0.3D1
      D(1) = -(xj+dx(3)+dx(4))**2*num4
      D(2) = (xj+dx(3)+dx(4))*num4
      D(3) = -num4/0.3D1

    ELSEIF((this%extent(1).eq.2).and.(this%extent(2).eq.0).and.Rend0) THEN
      ! This is a right hand element
      ! These are the coefficients for x^n, n=0..3    
      ! These are for val(spline,xj)=0, diff(spline,x)@x=xj.ne.0, diff(spline,x,x)@x=xj.ne.0  
      num1 = 2/dx(1)**3
      A(0) = -(xj-dx(1)-dx(2))**3*num1/0.3D1
      A(1) = (xj-dx(1)-dx(2))**2*num1
      A(2) = (-xj+dx(1)+dx(2))*num1
      A(3) = num1/0.3D1
      
      num2 = 2 / dx(2) ** 3 / dx(1) ** 2
      B(0) = dble((-3*xj*dx(2)*dx(1)**2+xj**2*dx(1)**2+&
      & 3*dx(1)**2*dx(2)**2+6*dx(2)**3*dx(1)+3*xj**2*dx(2)*dx(1)-&
      & 9*xj*dx(2)**2*dx(1)-6*xj*dx(2)**3+&
      & 3*xj**2*dx(2)**2+3*dx(2)**4)*xj*num2)/0.3D1
      B(1) = (2*xj*dx(2)*dx(1)**2-xj**2*dx(1)**2-dx(1)**2*dx(2)**2+&
      & 6*xj*dx(2)**2*dx(1)-2*dx(2)**3*dx(1)-&
      & 3*xj**2*dx(2)*dx(1)+4*xj*dx(2)**3-&
      & 3*xj**2*dx(2)**2-dx(2)**4)*num2
      B(2) = (xj*dx(1)**2-dx(2)*dx(1)**2-&
      & 3*dx(1)*dx(2)**2+3*dx(2)*xj*dx(1)+&
      & 3*xj*dx(2)**2-2*dx(2)**3)*num2
      B(3) = (-dble(dx(1)**2)/0.3D1-&
      & dble(dx(1)*dx(2))-dble(dx(2)**2))*dble(num2)
      
      C = 0.;  D = 0.
    ELSEIF((this%extent(1).eq.1).and.(this%extent(2).eq.0).and.Rend0) THEN
      ! This is a right hand element
      num2 = 2 / dx(2) ** 3
      B(0) = -(xj - dx(2)) ** 3 * num2 / 0.6D1
      B(1) = (xj - dx(2)) ** 2 * num2 / 0.2D1
      B(2) = (-xj / 0.2D1 + dx(2) / 0.2D1) * num2
      B(3) = num2 / 0.6D1

      A = 0.; C = 0.; D = 0.
    ELSEIF((this%extent(1).eq.1).and.(this%extent(2).eq.1))  THEN
      A = 0.
      D = 0.
      IF(Lend0) THEN
        num2 = 2/dx(2)**3
        B(0) = dble((2*xj+dx(2))*(xj-dx(2))**2*num2)/0.3D1
        B(1) = -2*xj*(xj-dx(2))*num2
        B(2) = (-dx(2)+2*xj)*num2
        B(3) = -0.2D1/0.3D1*dble(num2)
        
        
        num3 = 2/dx(3)**3
        C(0) = -dble((2*xj-dx(3))*(xj+dx(3))**2*num3)/0.3D1
        C(1) = 2*xj*(xj+dx(3))*num3
        C(2) = (-dx(3)-2*xj)*num3
        C(3) = 0.2D1/0.3D1*dble(num3)
        IF((this%deriv)) THEN
          B = 0.
          dx(2) = 0
          this%extent(1) = 0
        ENDIF
      ELSEIF(Rend0) THEN
        
        num2 = 0.2D1/0.3D1/dx(2)**2
        B(0) = -xj*(xj-dx(2))**2*num2
        B(1) = (3*xj-dx(2))*(xj-dx(2))*num2
        B(2) = (-3*xj+2*dx(2))*num2
        B(3) = num2
        
        num3 = 0.2D1/0.3D1/dx(3)**2
        C(0) = -xj*(xj+dx(3))**2*num3
        C(1) = (xj+dx(3))*(3*xj+dx(3))*num3
        C(2) = (-3*xj-2*dx(3))*num3
        C(3) = num3
      ENDIF
    ELSEIF((this%extent(1).eq.0).and.(this%extent(2).eq.1).and.Lend0) THEN
      A = 0.
      B = 0.
      D = 0.
      
      ! This is for val(spline,xj=0).ne.0
      num3 = 2 / dx(3) ** 3
      C(0) = (xj + dx(3)) ** 3 * num3 / 0.3D1
      C(1) = -(xj + dx(3)) ** 2 * num3
      C(2) = (xj + dx(3)) * num3
      C(3) = -num3 / 0.3D1

    ENDIF   

     
    this%A = A
    this%B = B
    this%C = C
    this%D = D
    IF(this%extent(1).eq.1) this%dx(1) = 0.
    IF(this%extent(2).eq.1) this%dx(4) = 0.
  END SUBROUTINE bspline_init
  
  ELEMENTAL FUNCTION bspline_val(this,x)
    TYPE(bspline), INTENT(IN) :: this
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: bspline_val, xj
    REAL(r8), DIMENSION(4) :: dx
    REAL(r8), DIMENSION(0:3) :: coeff
    INTEGER :: j
    IF (this%deriv) THEN
      bspline_val = bspline_val_prime(this,x)
      RETURN
    ELSE
      dx = this%dx
      xj = this%xj
      coeff = 0.
      bspline_val = 0.
      IF ((x > xj-dx(1)-dx(2)).and.(x <= xj-dx(2)))   THEN 
        coeff = this%A
      ELSEIF ((x > xj-dx(2)).and.(x < xj))          THEN 
        coeff = this%B
      ELSEIF ((x.eq.xj))                            THEN
        IF ((this%extent(1).gt.0)) THEN
          coeff = this%B
        ELSE
          coeff = this%C
        ENDIF
      ELSEIF ((x >= xj).and.(x < xj+dx(3)))           THEN 
        coeff = this%C
      ELSEIF ((x >= xj+dx(3)).and.(x < xj+dx(4)+dx(3))) THEN 
        coeff = this%D
      ELSE
        coeff = 0.
      ENDIF
      bspline_val = sum((/ (coeff(j)*(x-xj)**j, j=0,3)/))
      ENDIF
    RETURN
  END FUNCTION bspline_val
  
  ELEMENTAL FUNCTION bspline_val_prime(this,x)
    TYPE(bspline), INTENT(IN) :: this
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: bspline_val_prime, xj
    REAL(r8), DIMENSION(4) :: dx
    REAL(r8), DIMENSION(0:3) :: coeff
    INTEGER :: j
    dx = this%dx
    xj = this%xj
    coeff = 0.
    bspline_val_prime = 0.
    IF ((x > xj-dx(1)-dx(2)).and.(x <= xj-dx(2)))   THEN 
      coeff = this%A
    ELSEIF ((x > xj-dx(2)).and.(x < xj))          THEN 
      coeff = this%B
    ELSEIF ((x.eq.xj))                            THEN
      IF ((this%extent(1).gt.0)) THEN
        coeff = this%B
      ELSE
        coeff = this%C
      ENDIF
    ELSEIF ((x > xj).and.(x < xj+dx(3)))           THEN 
      coeff = this%C
    ELSEIF ((x >= xj+dx(3)).and.(x < xj+dx(4)+dx(3))) THEN 
      coeff = this%D
    ELSE
      coeff = 0.
    ENDIF
    bspline_val_prime = sum((/ (j*coeff(j)*(x-xj)**(j-1), j=1,3)/))
    RETURN
  END FUNCTION bspline_val_prime
  
  FUNCTION int_spline_spline_func(spline1,spline2,func,deriv1,deriv2)
    REAL(r8) :: int_spline_spline_func
    TYPE(bspline) :: spline1,spline2
    INTERFACE
      FUNCTION func(r)
        USE local
        REAL(r8), INTENT(IN), DIMENSION(:) :: r
        REAL(r8), DIMENSION(size(r)) :: func
      END FUNCTION func
    END INTERFACE
    REAL(r8) :: a, b
    INTEGER(i4), PARAMETER :: max_div = 16        ! maximum divisions of [a,b]
    INTEGER(i8) :: j, k, min_k, i, i_max, ind, stat
    REAL(r8), DIMENSION(max_div,max_div) ::  int_arr
    REAL(r8) :: error , r_int, error2, inte, h_k, temp_int(4), xj, dx(4)
    REAL(r8), DIMENSION(:), ALLOCATABLE :: tempa
    LOGICAL, OPTIONAL, INTENT(IN) :: deriv1, deriv2
    LOGICAL :: prime1, prime2
    CALL cpu_time(it1)
    prime1 = .false.
    prime2 = .false.
    IF (present(deriv1)) prime1 = deriv1
    IF (present(deriv2)) prime2 = deriv2
    int_spline_spline_func = 0.
    temp_int=0.
    xj = spline1%xj
    dx = spline1%dx
    IF(prime2.or.spline2%deriv) THEN
      xj = spline2%xj
      dx = spline2%dx
    ENDIF
    DO ind=1,4
      IF(ind.eq.1) THEN
        a = maxval((/spline1%xj-spline1%dx(1)-spline1%dx(2),0./))
        b = maxval((/spline1%xj-spline1%dx(2),0./))
      ELSEIF(ind.eq.2) THEN
        a = maxval((/spline1%xj-spline1%dx(2),0./))
        b = maxval((/spline1%xj,0./))
      ELSEIF(ind.eq.3) THEN
        a = maxval((/spline1%xj,0./))
        b = maxval((/spline1%xj+spline1%dx(3),0./))
      ELSEIF(ind.eq.4) THEN
        a = maxval((/spline1%xj+spline1%dx(3),0./))
        b = maxval((/spline1%xj+spline1%dx(4)+spline1%dx(3),0./))
      ENDIF
      IF ((b<=a).or.(a>=spline2%xj+spline2%dx(3)+spline2%dx(4)).or.(b<=spline2%xj-spline2%dx(1)-spline2%dx(2))) THEN 
        temp_int(ind) = 0.
        !RETURN
      ELSE
        IF (a.eq.0) a = a+b*epsilon(b)
        min_k = 1
        inte = 0.
        ! Outer loop yields next composite trapezoidal rule estimate
        DO k = min_k, max_div
          h_k = (b - a)/(2**(k - 1))
          i_max = 2**(k-2)
          IF(k.eq.min_k) THEN
            stat = 0
            ALLOCATE(tempa(1:2),stat=stat)
            IF(stat.gt.0) WRITE (*,*) 'The problem is here'
            tempa = (/a,b/)
            IF ((prime1).and.(.not.prime2)) THEN
              inte = 0.5 * (b-a) *(sum(func(tempa)*bspline_val_prime(spline1,tempa)*bspline_val(spline2,tempa)))
            ELSEIF ((prime2).and.(.not.prime1)) THEN
              inte = 0.5 * (b-a) *(sum(func(tempa)*bspline_val(spline1,tempa)*bspline_val_prime(spline2,tempa)))
            ELSEIF ((prime1).and.(prime2)) THEN
              inte = 0.5 * (b-a) *(sum(func(tempa)*bspline_val_prime(spline1,tempa)*bspline_val_prime(spline2,tempa)))
            ELSE
              inte = 0.5 * (b-a) *(sum(func(tempa)*bspline_val(spline1,tempa)*bspline_val(spline2,tempa)))
            ENDIF
            DEALLOCATE(tempa)
            IF((i_max*2-1)>0) THEN
              !WRITE (*,'(a,i)') 'i_max*2-1=',i_max*2-1
              ALLOCATE(tempa(i_max*2-1))!,stat=stat)
              !IF(stat.gt.0) WRITE (*,'(a,i)') 'i_max=',i_max
              tempa = (/ (a+h_k*i, i=1,max(i_max*2-1,1))/)
            ENDIF
          ELSE
            !WRITE (*,'(a,i)') 'i_max=',i_max
            ALLOCATE(tempa(i_max))!,stat=stat)
            !IF(stat.gt.0) WRITE (*,'(a,i)') 'i_max=',i_max
            tempa = (/(a + (2.*i-1.)*h_k, i=1,i_max)/)
          ENDIF
          IF (k >= 2) THEN
            IF ((prime1).and.(.not.prime2)) THEN
              inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*bspline_val_prime(spline1,tempa)*bspline_val(spline2,tempa)))
            ELSEIF ((prime2).and.(.not.prime1)) THEN
              inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*bspline_val(spline1,tempa)*bspline_val_prime(spline2,tempa)))
            ELSEIF ((prime1).and.(prime2)) THEN            
              IF (size(tempa)==1) tempa=tempa+tempa*epsilon(tempa(1))
              inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*bspline_val_prime(spline1,tempa)*bspline_val_prime(spline2,tempa)))
            ELSE
              inte = 0.5 *(inte + h_k*2.*sum(func(tempa)*bspline_val(spline1,tempa)*bspline_val(spline2,tempa)))
            ENDIF
          ENDIF
          IF (allocated(tempa)) DEALLOCATE(tempa)
          int_arr(k,min_k)  = inte
          ! Inner loop produces entire row of Romberg tableau
          DO j= (min_k + 1),k
            int_arr(k,j) = int_arr(k,j-1) + (int_arr(k,j-1) - int_arr(k-1,j-1))/(4.**(j-min_k)-1.)
          END DO
          IF (k.ge.min_k+2) THEN
            r_int = int_arr(k,k)
            error = r_int - int_arr(k-1,k-1)
            IF(abs(error).le.epsilo*abs(r_int).or. (abs(r_int).lt.epsilo.and.k.gt.6)) THEN
              error2 = r_int - int_arr(k-2,k-2)
              !test the next nearest neighbor in the diagonal
              IF((abs(error2).le.epsilo*abs(r_int)).or. (abs(r_int).lt.epsilo.and.k.gt.6)) THEN
                temp_int(ind) = r_int
                !WRITE (*,*) 'Good:  int_bspline_bspline_func at grid points ', spline1%xj,' ', spline2%xj,&
                !  &' converged well.   prime1 = ', prime1, '  prime2 = ', prime2
                !WRITE (*,'(a,2g)') '(a,b)=',(/a,b/) , ' func(a,b)=', func((/a,b/)), 'spline1(a,b) = ',bspline_val(spline1,(/a,b/))
                !   'spline2(a,b) = ', bspline_val(spline1,(/a,b/)), 'spline1prime(a,b) = ',bspline_val_prime(spline1,(/a,b/)), 
                !   'spline2prime(a,b) = ',bspline_val_prime(spline2,(/a,b/))
                !WRITE (*,*) 'inte(k) = ', (/ (int_arr(k1,min_k), k1=min_k,k) /)
                !WRITE (*,*) 'int_arr(k,k) = ', (/ (int_arr(k1,k1), k1=min_k,k) /)
                !WRITE (*,*) ' '
                EXIT
              END IF
            END IF
          END IF
        END DO
        IF (k.ge.max_div) THEN
          IF (verbose) THEN
            WRITE (*,*) 'Warning:  int_spline_spline_func at grid points ', spline1%xj,' ', spline2%xj,&
              &' did not converge well.   prime1 = ', prime1, '  prime2 = ', prime2
            WRITE (*,'(a20,2g16.7)') '(a,b)=',(/a,b/) , &
              & ' func(a,b)=', func((/a,b/)), &
              & 'spline1(a,b) = ', bspline_val(spline1,(/a,b/)),&
              & 'spline2(a,b) = ', bspline_val(spline1,(/a,b/)), &
              & 'spline1prime(a,b) = ', bspline_val_prime(spline1,(/a,b/)), &
              & 'spline2prime(a,b) = ',bspline_val_prime(spline2,(/a,b/))
            WRITE (*,*) 'inte(k) = ', (/ (int_arr(k,min_k), k=min_k,max_div) /)
            WRITE (*,*) 'int_arr(k,k) = ', (/ (int_arr(k,k), k=min_k,max_div) /)
            WRITE (*,*) ' '
          ENDIF
          temp_int(ind) = r_int
        ENDIF
      ENDIF
    ENDDO
    int_spline_spline_func = sum(temp_int)
    CALL cpu_time(it2)
    it = it+it2-it1
  END FUNCTION int_spline_spline_func
  
  FUNCTION int_constant_constant_func(const1,const2,func)
    IMPLICIT NONE
    REAL(r8) :: int_constant_constant_func
    TYPE(constant), INTENT(IN) :: const1, const2
    INTERFACE
      FUNCTION func(r)
        USE local
        REAL(r8), INTENT(IN), DIMENSION(:) :: r
        REAL(r8), DIMENSION(size(r)) :: func
      END FUNCTION func
    END INTERFACE
    REAL(r8) :: a, b
    INTEGER(i4), PARAMETER :: max_div = 16        ! maximum divisions of [a,b]
    INTEGER(i4) :: j, k, min_k, i, i_max
    REAL(r8), DIMENSION(max_div,max_div) ::  int_arr
    REAL(r8) :: error , r_int, error2, inte, h_k
    REAL(r8), DIMENSION(:), ALLOCATABLE :: tempa
    !
    a = maxval((/const1%xj,const2%xj/))
    b = minval((/const1%xj+const1%dx,const2%xj+const2%dx/))
    a = a + b*epsilon(b)
    b = b - b*epsilon(b)
    IF (b<=a) THEN 
      int_constant_constant_func = 0.
    ELSE
      min_k = 1
      inte = 0.
      r_int = 0.
      ! Outer loop yields next composite trapezoidal rule estimate
      DO k = min_k, max_div
        h_k = (b - a)/(2.**(k - 1.))
        i_max = 2.**(k-2.)
        IF(k.eq.1) THEN
          inte = 0.5 * (b - a) *(sum(func((/a, b/))*constant_val(const1,(/a,b/))*constant_val(const2,(/a,b/))))
        ELSE
          ALLOCATE(tempa(i_max))
          tempa = (/(a + (2.*i-1.)*h_k, i=1,i_max)/)
          inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*constant_val(const1,tempa)*constant_val(const2,tempa)))
          DEALLOCATE(tempa)
        END IF
        int_arr(k,1)  = inte
        ! Inner loop produces entire row of Romberg tableau
        DO j= (min_k + 1),k
          int_arr(k,j) = int_arr(k,j-1) + (int_arr(k,j-1) - int_arr(k-1,j-1))/(4.**(j-1.)-1.)
        END DO
        IF (k.gt.2) THEN
          r_int = int_arr(k,k)
          error = r_int - int_arr(k-1,k-1)
          IF(abs(error).le.epsilo*abs(r_int).or. (abs(r_int).lt.epsilo.and.k.gt.15)) THEN
            error2 = r_int - int_arr(k-2,k-2)
            !test the next nearest neighbor in the diagonal
            IF((abs(error2).le.epsilo*abs(r_int).and.k.gt.4).or. (abs(r_int).lt.epsilo.and.k.gt.15)) THEN
              int_constant_constant_func = r_int
              RETURN
            END IF
          END IF
        END IF
      END DO
      IF(verbose) THEN
        WRITE (*,*) 'Warning:  int_constant_constant_func at grid points ', const1%xj,' ', const2%xj,&
          &' did not converge well (a,b)=',(/a,b/) , ' for f(a,b)*(a,b)=', func((/a,b/))*(/a,b/)
      ENDIF
      int_constant_constant_func = r_int
    ENDIF
  END FUNCTION int_constant_constant_func
  
  FUNCTION int_constant_linear_func(const1,lin1,func,deriv2)
    IMPLICIT NONE
    REAL(r8) :: int_constant_linear_func
    TYPE(constant), INTENT(IN) :: const1
    TYPE(linear), INTENT(IN) :: lin1
    INTERFACE
      FUNCTION func(r)
        USE local
        REAL(r8), INTENT(IN), DIMENSION(:) :: r
        REAL(r8), DIMENSION(size(r)) :: func
      END FUNCTION func
    END INTERFACE
    REAL(r8) :: a, b
    INTEGER(i4), PARAMETER :: max_div = 20        ! maximum divisions of [a,b]
    INTEGER(i4) :: j, k, min_k, i, i_max
    REAL(r8), DIMENSION(max_div,max_div) ::  int_arr
    REAL(r8) :: error , r_int, error2, inte, h_k
    REAL(r8), DIMENSION(:), ALLOCATABLE :: tempa
    LOGICAL, INTENT(IN), OPTIONAL :: deriv2
    LOGICAL :: prime
    !
    a = maxval((/const1%xj,lin1%xj-lin1%extent(1)*lin1%dx(1)/))
    b = minval((/const1%xj+const1%dx,lin1%xj+lin1%extent(2)*lin1%dx(2)/))
    a = a + b*epsilon(b)
    b = b - b*epsilon(b)
    IF (b<=a) THEN 
      int_constant_linear_func = 0.
    ELSE
      r_int = 0
      prime = .false.
      IF (present(deriv2)) prime = deriv2
      min_k = 1
      inte = 0.
      ! Outer loop yields next composite trapezoidal rule estimate
      DO k = min_k, max_div
        h_k = (b - a)/(2.**(k - 1.))
        i_max = 2.**(k-2.)
        IF(k.eq.1) THEN
          IF (prime) THEN
            inte = 0.5 * (b - a) *(sum(func((/a, b/))*constant_val(const1,(/a,b/))*linear_val_prime(lin1,(/a,b/))))
          ELSE
            inte = 0.5 * (b - a) *(sum(func((/a, b/))*constant_val(const1,(/a,b/))*linear_val(lin1,(/a,b/))))
          ENDIF
        ELSE
          ALLOCATE(tempa(i_max))
          tempa = (/(a + (2.*i-1.)*h_k, i=1,i_max)/)
          IF (prime) THEN
            inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*constant_val(const1,tempa)*linear_val_prime(lin1,tempa)) )
          ELSE
            inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*constant_val(const1,tempa)*linear_val(lin1,tempa)) )
          ENDIF
          DEALLOCATE(tempa)
        END IF
        int_arr(k,1)  = inte
        ! Inner loop produces entire row of Romberg tableau
        DO j= (min_k + 1),k
          int_arr(k,j) = int_arr(k,j-1) + (int_arr(k,j-1) - int_arr(k-1,j-1))/(4.**(j-1.)-1.)
        END DO
        IF (k.gt.2) THEN
          r_int = int_arr(k,k)
          error = r_int - int_arr(k-1,k-1)
          IF(abs(error).le.epsilo*abs(r_int).or. (abs(r_int).lt.epsilo.and.k.gt.15)) THEN
            error2 = r_int - int_arr(k-2,k-2)
            !test the next nearest neighbor in the diagonal
            IF((abs(error2).le.epsilo*abs(r_int).and.k.gt.4).or. (abs(r_int).lt.epsilo.and.k.gt.15)) THEN
              int_constant_linear_func = r_int
              RETURN
            END IF
          END IF
        END IF
      END DO
      int_constant_linear_func = r_int
      IF(verbose) THEN
        WRITE (*,*) 'Warning:  int_constant_linear_func at grid points ', const1%xj,' ', lin1%xj,&
          &' did not converge well (a,b)=',(/a,b/) , ' for f(a,b)=', func((/a,b/)), &
          &' prime = ', prime
      ENDIF
    ENDIF
  END FUNCTION int_constant_linear_func
  
  FUNCTION int_linear_linear_func(lin1,lin2,func,deriv1,deriv2)
    IMPLICIT NONE
    REAL(r8) :: int_linear_linear_func
    TYPE(linear), INTENT(IN) :: lin1, lin2
    INTERFACE
      FUNCTION func(r)
        USE local
        REAL(r8), INTENT(IN), DIMENSION(:) :: r
        REAL(r8), DIMENSION(size(r)) :: func
      END FUNCTION func
    END INTERFACE
    REAL(r8) :: a, b
    INTEGER(i4), PARAMETER :: max_div = 16        ! maximum divisions of [a,b]
    INTEGER(i4) :: j, k, min_k, i, i_max, k1, ind
    REAL(r8), DIMENSION(max_div,max_div) ::  int_arr
    REAL(r8) :: error , r_int, error2, inte, h_k, temp_int(2)
    REAL(r8), DIMENSION(:), ALLOCATABLE :: tempa
    LOGICAL, OPTIONAL, INTENT(IN) :: deriv1, deriv2
    LOGICAL :: prime1, prime2
    DO ind=1,2
      temp_int(ind) = 0.
      !IF (lin2%xj == lin1%xj) THEN  
      IF (ind==1) THEN
        a = max(lin2%xj-lin2%extent(1)*lin2%dx(1),lin1%xj-lin1%extent(1)*lin1%dx(1))
        b = lin1%xj
      ELSEIF (ind==2) THEN
        a = lin1%xj
        b = min(lin2%xj+lin2%extent(2)*lin2%dx(2),lin1%xj+lin1%extent(2)*lin1%dx(2))
      ENDIF
      !ELSE
        !a = maxval((/lin2%xj-lin2%extent(1)*lin2%dx(1),lin1%xj-lin1%extent(1)*lin1%dx(1)/))
        !b = minval((/lin2%xj+lin2%extent(2)*lin2%dx(2),lin1%xj+lin1%extent(2)*lin1%dx(2)/))
      !ENDIF
      a = a+b*epsilon(b)
      b = b-b*epsilon(b)
      IF (b<=a) THEN 
        temp_int(ind) = 0.
      ELSE
        prime1 = .false.
        prime2 = .false.
        IF (present(deriv1)) prime1 = deriv1
        IF (present(deriv2)) prime2 = deriv2
        min_k = 1
        inte = 0.
        ! Outer loop yields next composite trapezoidal rule estimate
        DO k = min_k, max_div
          h_k = (b - a)/(2**(k - 1))
          i_max = 2**(k-2)
          IF(k.eq.min_k) THEN
            ALLOCATE(tempa(2))
            tempa = (/a,b/)
            IF ((prime1).and.(.not.prime2)) THEN
              inte = 0.5 * (b-a) *(sum(func(tempa)*linear_val_prime(lin1,tempa)*linear_val(lin2,tempa)))
            ELSEIF ((prime2).and.(.not.prime1)) THEN
              inte = 0.5 * (b-a) *(sum(func(tempa)*linear_val(lin1,tempa)*linear_val_prime(lin2,tempa)))
            ELSEIF ((prime1).and.(prime2)) THEN
              inte = 0.5 * (b-a) *(sum(func(tempa)*linear_val_prime(lin1,tempa)*linear_val_prime(lin2,tempa)))
            ELSE
              inte = 0.5 * (b-a) *(sum(func(tempa)*linear_val(lin1,tempa)*linear_val(lin2,tempa)))
            ENDIF
            DEALLOCATE(tempa)
            ALLOCATE(tempa(i_max*2-1))
            IF (i_max*2-1.gt.0) THEN 
              tempa = (/ (a+h_k*i, i=1,max(i_max*2-1,1))/)
            ENDIF
          ELSE
            ALLOCATE(tempa(i_max))
            tempa = (/(a + (2.*i-1.)*h_k, i=1,i_max)/)
          ENDIF
          !WRITE (*,*) 'tempa=',tempa,'k=',k,'i_max = ', i_max
          IF (k >= 2) THEN
            IF ((prime1).and.(.not.prime2)) THEN
              !IF (size(tempa)==1) tempa=tempa+tempa*epsilon(tempa(1))
              inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*linear_val_prime(lin1,tempa)*linear_val(lin2,tempa)))
            ELSEIF ((prime2).and.(.not.prime1)) THEN
              !IF (size(tempa)==1) tempa=tempa+epsilon(tempa(1))
              inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*linear_val(lin1,tempa)*linear_val_prime(lin2,tempa)))
            ELSEIF ((prime1).and.(prime2)) THEN            
              IF (size(tempa)==1) tempa=tempa+tempa*epsilon(tempa(1))
              inte = 0.5 * (inte + h_k*2.*sum(func(tempa)*linear_val_prime(lin1,tempa)*linear_val_prime(lin2,tempa)))
            ELSE
              inte = 0.5 *(inte + h_k*2.*sum(func(tempa)*linear_val(lin1,tempa)*linear_val(lin2,tempa)*tempa))
            ENDIF
          ENDIF
          DEALLOCATE(tempa)
          int_arr(k,min_k)  = inte
          ! Inner loop produces entire row of Romberg tableau
          DO j= (min_k + 1),k
            int_arr(k,j) = int_arr(k,j-1) + (int_arr(k,j-1) - int_arr(k-1,j-1))/(4.**(j-min_k)-1.)
          END DO
          IF (k.gt.min_k+3) THEN
            r_int = int_arr(k,k)
            error = r_int - int_arr(k-1,k-1)
            IF(abs(error).le.epsilo*abs(r_int).or. (abs(r_int).lt.epsilo.and.k.gt.15)) THEN
              error2 = r_int - int_arr(k-2,k-2)
              !test the next nearest neighbor in the diagonal
              IF((abs(error2).le.epsilo*abs(r_int).and.k.gt.4).or. (abs(r_int).lt.epsilo.and.k.gt.15)) THEN
                temp_int(ind) = r_int
                !WRITE (*,*) 'Good:  int_linear_linear_func at grid points ', lin1%xj,' ', lin2%xj,&
                !  &' converged well.   prime1 = ', prime1, '  prime2 = ', prime2
                !WRITE (*,'(a,2g)') '(a,b)=',(/a,b/) , ' func(a,b)=', func((/a,b/)), 'lin1(a,b) = ',linear_val(lin1,(/a,b/)),&
                !  & 'lin2(a,b) = ', linear_val(lin1,(/a,b/)), 'lin1prime(a,b) = ',linear_val_prime(lin1,(/a,b/)), &
                !  & 'lin2prime(a,b) = ',linear_val_prime(lin2,(/a,b/))
                !WRITE (*,*) 'inte(k) = ', (/ (int_arr(k1,min_k), k1=min_k,k) /)
                !WRITE (*,*) 'int_arr(k,k) = ', (/ (int_arr(k1,k1), k1=min_k,k) /)
                !WRITE (*,*) ' '
                EXIT
              END IF
            END IF
          END IF
        END DO
        IF (k >= max_div) THEN
          IF (verbose) THEN
            WRITE (*,*) 'Warning:  int_linear_linear_func at grid points ', lin1%xj,' ', lin2%xj,&
              &' did not converge well.   prime1 = ', prime1, '  prime2 = ', prime2
            WRITE (*,'(a20,2g16.7)') '(a,b)=',(/a,b/) , ' func(a,b)=', func((/a,b/)), 'lin1(a,b) = ',linear_val(lin1,(/a,b/)),&
              & 'lin2(a,b) = ', linear_val(lin1,(/a,b/)), 'lin1prime(a,b) = ',linear_val_prime(lin1,(/a,b/)), &
              & 'lin2prime(a,b) = ',linear_val_prime(lin2,(/a,b/))
            WRITE (*,*) 'inte(k) = ', (/ (int_arr(k,min_k), k=min_k,max_div) /)
            WRITE (*,*) 'int_arr(k,k) = ', (/ (int_arr(k,k), k=min_k,max_div) /)
            WRITE (*,*) ' '
          ENDIF
          temp_int(ind) = r_int
        ENDIF
      ENDIF
      !WRITE (*,*) 'temp_int(',ind,') = ', temp_int(ind), 'lin1%xj = ', lin1%xj, 'lin2%xj = ', lin2%xj
    ENDDO
    int_linear_linear_func = sum(temp_int)!/real(size(temp_int))
  END FUNCTION int_linear_linear_func
END MODULE finite_elements_module
