MODULE cyl_funcs_module
  USE local
  IMPLICIT NONE
  REAL(r8), EXTERNAL :: s17aef, s17aff
  INTEGER :: mt, equilib, ifail, nz  !m_theta, choice of equilibrium configuration, 
  REAL(r8) :: kz, gamma, ar, br, rho0, eps, Bz0, Bt0, s2, P0, P1,lambd !ar is the radius of the plasma, br is the radius of the wall
  REAL(r8), PARAMETER, DIMENSION(2) :: gamma_mt_1 = (/1.841183781,3.054236928 /)
  LOGICAL :: homo
CONTAINS
  FUNCTION alfven_range(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(2) :: alfven_range
    alfven_range(1) = minval(1/(rho(r))* kz**2*(Bmag(r))**2 )
    alfven_range(2) = maxval(1/(rho(r))* kz**2*(Bmag(r))**2 )
    RETURN    
  END FUNCTION alfven_range
  FUNCTION slow_inf_range(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(2) :: slow_inf_range
    slow_inf_range(1) = minval(gamma*P(r)*kz**2/(rho(r)*(1.+gamma*P(r)/(Bz(r))**2)))
    slow_inf_range(2) = maxval(gamma*P(r)*kz**2/(rho(r)*(1.+gamma*P(r)/(Bz(r))**2)))
    RETURN 
  END FUNCTION slow_inf_range
  FUNCTION max_slow()
    IMPLICIT NONE
    REAL(r8) :: max_slow
    max_slow = Bz0**2*(ar**2*kz**2+gamma_mt_1(mt)**2)*(1+s2)/(2.*maxval(rho((/ar/)))*ar**2)*&
    & (1.-sqrt(1-4*s2*(kz*ar)**2/((1+s2)**2*((ar*kz)**2+gamma_mt_1(mt)**2))))
  END FUNCTION
  FUNCTION rho(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: rho
    SELECT CASE (equilib)
      CASE(1)
        rho = P(r)*rho0 !for Appert Homogeneous Plasma
      CASE(2)
        rho = rho0 * (1-eps*(r**2/ar**2))
      CASE(3:4)
        rho = rho0
      CASE(5)
        rho = P(r)*rho0
    ENDSELECT
    RETURN
  END FUNCTION rho
  FUNCTION Bz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bz
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:3)
        Bz = Bz0
      CASE(4)
        DO i=1,size(r)
          Bz(i) = sqrt(1.-P1)* s17aef(lambd*r(i),ifail)
        ENDDO
      CASE(5)
        Bz = Bz0
    ENDSELECT
    RETURN
  END FUNCTION Bz
  FUNCTION Bt(r) !for B_theta
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bt
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:3)
        Bt = Bt0*r/ar
      CASE(4)
        DO i=1,size(r)
          Bt(i) = s17aff(lambd*r(i),ifail)
        ENDDO
      CASE(5)
        Bt = Bt0*r/ar
    ENDSELECT
    RETURN
  END FUNCTION Bt
  FUNCTION Bt_prime(r) !for B_theta
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bt_prime
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:3)
        Bt_prime = Bt0/ar
      CASE(4)
        WRITE(*,*) 'Cannot do equilib 4 for Bt_prime'
        STOP
      CASE(5)
        Bt_prime = Bt0/ar
    ENDSELECT
    RETURN
  END FUNCTION Bt_prime
  FUNCTION Bmag(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bmag
    Bmag = sqrt(Bt(r)**2+Bz(r)**2)
  END FUNCTION Bmag
  FUNCTION P(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: P
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:2)
        P = s2/gamma*Bz0**2  
      CASE(3)
        P = Bt0**2*(1-r**2/ar**2)
      CASE(4)
        DO i=1,size(r)
          P(i) = P0+P1/2.*s17aef(lambd*r(i),ifail)**2
        ENDDO
      CASE(5)
        P = Bt0**2*(1-r**2/ar**2)
    ENDSELECT
    RETURN
  END FUNCTION P
  
  FUNCTION q(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: q
    IF (Bt0.ne.0) THEN
      q = kz*ar*Bz(r)/Bt0
    ELSE
      q = 0
    ENDIF
    RETURN
  END FUNCTION q
  FUNCTION equilibrium(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: equilibrium
    equilibrium = diff(P,r) + Bz(r)*diff(Bz,r)+Bt(r)*diff(Bt,r)+Bt(r)**2/r
  END FUNCTION equilibrium
  FUNCTION diff(func,r)
    REAL(r8), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: diff
    INTEGER :: i
    REAL(r8) :: epsr
    INTERFACE
      FUNCTION func(r)
        USE local
        REAL(r8), INTENT(IN), DIMENSION(:) :: r
        REAL(r8), DIMENSION(size(r)) :: func
      END FUNCTION func
    END INTERFACE
    DO i=1,size(r)
      epsr = 10*epsilon(r(i))
      diff(i) = maxval((func((/r(i)+epsr/))-func((/r(i)-epsr/)))/(2*epsr))
    ENDDO
  END FUNCTION diff
  FUNCTION K12(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: K12
    K12 = rho(r)*r
  END FUNCTION K12
  FUNCTION K22(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: K22
    K22 = rho(r)*r**2
    RETURN
  END FUNCTION K22
  FUNCTION W11A(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W11A
    W11A = gamma*P(r)+(Bz(r))**2+(Bt(r))**2
    RETURN
  END FUNCTION W11A
  FUNCTION W11B(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W11B, Bt1, Bz1, kr2
    Bt1 = Bt(r)
    Bz1 = Bz(r)
    kr2 = (kz*r)**2
    W11B = (kz*r*Bz(r)-mt*Bt1)*Bt1/(mt*r)!-(Bt1)**2/r + (kz*mt*Bz1*Bt1 + kz*kr2*Bz1*Bt1/mt)/(kr2+mt**2)
    RETURN
  END FUNCTION W11B
  FUNCTION W11C(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W11C, krBz, mBt, Bt1
    Bt1 = Bt(r)
    krBz = kz*r*Bz(r)
    mBt = mt*Bt1
    W11C = ((krBz-mBt)**2/mt**2-2*Bt1**2+(krBz+mBt)**2-2*Bt1*r*Bt_prime(r))/r**2
    RETURN
  END FUNCTION W11C
  FUNCTION W12A(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W12A, Bz1
    Bz1 = Bz(r)
    W12A = gamma*P(r)+Bz1**2-kz*r/mt*Bt(r)*Bz1
    RETURN
  END FUNCTION W12A
  FUNCTION W12B(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W12B, Bz1
    Bz1 = Bz(r)
    W12B = kz*Bz1/mt**2*(mt*Bt(r)-kz*r*Bz1)
    RETURN
  END FUNCTION W12B
  FUNCTION W13A(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W13A, Bt1
    Bt1 = Bt(r)
    W13A = gamma*P(r)-Bt1*mt*Bz(r)/(kz*r)+Bt1**2
    RETURN
  END FUNCTION W13A
  FUNCTION W13B(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W13B, Bt1
    Bt1 = Bt(r)
    W13B = kz*Bz(r)*Bt1/mt - Bt1**2/r
    RETURN
  END FUNCTION W13B
  FUNCTION W22A(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W22A
    W22A = gamma*P(r)+(Bz(r))**2 *(1+(kz*r)**2/mt**2)
    RETURN
  END FUNCTION W22A
  FUNCTION W23A(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W23A, kr
    kr = kz*r
    W23A = gamma*P(r)-(kr**2+mt**2)/(kr*mt)*Bz(r)*Bt(r)
    RETURN
  END FUNCTION W23A
  FUNCTION W33A(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: W33A
    W33A = gamma*P(r) + (Bt(r)/(kz*r))**2 * ((kz*r)**2+mt**2)
    RETURN
  END FUNCTION W33A
END MODULE cyl_funcs_module
