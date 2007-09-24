MODULE cyl_funcs_module
  USE local
  IMPLICIT NONE
  REAL(r8), EXTERNAL :: s17aef, s17aff
  INTEGER :: mt, equilib, ifail, nz  !m_theta, choice of equilibrium configuration, 
  REAL(r8) :: kz, gamma, ar, br, rho0, eps, Bz0, Bt0, s2, P0, P1,lambd !ar is the radius of the plasma, br is the radius of the wall
  REAL(r8), PARAMETER, DIMENSION(2) :: gamma_mt_1 = (/1.841183781,3.054236928 /)
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
      CASE(3:4,6:8)
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
      CASE(5:6)
        Bz = Bz0
      CASE(7)
        Bz = sqrt(Bz0**2-2*p0*exp(-r**2/ar**2))
      CASE(8)
        Bz = 0
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
      CASE(6)
        Bt = sqrt(P0*eps)*r/ar
      CASE(7)
        Bt = 0
      CASE(8)
        Bt = 2*Bt0*r/(r**2+P0**2)
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
      CASE(6)
        Bt_prime = sqrt(P0*eps)/ar
      CASE(7)
        Bt_prime = 0
      CASE(8)
        Bt_prime = 2.*Bt0/(r**2+P0**2)-4.*Bt0*r**2/(r**2+P0**2)**2
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
      CASE(6)
        P = P0*(1-eps*r**2/ar**2)
      CASE(7)
        P = P0*exp(-r**2/ar**2)
      CASE(8)
        P = 2*Bt0**2*P0**2/(r**2+P0**2)**2
    ENDSELECT
    RETURN
  END FUNCTION P
  FUNCTION q(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: q
    IF (Bt0.ne.0) THEN
      q = kz*ar*Bz(r)/Bt0
    ELSE IF (P0.ne.0.and.eps.ne.0) THEN
      q = kz*ar*Bz(r)/sqrt(P0*eps)
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
END MODULE cyl_funcs_module
