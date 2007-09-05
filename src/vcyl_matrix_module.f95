MODULE vcyl_matrix_module
  USE local
  IMPLICIT NONE
  REAL(r8), EXTERNAL :: s17aef, s17aff
  INTEGER :: mt, equilib, ifail, nz  !m_theta, choice of equilibrium configuration, 
  REAL(r8) :: kz, gamma, ar, br, rho0, eps, Bz0, Bt0, Vz0, epsVz, Vp0, epsVp, s2, P0, P1,lambd !ar is the radius of the plasma, br is the radius of the wall
  REAL(r8), PARAMETER, DIMENSION(2) :: gamma_mt_1 = (/1.841183781,3.054236928 /) ! These are the 1st zeros of d J_m/ dx
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
  ELEMENTAL FUNCTION Vz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Vz
    SELECT CASE (equilib)
      CASE(1:5)
        Vz = Vz0*(1-epsVz*r**2/ar**2)        
    ENDSELECT
  END FUNCTION Vz
  FUNCTION Vp(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Vp
    REAL(r8) :: epsr
    epsr = epsilon(r)
    SELECT CASE (equilib)
      CASE(1:2)
        ! In making the following change to ensure equilibrium, there was a noticable slowdown
        Vp = sqrt(((P(r+epsr)-P(r))/epsr + Bmag(r)*(Bmag(r+epsr)-Bmag(r))/epsr+Bt(r)**2/r)/(rho(r)*r)) 
      CASE(3:5)
        Vp = Vp0*(1-epsVp*r**2/ar**2)
    ENDSELECT
    RETURN
  END FUNCTION Vp
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
    equilibrium = diff(P,r) + Bz(r)*diff(Bz,r)+Bt(r)*diff(Bt,r)+Bt(r)**2/r-rho(r)*r*Vp(r)**2
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
  FUNCTION A11(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A11
    A11 = rho(r)*(1+1./mt**2)
  END FUNCTION A11
  FUNCTION A12(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A12
    A12 = -rho(r)*r/mt**2
  END FUNCTION A12
  FUNCTION B11(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B11
    B11 = (rho(r))*(mt*Vp(r)/2.*(1-1/mt**2)+kz*Vz(r)*(1+1./mt**2))
  END FUNCTION B11
  FUNCTION B12(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B12
    B12 = -rho(r)*r/mt**2*kz*Vz(r)
  END FUNCTION B12
  FUNCTION A22(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A22
    A22 = rho(r)*r**2/mt**2
  END FUNCTION A22
  FUNCTION B22(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B22
    B22 = rho(r)*r**2/mt**2*(kz*Vz(r)+mt*Vp(r)/2.)
  END FUNCTION B22
  FUNCTION A33(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A33
    A33 = rho(r)/kz**2
  END FUNCTION A33
  FUNCTION B33(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B33
    B33 = rho(r)/kz**2*(kz*Vz(r)+mt*Vp(r)/2.)
  END FUNCTION B33
  FUNCTION B41a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B41a
    B41a = (gamma*P(r)+Bmag(r)**2)
  END FUNCTION B41a
  FUNCTION B41b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B41b
    B41b = Bt(r)*Bz(r)*kz/mt
  END FUNCTION B41b
  FUNCTION B41c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B41c
    B41c = rho(r)*Vp(r)**2/r*(1./mt**2-1) + Bz(r)**2*kz**2*(1+1./mt**2) - &
    & 2*Bt(r)*Bz(r)/(mt*r)*kz*(1-mt**2)+(mt**2-1)*Bt(r)**2/r**2
  END FUNCTION B41c
  FUNCTION B42b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42b
    B42b = gamma*P(r)+Bz(r)**2-Bt(r)*Bz(r)*r*kz/mt
  END FUNCTION B42b
  FUNCTION B42c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42c
    B42c = rho(r)*Vp(r)**2*(1-1./mt**2)-Bz(r)**2*kz**2/mt**2*r-Bt(r)*Bz(r)*kz/mt
  END FUNCTION B42c
  FUNCTION B43b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B43b
    B43b = gamma*P(r)-Bt(r)*Bz(r)*mt/(r*kz)+Bt(r)**2
  END FUNCTION B43b
  FUNCTION B43c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B43c
    B43c = rho(r)*Vp(r)**2+Bt(r)*Bz(r)*kz/mt-Bt(r)**2/r
  END FUNCTION B43c
  FUNCTION B52(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B52
    B52 = gamma*P(r)+rho(r)*r*Vp(r)**2/mt**2+Bz(r)**2*(1+kz**2*r**2/mt**2)
  END FUNCTION B52
  FUNCTION B53(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B53
    B53 = gamma*P(r)-Bz(r)*Bt(r)/(mt*r*kz)*(mt**2+r**2*kz**2)
  END FUNCTION B53
  FUNCTION B63(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B63
    B63 = gamma*P(r)+(mt**2/(kz**2*r**2)+1)*Bt(r)**2
  END FUNCTION B63
END MODULE vcyl_matrix_module
