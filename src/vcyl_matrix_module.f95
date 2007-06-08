MODULE vcyl_matrix_module
  USE local
  IMPLICIT NONE
  INTEGER :: mt  !m_theta
  REAL(r8) :: kz, gamma, ar, br, rho0, eps, Bz0, Bt0, Vz0, epsVz, Vp0, epsVp !ar is the radius of the plasma, br is the radius of the wall
  LOGICAL :: homo
CONTAINS
  FUNCTION alfven_range(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(2) :: alfven_range
    alfven_range(1) = minval(1/(rho(r))* kz**2*(Bz(r))**2 )
    alfven_range(2) = maxval(1/(rho(r))* kz**2*(Bz(r))**2 )
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
  FUNCTION rho(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: rho
    IF(homo) THEN
      rho = P(r)*rho0 !for Appert Homogeneous Plasma
    ELSE
      rho = rho0 * (1-eps*(r**2/ar**2))
    ENDIF
    RETURN
  END FUNCTION rho
  ELEMENTAL FUNCTION Bz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Bz
    Bz = Bz0
    RETURN
  END FUNCTION Bz
  ELEMENTAL FUNCTION Bt(r) !for B_theta
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Bt
    Bt = Bt0*r/ar
    RETURN
  END FUNCTION Bt
  ELEMENTAL FUNCTION Bt_prime(r) !for B_theta
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Bt_prime
    Bt_prime = Bt0
    RETURN
  END FUNCTION Bt_prime
  ELEMENTAL FUNCTION Bmag(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Bmag
    Bmag = sqrt(Bt(r)**2+Bz(r)**2)
  END FUNCTION Bmag
  FUNCTION P(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: P
    IF(homo) THEN
      P = 1./12./gamma*Bz0**2
    ELSE
      P = rho(r) ! This can't be true or it would violate the equilibrium condition p'=0
      P = 1./12./gamma*Bz0**2
    ENDIF
    RETURN
  END FUNCTION P
  ELEMENTAL FUNCTION Vz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Vz
    Vz = Vz0*(1-epsVz*r**2/ar**2)
  END FUNCTION Vz
  FUNCTION Vp(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Vp
    REAL(r8) :: epsr
    epsr = epsilon(r)
    ! In making the following change to ensure equilibrium, there was a noticable slowdown
    Vp = sqrt(((P(r+epsr)-P(r))/epsr + Bmag(r)*(Bmag(r+epsr)-Bmag(r))/epsr+Bt(r)**2/r)/(rho(r)*r))
    !Vp0*(1-epsVp*r**2/ar**2)
  END FUNCTION Vp
  FUNCTION equilibrium(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: equilibrium
    REAL(r8) :: epsr
    epsr = epsilon(r)
    equilibrium = (P(r+epsr)-P(r))/epsr + Bmag(r)*(Bmag(r+epsr)-Bmag(r))/epsr+Bt(r)**2/r-rho(r)*r*Vp(r)**2
  END FUNCTION equilibrium
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
