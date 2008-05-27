MODULE vcyl_matrix_module
  USE local
  USE cyl_funcs_module
  USE vcyl_funcs_module
  IMPLICIT NONE
CONTAINS
  FUNCTION A11(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A11
    A11 = rho(r)/r
  END FUNCTION A11
  FUNCTION A22(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A22
    A22 = rho(r)/r
  END FUNCTION A22
  FUNCTION A33(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A33
    A33 = rho(r)/r
  END FUNCTION A33
  FUNCTION B11(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B11
    B11 = rho(r)/r*kV(r)
  END FUNCTION B11
  FUNCTION B12(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B12
    B12 = rho(r)*Vp(r)*Bz(r)/Bmag(r)/r
  END FUNCTION B12
  FUNCTION B13(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B13
    B13 = rho(r)*Vp(r)*Bt(r)/Bmag(r)/r
  END FUNCTION B13
  FUNCTION B22(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B22
    B22 = rho(r)*kV(r)/r
  END FUNCTION B22
  FUNCTION B33(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B33
    B33 = rho(r)*kV(r)/r
  END FUNCTION B33
  FUNCTION B41a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B41a
    B41a = (gamma*P(r)+Bz(r)**2+Bt(r)**2)/r
  END FUNCTION B41a
  FUNCTION B41b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B41b
    B41b = -Bt(r)**2/r**2
  END FUNCTION B41b
  FUNCTION B41c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B41c
    B41c = kdotB(r)**2/r+rho(r)*Vp(r)**2/r
  END FUNCTION B41c
  FUNCTION B42b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42b
    B42b = (gamma*P(r)/Bmag(r)+Bmag(r))*(mt*Bz(r)/r-kz*Bt(r))/r
  END FUNCTION B42b
  FUNCTION B42c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42c
    B42c = 2*Bt(r)*kz*Bmag(r)/r**2+Vp(r)**2*rho(r)*(mt*Bz(r)/r-kz*Bt(r))/Bmag(r)
  END FUNCTION B42c
  FUNCTION B43b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B43b
    B43b = gamma*P(r)*kdotB(r)/Bmag(r)/r
  END FUNCTION B43b
  FUNCTION B43c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B43c
    B43c = rho(r)*Vp(r)**2*kdotB(r)/Bmag(r)
  END FUNCTION B43c
  FUNCTION B52(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B52
    B52 = gamma*P(r)/r**3/Bmag(r)**2*(kz*Bt(r)*r-mt*Bz(r))**2+&
    & rho(r)/r*Vp(r)**2*Bz(r)**2/Bmag(r)**2+Bmag(r)**2*(mt**2+kz**2*r**2)/r**3
  END FUNCTION B52
  FUNCTION B53(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B53
    B53 = (gamma*P(r)*kdotB(r)*(mt*Bz(r)-kz*Bt(r)*r)/r**2+&
    & Bz(r)*Bt(r)*rho(r)*Vp(r)**2/r)/Bmag(r)**2
  END FUNCTION B53
  FUNCTION B63(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B63
    B63 = (gamma*P(r)*kdotB(r)**2/r+rho(r)/r*Vp(r)**2*Bt(r)**2)/Bmag(r)**2
  END FUNCTION B63
  ! These are the surface (boundary) terms
  FUNCTION BCB1(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCB1
    BCB1 = -Bt(r)**2/r**2+rho(r)*Vp(r)**2
  END FUNCTION BCB1
  FUNCTION BCB1nw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCB1nw
    BCB1nw = -Ka*kdotB(r)**2/(r*kz*Kadot)
  END FUNCTION BCB1nw
  FUNCTION BCB1cw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCB1cw
    BCB1cw = -kdotB(r)**2*(Kbdot*La-Lbdot*Ka)/(r*kz*(Ladot*Kbdot-Kadot*Lbdot))
  END FUNCTION BCB1cw
  FUNCTION BJCP1(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BJCP1
    BJCP1 = -kz*Bz(r)*Bt(r)/mt+Bt(r)**2/r-rho(r)*Vp(r)**2*r
  END FUNCTION BJCP1
  FUNCTION BJCP1a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BJCP1a
    BJCP1a = -gamma*P(r)-Bt(r)**2-Bz(r)**2
  END FUNCTION BJCP1a
  FUNCTION BJCP2(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BJCP2
    BJCP2 = -gamma*P(r)+kz*r*Bz(r)*Bt(r)/mt-Bz(r)**2
  END FUNCTION BJCP2
  FUNCTION BJCP3(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BJCP3
    BJCP3 = -gamma*P(r)-Bt(r)**2+mt*Bz(r)*Bt(r)/(kz*r)
  END FUNCTION BJCP3
  FUNCTION n0(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: n0
    n0 = Ka*(mt**2+kz**2*br**2)*(mt*Bt(r)/r+kz*Bz(r))**2*(Lb*Kbdot-Lbdot*Kb)
  END FUNCTION n0
  FUNCTION n1(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: n1
    n1 = -kz*br*tw*Kbdot*(mt*Bt(r)/r+kz*Bz(r))**2*(Kbdot*La-Ka*Lbdot)
  END FUNCTION n1
  FUNCTION d0(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: d0
    d0 = kz*Kadot*(mt**2+kz**2*br**2)*(Lbdot*Kb-Lb*Kbdot)
  END FUNCTION d0
  FUNCTION d1(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: d1
    d1 = -kz**2*br*tw*Kbdot*(Kadot*Lbdot-Kbdot*Ladot)
  END FUNCTION d1
  FUNCTION BCB1rw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCB1rw
    BCB1rw = n0(r)/d0(r)
  END FUNCTION BCB1rw
  FUNCTION BCA1rw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCA1rw
    BCA1rw = d1(r)/d0(r)*BJCP1(r)-n1(r)/d0(r)
  END FUNCTION BCA1rw
  FUNCTION BCA1arw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCA1arw
    BCA1arw = BJCP1a(r)*d1(r)/d0(r)
  END FUNCTION BCA1arw
  FUNCTION BCA2rw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCA2rw
    BCA2rw = BJCP2(r)*d1(r)/d0(r)
  END FUNCTION BCA2rw
  FUNCTION BCA3rw(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: BCA3rw
    BCA3rw = BJCP3(r)*d1(r)/d0(r)
  END FUNCTION BCA3rw
END MODULE vcyl_matrix_module
