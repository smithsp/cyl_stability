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
  FUNCTION A12(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A12
    A12 = rho(r)/r*mt
  END FUNCTION A12
  FUNCTION A22a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A22a
    A22a = r*rho(r)*Bmag(r)**2/Bz(r)**2
  END FUNCTION A22a
  FUNCTION A22c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A22c
    A22c = rho(r)/r*mt**2
  END FUNCTION A22c
  FUNCTION A33(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: A33
    A33 = rho(r)*r*Bmag(r)**2
  END FUNCTION A33
  ! Note that B11, B12c, B22a, B22c, and B33 are shifted by kVa which is reshifted after solving for omega in vcyl.f95
  FUNCTION B11(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B11
    B11 = rho(r)/r*(kV(r)-kVa)
  END FUNCTION B11
  FUNCTION B12b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B12b
    B12b = -rho(r)*Vp(r)
  END FUNCTION B12b
  FUNCTION B12c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B12c
    B12c = rho(r)/r*mt*(kV(r)-kVa)
  END FUNCTION B12c
  FUNCTION B13(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B13
    B13 = -rho(r)*Vp(r)*Bt(r)
  END FUNCTION B13
  FUNCTION B22a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B22a
    B22a = rho(r)*r*(kV(r)-kVa)*Bmag(r)**2/Bz(r)**2
  END FUNCTION B22a
  FUNCTION B22b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B22b
    B22b = -rho(r)*Vp(r)*mt
  END FUNCTION B22b
  FUNCTION B22c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B22c
    B22c = rho(r)/r*(kV(r)-kVa)*mt**2
  END FUNCTION B22c
  FUNCTION B23(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B23
    B23 = -rho(r)*mt*Vp(r)*Bt(r)
  END FUNCTION B23
  FUNCTION B33(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B33
    B33 = rho(r)*r*(kV(r)-kVa)*Bmag(r)**2
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
  FUNCTION B42a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42a
    B42a = gamma*P(r)*kz*Bt(r)/Bz(r)+Bmag(r)**2*kz*Bt(r)/Bz(r)
  END FUNCTION B42a
  FUNCTION B42b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42b
    B42b = -Bt(r)**2*mt/r**2-2*Bt(r)**3*kz/(Bz(r)*r)-2*Bt(r)*kz*Bz(r)/r+&
    & (r*kz*Bt(r)/Bz(r)-mt)*rho(r)*Vp(r)**2
  END FUNCTION B42b
  FUNCTION B42bb(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42bb
    B42bb = -Bt(r)**2*mt/r**2
  END FUNCTION B42bb
  FUNCTION B42c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B42c
    B42c = mt*kdotB(r)**2/r+Vp(r)**2*rho(r)*mt/r
  END FUNCTION B42c
  FUNCTION B43b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B43b
    B43b = -gamma*P(r)*kdotB(r)
  END FUNCTION B43b
  FUNCTION B43c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B43c
    B43c = -rho(r)*Vp(r)**2*kdotB(r)*r
  END FUNCTION B43c
  FUNCTION B44i(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B44i
    B44i = 1./r
  END FUNCTION B44i
  FUNCTION B45i(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B45i
    B45i = kz*Bt(r)/Bz(r)
  END FUNCTION B45i
  FUNCTION B52a(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B52a
    B52a = gamma*P(r)*r*kz**2*Bt(r)**2/Bz(r)**2+rho(r)*r*Vp(r)**2+&
    & Bt(r)**2*Bmag(r)**2*(mt**2+kz**2*r**2)/r/Bz(r)**2+&
    & 2*mt*Bt(r)*kz*Bmag(r)**2/Bz(r)+r*kz**2*Bmag(r)**2
  END FUNCTION B52a
  FUNCTION B52b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B52b
    B52b = (r*kz*Bt(r)/Bz(r)-mt)*rho(r)*Vp(r)**2*mt-&
    & 2*Bt(r)*mt*kz*Bmag(r)**2/r/Bz(r)-Bt(r)**2*mt**2/r**2
  END FUNCTION B52b
  FUNCTION B52c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B52c
    B52c = mt*B42c(r)
  END FUNCTION B52c
  FUNCTION B53b(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B53b
    B53b = -gamma*P(r)*kdotB(r)*r*kz*Bt(r)/Bz(r)+Bt(r)*rho(r)*Vp(r)**2*r
  END FUNCTION B53b
  FUNCTION B53c(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B53c
    B53c = mt*B43c(r)
  END FUNCTION B53c
  FUNCTION B55i(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B55i
    B55i = r*(kz*Bt(r)/Bz(r))**2
  END FUNCTION B55i
  FUNCTION B63(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B63
    B63 = gamma*P(r)*kdotB(r)**2*r+r*rho(r)*Vp(r)**2*Bt(r)**2
  END FUNCTION B63
  FUNCTION B64i(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B64i
    B64i = -(mt*Bt(r)/r+kz*Bz(r))/r*0.
  END FUNCTION B64i
  FUNCTION B65i(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B65i
    B65i = -(mt*Bt(r)/r+kz*Bz(r))*kz*Bz(r)/Bt(r)*0.
  END FUNCTION B65i
  FUNCTION B66i(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B66i
    B66i = -rho(r)*Bmag(r)**2*&
    & abs((mt+kz*Bz(r)*r/Bt(r))*sqrt(p(r)/rho(r)))
  END FUNCTION B66i
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
