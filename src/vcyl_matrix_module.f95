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
    B11 = (rho(r))*(mt*Vp(r)*(1-1/mt**2)+kz*Vz(r)*(1+1./mt**2))
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
    B22 = rho(r)*r**2/mt**2*(kz*Vz(r)+mt*Vp(r))
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
    B33 = rho(r)/kz**2*(kz*Vz(r)+mt*Vp(r))
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
    B41c = rho(r)*Vp(r)**2*(1./mt**2-1) + Bz(r)**2*kz**2*(1+1./mt**2) - &
    & 2*Bt(r)*Bz(r)/(mt*r)*kz*(1-mt**2)+(mt**2-1)*Bt(r)**2/r**2  !Units?
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
    B42c = rho(r)*Vp(r)**2*r*(1-1./mt**2)-Bz(r)**2*kz**2/mt**2*r+Bt(r)*Bz(r)*kz/mt !Units?
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
    B43c = rho(r)*Vp(r)**2*r+Bt(r)*Bz(r)*kz/mt-Bt(r)**2/r !Units?
  END FUNCTION B43c
  FUNCTION B52(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: B52
    B52 = gamma*P(r)+rho(r)*r**2*Vp(r)**2/mt**2+Bz(r)**2*(1+kz**2*r**2/mt**2) !Units?
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
