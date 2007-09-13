MODULE cyl_matrix_module
  USE local
  USE cyl_funcs_module
  IMPLICIT NONE
CONTAINS
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
END MODULE cyl_matrix_module
