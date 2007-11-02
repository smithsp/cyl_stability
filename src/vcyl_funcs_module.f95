MODULE vcyl_funcs_module
  USE local
  USE cyl_funcs_module
  IMPLICIT NONE
CONTAINS
  ELEMENTAL FUNCTION Vz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Vz
    SELECT CASE (equilib)
      CASE(1:10)
        Vz = Vz0*(1-epsVz*r**2/ar**2)
      CASE(11)
        Vz = 0
    ENDSELECT
  END FUNCTION Vz
  FUNCTION Vp(r)
    ! Note that this corresponds to Omega in the formulation
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Vp
    REAL(r8) :: epsr
    epsr = epsilon(r)
    SELECT CASE (equilib)
      CASE(1:2)
        ! In making the following change to ensure equilibrium, there was a noticable slowdown
        Vp = 0! sqrt(((P(r+epsr)-P(r))/epsr + Bmag(r)*(Bmag(r+epsr)-Bmag(r))/epsr+Bt(r)**2/r)/(rho(r)*r)) 
      CASE(3:10)
        Vp = 0! Vp0*(1-epsVp*r**2/ar**2)
      CASE(11)
        Vp = Vp0+epsVp*r+eps*r**2
    ENDSELECT
    RETURN
  END FUNCTION Vp
  FUNCTION equilibrium_V(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: equilibrium_V
    equilibrium_V = P_prime(r) + Bz(r)*Bz_prime(r)+Bt(r)*Bt_prime(r)+Bt(r)**2/r-rho(r)*r*Vp(r)**2
  END FUNCTION equilibrium_V
END MODULE vcyl_funcs_module
