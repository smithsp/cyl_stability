MODULE vcyl_funcs_module
  USE local
  USE cyl_funcs_module
  IMPLICIT NONE
  REAL(r8) :: Vz0, epsVz, Vp0, epsVp
CONTAINS
  ELEMENTAL FUNCTION Vz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Vz
    SELECT CASE (equilib)
      CASE(1:9)
        Vz = Vz0*(1-epsVz*r**2/ar**2)        
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
        Vp = sqrt(((P(r+epsr)-P(r))/epsr + Bmag(r)*(Bmag(r+epsr)-Bmag(r))/epsr+Bt(r)**2/r)/(rho(r)*r)) 
      CASE(3:9)
        Vp = Vp0*(1-epsVp*r**2/ar**2)
    ENDSELECT
    RETURN
  END FUNCTION Vp
  FUNCTION equilibrium_V(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: equilibrium_V
    equilibrium_V = diff(P,r) + Bz(r)*diff(Bz,r)+Bt(r)*diff(Bt,r)+Bt(r)**2/r-rho(r)*r*Vp(r)**2
  END FUNCTION equilibrium_V
END MODULE vcyl_funcs_module
