MODULE vcyl_funcs_module
  USE local
  USE cyl_funcs_module
  IMPLICIT NONE
  REAL(r8) :: Kbdot, Lbdot, Kb, Lb, Kadot, Ladot, Ka, La, tw
CONTAINS
  ELEMENTAL FUNCTION Vz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN) :: r
    REAL(r8) :: Vz
    SELECT CASE (equilib)
      CASE(1:10)
        Vz = Vz0*(1-epsVz*r/pi/ar**2)
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
        Vp = 0
      CASE(3:10)
        Vp = 0
      CASE(11)
        Vp = Vp0+epsVp*sqrt(r/pi)+eps*r/pi
    ENDSELECT
    RETURN
  END FUNCTION Vp
  FUNCTION kV(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: kV
    kV = mt*Vp(r)+kz*Vz(r)
  END FUNCTION kV
  FUNCTION equilibrium_V(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: equilibrium_V
    equilibrium_V = P_prime(r) + Bz(r)*Bz_prime(r)+Bt(r)*Bt_prime(r)+2* Bt(r)**2/r-rho(r)*Vp(r)**2/2./pi
  END FUNCTION equilibrium_V
  SUBROUTINE init_bc
    COMPLEX(r8), DIMENSION(1) :: res
    INTEGER :: nz, ifail
    ifail = 0
    CALL s18dcf(real(abs(mt)  ,r8),cmplx(kz*ar,0,r8),1,'U',res,nz,ifail)
    Ka = REAL(res(1))
    CALL s18dcf(real(abs(mt)+1,r8),cmplx(kz*ar,0,r8),1,'U',res,nz,ifail)
    Kadot = -REAL(res(1))+abs(mt)/kz/ar*Ka
    CALL s18def(real(abs(mt)  ,r8),cmplx(kz*ar,0,r8),1,'U',res,nz,ifail)
    La = REAL(res(1))
    CALL s18def(real(abs(mt)+1,r8),cmplx(kz*ar,0,r8),1,'U',res,nz,ifail)
    Ladot = REAL(res(1))+abs(mt)/kz/ar*La
    IF (br.le.1000.0) THEN
      CALL s18dcf(real(abs(mt)  ,r8),cmplx(kz*br,0,r8),1,'U',res,nz,ifail)
      Kb = REAL(res(1))
      CALL s18dcf(real(abs(mt)+1,r8),cmplx(kz*br,0,r8),1,'U',res,nz,ifail)
      Kbdot = -REAL(res(1))+abs(mt)/kz/br*Kb
      CALL s18def(real(abs(mt)  ,r8),cmplx(kz*br,0,r8),1,'U',res,nz,ifail)
      Lb = REAL(res(1))
      CALL s18def(real(abs(mt)+1,r8),cmplx(kz*br,0,r8),1,'U',res,nz,ifail)
      Lbdot = REAL(res(1))+abs(mt)/kz/br*Lb
    ELSE
      Kb = 0; Kbdot = 0; Lb = 0; Lbdot = 0;
    ENDIF
  END SUBROUTINE init_bc
END MODULE vcyl_funcs_module
