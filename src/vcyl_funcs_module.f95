MODULE vcyl_funcs_module
  USE local
  USE cyl_funcs_module
  IMPLICIT NONE
  REAL(r8) :: Kbdot, Lbdot, Kb, Lb, Kadot, Ladot, Ka, La, tw
  REAL(r8) :: BCB11a, BCB11nw, BCB11cw
  REAL(r8) :: BCA25aa, BCB21aa, BCB21ac, BCB22aa, BCB22ac, BCB23a, BCB24a, BCB25aa, BCB25ac
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
      CASE(12)
        Vz =  Vp0/Bt0*sqrt(ar**2*(Bz0**2-2.*P0+2.*Bt0**2)+2.*r**2*(P0-Bt0**2))
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
        Vp = Vp0+epsVp*r+eps*r**2
      CASE(12)
        Vp = Vp0
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
    equilibrium_V = P_prime(r) + Bz(r)*Bz_prime(r)+Bt(r)*Bt_prime(r)+Bt(r)**2/r-rho(r)*r*Vp(r)**2
  END FUNCTION equilibrium_V
  SUBROUTINE init_bc
    COMPLEX(r8), DIMENSION(1) :: res
    INTEGER :: nz, ifail
    REAL(r8), DIMENSION(1) :: a
    a = ar
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
      Kb = 0; Kbdot = 0; Lb = 0; Lbdot = 0; tw = 0
    ENDIF
    ! The form of the boundary terms of the 4th row is
    ! phi_1i(a)*[(BCB11nw+BCB11a)*xi1(a)+(BCB11nw+BCB11a)*mt*xi2(a)]   or
    ! phi_1i(a)*[(BCB11cw+BCB11a)*xi1(a)+(BCB11cw+BCB11a)*mt*xi2(a)]
    BCB11a = minval(-Bt(a)**2/a**2+rho(a)*Vp(a)**2)
    BCB11nw = minval(-Ka*(mt*Bt(a)/a+kz*Bz(a))**2/(a*kz*Kadot))
    BCB11cw = minval(-(mt*Bt(a)/a+kz*Bz(a))**2*(Kbdot*La-Lbdot*Ka)/&
      & (a*kz*(Ladot*Kbdot-Kadot*Lbdot)))
    ! The form of the boundary terms of the 5th row is
    ! phi_2i(a)*{mt*[(BCB11nw+BCB11a)*xi1(a)+(BCB11nw+BCB11a)*mt*xi2(a)]+BCB21ac*xi1(a)+BCB21aa*xi1'(a)+BCB22ac*xi2(a)+BCB22aa*xi2'(a)+BCB23a*xi3(a)+
    !            BCB24a*u1(a)+BCB25aa*u2(a)+BCB25ac*u2'(a)} or
    ! phi_2i(a)*{mt*[(BCB11cw+BCB11a)*xi1(a)+(BCB11cw+BCB11a)*mt*xi2(a)]+BCB21ac*xi1(a)+BCB21aa*xi1'(a)+BCB22ac*xi2(a)+BCB22aa*xi2'(a)+BCB23a*xi3(a)+
    !            BCB24a*u1(a)+BCB25aa*u2(a)+BCB25ac*u2'(a)} 
    ! The form of the boundary terms of the 2nd row is
    ! phi_5i(a)*{BCB24a*xi1(a)+BCB25ac*xi2(a)+BCB25aa*xi2'(a)+BCA25a*u2'(a)}
    BCA25aa = minval(-rho(a)*Bmag(a)**2/Bz(a)**2*a)
    !This is for xi1
    BCB21ac = minval((mt/a-kz*Bt(a)/Bz(a))*(rho(a)*Vp(a)**2)+&
      & 2*kz*Bt(a)/Bz(a)*Bmag(a)**2/a**2)*ar*0.
    !This is for xi1_prime
    BCB21aa = minval((mt/a-kz*Bt(a)/Bz(a))*(Bmag(a)**2+gamma*p(a))/a)*ar
    !This is for xi2
    BCB22ac = mt*BCB21ac
    !This is for xi2_prime
    BCB22aa = minval((mt/a-kz*Bt(a)/Bz(a))*(Bmag(a)**2/a*mt+&
      & gamma*p(a)*kz*Bt(a)/Bz(a))-&
      & rho(a)*Vp(a)**2-Bmag(a)**4/Bz(a)**2*(mt**2/a**2+kz**2))*ar
    BCB23a = minval(-gamma*p(a)*(mt*Bt(a)/a+kz*Bz(a))*(mt/a-kz*Bt(a)/Bz(a))-Bt(a)*rho(a)*Vp(a)**2)*ar
    BCB24a = minval(rho(a)*Vp(a))*ar*0.
    !This is for u2_prime
    BCB25aa = minval(-rho(a)*Bmag(a)**2/Bz(a)**2*kV(a))*ar
    !This is for u2
    BCB25ac = minval(rho(a)*Vp(a)*mt)*ar*0.
  END SUBROUTINE init_bc
END MODULE vcyl_funcs_module
