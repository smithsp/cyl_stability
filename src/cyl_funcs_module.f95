MODULE cyl_funcs_module
  USE local
  IMPLICIT NONE
  REAL(r8), EXTERNAL :: s17aef, s17aff
  INTEGER :: mt, equilib, ifail  !m_theta, choice of equilibrium configuration, 
  REAL(r8) :: kz, gamma, ar, br, rho0, eps, Bz0, Bt0, s2, P0, P1,lambd !ar is the radius of the plasma, br is the radius of the wall
  REAL(r8) :: Vz0, epsVz, Vp0, epsVp, nu
  REAL(r8) :: rs, alpha,w, maxw, minw, h ! These are parameters for localizing the grid around rs
  REAL(r8), PARAMETER, DIMENSION(2) :: gamma_mt_1 = (/1.841183781,3.054236928 /)

CONTAINS
  FUNCTION kdotB(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: kdotB
    kdotB = mt*Bt(r)/r+kz*Bz(r)
  END FUNCTION kdotB
  FUNCTION alfven_range(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(2) :: alfven_range
    alfven_range(1) = minval(1/(rho(r))* kdotB(r)**2)
    alfven_range(2) = maxval(1/(rho(r))* kdotB(r)**2)
    RETURN    
  END FUNCTION alfven_range
  FUNCTION slow_inf_range(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(2) :: slow_inf_range
    slow_inf_range(1) = minval(gamma*P(r)*kdotB(r)**2/(rho(r)*(Bmag(r)**2+gamma*P(r))))
    slow_inf_range(2) = maxval(gamma*P(r)*kdotB(r)**2/(rho(r)*(Bmag(r)**2+gamma*P(r))))
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
      CASE(2,10)
        rho = rho0 * (1-eps*(r**2/ar**2))
      CASE(3:4,6:9,12,13)
        rho = rho0
      CASE(5)
        rho = P(r)*rho0
      CASE(11)
        rho = 1.+(rho0-1.)*r**2
    ENDSELECT
    RETURN
  END FUNCTION rho
  FUNCTION Bz(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bz
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:3,5:6)
        Bz = Bz0
      CASE(4)
        DO i=1,size(r)
          Bz(i) = sqrt(1.-P1)* s17aef(lambd*r(i),ifail)
        ENDDO
      CASE(7)
        Bz = sqrt(Bz0**2-2*p0*exp(-r**2/ar**2))
      CASE(8)
        Bz = 0
      CASE(9)
        Bz = sqrt(2.0)*sqrt(Bt0**2*(2*ar**2-r**2))/ar
      CASE(10,12)
        Bz = sqrt((Bz0**2-2*p0+2*Bt0**2)*ar**2+2*(p0-Bt0**2)*r**2)/ar
      CASE(11)
        Bz = sqrt(2.0)*sqrt(0.1D1/gamma*P0)*(1.-Bz0*r**2)
      CASE(13)
        Bz = (1./6.)*sqrt(-15.*Bt0**2*eps/ar**8&
                         &*(-25.*ar**8+48.*r**2*ar**6-36.*r**4*ar**4+16.*r**6*ar**2-3.*r**8&
                           & +8.*ar**8*((1.-r**2/ar**2))**(3/2.)+24*ar**8*sqrt((1.-r**2/ar**2))&
                           & -12*ar**8*(log(2*ar**2+2*sqrt((1.-r**2/ar**2))*ar**2-r**2)-log(ar**2)))+36*Bz0**2)
    ENDSELECT
    RETURN
  END FUNCTION Bz
  FUNCTION Bz_prime(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bz_prime
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:3,5:6,8)
        Bz_prime = 0
      CASE(4)
        DO i=1,size(r)
          Bz_prime(i) = -sqrt(1.-P1)* s17aff(lambd*r(i),ifail)*lambd
        ENDDO
      CASE(7)
        Bz_prime = 2*(Bz0**2-2*p0*exp(-r**2/ar**2))**(-.5)*p0*r&
        & /ar**2*exp(-r**2/ar**2)
      CASE(9)
        Bz_prime = -sqrt(2.0)*(Bt0**2*(2*ar**2-r**2))**(-.5)/ar*Bt0**2*r
      CASE(10,12)
        Bz_prime = 2*r*(p0-Bt0**2)*(ar**2*Bz0**2-2*ar**2*p0+2*ar**2*Bt0**2+&
        & 2*r**2*p0-2*r**2*Bt0**2)**(-0.5)/ar
      CASE(11)
        Bz_prime = -2*sqrt(2.0)*sqrt(0.1D1/gamma*P0)*Bz0*r
      CASE(13)
        Bz_prime = -0.5D1/0.4D1*Bt0**2/ar**8*eps&
                    &*(-15*Bt0**2/ar**8*eps&
                      &*(-25*ar**8+48*r**2*ar**6-36*r**4*ar**4+16*r**6*ar**2&
                        &-3*r**8+8*ar**8*(1.-r**2/ar**2)**(0.3D1/0.2D1)&
                        &+24*ar**8*sqrt(1.-r**2/ar**2)&
                        &-12*ar**8*log(2+2*sqrt(1.-r**2/ar**2)-r**2/ar**2))&
                      &+36*Bz0**2)**(-0.1D1/0.2D1)&
                    &*(96*r*ar**6-144*r**3*ar**4+96*r**5*ar**2-24*r**7&
                      &-24*r*ar**4*(2*ar**2-r**2)/sqrt(1-r**2/ar**2)&
                      &+24*r*ar**6*(1+ar/sqrt(ar**2-r**2))&
                        &/(2+2*sqrt((1.-r**2/ar**2))-r**2/ar**2))
    ENDSELECT
    RETURN
  END FUNCTION Bz_prime
  FUNCTION Bt(r) !for B_theta
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bt
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:2,7)
        Bt = 0
      CASE(3,5,9,10,12)
        Bt = Bt0*r/ar
      CASE(4)
        DO i=1,size(r)
          Bt(i) = s17aff(lambd*r(i),ifail)
        ENDDO
      CASE(6)
        Bt = sqrt(P0*eps)*r/ar
      CASE(8)
        Bt = 2*Bt0*r/(r**2+P0**2)
      CASE(11)
        Bt = Bt0*r*(1.-lambd*r**2/2.)/2.
      CASE(13)
        Bt = Bt0*(ar-(ar**2-r**2)**2*sqrt(1-r**2/ar**2)/ar**3)/r
    ENDSELECT
    RETURN
  END FUNCTION Bt
  FUNCTION Bt_prime(r) !for B_theta
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: Bt_prime
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:2)
        Bt_prime = 0
      CASE(3,5,9,10,12)
        Bt_prime = Bt0/ar
      CASE(4)
        DO i=1,size(r)
          Bt_prime(i) = (s17aef(lambd*r(i),ifail)-s17aff(lambd*r(i),ifail)/(lambd*r(i)))*lambd
        ENDDO
      CASE(6)
        Bt_prime = sqrt(P0*eps)/ar
      CASE(7)
        Bt_prime = 0
      CASE(8)
        Bt_prime = 2.*Bt0/(r**2+P0**2)-4.*Bt0*r**2/(r**2+P0**2)**2
      CASE(11)
        Bt_prime = Bt0*(1-lambd*r**2/2.)/2.-Bt0*r**2*lambd/2.
      CASE(13)
        Bt_prime = Bt0*(-ar+sqrt(1-r**2/ar**2)*(ar**4+3*r**2*ar**2-4*r**4)/ar**3)/r**2
    ENDSELECT
    RETURN
  END FUNCTION Bt_prime
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
      CASE(3,5)
        P = Bt0**2*(1-r**2/ar**2)
      CASE(4)
        DO i=1,size(r)
          P(i) = P0+P1/2.*s17aef(lambd*r(i),ifail)**2
        ENDDO
      CASE(6)
        P = P0*(1-eps*r**2/ar**2)
      CASE(7)
        P = P0*exp(-r**2/ar**2)
      CASE(8)
        P = 2*Bt0**2*P0**2/(r**2+P0**2)**2
      CASE(9)
        P = P0
      CASE(10)
        P = P0*(1-r**2/ar**2)
      CASE(11)
        P = eps**2*(rho0-1)*r**8/8&
        & +0.2D1/0.7D1*epsVp*eps*(rho0-1)*r**7&
        & +(-Bt0**2*lambd**2/24+(2*Vp0*eps+epsVp**2)*(rho0-1)/6+eps**2/6)*r**6&
        & +(0.2D1/0.5D1*Vp0*epsVp*(rho0-1)+0.2D1/0.5D1*epsVp*eps)*r**5&
        & +(Vp0**2*(rho0-1)/4+Vp0*eps/2+epsVp**2/4&
        & +0.3D1/0.16D2*Bt0**2*lambd-0.1D1/gamma*P0*Bz0**2)*r**4&
        & +0.2D1/0.3D1*Vp0*epsVp*r**3&
        & +(-Bt0**2/4+Vp0**2/2+2/gamma*P0*Bz0)*r**2&
        & -0.1D1/gamma*P0+P1
      CASE(12)
        P = -(1-r**2/ar**2)/2.*(ar**2*Vp0**2-2.*P0)
      CASE(13)
        P = -0.5D1/0.24D2*(1-eps)*Bt0**2/ar**8&
        &*(-25*ar**8+48*r**2*ar**6-36*r**4*ar**4+16*r**6*ar**2-3*r**8&
          &+32*ar**8*sqrt(1-r**2/ar**2)-8*ar**6*sqrt(1-r**2/ar**2)*r**2&
          &-12*ar**8*log(2+2*sqrt(1/ar**2*(ar**2-r**2))-r**2/ar**2))
    ENDSELECT
    RETURN
  END FUNCTION P
  FUNCTION P_prime(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: P_prime
    INTEGER :: i
    SELECT CASE (equilib)
      CASE(1:2)
        P_prime = 0  
      CASE(3,5)
        P_prime = -Bt0**2*(2*r/ar**2)
      CASE(4)
        DO i=1,size(r)
          P_prime(i) = -P1*lambd*s17aef(lambd*r(i),ifail)*s17aff(lambd*r(i),ifail)
        ENDDO
      CASE(6)
        P_prime = P0*(-2*eps*r/ar**2)
      CASE(7)
        P_prime = -2*P0*exp(-r**2/ar**2)*r/ar**2
      CASE(8)
        P_prime = -8*Bt0**2*P0**2/(r**2+P0**2)**3 *r
      CASE(9)
        P_prime = 0
      CASE(10)
        P_prime = -2*p0*r/ar**2
      CASE(11)
        P_prime = eps**2*(rho0-1)*r**7&
        & +2*epsVp*eps*(rho0-1)*r**6&
        & +6*(-Bt0**2*lambd**2/24+(2*Vp0*eps+epsVp**2)*(rho0-1)/6+eps**2/6)*r**5&
        & +5*(0.2D1/0.5D1*Vp0*epsVp*(rho0-1)+0.2D1/0.5D1*epsVp*eps)*r**4&
        & +4*(Vp0**2*(rho0-1)/4+Vp0*eps/2+epsVp**2/4+0.3D1/0.16D2*Bt0**2*lambd-0.1D1/Gamma*P0*Bz0**2)*r**3&
        & +2*Vp0*epsVp*r**2+2*(-Bt0**2/4+Vp0**2/2+2/Gamma*P0*Bz0)*r
      CASE(12)
        P_prime = (ar**2*Vp0**2-2*P0)*r/ar**2
      CASE(13)
        P_prime = -0.5D1/0.24D2*(1-eps)*Bt0**2/ar**8&
                    &*(96*r*ar**6-144*r**3*ar**4+96*r**5*ar**2-24*r**7&
                      &-24*r*ar**4*(2*ar**2-r**2)/sqrt(1-r**2/ar**2)&
                      &+24*r*ar**6*(1+ar/sqrt(ar**2-r**2))&
                        &/(2+2*sqrt((1.-r**2/ar**2))-r**2/ar**2))

    ENDSELECT
    RETURN
  END FUNCTION P_prime
  FUNCTION q(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: q    
    q = kz*r*Bz(r)/Bt(r)
    RETURN
  END FUNCTION q
  FUNCTION equilibrium(r)
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(size(r)) :: equilibrium
    equilibrium = P_prime(r) + Bz(r)*Bz_prime(r)+Bt(r)*Bt_prime(r)+Bt(r)**2/r
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
  FUNCTION new_grid(grid)
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(size(grid)) :: new_grid
    INTEGER :: rs_ind, ng
    ng = size(grid)  
    IF (verbose) THEN 
      WRITE(*,*) 'The value of nu used to determine whether or not to modify it was', nu
    ENDIF
    IF (rs.lt.0 .or. rs.gt.ar.or.alpha.eq.0 ) THEN
      new_grid(1) = new_grid(2)/1000.0
      new_grid = grid
      h = minval(new_grid(3:)-new_grid(2:ng-1))
      IF (abs(nu)>1) THEN ! This makes it possible to specify nu directly
        nu = -abs(nu)*h**(-2.7)*exp(-23.95)
      ENDIF
      RETURN
    ENDIF
    minw=0
    maxw=.04/alpha
    new_grid = findw(ng)
    h = minval(new_grid(3:)-new_grid(2:ng-1))    
    IF (abs(nu)>1) THEN ! This makes it possible to specify nu directly
      nu = -abs(nu)*h**(-2.7)*exp(-23.95)
    ENDIF
    RETURN    
  END FUNCTION new_grid
  FUNCTION localize(x)
    REAL(r8), INTENT(IN) :: x
    REAL(r8) :: localize
    localize = alpha*w/((x-rs)**2+w)
    RETURN 
  END FUNCTION localize
  RECURSIVE FUNCTION findw(ngrid) RESULT(new_grid)
    INTEGER, INTENT(IN) :: ngrid
    REAL(r8), DIMENSION(ngrid) :: new_grid
    REAL(r8), DIMENSION(ngrid+10) :: grid
    INTEGER :: i, ng
    REAL(r8), DIMENSION(1) :: temp
    w = (minw+maxw)/2.
    grid(1) = 0.
    DO i=2,size(grid)
      grid(i) = grid(i-1)+1./((localize(grid(i-1))+1./(alpha+1.))*ngrid)
    ENDDO
    ng = count(grid<=ar)
    IF (verbose) THEN
      WRITE (*,*) grid
    ENDIF
    IF (ng==ngrid) THEN
      IF (ar-grid(ngrid)<1e-4) THEN
        new_grid = grid(1:ngrid)
        temp = (grid(minloc(abs(grid-rs))))
        new_grid = grid(1:ngrid)-temp(1)+rs
        new_grid(1) = new_grid(2)/1000.0
        new_grid(ngrid) = ar
        RETURN
      ELSE
        maxw = w
        new_grid = findw(ngrid)
      ENDIF
    ELSEIF (ng<ngrid) THEN
      minw = w
      new_grid = findw(ngrid)
    ELSEIF (ng>ngrid) THEN
      maxw = w
      new_grid = findw(ngrid)
    ENDIF
  END FUNCTION findw
  SUBROUTINE calc_rs(grid)
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(size(grid)-1) :: y
    INTEGER :: rs_ind, ng
    ng = size(grid)
    y = (mt*Bt(grid(2:))+kz*Bz(grid(2:))*grid(2:))
    rs_ind = minval(minloc(y(2:ng-1)*y(3:)))
    rs = (grid(rs_ind+1)*y(rs_ind+1)-y(rs_ind)*grid(rs_ind+2))/(y(rs_ind+1)-y(rs_ind)) !Linear interpolation
  END SUBROUTINE calc_rs
END MODULE cyl_funcs_module
