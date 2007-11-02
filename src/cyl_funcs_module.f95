MODULE cyl_funcs_module
  USE local
  IMPLICIT NONE
  REAL(r8), EXTERNAL :: s17aef, s17aff
  INTEGER :: mt, equilib, ifail, nz  !m_theta, choice of equilibrium configuration, 
  REAL(r8) :: kz, gamma, ar, br, rho0, eps, Bz0, Bt0, s2, P0, P1,lambd !ar is the radius of the plasma, br is the radius of the wall
  REAL(r8) :: Vz0, epsVz, Vp0, epsVp
  REAL(r8) :: rs, alpha ! These are parameters for localizing the grid around rs
  REAL(r8), PARAMETER, DIMENSION(2) :: gamma_mt_1 = (/1.841183781,3.054236928 /)
  LOGICAL :: kB  !This is a switch to calculate the resonant surface at which kB=0 (instead of where k_\parallel=0)
CONTAINS
  FUNCTION alfven_range(r)
    IMPLICIT NONE
    REAL(r8), INTENT(IN), DIMENSION(:) :: r
    REAL(r8), DIMENSION(2) :: alfven_range
    alfven_range(1) = minval(1/(rho(r))* kz**2*(Bmag(r))**2 )
    alfven_range(2) = maxval(1/(rho(r))* kz**2*(Bmag(r))**2 )
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
      CASE(2)
        rho = rho0 * (1-eps*(r**2/ar**2))
      CASE(3:4,6:10)
        rho = rho0
      CASE(5)
        rho = P(r)*rho0
      CASE(11)
        rho = 1+(rho0-1)*r**2
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
      CASE(10)
        Bz = sqrt((Bz0**2-2*p0+2*Bt0**2)*ar**2+2*(p0-Bt0**2)*r**2)/ar
      CASE(11)
        Bz = sqrt(2.0)*sqrt(0.1D1/gamma*P0)*(1-Bz0*r**2)
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
      CASE(10)
        Bz_prime = 2*r*(p0-Bt0**2)*(ar**2*Bz0**2-2*ar**2*p0+2*ar**2*Bt0**2+&
        & 2*r**2*p0-2*r**2*Bt0**2)**(-0.5)/ar
      CASE(11)
        Bz_prime = -2*sqrt(2.0)*sqrt(0.1D1/gamma*P0)*Bz0*r
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
      CASE(3,5,9,10)
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
      CASE(3,5,9,10)
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
        Bt_prime = Bt0*(1-lambd*r**2/2)/2-Bt0*r**2*lambd/2
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
        & -0.1D1/gamma*P0
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
          P_prime(i) = -P1*lambd*s17aef(lambd*r(i),ifail)*s17aff(lambd*r(i))
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
    IF (rs.le.0 .or. rs.ge.ar ) THEN
      new_grid = grid
      RETURN
    ENDIF
    rs_ind = minval(minloc(abs(grid-rs)))
    IF (rs-grid(rs_ind)<0) THEN
      rs_ind=  rs_ind-1
    ENDIF
    IF (rs.ne.0) THEN
      new_grid(1:rs_ind) = -rs*((rs-grid(1:rs_ind))/rs)**alpha+rs
    ENDIF
    IF (rs.ne.ar) THEN
      new_grid(rs_ind+1:) = (ar-rs)*((grid(rs_ind+1:)-rs)/(ar-rs))**alpha+rs
    ENDIF
    write(*,*) 'new_grid = ', new_grid
    RETURN    
  END FUNCTION new_grid
  SUBROUTINE calc_rs(grid)
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(size(grid)-1) :: y
    INTEGER :: rs_ind, ng
    ng = size(grid)
    rs = -mt/kz!minval(grid(minloc(abs(mt+kz*grid(2:)))+1))
    IF (kB) THEN
      y = (mt*Bt(grid(2:))+kz*Bz(grid(2:))*grid(2:))
      rs_ind = minval(minloc(y(2:ng-1)*y(3:)))
      rs = (grid(rs_ind+1)*y(rs_ind+1)-y(rs_ind)*grid(rs_ind+2))/(y(rs_ind+1)-y(rs_ind)) !Linear interpolation
    ENDIF
  END SUBROUTINE calc_rs
END MODULE cyl_funcs_module
