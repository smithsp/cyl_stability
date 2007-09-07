PROGRAM cyl
  USE local
  USE cyl_funcs_module
  USE finite_elements_module
  USE sort_module
  INTEGER :: epsi, ref, start, fin
  LOGICAL :: spline, psi_deriv, chi_deriv
  CHARACTER(LEN=40) :: filename, date, time
!The following are all used in subroutines below and should not be used in the main (cyl) program
  INTEGER :: N, i, j, k, l, nphi, npsi, nchi, INFO, num
  REAL(r8) :: temp, tempA, tempB, tempC, at1, at2, st1, st2, t1, t2
  CHARACTER(LEN=30) :: FMT, FMTR
  INTEGER :: NN
  INTEGER :: LDVL=1, LWORK, LDVR, lower(3), upper(3)
  
  NAMELIST /control_params/  ref, start, fin, verbose
  NAMELIST /cyl_params/ ar, kz, gamma, mt, rho0, eps, Bz0, Bt0, nz, s2, equilib, N,  psi_deriv, chi_deriv, spline, epsilo
  spline = .true.
  verbose = .false.
  psi_deriv = .true.
  chi_deriv = .true.
  CALL cpu_time(t1)
  OPEN(1,file='control_params.in',status='old',form='formatted')
  READ(1,nml=control_params)
  CLOSE(1)
  WRITE(*,nml=control_params) 
  IF (ref.ne.start) THEN
    num = ref
    WRITE(filename,'(a,i0,a)') 'input/',ref,'.in'
    WRITE(*,*)  'Inputting parameters from ', filename
    OPEN(1,file=filename,status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    OPEN(1,file=filename,status='replace',form='formatted')
    WRITE(1,nml=cyl_params)
    CLOSE(1)
    WRITE(*,nml=cyl_params)
    IF (spline) THEN
      CALL bspline_deriv()      
    ELSE 
      CALL linear_const()
    ENDIF
  ENDIF
  DO num = start,fin
    WRITE(filename,'(a,i0,a)') 'input/',num,'.in'
    WRITE(*,*)  'Inputting parameters from ', filename
    OPEN(1,file=filename,status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    OPEN(1,file=filename,status='replace',form='formatted')
    WRITE(1,nml=cyl_params)
    CLOSE(1)
    WRITE(*,nml=cyl_params)
    IF (spline) THEN
      CALL bspline_deriv()      
    ELSE 
      CALL linear_const()
    ENDIF
  ENDDO
  CALL cpu_time(t2)
  WRITE (0,*) 'Time taken: ', t2-t1, ' seconds.'
  WRITE (0,*) 'Time taken by EV solver: ',st, 'seconds.'
  WRITE (0,*) 'Time taken to assemble matrices: ',at, 'seconds.'
  WRITE (0,*) 'Time taken to integrate matrix elements: ',it, 'seconds.'
CONTAINS
  SUBROUTINE linear_const()
    IMPLICIT NONE
    REAL(r8), DIMENSION(N):: grid
    REAL(r8), DIMENSION(N-1) :: xi_1, xi_2, xi_3
    TYPE(linear), DIMENSION(N-1) :: phi 
    TYPE(constant), DIMENSION(N-1) :: psi, chi
    REAL(r8), DIMENSION(3*(N-1),3*(N-1)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,3*(N-1)) :: VL
    INTEGER inds(3*size(phi))
    REAL(r8), DIMENSION(3*(N-1)) :: ALPHAR, ALPHAI, BETA, RWORK(8*3*(N-1)), WORK(10*3*(N-1)), slow, lambda
    REAL(r8) :: xgrid(N*30)
    
  !Initialize the finite elements
    lower = (/lbound(phi,1), lbound(psi,1), lbound(chi,1)/)
    upper = (/ubound(phi,1), ubound(psi,1), ubound(chi,1)/)
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    nphi =  size(phi); npsi = size(psi); nchi = size(chi)

    CALL init(phi(1),grid(1),p3=grid(2))
    CALL init(phi(2:N-1),grid(2:N-1),p2=grid(1:N-2),p3=grid(3:N))
    CALL init(psi(1:N-1),grid(1:N-1),p3=grid(2:N))
    CALL init(chi(1:N-1),grid(1:N-1),p3=grid(2:N))


  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    NN = size(phi)+size(psi)+size(chi)
    LWORK=10*NN
    LDVR = NN
    WRITE(FMT,'(a1,I,a)') '(',NN,'g14.5)'
    WRITE(FMTR,'(a1,I,a)') '(',nphi,'g20.12)'

    A = 0.
    B = 0.
    k = 1
    CALL cpu_time(at1)
    DO i=lower(k),upper(k)
      l = 1
      DO j=i,min(i+phi(i)%extent(2),upper(l))
        temp = (1.+1./mt**2)*int_func(phi(i),phi(j),rho)
        A(3*(i-lower(k))+k,3*(j-lower(l))+l)  = temp
        IF ((i.ne.j).and.(temp.ne.0.))  A(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
        tempA = int_func(phi(i),phi(j),W11A,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi(i),phi(j),W11B,deriv1=.true.)+int_func(phi(i),phi(j),W11B,deriv2=.true.)
        tempC = int_func(phi(i),phi(j),W11C)
        temp = tempA+tempB+tempC
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = temp
        IF ((i.ne.j).and.(temp.ne.0.))  B(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
      ENDDO
      l = 2
      DO j=max(i-phi(i)%extent(1),lower(l)),min(i+phi(i)%extent(2),upper(l))
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = -1./mt**2 * int_func(psi(j),phi(i),K12)
        tempA = int_func(psi(j),phi(i),W12A,deriv2=.true.)
        tempB = int_func(psi(j),phi(i),W12B)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA+tempB
      ENDDO
      l = 3
      DO j=max(i-phi(i)%extent(1),lower(l)),min(i+phi(i)%extent(2),upper(l))
        tempA = int_func(chi(j),phi(i),W13A,deriv2=.true.)
        tempB = int_func(chi(j),phi(i),W13B)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA+tempB
      ENDDO
    ENDDO
    k = 2
    DO i=lower(k),upper(k)
      l = 1
      DO j=max(i-psi(i)%extent(1),lower(l)),min(i+psi(i)%extent(2),upper(l))
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = A(3*(j-lower(l))+l,3*(i-lower(k))+k)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l = 2
      DO j=i,min(i+psi(i)%extent(2),upper(l))
        temp = 1./mt**2 * int_func(psi(i),psi(j),K22)
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = temp
        IF ((i.ne.j).and.(temp.ne.0.))  A(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
        tempA = int_func(psi(i),psi(j),W22A)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA
        IF ((i.ne.j).and.(tempA.ne.0.)) B(3*(j-lower(k))+k,3*(i-lower(l))+l) = tempA
      ENDDO
      l = 3
      DO j=max(i-psi(i)%extent(1),lower(l)),min(i+psi(i)%extent(2),upper(l))
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = int_func(psi(i),chi(j),W23A)
      ENDDO
    ENDDO
    k = 3
    DO i=lower(k),upper(k)
      l = 1
      DO j=max(i-chi(i)%extent(1),lower(l)),min(i+chi(i)%extent(2),upper(l))
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l=2
      DO j=max(i-chi(i)%extent(1),lower(l)),min(i+chi(i)%extent(2),upper(l))
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l = 3
      DO j=i,min(i+chi(i)%extent(2),upper(l))
        temp = 1/kz**2 * int_func(chi(i),chi(j),rho)
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = temp
        IF ((i.ne.j).and.(temp.ne.0.))  A(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
        tempA = int_func(chi(i),chi(j),W33A)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA
        IF ((i.ne.j).and.(tempA.ne.0.))  B(3*(j-lower(k))+k,3*(i-lower(l))+l) = tempA
      ENDDO
    ENDDO
    CALL cpu_time(at2)
    at = at+at2-at1
    IF(verbose) THEN
      OPEN(1,status='replace',file='lin_A.txt')
      WRITE (1,*) 'A='
      WRITE (1,FMT) transpose(A)
      CLOSE(1)
      OPEN(1,status='replace',file='lin_B.txt')
      WRITE (1,*) 'B='
      WRITE (1,FMT) transpose(B)
      CLOSE(1)
    ENDIF
    C=A(:,:)
    D=B(:,:)
    CALL cpu_time(st1)
    CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs_linconst(phi,psi,chi,grid,VR)
  END SUBROUTINE linear_const
  
  SUBROUTINE bspline_deriv()
    IMPLICIT NONE
    REAL(r8), DIMENSION(N):: grid!, slow!, slow_sort
    REAL(r8), DIMENSION(N-1) :: xi_1, xi_2, xi_3
    TYPE(bspline), DIMENSION(0:N) :: phi
    TYPE(bspline), DIMENSION(0:N) :: psi, chi 
    REAL(r8), DIMENSION(size(phi)+size(psi)+size(chi),size(phi)+size(psi)+size(chi)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,size(phi)+size(psi)+size(chi)) :: VL
    INTEGER ::  inds(size(phi)+size(psi)+size(chi))
    REAL(r8), DIMENSION(size(phi)+size(psi)+size(chi)) :: ALPHAR, ALPHAI, BETA, slow,lambda,&
    & RWORK(8*(size(phi)+size(psi)+size(chi))), WORK(10*(size(phi)+size(psi)+size(chi)))
    LOGICAL :: Lend0
    REAL(r8) :: xgrid(N*30)
 
  !Initialize the finite elements
    Lend0 = .true.
    lower = (/lbound(phi,1), lbound(psi,1), lbound(chi,1)/)
    upper = (/ubound(phi,1), ubound(psi,1), ubound(chi,1)/) 
    grid(1:N) = (/ (i*ar/(N-1), i=0,N-1) /)
    xgrid = (/ (i*ar/real(size(xgrid)-1), i=0,size(xgrid)-1) /)
    nphi =  size(phi); npsi = size(psi); nchi = size(chi)

    CALL init(phi(0),grid(1),p3=grid(2),LendZero=.true.)
    CALL init(phi(1),grid(1),p3=grid(2),p4=grid(3),LendZero=.true.)
    CALL init(phi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=.true.)
    CALL init(phi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    CALL init(phi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
  
    psi = phi
    psi%deriv = psi_deriv
    chi=psi
    chi%deriv = chi_deriv

  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    NN = size(phi)+size(psi)+size(chi)
    LWORK=10*NN
    LDVR = NN
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi,'G20.11)'

    A = 0.
    B = 0.
    k = 1
    CALL cpu_time(at1)
    DO i=lower(k),upper(k)
      l = 1
      DO j=i,min(i+3,upper(l))
        temp = (1.+1./mt**2)*int_func(phi(i),phi(j),rho)
        A(3*(i-lower(k))+k,3*(j-lower(l))+l)  = temp
        IF ((i.ne.j))  A(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
        tempA = int_func(phi(i),phi(j),W11A,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi(i),phi(j),W11B,deriv1=.true.)+int_func(phi(i),phi(j),W11B,deriv2=.true.)
        tempC = int_func(phi(i),phi(j),W11C)
        temp = tempA+tempB+tempC
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = temp
        IF ((i.ne.j))  B(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = -1./mt**2 * int_func(psi(j),phi(i),K12)
        tempA = int_func(psi(j),phi(i),W12A,deriv2=.true.)
        tempB = int_func(psi(j),phi(i),W12B)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA+tempB
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(chi(j),phi(i),W13A,deriv2=.true.)
        tempB = int_func(chi(j),phi(i),W13B)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA+tempB
      ENDDO
    ENDDO
    k = 2
    DO i=lower(k),upper(k)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = A(3*(j-lower(l))+l,3*(i-lower(k))+k)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l = 2
      DO j=i,min(i+3,upper(l))
        temp = 1./mt**2 * int_func(psi(i),psi(j),K22)
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = temp
        IF ((i.ne.j))  A(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
        tempA = int_func(psi(i),psi(j),W22A)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA
        IF ((i.ne.j)) B(3*(j-lower(k))+k,3*(i-lower(l))+l) = tempA
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = int_func(psi(i),chi(j),W23A)
      ENDDO
    ENDDO
    k = 3
    DO i=lower(k),upper(k)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l=2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l = 3
      DO j=i,min(i+3,upper(l))
        temp = 1/kz**2 * int_func(chi(i),chi(j),rho)
        A(3*(i-lower(k))+k,3*(j-lower(l))+l) = temp
        IF ((i.ne.j))  A(3*(j-lower(k))+k,3*(i-lower(l))+l) = temp
        tempA = int_func(chi(i),chi(j),W33A)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = tempA
        IF ((i.ne.j))  B(3*(j-lower(k))+k,3*(i-lower(l))+l) = tempA
      ENDDO
    ENDDO
    CALL cpu_time(at2)
    at = at+at2-at1
    IF(verbose) THEN
      OPEN(1,status='replace',file='A.txt')
      WRITE (1,*) 'A='
      WRITE (1,FMT) transpose(A)
      CLOSE(1)
      OPEN(1,status='replace',file='B.txt')
      WRITE (1,*) 'B='
      WRITE (1,FMT) transpose(B)
      CLOSE(1)
    ENDIF
    C=A(:,:)
    D=B(:,:)
    CALL cpu_time(st1)
    CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs_spline(phi,psi,chi,grid,VR)    
  END SUBROUTINE bspline_deriv
  
  SUBROUTINE plot_equilibrium
    INTEGER, PARAMETER :: num_pts = 300
    INTEGER :: i
    REAL(r8), DIMENSION(num_pts) :: grid, pressure, zb, tb, density, pv, zv, p_prime, Bt_prime, Bz_prime, safety
    grid = (/ (i*ar/(num_pts-1), i=0,num_pts-1) /)
    WRITE (*,nml=cyl_params)
    WRITE (*,'(a,i)') 'Plotting equilib = ',equilib
    WRITE (*,*) 'Eqilibrium condition at grid points='
    WRITE (*,'(g)') equilibrium(grid(2:))
    pressure = P(grid)
    tb = Bt(grid)
    zb = Bz(grid)
    density = rho(grid)
    pv = 0
    zv = 0
    safety = q(grid)
    p_prime = diff(P,grid)
    Bz_prime = diff(Bz,grid)
    Bt_prime = diff(Bt,grid)
    WRITE (*,'(5(a20))') 'pressure','p_prime', 'Bz*Bz_prime', 'Bt*Bt_prime', 'Bt**2/r'
    DO i=2,size(grid)
      WRITE (*,'(5(g20.10))') pressure(i), p_prime(i), zb(i)*Bz_prime(i), tb(i)*Bt_prime(i), tb(i)**2/grid(i)
    ENDDO
    WRITE(filename,'(a,i0,a)') 'equilib',equilib,'.txt'
    OPEN(1,file=filename,status='replace')
    DO i=1,num_pts
      WRITE (1,'(g,g,g,g,g,g,g,g)') grid(i), pressure(i),zb(i),tb(i),density(i),pv(i),zv(i), safety(i)
    ENDDO
    CLOSE(1)
  END SUBROUTINE plot_equilibrium

  SUBROUTINE output_params(assemble_t,solve_t)
    REAL(r8), INTENT(IN) :: assemble_t,solve_t
    WRITE(filename,'(a,i0,a)') 'output_cyl/',num,'.dat'
    OPEN(1,file=filename,status='replace')
    IF (spline) THEN
      WRITE(1,'(a)') 'spline'
    ELSE 
      WRITE(1,'(a)') 'linconst'
    ENDIF
    WRITE(1,'(i)') N, NN, mt, equilib, num
    WRITE(1,'(g)') assemble_t, solve_t, kz, ar, rho0
    SELECT CASE (equilib)
      CASE(1)
        WRITE(1,'(g)') Bz0, Bt0, s2
      CASE(2)
        WRITE(1,'(g)') Bz0, Bt0, s2, eps
      CASE(3)
        WRITE(1,'(g)') Bz0, Bt0
        WRITE(1,'(i)') nz
      CASE(4)
        WRITE(1,'(g)') P0, P1, lambd
      CASE(5)
        WRITE(1,'(g)') Bz0, Bt0
        WRITE(1,'(i)') nz
    ENDSELECT
    CLOSE(1)
  END SUBROUTINE output_params
  SUBROUTINE output_evals(evalsr,evalsi)
    REAL(r8), DIMENSION(:), INTENT(IN) :: evalsr, evalsi
    WRITE(filename,'(a,i0,a)') 'output_cyl/',num,'.evalsr'
    OPEN(1,file=filename,status='replace')
    WRITE(1,'(g)') evalsr
    CLOSE(1)
    WRITE(filename,'(a,i0,a)') 'output_cyl/',num,'.evalsi'
    OPEN(1,file=filename,status='replace')
    WRITE(1,'(g)') evalsi
    CLOSE(1)    
  END SUBROUTINE output_evals
  SUBROUTINE output_evecs_spline(phi1,phi2,phi3,grid,VR)
    TYPE(bspline), DIMENSION(:), INTENT(IN) :: phi1, phi2, phi3
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: VR
    TYPE(bspline), DIMENSION(size(phi1),3) :: phi
    INTEGER :: N , i, j, k
    CHARACTER(LEN=30) FMT
    phi = reshape((/phi1,phi2,phi3/),(/size(phi1),3/))
    N = size(grid)
    WRITE(FMT,'(a,i0,a)') '(',N,'(g,a))'
    WRITE(filename,'(a,i0,a)') 'output_cyl/',num,'.grid'
    OPEN(1,file=filename,status='replace')
    WRITE(1,FMT) (grid(i),',',i=1,size(grid))
    CLOSE(1)
    DO j=1,3
      WRITE(filename,'(2(a,i0))') 'output_cyl/',num,'.evecs',j
      OPEN(1,file=filename,status='replace')
      DO k=1,size(VR(1,:))
        WRITE(1,FMT) ( sum(VR(j::3,k)*val(phi(:,j),grid(i))),',' , i=1,size(grid) )
      ENDDO
      CLOSE(1)
    ENDDO 
  END SUBROUTINE output_evecs_spline
  SUBROUTINE output_evecs_linconst(phi1,phi2,phi3,grid,VR)
    TYPE(linear), DIMENSION(:), INTENT(IN) :: phi1
    TYPE(constant), DIMENSION(:), INTENT(IN) :: phi2, phi3
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: VR
    TYPE(constant), DIMENSION(size(phi2),2:3) :: phi
    INTEGER :: N , i, j, k
    CHARACTER(LEN=30) FMT
    phi = reshape((/phi2,phi3/),(/size(phi2),2/))
    N = size(grid)
    WRITE(FMT,'(a,i0,a)') '(',N,'(g,a))'
    WRITE(filename,'(a,i0,a)') 'output_cyl/', num,'.grid'
    OPEN(1,file=filename,status='replace')
    WRITE(1,FMT) (grid(i),',',i=1,size(grid))
    CLOSE(1)
    j=1
    WRITE(filename,'(2(a,i0))') 'output_cyl/', num,'.evecs',j
    OPEN(1,file=filename,status='replace')
    DO k=1,size(VR(1,:))
      WRITE(1,FMT) ( sum(VR(j::3,k)*val(phi1,grid(i))),',' , i=1,size(grid) )
    ENDDO
    CLOSE(1)
    DO j=2,3
      WRITE(filename,'(2(a,i0))') 'output_cyl/', num,'.evecs',j
      OPEN(1,file=filename,status='replace')
      DO k=1,size(VR(1,:))
        WRITE(1,FMT) ( sum(VR(j::3,k)*val(phi(:,j),grid(i))),',' , i=1,size(grid) )
      ENDDO
      CLOSE(1)
    ENDDO 
  END SUBROUTINE output_evecs_linconst
END PROGRAM cyl
