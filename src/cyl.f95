PROGRAM cyl
  USE local
  USE cyl_funcs_module
  USE finite_elements_module
  USE sort_module
  INTEGER :: N, min_N, max_N, epsi
  LOGICAL :: linconst, spline, hermite, verbose, slow_evals, slow_evecs, alfven_evecs, psi_deriv, chi_deriv, homo_plasma, inhomo
!The following are all used in subroutines below and should not be used in the main (cyl) program
  INTEGER :: i, j, k, l, nphi, npsi, nchi, INFO, pick_val, ind4(1), ind8(1), stat, ind(1)
  REAL(r8) :: temp, tempA, tempB, tempC, areps, slow_inf
  CHARACTER(LEN=30) :: FMT, FMTR
  INTEGER :: NN
  INTEGER :: LDVL=1, LWORK, LDVR, lower(3), upper(3)
  
  NAMELIST /control_params/ min_N, max_N, linconst, spline, hermite, verbose, slow_evals, slow_evecs, alfven_evecs,&
  & psi_deriv, chi_deriv, homo_plasma, inhomo
  NAMELIST /cyl_params/ ar, kz, gamma, mt, rho0, eps, homo, Bz0, Bt0
  min_N = 5
  max_N = 6
  linconst = .true.
  spline = .true.
  hermite = .true.
  verbose = .false.
  slow_evals = .true.
  slow_evecs = .false.
  psi_deriv = .true.
  chi_deriv = .true.
  homo_plasma = .true.
  inhomo = .false.
  OPEN(1,file='control_params.in',status='old',form='formatted')
  READ(1,nml=control_params)
  CLOSE(1)
  WRITE(*,nml=control_params) 
  IF(homo_plasma) THEN
    OPEN(1,file='homogeneous.in',status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    WRITE (*,nml=cyl_params)
    OPEN(2,file='lin_slow_evals.txt',status='replace')
    DO N = min_N,max_N
      IF (linconst) THEN
        WRITE (*,*) 'N = ', N
        WRITE (*,*) 'With linear and constant elements.'
        CALL linear_const()
      ENDIF
    ENDDO
    CLOSE(2)
    IF (linconst.and.slow_evecs) THEN
      slow_evals = .false.
      CALL linear_const()
    ENDIF
    OPEN(1,file='control_params.in',status='old',form='formatted')
    READ(1,nml=control_params)
    CLOSE(1)
    WRITE(*,nml=control_params) 
    OPEN(2,file='spline_slow_evals.txt',status='replace')
    DO N = min_N,max_N
      IF (spline) THEN
        WRITE (*,*) 'N = ', N
        WRITE (*,*) 'With bspline elements.'
        CALL bspline_deriv()
      ENDIF
    ENDDO
    CLOSE(2)
    IF (spline.and.slow_evecs) THEN
      slow_evals = .false.
      CALL bspline_deriv()
    ENDIF
  ENDIF
  IF(inhomo) THEN
    OPEN(1,file='inhomogeneous.in',status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    WRITE (*,nml=cyl_params)
    N = min_N
    IF (alfven_evecs) THEN
      IF(linconst)  THEN
        OPEN (3, status='replace', file='lin_alfven_EVs.txt')
        DO N=min_N,max_N,10
          WRITE (*,*) 'N = ', N
          WRITE (*,*) 'With linear and constant elements.'
          CALL linear_const()
        ENDDO
        CLOSE(3)
      ENDIF
      IF(spline) THEN
        OPEN (3, status='replace',file='spline_alfven_EVs.txt')
        DO N=min_N,max_N,10
          WRITE (*,*) 'N = ', N
          WRITE (*,*) 'With bspline elements.'
          CALL bspline_deriv()
        ENDDO
        CLOSE(3)
      ENDIF
    ENDIF
    IF (slow_evals) THEN
      N = 20
      IF(linconst) THEN
      ENDIF
      IF(spline) THEN
        OPEN (2, status='replace',file='spline_slow_vareps.txt')
        eps = 0
        CALL bspline_deriv()
        DO epsi = -10,-4
          eps = 10**(epsi/2.0)
          CALL bspline_deriv()
        ENDDO
        CLOSE(3)
      ENDIF
    ENDIF
  ENDIF
  IF (hermite) THEN
    WRITE (*,*) 'With Hermite elements.'
    CALL hermite_elements()
  ENDIF
  WRITE (*,*) 'Alfven range (approx):', alfven_range((/0.,ar/))
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
    CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    !WRITE (*,*) 'INFO =', INFO
    !WRITE (*,*) 'real(ALPHA) ='
    !WRITE (*,'(g)') ALPHAR
    !WRITE (*,*) 'BETA ='
    !WRITE (*,'(g)') BETA
    !WRITE (*,*) 'real(LAMBDA) = '
    !WRITE (*,'(g)') (ALPHAR)/BETA
    !WRITE (*,*) 's^2 = ', gamma*P(grid)/(Bz(grid))**2
    lambda = sort(ALPHAR/BETA,rev=.true.)
    IF(.not.homo) THEN
      WRITE (*,*) 'alfven range = ', alfven_range(grid)
      WRITE (*,*) 'slow_inf range = ', slow_inf_range(grid)
      WRITE (*,*) 'Lambda = '
      WRITE (*,'(g)') lambda
    ENDIF
    slow_inf = minval(slow_inf_range(grid))
    !WRITE (*,*) '(real(LAMBDA) - slow_inf)/slow_inf= '
    slow = sort(((ALPHAR)/BETA-slow_inf)/slow_inf,rev=.true.,inds=inds)
    IF (slow_evals) WRITE(2,FMTR) slow(NN-nphi+1:NN)
    WRITE (*,*) '(lambda-slow_inf)/slow_inf = '
    DO i=1,NN
      !IF (slow(i) .lt. (((minval(kz**2 *Bz0**2/rho((/0./))))*0.9-slow_inf)/slow_inf)) 
      WRITE (*,'(g)')  slow(i)
    ENDDO
    WRITE (*,*) ' '
    !WRITE (*,*) 'k^2 Bz(grid)^2/rho(grid) = '
    !WRITE (*,'(g)') kz**2*(Bz(grid))**2/(rho(grid))
    !WRITE (*,*) 'B V1 - L1 A V1='
    pick_val = NN
    !WRITE (*,'(g,g)') matmul(D,VR(:,pick_val))-(ALPHAR(pick_val)+(0,1)*ALPHAI(pick_val))/BETA(pick_val)*matmul(C,VR(:,pick_val))

    IF(slow_evecs.and.(.not.slow_evals)) THEN
      OPEN (1, status='replace',file='lin_slow_EVs.txt')
      WRITE (1,'(i)') size(grid)-1
      WRITE (1,'(i)') size(phi)
      DO i=1,size(grid)-1
        WRITE (1,'(g)') (grid(i)+grid(i+1))/2.
      ENDDO
      DO j=size(phi)-1,0,-1
        WRITE(1,*) ''  
        DO i=1,size(grid)-1
          WRITE (1,'(g)') sum(VR(3::3,inds(NN-j))*val(chi,(grid(i)+grid(i+1))/2))
        ENDDO    
      ENDDO
      CLOSE(1)
    ENDIF
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    pick_val = NN
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    ind = minloc(abs(ALPHAR/BETA-0.513e-04))
    !WRITE (*,*) 'Lambda(',ind(1),') = ',ALPHAR(ind(1))/BETA(ind(1)) 
    !WRITE (*,'(g)') (VR(2::3,ind(1)))
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
    !CALL init(phi(0),grid(1)-grid(2)+grid(1),p1=grid(1)-3*(grid(2)+grid(1)),p2=grid(1)-2*(grid(2)+grid(1)),p3=grid(1),p4=grid(2))
    !CALL init(phi(2),grid(3),p1=grid(1),p2=grid(2),p3=grid(4),RendZero=.true.)
    !CALL init(phi(0),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),RendZero=.true.)
    !CALL init(phi(1),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),LendZero=.true.)
    !CALL init(phi(0),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.)
    !CALL init(phi(2),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.)    
    CALL init(phi(0),grid(1),p3=grid(2),LendZero=.true.)
    CALL init(phi(1),grid(1),p3=grid(2),p4=grid(3),LendZero=.true.)
    CALL init(phi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=.true.)
    CALL init(phi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    CALL init(phi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
    !CALL init(phi(N+1),grid(N),p2=grid(N-1),RendZero=.true.)
    !CALL init(psi(0),grid(1)-grid(2)+grid(1),p1=grid(1)-3*(grid(2)+grid(1)),p2=grid(1)-2*(grid(2)+grid(1)),&
    !& p3=grid(1),p4=grid(2),deriv=psi_deriv)
    !CALL init(psi(0),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=psi_deriv,RendZero=.true.)
    !CALL init(psi(1),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=psi_deriv,LendZero=.true.)
    !CALL init(psi(0),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(2),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(2),grid(3),p1=grid(1),p2=grid(2),p3=grid(4),deriv=psi_deriv,RendZero=.true.)    
    !CALL init(psi(2),grid(1),p3=grid(2),p4=grid(3),deriv=psi_deriv)
    CALL init(psi(0),grid(1),p3=grid(2),deriv=psi_deriv,LendZero=.true.)
    CALL init(psi(1),grid(1),p3=grid(2),p4=grid(3),deriv=psi_deriv,LendZero=.true.)
    CALL init(psi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),deriv=psi_deriv,LendZero=.true.)
    CALL init(psi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N),deriv=psi_deriv)
    CALL init(psi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),deriv=psi_deriv,RendZero=.true.)!)
    CALL init(psi(N),grid(N),p1=grid(N-2),p2=grid(N-1),deriv=psi_deriv,RendZero=.true.)
    !CALL init(psi(N+1),grid(N),p2=grid(N-1),RendZero=.true.)
    !CALL init(chi(0),grid(1)-grid(2)+grid(1),p1=grid(1)-3*(grid(2)+grid(1)),p2=grid(1)-2*(grid(2)+grid(1)),&
    !& p3=grid(1),p4=grid(2),deriv=chi_deriv)
    !CALL init(chi(0),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=chi_deriv,RendZero=.true.)
    !CALL init(chi(1),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=chi_deriv,LendZero=.true.)
    !CALL init(chi(0),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(2),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(2),grid(3),p1=grid(1),p2=grid(2),p3=grid(4),deriv=chi_deriv,RendZero=.true.)    
    !CALL init(chi(2),grid(1),p3=grid(2),p4=grid(3),deriv=chi_deriv)
    CALL init(chi(0),grid(1),p3=grid(2),deriv=chi_deriv,LendZero=.true.)
    CALL init(chi(1),grid(1),p3=grid(2),p4=grid(3),deriv=chi_deriv,LendZero=.true.)
    CALL init(chi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),deriv=chi_deriv,LendZero=.true.)
    CALL init(chi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N),deriv=chi_deriv)
    CALL init(chi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),deriv=chi_deriv,RendZero=.true.)
    CALL init(chi(N),grid(N),p1=grid(N-2),p2=grid(N-1),deriv=chi_deriv,RendZero=.true.)
    !CALL init(chi(N+1),grid(N),p2=grid(N-1),RendZero=.true.)

  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    NN = size(phi)+size(psi)+size(chi)
    LWORK=10*NN
    LDVR = NN
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi,'G20.11)'

    A = 0.
    B = 0.
    k = 1
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
    CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    !WRITE (*,*) 'INFO =', INFO
    !WRITE (*,*) 'real(ALPHA) ='
    !WRITE (*,'(g)') ALPHAR
    !WRITE (*,*) 'BETA ='
    !WRITE (*,'(g)') BETA
    !WRITE (*,*) 'real(LAMBDA) = '
    !WRITE (*,'(g)') (ALPHAR)/BETA
    !WRITE (*,*) 's^2 = ', gamma*P(grid)/(Bz(grid))**2
    slow_inf = minval(slow_inf_range(grid))
    !WRITE (*,*) '(real(LAMBDA) - slow_inf)/slow_inf= '
    lambda = sort(ALPHAR/BETA,rev=.true.)
    IF(.not.homo) THEN
      WRITE (*,*) 'alfven range = ', alfven_range(grid)
      WRITE (*,*) 'slow_inf range = ', slow_inf_range(grid)
      WRITE (*,*) 'Lambda = '
      WRITE (*,'(g)') lambda
    ENDIF
    slow = sort(((ALPHAR)/BETA-slow_inf)/slow_inf,rev=.true.,inds=inds)
    IF (slow_evals) THEN
      IF(homo) WRITE(2,FMTR) slow(NN-nphi+1:NN)
      IF(inhomo) THEN
        WRITE (2,'(g)') eps
        WRITE (2,'(g)') (maxval(slow_inf_range(grid))-slow_inf)/slow_inf
        WRITE (2,'(i)') nphi
        WRITE (2,'(g)') slow(NN-nphi+1:NN)
        WRITE (2,*) ''
      ENDIF
    ENDIF
    WRITE (*,*) 'inds =', inds
    DO i=1,NN
      !IF (slow(i) .lt. (((minval(kz**2 *Bz0**2/rho((/0./))))*0.9-slow_inf)/slow_inf)) 

      WRITE (*,'(g,a,g)')  slow(i), ' = ', (ALPHAR(inds(i))/BETA(inds(i))-slow_inf)/slow_inf
    ENDDO
    WRITE (*,*) ' '
    !WRITE (*,*) 'k^2 Bz(grid)^2/rho(grid) = '
    !WRITE (*,'(g)') kz**2*(Bz(grid))**2/(rho(grid))
    !WRITE (*,*) 'B V1 - L1 A V1='
    pick_val = NN
    !WRITE (*,'(g,g)') matmul(D,VR(:,pick_val))-(ALPHAR(pick_val)+(0,1)*ALPHAI(pick_val))/BETA(pick_val)*matmul(C,VR(:,pick_val))
    IF(slow_evecs.and.(.not.slow_evals)) THEN
      IF(homo) THEN
        OPEN (1, status='replace',file='spline_slow_EVs.txt')
        WRITE (1,'(i)') size(xgrid)
        WRITE (1,'(i)') size(phi)
        DO i=1,size(xgrid)
          WRITE (1,'(g)') xgrid(i)
        ENDDO
        DO j=size(phi)-1,0,-1
          WRITE(1,*) ''  
          DO i=1,size(xgrid)
            WRITE (1,'(g)') sum(VR(3::3,inds(NN-j))*val(chi,xgrid(i)))
          ENDDO    
        ENDDO
        CLOSE(1)   
      ENDIF
    ENDIF
    
    IF(alfven_evecs) THEN
      ind = inds(nphi+nphi/2)!minloc(abs(ALPHAR/BETA-0.513e-04))
      WRITE (*,*) 'Lambda(',ind(1),') = ',ALPHAR(ind(1))/BETA(ind(1)) 
      WRITE (*,'(g)') (VR(2::3,ind(1)))
      WRITE (3,'(i)') size(xgrid)
      WRITE (3,'(i)') 1
      DO i=1,size(xgrid)
        WRITE (3,'(g)') xgrid(i)
      ENDDO
      WRITE(3,*) ''  
      DO i=1,size(xgrid)
        WRITE (3,'(g)') sum(VR(2::3,ind(1))*val(psi,xgrid(i)))
      ENDDO 
      WRITE (3,*) ''
      WRITE (3,'(i)') size(grid)
      WRITE (3,'(g)') 0.1D1/rho0/lambda(ind(1))/eps*sqrt(rho0*lambda(ind(1))*eps*(rho0*lambda(ind(1))-kz**2*Bz0**2))*ar

      DO i=1,size(grid)
        WRITE (3,'(g)') grid(i)
      ENDDO
      WRITE(3,*) ''  
      DO i=1,size(grid)
        WRITE (3,'(g)') sum(VR(2::3,ind(1))*val(psi,grid(i)))
      ENDDO 
      WRITE (3,*) ''
    ENDIF
    !xi_1 = sum(VR(1::3,pick_val)*val(phi,areps))
    !xi_2 = sum(VR(2::3,pick_val)*val(psi,areps))
    !xi_3 = sum(VR(3::3,pick_val)*val(chi,areps))
    !WRITE (*,*) 'xi_r(ar) = ', xi_1, 'xi_theta(ar) = ', (0,1)/mt*(xi_1-ar*xi_2), 'xi_z(ar) = ', -(0,1)*xi_3/kz
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    pick_val = NN
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    
  END SUBROUTINE bspline_deriv
  
  SUBROUTINE assemble(A,B)
    USE local
    REAL(r8), DIMENSION(:,:), INTENT(OUT) :: A,B
    A = 0.
    B = 0.
  END SUBROUTINE assemble
  
  SUBROUTINE hermite_elements()
    IMPLICIT NONE
    REAL(r8), DIMENSION(N):: grid!, slow!, slow_sort
    REAL(r8), DIMENSION(N-1) :: xi_1, xi_2, xi_3
    TYPE(bspline), DIMENSION(2:2*(N)-1) :: phi
    TYPE(bspline), DIMENSION(2:2*(N)-1) :: psi, chi 
    REAL(r8), DIMENSION(size(phi)+size(psi)+size(chi),size(phi)+size(psi)+size(chi)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,size(phi)+size(psi)+size(chi)) :: VL
    INTEGER ::  inds(size(phi)+size(psi)+size(chi))
    REAL(r8), DIMENSION(size(phi)+size(psi)+size(chi)) :: ALPHAR, ALPHAI, BETA, slow, RWORK(8*(size(phi)+size(psi)+size(chi))), WORK(10*(size(phi)+size(psi)+size(chi)))
    LOGICAL:: psi_deriv, chi_deriv
    LOGICAL :: Lend0
    REAL(r8) :: xgrid(N*30)
    
  !Initialize the finite elements
    Lend0 = .true.
    lower = (/lbound(phi,1), lbound(psi,1), lbound(chi,1)/)
    upper = (/ubound(phi,1), ubound(psi,1), ubound(chi,1)/) 
    grid(1:N) = (/ (i*ar/(N-1), i=0,N-1) /)
    !grid(1) = grid(2)/2.
    !grid(N) = (grid(1)+grid(2))/2.
    !grid = sort(grid)
    xgrid = (/ (i*ar/real(size(xgrid)-1), i=0,size(xgrid)-1) /)
    !grid(1) = epsilon(grid(2))
    nphi =  size(phi); npsi = size(psi); nchi = size(chi)
    
    !CALL init(phi(0),grid(1)-grid(2)+grid(1),p1=grid(1)-3*(grid(2)+grid(1)),p2=grid(1)-2*(grid(2)+grid(1)),p3=grid(1),p4=grid(2))
    !CALL init(phi(2),grid(3),p1=grid(1),p2=grid(2),p3=grid(4),RendZero=.true.)
    !CALL init(phi(0),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),RendZero=.true.)
    !CALL init(phi(1),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),LendZero=.true.)
    !CALL init(phi(0),grid(1),p3=grid(2),LendZero=.true.)
    !CALL init(phi(0),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.)
    !CALL init(phi(2),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.)    
    !CALL init(phi(1),grid(1),p3=grid(2),p4=grid(3),LendZero = Lend0)
    !CALL init(phi(1),grid(1))
    !CALL init(phi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero = Lend0)
    !CALL init(phi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    !CALL init(phi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi(2),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),LendZero=.true.)
    !CALL init(phi(2),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),RendZero=.true.)
    CALL init(phi(3:2*(N-1):2),grid(2:N-1), p2=grid(1:N-2),p3=grid(3:N), LendZero=.true.)
    CALL init(phi(4:2*(N-1):2),grid(2:N-1), p2=grid(1:N-2),p3=grid(3:N), RendZero=.true.)
    CALL init(phi(2*(N)-1),grid(N), p2=grid(N-1), p3=grid(N)+grid(N)-grid(N-1), RendZero=.true.)
    !CALL init(phi(2*(N)),grid(N), p2=grid(N-1), p3=grid(N)+grid(N)-grid(N-1), RendZero=.true.)
    
    
    !CALL init(phi(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
    !CALL init(phi(N+1),grid(N),p2=grid(N-1),RendZero=.true.)
    !CALL init(psi(0),grid(1)-grid(2)+grid(1),p1=grid(1)-3*(grid(2)+grid(1)),p2=grid(1)-2*(grid(2)+grid(1)),&
    !& p3=grid(1),p4=grid(2),deriv=psi_deriv)
    !CALL init(psi(0),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=psi_deriv,RendZero=.true.)
    !CALL init(psi(1),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=psi_deriv,LendZero=.true.)
    !CALL init(psi(0),grid(1),p3=grid(2),deriv=psi_deriv,LendZero=.true.)
    !CALL init(psi(0),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(2),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(2),grid(3),p1=grid(1),p2=grid(2),p3=grid(4),deriv=psi_deriv,RendZero=.true.)    
    !CALL init(psi(2),grid(1),p3=grid(2),p4=grid(3),deriv=psi_deriv)
    !CALL init(psi(1),grid(1),p3=grid(2),p4=grid(3),deriv=psi_deriv,LendZero = Lend0)
    !CALL init(psi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),deriv=psi_deriv,LendZero = Lend0)
    !CALL init(psi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N),deriv=psi_deriv)
    !CALL init(psi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),deriv=psi_deriv,RendZero=.true.)!)
    CALL init(psi(2),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),LendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(2),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),RendZero=.true.,deriv=psi_deriv)
    CALL init(psi(3:2*(N-1):2),grid(2:N-1), p2=grid(1:N-2),p3=grid(3:N), LendZero=.true.,deriv=psi_deriv)
    CALL init(psi(4:2*(N-1):2),grid(2:N-1), p2=grid(1:N-2),p3=grid(3:N), RendZero=.true.,deriv=psi_deriv)
    CALL init(psi(2*(N)-1),grid(N), p2=grid(N-1), p3=grid(N)+grid(N)-grid(N-1), RendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(2*(N)),grid(N), p2=grid(N-1), p3=grid(N)+grid(N)-grid(N-1), RendZero=.true.,deriv=psi_deriv)
    !CALL init(psi(N),grid(N),p1=grid(N-2),p2=grid(N-1),deriv=psi_deriv,RendZero=.true.)
    !CALL init(psi(N+1),grid(N),p2=grid(N-1),RendZero=.true.)
    !CALL init(chi(0),grid(1)-grid(2)+grid(1),p1=grid(1)-3*(grid(2)+grid(1)),p2=grid(1)-2*(grid(2)+grid(1)),&
    !& p3=grid(1),p4=grid(2),deriv=chi_deriv)
    !CALL init(chi(0),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=chi_deriv,RendZero=.true.)
    !CALL init(chi(1),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),deriv=chi_deriv,LendZero=.true.)
    !CALL init(chi(0),grid(1),p3=grid(2),deriv=chi_deriv,LendZero=.true.)
    !CALL init(chi(0),grid(2),p2=grid(1),p3=grid(3),RendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(2),grid(2),p2=grid(1),p3=grid(3),LendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(2),grid(3),p1=grid(1),p2=grid(2),p3=grid(4),deriv=chi_deriv,RendZero=.true.)    
    !CALL init(chi(2),grid(1),p3=grid(2),p4=grid(3),deriv=chi_deriv)
    !CALL init(chi(1),grid(1),p3=grid(2),p4=grid(3),deriv=chi_deriv,LendZero = Lend0)
    !CALL init(chi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),deriv=chi_deriv,LendZero = Lend0)
    !CALL init(chi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N),deriv=chi_deriv)
    !CALL init(chi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),deriv=chi_deriv,RendZero=.true.)
    CALL init(chi(2),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),LendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(2),grid(1),p2=grid(1)-grid(2)+grid(1),p3=grid(2),RendZero=.true.,deriv=chi_deriv)
    CALL init(chi(3:2*(N-1):2),grid(2:N-1), p2=grid(1:N-2),p3=grid(3:N), LendZero=.true.,deriv=chi_deriv)
    CALL init(chi(4:2*(N-1):2),grid(2:N-1), p2=grid(1:N-2),p3=grid(3:N), RendZero=.true.,deriv=chi_deriv)
    CALL init(chi(2*(N)-1),grid(N), p2=grid(N-1), p3=grid(N)+grid(N)-grid(N-1), RendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(2*(N)),grid(N), p2=grid(N-1), p3=grid(N)+grid(N)-grid(N-1), RendZero=.true.,deriv=chi_deriv)
    !CALL init(chi(N),grid(N),p1=grid(N-2),p2=grid(N-1),deriv=chi_deriv,RendZero=.true.)
    !CALL init(chi(N+1),grid(N),p2=grid(N-1),RendZero=.true.)

  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    NN = size(phi)+size(psi)+size(chi)
    LWORK=10*NN
    LDVR = NN
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a1,I,a)') '(',nphi,'g20.12)'

    A = 0.
    B = 0.
    k = 1
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
        !A(3*(i-lower(k))+k,3*(j-lower(l))+l) = -1./(kz*mt) * int_func(chi(j),phi(i),rho)
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
        !A(3*(i-lower(k))+k,3*(j-lower(l))+l) = 1/(kz*mt) * int_func(psi(i),chi(j),K12)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = int_func(psi(i),chi(j),W23A)
      ENDDO
    ENDDO
    k = 3
    DO i=lower(k),upper(k)
      l = 1
      DO j=max(i-chi(i)%extent(1),lower(l)),min(i+chi(i)%extent(2),upper(l))
        !A(3*(i-lower(k))+k,3*(j-lower(l))+l) = A(3*(j-lower(l))+l,3*(i-lower(k))+k)
        B(3*(i-lower(k))+k,3*(j-lower(l))+l) = B(3*(j-lower(l))+l,3*(i-lower(k))+k)
      ENDDO
      l=2
      DO j=max(i-chi(i)%extent(1),lower(l)),min(i+chi(i)%extent(2),upper(l))
        !A(3*(i-lower(k))+k,3*(j-lower(l))+l) = A(3*(j-lower(l))+l,3*(i-lower(k))+k)
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
    !B(NN,NN) = 1.
    !WRITE (*,*) 'A='
    !WRITE (*,FMT) transpose(A)
    !WRITE (*,*) 'B='
    !WRITE (*,FMT) transpose(B)
    !WRITE (*,*) 'W11B(0,1) = ', W11B((/0.,1./))
    C=A(:,:)
    D=B(:,:)
    CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    !WRITE (*,*) 'INFO =', INFO
    !WRITE (*,*) 'real(ALPHA) ='
    !WRITE (*,'(g)') ALPHAR
    !WRITE (*,*) 'BETA ='
    !WRITE (*,'(g)') BETA
    !WRITE (*,*) 'real(LAMBDA) = '
    !WRITE (*,'(g)') (ALPHAR)/BETA
    !WRITE (*,*) 's^2 = ', gamma*P(grid)/(Bz(grid))**2
    slow_inf = minval(slow_inf_range(grid))
    !WRITE (*,*) '(real(LAMBDA) - slow_inf)/slow_inf= '
    slow = sort(((ALPHAR)/BETA-slow_inf)/slow_inf,rev=.true.,inds=inds)
    WRITE (*,*) 'inds =', inds
    DO i=1,NN
      !IF (slow(i) .lt. (((minval(kz**2 *Bz0**2/rho((/0./))))*0.9-slow_inf)/slow_inf)) 

      WRITE (*,'(g,a,g)')  slow(i), ' = ', (ALPHAR(inds(i))/BETA(inds(i))-slow_inf)/slow_inf
    ENDDO
    WRITE (*,*) ' '
    !WRITE (*,*) 'k^2 Bz(grid)^2/rho(grid) = '
    !WRITE (*,'(g)') kz**2*(Bz(grid))**2/(rho(grid))
    !WRITE (*,*) 'B V1 - L1 A V1='
    pick_val = NN
    !WRITE (*,'(g,g)') matmul(D,VR(:,pick_val))-&
    !&  (ALPHAR(pick_val)+(0,1)*ALPHAI(pick_val))/BETA(pick_val)*matmul(C,VR(:,pick_val))
    IF(slow_evecs) THEN
      OPEN (1, status='replace',file='hermite_slow_EVs.txt')
      WRITE (1,'(i)') size(xgrid)
      WRITE (1,'(i)') size(phi)
      DO i=1,size(xgrid)
        WRITE (1,'(g)') xgrid(i)
      ENDDO
      DO j=size(phi)-1,0,-1
        WRITE(1,*) ''  
        DO i=1,size(xgrid)
          WRITE (1,'(g)') sum(VR(3::3,inds(NN-j))*val(chi,xgrid(i)))
        ENDDO    
      ENDDO
      CLOSE(1)
    ENDIF
    !xi_1 = sum(VR(1::3,pick_val)*val(phi,areps))
    !xi_2 = sum(VR(2::3,pick_val)*val(psi,areps))
    !xi_3 = sum(VR(3::3,pick_val)*val(chi,areps))
    !WRITE (*,*) 'xi_r(ar) = ', xi_1, 'xi_theta(ar) = ', (0,1)/mt*(xi_1-ar*xi_2), 'xi_z(ar) = ', -(0,1)*xi_3/kz
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    pick_val = NN
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    ind = minloc(abs(ALPHAR/BETA-0.513e-04))
    !WRITE (*,*) 'Lambda(',ind(1),') = ',ALPHAR(ind(1))/BETA(ind(1)) 
    !WRITE (*,'(g)') (VR(2::3,ind(1)))
  END SUBROUTINE hermite_elements
END PROGRAM cyl
