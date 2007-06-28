PROGRAM cyl
  USE local
  USE vcyl_matrix_module
  USE finite_elements_module
  USE sort_module
  INTEGER :: N, min_N, max_N, delta_N, ephi2, max_Vz0, delta_Vz0, max_Bt0, delta_Bt0, Bt_ind, k_ind, delta_k, max_k
  LOGICAL :: linconst, spline, hermite, slow_evals, slow_evecs, alfven_evecs, phi2_deriv, phi3_deriv, homo_plasma, inhomo, axial_flow, azi_flow, slow_evecs1, alfven_evecs1, slow_evals1, var_Bt0, var_k, converge
!The following are all used in subroutines below and should not be used in the main (cyl) program
  INTEGER :: i, j, k, l, m, nphi1, nphi2, nphi3, nphi4, nphi5, nphi6, INFO, pick_val, ind4(1), ind8(1), stat, ind(1)
  REAL(r8) :: temp, tempA, tempB, tempC, areps, slow_inf, t1, t2, log_Bt0, TT, log_k, st1, st2, st, at1, at2, at
  CHARACTER(LEN=30) :: FMT, FMTR
  CHARACTER(LEN=40) :: filename
  INTEGER :: NN, ifail
  INTEGER :: LDVL=1, LWORK, LDVR, lower(9), upper(9)
  
  NAMELIST /control_params/ min_N, max_N, delta_N, linconst, spline, hermite, verbose, slow_evals, slow_evecs, alfven_evecs,phi2_deriv, phi3_deriv, homo_plasma, inhomo, axial_flow, max_Vz0, delta_Vz0, &
  & max_Bt0, delta_Bt0, azi_flow, log_Bt0, epsilo, converge
  NAMELIST /cyl_params/ ar, br, kz, gamma, mt, rho0, eps, homo, Bz0, Bt0, Vz0, epsVz, Vp0, epsVp, s2
  NAMELIST /azi_flow_params/ Bt0, max_Bt0, delta_Bt0, log_Bt0, max_k, delta_k, log_k, var_Bt0, var_k, kz
  CALL cpu_time(t1)
  converge = .false.
  it = 0.
  st = 0.
  at = 0.
  st1 = 0.
  st2 = 0.
  min_N = 5
  max_N = 6
  linconst = .true.
  spline = .true.
  hermite = .true.
  verbose = .false.
  slow_evals = .true.
  slow_evecs = .false.
  phi2_deriv = .true.
  phi3_deriv = .true.
  homo_plasma = .true.
  inhomo = .false.
  axial_flow = .false.
  max_Bt0 = 0
  delta_Bt0 = 0
  Bt0 = 0
  s2 = 1./12.
  ifail = 1
  TT = s17aef(1.2D0,ifail)
  epsilo = 1.D-10
  WRITE (*,*) 'bessel_j0(1.2) = ', TT, ifail
  OPEN(1,file='src/vcontrol_params.in',status='old',form='formatted')
  READ(1,nml=control_params)
  CLOSE(1)
  WRITE(*,nml=control_params) 
  slow_evecs1 = slow_evecs
  alfven_evecs1 = alfven_evecs
  slow_evals1 = slow_evals
  IF(spline.and.converge) THEN
    slow_evals = .false.
    slow_evecs = .false.
    alfven_evecs = .false.
    OPEN(1,file='src/converge.in',status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    WRITE(filename,'(a,i0,a,i0,a,i0,a)') 'spline_converge_',min_N,'_',max_N,'_', delta_N,'.txt'
    OPEN(10,file=filename,status='replace')
    WRITE (0,*) 'Output written to ', filename
    DO N=min_N, max_N, delta_N
      CALL bspline_deriv()
    ENDDO
    CLOSE(10)
  ENDIF
  IF(homo_plasma) THEN
    OPEN(1,file='src/vhomogeneous.in',status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    WRITE (*,nml=cyl_params)
    IF (linconst.and.slow_evals1) THEN
      slow_evecs=.false.
      alfven_evecs=.false.
      WRITE(filename,'(a,i0,a,i0,a)') 'lin_slow_evals_',min_N,'_',max_N,'.txt'
      WRITE(0,*) 'Output written to ', filename
      OPEN(2,file=filename,status='replace')
      DO N = min_N,max_N
        WRITE (*,*) 'N = ', N
        WRITE (*,*) 'With linear and constant elements.'
        CALL linear_const()
      ENDDO
      CLOSE(2)
    ENDIF
    IF (linconst.and.slow_evecs1) THEN
      slow_evecs = .true.
      slow_evals = .false.
      CALL linear_const()
    ENDIF
    IF(slow_evals1.and.spline) THEN
      slow_evals = .true.
      slow_evecs = .false.
      WRITE (filename,'(a,i0,a,i0,a)') 'spline_slow_evals_',min_N,'_',max_N,'.txt'
      OPEN(2,file=filename,status='replace')
      WRITE (0,*) 'Output written to ', filename
      DO N = min_N,max_N
        WRITE (*,*) 'N = ', N
        WRITE (*,*) 'With bspline elements.'
        CALL bspline_deriv()
      ENDDO
      CLOSE(2)
    ENDIF
    IF (spline.and.slow_evecs1) THEN
      slow_evals = .false.
      slow_evecs = .true.
      N = max_N
      CALL bspline_deriv()
    ENDIF
    IF (axial_flow) THEN
      slow_evals = .false.
      slow_evecs = .false.
      alfven_evecs = .false.
      N = min_N
      WRITE(filename,'(a,i0,a,i0,a,i0,a)') 'spline',N,'_var_Vz0_',delta_Vz0,'_',max_Vz0,'.txt'
      OPEN(4,file=filename,status='replace')
      WRITE (0,*) 'Output written to ', filename
      DO Vz0 = 0,max_Vz0, delta_Vz0
        IF (spline) THEN
          WRITE(*,'(a,i,a,g)') 'With ',N,'bspline elements    Vz0 = ', Vz0
          CALL bspline_deriv()
        ENDIF
      ENDDO
      CLOSE(4)
    ENDIF
    IF (azi_flow) THEN
      slow_evals = .false.
      slow_evecs = .false.
      alfven_evecs = .false.
      N = min_N
      Vz0 = 0
      OPEN(1,file='src/azi_flow.in',status='old',form='formatted')
      READ(1,nml=azi_flow_params)
      CLOSE(1)
      IF(var_Bt0) THEN
        WRITE(filename,'(a,i0,a,i0,a,i0,a,f0.1,a)') 'spline',N,'_var_Bt0_',delta_Bt0,'_',max_Bt0,'_',log_Bt0,'.txt'

        WRITE (0,*) 'Output written to ', filename
        OPEN(7,file=filename,status='replace')
        DO Bt_ind = 0, max_Bt0, delta_Bt0
          Bt0 = Bt_ind*10**log_Bt0
          IF (spline) THEN
            WRITE(*,'(a,i3,a)') 'With ',N,' bspline elements'
            CALL bspline_deriv()
          ENDIF
        ENDDO
        CLOSE(7)
      ENDIF
      OPEN(1,file='src/azi_flow.in',status='old',form='formatted')
      READ(1,nml=azi_flow_params)
      CLOSE(1)
      IF(var_k) THEN
        WRITE(filename,'(a,i0,a,i0,a,i0,a,f0.1,a)') 'spline',N,'_var_k_',delta_k,'_',max_k,'_',log_k,'.txt'
        WRITE (0,*) 'Output written to ', filename
        OPEN(7,file=filename,status='replace')
        DO k_ind = delta_k, max_k, delta_k
          kz = k_ind*10**log_k
          IF (spline) THEN
            WRITE(*,'(a,i3,a)') 'With ',N,' bspline elements'
            CALL bspline_deriv()
          ENDIF
        ENDDO
        CLOSE(7)
      ENDIF
    ENDIF
  ENDIF
  IF(inhomo) THEN
    OPEN(1,file='src/vinhomogeneous.in',status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    WRITE (*,nml=cyl_params)
    N = min_N
    IF (alfven_evecs1) THEN
      alfven_evecs = .true.
      slow_evecs = .false.
      slow_evals = .false.
      IF(linconst)  THEN
        WRITE (filename,'(a,3(i0,a),es6.0,a)') 'lin',min_N,'_',max_N,'_',delta_N,'_alfven_EVs_rho_eps',eps,'.txt'
        WRITE (*,*) 'Output written to ', filename
        OPEN (3, status='replace', file=filename)
        DO N=min_N,max_N,delta_N
          WRITE (*,*) 'N = ', N
          WRITE (*,*) 'With linear and constant elements.'
          CALL linear_const()
        ENDDO
        CLOSE(3)
      ENDIF
      IF(spline) THEN
        WRITE (filename,'(a,3(i0,a),es6.0,a)') 'spline',min_N,'_',max_N,'_',delta_N,'_alfven_EVs_rho_eps',eps,'.txt'
        WRITE (*,*) 'Output written to ', filename
        OPEN (3, status='replace',file=filename)
        DO N=min_N,max_N,delta_N
          WRITE (*,*) 'N = ', N
          WRITE (*,*) 'With bspline elements.'
          CALL bspline_deriv()
        ENDDO
        CLOSE(3)
      ENDIF
    ENDIF
    IF (slow_evals1) THEN
      slow_evals = .true.
      slow_evecs = .false.
      alfven_evecs = .false.
      N = 20
      IF(linconst) THEN
      ENDIF
      IF(spline) THEN
        WRITE (filename,'(a,i0,a,i0,a,i0,a)') 'spline',N,'_slow_vareps_',-10,'_',-4,'.txt'
        OPEN (2, status='replace',file=filename)
        WRITE (*,*) 'Output written to ', filename
        eps = 0
        CALL bspline_deriv()
        DO ephi2 = -10,-4
          eps = 10**(ephi2/2.0)
          CALL bspline_deriv()
        ENDDO
        CLOSE(3)
      ENDIF
    ENDIF
  ENDIF
  IF (hermite) THEN
    WRITE (*,*) 'With Hermite elements.'
    !CALL hermite_elements()
  ENDIF
  WRITE (*,*) 'Alfven range (approx):', alfven_range((/0.,ar/))
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
    TYPE(linear), DIMENSION(N-1) :: phi1 , phi4
    TYPE(constant), DIMENSION(N-1) :: phi2, phi3, phi5, phi6
    REAL(r8), DIMENSION(3*(N-1),3*(N-1)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,3*(N-1)) :: VL
    INTEGER inds(3*size(phi1))
    REAL(r8), DIMENSION(3*(N-1)) :: ALPHAR, ALPHAI, BETA, RWORK(8*3*(N-1)), WORK(10*3*(N-1)), slow, lambda
    REAL(r8) :: xgrid(N*30)
    
  !Initialize the finite elements
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3)

    CALL init(phi1(1),grid(1),p3=grid(2))
    CALL init(phi1(2:N-1),grid(2:N-1),p2=grid(1:N-2),p3=grid(3:N))
    CALL init(phi2(1:N-1),grid(1:N-1),p3=grid(2:N))
    CALL init(phi3(1:N-1),grid(1:N-1),p3=grid(2:N))


  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    NN = size(phi1)+size(phi2)+size(phi3)
    LWORK=10*NN
    LDVR = NN
    WRITE(FMT,'(a1,I,a)') '(',NN,'g14.5)'
    WRITE(FMTR,'(a1,I,a)') '(',nphi1,'g20.12)'

    A = 0.
    B = 0.
    k = 1

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
    IF (slow_evals) WRITE(2,FMTR) slow(NN-nphi1+1:NN)
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
      WRITE (1,'(i)') size(phi1)
      DO i=1,size(grid)-1
        WRITE (1,'(g)') (grid(i)+grid(i+1))/2.
      ENDDO
      DO j=size(phi1)-1,0,-1
        WRITE(1,*) ''  
        DO i=1,size(grid)-1
          WRITE (1,'(g)') sum(VR(3::3,inds(NN-j))*val(phi3,(grid(i)+grid(i+1))/2))
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
    TYPE(bspline), DIMENSION(0:N) :: phi1, phi2, phi3, phi4, phi5, phi6
    REAL(r8), DIMENSION(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6),&
    & size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)) :: VL
    INTEGER ::  inds(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6))
    REAL(r8), DIMENSION(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)) :: ALPHAR, ALPHAI, BETA, slow,lambda,&
    & RWORK(8*(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6))), &
    & WORK(10*(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)))
    LOGICAL :: Lend0
    REAL(r8) :: xgrid(N*30)
 
  !Initialize the finite elements
    Lend0 = .true.
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
    grid(1:N) = (/ (i*ar/(N-1), i=0,N-1) /)
    IF(verbose) THEN
      WRITE (*,*) 'Eqilibrium condition at grid points='
      WRITE (*,'(g)') equilibrium(grid(2:))
    ENDIF
    xgrid = (/ (i*ar/real(size(xgrid)-1), i=0,size(xgrid)-1) /)
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3); nphi4 = size(phi4); nphi5 = size(phi5); nphi6 = size(phi6)
    IF(lower(1).eq.0)   CALL init(phi1(0),grid(1),p3=grid(2),LendZero=.true.)
    IF(lower(1).le.1)   CALL init(phi1(1),grid(1),p3=grid(2),p4=grid(3),LendZero=.true.)
    IF(lower(1).le.2)   CALL init(phi1(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=.true.)
    CALL init(phi1(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    CALL init(phi1(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi1(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
    phi2 = phi1
    phi3 = phi1
    phi4 = phi1
    phi5 = phi1
    phi6 = phi1
    phi2%deriv = phi2_deriv
    phi3%deriv = phi3_deriv
    phi5%deriv = phi2_deriv
    phi6%deriv = phi3_deriv
    
  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)
    LWORK=10*NN
    LDVR = NN
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi1,'G20.11)'

    A = 0.
    B = 0.
    k = 1
    m = k+3
    CALL cpu_time(at1)
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l)  = int_func(phi4(i),phi1(j),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),B11)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l)  = int_func(phi4(i),phi2(j),A12)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),B12)
      ENDDO
      l = 4
      DO j=i,min(i+3,upper(l))
        temp = int_func(phi4(j),phi4(i),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = temp
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = temp
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi5(j),phi4(i),A12)
      ENDDO
    ENDDO
    k = 2
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi5(i),phi1(j),A12)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi5(i),phi1(j),B12)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi5(i),phi2(j),A22)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi5(i),phi2(j),B22)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+1,6*(i-lower(m))+5)
      ENDDO
      l = 5
      DO j=i,min(i+3,upper(l))
        temp = int_func(phi5(i),phi5(j),A22)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = temp
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = temp
      ENDDO
    ENDDO
    k = 3
    m = k+3
    DO i=lower(m),upper(m)
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l)  = int_func(phi6(i),phi3(j),A33)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),B33)
      ENDDO
      l = 6
      DO j=i,min(i+3,upper(l))
        temp = int_func(phi6(i),phi6(j),A33)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = temp
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = temp
      ENDDO
    ENDDO
    k = 4
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=i,min(i+3,upper(l))
        tempA = int_func(phi1(i),phi1(j),B41a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi1(i),phi1(j),B41b,deriv2=.true.)+int_func(phi1(i),phi1(j),B41b,deriv1=.true.)
        tempC = int_func(phi1(i),phi1(j),B41c)
        temp = tempA+tempB+tempC-&
        & 1./2.*( 2*Bz(ar)*Bz(ar)*kz/mt*val(phi1(i),ar)*val(phi1(j),ar) + &
        &         Bmag(ar)**2*(val(phi1(i),ar)*val_prime(phi1(j),ar)+val(phi1(j),ar)*val_prime(phi1(i),ar)))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = temp
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = temp
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi1(i),phi2(j),B42b,deriv1=.true.)
        tempC = int_func(phi1(i),phi2(j),B42c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = &
        & tempB+tempC-1./2.*val(phi1(i),ar)*val(phi2(j),ar)*(Bz(ar)**2-ar*Bz(ar)*Bt(ar)*kz/mt)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi1(i),phi3(j),B43b,deriv1=.true.)
        tempC = int_func(phi1(i),phi3(j),B43c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC - &
        & 1./2.*val(phi1(i),ar)*val(phi3(j),ar)*(Bt(ar)**2-Bz(ar)*Bt(ar)*mt/(kz*ar))
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l)  = int_func(phi1(i),phi4(j),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+m,6*(i-lower(m))+m)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(l))+2,6*(i-lower(m))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+2,6*(i-lower(m))+1)
      ENDDO
    ENDDO
    k = 5
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+4,6*(i-lower(m))+2)
      ENDDO
      l = 2
      DO j=i,min(i+3,upper(l))
        tempB = int_func(phi2(i),phi2(j),B52)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = tempB
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(i),phi3(j),B53)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(l))+1,6*(i-lower(m))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+1,6*(i-lower(m))+2)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(l))+2,6*(i-lower(m))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+2,6*(i-lower(m))+2)
      ENDDO
    ENDDO
    k = 6
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+4,6*(i-lower(m))+3)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+5,6*(i-lower(m))+3)
      ENDDO
      l = 3
      DO j=i,min(i+3,upper(l))
        tempB = int_func(phi3(i),phi3(j),B63)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = tempB
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(l))+3,6*(i-lower(m))+3)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+3,6*(i-lower(m))+3)
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
    IF(slow_evecs.or.alfven_evecs) THEN
      CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ELSE
      CALL DGGEV('N','N',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    !WRITE (*,*) 'INFO =', INFO
    !WRITE (*,*) 'real(ALPHA) ='
    !WRITE (*,'(g)') ALPHAR
    !WRITE (*,*) 'BETA ='
    !WRITE (*,'(g)') BETA
    IF(axial_flow) THEN
      WRITE (4,'(g)') Vz0
      WRITE (4,'(g)') kz
      WRITE (4,'(i)') NN
      WRITE (4,'(g)') (ALPHAR/BETA)
      WRITE (4,'(g)') (ALPHAI/BETA)
      WRITE (4,*) ''
    ENDIF
    IF(azi_flow) THEN
      WRITE (7,'(g)') Bt0
      WRITE (7,'(g)') rho0
      WRITE (7,'(g)') kz
      WRITE (7,'(i)') mt
      WRITE (7,'(g)') ar
      WRITE (7,'(i)') NN
      WRITE (7,'(g)') (ALPHAR/BETA)
      WRITE (7,'(g)') (ALPHAI/BETA)
      WRITE (7,*) ''
    ENDIF
    slow_inf = minval(slow_inf_range(grid))
    IF(converge) THEN
      WRITE (10,'(i)') N
      WRITE (10,'(g)') sqrt(alfven_range((/0.,ar/)))
      WRITE (10,'(g)') sqrt(max_slow())
      WRITE (10,'(g)') sqrt(slow_inf_range((/0.,ar/)))
      WRITE (10,'(i)') NN
      WRITE (10,'(g)') (ALPHAR/BETA)
      WRITE (10,'(g)') (ALPHAI/BETA)
      WRITE (10,*) ''
    ENDIF
    !WRITE (*,*) 's^2 = ', gamma*P(grid)/(Bz(grid))**2
    !WRITE (*,*) '(real(LAMBDA) - slow_inf)/slow_inf= '
    lambda = ALPHAR**2/BETA**2
    IF(.not.homo) THEN
      WRITE (*,*) 'alfven range = ', alfven_range(grid)
      WRITE (*,*) 'slow_inf range = ', slow_inf_range(grid)
      WRITE (*,*) 'Lambda = '
      WRITE (*,'(g)') sort(lambda,rev=.true.)
    ENDIF
    slow = sort(((ALPHAR)**2/BETA**2-slow_inf)/slow_inf,rev=.true.,inds=inds)
    IF (slow_evals) THEN
      IF(homo) WRITE(2,FMTR) slow(NN-2*nphi1+1:NN:2)
      IF(inhomo) THEN
        WRITE (2,'(g)') eps
        WRITE (2,'(g)') (maxval(slow_inf_range(grid))-slow_inf)/slow_inf
        WRITE (2,'(i)') nphi1
        WRITE (2,'(g)') slow(NN-2*nphi1+1:NN:2)
        WRITE (2,*) ''
      ENDIF
    ENDIF
    WRITE (*,*) 'inds =', inds
    WRITE (*,*) 'max_slow = ', max_slow(), '(max_slow-slow_inf)/slow_inf = ', (max_slow()-slow_inf)/slow_inf
    DO i=1,NN
      !IF (slow(i) .lt. (((minval(kz**2 *Bz0**2/rho((/0./))))*0.9-slow_inf)/slow_inf)) 

      WRITE (*,'(g,a,g)')  slow(i), ' = ', (ALPHAR(inds(i))**2/BETA(inds(i))**2-slow_inf)/slow_inf
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
        WRITE (1,'(i)') size(phi1)
        DO i=1,size(xgrid)
          WRITE (1,'(g)') xgrid(i)
        ENDDO
        DO j=2*size(phi1)-1,0,-2
          WRITE(1,*) ''  
          IF(ALPHAI(inds(NN-j)).ne.0.) THEN
            WRITE (0,*) 'Warning: Eigenmode for complex eigenvalue not computed correctly'
          ENDIF
          DO i=1,size(xgrid)
            WRITE (1,'(g)') sum(VR(3::6,inds(NN-j))*val(phi3,xgrid(i)))
          ENDDO    
        ENDDO
        CLOSE(1)   
      ENDIF
    ENDIF
    
    IF(alfven_evecs) THEN
      !WRITE (*,'(a20,a20)') 'Lambda','Singular r'
      !DO i=1,2*nphi1
        !WRITE (*,'(g20,g20)') lambda(2*(nphi1+1)+i)
      !ENDDO
      !WRITE (*,
      ind = inds(2*(nphi1+1)+nphi1+1)!minloc(abs(ALPHAR/BETA-0.513e-04))
      WRITE (*,*) 'Lambda(',ind(1),') = ',ALPHAR(ind(1))**2/BETA(ind(1))**2
      WRITE (*,'(g)') (VR(2::6,ind(1)))
      WRITE (3,'(i)') size(xgrid)
      WRITE (3,'(i)') 1
      DO i=1,size(xgrid)
        WRITE (3,'(g)') xgrid(i)
      ENDDO
      WRITE(3,*) ''  
      DO i=1,size(xgrid)
        WRITE (3,'(g)') sum(VR(2::6,ind(1))*val(phi2,xgrid(i)))
      ENDDO 
      WRITE (3,*) ''
      WRITE (3,'(i)') size(grid)
      WRITE (3,'(g)') 0.1D1/rho0/lambda(ind(1))/eps*sqrt(rho0*lambda(ind(1))*eps*(rho0*lambda(ind(1))-kz**2*Bz0**2))*ar

      DO i=1,size(grid)
        WRITE (3,'(g)') grid(i)
      ENDDO
      WRITE(3,*) ''  
      DO i=1,size(grid)
        WRITE (3,'(g)') sum(VR(2::6,ind(1))*val(phi2,grid(i)))
      ENDDO 
      WRITE (3,*) ''
      
    ENDIF
    !WRITE (*,*) 'xi_r(ar) = ', xi_1, 'xi_theta(ar) = ', (0,1)/mt*(xi_1-ar*xi_2), 'xi_z(ar) = ', -(0,1)*xi_3/kz
    !WRITE (*,*) 'Lambda(',pick_val,') - slow_inf = ', ALPHAR(pick_val)/BETA(pick_val) - slow_inf
    !WRITE (*,'(g)') VR(3::3,pick_val)
    
  END SUBROUTINE bspline_deriv  
END PROGRAM cyl
