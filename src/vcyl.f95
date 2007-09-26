PROGRAM cyl
  USE local
  USE cyl_matrix_module
  USE vcyl_matrix_module
  USE finite_elements_module
  USE sort_module
  INTEGER :: ref, start, fin, N, num
  LOGICAL :: spline, phi2_deriv, phi3_deriv, vcyl
  REAL(r8) :: log_Vz0,  nq
!The following are all used in subroutines below and should not be used in the main (cyl) program
  INTEGER :: NN, ifail1, i, j, k, l, m, nphi1, nphi2, nphi3, nphi4, nphi5, nphi6, nphi, npsi, nchi,INFO
  REAL(r8) :: temp, tempA, tempB, tempC, t1, t2, st1, st2, st, at1, at2, at
  CHARACTER(LEN=30) :: FMT, FMTR
  CHARACTER(LEN=40) :: filename, date, time
  INTEGER :: LDVL=1, LWORK, LDVR, lower(9), upper(9)
  
  NAMELIST /control_params/  ref, start, fin, verbose
  NAMELIST /cyl_params/ ar, br, kz, gamma, mt, rho0, eps, Bz0, Bt0, nz, s2, equilib, N, phi2_deriv, phi3_deriv, spline, epsilo, lambd, P0, P1,  Vz0, epsVz, Vp0, epsVp
  spline = .true.
  verbose = .false.
  phi2_deriv = .true.
  phi3_deriv = .true.
  CALL cpu_time(t1)
  it = 0.
  st = 0.
  at = 0.
  st1 = 0.
  st2 = 0.
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
    CALL plot_equilibrium
    IF (spline) THEN
      CALL bspline_deriv_sa()
      CALL plot_equilibrium
      CALL bspline_deriv()   
      CALL plot_equilibrium   
    ELSE 
      CALL linear_const_sa()
      CALL plot_equilibrium  
      CALL linear_const()
      CALL plot_equilibrium  
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
      CALL bspline_deriv_sa()
      CALL plot_equilibrium  
      CALL bspline_deriv()   
      CALL plot_equilibrium     
    ELSE 
      CALL linear_const_sa()
      CALL plot_equilibrium  
      CALL linear_const()
      CALL plot_equilibrium  
    ENDIF
  ENDDO
  CALL cpu_time(t2)
  WRITE (0,*) 'Time taken: ', t2-t1, ' seconds.'
  WRITE (0,*) 'Time taken by EV solver: ',st, 'seconds.'
  WRITE (0,*) 'Time taken to assemble matrices: ',at, 'seconds.'
  WRITE (0,*) 'Time taken to integrate matrix elements: ',it, 'seconds.'
CONTAINS
SUBROUTINE linear_const_sa() !sa stands for self adjoint
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
    
    vcyl = .false.
    
  !Initialize the finite elements
    lower(1:3) = (/lbound(phi,1), lbound(psi,1), lbound(chi,1)/)
    upper(1:3) = (/ubound(phi,1), ubound(psi,1), ubound(chi,1)/)
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
  END SUBROUTINE linear_const_sa
  
  SUBROUTINE bspline_deriv_sa() !sa stands for self-adjoint
    IMPLICIT NONE
    REAL(r8), DIMENSION(N):: grid
    TYPE(bspline), DIMENSION(0:N) :: phi, psi, chi 
    REAL(r8), DIMENSION(size(phi)+size(psi)+size(chi),size(phi)+size(psi)+size(chi)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,size(phi)+size(psi)+size(chi)) :: VL
    REAL(r8), DIMENSION(size(phi)+size(psi)+size(chi)) :: ALPHAR, ALPHAI, BETA, slow,lambda,&
    & RWORK(8*(size(phi)+size(psi)+size(chi))), WORK(10*(size(phi)+size(psi)+size(chi)))
    LOGICAL :: Lend0
 
    vcyl = .false.
    
  !Initialize the finite elements
    Lend0 = .true.
    lower(1:3) = (/lbound(phi,1), lbound(psi,1), lbound(chi,1)/)
    upper(1:3) = (/ubound(phi,1), ubound(psi,1), ubound(chi,1)/) 
    grid(1:N) = (/ (i*ar/(N-1), i=0,N-1) /)
    nphi =  size(phi); npsi = size(psi); nchi = size(chi)

    CALL init(phi(0),grid(1),p3=grid(2),LendZero=.true.)
    CALL init(phi(1),grid(1),p3=grid(2),p4=grid(3),LendZero=.true.)
    CALL init(phi(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=.true.)
    CALL init(phi(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    CALL init(phi(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
  
    psi = phi
    psi%deriv = phi2_deriv
    chi=psi
    chi%deriv = phi3_deriv

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
  END SUBROUTINE bspline_deriv_sa
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
    
    vcyl = .true.
    
  !Initialize the finite elements
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3)

    CALL init(phi1(1),grid(1),p3=grid(2))
    CALL init(phi1(2:N-1),grid(2:N-1),p2=grid(1:N-2),p3=grid(3:N))
    CALL init(phi2(1:N-1),grid(1:N-1),p3=grid(2:N))
    CALL init(phi3(1:N-1),grid(1:N-1),p3=grid(2:N))
    WRITE (*,*) 'The linear finite elements are not implemented.'
  END SUBROUTINE linear_const
  
  
  SUBROUTINE bspline_deriv()
    IMPLICIT NONE
    REAL(r8), DIMENSION(N):: grid
    TYPE(bspline), DIMENSION(0:N) :: phi1, phi2, phi3, phi4, phi5, phi6
    REAL(r8), DIMENSION(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6),&
    & size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)):: A, B, C, D, VR
    REAL(r8), DIMENSION(1,size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)) :: VL
    INTEGER ::  inds(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6))
    REAL(r8), DIMENSION(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)) :: ALPHAR, ALPHAI, BETA, slow,lambda,&
    & RWORK(8*(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6))), &
    & WORK(10*(size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)))
    LOGICAL :: Lend0
 
    vcyl = .true.
    
  !Initialize the finite elements
    Lend0 = .true.
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
    grid(1:N) = (/ (i*ar/(N-1), i=0,N-1) /)
    
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
        temp = minval(tempA+tempB+tempC-&
        & 1./2.*( 2*Bz((/ar/))*Bt((/ar/))*kz/mt*val(phi1(i),ar)*val(phi1(j),ar) + &
        &         Bmag((/ar/))**2*(val(phi1(i),ar)*val_prime(phi1(j),ar)+val(phi1(j),ar)*val_prime(phi1(i),ar))))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = temp
        IF((i.ne.j).and.(temp.ne.0)) B(6*(j-lower(m))+k,6*(i-lower(l))+l) = temp
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi1(i),phi2(j),B42b,deriv1=.true.)
        tempC = int_func(phi1(i),phi2(j),B42c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = &
        & minval(tempB+tempC-1./2.*val(phi1(i),ar)*val(phi2(j),ar)*(Bz((/ar/))**2-ar*Bz((/ar/))*Bt((/ar/))*kz/mt))
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi1(i),phi3(j),B43b,deriv1=.true.)
        tempC = int_func(phi1(i),phi3(j),B43c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC - &
        & minval(1./2.*val(phi1(i),ar)*val(phi3(j),ar)*(Bt((/ar/))**2-Bz((/ar/))*Bt((/ar/))*mt/(kz*ar)))
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
    CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs(phi1,phi2,phi3,phi4,phi5,phi6,grid,VR) 
  END SUBROUTINE bspline_deriv  
  SUBROUTINE plot_equilibrium
    INTEGER num_pts
    INTEGER :: i
    REAL(r8), DIMENSION(N) :: grid, pressure, zb, tb, density, pv, zv, p_prime, Bt_prime, Bz_prime, safety
    num_pts = N
    grid = (/ (i*ar/(num_pts-1), i=0,num_pts-1) /)
    pressure = P(grid)
    tb = Bt(grid)
    zb = Bz(grid)
    density = rho(grid)
    pv = Vp(grid)*grid
    zv = Vz(grid)
    safety = q(grid)
    p_prime = diff(P,grid)
    Bz_prime = diff(Bz,grid)
    Bt_prime = diff(Bt,grid)
    IF(vcyl) THEN
      WRITE(filename,'(a,i0,a)') 'equilibria_vcyl/',num,'.txt'
    ELSE
      pv = 0.
      zv = 0.
      WRITE(filename,'(a,i0,a)') 'equilibria_cyl/',num,'.txt'
    ENDIF
    WRITE (*,'(a,i0,a,a)') 'Writing equilibrium data for run number ',num, ' to ',filename
    OPEN(1,file=filename,status='replace')
    DO i=1,num_pts
      WRITE (1,'(g,g,g,g,g,g,g,g)') grid(i), pressure(i),zb(i),tb(i),density(i),pv(i),zv(i), safety(i)
    ENDDO
    CLOSE(1)
    IF (verbose) THEN
      WRITE (*,*) 'Eqilibrium condition at grid points='
      WRITE (*,'(g)') equilibrium_v(grid(2:))
      WRITE (*,'(5(a20))') 'pressure','p_prime', 'Bz*Bz_prime', 'Bt*Bt_prime', 'Bt**2/r'
      DO i=2,size(grid)
        WRITE (*,'(5(g20.10))') pressure(i), p_prime(i), zb(i)*Bz_prime(i), tb(i)*Bt_prime(i), tb(i)**2/grid(i)
      ENDDO
    ENDIF
  END SUBROUTINE plot_equilibrium

  SUBROUTINE output_params(assemble_t,solve_t)
    REAL(r8), INTENT(IN) :: assemble_t,solve_t
    IF(vcyl) THEN
      WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.dat'
    ELSE
      WRITE(filename,'((a,i0),a)') 'output_cyl/', num,'.dat'
    ENDIF
    OPEN(1,file=filename,status='replace')
    IF (spline) THEN
      WRITE(1,'(a)') 'spline'
    ELSE 
      WRITE(1,'(a)') 'linconst'
    ENDIF
    WRITE(1,'(i)') N, NN, mt, equilib, num, nz
    WRITE(1,'(g)') epsilo, assemble_t, solve_t, kz, ar, rho0, Bz0, Bt0, s2, eps, P0, P1, lambd, Vz0, epsVz, Vp0, epsVp
    CLOSE(1)
  END SUBROUTINE output_params
  SUBROUTINE output_evals(evalsr,evalsi)
    REAL(r8), DIMENSION(:), INTENT(IN) :: evalsr, evalsi
    IF(vcyl) THEN
      WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.evalsr'
    ELSE
      WRITE(filename,'((a,i0),a)') 'output_cyl/', num,'.evalsr'
    ENDIF
    OPEN(1,file=filename,status='replace')
    WRITE(1,'(g)') evalsr
    CLOSE(1)
    IF(vcyl) THEN
      WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.evalsi'
    ELSE
      WRITE(filename,'((a,i0),a)') 'output_cyl/', num,'.evalsi'
    ENDIF
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
  SUBROUTINE output_evecs(phi1,phi2,phi3,phi4,phi5,phi6,grid,VR)
    TYPE(bspline), DIMENSION(:), INTENT(IN) :: phi1, phi2, phi3,phi4,phi5,phi6
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: VR
    TYPE(bspline), DIMENSION(size(phi1),6) :: phi
    INTEGER :: N , i, j, k
    CHARACTER(LEN=30) FMT
    phi = reshape((/phi1,phi2,phi3,phi4,phi5,phi6/),(/size(phi1),6/))
    N = size(grid)
    WRITE(FMT,'(a,i0,a)') '(',N,'(g,a))'
    WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.grid'
    OPEN(1,file=filename,status='replace')
    WRITE(1,FMT) (grid(i),',',i=1,size(grid))
    CLOSE(1)
    DO j=1,6
      WRITE(filename,'(2(a,i0))') 'output_vcyl/', num,'.evecs',j
      OPEN(1,file=filename,status='replace')
      DO k=1,size(VR(1,:))
        WRITE(1,FMT) ( sum(VR(j::6,k)*val(phi(:,j),grid(i))),',' , i=1,size(grid) )
      ENDDO
      CLOSE(1)
    ENDDO 
  END SUBROUTINE output_evecs
END PROGRAM cyl
