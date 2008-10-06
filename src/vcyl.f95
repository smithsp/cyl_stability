PROGRAM cyl
  USE local
  USE cyl_matrix_module
  USE vcyl_matrix_module
  USE finite_elements_module
  USE sort_module
  INTEGER :: ref, start, fin, N, num, BCrow
  LOGICAL :: spline, phi2_deriv, phi3_deriv, vcyl
  REAL(r8) :: log_Vz0,  nq
  CHARACTER(LEN=40) :: filename, date, time
!The following are all used in subroutines below and should not be used in the main (cyl) program
  INTEGER :: NN, ifail1, i, j, k, l, m, nphi1, nphi2, nphi3, nphi4, nphi5, nphi6, nphi, npsi, nchi,INFO
  REAL(r8) :: temp, tempA, tempB, tempC, t1, t2, st1, st2, st, at1, at2, at
  REAL(r8), DIMENSION(:), ALLOCATABLE :: grid, ALPHAR, ALPHAI, BETA, RWORK, WORK
  REAL(r8), DIMENSION(:,:), ALLOCATABLE :: A, B, C, D, VR
  REAL(r8), DIMENSION(:,:), ALLOCATABLE :: VL
  CHARACTER(LEN=30) :: FMT, FMTR
  INTEGER :: LDVL=1, LWORK, LDVR, lower(9), upper(9)
  LOGICAL :: Lend0, Evecs
  
  NAMELIST /control_params/  ref, start, fin, verbose
  NAMELIST /cyl_params/ equilib, N, BCrow, kz, mt, ar, br, tw, gamma, Lend0, Evecs, spline, alpha, rs, epsilo, rho0, eps, P0, P1, s2, Bz0, Bt0, lambd, Vz0, epsVz, Vp0, epsVp
  spline = .false.
  verbose = .false.
  phi2_deriv = .true.
  phi3_deriv = .true.
  Evecs = .false.
  CALL cpu_time(t1)
  it = 0.;  st = 0.;  at = 0.;  st1 = 0.;  st2 = 0.
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
    IF (br.ge.(1.e4)) THEN
      tw=0
      br=1.e4
    ENDIF
    OPEN(1,file=filename,status='replace',form='formatted')
    WRITE(1,nml=cyl_params)
    CLOSE(1)
    WRITE(*,nml=cyl_params)
    IF (spline) THEN
      IF((Vz0.eq.0).and.(Vp0.eq.0).and.(epsVp.eq.0).and.(epsVz.eq.0).and.(br.le.ar)) THEN
        CALL bspline_deriv_sa
      ELSEIF ((tw.gt.0).and.(br.gt.ar)) THEN
        CALL bspline_deriv_rw
      ELSE
        CALL bspline_deriv
      ENDIF
    ELSE 
      IF((Vz0.eq.0).and.(Vp0.eq.0).and.(epsVp.eq.0).and.(epsVz.eq.0).and.(br.le.ar)) THEN
        CALL linear_const_sa
      ELSEIF ((tw.gt.0).and.(br.gt.ar)) THEN
        CALL linear_const_rw
      ELSE
        CALL linear_const
      ENDIF
    ENDIF
    CALL plot_equilibrium  
  ENDIF
  DO num = start,fin
    WRITE(filename,'(a,i0,a)') 'input/',num,'.in'
    WRITE(*,*)  'Inputting parameters from ', filename
    OPEN(1,file=filename,status='old',form='formatted')
    READ(1,nml=cyl_params)
    CLOSE(1)
    IF (br.ge.(1.e4)) THEN
      tw=0
      br=1.e4
    ENDIF
    OPEN(1,file=filename,status='replace',form='formatted')
    WRITE(1,nml=cyl_params)
    CLOSE(1)
    WRITE(*,nml=cyl_params)
    IF (spline) THEN
      IF((Vz0.eq.0).and.(Vp0.eq.0).and.(epsVp.eq.0).and.(epsVz.eq.0).and.(br.le.ar)) THEN
        CALL bspline_deriv_sa
      ELSEIF ((tw.gt.0).and.(br.gt.ar)) THEN
        CALL bspline_deriv_rw
      ELSE
        CALL bspline_deriv
      ENDIF
    ELSE 
      IF((Vz0.eq.0).and.(Vp0.eq.0).and.(epsVp.eq.0).and.(epsVz.eq.0).and.(br.le.ar)) THEN
        CALL linear_const_sa
      ELSEIF ((tw.gt.0).and.(br.gt.ar)) THEN
        CALL linear_const_rw
      ELSE
        CALL linear_const
      ENDIF
    ENDIF
    CALL plot_equilibrium  
  ENDDO
  CALL cpu_time(t2)
  WRITE (0,*) 'Time taken: ', t2-t1, ' seconds.'
  WRITE (0,*) 'Time taken by EV solver: ',st, 'seconds.'
  WRITE (0,*) 'Time taken to assemble matrices: ',at, 'seconds.'
  WRITE (0,*) 'Time taken to integrate matrix elements: ',it, 'seconds.'
CONTAINS
SUBROUTINE linear_const_sa() !sa stands for self adjoint
    IMPLICIT NONE
    TYPE(linear), DIMENSION(:), ALLOCATABLE :: phi 
    TYPE(constant), DIMENSION(:), ALLOCATABLE :: psi, chi
    
    ALLOCATE(grid(N),psi(N-1),chi(N-1))
    IF ((abs(mt).eq.1).and.(.not.Lend0)) THEN
      ALLOCATE(phi(N-1))
    ELSE
      ALLOCATE(phi(2:N-1))
    ENDIF
    NN = size(phi)+size(psi)+size(chi)
    ALLOCATE(A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),VR(NN,NN),VL(1,NN),ALPHAI(NN),ALPHAR(NN),BETA(NN),RWORK(8*NN),WORK(10*NN))
    
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
    IF (Evecs) THEN
      CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ELSE
      CALL DGGEV('N','N',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs(nphi,grid,VR)
    DEALLOCATE(grid,phi,psi,chi,A,B,C,D,VR,VL,ALPHAI,ALPHAR,BETA,RWORK,WORK) 
  END SUBROUTINE linear_const_sa
  
  SUBROUTINE bspline_deriv_sa() !sa stands for self-adjoint
    IMPLICIT NONE
    TYPE(bspline), DIMENSION(:), ALLOCATABLE :: phi, psi, chi 

    vcyl = .false.
    ALLOCATE(grid(N))
    IF (abs(mt).eq.1.and.(.not.Lend0)) THEN
      ALLOCATE(phi(0:N),psi(0:N),chi(0:N))
    ELSE
      ALLOCATE(phi(1:N),psi(1:N),chi(1:N))
    ENDIF
    NN = size(phi)+size(psi)+size(chi)
    ALLOCATE(A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),VR(NN,NN),VL(1,NN),ALPHAI(NN),ALPHAR(NN),BETA(NN),RWORK(8*NN),WORK(10*NN))
    
  !Initialize the finite elements
    lower(1:3) = (/lbound(phi,1), lbound(psi,1), lbound(chi,1)/)
    upper(1:3) = (/ubound(phi,1), ubound(psi,1), ubound(chi,1)/) 

    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    grid = new_grid(grid)
    nphi =  size(phi); npsi = size(psi); nchi = size(chi)
    IF (lower(1).eq.0) THEN
      CALL init(phi(0),grid(1),p3=grid(2),LendZero=.true.)
    ENDIF
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
    IF (Evecs) THEN
      CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ELSE
      CALL DGGEV('N','N',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs(nphi,grid,VR)   
    DEALLOCATE(grid,phi,psi,chi,A,B,C,D,VR,VL,ALPHAI,ALPHAR,BETA,RWORK,WORK) 
  END SUBROUTINE bspline_deriv_sa
  SUBROUTINE linear_const()
    IMPLICIT NONE
    TYPE(linear),   DIMENSION(:), ALLOCATABLE :: phi1, phi4
    TYPE(linear),   DIMENSION(:), ALLOCATABLE :: phi2, phi5
    TYPE(constant), DIMENSION(:), ALLOCATABLE :: phi3, phi6
    REAL(r8), DIMENSION(:,:), ALLOCATABLE :: A, B, C, D, VR, VL
    REAL(r8), DIMENSION(:), ALLOCATABLE :: ALPHAR, ALPHAI, BETA, RWORK, WORK
    REAL(r8) :: xi1pNa, xi2pNa, xi1Na, xi2Na, xi3Na,xi1pN1a,xi2pN1a
    INTEGER :: lower1, upper1
    WRITE (*,*) '************Using linear_const**************'
    vcyl = .true.
    lower1 = 1
    upper1 = N-1
  ! The wall is not at the plasma surface  
    IF (br.gt.ar) THEN
      CALL init_bc
      upper1 = N
    ENDIF
  ! Allocate the grid, finite elements, vectors, and matrices
    ALLOCATE(grid(N))    
    ALLOCATE(phi1(lower1:upper1),phi2(lower1:upper1),phi3(lower1:upper1),phi4(lower1:upper1),phi5(lower1:upper1),phi6(lower1:upper1))
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3); nphi4 = size(phi4); nphi5 = size(phi5); nphi6 = size(phi6)
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)
    LWORK=10*NN
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)
    ALLOCATE(A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),VR(NN,NN),VL(1,NN),ALPHAI(NN),ALPHAR(NN),BETA(NN),RWORK(8*NN),WORK(10*NN))
    LDVR = NN
  ! Initialize the grid
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    grid = new_grid(grid)
  ! Set the arrays containing the lower and upper bounds
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
  ! Initialize the finite elements  
    CALL init(phi1(1),grid(1),p3=grid(2))
    CALL init(phi1(2:N-1),grid(2:N-1),p2=grid(1:N-2),p3=grid(3:N))
    CALL init(phi3(2:N),grid(1:N-1),p3=grid(2:N))
    IF(upper(1).eq.N) THEN
      CALL init(phi1(N),grid(N),p2=grid(N-1))
      !CALL init(phi3(N),grid(N-1),p3=grid(N))
    ENDIF
    phi2 = phi1
    phi4 = phi1
  ! These 4 lines modify the linear elements to properly
  ! allow for div(xi)=0
    phi1%mod_lin = .true.
    phi4%mod_lin = .true.
    CALL linear_int_fac(phi1)
    CALL linear_int_fac(phi4)
    
    phi5 = phi2
    phi6 = phi3
    IF (verbose) THEN
      WRITE (*,'(a,g)') 'Ka=',Ka,'La=',La,'Kb=',Kb,'Lb=',Lb,'Kadot=',Kadot,'Ladot=',Ladot,'Kbdot=',Kbdot,'Lbdot=',Lbdot
      WRITE (*,'(a,g)') 'BCB1=',BCB1((/ar/)),'BCB1nw=',BCB1nw((/ar/)),'BCB1cw=',BCB1cw((/ar/))
      WRITE (*,'(a,g)') 'a=',ar,'val(phi1(N),a)=',val(phi1(N),(/ar/)), 'phi1(N)%xj-a=',phi1(N)%xj-ar
    ENDIF
  ! Set some format specifiers for outputting matrices
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi1,'G20.11)'
    
    CALL cpu_time(at1)
  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    A = 0.;    B = 0.;    C = 0.;    D = 0.
  !k is the row of the submatrix
    k = 1
    m = k+3
    ! Note that I've made the assumption that phi1=phi4, phi2=phi5, phi3=phi6
    DO i=lower(m),upper(m)
    ! l is the column of the submatrix
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),B11)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),A12)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),B12b,deriv1=.true.)+int_func(phi2(j),phi4(i),B12c)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi4(i),B13)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 2
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(m))+l,6*(i-lower(l))+k)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi5(i),phi2(j),A22a,deriv1=.true.,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),A22c)
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempC
        tempA = int_func(phi5(i),phi2(j),B22a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi5(i),phi2(j),B22b,deriv1=.true.)+int_func(phi5(i),phi2(j),B22b,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),B22c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi5(i),B23)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 3
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),A33)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),B33)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+3,6*(j-lower(l))+3)
      ENDDO
    ENDDO
    k = 4
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi1(i),phi1(j),B41a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi1(i),phi1(j),B41b,deriv1=.true.)+int_func(phi1(i),phi1(j),B41b,deriv2=.true.)
        tempC = int_func(phi1(i),phi1(j),B41c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(j),phi1(i),B42a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(j),phi1(i),B42bb,deriv2=.true.)+int_func(phi2(j),phi1(i),B42b,deriv1=.true.)
        tempC = int_func(phi2(j),phi1(i),B42c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi1(i),B43b,deriv2=.true.)
        tempC = int_func(phi3(j),phi1(i),B43c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC 
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+m,6*(j-lower(m))+m)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 5
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+2)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(i),phi2(j),B52a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(i),phi2(j),B52b,deriv1=.true.)+int_func(phi2(i),phi2(j),B52b,deriv2=.true.)
        tempC = int_func(phi2(i),phi2(j),B52c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi2(i),B53b,deriv2=.true.)
        tempC = int_func(phi3(j),phi2(i),B53c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 6
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+3)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+5,6*(i-lower(m))+3)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(i),phi3(j),B63)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+3,6*(j-lower(m))+3)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    ! This is for a wall not at the plasma surface
    IF (br.gt.ar) THEN
      xi1Na = 1.
      xi1pNa = minval(mod_lin_func((/ar/)))/phi1(N)%int_fac(1)*0.
      xi1pN1a = -minval(mod_lin_func((/ar/)))/phi1(N-1)%int_fac(2)*0.
      xi2Na = 1.
      xi2pNa = 1./phi2(N)%dx(1)*0.
      xi2pN1a = -1./phi2(N-1)%dx(2)*0.
      xi3Na = 1.*0.
      
      ! These are the extra terms after integrating by parts, some replaced by boundary conditions
      k = 2
      m = k+3
      i = N
        l = 1
        j = N
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB24a*xi1Na
        l = 2
        j = N-1
          C(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCA25aa*xi2pN1a
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB25aa*xi2pN1a
        j = N
          C(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCA25aa*xi2pNa
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB25aa*xi2pNa+BCB25ac*xi2Na
        l = 5
        j = N-1
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCA25aa*xi2pN1a
        j = N
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCA25aa*xi2pNa
      k = 4
      m = k-3
      !DO i=lower(m),upper(m)
      i = N
        l = 1
        !DO j=max(i-3,lower(l)),min(i+3,upper(l))
        j = N
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  (BCB11a+BCB11nw)*xi1Na
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  (BCB11a+BCB11cw)*xi1Na
          ENDIF
        !ENDDO
        l = 2
        !DO j=max(i-3,lower(l)),min(i+3,upper(l))
        j = N
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*(BCB11a+BCB11nw)*xi2Na
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*(BCB11a+BCB11cw)*xi2Na
          ENDIF
        !ENDDO
      !ENDDO
      k = 5
      m = k-3
      !DO i=lower(m),upper(m)
      i = N
        l = 1
        !DO j=max(i-3,lower(l)),min(i+3,upper(l))
        j = N-1
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB21aa*xi1pN1a
        j = N
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  (mt*(BCB11a+BCB11nw)+BCB21ac)*xi1Na+BCB21aa*xi1pNa
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  (mt*(BCB11a+BCB11cw)+BCB21ac)*xi1Na+BCB21aa*xi1pNa
          ENDIF
        !ENDDO
        l = 2
        !DO j=max(i-3,lower(l)),min(i+3,upper(l))
        j = N-1
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB22aa*xi2pN1a
        j = N
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  (mt**2*(BCB11a+BCB11nw)+BCB22ac)*xi2Na+BCB22aa*xi2pNa
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  (mt**2*(BCB11a+BCB11cw)+BCB22ac)*xi2Na+BCB22aa*xi2pNa
          ENDIF
        !ENDDO
        l = 3
        j = N
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB23a*xi3Na
        l = 4
        j = N
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB24a*xi1Na
        l = 5
        j = N-1
          C(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCA25aa*xi2pN1a
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB25aa*xi2pN1a
        j = N
          C(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCA25aa*xi2pNa
          D(6*(i-lower(m))+k,6*(j-lower(l))+l) = BCB25ac*xi2Na+BCB25aa*xi2pNa
      !ENDDO
    ENDIF
    B=B+D !Sign?
    A=A+C
  ! This is to account for regularity at the center
    !IF (mt.ne.1) THEN
      IF (lower(1).le.1) THEN
        DO i = 1,1
          m = 1
          DO k=1,4,3
            A(6*(i-lower(m))+k,:) = 0
            B(6*(i-lower(m))+k,:) = 0
            A(:,6*(i-lower(m))+k) = 0
            B(:,6*(i-lower(m))+k) = 0
          ENDDO
          k=1;  l=1;  j=i
          A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
          k=4;  l=4;  j=i
          A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        ENDDO
      ENDIF
    !ENDIF
    IF (lower(2).le.1) THEN
      DO i = 1,1
        j = i
        m = 2
        DO k=2,5,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=2;  l=2
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=5;  l=5
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
    IF (lower(3).le.1) THEN
      DO i = 1,1
        j = i
        m = 3
        DO k=3,6,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=3;  l=3
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=6;  l=6
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
      
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
      OPEN(1,status='replace',file='C.txt')
      WRITE (1,*) 'C='
      WRITE (1,FMT) transpose(C)
      CLOSE(1)
      OPEN(1,status='replace',file='D.txt')
      WRITE (1,*) 'D='
      WRITE (1,FMT) transpose(D)
      CLOSE(1)
    ENDIF
    C=A(:,:)
    D=B(:,:)
    CALL cpu_time(st1)
    IF (Evecs) THEN
      CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ELSE
      CALL DGGEV('N','N',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs(nphi1,grid,VR)
    IF (verbose) THEN
      WRITE (*,*) ALLOCATED(grid),ALLOCATED(phi1),ALLOCATED(phi2),ALLOCATED(phi3),&
      & ALLOCATED(phi4),ALLOCATED(phi5),ALLOCATED(phi6),&
      & ALLOCATED(A),ALLOCATED(B),ALLOCATED(C),ALLOCATED(D),&
      & ALLOCATED(VR),ALLOCATED(VL),ALLOCATED(ALPHAI),ALLOCATED(ALPHAR),ALLOCATED(BETA),&
      & ALLOCATED(RWORK),ALLOCATED(WORK)
    ENDIF
    DEALLOCATE(grid,phi1,phi2,phi3,phi4,phi5,phi6,A,B,C,D,VR,VL,ALPHAI,ALPHAR,BETA,RWORK,WORK,stat=k) 
    IF (VERBOSE) THEN
      WRITE(*,*) 'DEALLOCATION STATUS IS', k
    ENDIF
  END SUBROUTINE linear_const
  SUBROUTINE linear_const_rw()
    IMPLICIT NONE
    TYPE(linear),   DIMENSION(:), ALLOCATABLE :: phi1, phi4
    TYPE(linear),   DIMENSION(:), ALLOCATABLE :: phi2, phi5
    TYPE(constant), DIMENSION(:), ALLOCATABLE :: phi3, phi6
    COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: A, B, C, D, VR, VL
    COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: ALPHA, BETA,  WORK
    INTEGER :: lower1, upper1
    WRITE (*,*) '************Using linear_const_rw**************'
    vcyl = .true.
    lower1 = 1
    upper1 = N-1
  ! The wall is not at the plasma surface  
    IF (br.gt.ar) THEN
      CALL init_bc
      upper1 = N
    ENDIF
  ! Allocate the grid, finite elements, vectors, and matrices
    ALLOCATE(grid(N))    
    ALLOCATE(phi1(lower1:upper1),phi2(lower1:upper1),phi3(lower1:upper1),phi4(lower1:upper1),phi5(lower1:upper1),phi6(lower1:upper1))
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3); nphi4 = size(phi4); nphi5 = size(phi5); nphi6 = size(phi6)
  ! We're going to add 1 to NN for the coefficient C2 of the vacuum solution
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)+1
    LWORK=10*NN
    ALLOCATE(A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),VR(NN,NN),VL(1,NN),ALPHA(NN),BETA(NN),RWORK(8*NN),WORK(LWORK))
    LDVR = NN
  ! Initialize the grid
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    grid = new_grid(grid)
  ! Set the arrays containing the lower and upper bounds
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
  ! Initialize the finite elements  
    CALL init(phi1(1),grid(1),p3=grid(2))
    CALL init(phi1(2:N-1),grid(2:N-1),p2=grid(1:N-2),p3=grid(3:N))
    CALL init(phi3(2:N),grid(1:N-1),p3=grid(2:N))
    IF(upper(1).eq.N) THEN
      CALL init(phi1(N),grid(N),p2=grid(N-1))
      !CALL init(phi3(N),grid(N-1),p3=grid(N))
    ENDIF
    phi2 = phi1
    phi4 = phi1
  ! These 4 lines modify the linear elements to properly
  ! allow for div(xi)=0
    phi1%mod_lin = .true.
    phi4%mod_lin = .true.
    CALL linear_int_fac(phi1)
    CALL linear_int_fac(phi4)
    
    phi5 = phi2
    phi6 = phi3
    IF (verbose) THEN
      WRITE (*,'(a,g)') 'Ka=',Ka,'La=',La,'Kb=',Kb,'Lb=',Lb,'Kadot=',Kadot,'Ladot=',Ladot,'Kbdot=',Kbdot,'Lbdot=',Lbdot
      WRITE (*,'(a,g)') 'BCB1=',BCB1((/ar/)),'BCB1nw=',BCB1nw((/ar/)),'BCB1cw=',BCB1cw((/ar/))
      WRITE (*,'(a,g)') 'a=',ar,'val(phi1(N),a)=',val(phi1(N),(/ar/)), 'phi1(N)%xj-a=',phi1(N)%xj-ar
      write (*,'(a,g)') 'NN=',NN
    ENDIF
  ! Set some format specifiers for outputting matrices
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi1,'G20.11)'
    
    CALL cpu_time(at1)
  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    A = 0.;    B = 0.;    C = 0.;    D = 0.
  !k is the row of the submatrix
    k = 1
    m = k+3
    ! Note that I've made the assumption that phi1=phi4, phi2=phi5, phi3=phi6
    DO i=lower(m),upper(m)
    ! l is the column of the submatrix
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),B11)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),A12)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),B12b,deriv1=.true.)+int_func(phi2(j),phi4(i),B12c)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi4(i),B13)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 2
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(m))+l,6*(i-lower(l))+k)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi5(i),phi2(j),A22a,deriv1=.true.,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),A22c)
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempC
        tempA = int_func(phi5(i),phi2(j),B22a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi5(i),phi2(j),B22b,deriv1=.true.)+int_func(phi5(i),phi2(j),B22b,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),B22c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi5(i),B23)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 3
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),A33)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),B33)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+3,6*(j-lower(l))+3)
      ENDDO
    ENDDO
    k = 4
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi1(i),phi1(j),B41a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi1(i),phi1(j),B41b,deriv1=.true.)+int_func(phi1(i),phi1(j),B41b,deriv2=.true.)
        tempC = int_func(phi1(i),phi1(j),B41c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(j),phi1(i),B42a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(j),phi1(i),B42bb,deriv2=.true.)+int_func(phi2(j),phi1(i),B42b,deriv1=.true.)
        tempC = int_func(phi2(j),phi1(i),B42c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi1(i),B43b,deriv2=.true.)
        tempC = int_func(phi3(j),phi1(i),B43c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC 
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+m,6*(j-lower(m))+m)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 5
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+2)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(i),phi2(j),B52a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(i),phi2(j),B52b,deriv1=.true.)+int_func(phi2(i),phi2(j),B52b,deriv2=.true.)
        tempC = int_func(phi2(i),phi2(j),B52c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi2(i),B53b,deriv2=.true.)
        tempC = int_func(phi3(j),phi2(i),B53c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 6
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+3)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+5,6*(i-lower(m))+3)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(i),phi3(j),B63)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+3,6*(j-lower(m))+3)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    ! This is for a wall not at the plasma surface
    IF (br.gt.ar) THEN
      ! These are the extra terms after integrating by parts, replaced by boundary conditions
      k = 4
      m = k-3
      DO i=upper(m)-2,upper(m)
        l = 1
        DO j=max(i-3,lower(l)),min(i+3,upper(l))
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  minval(val(phi1(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi1(j),(/ar/)))
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  minval(val(phi1(i),(/ar/))*(BCB1((/ar/))+BCB1cw((/ar/)))*val(phi1(j),(/ar/)))
          ELSE! This is for a resistive wall
            ! This is for unknown C2
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  minval(val(phi1(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi1(j),(/ar/)))
          ENDIF
        ENDDO
        l = 2
        DO j=max(i-3,lower(l)),min(i+3,upper(l))
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*minval(val(phi1(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi2(j),(/ar/)))
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*minval(val(phi1(i),(/ar/))*(BCB1((/ar/))+BCB1cw((/ar/)))*val(phi2(j),(/ar/)))
          ELSE! This is for a resistive wall
            ! This is for unknown C2
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*minval(val(phi1(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi2(j),(/ar/)))
          ENDIF
        ENDDO
      ! This is for unknown C2
        D(6*(i-lower(m))+k,NN) = cmplx(0,1,r8)*minval(val(phi1(i),(/ar/))*(mt*Bt((/ar/))+kz*ar*Bz((/ar/))))*(La*Kadot-Ladot*Ka)/Kadot
      ENDDO
      k = 5
      m = k-3
      DO i=upper(m)-2,upper(m)
        l = 1
        DO j=max(i-3,lower(l)),min(i+3,upper(l))
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*minval(val(phi2(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi1(j),(/ar/)))
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*minval(val(phi2(i),(/ar/))*(BCB1((/ar/))+BCB1cw((/ar/)))*val(phi1(j),(/ar/)))
          ELSE! This is for a resistive wall
            ! This is for unknown C2
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt*minval(val(phi2(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi1(j),(/ar/)))
          ENDIF
        ENDDO
        l = 2
        DO j=max(i-3,lower(l)),min(i+3,upper(l))
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt**2*minval(val(phi2(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi2(j),(/ar/)))
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt**2*minval(val(phi2(i),(/ar/))*(BCB1((/ar/))+BCB1cw((/ar/)))*val(phi2(j),(/ar/)))
          ELSE! This is for a resistive wall
            ! This is for unknown C2
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  mt**2*minval(val(phi2(i),(/ar/))*(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi2(j),(/ar/)))
          ENDIF
        ENDDO
      ! This is for unknown C2
        D(6*(i-lower(m))+k,NN) = mt*cmplx(0,1,r8)*minval(val(phi2(i),(/ar/))*(mt*Bt((/ar/))+kz*ar*Bz((/ar/))))*(La*Kadot-Ladot*Ka)/Kadot/ar
      ENDDO
    ENDIF
    l=1
    ! This is for unknown C2
    DO j=lower(l),upper(l)
      C(NN,6*(j-lower(l))+l) = br*Kbdot*minval((mt*Bt((/ar/))+kz*ar*Bz((/ar/)))*val(phi1(j),(/ar/)))&
      & /ar**2/(mt**2+kz**2*br**2)/Kadot
      D(NN,6*(j-lower(l))+l) = 0
    ENDDO
    l=2
    ! This is for unknown C2
    DO j=lower(l),upper(l)
      C(NN,6*(j-lower(l))+l) = mt*br*Kbdot*minval((mt*Bt((/ar/))+kz*ar*Bz((/ar/)))*val(phi2(j),(/ar/)))&
      & /ar**2/(mt**2+kz**2*br**2)/Kadot
      D(NN,6*(j-lower(l))+l) = 0
    ENDDO
    C(NN,NN) = cmplx(0,1,r8)*kz*br*(Kbdot*Ladot-Kadot*Lbdot)/(mt**2+kz**2*br**2)/Kadot
    D(NN,NN) = -(Kbdot*Lb-Kb*Lbdot)/Kbdot/tw
    B=B+D!Sign?
    A=A+C
  ! This is to account for regularity at the center
    !IF (mt.ne.1) THEN
      IF (lower(1).le.1) THEN
        DO i = 1,1
          m = 1
          DO k=1,4,3
            A(6*(i-lower(m))+k,:) = 0
            B(6*(i-lower(m))+k,:) = 0
            A(:,6*(i-lower(m))+k) = 0
            B(:,6*(i-lower(m))+k) = 0
          ENDDO
          k=1;  l=1;  j=i
          A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
          k=4;  l=4;  j=i
          A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        ENDDO
      ENDIF
    !ENDIF
    IF (lower(2).le.1) THEN
      DO i = 1,1
        j = i
        m = 2
        DO k=2,5,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=2;  l=2
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=5;  l=5
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
    IF (lower(3).le.1) THEN
      DO i = 1,1
        j = i
        m = 3
        DO k=3,6,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=3;  l=3
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=6;  l=6
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
      
    CALL cpu_time(at2)
    at = at+at2-at1
    IF(verbose.and.num.eq.fin) THEN
      OPEN(1,status='replace',file='A.txt')
      WRITE (1,*) 'A='
      WRITE (1,FMT) REAL(transpose(A))
      CLOSE(1)
      OPEN(1,status='replace',file='AI.txt')
      WRITE (1,*) 'AI='
      WRITE (1,FMT) AIMAG(transpose(A))
      CLOSE(1)
      OPEN(1,status='replace',file='B.txt')
      WRITE (1,*) 'B='
      WRITE (1,FMT) REAL(transpose(B))
      CLOSE(1)
      OPEN(1,status='replace',file='BI.txt')
      WRITE (1,*) 'BI='
      WRITE (1,FMT) AIMAG(transpose(B))
      CLOSE(1)
      OPEN(1,status='replace',file='D.txt')
      WRITE (1,*) 'D='
      WRITE (1,FMT) REAL(transpose(D))
      CLOSE(1)
      WRITE(0,*) 'Alfven range=',alfven_range(grid)
    ENDIF
    C=A
    D=B
    CALL cpu_time(st1)
    IF (Evecs) THEN
      CALL ZGGEV('N','V',NN,B,NN,A,NN,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)
    ELSE
      CALL ZGGEV('N','N',NN,B,NN,A,NN,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    IF (maxval(abs(AIMAG(BETA))).EQ.0) THEN
      CALL output_evals(REAL(ALPHA)/REAL(BETA),AIMAG(ALPHA)/REAL(BETA))
    ELSE
      CALL output_evals(REAL(ALPHA/BETA),AIMAG(ALPHA/BETA))
    ENDIF
    !CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    IF (Evecs) THEN
      CALL output_evecs(nphi1,grid,REAL(VR))
      CALL output_evecs(nphi1,grid,AIMAG(VR),imag=.true.)
    ENDIF
    !CALL output_evecs(nphi1,grid,VR)
    !IF (verbose) THEN
      !WRITE (*,*) ALLOCATED(grid),ALLOCATED(phi1),ALLOCATED(phi2),ALLOCATED(phi3),&
      !& ALLOCATED(phi4),ALLOCATED(phi5),ALLOCATED(phi6),&
      !& ALLOCATED(A),ALLOCATED(B),ALLOCATED(C),ALLOCATED(D),&
      !& ALLOCATED(VR),ALLOCATED(VL),ALLOCATED(ALPHAI),ALLOCATED(ALPHAR),ALLOCATED(BETA),&
      !& ALLOCATED(RWORK),ALLOCATED(WORK)
    !ENDIF
    DEALLOCATE(grid,phi1,phi2,phi3,phi4,phi5,phi6,A,B,C,D,VR,VL,ALPHA,BETA,RWORK,WORK,stat=k) 
    !IF (VERBOSE) THEN
      !WRITE(*,*) 'DEALLOCATION STATUS IS', k
    !ENDIF
  END SUBROUTINE linear_const_rw
  
  SUBROUTINE bspline_deriv()
    IMPLICIT NONE
    TYPE(bspline), DIMENSION(:), ALLOCATABLE :: phi1, phi2, phi3, phi4, phi5, phi6
    INTEGER :: lower1, upper1
    
    vcyl = .true.
    WRITE (*,*) '************Using bspline_deriv**************'
    upper1 = N
  ! The wall is not at the plasma surface  
    IF (br.gt.ar) THEN
      CALL init_bc
      upper1 = N+1
    ENDIF
    IF (verbose) THEN
      WRITE (*,'(a,g)') 'Ka=',Ka,'La=',La,'Kb=',Kb,'Lb=',Lb,'Kadot=',Kadot,'Ladot=',Ladot,'Kbdot=',Kbdot,'Lbdot=',Lbdot
      WRITE (*,'(a,g)') 'BCB1nw=',BCB1nw((/ar/)),'BCB1cw=',BCB1cw((/ar/))
    ENDIF

    lower1=0
    
    IF (verbose) THEN
      WRITE (*,*) 'lower1=',lower1,'upper1=',upper1
    ENDIF
    
  ! Allocate the grid, finite elements, vectors, and matrices
    ALLOCATE(grid(N))    
    ALLOCATE(phi1(lower1:upper1),phi2(lower1:upper1),phi3(lower1:upper1),phi4(lower1:upper1),phi5(lower1:upper1),phi6(lower1:upper1))
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3); nphi4 = size(phi4); nphi5 = size(phi5); nphi6 = size(phi6)
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)
    LWORK=10*NN
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)
    ALLOCATE(A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),VR(NN,NN),VL(1,NN),ALPHAI(NN),ALPHAR(NN),BETA(NN),RWORK(8*NN),WORK(10*NN))
    LDVR = NN
    
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    grid = new_grid(grid)
    grid(1) = grid(2)/100.0
  !Initialize the finite elements
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
    
    IF(lower(1).eq.0)   CALL init(phi1(0),grid(1),p3=grid(2),LendZero=.true.)
    IF(lower(1).le.1)   CALL init(phi1(1),grid(1),p3=grid(2),p4=grid(3),LendZero=.true.)
    IF(lower(1).le.2)   CALL init(phi1(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=.true.)
    CALL init(phi1(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    CALL init(phi1(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi1(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
    IF(upper(1).eq.N+1) CALL init(phi1(N+1),grid(N),p2=grid(N-1),RendZero=.true.)
    phi2 = phi1
    phi3 = phi1
    phi4 = phi1
    phi5 = phi1
    phi6 = phi1
    phi2%deriv = phi2_deriv
    phi3%deriv = phi3_deriv
    phi5%deriv = phi2_deriv
    phi6%deriv = phi3_deriv
  ! These are unique elements that enforce behavior at r=0
    !phi2(1)=phi2(1)*phi2(0)%C(1)/phi2(1)%C(1)-phi2(0)
    phi3(1) = phi3(1)*phi3(0)%C(2)/phi3(1)%C(2)-phi3(0)
    phi3(0) = phi3(3)+shift_spline(phi3(2))
    !phi3(0)%C=(/0.,0.,-3*phi3(0)%dx(3),2./)
    !CALL shift_spline(phi3(0))
    !WRITE(0,*) phi3(0)
    !phi3(2)=phi3(2)*(2*phi3(0)%B(2)-6*phi3(0)%B(3)*phi3(0)%xj)/&
    !& (2*phi3(2)%B(2)-6*phi3(2)%B(3)*phi3(2)%xj)-phi3(0)
    !WRITE(0,*) phi2(1)
    !WRITE(0,*) phi3(1)
    !WRITE(0,*) phi3(2)
    
  ! Set some format specifiers for outputting matrices
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi1,'G20.11)'

    CALL cpu_time(at1)
  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    A = 0.;    B = 0.;    C = 0.;    D = 0.
  !k is the row of the submatrix
    k = 1
    m = k+3
    ! Note that I've made the assumption that phi1=phi4, phi2=phi5, phi3=phi6
    DO i=lower(m),upper(m)
    ! l is the column of the submatrix
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),B11)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),A12)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),B12b,deriv1=.true.)+int_func(phi2(j),phi4(i),B12c)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi4(i),B13)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 2
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(m))+k,6*(i-lower(l))+l)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi5(i),phi2(j),A22a,deriv1=.true.,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),A22c)
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempC
        tempA = int_func(phi5(i),phi2(j),B22a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi5(i),phi2(j),B22b,deriv1=.true.)+int_func(phi5(i),phi2(j),B22b,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),B22c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi5(i),B23)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 3
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),A33)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),B33)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+3,6*(j-lower(l))+3)
      ENDDO
    ENDDO
    k = 4
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi1(i),phi1(j),B41a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi1(i),phi1(j),B41b,deriv1=.true.)+int_func(phi1(i),phi1(j),B41b,deriv2=.true.)
        tempC = int_func(phi1(i),phi1(j),B41c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(j),phi1(i),B42a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(j),phi1(i),B42bb,deriv2=.true.)+int_func(phi2(j),phi1(i),B42b,deriv1=.true.)
        tempC = int_func(phi2(j),phi1(i),B42c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi1(i),B43b,deriv2=.true.)
        tempC = int_func(phi3(j),phi1(i),B43c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC 
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+m,6*(j-lower(m))+m)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 5
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+2)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(i),phi2(j),B52a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(i),phi2(j),B52b,deriv1=.true.)+int_func(phi2(i),phi2(j),B52b,deriv2=.true.)
        tempC = int_func(phi2(i),phi2(j),B52c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi2(i),B53b,deriv2=.true.)
        tempC = int_func(phi3(j),phi2(i),B53c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 6
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+3)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+5,6*(i-lower(m))+3)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(i),phi3(j),B63)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+3,6*(j-lower(m))+3)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    ! This is for a wall not at the plasma surface
    IF (br.gt.ar) THEN
      WRITE (0,*) 'The section for a free boundary has not been updated.'
      ! These are the extra terms after integrating by parts, replaced by boundary conditions
      k = 4
      m = k-3
      DO i=lower(m),upper(m)
        l = 1
        DO j=max(i-3,lower(l)),min(i+3,upper(l))
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  val(phi1(i),ar)*minval((BCB1((/ar/))+BCB1nw((/ar/))))*val(phi1(j),ar)
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  val(phi1(i),ar)*minval((BCB1((/ar/))+BCB1cw((/ar/))))*val(phi1(j),ar)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    B=B+D !Sign?
  ! This is to account for regularity at the center
    !IF (mt.ne.1) THEN
    IF (lower(1).lt.1) THEN
      DO i = 1,0
        m = 1
        DO k=1,4,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=1;  l=1;  j=i
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=4;  l=4;  j=i
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
    !ENDIF
    IF (lower(2).le.0) THEN
      DO  i = 1,0
        j = i
        m = 2
        DO k=2,5,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=2;  l=2
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=5;  l=5
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
    IF (lower(3).lt.1) THEN
      DO i = 0,1
        j = i
        m = 3
        DO k=3,6,3
          A(6*(i-lower(m))+k,:) = 0
          B(6*(i-lower(m))+k,:) = 0
          A(:,6*(i-lower(m))+k) = 0
          B(:,6*(i-lower(m))+k) = 0
        ENDDO
        k=3;  l=3
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
        k=6;  l=6
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
      ENDDO
    ENDIF
      
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
      OPEN(1,status='replace',file='D.txt')
      WRITE (1,*) 'D='
      WRITE (1,FMT) transpose(D)
      CLOSE(1)
      WRITE(0,*) 'Alvfen range=',alfven_range(grid)
    ENDIF
    C=A(:,:)
    D=B(:,:)
    CALL cpu_time(st1)
    IF (Evecs) THEN
      CALL DGGEV('N','V',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ELSE
      CALL DGGEV('N','N',NN,B,NN,A,NN,ALPHAR,ALPHAI,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    CALL output_evals(ALPHAR/BETA,ALPHAI/BETA)
    CALL output_evecs(nphi1,grid,VR) 
    CALL output_error_real(ALPHAR,ALPHAI,BETA,VR,phi1,phi2,phi3)
    IF (verbose) THEN
      WRITE (*,*) ALLOCATED(grid),ALLOCATED(phi1),ALLOCATED(phi2),ALLOCATED(phi3),&
      & ALLOCATED(phi4),ALLOCATED(phi5),ALLOCATED(phi6),&
      & ALLOCATED(A),ALLOCATED(B),ALLOCATED(C),ALLOCATED(D),&
      & ALLOCATED(VR),ALLOCATED(VL),ALLOCATED(ALPHAI),ALLOCATED(ALPHAR),ALLOCATED(BETA),&
      & ALLOCATED(RWORK),ALLOCATED(WORK)
    ENDIF
    DEALLOCATE(grid,phi1,phi2,phi3,phi4,phi5,phi6,A,B,C,D,VR,VL,ALPHAI,ALPHAR,BETA,RWORK,WORK,stat=k) 
    IF (VERBOSE) THEN
      WRITE(*,*) 'DEALLOCATION STATUS IS', k
    ENDIF
  END SUBROUTINE bspline_deriv
  SUBROUTINE bspline_deriv_rw()
    IMPLICIT NONE
    COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: ALPHA, BETA, WORK
    COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: A, B, C, D, VR, VL
    TYPE(bspline), DIMENSION(:), ALLOCATABLE :: phi1, phi2, phi3, phi4, phi5, phi6
    INTEGER :: lower1, upper1
    WRITE (*,*) '************Using bspline_deriv_rw**************'
    vcyl = .true.
  ! Decide if there should be an element that allows phi1(r=0).ne.0
    IF (abs(mt).eq.1.and.(.not.Lend0)) THEN
      lower1=0
    ELSE
      lower1=1
    ENDIF
    lower1=0
  ! The wall is not at the plasma surface
    IF (br.gt.ar) THEN
      CALL init_bc
      upper1 = N+1
    ENDIF
    IF (verbose) THEN
      WRITE (*,*) 'lower1=',lower1,'upper1=',upper1
    ENDIF
  ! Allocate the grid, finite elements, vectors, and matrices
    ALLOCATE(grid(N))    
    ALLOCATE(phi1(lower1:upper1),phi2(lower1:upper1),phi3(lower1:upper1),phi4(lower1:upper1),phi5(lower1:upper1),phi6(lower1:upper1))
    nphi1 =  size(phi1); nphi2 = size(phi2); nphi3 = size(phi3); nphi4 = size(phi4); nphi5 = size(phi5); nphi6 = size(phi6)
  ! We're going to add 1 to NN for the coefficient C2 of the vacuum solution
    NN = size(phi1)+size(phi2)+size(phi3)+size(phi4)+size(phi5)+size(phi6)+1
    LWORK=10*NN
    ALLOCATE(A(NN,NN),B(NN,NN),C(NN,NN),D(NN,NN),VR(NN,NN),VL(1,NN),ALPHA(NN),BETA(NN),RWORK(8*NN),WORK(LWORK))
    LDVR = NN

  ! Initialize the grid
    grid = (/ (i*ar/(N-1), i=0,N-1) /)
    grid = new_grid(grid)
    grid(1) = grid(2)/1000.0
  
  ! Set the upper and lower bounds
    lower = (/lbound(phi1,1), lbound(phi2,1), lbound(phi3,1), lbound(phi4,1), lbound(phi5,1), lbound(phi6,1), lbound(phi1,1), lbound(phi2,1), lbound(phi3,1)/)
    upper = (/ubound(phi1,1), ubound(phi2,1), ubound(phi3,1), ubound(phi4,1), ubound(phi5,1), ubound(phi6,1), ubound(phi1,1), ubound(phi2,1), ubound(phi3,1)/) 
  ! Initialize the finite elements  
    IF(lower(1).eq.0)   CALL init(phi1(0),grid(1),p3=grid(2),LendZero=.true.)
    IF(lower(1).le.1)   CALL init(phi1(1),grid(1),p3=grid(2),p4=grid(3),LendZero=.true.)
    IF(lower(1).le.2)   CALL init(phi1(2),grid(2),p2=grid(1),p3=grid(3),p4=grid(4),LendZero=.true.)
    CALL init(phi1(3:N-2),grid(3:N-2),p1=grid(1:N-4),p2=grid(2:N-3),p3=grid(4:N-1),p4=grid(5:N))
    CALL init(phi1(N-1),grid(N-1),p1=grid(N-3),p2=grid(N-2),p3=grid(N),RendZero=.true.)
    CALL init(phi1(N),grid(N),p1=grid(N-2),p2=grid(N-1),RendZero=.true.)
    IF(upper(1).eq.N+1) CALL init(phi1(N+1),grid(N),p2=grid(N-1),RendZero=.true.)
    phi2 = phi1
    phi3 = phi1
    phi4 = phi1
    phi5 = phi1
    phi6 = phi1
    phi2%deriv = phi2_deriv
    phi3%deriv = phi3_deriv
    phi5%deriv = phi2_deriv
    phi6%deriv = phi3_deriv
    
  ! Set some format specifiers for outputting matrices
    WRITE(FMT,'(a,I,a)') '(',NN,'G13.5)'
    WRITE(FMTR,'(a,I,a)') '(',nphi1,'G20.11)'
    CALL cpu_time(at1)
  !Initialize the matrices and vectors needed for the eigenvalue/eigenvector solver
    A = 0.;    B = 0.;    C = 0.;    D = 0.
  !k is the row of the submatrix
    k = 1
    m = k+3
    ! Note that I've made the assumption that phi1=phi4, phi2=phi5, phi3=phi6
    DO i=lower(m),upper(m)
    ! l is the column of the submatrix
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),A11)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi4(i),phi1(j),B11)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),A12)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi2(j),phi4(i),B12b,deriv1=.true.)+int_func(phi2(j),phi4(i),B12c)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi4(i),B13)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 2
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(j-lower(m))+k,6*(i-lower(l))+l)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi5(i),phi2(j),A22a,deriv1=.true.,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),A22c)
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempC
        tempA = int_func(phi5(i),phi2(j),B22a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi5(i),phi2(j),B22b,deriv1=.true.)+int_func(phi5(i),phi2(j),B22b,deriv2=.true.)
        tempC = int_func(phi5(i),phi2(j),B22c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(j),phi5(i),B23)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+2,6*(j-lower(l))+2)
      ENDDO
    ENDDO
    k = 3
    m = k+3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+l,6*(i-lower(l))+k)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),A33)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi6(i),phi3(j),B33)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+3,6*(j-lower(l))+3)
      ENDDO
    ENDDO
    k = 4
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi1(i),phi1(j),B41a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi1(i),phi1(j),B41b,deriv1=.true.)+int_func(phi1(i),phi1(j),B41b,deriv2=.true.)
        tempC = int_func(phi1(i),phi1(j),B41c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(j),phi1(i),B42a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(j),phi1(i),B42bb,deriv2=.true.)+int_func(phi2(j),phi1(i),B42b,deriv1=.true.)
        tempC = int_func(phi2(j),phi1(i),B42c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi1(i),B43b,deriv2=.true.)
        tempC = int_func(phi3(j),phi1(i),B43c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC 
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+m,6*(j-lower(m))+m)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(m))+1,6*(j-lower(l))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+1,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 5
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+2)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempA = int_func(phi2(i),phi2(j),B52a,deriv1=.true.,deriv2=.true.)
        tempB = int_func(phi2(i),phi2(j),B52b,deriv1=.true.)+int_func(phi2(i),phi2(j),B52b,deriv2=.true.)
        tempC = int_func(phi2(i),phi2(j),B52c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempA + tempB + tempC
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        tempB = int_func(phi3(j),phi2(i),B53b,deriv2=.true.)
        tempC = int_func(phi3(j),phi2(i),B53c)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = tempB + tempC
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+1)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+2,6*(j-lower(m))+2)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+2,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    k = 6
    m = k-3
    DO i=lower(m),upper(m)
      l = 1
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(m))+4,6*(i-lower(l))+3)
      ENDDO
      l = 2
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(j-lower(l))+5,6*(i-lower(m))+3)
      ENDDO
      l = 3
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = int_func(phi3(i),phi3(j),B63)
      ENDDO
      l = 4
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+1)
      ENDDO
      l = 5
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+2)
      ENDDO
      l = 6
      DO j=max(i-3,lower(l)),min(i+3,upper(l))
        A(6*(i-lower(m))+k,6*(j-lower(l))+l) = A(6*(i-lower(l))+3,6*(j-lower(m))+3)
        B(6*(i-lower(m))+k,6*(j-lower(l))+l) = B(6*(i-lower(l))+3,6*(j-lower(m))+3)
      ENDDO
    ENDDO
    ! This is for a wall not at the plasma surface
    IF (br.gt.ar) THEN
      WRITE (0,*) 'The section for a free boundary has not been updated.'
    ! These are the extra terms from integrating by parts, replaced by perturbed pressure balance
      k = 4
      m = k-3
      DO i=lower(m),upper(m)
        l = 1
        DO j=max(i-3,lower(l)),min(i+3,upper(l))
          IF (br.gt.1000.0.or.tw.eq.0) THEN ! This for no wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  val(phi1(i),ar)*minval(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi1(j),ar)
          ELSE IF (tw.lt.0) THEN ! This is for a conducting wall
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  val(phi1(i),ar)*minval(BCB1((/ar/))+BCB1cw((/ar/)))*val(phi1(j),ar)
          ELSE ! This is for a resistive wall
            ! This is for unknown C2
            D(6*(i-lower(m))+k,6*(j-lower(l))+l) =  val(phi1(i),ar)*minval(BCB1((/ar/))+BCB1nw((/ar/)))*val(phi1(j),ar)
          ENDIF
        ENDDO
        ! This is for unknown C2
        D(6*(i-lower(m))+k,NN) = val(phi1(i),ar)*cmplx(0,1,r8)*minval(mt*Bt((/ar/))+kz*ar*Bz((/ar/)))*(La*Kadot-Ladot*Ka)/Kadot/ar
      ENDDO
    ENDIF
    l=1
    ! This is for unknown C2
    DO j=lower(l),upper(l)
 ! There was a k,kz, m,mt mistake in the following line
      C(NN,6*(j-lower(l))+l) = br*Kbdot*minval(mt*Bt((/ar/))+kz*ar*Bz((/ar/)))&
      & /ar**2/(mt**2+kz**2*br**2)/Kadot*val(phi1(j),ar)
      D(NN,6*(j-lower(l))+l) = 0
    ENDDO
    C(NN,NN) = cmplx(0,1,r8)*kz*br*(Kbdot*Ladot-Kadot*Lbdot)/(mt**2+kz**2*br**2)/Kadot
    D(NN,NN) = -(Kbdot*Lb-Kb*Lbdot)/Kbdot/tw
    B=B+D
    A=A+C
    ! This is to account for regularity at the center
    !IF (mt.ne.1) THEN
      IF (lower(1).lt.1) THEN
        DO i = 0,0
          m = 1
          DO k=1,4,3
            A(6*(i-lower(m))+k,:) = 0
            B(6*(i-lower(m))+k,:) = 0
            A(:,6*(i-lower(m))+k) = 0
            B(:,6*(i-lower(m))+k) = 0
          ENDDO
          k=1;  l=1;  j=i
          A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
          k=4;  l=4;  j=i
          A(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1.0
          !k=4;  l=1;  j=1
          !B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
          !k=1;  l=4;  j=1
          !B(6*(i-lower(m))+k,6*(j-lower(l))+l) = 1e-20
        ENDDO
      ENDIF
    CALL cpu_time(at2)
    at = at+at2-at1
    IF(verbose.and.num.eq.fin) THEN
      OPEN(1,status='replace',file='A.txt')
      WRITE (1,*) 'A='
      WRITE (1,FMT) REAL(transpose(A))
      CLOSE(1)
      OPEN(1,status='replace',file='AI.txt')
      WRITE (1,*) 'AI='
      WRITE (1,FMT) AIMAG(transpose(A))
      CLOSE(1)
      OPEN(1,status='replace',file='B.txt')
      WRITE (1,*) 'B='
      WRITE (1,FMT) REAL(transpose(B))
      CLOSE(1)
      OPEN(1,status='replace',file='BI.txt')
      WRITE (1,*) 'BI='
      WRITE (1,FMT) AIMAG(transpose(B))
      CLOSE(1)
      OPEN(1,status='replace',file='D.txt')
      WRITE (1,*) 'D='
      WRITE (1,FMT) REAL(transpose(D))
      CLOSE(1)
      WRITE(0,*) 'Alfven range=',alfven_range(grid)
    ENDIF
    WRITE (*,'(a,g)') 'Ka=',Ka,'La=',La,'Kb=',Kb,'Lb=',Lb,'Kadot=',Kadot,'Ladot=',Ladot,'Kbdot=',Kbdot,'Lbdot=',Lbdot
    C = A
    D = B
    CALL cpu_time(st1)
    IF (Evecs) THEN
      CALL ZGGEV('N','V',NN,B,NN,A,NN,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)
    ELSE
      CALL ZGGEV('N','N',NN,B,NN,A,NN,ALPHA,BETA,VL,LDVL,VR,LDVR,WORK,LWORK,RWORK,INFO)
    ENDIF
    CALL cpu_time(st2)
    st = st+st2-st1
    CALL output_params(at2-at1,st2-st1)
    IF (VERBOSE) THEN
      WRITE (*,*) 'ALPHA = '
      WRITE (*,'(2G)') ALPHA
      WRITE (*,*) 'BETA = '
      WRITE (*,'(2G)') BETA
    ENDIF
    IF (maxval(abs(AIMAG(BETA))).EQ.0) THEN
      CALL output_evals(REAL(ALPHA)/REAL(BETA),AIMAG(ALPHA)/REAL(BETA))
    ELSE
      CALL output_evals(REAL(ALPHA/BETA),AIMAG(ALPHA/BETA))
    ENDIF
    IF (Evecs) THEN
      CALL output_evecs(nphi1,grid,REAL(VR))
      CALL output_evecs(nphi1,grid,AIMAG(VR),imag=.true.)
      CALL output_error(ALPHA,BETA,VR,phi1,phi2,phi3)
    ENDIF
    DEALLOCATE(grid,phi1,phi2,phi3,phi4,phi5,phi6,A,B,C,D,VR,VL,ALPHA,BETA,RWORK,WORK) 
  END SUBROUTINE bspline_deriv_rw 
  
  SUBROUTINE plot_equilibrium
    INTEGER num_pts
    INTEGER :: i
    REAL(r8), DIMENSION(N) :: grid, pressure, zb, tb, density, pv, zv, safety, pp, bzp, btp
    num_pts = N
    grid = (/ (i*ar/real(num_pts-1), i=0,num_pts-1) /)
    pressure = P(grid)
    tb = Bt(grid)
    zb = Bz(grid)
    density = rho(grid)
    pv = Vp(grid)*grid
    zv = Vz(grid)
    pp = P_prime(grid)
    bzp = Bz_prime(grid)
    btp = Bt_prime(grid)
    safety(1) = 0.
    safety(2:N) = q(grid(2:N))
    WRITE(filename,'(a,i0,a)') 'equilibria_vcyl/',num,'.txt'
    IF(.not.vcyl) THEN
      pv = 0.
      zv = 0.
    ENDIF
    WRITE (*,'(a,i0,a,a)') 'Writing equilibrium data for run number ',num, ' to ',filename
    OPEN(1,file=filename,status='replace')
    DO i=1,num_pts
      WRITE (1,'(8g25.16)') grid(i), pressure(i),zb(i),tb(i),density(i),pv(i),zv(i), safety(i)
    ENDDO
    CLOSE(1)
    IF (verbose) THEN
      WRITE (*,*) 'Eqilibrium condition at grid points='
      WRITE (*,'(g)') equilibrium_v(grid(2:))
      WRITE (*,'(5(a20))') 'pressure','p_prime', 'Bz*Bz_prime', 'Bt*Bt_prime', 'Bt**2/r'
      DO i=2,size(grid)
        WRITE (*,'(5(g20.10))') pressure(i), pp(i), zb(i)*bzp(i), tb(i)*btp(i), tb(i)**2/grid(i)
      ENDDO
    ENDIF
  END SUBROUTINE plot_equilibrium

  SUBROUTINE output_params(assemble_t,solve_t)
    REAL(r8), INTENT(IN) :: assemble_t,solve_t
    WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.dat'
    OPEN(1,file=filename,status='replace')
    IF (spline) THEN
      WRITE(1,'(a)') 'spline'
    ELSE 
      WRITE(1,'(a)') 'linconst'
    ENDIF
    WRITE(1,'(i)') N, NN, mt, equilib, num, BCrow
    WRITE(1,'(g)') epsilo, assemble_t, solve_t, kz, ar, br, tw, rho0, Bz0, Bt0, s2, eps, P0, P1, lambd, Vz0, epsVz, Vp0, epsVp, rs, alpha
    WRITE(1,'(l)') Lend0, vcyl
    CLOSE(1)
  END SUBROUTINE output_params
  SUBROUTINE output_evals(evalsr,evalsi)
    REAL(r8), DIMENSION(:), INTENT(IN) :: evalsr, evalsi
    REAL(r8), DIMENSION(size(evalsr)) :: evalsr1, evalsi1
    evalsr1 = evalsr
    evalsi1 = evalsi
    IF (.not.vcyl) THEN
      evalsr1 = real((evalsr+(0,1.0)*evalsi)**0.5)
      evalsi1 = aimag((evalsr+(0,1.0)*evalsi)**0.5)
    ENDIF
    WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.evalsr'
    OPEN(1,file=filename,status='replace')
    WRITE(1,'(g)') evalsr1
    CLOSE(1)
    WRITE(filename,'((a,i0),a)') 'output_vcyl/', num,'.evalsi'
    OPEN(1,file=filename,status='replace')
    WRITE(1,'(g)') evalsi1
    CLOSE(1)    
  END SUBROUTINE output_evals
  SUBROUTINE output_evecs(nphi,grid,VR,imag)
    INTEGER, INTENT(IN) :: nphi
    REAL(r8), DIMENSION(:), INTENT(IN) :: grid
    REAL(r8), DIMENSION(:,:), INTENT(IN) :: VR
    LOGICAL, INTENT(IN), OPTIONAL :: imag
    LOGICAL :: imag1
    INTEGER :: N , i, j, k, i_skip
    CHARACTER(LEN=30) FMT
    imag1 = .false.
    IF(present(imag)) imag1 = imag
    IF (Evecs) THEN
      i_skip = 3
      IF (vcyl) THEN
        i_skip=6
      ENDIF
      N = size(grid)
      WRITE(FMT,'(a,i0,a)') '(',N,'(g,a))'
      WRITE(filename,'(a,i0,a)') 'output_vcyl/',num,'.grid'
      OPEN(1,file=filename,status='replace')
      WRITE(1,FMT) (grid(i),',',i=1,size(grid))
      CLOSE(1)
      WRITE(FMT,'(a,i0,a)') '(',nphi,'(g,a))'
      DO j=1,3
        WRITE(filename,'(2(a,i0))') 'output_vcyl/',num,'.evecs',j
        IF(imag1) WRITE(filename,'(2(a,i0))') 'output_vcyl/',num,'.evecs_imag',j
        OPEN(1,file=filename,status='replace')
        DO k=1,size(VR(1,:))
          WRITE(1,FMT) ( VR(i,k),',' , i=j,nphi*i_skip,i_skip )
        ENDDO
        CLOSE(1)
      ENDDO 
    ENDIF
  END SUBROUTINE output_evecs
  SUBROUTINE output_error_real(ALPHAR,ALPHAI,BETA,VR,phi1,phi2,phi3)
    REAL(r8), INTENT(IN), DIMENSION(:,:) :: VR
    REAL(r8), INTENT(IN), DIMENSION(:) :: ALPHAR, ALPHAI, BETA
    TYPE(bspline), DIMENSION(:), INTENT(IN) :: phi1, phi2, phi3
    COMPLEX(r8), DIMENSION(size(ALPHAR)) :: tempR, tempL, temp1, temp2, temp3
    WRITE(filename,'((a,i0,a))') 'output_vcyl/',num,'.BCerr'
    OPEN(1,file=filename,status='replace') 
    DO i=1,size(ALPHAR)
      temp1(i) = minval(BJCP1((/ar/)))*sum(VR(1::6,i)*val(phi1,ar)) + minval(BJCP1a((/ar/)))*sum(VR(1::6,i)*val_prime(phi1,ar)) 
      temp2(i) = minval(BJCP2((/ar/)))*sum(VR(2::6,i)*val(phi2,ar))
      temp3(i) = minval(BJCP3((/ar/)))*sum(VR(3::6,i)*val(phi3,ar))
      tempL(i) = temp1(i)+temp2(i)+temp3(i)
      IF (br.gt.1000.0.or.tw.eq.0) THEN !No Wall
        WRITE (*,*) 'There is no wall'
        tempR(i) = sum(val(phi1,ar)*VR(1::6,i))*minval(BCB1nw((/ar/)))
      ELSEIF (tw.lt.0) THEN !Conducting Wall
        WRITE (*,*) 'There is a conducting wall'
        tempR(i) = sum(val(phi1,ar)*VR(1::6,i))*minval(BCB1cw((/ar/)))
      ENDIF
      IF (ALPHAI(i).gt.0) THEN
        temp1(i) = temp1(i)+cmplx(0,1,r8)*(minval(BJCP1((/ar/)))*sum(VR(1::6,i+1)*val(phi1,ar)) + minval(BJCP1a((/ar/)))*sum(VR(1::6,i+1)*val_prime(phi1,ar)))
        temp2(i) = temp2(i)+cmplx(0,1,r8)*minval(BJCP2((/ar/)))*sum(VR(2::6,i+1)*val(phi2,ar))
        temp3(i) = temp3(i)+cmplx(0,1,r8)*minval(BJCP3((/ar/)))*sum(VR(3::6,i+1)*val(phi3,ar))
        tempL(i) = temp1(i)+temp2(i)+temp3(i)
        IF (br.gt.1000.0.or.tw.eq.0) THEN !No Wall
          tempR(i) = tempR(i)+cmplx(0,1,r8)*minval(BCB1nw((/ar/)))*sum(VR(1::6,i+1)*val(phi1,ar))
        ELSEIF (tw.lt.0) THEN !Conducting Wall
          tempR(i) = tempR(i)+cmplx(0,1,r8)*minval(BCB1cw((/ar/)))*sum(VR(1::6,i+1)*val(phi1,ar))
        ENDIF
      ELSEIF (ALPHAI(i).lt.0) THEN
        temp1(i) = conjg(temp1(i-1))
        temp2(i) = conjg(temp2(i-1))
        temp3(i) = conjg(temp3(i-1))
        tempL(i) = temp1(i)+temp2(i)+temp3(i)
        tempR(i) = conjg(tempR(i-1))
       ENDIF
    ENDDO
    WRITE(1,'(f)') abs(tempL-tempR)
    CLOSE(1)
    IF(VERBOSE) THEN
      WRITE(*,'(4A40)') 'temp1', 'temp2', 'temp3', 'tempR'
      WRITE(*,'(8G20.10)') ((/temp1(i),temp2(i),temp3(i),tempR(i)/),i=1,size(ALPHAI))
    ENDIF
  END SUBROUTINE output_error_real
  SUBROUTINE output_error(ALPHA,BETA,VR,phi1,phi2,phi3)
    COMPLEX(r8), INTENT(IN), DIMENSION(:,:) :: VR
    COMPLEX(r8), INTENT(IN), DIMENSION(:) :: ALPHA, BETA
    TYPE(bspline), DIMENSION(:), INTENT(IN) :: phi1, phi2, phi3
    COMPLEX(r8), DIMENSION(size(ALPHA)) :: tempR, tempL, temp2, temp3, temp1
    WRITE(filename,'((a,i0,a))') 'output_vcyl/',num,'.BCerr'
    OPEN(1,file=filename,status='replace') 
    temp1 = (/ (minval(BJCP1((/ar/)))*sum(VR(1::6,i)*val(phi1,ar)) + minval(BJCP1a((/ar/)))*sum(VR(1::6,i)*val_prime(phi1,ar)) , i=1,size(ALPHA) ) /)
    temp2 = (/ (minval(BJCP2((/ar/)))*sum(VR(2::6,i)*val(phi2,ar)), i=1,size(ALPHA))/)
    temp3 = (/ (minval(BJCP3((/ar/)))*sum(VR(3::6,i)*val(phi3,ar)), i=1,size(ALPHA))/)
    tempL = (/ ((minval(BJCP1((/ar/)))*sum(VR(1::6,i)*val(phi1,ar)) + minval(BJCP1a((/ar/)))*sum(VR(1::6,i)*val_prime(phi1,ar)) &
    & +minval(BJCP2((/ar/)))*sum(VR(2::6,i)*val(phi2,ar)) + minval(BJCP3((/ar/)))*sum(VR(3::6,i)*val(phi3,ar))), i=1,size(ALPHA) ) /)
    IF (br.gt.1000.0.or.tw.eq.0) THEN !No Wall
      tempR = (/ (sum(val(phi1,ar)*VR(1::6,i))*minval(BCB1nw((/ar/))),i=1,size(ALPHA) ) /)
    ELSEIF (tw.lt.0) THEN !Conducting Wall
      tempR = (/ (sum(val(phi1,ar)*VR(1::6,i))*minval(BCB1cw((/ar/))),i=1,size(ALPHA) ) /)
    ELSE
      tempR = (/ (sum(val(phi1,ar)*VR(1::6,i))/2.*&
      & (minval(n0((/ar/)))*CMPLX(0,1,r8) + minval(n1((/ar/)))*ALPHA(i)/BETA(i))/ &
      & (minval(d0((/ar/)))*CMPLX(0,1,r8) + minval(d1((/ar/)))*ALPHA(i)/BETA(i)),i=1,size(ALPHA) ) /)
    ENDIF
    WRITE(1,'(f)') abs(tempL-tempR)/abs(temp1)
    CLOSE(1)
    IF(VERBOSE) THEN
      WRITE(*,'(4A40)') 'temp1', 'temp2', 'temp3', 'tempR'
      WRITE(*,'(8G20.10)') ((/temp1(i),temp2(i),temp3(i),tempR(i)/),i=1,size(ALPHA))
    ENDIF
  END SUBROUTINE output_error
END PROGRAM cyl
