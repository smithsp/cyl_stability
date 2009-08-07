PROGRAM test
  USE local
  !USE nag_lib_support, ONLY: nag_lib_ident
  !USE nag_inv_hyp_fun, ONLY: nag_arccosh
  USE nag_bessel_fun, ONLY: nag_bessel_j0, nag_bessel_k
  INTEGER ifail, i, j, N
  REAL(r8) :: x(2),y(2)
  !CALL nag_lib_ident
  x= 1.
  write(*,*) nag_bessel_j0(x(1))
  N=10
  DO j=1,10
    DO i=1,N
      WRITE(*,*) i/real(N),nag_bessel_k(cmplx(i/real(N),0,r8),real(j,r8))
    ENDDO
  ENDDO
  STOP
END PROGRAM test
