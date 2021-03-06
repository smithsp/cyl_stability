USE cyl_funcs_module
USE vcyl_funcs_module
real :: a,b,c
real :: ddd(3,5) , d1(3), d2(3), d3(3), d4(3), d5(3)
integer, parameter :: N=41
real, dimension(N) :: grid
integer :: i
equilib = 4
kz=1.2
ar=1.0
mt=-1
P0=.0
P1=.1
lambd=2
Bz0=10.0
Bt0=.5
gamma=5./3.
eps=0.5
s2=1./12.
rho0=1.0
Vp0=1.0
epsVp=.5
Vz0=1.0
epsVz=.5

grid = (/ (i*ar/real(N-1), i=0,N-1) /)
grid(1) = grid(2)/10000.
write (*,*) '0=',2.D0/real(1.D16)
write (*,*) 'Inf=',Inf

DO equilib=1,13
  WRITE (*,*) 'equilib = ', equilib
  write (*,'(g)') equilibrium_v(grid)
ENDDO
equilib = 13
write (*,*) 'grid'
write (*,'(g)') grid
write (*,*) 'Bz'
write (*,'(g)') Bz(grid)
write (*,*) 'Bz_prime'
write (*,'(g)') Bz_prime(grid)
write (*,*) 'Bt'
write (*,'(g)') Bt(grid)
write (*,*) 'Bt_prime'
write (*,'(g)') Bt_prime(grid)
write (*,*) 'Pressure'
write (*,'(g)') p(grid)
write (*,*) 'P_prime'
write (*,'(g)') p_prime(grid)
d1 = (/1,2,3/)
d2 = (/4,5,6/)
d3 = (/7,8,9/)
d4 = (/10,11,12/)
d5 = (/13,14,15/)
ddd = reshape((/d1,d2,d3,d4,d5/),(/3,5/))
!write(*,'(3(g))') ddd
!write(*,*) ddd(:,1)
a = 1e-10
!READ(*,'(g,g,g)') a,b,c
!WRITE (*,'(g,g,g)') a,b,c

end
