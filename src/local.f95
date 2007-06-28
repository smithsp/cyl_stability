
MODULE local
 
 IMPLICIT NONE

 INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND (4)
 INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND (8)
 INTEGER, PARAMETER :: r8 = KIND(0.d0) !SELECTED_REAL_KIND (15,50) !
 REAL(r8), PARAMETER  :: pi = 3.14159265358979323846_r8, two_pi = 2*pi
 LOGICAL :: verbose
END MODULE local         
