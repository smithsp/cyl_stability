MODULE sort_module
  USE local
CONTAINS
  FUNCTION sort(vec,rev,inds)
    REAL(r8), DIMENSION(:), INTENT(IN) :: vec
    REAL(r8), DIMENSION(size(vec)) :: sort, tempvec
    REAL(r8) :: temp
    INTEGER :: ind, i, j, indices(size(vec)), tempi
    LOGICAL, OPTIONAL :: rev
    LOGICAL :: reverse, swapped
    INTEGER, OPTIONAL, DIMENSION(size(vec)), INTENT(OUT) :: inds
    reverse = .false.
    indices = (/ (i,i=1,size(vec)) /)
    IF (present(rev)) reverse = rev
    tempvec(:) = vec(:)
    ind=1
    DO j=1,size(sort)
      swapped = .false.
      DO i=1,size(sort)-ind
        IF (tempvec(i)>tempvec(i+1)) THEN
          temp = tempvec(i)
          tempvec(i) = tempvec(i+1)
          tempvec(i+1) = temp
          IF(present(inds)) THEN
            tempi = indices(i)
            indices(i) = indices(i+1)
            indices(i+1) =tempi
          ENDIF
          swapped = .true.
        ENDIF
      ENDDO
      IF(.not.swapped) EXIT
      ind = ind+1
    ENDDO
    
    IF(reverse) THEN
      sort = tempvec(size(vec):1:-1)
      IF(present(inds)) inds = indices(size(vec):1:-1)
    ELSE
      sort = tempvec(:)
      IF(present(inds)) inds = indices(:)
    ENDIF
  END FUNCTION sort
END MODULE sort_module
