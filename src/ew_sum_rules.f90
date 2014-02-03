SUBROUTINE EW_SUM_RULE(nstates, wf_matrix)
  !
  ! by LauPK
  !
  USE nrtype
  USE phys_constants
  USE global_vars
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) ::  nstates
  REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: wf_matrix
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i
  REAL(KIND = DP) :: error,  Ews_temp
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: ew_vec
  !
  ALLOCATE(ew_vec(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "ew_vec allocation request denied."
     STOP
  ENDIF
  !
  ! ENERGY  WEIGHTED SUM RULES --- X2 operator
  WRITE(*,*) "   "
  WRITE(*,*) "E-W SUM RULE test value (x^2 operator)"
  !
  ! 
  DO i = 2, nstates !column 1 is the x_grid
     !
     ew_vec = 0.0_DP
     ew_vec = wf_matrix(:, i)*X_Grid(:)*X_Grid(:)*wf_matrix(:, i)
     !
     Ifail = 0
     Ews_temp = 0.0_DP
     CALL D01GAF(X_Grid, ew_vec, npoints, Ews_temp, Error, Ifail)
     !
     PRINT*, "from ", i-2, "-th state: ", (Ews_temp*2.0_DP)
     !
  ENDDO
  !
  DEALLOCATE(ew_vec, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "ew_vec deallocation request denied."
     STOP
  ENDIF
  !
  RETURN
  !
END SUBROUTINE EW_SUM_RULE
