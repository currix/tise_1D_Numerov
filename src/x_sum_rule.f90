SUBROUTINE X_SUM_RULE(nstates, wf_matrix)
  !
  ! SUM_i/=2 |<AVEC(2)|X|AVEC(i)>|^2 = <AVEC(2)|X^2|AVEC(2)> - <AVEC(2)|X|AVEC(2)>^2 
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
  REAL(KIND = DP) :: error, Test_value, res_1, res_2
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: TEST_element_1, TEST_element_2
  !
  ALLOCATE(TEST_element_1(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "TEST_element_1 allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(TEST_element_2(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "TEST_element_2 allocation request denied."
     STOP
  ENDIF
  !
  ! SUM RULES --- X operator
  WRITE(*,*) "   "
  WRITE(*,*) "SUM RULE test value (x operator)"
  ! 
  DO i = 2, nstates !column 1 is the x_grid
     !
     TEST_element_1 = 0.0_DP
     TEST_element_2 = 0.0_DP
     TEST_element_1 = wf_matrix(:, i)*X_grid(:)*X_grid(:)*wf_matrix(:, i)
     TEST_element_2 = wf_matrix(:, i)*X_grid(:)*wf_matrix(:, i) 
     !
     Test_value = 0.0_DP
     res_1 = 0.0_DP
     res_2 = 0.0_DP
     !
     Ifail = 0
     CALL D01GAF(X_Grid, TEST_element_1, npoints, res_1, error, Ifail)
     Ifail = 0
     CALL D01GAF(X_Grid, TEST_element_2, npoints, res_2, error, Ifail)
     !
     Test_value = res_1 - res_2**2
     !
     PRINT*, "from ", i-2, "-th state: ", Test_value
     !
  ENDDO
  !
  DEALLOCATE(TEST_element_1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_1 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(TEST_element_2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_2 deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE X_SUM_RULE
