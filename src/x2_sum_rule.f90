SUBROUTINE X2_SUM_RULE(nstates, wf_matrix)
  !
  ! SUM_i/=2 |<AVEC(2)|X^2|AVEC(i)>|^2 = <AVEC(2)|Y^2|AVEC(2)> - <AVEC(2)|Y|AVEC(2)>^2 
  ! where Y = X^2 - <AVEC(2)|X^2|AVEC(2)>
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
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, n
  REAL(KIND = DP) :: error, Test_value, res_scale
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: integrand, Rescaled_X_Grid, scale_vector
  !
  ALLOCATE(integrand(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "integrand allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Rescaled_X_Grid(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "Rescaled_X_Grid allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(scale_vector(1:npoints), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "scale_vector allocation request denied."
     STOP
  ENDIF
  !
  ! SUM RULES  --- X2 operator
  WRITE(*,*) "   "
  WRITE(*,*) "SUM RULE test value (x^2 operator)"
  !
  DO i = 2, nstates !column 1 is the x_grid
     !
     integrand = 0.0_DP
     Rescaled_X_Grid = 0.0_DP
     scale_vector = 0.0_DP
     Test_value = 0.0_DP
     res_scale = 0.0_DP
     ! 
     !OPERATOR SCALING
     scale_vector = wf_matrix(:,i)*(X_grid(:)**2)*wf_matrix(:,i)
     Ifail = 0
     CALL D01GAF(X_Grid, scale_vector, npoints, res_scale, error, Ifail)
     !
     DO n = 1, npoints 
        Rescaled_X_Grid(n) = X_grid(n)**2 - res_scale
     ENDDO
     !
     ! m_0 calculation
     integrand = wf_matrix(:,i)*(Rescaled_X_Grid(:)**2)*wf_matrix(:,i)
     Ifail = 0
     CALL D01GAF(X_Grid, integrand, npoints, Test_value, error, Ifail)
     !
     PRINT*, "from ", i-2, "-th state: ", Test_value
     !
  ENDDO
  !
  DEALLOCATE(integrand, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "integrand deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Rescaled_X_Grid, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Rescaled_X_Grid deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(scale_vector, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "scale_vector deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE X2_SUM_RULE
