PROGRAM tise_1D_numerov
  !
  ! Main program to solve the time independent Schroedinger equation in 
  ! one dimension using a Numerov approach
  !
  ! Bound states calculation
  !
  ! by Currix TM
  !
  USE nrtype
  USE phys_constants
  USE global_vars
  USE pot_param
  USE numerov_alg_box
  !
  IMPLICIT NONE
  !
  !
  REAL(KIND = DP) :: energy, eigenvalue, e_min, x_e_min, delta_engy, engy_threshold = 0.0_DP 
  REAL(KIND = DP) :: tol_minpack, tol_wf
  REAL(KIND = DP) :: wdiff0, wdiff1 
  REAL(KIND = DP), DIMENSION(4) :: eig_result
  REAL(KIND = DP), DIMENSION(1) :: eigenval, diff_wf_val
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: en_matrix, wf_vector
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: wf_matrix
  REAL(KIND = DP) :: norm_val, right_val
  INTEGER(KIND = I4B) :: quanta = 0, n_states = 0, index, ierr, i_sum_rules
  !
  ! Pointer for eigenvalues
  TYPE element
     REAL(KIND = DP) :: energy_val
     TYPE(element), POINTER :: next
  END TYPE element
  !
  TYPE(element), TARGET :: head_eigenval
  TYPE(element), POINTER :: current_eigenval, temp_eigenval
  !
  !
  INTERFACE X_SUM_RULE
     SUBROUTINE X_SUM_RULE(nstates, wf_matrix)
       !
       !
       USE nrtype
       USE phys_constants
       USE global_vars
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: nstates
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: wf_matrix
       !
     END SUBROUTINE X_SUM_RULE
  END INTERFACE X_SUM_RULE
  !
  INTERFACE X2_SUM_RULE
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
       INTEGER(KIND = I4B), INTENT(IN) :: nstates
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: wf_matrix
       !
     END SUBROUTINE X2_SUM_RULE
  END INTERFACE X2_SUM_RULE
  !
  INTERFACE EW_SUM_RULE
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
       INTEGER(KIND = I4B), INTENT(IN) :: nstates
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: wf_matrix
       !
     END SUBROUTINE EW_SUM_RULE
  END INTERFACE EW_SUM_RULE
  !
  ! READING INPUT
  !
  NAMELIST/INP_DEBUG/ iprint
  NAMELIST/INP_X/ x_min, x_max, npoints ! Box walls at xmin and xmax
  NAMELIST/INP_POTENTIAL/ param_pot
  NAMELIST/INP_SYSTEM/ red_mass
  NAMELIST/INP_NUMEROV/ match_p, tol_minpack, tol_wf, delta_engy
  NAMELIST/INP_BOX/ engy_threshold
  NAMELIST/INP_OUTPUT/  i_wf_save, wf_filename, i_e_save, e_filename, i_sum_rules
  !
  READ(UNIT=*,NML=INP_DEBUG)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "iprint = ", iprint
     WRITE(*,*)
  ENDIF
  !
  READ(UNIT=*,NML=INP_X)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "x_min = ", x_min, ", x_max = ", x_max, ", npoints = ", npoints
     WRITE(*,*)
  ENDIF
  !
  READ(UNIT=*,NML=INP_POTENTIAL)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "potential params = ", param_pot
     WRITE(*,*)
  ENDIF
  !
  READ(UNIT=*,NML=INP_SYSTEM)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "reduced mass (amu) = ", red_mass
     WRITE(*,*)
  ENDIF
  !
  READ(UNIT=*,NML=INP_NUMEROV)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "match point index = ", match_p
     WRITE(*,*)
     WRITE(*,*) "tol minpack = ", tol_minpack, ", tol wf = ", tol_wf
     WRITE(*,*)
     WRITE(*,*) "delta energy = ", delta_engy
     WRITE(*,*)
  ENDIF
  !
  READ(UNIT=*,NML=INP_BOX)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "Energy threshold = ", engy_threshold
     WRITE(*,*)
  ENDIF
  !
  READ(UNIT=*,NML=INP_OUTPUT)
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "iprint = ", iprint
     WRITE(*,*)
     WRITE(*,*) "wf save = ", i_wf_save, ", wf filename = ", wf_filename, "sum rules = ",  i_sum_rules
     WRITE(*,*)
     WRITE(*,*) "e save = ", i_e_save, ", e filename = ", e_filename
  ENDIF
  ! 
  !     DEFINE X GRID
  !
  ALLOCATE(x_grid(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "x_grid allocation request denied."
     STOP
  ENDIF
  !
  x_step = (x_max - x_min)/REAL(npoints - 1, DP)  
  x_grid = x_min + x_step*REAL((/ (index, index = 0, npoints - 1) /),DP)
  !
  IF (iprint > 1) WRITE(*,*) ' X grid step = ', x_step, 'fm'
  !
  !    EVALUATE POTENTIAL
  !
  ALLOCATE(vpot(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*), "vpot allocation request denied."
     STOP
  ENDIF
  !
  vpot = Potf(x_grid)
  !
  !    WAVE FUNCTION VECTOR
  !
  ALLOCATE(wf_vector(1:npoints), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*), "wf_vector allocation request denied."
     STOP
  ENDIF
  !
  !
  ! Main Loop
  !
  e_min = MINVAL(vpot)
  x_e_min = x_grid(MINLOC(ARRAY = vpot, DIM = 1))
  !
  IF (iprint > 1) THEN 
     WRITE(*,*) "E min = ", e_min,", x min = ", x_e_min, MINLOC(ARRAY = vpot, DIM = 1)
  ENDIF
  !
  !
  energy = e_min
  !
  CALL wf_numerov_diff(energy, wdiff0)
  !
  DO WHILE (energy < engy_threshold)
     !
     CALL wf_numerov_diff(energy + delta_engy, wdiff1) 
     ! 
     IF (iprint > 1) WRITE(*,*) "Energy", energy, "wdiff1", wdiff1
     !
     IF (wdiff0*wdiff1 < 0) THEN
        !
        IF (iprint > 1) THEN 
           WRITE(*,*) "Possible eigenvalue in [", energy, "; ", energy + delta_engy, "] -> [ ", wdiff0, "; ", wdiff1, "]"
        ENDIF
        !
        !
        ! initial x value
        eigenval(1) = energy + delta_engy/2.0_DP
        !
        ! initial diff value
        CALL wf_numerov_diff_min(1, eigenval, diff_wf_val, ierr)
        !
        IF (iprint > 1) THEN 
           WRITE(*,*) "Initial value for minimization", eigenval, diff_wf_val
        ENDIF
        !
        CALL hybrd1(wf_numerov_diff_min, 1, eigenval, diff_wf_val, tol_minpack, ierr)
        !
        IF (iprint > 1) THEN
           WRITE ( *, * ) ' '
           WRITE ( *, '(a,i6)' ) '  Returned value of INFO = ', ierr
           WRITE ( *, *) eigenval, diff_wf_val
        ENDIF
        !
        eigenvalue = eigenval(1)
        ! 
        ! Check if the eigenvalue found is in the interval under study
        !
        IF (eigenvalue < energy .OR. eigenvalue > energy + delta_engy ) THEN ! Check whether the eigenvalue found is out of the search interval
           !
           IF ( iprint > 1 ) WRITE(*,*) " Wrong eigenvalue buddy... Keep on searching!"
           !
        ELSE
           !
           ! Compute the relative difference of the right and left wave functions
           !
           eig_result = wf_Numerov_lr(eigenvalue)
           !
           ! Renormalize output
           norm_val = 1.0_DP/eig_result(3)
           right_val = 1.0_DP/eig_result(4)
           eig_result(1) = eig_result(1)*norm_val  ! wfl
           eig_result(2) = eig_result(2)*right_val ! wfr
           eig_result(3) = 1.0_DP ! wfpl
           eig_result(4) = 1.0_DP ! wfpr
           !
           IF ( iprint > 1 ) THEN
              !
              WRITE(*,*) "Possible eigenvalue = ", eigenvalue
              !
              WRITE(*,*) "Matching data  Ia ", eig_result
              WRITE(*,*) "Matching data IIa ",eig_result(1)-eig_result(2), &
                   2*ABS((eig_result(1)-eig_result(2))/(eig_result(1)+eig_result(2)))
              !
           ENDIF
           !
           ! Build normalized wave function
           !
           CALL wf_numerov(eigenvalue, wf_vector)
           !
           norm_val = wf_vector(match_p)/eig_result(1)
           eig_result(1) = eig_result(1)*norm_val
           eig_result(2) = eig_result(2)*norm_val
           !
           WRITE(*,*) "Possible eigenvalue : ", eigenvalue, norm_val, wf_vector(match_p), eig_result(1)
           !
           IF ( iprint >= 1 ) THEN
              !
              WRITE(*,*) "Possible eigenvalue = ", eigenvalue
              !
              WRITE(*,*) "Matching data  Ib ", eig_result
              WRITE(*,*) "Matching data IIb ",eig_result(1)-eig_result(2), &
                   2*ABS((eig_result(1)-eig_result(2))/(eig_result(1)+eig_result(2)))
              !
           ENDIF
           !
           IF ( ABS(eig_result(1)-eig_result(2)) < tol_wf ) THEN
              !
              IF ( iprint >= 1 ) THEN
                 !
                 WRITE(*,*) n_states,"-th eigenvalue = ", eigenval
                 !
              ENDIF
              ! DEFINE HEAD OF EIGENVALUES POINTER
              IF (n_states == 0) THEN 
                 !
                 head_eigenval%energy_val = eigenvalue
                 NULLIFY(head_eigenval%next)
                 current_eigenval => head_eigenval
                 !
              ELSE
                 !
                 ! FILL LIST EIGENVALUE LIST
                 !
                 ALLOCATE(temp_eigenval, STAT = ierr)
                 IF (ierr /= 0) STOP 'STOP :: POINTER ALLOCATION FAILED.' 
                 temp_eigenval%energy_val = eigenvalue
                 NULLIFY(temp_eigenval%next)
                 current_eigenval%next => temp_eigenval 
                 current_eigenval => temp_eigenval 
                 !
              ENDIF
              !
              !
              n_states = n_states + 1
              nodes_num = nodes_num + 1
              CALL wf_numerov_diff(energy + delta_engy, wdiff1)
              !
           ENDIF
           !
        ENDIF
        !
     ENDIF
     !
     energy = energy + delta_engy
     !
     wdiff0 = wdiff1
     ! 
  END DO
  !
  !
  ! PRINT THE EIGENVALUE LIST
  !
  ALLOCATE(en_matrix(1:n_states), STAT = ierr)    
  IF (ierr /= 0) THEN
     WRITE(*,*) "en_matrix allocation request denied."
     STOP
  ENDIF
  !
  current_eigenval => head_eigenval
  quanta = 0
  DO
     IF (.NOT.ASSOCIATED(current_eigenval)) EXIT
     WRITE(*,*) quanta, current_eigenval%energy_val
     en_matrix(quanta + 1) =  current_eigenval%energy_val
     quanta = quanta + 1 
     current_eigenval => current_eigenval%next
  ENDDO
  !
  ! SAVE THE EIGENVALUES
  IF (i_e_save /= 0) THEN
     !
     OPEN(UNIT = 20, FILE = e_filename, STATUS = 'UNKNOWN')
     !
     DO quanta = 0, n_states - 1
        WRITE(20,*) quanta, en_matrix(quanta + 1)
     ENDDO
     !
  ENDIF
  ! 
  ! SAVE THE WAVE FUNCTION
  IF (i_wf_save /= 0) THEN
     !
     OPEN(UNIT = 30, FILE = wf_filename, STATUS = 'UNKNOWN')
     !
     ALLOCATE(wf_matrix(1:npoints,1:n_states + 1), STAT = ierr) ! x_grid and wavefunctions    
     IF (ierr /= 0) THEN
        WRITE(*,*) "wf_matrix allocation request denied."
        STOP
     ENDIF
     !
     ! First column x_grid
     wf_matrix(1:npoints,1) = x_grid(1:npoints)
     !
     DO quanta = 0, n_states - 1
        CALL wf_numerov(en_matrix(quanta + 1), wf_vector)
        wf_matrix(1:npoints,quanta + 2) = wf_vector(1:npoints)
     ENDDO
     !
     DO index = 1, npoints
        WRITE(30,*) wf_matrix(index, 1:n_states + 1) ! First column is xgrid
     ENDDO
     !
  ENDIF
  !
  ! SUM RULES CALCULATION 
  !
  IF(i_sum_rules /= 0) THEN
     !
     ! SUM RULES --- X operator
     CALL X_SUM_RULE(n_states+1, wf_matrix)
     !
     ! SUM RULES  --- X2 operator
     CALL X2_SUM_RULE(n_states+1, wf_matrix)
     !
     ! ENERGY  WEIGHTED SUM RULES --- X2 operator
     CALL EW_SUM_RULE(n_states+1, wf_matrix)
     !
  ENDIF
  !
  IF (iprint >= 1) WRITE(*,*) 'Sayonara baby...'
  !
END PROGRAM tise_1D_numerov
