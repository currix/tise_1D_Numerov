MODULE numerov_alg
  !
  USE nrtype
  USE phys_constants
  USE global_vars
  USE pot_param
  !
  IMPLICIT NONE
  !
  ! DEFINITION OF GLOBAL VARIABLES
  !
  INTEGER(KIND = I4B) :: match_p ! index of match point
  !
  INTEGER(KIND = I4B) :: nodes_num ! number of nodes flag: even/odd wf
  !
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: vpot ! MeV potential values
  !
  ! CONTROLLING OUTPUT
  !
  INTEGER(KIND = I4B) :: i_wf_save ! save wave functions flag
  !
  CHARACTER(LEN = 64) :: wf_filename ! wave functions file name
  !
  INTEGER(KIND = I4B) :: i_e_save ! save energies flag
  !
  CHARACTER(LEN = 64) :: e_filename ! energies file name
  !
  ! Potential external function
  !
  INTERFACE POTF
     ELEMENTAL FUNCTION Potf(x)
       USE nrtype
       USE global_vars
       USE pot_param
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       REAL(KIND = DP), INTENT(IN) :: x
       !
       REAL(KIND = DP) ::  Potf
     END FUNCTION Potf
  END INTERFACE POTF  
  !
CONTAINS
  !
  FUNCTION wf_numerov_lr(energy)
    !
    ! energy input 
    REAL(KIND = DP), INTENT(IN) :: energy
    !
    ! wf_left wf_right wf_left' wf_right'
    REAL(KIND = DP), DIMENSION(4) :: wf_numerov_lr
    !
    ! Local variables
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: gn, aux_vec, y_left, y_right
    REAL(KIND = DP)  :: gn_min, gn_max, delta, nu_min, nu_min_plus, nu_min_minus
    INTEGER(KIND = I4B) :: index, ierr
    ! 
    ! Allocate vectors
    ALLOCATE(gn(1:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "gn allocation request denied."
       STOP
    ENDIF
    ALLOCATE(aux_vec(1:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "aux_vec allocation request denied."
       STOP
    ENDIF
    ALLOCATE(y_left(1:match_p + 2), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "y_left allocation request denied."
       STOP
    ENDIF
    ALLOCATE(y_right(match_p - 2:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "y_right allocation request denied."
       STOP
    ENDIF
    !
    ! Preliminary calculations
    !
    gn = (2.0_DP*amu*red_mass/(hbarc**2))*(energy - vpot)
    !
    aux_vec = 1.0_DP + ((x_step**2)/12.0_DP)*gn
    !
    gn = gn/aux_vec ! rescaled gn
    !
    gn_min = gn(1)
    gn_max = gn(npoints)
    !
    !delta = SQRT(x_step**2*gn_min*(x_step**2*gn_min-4.0_DP))
    delta = x_step*SQRT(gn_min*(x_step**2*gn_min-4.0_DP))
    !
    nu_min_plus = 0.5_DP*(x_step**2*gn_min + delta)
    nu_min_minus = 0.5_DP*(x_step**2*gn_min - delta)
    !
    ! Left branch  
    IF ( ABS(1.0_DP/(1.0_DP - nu_min_plus)) > 1 ) THEN
       nu_min = nu_min_plus
    ELSE
       nu_min = nu_min_minus
    ENDIF
    !
    y_left(1) = 1.0_DP
    y_left(2) = y_left(1)/(1.0_DP - nu_min)
    !
    y_left(1:2) = y_left(1:2)*aux_vec(1:2)  ! Renormalization
    !
    DO index = 3, match_p + 2
       !
       y_left(index) = 2.0_DP*y_left(index - 1) - y_left(index - 2) - y_left(index - 1)*gn(index-1)*x_step**2
       !
    ENDDO
    !
    y_left = y_left/aux_vec(1:match_p + 2)
    !
    !
    ! Right branch 
    IF ( ABS(1.0_DP - nu_min_plus) > 1 ) THEN
       nu_min = nu_min_plus
    ELSE
       nu_min = nu_min_minus
    ENDIF
    !
    IF (MOD(nodes_num,2) == 0) THEN
       y_right(npoints) = 1.0_DP ! even case
    ELSE
       y_right(npoints) = -ABS(1.0_DP) ! odd case
    ENDIF
    !
    y_right(npoints - 1) = y_right(npoints)*(1 - nu_min);
    !
    ! Rescaled wave function
    y_right(npoints - 1:npoints) = y_right(npoints - 1:npoints)*aux_vec(npoints - 1:npoints)  ! Renormalization
    !
    DO index = npoints - 2, match_p - 2, -1 ! Minus 2 for slope calculation
       !
       y_right(index) = 2*y_right(index + 1) - y_right(index + 2) - y_right(index + 1)*gn(index + 1)*x_step**2;
       !
    ENDDO
    !
    y_right = y_right / aux_vec(match_p - 2: npoints);
    !
    ! Prepare output
    !
    wf_numerov_lr(1) = y_left(match_p)   ! Left wf value
    wf_numerov_lr(2) = y_right(match_p)  ! Right wf value
    !
    wf_numerov_lr(3) = ( y_left(match_p + 1) - y_left(match_p - 1)  -&
         (1.0_DP/8.0_DP)*y_left(match_p + 2) + (1.0_DP/8.0_DP)*y_left(match_p - 2) ) /&
         (3.0_DP*x_step/2.0_DP) ! Left wf derivative
    !
    wf_numerov_lr(4) = ( y_right(match_p + 1) - y_right(match_p - 1)  -&
         (1.0_DP/8.0_DP)*y_right(match_p + 2) + (1.0_DP/8.0_DP)*y_right(match_p - 2) ) /&
         (3.0_DP*x_step/2.0_DP) ! Right wf derivative
    !
    !
    !    
  END FUNCTION wf_numerov_lr
  !
  SUBROUTINE wf_numerov_diff(energy, diff_val)
    !
    ! energy input 
    REAL(KIND = DP), INTENT(IN) :: energy
    !
    ! 
    REAL(KIND = DP), INTENT(OUT) :: diff_val
    !
    ! local variables
    REAL(KIND = DP), DIMENSION(4) :: temp_result
    REAL(KIND = DP) ::  norm_val, right_val
    !
    !
    temp_result = wf_numerov_lr(energy)
    !
    !
    ! Renormalize output
    norm_val = 1.0_DP/temp_result(1)
    right_val = temp_result(3)/temp_result(4)
    temp_result(1) = 1.0_DP ! wfl
    temp_result(2) = temp_result(2)*right_val*norm_val ! wfr
    temp_result(3) = temp_result(3)*norm_val ! wfpl
    temp_result(4) = temp_result(3) ! wfpr
    !
    diff_val = SIGN(1.0_DP, right_val)*(temp_result(1) - temp_result(2)) !  To have a definite criteria with function signs
    !
    !
  END SUBROUTINE wf_numerov_diff
  !
  SUBROUTINE wf_numerov_diff_min(n, energy, diff_val, iflag)
    !
    ! Number of equations
    INTEGER(KIND = I4B), INTENT(IN) :: n
    !
    ! energy input 
    REAL(KIND = DP), DIMENSION(n), INTENT(IN) :: energy
    !
    ! difference left/right output 
    REAL(KIND = DP), DIMENSION(n), INTENT(OUT) :: diff_val
    !
    ! flag
    INTEGER(KIND = I4B), INTENT(OUT) :: iflag
    !
    ! local variables
    REAL(KIND = DP), DIMENSION(4) :: temp_result
    REAL(KIND = DP) ::  norm_val, right_val 
    !
    !
    temp_result = wf_numerov_lr(energy(1))
    !
    !
    ! Renormalize output
    norm_val = 1.0_DP/temp_result(1)
    right_val = temp_result(3)/temp_result(4)
    temp_result(1) = 1.0_DP ! wfl
    temp_result(2) = temp_result(2)*right_val*norm_val ! wfr
    temp_result(3) = temp_result(3)*norm_val ! wfpl
    temp_result(4) = temp_result(3) ! wfpr
    !
    diff_val = SIGN(1.0_DP, right_val)*(temp_result(1) - temp_result(2)) !  To have a definite criteria with function signs
    !
    iflag = 0 !  Added to remove compilation warning. No use apart of this.
    !
  END SUBROUTINE wf_numerov_diff_min
  !
  SUBROUTINE wf_numerov(quanta, energy, wf_matrix)
    !
    ! number of quanta
    INTEGER(KIND = I4B), INTENT(IN) :: quanta
    !
    ! eigenvalue
    REAL(KIND = DP), INTENT(IN) :: energy
    !
    ! eigenfunctions
    REAL(KIND = DP), DIMENSION(:,:), INTENT(OUT) :: wf_matrix
    !
    ! Local variables
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: gn, aux_vec, y_left, y_right, y_wf
    REAL(KIND = DP)  :: gn_min, gn_max, delta, nu_min, nu_min_plus, nu_min_minus, norm_con
    REAL(KIND = DP)  :: wf_der_right, wf_der_left
    INTEGER(KIND = I4B) :: index, ierr
    ! 
    ! Allocate vectors
    ALLOCATE(gn(1:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "gn allocation request denied."
       STOP
    ENDIF
    ALLOCATE(aux_vec(1:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "aux_vec allocation request denied."
       STOP
    ENDIF
    ALLOCATE(y_left(1:match_p + 2), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "y_left allocation request denied."
       STOP
    ENDIF
    ALLOCATE(y_right(match_p - 2:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "y_right allocation request denied."
       STOP
    ENDIF
    ALLOCATE(y_wf(1:npoints), STAT = ierr)    
    IF (ierr /= 0) THEN
       WRITE(*,*) "y_right allocation request denied."
       STOP
    ENDIF
    !
    ! Preliminary calculations
    !
    gn = (2.0_DP*amu*red_mass/(hbarc**2))*(energy - vpot)
    !
    aux_vec = 1.0_DP + ((x_step**2)/12.0_DP)*gn
    !
    gn = gn/aux_vec ! rescaled gn
    !
    gn_min = gn(1)
    gn_max = gn(npoints)
    !
    !delta = SQRT(x_step**2*gn_min*(x_step**2*gn_min-4.0_DP))
    delta = x_step*SQRT(gn_min*(x_step**2*gn_min-4.0_DP))
    !
    nu_min_plus = 0.5_DP*(x_step**2*gn_min + delta)
    nu_min_minus = 0.5_DP*(x_step**2*gn_min - delta)
    !
    ! Left branch  
    IF ( ABS(1.0_DP/(1.0_DP - nu_min_plus)) > 1 ) THEN
       nu_min = nu_min_plus
    ELSE
       nu_min = nu_min_minus
    ENDIF
    !
    y_left(1) = 1.0_DP
    y_left(2) = y_left(1)/(1.0_DP - nu_min)
    !
    y_left(1:2) = y_left(1:2)*aux_vec(1:2)  ! Renormalization
    !
    DO index = 3, match_p + 2
       !
       y_left(index) = 2.0_DP*y_left(index - 1) - y_left(index - 2) - y_left(index - 1)*gn(index-1)*x_step**2
       !
    ENDDO
    !
    y_left = y_left/aux_vec(1:match_p + 2)
    !
    !
    ! Right branch 
    IF ( ABS(1.0_DP - nu_min_plus) > 1 ) THEN
       nu_min = nu_min_plus
    ELSE
       nu_min = nu_min_minus
    ENDIF
    !
    IF (MOD(nodes_num,2) == 0) THEN
       y_right(npoints) = 1.0_DP ! even case
    ELSE
       y_right(npoints) = -ABS(1.0_DP) ! odd case
    ENDIF
    !
    y_right(npoints - 1) = y_right(npoints)*(1 - nu_min);
    !
    ! Rescaled wave function
    y_right(npoints - 1:npoints) = y_right(npoints - 1:npoints)*aux_vec(npoints - 1:npoints)  ! Renormalization
    !
    DO index = npoints - 2, match_p - 2, -1 ! Minus 2 for slope calculation
       !
       y_right(index) = 2*y_right(index + 1) - y_right(index + 2) - y_right(index + 1)*gn(index + 1)*x_step**2;
       !
    ENDDO
    !
    y_right = y_right / aux_vec(match_p - 2: npoints);
    !
    !
    ! Left and right derivatives at the matching point
    !
    wf_der_left = ( y_left(match_p + 1) - y_left(match_p - 1)  -&
         (1.0_DP/8.0_DP)*y_left(match_p + 2) + (1.0_DP/8.0_DP)*y_left(match_p - 2) ) /&
         (3.0_DP*x_step/2.0_DP) ! Left wf derivative
    !
    wf_der_right = ( y_right(match_p + 1) - y_right(match_p - 1)  -&
         (1.0_DP/8.0_DP)*y_right(match_p + 2) + (1.0_DP/8.0_DP)*y_right(match_p - 2) ) /&
         (3.0_DP*x_step/2.0_DP) ! Right wf derivative
    !
    ! Equate derivatives at the matching point
    y_right = (wf_der_left/wf_der_right)*y_right
    !
    ! Unnormalized wf
    y_wf(1:match_p) = y_left(1:match_p)
    y_wf(match_p+1:npoints) = y_right(match_p+1:npoints)
    !
    ! Renormalize
    y_wf = y_wf/y_wf(1)
    !
    ! Integrate with simpsn intlib subroutine
    CALL simpsn(npoints, x_step, y_wf**2, norm_con)
    y_wf = y_wf/SQRT(norm_con)
    !
    IF (i_wf_save > 0) THEN
       !
       IF (quanta == 0) THEN
          wf_matrix(1:npoints,1) = x_grid(1:npoints)
          wf_matrix(1:npoints,2) = y_wf(1:npoints)
       ELSE
          wf_matrix(1:npoints,quanta + 2) = y_wf(1:npoints)
       ENDIF
       !
    ENDIF
    !
  END SUBROUTINE wf_numerov
  !
  !
END MODULE numerov_alg
