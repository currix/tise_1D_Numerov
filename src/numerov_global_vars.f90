  MODULE global_vars
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! DEFINITION OF GLOBAL VARIABLES
    !
    REAL(KIND = DP) :: x_min, x_max, x_step  ! fm
    !
    INTEGER(KIND = I4B) :: npoints ! x_grid dimension
    !
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: x_grid ! fm
    !
    INTEGER(KIND = I4B) :: iprint ! verbosity flag
    !
    !
    ! DEFINITION OF SYSTEM VARIABLES
    !
    REAL(KIND = DP) :: red_mass
    !
  END MODULE global_vars
