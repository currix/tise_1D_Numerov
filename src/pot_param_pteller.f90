MODULE pot_param
  !
  USE nrtype
  USE global_vars
  !
  IMPLICIT NONE
  !
  ! MPARAM -> MAXIMUM NUMBER OF POTENTIAL PARAMETERS
  INTEGER(KIND = I4B), PARAMETER :: MPARAM = 5
  ! POTENTIAL PARAMETER VALUES AND UNITS
  REAL(KIND = DP), DIMENSION(1:MPARAM) :: param_pot
  !
CONTAINS
  !
  ELEMENTAL FUNCTION Potf(x)
    !
    !     POESCHL-TELLER 1D POTENTIAL
    !
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  Potf
    !
    !     POTENTIAL PARAMETERS
    !     PARAM_POT(1) --> VP0 DEPTH
    !     PARAM_POT(2) --> AP0 DIFFUSIVITY
    !
    REAL(KIND = DP) :: AUX
    !
    AUX = COSH(PARAM_pot(2)*X)
    POTf = PARAM_pot(1)/(AUX*AUX)
    !
  END FUNCTION Potf
  !  
  !
END MODULE pot_param

