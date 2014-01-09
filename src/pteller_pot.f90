ELEMENTAL FUNCTION Potf(x)
  !
  !     POESCHL-TELLER 1D POTENTIAL
  !
  USE nrtype
  USE constants
  USE pot_param
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
  !     PARAM_POT(2) --> RP  RADIUS
  !     PARAM_POT(3) --> AP0 DIFFUSIVITY
  !
  REAL(KIND = DP) :: AUX
  !
  AUX = COSH(PARAM_pot(2)*X)
  POTf = PARAM_pot(1)/(AUX*AUX)
  !
END FUNCTION Potf
