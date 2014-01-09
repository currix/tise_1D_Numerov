ELEMENTAL FUNCTION Potf(x)
  !
  !     MORSE 1D POTENTIAL
  !
  !     Note that x = R - R_e
  !
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
  !
  !     POTENTIAL PARAMETERS
  !     PARAM_POT(1) --> VM DEPTH
  !     PARAM_POT(2) --> AM DIFFUSIVITY
  !
  REAL(KIND = DP) :: AUX
  !
  AUX = EXP(-PARAM_pot(2)*X)
  POTf = PARAM_pot(1)*(AUX*AUX - 2.0_DP*AUX)
  !
END FUNCTION Potf
