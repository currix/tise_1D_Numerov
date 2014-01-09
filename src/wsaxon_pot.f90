ELEMENTAL FUNCTION Potf(x)
  !
  !     WOODS-SAXON 1D POTENTIAL
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
  !     PARAM_POT(1) --> VP0 DEPTH (MeV)
  !     PARAM_POT(2) --> RP  RADIUS (fm)
  !     PARAM_POT(3) --> AP0 DIFFUSIVITY (fm)
  !
  REAL(KIND = DP) :: AUX
  !
  IF (x > 0.0_DP) THEN
     AUX = EXP((PARAM_POT(2)-x)/PARAM_POT(3))
     POTf = PARAM_POT(1)*AUX/(1.0_DP+AUX)
  ELSE
     AUX = EXP((PARAM_POT(2)+x)/PARAM_POT(3))         
     POTf = PARAM_POT(1)*AUX/(1.0_DP+AUX)
  ENDIF
  !
  !
END FUNCTION Potf
