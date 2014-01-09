  MODULE phys_constants
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! PHYSICAL CONSTANTS
    REAL(KIND = DP), PARAMETER :: hbarc = 197.32858E0_DP  ! MeV fm
    REAL(KIND = DP), PARAMETER :: amu = 938.92635E0_DP ! MeV/c^2 / amu
    REAL(KIND = DP), PARAMETER :: H2OM = 41.4713768_DP ! MeV/fm^2 (hbar^2/1 amu)
    !
  END MODULE phys_constants
