module lshell_setting

  use constant_parameter
  implicit none

  DOUBLE PRECISION, PARAMETER :: mu_0    = 4d0 * pi * 1d-7 !SI
  DOUBLE PRECISION, PARAMETER :: moment  = 7.75d22 ! the Earth's dipole moment model [Am^2]
  DOUBLE PRECISION, PARAMETER :: R_E     = 6371d3  ! radius of the Earth [m]

  DOUBLE PRECISION, PARAMETER :: L         = 9d0 ! L-shell
  DOUBLE PRECISION, PARAMETER :: B0_eq     = (1d-7 * moment) / (L * R_E)**3 * 1d4 ![T]→[G]
  DOUBLE PRECISION, PARAMETER :: Omega0_eq = q * B0_eq / m / c ![rad/s]
  DOUBLE PRECISION, PARAMETER :: fce_eq    = Omega0_eq / (2d0 * pi) ![s^-1]

  DOUBLE PRECISION, PARAMETER :: z_unit = c / Omega0_eq ![cm]
  DOUBLE PRECISION, PARAMETER :: t_unit = 1d0 / Omega0_eq ![s]
  DOUBLE PRECISION, PARAMETER :: r_eq   = (R_E * L) * 1d2 / z_unit ![cm]→[]
  DOUBLE PRECISION, PARAMETER :: e_unit = q ![statC]
  DOUBLE PRECISION, PARAMETER :: J_unit = m * c**2d0 ![erg]
  DOUBLE PRECISION, PARAMETER :: V_unit = m * c**2d0 / q ![statV]
  DOUBLE PRECISION, PARAMETER :: charge = q / e_unit ![statC]→[]
  DOUBLE PRECISION, PARAMETER :: c_normal = c / (z_unit/t_unit)
  DOUBLE PRECISION, PARAMETER :: ion_mass = 1.672621898d-27 * 1d3 / m ![kg]→[g]→[]
  DOUBLE PRECISION, PARAMETER :: electron_mass = m / m ![g]→[]

  DOUBLE PRECISION, PARAMETER :: number_density_eq = 1d0 * z_unit**3d0 ![cm^-3]→[]
  DOUBLE PRECISION, PARAMETER :: Temperature_ion = 1000 * (q*c/1d1) * 1d7 / J_unit ![eV]→[erg]→[]
  DOUBLE PRECISION, PARAMETER :: Temperature_electron = 100 * (q*c/1d1) * 1d7 / J_unit ![eV]→[erg]→[]



end module lshell_setting