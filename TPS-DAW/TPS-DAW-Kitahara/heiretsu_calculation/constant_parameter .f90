module constant_parameter
  implicit none

  !-------------------------------------
  ! mathematical and physical constants
  !-------------------------------------
  DOUBLE PRECISION, PARAMETER :: pi = 4d0*DATAN(1d0)
  DOUBLE PRECISION, PARAMETER :: c  = 299792458d0
  DOUBLE PRECISION, PARAMETER :: q  = 1.6021766208d-19
  DOUBLE PRECISION, PARAMETER :: m  = 9.10938356d-31

  DOUBLE PRECISION, PARAMETER :: m_e = m * c**2 / (q * 10d0**3) ! [keV] , 1 eV = 1 C

  DOUBLE PRECISION, PARAMETER :: rad2deg = 180d0 / pi
  DOUBLE PRECISION, PARAMETER :: deg2rad = pi / 180d0
  
end module constant_parameter
