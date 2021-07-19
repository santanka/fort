module constant_parameter

  implicit none
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  double precision, parameter :: c  = 2.99792458d8
  double precision, parameter :: q  = 1.6021766208d-19
  double precision, parameter :: m  = 9.10938356d-31
  double precision, parameter :: m_e= m*c**2/(q*1.d3) ![keV] 1eV=1c
  double precision, parameter :: rad2deg = 180.d0/pi
  double precision, parameter :: deg2rad = pi/180.d0
   
end module constant_parameter
