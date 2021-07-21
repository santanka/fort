module lshell_setting

  use constant_parameter
  implicit none

  double precision, parameter :: mu_0 =4.d0*pi*1.d-7
  double precision, parameter :: Earths_dipole_moment =8.1d22 
  double precision, parameter :: R_E =6.371d6
  double precision, parameter :: L =6.d0 !L-shell
  double precision, parameter :: B0_eq =(1.d-7*Earths_dipole_moment)/(L*R_E)**3
  double precision, parameter :: Omega0_eq =q*B0_eq/m
  double precision, parameter :: z_unit =c/Omega0_eq
  double precision, parameter :: r_eq =R_E*L/z_unit
  

end module lshell_setting