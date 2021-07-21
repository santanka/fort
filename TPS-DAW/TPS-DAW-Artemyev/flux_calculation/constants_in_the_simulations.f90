module constants_in_the_simulations
  use lshell_setting

  implicit none

  !-------------------------------------
  ! simulation variation
  !-------------------------------------
  LOGICAL, PARAMETER :: frequency_sweep  = .false. ! .true. 
  LOGICAL, PARAMETER :: wave_existance   = .true. ! .true. 
  LOGICAL, PARAMETER :: categorization   = .false. ! .false.
  
  !-------------------------------------
  ! initial setting of simulation system
  !-------------------------------------
  INTEGER, PARAMETER          :: n_time = 80000 ! (10.9932 [s])
  INTEGER, PARAMETER          :: n_z = 3600 ! (n + 1) for dB_dz
  DOUBLE PRECISION, PARAMETER :: d_t = 1.0d0
  DOUBLE PRECISION, PARAMETER :: d_z = 0.5d0
  DOUBLE PRECISION, PARAMETER :: L_t = DBLE(n_time) * d_t
  DOUBLE PRECISION, PARAMETER :: L_z = DBLE(n_z) * d_z ! (n - 1) for dB_dz

  !-------------------------------------
  ! initial setting of wave
  !-------------------------------------
  DOUBLE PRECISION, PARAMETER :: wavelength_perp_eq = 150d3 / z_unit ![m]
  DOUBLE PRECISION, PARAMETER :: electrostatic_potential_0 = 200 / V_unit ![V]

  !!�s�g�p
  !DOUBLE PRECISION, PARAMETER :: w_p = 4d0
  !DOUBLE PRECISION, PARAMETER :: B_w = 10d-4
  !DOUBLE PRECISION, PARAMETER :: sweep_rate = 7.6d3 !fce_eq [kHz] 
  !DOUBLE PRECISION, PARAMETER :: freq_start = 0.3d0 ! [fce]
  !DOUBLE PRECISION, PARAMETER :: freq_end   = 0.45d0 ! [fce]
  !DOUBLE PRECISION, PARAMETER :: freq_sweep = freq_end - freq_start ! [fce]
  !DOUBLE PRECISION, PARAMETER :: t_element  = freq_sweep * fce_eq / sweep_rate / t_unit ! [t_unit]
  !DOUBLE PRECISION, PARAMETER :: dfdt = sweep_rate * t_unit / fce_eq
  !DOUBLE PRECISION, PARAMETER :: damping_region = 2000d0 ![c��^-1]
  !���Ԍo�߂ƂƂ���z>6000��freq��NaN�ɂȂ�B   
  
  !-------------------------------------
  ! initial setting of particle
  !-------------------------------------
  INTEGER            :: N_particle
  INTEGER, PARAMETER :: n_thread = 2

end module constants_in_the_simulations