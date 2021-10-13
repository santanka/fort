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
  INTEGER, PARAMETER          :: n_time = 20000  ! !80000 (10.9932 [s])
  INTEGER, PARAMETER          :: n_z = 3000 ! (n + 1) for dB_dz
  DOUBLE PRECISION, PARAMETER :: d_t = 1.0d0 / Omega0_eq
  DOUBLE PRECISION, PARAMETER :: d_z = 0.5d0 / Omega0_eq * c_normal
  DOUBLE PRECISION, PARAMETER :: L_t = DBLE(n_time) * d_t
  DOUBLE PRECISION, PARAMETER :: L_z = DBLE(n_z) * d_z ! (n - 1) for dB_dz

  !-------------------------------------
  ! initial setting of wave
  !-------------------------------------
  DOUBLE PRECISION, PARAMETER :: electrostatic_potential_0 = 0d0 * 1d8 / c / V_unit ![V]?[statV]?[]

  !-------------------------------------
  ! initial setting of particle
  !-------------------------------------
  INTEGER            :: N_particle
  INTEGER, PARAMETER :: n_thread = 2

end module constants_in_the_simulations
