module constants_in_the_simulations
  use lshell_setting

  implicit none

  !-------------------------------------
  ! simulation variation
  !-------------------------------------
  LOGICAL, PARAMETER :: wave_existence  = .true. ! .false. ! true or false
  LOGICAL, PARAMETER :: frequency_sweep = .true. ! .true. !
  LOGICAL, PARAMETER :: wave_packet     = .true. ! .true. !

  !-------------------------------------
  ! initial setting of simulation system
  !-------------------------------------
  INTEGER, PARAMETER          :: n_t = 90000!30000
  INTEGER, PARAMETER          :: n_z = 4000 ! (n + 1) for dB_dz
  DOUBLE PRECISION, PARAMETER :: d_t = 0.5d0
  DOUBLE PRECISION, PARAMETER :: d_z = 0.5d0

!  DOUBLE PRECISION, PARAMETER :: L_t = DBLE(n_t) * d_t
  DOUBLE PRECISION, PARAMETER :: L_z = DBLE(n_z) * d_z ! (n - 1) for dB_dz

  !-------------------------------------
  ! initial setting of wave
  !-------------------------------------
  DOUBLE PRECISION, PARAMETER :: w_p = 4d0
  DOUBLE PRECISION, PARAMETER :: B_w = 5d-4

  DOUBLE PRECISION, PARAMETER :: sweep_rate = fce_eq ! 7.6d3 ! [kHz]
  DOUBLE PRECISION, PARAMETER :: freq_start = 0.30d0 ! [fce]
  DOUBLE PRECISION, PARAMETER :: freq_end   = 0.45d0 ! [fce]
  DOUBLE PRECISION, PARAMETER :: freq_sweep = freq_end - freq_start ! [fce]
  DOUBLE PRECISION, PARAMETER :: t_element  = freq_sweep * fce_eq / sweep_rate / t_unit ! [t_unit]
  DOUBLE PRECISION, PARAMETER :: dfdt = sweep_rate * t_unit / fce_eq

  !-------------------------------------
  ! initial setting of particle
  !-------------------------------------
  DOUBLE PRECISION, PARAMETER :: gamma_input = 20.d0 ! [keV]

  INTEGER, PARAMETER :: N_p   = 100 ! (about)  
!  INTEGER, PARAMETER :: N_rdm = 100!19900 !  1720 -> 1M, 3450 -> 2M  
 
  DOUBLE PRECISION, PARAMETER :: alpha_max0 = 10d0 * deg2rad
  DOUBLE PRECISION, PARAMETER :: alpha_min0 = 0d0 * deg2rad 
  
  INTEGER, PARAMETER          :: n_alpha_output = 90 
  DOUBLE PRECISION, PARAMETER :: d_alpha_output = 90d0 / DBLE(n_alpha_output)

  INTEGER, PARAMETER          :: n_thread = 1
end module constants_in_the_simulations
