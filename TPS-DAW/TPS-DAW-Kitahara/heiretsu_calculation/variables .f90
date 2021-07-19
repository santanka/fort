module variables
  use lshell_setting
  use constants_in_the_simulations

  implicit none

  !-------------------------------------
  ! variables
  !-------------------------------------

  INTEGER          :: i_t, i_z, i_p, i_phase, i_alpha, i_grid
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: z(0 : n_z)
  DOUBLE PRECISION :: B(0 : n_z)
  DOUBLE PRECISION :: k_init(0 : n_z)
  DOUBLE PRECISION :: phase(0 : n_z)
  DOUBLE PRECISION :: rnd
  
  INTEGER          :: count_p(0 : n_z), count, total, progress
  DOUBLE PRECISION :: loss_cone(0 : n_z), r0, lambda0

  CHARACTER(64)    :: file_particle, file_wave, file_wave_front, file_equator, file_position, file_losscone
  CHARACTER(20)    :: string

  !-------------------------------------
  ! for wave
  !-------------------------------------

  DOUBLE PRECISION :: freq_0(-n_z : n_z), freq_1(-n_z : n_z)
  DOUBLE PRECISION :: ampl_0(-n_z : n_z), ampl_1(-n_z : n_z)
  DOUBLE PRECISION :: V_g(-n_z : n_z), V_g_0
  DOUBLE PRECISION :: z_front, B_front, V_g_front
  DOUBLE PRECISION :: z_edge,  B_edge,  V_g_edge

  !-------------------------------------
  ! for particle
  !-------------------------------------

  INTEGER          :: clock!, N_p_real, i_rdm

  DOUBLE PRECISION :: v0, v_para0, v_perp0
  DOUBLE PRECISION :: alpha0, gamma0, phase0
  DOUBLE PRECISION :: alpha, gamma, energy, B0_p

  DOUBLE PRECISION :: z_p(1:n_p), u_p(0:2, 1:n_p), zeta
  DOUBLE PRECISION :: z_p_sim, u_p_sim(0:2)
  LOGICAL          :: equator_flag

  INTEGER, PARAMETER :: n_thr = 20
  INTEGER          :: n_file, i_thr

  DOUBLE PRECISION :: alpha_min, alpha_max

  !-------------------------------------
  ! for resonance condition
  !-------------------------------------

  INTEGER          :: i_w
  DOUBLE PRECISION :: w, alpha_c1, alpha_c2, alpha_R, V_R, V_para

  !-------------------------------------
  ! output
  !-------------------------------------

  INTEGER          :: i_alpha_output, i_input
  DOUBLE PRECISION :: B_output, alpha_output, count_output(1:n_alpha_output), alpha_p
  
  !------------------------------
  ! decide z
  !------------------------------  !now
  INTEGER          :: n_z_mirror 
  DOUBLE PRECISION :: B_mirror, z_mirror, norm, cumulative(0:n_z), B_p
  
contains

  subroutine initial_display
    !---------------------------------
    ! initial setting of wave (frequency, and sweep rate)
    !---------------------------------

    write(*, *) '!---------------------------------'
    write(*, *) '! initial setting of wave'
    write(*, *) '!---------------------------------'
    write(*, *) ' '
    write(*, *) 'We_eq  [rad/s] = ', Omega0_eq
    write(*, *) 'fce_eq    [Hz] = ', fce_eq
    write(*, *) 'freq_str [fce] = ', freq_start
    write(*, *) 'freq_end [fce] = ', freq_end
    write(*, *) 'delta_t [t_unit] = ', t_element
    write(*, *) 'one_sec [t_step] = ', 1/t_unit
    write(*, *) ' '
    write(*, *) 'sweeprate [Hz/s] = ', sweep_rate
    write(*, *) 'sweeprate [unit] = ', dfdt
    write(*, *) ' '
    write(*, *) ' '

    !-------------------------------------
    ! initial setting of particles
    !-------------------------------------

    gamma0 = gamma_input / m_e + 1d0
    call gamma_to_vabs(gamma0, v0)
    v0      = SQRT(1d0 - 1d0 / gamma0**2)

    write(*, *) '!---------------------------------'
    write(*, *) '! initial setting of particles'
    write(*, *) '!---------------------------------'
    write(*, *) ' '
    write(*, *) 'energy  = ', gamma_input, '[keV]'
    write(*, *) ' '
    write(*, *) 'gamma0  = ', gamma0
    write(*, *) 'v0      = ', v0
    write(*, *) ' '

  end subroutine initial_display
end module variables
