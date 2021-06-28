module variables
  use lshell_setting
  use constants_in_the_simulations

  implicit none

  !-------------------------------------
  ! variables
  !-------------------------------------

  INTEGER          :: i_t, i_z, i_p, i_phase, i_alpha, i_grid, i_thr
  INTEGER          :: i, j, ios, a, E, i_v_para, i_v_perp
  INTEGER          :: n_file
  INTEGER, PARAMETER :: n_thr = 20
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: z(-n_z : n_z)
  DOUBLE PRECISION :: B(-n_z : n_z)
  DOUBLE PRECISION :: k_init(-n_z : n_z)
  DOUBLE PRECISION :: phase(-n_z : n_z)
  DOUBLE PRECISION :: rnd
  CHARACTER(64)    :: file_output, file_particle, file_wave, file_data
  CHARACTER(64)    :: file_energy, file_alpha, file_distribution, file_phase_space
  CHARACTER(64)    :: file_equator
  CHARACTER(20)    :: string
  CHARACTER(64) :: command  
  
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

  INTEGER          :: clock  !, N_p_real, i_rdm
  DOUBLE PRECISION, allocatable :: alpha0(:), gamma0(:), energy0(:), alpha_eq(:)
  DOUBLE PRECISION :: alpha, gamma, energy, B0_p
  INTEGER          :: alpha_loop, energy_loop, phi_loop
  DOUBLE PRECISION :: v, zeta
  DOUBLE PRECISION :: v_para, v_perp
  DOUBLE PRECISION :: v_0, v_1
  DOUBLE PRECISION, allocatable :: z_p(:), u_p(:,:), u_p_eq(:,:), v_eq(:,:)
  DOUBLE PRECISION :: z_p_sim, u_p_sim(0:2)
  DOUBLE PRECISION,allocatable :: equator_time(:)
  INTEGER,allocatable          :: equator_flag(:), wave_flag(:), edge_flag(:)
  INTEGER :: equator_flag_sim, wave_flag_sim, edge_flag_sim 

  !------------------------------
  ! decide z
  !------------------------------  !now
  INTEGER          :: n_z_mirror 
  DOUBLE PRECISION :: B_mirror, z_mirror, norm, cumulative(0:n_z)

  
  !-------------------------------------
  ! categorization 
  !-------------------------------------

  DOUBLE PRECISION, allocatable :: sign_theta0(:), sign_theta1(:)
  DOUBLE PRECISION :: B_p, alpha_p, energy_p, zeta_p
  DOUBLE PRECISION :: freq_p, ampl_p, k_p, V_g_p
  DOUBLE PRECISION :: gamma_p, V_R_p, theta_p, Cw_p, w_tr_p, dB_dz_p, dk_dB_p, S_p
  INTEGER, allocatable :: cross_theta_0(:)
  INTEGER, allocatable :: Cw_flag(:), S_flag(:)
  INTEGER :: Cw_flag_sim, S_flag_sim, cross_theta_0_sim
  DOUBLE PRECISION :: sign_theta0_sim, sign_theta1_sim
  
end module variables
