subroutine z_position_to_radius_MLAT(z_position, radius, MLAT)
    use lshell_setting

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position
    DOUBLE PRECISION, INTENT(OUT) :: radius, MLAT
    DOUBLE PRECISION :: ff, gg, MLAT0, MLAT1
    INTEGER :: ii

    MLAT0 = 1d0
    do ii = 1, 1000000
        if (ii == 1000000) then
            print *, "Error!: solution is not found. z_position = ", z_position
        endif

        ff = r_eq * ((1d0 / 2d0) * DSIN(MLAT0) * DSQRT(3d0 * DSIN(MLAT0)**2 + 1d0) &
            & + (1d0 / (2d0 * DSQRT(3d0))) * DLOG(DSQRT(3d0) * DSIN(MLAT0) &
            & + DSQRT(3d0 * DSIN(MLAT0)**2 + 1d0))) &
            & - z_position
        gg = r_eq * DCOS(MLAT0) * DSQRT(3d0 * DSIN(MLAT0)**2 + 1d0)

        MLAT1 = MLAT0 - ff / gg
        if (DABS(MLAT1 - MLAT0) <= 1d-5) exit
        MLAT0 = MLAT1
    end do !ii

    MLAT = MLAT1
    radius = r_eq * DCOS(MLAT)**2d0

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_BB(z_position, BB)

    implicit none
      
    DOUBLE PRECISION, INTENT(IN)  :: z_position
    DOUBLE PRECISION, INTENT(OUT) :: BB
    DOUBLE PRECISION :: radius, MLAT
    
    call z_position_to_radius_MLAT(z_position, radius, MLAT)
  
    BB = DSQRT(1d0 + 3d0 * DSIN(MLAT)**2) / DCOS(MLAT)**6

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_number_density(z_position, number_density)
    use constants_in_the_simulations

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position
    DOUBLE PRECISION, INTENT(OUT) :: number_density
    DOUBLE PRECISION :: radius, MLAT

    call z_position_to_radius_MLAT(z_position, radius, MLAT)

    number_density = number_density_eq * DCOS(MLAT)**(-5d0)

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_alfven_velocity(z_position, BB, alfven_velocity)
    use constants_in_the_simulations
    use lshell_setting

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position, BB
    DOUBLE PRECISION, INTENT(OUT) :: alfven_velocity
    DOUBLE PRECISION :: number_density

    CALL z_position_to_number_density(z_position, number_density)

    alfven_velocity = BB / DSQRT(mu_0 * number_density * ion_mass)

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_ion_acoustic_gyroradius(z_position, BB, ion_acoustic_gyroradius)
    use constant_parameter
    use lshell_setting

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position, BB
    DOUBLE PRECISION, INTENT(OUT) :: ion_acoustic_gyroradius
    DOUBLE PRECISION :: number_density
    
    CALL z_position_to_number_density(z_position, number_density)

    ion_acoustic_gyroradius = DSQRT((1.5d0 * Temperature_ion + 2d0 * Temperature_electron) * ion_mass) &
        & / charge / BB

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_wave_frequency(z_position, wave_frequency)
    use constant_parameter
    use lshell_setting

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position
    DOUBLE PRECISION, INTENT(OUT) :: wave_frequency
    DOUBLE PRECISION :: ss, alfven_velocity_eq

    ss = z_position / r_eq
    CALL z_position_to_alfven_velocity(0d0, 1d0, alfven_velocity_eq)
    wave_frequency = pi / z_position * alfven_velocity_eq

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine BB_to_wave_number_perp(BB, wave_number_perp)
    use constant_parameter
    use constants_in_the_simulations

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: BB
    DOUBLE PRECISION, INTENT(OUT) :: wave_number_perp

    wave_number_perp = 2d0 * pi / wavelength_perp_eq * DSQRT(BB)

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_wave_number_para(z_position, BB, wave_number_perp, wave_number_para)
    use lshell_setting

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position, BB, wave_number_perp
    DOUBLE PRECISION, INTENT(OUT) :: wave_number_para
    DOUBLE PRECISION :: wave_frequency, alfven_velocity, ion_acoustic_gyroradius

    CALL z_position_to_wave_frequency(z_position, wave_frequency)
    CALL z_position_to_alfven_velocity(z_position, BB, alfven_velocity)
    CALL z_position_to_ion_acoustic_gyroradius(z_position, BB, ion_acoustic_gyroradius)

    wave_number_para = wave_frequency / alfven_velocity / &
        & DSQRT(1d0 + wave_number_perp**2d0 * ion_acoustic_gyroradius**2d0 * (1d0 + Temperature_ion / Temperature_electron))

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine wave_number_para_to_wave_phase_initial(wave_number_para_pre, wave_number_para, wave_phase_pre, wave_phase)
    use constants_in_the_simulations

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: wave_number_para_pre, wave_number_para, wave_phase_pre
    DOUBLE PRECISION, INTENT(OUT) :: wave_phase

    wave_phase = wave_phase_pre + (wave_number_para_pre + wave_number_para) / 2d0 * d_z

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_position_to_electrostatic_potential(z_position, BB, electrostatic_potential)
    use constant_parameter
    use lshell_setting
    use constants_in_the_simulations

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position, BB
    DOUBLE PRECISION, INTENT(OUT) :: electrostatic_potential
    DOUBLE PRECISION :: radius, MLAT, wave_number_perp, ion_acoustic_gyroradius, g_function

    CALL z_position_to_radius_MLAT(z_position, radius, MLAT)
    CALL BB_to_wave_number_perp(BB, wave_number_perp)
    CALL z_position_to_ion_acoustic_gyroradius(z_position, BB, ion_acoustic_gyroradius)

    g_function = 1d0 / 2d0 * (DTANH(360d0 / pi * DABS(MLAT) - 5d0) - DTANH(-5d0))

    electrostatic_potential = - wave_number_perp**2d0 * ion_acoustic_gyroradius**2d0 &
        & * (1d0 + Temperature_ion / Temperature_electron) * electrostatic_potential_0 * g_function
    
end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine time_to_wave_phase_update(wave_phase, wave_frequency)
    use constants_in_the_simulations

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: wave_frequency
    DOUBLE PRECISION, INTENT(INOUT) :: wave_phase

    wave_phase = wave_phase - wave_frequency * d_t

end subroutine

!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_particle_to_position(z_particle, z_position, i_z_left, i_z_right, ratio)
    use constants_in_the_simulations, only: n_z

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position(-n_z:n_z), z_particle
    INTEGER, INTENT(OUT) :: i_z_left, i_z_right
    DOUBLE PRECISION, INTENT(OUT) :: ratio
    DOUBLE PRECISION :: difference(-n_z:n_z), difference_min
    INTEGER :: i_min(1)

    difference = DABS(z_position - z_particle)
    i_min = MINLOC(difference) - (n_z + 1) !array_difference : 0~2*n_z+1

    if (i_min(1) >= n_z) then
        i_z_left = n_z - 1
        i_z_right = n_z
        ratio = 1d0

    else if (i_min(1) <= -n_z) then
        i_z_left = n_z - 1
        i_z_right = n_z
        ratio = 1d0

    else
        difference_min = z_position(i_min(1)) - z_particle
        if (difference_min > 0) then
            i_z_left = i_min(1) - 1
            i_z_right = i_min(1)
        else if (difference_min <= 0) then
            i_z_left = i_min(1)
            i_z_right = i_min(1) + 1
        end if

        ratio = (z_particle - z_position(i_z_left)) / (z_position(i_z_right) - z_position(i_z_left))

    end if

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine u_particle_to_gamma(u_particle, gamma)

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: u_particle(0:2)
    DOUBLE PRECISION, INTENT(OUT) :: gamma

    gamma = DSQRT(1d0 + u_particle(0)**2d0 + u_particle(1)**2d0)

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine u_particle_to_v_particle_para(u_particle, v_particle_para)

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: u_particle(0:2)
    DOUBLE PRECISION, INTENT(OUT) :: v_particle_para
    DOUBLE PRECISION :: gamma

    CALL u_particle_to_gamma(u_particle, gamma)

    v_particle_para = u_particle(0) / gamma

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine z_particle_to_dB_dz(z_particle, dB_dz)

    use lshell_setting, only: r_eq

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_particle
    DOUBLE PRECISION, INTENT(OUT) :: dB_dz
    DOUBLE PRECISION :: radius, MLAT

    CALL z_position_to_radius_MLAT(z_particle, radius, MLAT)

    dB_dz = 3d0 * DSIN(MLAT) * (5d0 * DSIN(MLAT)**2d0 + 3d0) / DCOS(MLAT)**8d0 / (3d0 * DSIN(MLAT)**2d0 + 1d0) / r_eq

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine electrostatic_potential_to_EE_wave_para(electrostatic_potential, wave_number_para, wave_phase, EE_wave_para)

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: electrostatic_potential, wave_number_para, wave_phase
    DOUBLE PRECISION, INTENT(OUT) :: EE_wave_para

    EE_wave_para = - wave_number_para * electrostatic_potential * DCOS(wave_phase)

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine Motion_of_Equation(z_position, wave_phase, z_p, u_p, force)
    !p -> particle

    use constant_parameter, only: pi
    use lshell_setting, only: charge
    use constants_in_the_simulations, only: n_z

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_position(-n_z:n_z), wave_phase(-n_z:n_z)
    DOUBLE PRECISION, INTENT(IN) :: z_p, u_p(0:2)
    DOUBLE PRECISION, INTENT(OUT) :: force(0:2)
    DOUBLE PRECISION :: gamma, ratio, BB_p, dB_dz_p, wave_number_perp_p, wave_number_para_p, force_wave(0:2)
    DOUBLE PRECISION :: electrostatic_potential_p, wave_phase_p, EE_wave_para_p
    INTEGER :: i_z_left, i_z_right

    CALL u_particle_to_gamma(u_p, gamma)
    CALL z_particle_to_position(z_p, z_position, i_z_left, i_z_right, ratio)
    CALL z_position_to_BB(z_p, BB_p)
    CALL z_particle_to_dB_dz(z_p, dB_dz_p)
    CALL BB_to_wave_number_perp(BB_p, wave_number_perp_p)
    CALL z_position_to_wave_number_para(z_p, BB_p, wave_number_perp_p, wave_number_para_p)
    CALL z_position_to_electrostatic_potential(z_p, BB_p, electrostatic_potential_p)

    wave_phase_p = (1d0 - ratio) * wave_phase(i_z_left) + ratio * wave_phase(i_z_right)

    CALL electrostatic_potential_to_EE_wave_para(electrostatic_potential_p, wave_number_para_p, wave_phase_p, EE_wave_para_p)

    !force(wave)
    force_wave(0) = - charge * EE_wave_para_p
    force_wave(1) = 0d0
    force_wave(2) = 0d0

    !force
    force(0) = - u_p(1)**2d0 / 2d0 / BB_p / gamma * dB_dz_p + force_wave(0)
    force(1) = u_p(0) * u_p(1) / 2d0 / BB_p / gamma * dB_dz_p + force_wave(1)
    force(2) = charge * BB_p / gamma

end subroutine
!
!!----------------------------------------------------------------------------------------------------------------------------------
!
subroutine particle_update_by_runge_kutta(z_in, wave_phase_in, z_particle, u_particle, equator_flag, edge_flag)

    use constants_in_the_simulations, only: d_t, n_z, L_z, d_z

    implicit none

    DOUBLE PRECISION, INTENT(IN) :: z_in(-n_z:n_z), wave_phase_in(-n_z:n_z)
    DOUBLE PRECISION, INTENT(INOUT) :: z_particle, u_particle(0:2)
    INTEGER, INTENT(OUT) :: equator_flag
    INTEGER, INTENT(OUT) :: edge_flag
    DOUBLE PRECISION :: ff_RK_1(0:2), ff_RK_2(0:2), ff_RK_3(0:2), ff_RK_4(0:2), u_particle_s(0:2)
    DOUBLE PRECISION :: kk_RK_1, kk_RK_2, kk_RK_3, kk_RK_4
    

    u_particle_s(:) = u_particle(:)

    !RK4
    CALL u_particle_to_v_particle_para(u_particle_s, kk_RK_1)
    CALL Motion_of_Equation(z_in, wave_phase_in, z_particle, u_particle, ff_RK_1)

    CALL u_particle_to_v_particle_para(u_particle_s + ff_RK_1 / 2d0 * d_t, kk_RK_2)
    CALL Motion_of_Equation(z_in, wave_phase_in, z_particle + kk_RK_1 / 2d0 * d_t, u_particle_s + ff_RK_1 / 2d0 * d_t, ff_RK_2)

    CALL u_particle_to_v_particle_para(u_particle_s + ff_RK_2 / 2d0 * d_t, kk_RK_3)
    CALL Motion_of_Equation(z_in, wave_phase_in, z_particle + kk_RK_2 / 2d0 * d_t, u_particle_s + ff_RK_2 / 2d0 * d_t, ff_RK_3)

    CALL u_particle_to_v_particle_para(u_particle_s + ff_RK_3 * d_t, kk_RK_4)
    CALL Motion_of_Equation(z_in, wave_phase_in, z_particle + kk_RK_3 * d_t, u_particle_s + ff_RK_3 * d_t, ff_RK_4)

    !particle update
    u_particle(:) = u_particle(:) + (ff_RK_1(:) + 2d0 * ff_RK_2(:) + 2d0 * ff_RK_3(:) + ff_RK_4(:)) * d_t / 6d0
    z_particle = z_particle + (kk_RK_1 + 2d0 * kk_RK_2 + 2d0 * kk_RK_3 + kk_RK_4) * d_t / 6d0
    
    if (z_particle <= 0d0) then
        z_particle = - z_particle
        u_particle(0) = DABS(u_particle(0))
        equator_flag = 1

    else if (z_particle >= L_z) then !mirror
        z_particle = L_z - (z_particle - L_z)
        u_particle(0) = - DABS(u_particle(0))
        edge_flag = 1

    end if

end subroutine
    




        

