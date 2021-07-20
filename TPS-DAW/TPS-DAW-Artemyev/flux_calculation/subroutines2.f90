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
    DOUBLE PRECISION :: ss

    ss = z_position / r_eq
    wave_frequency = 2d0 * pi / 2d0 / ss 

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