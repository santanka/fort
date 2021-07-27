program main
    !$ use omp_lib

    use mt19937
    use constant_parameter
    use lshell_setting
    use variables

    implicit none

!--MPI----------------------------------------------------------------------
    include 'mpif.h'
    INTEGER(KIND=4) ierr, nprocs, myrank
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
!--MPI----------------------------------------------------------------------

    call write_time(string)

    write(command, '(A13, I3.3)') 'mkdir results', myrank
    call system(command)

    !------------------------
    !initial setting of field
    !------------------------
    do i_z = -n_z, n_z
        z_position(i_z) = DBLE(i_z) * d_z
        call z_position_to_BB(z_position(i_z), BB(i_z))
        call z_position_to_wave_frequency(z_position(i_z), wave_frequency(i_z))
        call BB_to_wave_number_perp(BB(i_z), wave_number_perp(i_z))
        call z_position_to_wave_number_para(z_position(i_z), BB(i_z), wave_number_perp(i_z), wave_number_para(i_z))

        !initial wave_phase profile
        if (i_z == -n_z) then
            wave_phase(i_z) = 0d0
        else
            CALL wave_number_para_to_wave_phase_initial(wave_number_para(i_z - 1), wave_number_para(i_z), & 
                &wave_phase(i_z -1), wave_phase(i_z))
        end if
        
    end do !i_z


    !----------------------------
    !initial setting of particles
    !----------------------------
    CALL system_clock(count=clock)
    clock = 4267529
    CALL sgrnd(clock)
    print *, clock

    !count the quantity of the data
    WRITE(file_data, '(A17, I3.3, A4)') 'initial_condition', myrank, '.dat'
    OPEN (500, file = file_data)
    N_particle = 0
    do !s
        read(500, *, end = 99)
        N_particle = N_particle + 1
    end do !s

    99 N_particle = N_particle - 2
    print *, 'N_particle = ', N_particle
    close(500)

    !allocate the capacity
    allocate(alpha0(1:N_particle))
    allocate(gamma0(1:N_particle))
    allocate(energy0(1:N_particle))
    allocate(alpha_eq(1:N_particle))
    allocate(z_particle(1:N_particle))
    allocate(u_particle(0:2, 1:N_particle))
    allocate(equator_time(1:N_particle))
    allocate(equator_flag(1:N_particle))
    allocate(wave_flag(1:N_particle))
    allocate(edge_flag(1:N_particle))

    
    !get the values
    write(file_data,'(A17, I3.3, A4)') 'initial_condition', myrank, '.dat'
    open (500, file = file_data)
  
    do i = 1, 1
       read(500,*)
    end do !i
  
    do i = 1, N_particle
       read(500,*,iostat=ios) energy0(i), alpha0(i), z_particle(i), alpha_eq(i) 
       
       gamma0(i) = DBLE(energy0(i)) / m_e + 1d0     
       v_particle = DSQRT(1d0 - 1d0/gamma0(i)**2)
       v_particle_para = v_particle * DCOS(alpha0(i)*deg2rad)
       v_particle_perp = v_particle * DSIN(alpha0(i)*deg2rad)
       u_particle(0, i) = gamma0(i) * v_particle_para
       u_particle(1, i) = gamma0(i) * v_particle_perp
       u_particle(2, i) = 2d0 * pi * grnd() 
    end do !i

    !flag & sign reset
    equator_flag = 0
    equator_time = 0
    wave_flag = 0
    edge_flag = 0


    !-----------------------
    !initial setting of wave
    !-----------------------
    if (wave_existance == .true.) then
        do i_z = -n_z, n_z
            CALL z_position_to_electrostatic_potential(z_position(i_z), BB(i_z), electrostatic_potential(i_z))
        end do !i_z
    else
        electrostatic_potential = 0d0
    end if
    

    !---------
    !file open
    !---------
    do i_thr = 0, N_thr
        N_file = 20 + i_thr
        WRITE(file_equator, '(A7, I3.3, A17, I2.2, A4)') 'results', myrank, '/count_at_equator', i_thr, '.dat'
        OPEN(unit = N_file, file = file_equator)
    end do !i_thr


    !----------------
    !simulation start
    !----------------
    do i_time = 1, n_time
        time = DBLE(i_t) * d_t


        !-----------------
        !update wave_phase
        !-----------------
        do i_z = -n_z, n_z
            CALL time_to_wave_phase_update(wave_phase(i_z), wave_frequency(i_z))
        end do !i_z

        !$omp parallel num_threads(n_thread) &
        !$omp & private(i_particle, i_time, z_particle_sim, u_particle_sim, equator_flag_sim, wave_flag_sim, edge_flag_sim, &
        !$omp & N_file, alpha_particle_eq, energy_particle, BB_particle, alpha, v_particle, v_particle_para, v_particle_perp)


        do i_particle = 1, N_particle
            z_particle_sim = z_particle(i_particle)
            u_particle_sim(:) = u_particle(:, i_particle)
            equator_flag_sim = equator_flag(i_particle)
            wave_flag_sim = wave_flag(i_particle)
            edge_flag_sim = edge_flag(i_particle)


            N_file = 20 + omp_get_thread_num()

            CALL particle_update_by_runge_kutta(z_position, wave_phase, z_particle_sim, &
                & u_particle_sim, equator_flag_sim, edge_flag_sim)
            
            if (equator_flag_sim == 1) then
                CALL u_particle_to_energy(u_particle_sim, energy_particle)
                CALL u_particle_to_alpha(z_particle_sim, u_particle_sim, alpha_particle_eq)
                WRITE(unit = N_file, fmt = '(5E15.7, 4I3)') time, alpha_particle_eq, energy_particle, alpha_eq(i_particle), &
                    & energy0(i_particle), wave_flag(i_particle), cross_theta_0(i_particle), Cw_flag(i_particle), S_flag(i_particle)
            
                z_particle_sim = grnd() * ABS(z_particle_sim - L_z) + z_particle_sim
                CALL z_position_to_BB(z_particle_sim, BB_particle)
                alpha = ASIN(SQRT(BB_particle * SIN(alpha_eq(i_particle) * deg2rad)**2d0))
                v_particle = SQRT(1d0 - 1d0 / gamma0(i_particle)**2d0)
                v_particle_para = v_particle * COS(alpha)
                v_particle_perp = v_particle * SIN(alpha)
                u_particle_sim(0) = gamma0(i_particle) * v_particle_para
                u_particle_sim(1) = gamma0(i_particle) * v_particle_perp
                u_particle_sim(2) = 2d0 * pi * grnd()
                equator_flag_sim = 0
                wave_flag_sim = 0               
            end if

            if (edge_flag_sim == 1) then
                z_particle_sim = grnd() * ABS(z_particle_sim - L_z) + z_particle_sim
                CALL z_position_to_BB(z_particle_sim, BB_particle)
                alpha = ASIN(SQRT(BB_particle * SIN(alpha_eq(i_particle) * deg2rad)**2d0))
                v_particle = SQRT(1d0 - 1d0 / gamma0(i_particle)**2d0)
                v_particle_para = - v_particle * COS(alpha)
                v_particle_perp = v_particle * SIN(alpha)
                u_particle_sim(0) = gamma0(i_particle) * v_particle_para
                u_particle_sim(1) = gamma0(i_particle) * v_particle_perp
                u_particle_sim(2) = 2d0 * pi * grnd()
                edge_flag_sim = 0
            end if

            z_particle(i_particle) = z_particle_sim
            u_particle(:, i_particle) = u_particle_sim(:)
            wave_flag(i_particle) = wave_flag_sim
            edge_flag(i_particle) = edge_flag_sim

        end do !i_particle
        !$omp end do nowait
        !$omp end parallel

        !-------
        !out put
        !-------
        if (mod(i_time, 1000) == 0) then
            WRITE(string, '(F6.2, A2, A6, I2.2)') DBLE(i_time) / DBLE(n_time) * 100d0, ' %'
            CALL write_time(string)
        end if
            
    end do !i_time

    CLOSE(200)
    print *, "end"


    !-MPI------------------
    CALL MPI_FINALIZE(ierr)
    !-MPI------------------

    
end program main