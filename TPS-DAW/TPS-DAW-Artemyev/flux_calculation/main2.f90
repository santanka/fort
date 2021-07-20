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
            CALL wave_number_para_to_wave_phase_initial &
                & (wave_number_para(i_z - 1), wave_number_para(i_z), wave_phase(i_z -1), wave_phase(i_z))
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
    N_p = 0
    do !s
        read(500, *, end = 99)
        N_p = N_p + 1
    end do !s

    99 N_p = N_p - 2
    print *, 'N_p = ', N_p
    close(500)

    !allocate the capacity
    allocate(alpha0(1:n_p))
    allocate(gamma0(1:n_p))
    allocate(energy0(1:n_p))
    allocate(alpha_eq(1:n_p))
    allocate(z_p(1:n_p))
    allocate(u_p(0:2, 1:n_p))
    allocate(u_p_eq(0:2, 1:n_p))
    allocate(v_eq(0:2, 1:n_p))
    allocate(equator_time(1:n_p))
    allocate(equator_flag(1:n_p))
    allocate(wave_flag(1:n_p))
    allocate(edge_flag(1:n_p))
    allocate(sign_theta0(1:n_p))   
    allocate(sign_theta1(1:n_p))
    allocate(cross_theta_0(1:n_p))
    allocate(Cw_flag(1:n_p))
    allocate(S_flag(1:n_p))
    
    !get the values
    write(file_data,'(A17, I3.3, A4)') 'initial_condition', myrank, '.dat'
    open (500, file = file_data)
  
    do i = 1, 1
       read(500,*)
    end do !i
  
    do i = 1, n_p
       read(500,*,iostat=ios) energy0(i), alpha0(i), z_p(i), alpha_eq(i) 
       
       gamma0(i) = DBLE(energy0(i)) / m_e + 1d0     
       v         = DSQRT(1d0 - 1d0/gamma0(i)**2)
       v_para    = v * DCOS(alpha0(i)*deg2rad)
       v_perp    = v * DSIN(alpha0(i)*deg2rad)
       u_p(0, i) = gamma0(i) * v_para
       u_p(1, i) = gamma0(i) * v_perp
       u_p(2, i) = 2d0 * pi * grnd() 
    end do !i








    
end program main