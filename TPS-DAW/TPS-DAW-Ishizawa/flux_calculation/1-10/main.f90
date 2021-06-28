program main
  !$ use omp_lib
  
  use mt19937
  use constant_parameter
  use lshell_setting
  use variables

  implicit none
!--MPI-------------------------------------------------------------------
   include 'mpif.h'
   INTEGER(KIND=4)  ierr,nprocs,myrank
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!--MPI-------------------------------------------------------------------

  call initial_display
  call write_time(string)

  write(command, '(A13, I3.3)') 'mkdir results', myrank
  call system(command)  
    
  do i_z = -n_z, n_z   
     z(i_z) = DBLE(i_z) * d_z
     call z_to_B(z(i_z), B(i_z))
     call group_velocity(freq_start, B(i_z), w_p, V_g(i_z))   
     call wave_number(w_p, freq_start, B(i_z), k_init(i_z))
     if (i_z == -n_z) then     
        phase(i_z) = 0d0
     else
        phase(i_z) = phase(i_z - 1) - (k_init(i_z - 1) + k_init(i_z)) / 2d0 * d_z
     end if
  end do
  

  call group_velocity(freq_start, B(0), w_p, V_g_0)
  call group_velocity(freq_start, B(0), w_p, V_g_front)
  call group_velocity(freq_end,   B(0), w_p, V_g_edge)

  
  !------------------------------
  ! initial setting of particles
  !------------------------------

  call system_clock(count=clock) 
  clock = 4267529
  call sgrnd(clock)
  write(*,*) clock

    

  write(file_data,'(A17, I3.3, A4)') 'initial_condition', myrank, '.dat'
  open (500, file = file_data)
  ios = 1
  n_p = 0

  do 
     read(500,*,end = 99)
     n_p = n_p + 1
  end do
  
99 n_p = n_p-2
  write(*, *) "N_p = ", n_p
  close(500)

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
  
  write(file_data,'(A17, I3.3, A4)') 'initial_condition', myrank, '.dat'
  open (500, file = file_data)

  do i = 1, 1
     read(500,*)
  end do

  do i = 1, n_p
     read(500,*,iostat=ios) energy0(i), alpha0(i), z_p(i), alpha_eq(i) 
     
     gamma0(i) = DBLE(energy0(i)) / m_e + 1d0     
     v         = DSQRT(1d0 - 1d0/gamma0(i)**2)
     v_para    = v * DCOS(alpha0(i)*deg2rad)
     v_perp    = v * DSIN(alpha0(i)*deg2rad)
     u_p(0, i) = gamma0(i) * v_para
     u_p(1, i) = gamma0(i) * v_perp
     u_p(2, i) = 2d0 * pi * grnd() 

     if (ios < 0) exit
  end do  

  
  sign_theta0     = 0
  sign_theta1     = 0
  cross_theta_0   = 0
  Cw_flag         = 0
  S_flag          = 0
  equator_flag    = 0
  equator_time    = 0
  wave_flag       = 0
  edge_flag       = 0
  
  
  !------------------------------
  ! initial setting of wave
  !------------------------------
  if (frequency_sweep .eqv. .true.) then
     do i_z = -n_z, n_z   
        if (i_z <= 0) then
           freq_1(i_z) = freq_start + dfdt / V_g_0 * DBLE(abs(i_z)) * d_z
        else
           freq_1(i_z) = freq_start
        end if
     end do
  else
     freq_1(:) = freq_start
  end if
  
  
  if (wave_existance .eqv. .true.) then
     call wave_packet(ampl_1(-n_z:n_z))
  else
     ampl_1(:) = 0d0
  end if
  
   
  z_front         = 0d0
  z_edge          = 0d0


   !-------------------------------
   ! file open
   !-----------------------------
  
  do i_thr = 0, N_thr
     n_file = 20 + i_thr
     write(file_equator, '(A7, I3.3, A17, I2.2, A4)') 'results', myrank, '/count_at_equator', i_thr, '.dat'
     open(unit = n_file, file = file_equator)
  end do

  
  !------------------------------
  ! simulation start
  !------------------------------
   
  do i_t = 1, n_t
     time = DBLE(i_t) * d_t

      !------------------------------
      ! initialize freq. and V_g     !update_wave_at_all_space
      !------------------------------
     
      do i_z = -n_z, n_z
         freq_0(i_z) = freq_1(i_z)
         ampl_0(i_z) = ampl_1(i_z)

         if (i_z >= 0) then
            call group_velocity(freq_0(i_z), B(i_z), w_p, V_g(i_z))
         else
            V_g(i_z) = V_g_0
         end if
      end do
 
      !------------------------------
      ! update freq. and ampl.
      !------------------------------    

      call eno(freq_0, V_g * d_t / d_z, n_z, freq_1)
      call eno(ampl_0, V_g * d_t / d_z, n_z, ampl_1)
      
      do i_z = -n_z, n_z
         
         if (damping_region <= DBLE(i_z) * d_z) then
            freq_1(i_z) = freq_start
            ampl_1(i_z) = 0d0
         end if
         
         if (ampl_1(i_z) < 1d-20 ) then
            ampl_1(i_z) = 0d0            
         end if
      end do
      
      !------------------------------
      ! update phase
      !------------------------------
      
      do i_z = -n_z, n_z
         phase(i_z) = phase(i_z) + (freq_0(i_z) + freq_1(i_z)) / 2d0 * d_t
      end do
      
      !$omp parallel num_threads(n_thread) &
      !$omp & private(i_p, i_t, z_p_sim, u_p_sim, equator_flag_sim, wave_flag_sim, edge_flag_sim, & 
      !$omp & zeta, n_file, alpha_p, gamma_p, energy_p, B_p, alpha, v, v_para, v_perp, &
      !$omp & freq_p, ampl_p, k_p, zeta_p, V_R_p, Cw_p, theta_p, w_tr_p, V_g_p, dB_dz_p, dk_dB_p, S_p, &
      !$omp & sign_theta0_sim, sign_theta1_sim, cross_theta_0_sim, Cw_flag_sim, S_flag_sim)
      !$omp do
      do i_p = 1, N_p
         z_p_sim    = z_p(i_p)
         u_p_sim(:) = u_p(:, i_p)     
         equator_flag_sim  = equator_flag(i_p)
         wave_flag_sim     = wave_flag(i_p)
         edge_flag_sim     = edge_flag(i_p)
         sign_theta0_sim   = sign_theta0(i_p)
         sign_theta1_sim   = sign_theta1(i_p)
         cross_theta_0_sim = cross_theta_0(i_p)
         Cw_flag_sim       = Cw_flag(i_p)
         S_flag_sim        = S_flag(i_p)
         
         n_file = 20 + omp_get_thread_num()
         
         call runge_kutta (z, phase, freq_1(-n_z:n_z), ampl_1(-n_z:n_z), &
              & z_p_sim, u_p_sim, zeta, equator_flag_sim, wave_flag_sim, edge_flag_sim)
         
         if (equator_flag_sim == 1) then
            call p_to_energy(u_p_sim(0:2), gamma_p, energy_p)
            call z_to_B(z_p_sim, B_p)
            alpha_p = DASIN(DSIN(DATAN2(DABS(u_p_sim(1)), u_p_sim(0)))/ DSQRT(B_p)) * rad2deg   
            write(unit = n_file, fmt = '(5E15.7, 4I3)') time, alpha_p, energy_p, alpha_eq(i_p), energy0(i_p), &
                 & wave_flag(i_p), cross_theta_0(i_p), Cw_flag(i_p), S_flag(i_p) 
         end if
         
    !------------------------------------------
    ! categorization  
    !------------------------------------------
         if (categorization .eqv. .true.) then
            if (mod(i_t, 1) == 0) then 
               if ( u_p_sim(0) < 0d0 .and. z_p_sim <= damping_region) then                
                  call z_to_B(z_p_sim, B_p)
                  call p_to_energy(u_p_sim(0:2), gamma_p, energy_p)
                  call wave_at_z(z_p_sim, z, u_p_sim(0:2), phase, &
                       & freq_1(-n_z:n_z), ampl_1(-n_z:n_z), freq_p, ampl_p, k_p, zeta_p)   
                  call resonance_velocity(w_p, freq_p, B_p, gamma_p, V_R_p)
                  call to_theta(k_p, u_p_sim(0:2), gamma_p, V_R_p, theta_p)
                  call to_Cw(ampl_p, freq_p, k_p, gamma_p, u_p_sim(0:2), Cw_p)
                  call to_w_tr(freq_p, theta_p, k_p, u_p_sim(0:2), ampl_p, gamma_p, w_tr_p)
                  call group_velocity(freq_p, B_p, w_p, V_g_p)
                  call z_to_dB_dz(z_p_sim, dB_dz_p)
                  call to_dk_dB(freq_p, w_p, k_p, B_p, dk_dB_p)
                  call to_S(u_p_sim(0:2), gamma_p, V_g_p, k_p, B_p, dk_dB_p, dB_dz_p, w_tr_p, S_p)
                  
                  sign_theta1_sim = theta_p
                  
                  if( cross_theta_0_sim == 0 .and. ampl_p > 1d-20 ) then
                     
                     if (DABS(Cw_p / theta_p) > 1d0) then
                        Cw_flag_sim = 1
                     end if
                     
                     if (DABS(Cw_p / theta_p) > 1d0 .and. DABS(S_p) > 1d0) then
                        S_flag_sim = 1
                     end if
                     
                  end if
                  
                  
                  if ( sign_theta0_sim*sign_theta1_sim < 0d0 .and. ampl_p > 1d-20 ) then
                     cross_theta_0_sim = cross_theta_0_sim + 1                
                  end if
                  
                  sign_theta0_sim = sign_theta1_sim
                  
               end if
               
            end if
         end if
         
         
         if (equator_flag_sim == 1) then            
            z_p_sim    = grnd() * DABS(z_p_sim - 0d0) + 0d0
            call z_to_B(z_p_sim, B_p)
            alpha     = DASIN(DSQRT(B_p*DSIN(alpha_eq(i_p)*deg2rad)**2))
            v         = DSQRT(1d0 - 1d0/gamma0(i_p)**2)
            v_para    = v * DCOS(alpha)
            v_perp    = v * DSIN(alpha)
            u_p_sim(0) = gamma0(i_p) * v_para
            u_p_sim(1) = gamma0(i_p) * v_perp
            u_p_sim(2) = 2d0 * pi * grnd() 
            equator_flag_sim  = 0
            wave_flag_sim     = 0
            Cw_flag_sim       = 0
            S_flag_sim        = 0
            cross_theta_0_sim = 0
         end if
         
         if (edge_flag_sim == 1) then
            z_p_sim    = grnd() * DABS(z_p_sim - L_z) + z_p_sim
            call z_to_B(z_p_sim, B_p)
            alpha     = DASIN(DSQRT(B_p*DSIN(alpha_eq(i_p)*deg2rad)**2))
            v         = DSQRT(1d0 - 1d0/gamma0(i_p)**2)
            v_para    = - v * DCOS(alpha)
            v_perp    = v * DSIN(alpha)
            u_p_sim(0) = gamma0(i_p) * v_para
            u_p_sim(1) = gamma0(i_p) * v_perp
            u_p_sim(2) = 2d0 * pi * grnd() 
            edge_flag_sim = 0
         end if
         
         z_p(i_p)    = z_p_sim
         u_p(:, i_p) = u_p_sim(:)
         equator_flag(i_p)  = equator_flag_sim
         wave_flag(i_p)     = wave_flag_sim
         edge_flag(i_p)     = edge_flag_sim
         sign_theta0(i_p)   = sign_theta0_sim
         sign_theta1(i_p)   = sign_theta1_sim
         cross_theta_0(i_p) = cross_theta_0_sim
         Cw_flag(i_p)       = Cw_flag_sim
         S_flag(i_p)        = S_flag_sim
         
      end do !particle loop end
      !$omp end do nowait
      !$omp end parallel
   
      !------------------------------------------
      ! update wave outside the simulation box  
      !------------------------------------------

      if (frequency_sweep .eqv. .true.) then
         do i_z = -n_z, 0
            freq_1(i_z) = freq_0(i_z) + dfdt * d_t
         end do
      end if
      

      !------------------------------------------
      ! output
      !------------------------------------------
      if (mod(i_t, 1000) == 0) then
         write(string, '(F6.2, A2, A6, I2.2)') DBLE(i_t) / DBLE(n_t)*100.0, ' %'
         call write_time(string)
      end if
             
      
   end do !time_loop_end
   close(200)
   write(*,*) "end"

!-MPI-------------------------------------------------------------------
   call MPI_FINALIZE(ierr)      
!-MPI-------------------------------------------------------------------

 end program main
 
