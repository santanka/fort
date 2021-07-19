program wave_front_setting
  !$ use omp_lib

  use constant_parameter
  use lshell_setting
  use variables
  use mt19937

  implicit none

  call initial_display
  call write_time(string)

  !---------------------------------
  ! initial setting of field (grid)
  !---------------------------------
 
  write(command, '(A13)') 'mkdir results'
  call system(command)  
  
  write(file_losscone, '(A20)') 'results/losscone.dat'
  open(10, file = file_losscone)
       
  do i_z = 0, n_z
    z(i_z) = DBLE(i_z) * d_z
    call z_to_r_lambda(z(i_z), r0, lambda0)
    call z_to_B(z(i_z), B(i_z))
    
    loss_cone(i_z) = DASIN(DSQRT(1d0 / B(i_z)))
    write(10, '(14E15.7)') z(i_z), B(i_z), lambda0 * rad2deg, loss_cone(i_z) * rad2deg
  end do
  
  close(10)



  write(file_wave, '(A20)') 'results/wave_v_g.dat'
  open(10, file = file_wave)
       
  do i_z = 0, n_z
    z(i_z) = DBLE(i_z) * d_z
    call z_to_B(z(i_z), B(i_z))
    call z_to_theta(z(i_z), theta(i_z))
!    call group_velocity(freq_start, B(i_z), w_p, V_g(i_z))   
    call parallel_group_velocity(freq_start, B(i_z), w_p, theta(i_z), V_g(i_z))   
    call wave_number(w_p, freq_start, B(i_z), theta(i_z), k_init(i_z))
    
    k_init(i_z) = k_init(i_z) * DCOS(theta(i_z))

    if (i_z == 0) then
      phase(i_z) = 0d0
    else
      phase(i_z) = phase(i_z - 1) - (k_init(i_z - 1) + k_init(i_z)) / 2d0  * d_z 
    end if
    
    write(10, '(14E15.7)') z(i_z), B(i_z), V_g(i_z), k_init(i_z), phase(i_z), theta(i_z)
  end do
  
  close(10)

!  call group_velocity(freq_start, B(0), w_p, V_g_0)
!  call group_velocity(freq_start, B(0), w_p, V_g_front)
!  call group_velocity(freq_end,   B(0), w_p, V_g_edge)
  call parallel_group_velocity(freq_start, B(0), w_p, theta(0), V_g_0)
!  call parallel_group_velocity(freq_start, B(0), w_p, theta(0), V_g_front)
!  call parallel_group_velocity(freq_end,   B(0), w_p, theta(0), V_g_edge)
  
  do i_w = 10, 50, 5
    w = DBLE(i_w) * 0.01d0
    write(file_wave, '(A16, I3.3, A4)') 'results/wave_res', i_w, '.dat'
    open(11, file = file_wave)

    do i_z = 0, n_z
      call resonance_velocity(w_p, w, B(i_z), gamma0, theta(i_z), V_R)
      call vpara_and_vabs_to_alpha(V_R, v0, alpha_R)
      alpha_R = alpha_R * rad2deg

      write(11, '(14E15.7)') z(i_z), alpha_R, V_R
    end do
    close(11)
  end do


  

  !------------------------------
  ! initial setting of particles
  !------------------------------
  count_p(:) = 0
  i_p = 0
!    call system_clock(count=clock) 
    clock = 1987912
    call sgrnd(clock)
    write(*,*) clock
      
      
      

    do i_p = 1, N_p
    alpha0 = alpha_min0 + (alpha_max0 - alpha_min0) * grnd()
         
  !------------------------------------------------
  ! setting of initial position z_p(i_p)
  !------------------------------------------------
!  INTEGER          :: n_z_mirror 
!  DOUBLE PRECISION :: alpha_eq, B_mirror, z_mirror, norm, cumulative(1:n_z)    
      if (one_point_start .eqv. .false.) then
	    B_mirror = 1.0 / DSIN(alpha0)**2
	    do i_z = 1, n_z
	      z_mirror = z(i_z - 1)
	      n_z_mirror = i_z - 1
	      if (B(i_z) > B_mirror) exit
	    end do
	  	
	    norm = 0.d0
        cumulative(:) = 0.d0
   
	    do i_z = 1, n_z_mirror
          norm = norm + d_z * DSQRT(B(i_z))
          cumulative(i_z) = norm
        end do
        
        cumulative(:) = cumulative(:) / norm
	  
	    rnd = grnd()
	    do i_z = 1, n_z_mirror
	      if(rnd < cumulative(i_z)) then
            z_p(i_p) = (rnd - cumulative(i_z - 1))/(cumulative(i_z) - cumulative(i_z - 1)) * d_z + DBLE(i_z - 1) * d_z
            exit
	      end if 
	    end do
      else
        z_p(:) = start_point
      end if 

  !----------------------------------------------------------------------------
  ! setting of initial velocity from pitch angle and energy  u_p(i_p)
  !----------------------------------------------------------------------------

      if (one_point_start .eqv. .false.) then     
        call z_to_B(z_p(i_p), B_p)
        alpha_p = DASIN(DSQRT(B_p) * DSIN(alpha0))
        rnd = grnd()
        if (rnd > 0.5) then
          u_p(0, i_p) = gamma0 * v0 * DCOS(alpha_p)
        else
          u_p(0, i_p) = - gamma0 * v0 * DCOS(alpha_p)    
        end if
        u_p(1, i_p) = gamma0 * v0 * DSIN(alpha_p)
        u_p(2, i_p) = grnd() * 2d0 * pi 
      
      else
        alpha_p = alpha_point0
        u_p(0, i_p) = - gamma0 * v0 * DCOS(alpha_p)
		u_p(1, i_p) = gamma0 * v0 * DSIN(alpha_p)
        u_p(2, i_p) = 2d0 * pi / N_p * i_p 
      end if
      
    end do
    
  write(*, *) "N_p = ", N_p
    
  
!  do i_z = 0, n_z, 5
!
!    if (DSQRT(B(i_z)) * DSIN(alpha_max0) < 1.0) then
!      alpha_max = DASIN(DSQRT(B(i_z)) * DSIN(alpha_max0))
!      alpha_min = DASIN(DSQRT(B(i_z)) * DSIN(alpha_min0))   
!    else if (DSQRT(B(i_z)) * DSIN(alpha_min0) < 1.0 .and. DSQRT(B(i_z)) * DSIN(alpha_max0) >= 1.0) then
!      alpha_max = pi / 2d0
!      alpha_min = DASIN(DSQRT(B(i_z)) * DSIN(alpha_min0))      
!    else if (DSQRT(B(i_z)) * DSIN(alpha_min0) >= 1.0) then
!      alpha_min = pi / 2d0
!      alpha_max = pi / 2d0
!    end if 
!    
!    do i_rdm = 1, N_rdm 
!      alpha0 = grnd() * pi
!      rnd    = grnd() 
!    
!      if ( ( (alpha0 <= alpha_max .and. alpha0 > alpha_min) .or.          & 
!      &      (alpha0 <= (pi - alpha_min) .and. alpha0 > (pi - alpha_max)) )  .and. &
!      &       rnd <= loss_cone(i_z) * 2.d0 / pi) then
!!      &       rnd <= abs(cos(alpha0)) * loss_cone(i_z) * 2.d0 / pi) then
!!        write(10, '(14E15.7)') DBLE(i_rdm), alpha0 * rad2deg, rnd
!
!        i_p = i_p + 1
!        
!        z_p(i_p)    = z(i_z)
!        u_p(0, i_p) = gamma0 * v0 * DCOS(alpha0)
!        u_p(1, i_p) = gamma0 * v0 * DSIN(alpha0)
!        u_p(2, i_p) = grnd() * 2d0 * pi
!        
!        count_p(i_z) = count_p(i_z) + 1
!!        write(10, '(14E15.7)') DBLE(i_p), z_p(i_p), alpha0*rad2deg, u_p(0, i_p), u_p(1, i_p), u_p(2, i_p), rnd
!      end if
!      
!    end do
!!    close(10)
!
!
!    
!    if (mod(i_z, 100) == 0) then
!      write(*, '(A4, I4, A5, I6, A5, F7.2, A2, F7.2, A2, F7.2, A2, F7.2, A2, F7.2)')&
!      & 'z = ', i_z, ' N = ', count_p(i_z), ' % = ', loss_cone(i_z) * 2.d0 / pi * 100 ! &
!     ! &, '  ', alpha_max * rad2deg, '  ', alpha_min *rad2deg, '  ', (pi - alpha_min) * rad2deg, '  ', (pi - alpha_max) *rad2deg
!     end if
!  end do


!  N_p_real = sum(count_p)
  
!  write(*, *) 'total N = ', N_p_real
!  if (N_p_real > N_p) then
!    write(*, *)  'Error: N_p_real > N_p'
!  end if
  
  open(11, file = 'results/position_t00000.dat')
  do i_p = 1, N_p!N_p_real 
    write(11, '(6E15.7)') z_p(i_p), DATAN2(u_p(1, i_p), u_p(0, i_p)) * rad2deg
  end do
  close(11)
  
  call write_time(string)




  !------------------------------
  ! initial distribution
  !------------------------------    
  do i_p = 1, N_p  
    call z_to_B(z_p(i_p), B_output)
    alpha_output = DASIN(DSIN(DATAN2(u_p(1, i_p), u_p(0, i_p)))/ DSQRT(B_output)) * rad2deg
	
	i_alpha_output = INT(alpha_output / d_alpha_output) + 1
!	write(*,*)u_p(1, i_p), u_p(0, i_p)
	count_output(i_alpha_output) = count_output(i_alpha_output) + 1
  end do
  
  open(99, file = 'results/distribution_str.dat')
    do i_alpha_output = 1, n_alpha_output
      alpha_output = DBLE(i_alpha_output) * d_alpha_output - d_alpha_output/2d0
      write(99, '(25E15.7)') alpha_output, DBLE(count_output(i_alpha_output))
    end do
  close(99)
  
  write(*,*)"start"

  !------------------------------
  ! open files
  !------------------------------


  do i_thr = 0, N_thr
    n_file = 20 + i_thr
    write(file_equator, '(A24, I2.2, A4)') 'results/count_at_equator', i_thr, '.dat'
    open(unit = n_file, file = file_equator)
  end do


if (one_point_start .eqv. .true.) then
  write(command, '(A26)') 'mkdir results/trajectories'
  call system(command)  

  do i_p = 1, N_p
    write(file_particle, '(A31, I3.3, A4)') 'results/trajectories/trajectory', i_p, '.dat'
    open(unit = 100 + i_p, file = file_particle)
  end do

  
  write(command, '(A31)') 'mkdir results/wave_at_onr_point'
  call system(command)  

  do i_z = 0, N_z, N_z/4
    write(file_wave, '(A37, I4.4, A4)') 'results/wave_at_one_point/wavefield_z', INT(i_z * d_z), '.dat'
    open(unit = 200 + i_z, file = file_wave)
  end do
end if
  
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

  if ((wave_packet .eqv. .true.)) then
    do i_z = -n_z, n_z
      if (i_z <= 0 .and. V_g_0 * t_element >= DABS(DBLE(i_z) * d_z)) then
        ampl_1(i_z) = B_w
      else if (i_z <= 0 .and. V_g_0 * t_element < DABS(DBLE(i_z) * d_z) &
      &.and. V_g_0 * t_element > DABS(DBLE(i_z) * d_z + 10d0 * d_z)) then
        ampl_1(i_z) = B_w * DEXP(-2d0 * pi / (10d0 * d_z)* (-z(-i_z) - V_g_0 * t_element)**2)
      else 
        ampl_1(i_z) = 0d0!grnd() * 1d-14
      end if
    end do
  else
    ampl_1(:) = B_w
  end if
  
  write(file_wave, '(A27)')  'results/wave_amp_t00000.dat'
  open(12, file = file_wave)
  do i_z = 0, n_z
    write(12, '(6E15.7)') z(i_z), ampl_1(i_z), freq_1(i_z)
  end do
  close(12) 
	    
	     
	     
	     
	     
	     
  !------------------------------
  ! simulation start
  !------------------------------
  count = 0
  z_front = 0d0
  z_edge  = 0d0


  do i_t = 1, n_t
    time = DBLE(i_t) * d_t


      !------------------------------
      ! initialize freq. and V_g
      !------------------------------
      do i_z = -n_z, n_z
        freq_0(i_z) = freq_1(i_z)
        ampl_0(i_z) = ampl_1(i_z)
      end do
      
      do i_z = -n_z, n_z
        if (i_z > 0) then
          call parallel_group_velocity(freq_0(i_z), B(i_z), w_p, theta(i_z), V_g(i_z))
        else
          V_g(i_z) = V_g_0
        end if
      end do

      !------------------------------
      ! update freq. and ampl.
      !------------------------------    
      call eno(freq_0, V_g * d_t / d_z, n_z, freq_1)
      call eno(ampl_0, V_g * d_t / d_z, n_z, ampl_1)
            
      !------------------------------
      ! update wave front 
      !------------------------------    
!      if (z_front <= DBLE(n_z) * d_z) then
!        call z_to_B(z_front, B_front)
!        call group_velocity(freq_start, B_front, w_p, V_g_front)        
!        z_front = z_front + V_g_front * d_t
!      end if
!      
!      if (time >= t_element .and. z_edge <= DBLE(n_z) * d_z) then
!        call z_to_B(z_edge, B_edge)
!        call group_velocity(freq_end, B_edge, w_p, V_g_edge)        
!        z_edge = z_edge + V_g_edge * d_t
!      end if
!      
!      write(30, '(4E15.7)') time, z_front, z_edge

      !------------------------------
      ! update phase
      !------------------------------
      
      do i_z = 0, n_z
        phase(i_z) = phase(i_z) + (freq_0(i_z) + freq_1(i_z)) / 2d0 * d_t
      end do


if (one_point_start .eqv. .true.) then
  do i_z = 0, N_z, N_z/4
    call output_wave(w_p, freq_1(i_z), B(i_z), ampl_1(i_z), theta(i_z), phase(i_z), Bw_x, Bw_y, Bw_z, Ew_x, Ew_y, Ew_z)
    write(unit = 200 + i_z, fmt = '(14E15.7)') time, Bw_x, Bw_y, Bw_z, Ew_x, Ew_y, Ew_z, phase(i_z) 
  end do
end if
    !$omp parallel num_threads(n_thread) &
    !$omp & private(i_p, i_t, z_p_sim, u_p_sim, & !time, z, k_init, phase, freq_1, ampl_1, 
    !$omp & zeta, equator_flag, alpha, gamma, energy, string, n_file, count, total, progress)
    !$omp do
    do i_p = 1, N_p!N_p_real 
  
      z_p_sim    = z_p(i_p)
      u_p_sim(:) = u_p(:, i_p)

      
      n_file = 20 + omp_get_thread_num()
      count = count + 1
    
!     do i_t = 1, n_t
!       time = DBLE(i_t) * d_t
 

      call runge_kutta(z, phase, freq_1(0:n_z), ampl_1(0:n_z), z_p_sim, u_p_sim, zeta, equator_flag)


      !------------------------------
      ! output particle at equator
      !------------------------------

      alpha = DATAN2(u_p_sim(1), u_p_sim(0)) * rad2deg
      call p_to_energy(u_p_sim(0:2), gamma, energy)

      if (equator_flag .eqv. .true.) then
        write(unit = n_file, fmt = '(7E15.7)') time, z_p_sim, alpha, u_p_sim(0), u_p_sim(1), u_p_sim(2), energy
      end if


    !------------------------------
    ! timestamp 
    !------------------------------

!    if ( mod(count, 100) == 0 .and. omp_get_thread_num() == 4) then
!      total = N_p_real / omp_get_num_threads()
!      progress = i_p - omp_get_thread_num() * total
!      
!      write(string, '(F6.2, A2, A6, I2.2)') DBLE(progress) / DBLE(total)*100.0, ' %', &
!      & 'CPU =', omp_get_thread_num()
!      call write_time(string)
!      
!    end if


    z_p(i_p)    = z_p_sim
    u_p(:, i_p) = u_p_sim(:)
    !------------------------------------------
    ! particle loop end
    !------------------------------------------
    end do
   !$omp end do
   !$omp end parallel
   
      !------------------------------------------
      ! update wave outside the simulation box  
      !------------------------------------------

      if (frequency_sweep .eqv. .true.) then
        do i_z = -n_z, 0
          freq_1(i_z) = freq_0(i_z) + dfdt * d_t
        end do
      end if

!      if (wave_packet .eqv. .true.) then
!        if (time > t_element) then
!          ampl_1(i_z) = grnd() * 1d-7
!        end if
!      end if

    !------------------------------
    ! output position and wave 
    !------------------------------

	  if (mod(i_t, 300) == 0) then
	  
	    write(file_position, '(A18, I5.5, A4)')  'results/position_t', i_t, '.dat'
	    open(11, file = file_position)
	    do i_p = 1, N_p!N_p_real 
	      write(11, '(6E15.7)') z_p(i_p), DATAN2(u_p(1, i_p), u_p(0, i_p)) * rad2deg
	    end do
	    close(11)

	    write(file_wave, '(A18, I5.5, A4)')  'results/wave_amp_t', i_t, '.dat'
	    open(12, file = file_wave)
	    do i_z = 0, n_z
	      write(12, '(6E15.7)') z(i_z), ampl_1(i_z), freq_1(i_z)
	    end do
	    close(12) 
	   	    
        write(string, '(F6.2, A2, A6, I2.2)') DBLE(i_t) / DBLE(n_t)*100.0, ' %'
        call write_time(string)
	  end if
	  
if (one_point_start .eqv. .true.) then
	  do i_p = 1, N_p
!	    call force(z, phase, freq_1(0:n_z), ampl_1(0:n_z), z_p(i_p), u_p(:, i_p), force_p, zeta_p)
  	    call p_to_energy(u_p(0:2, i_p), gamma, energy)
	    write(unit = 100 + i_p, fmt = '(14E15.7)') time, z_p(i_p), DATAN2(u_p(1, i_p), u_p(0, i_p)) * rad2deg, &
	    & u_p(0, i_p), u_p(1, i_p), u_p(1, i_p), gamma, energy, zeta
	  end do
end if	  
  end do


  open(11, file = 'results/position_after.dat')
  do i_p = 1, N_p!N_p_real 
    write(11, '(6E15.7)') z_p(i_p), DATAN2(u_p(1, i_p), u_p(0, i_p)) * rad2deg
  end do
  close(11)
  
  !------------------------------
  ! final distribution
  !------------------------------    
  
  count_output(:) = 0d0
  do i_p = 1, N_p!N_p_real  
    call z_to_B(z_p(i_p), B_output)
    alpha_output = DASIN(DSIN(DATAN2(DABS(u_p(1, i_p)), u_p(0, i_p)))/ DSQRT(B_output)) * rad2deg
	
	i_alpha_output = INT(alpha_output / d_alpha_output) + 1
	count_output(i_alpha_output) = count_output(i_alpha_output) + 1
  end do
  
  open(99, file = 'results/distribution_end.dat')
    do i_alpha_output = 1, n_alpha_output
      alpha_output = DBLE(i_alpha_output) * d_alpha_output - d_alpha_output/2d0
      write(99, '(25E15.7)') alpha_output, DBLE(count_output(i_alpha_output))
    end do
  close(99)
  
  write(*,*) "end"
end program





