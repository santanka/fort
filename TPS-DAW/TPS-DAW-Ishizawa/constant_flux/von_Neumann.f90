program distribution_calculater

  use mt19937
  use constant_parameter
  use lshell_setting

  implicit none
  

  !--------------------------------------
  !  constants_in_the_simulation
  !--------------------------------------

  DOUBLE PRECISION, PARAMETER :: ave_p_para = 0.2d0
  DOUBLE PRECISION, PARAMETER :: ave_p_perp = 0.2d0 
  DOUBLE PRECISION, PARAMETER :: p_para_max = -0d0
  DOUBLE PRECISION, PARAMETER :: p_para_min = -1d0
  DOUBLE PRECISION, PARAMETER :: p_perp_max = 1d0
  DOUBLE PRECISION, PARAMETER :: p_perp_min = 0d0
  DOUBLE PRECISION, PARAMETER :: fmax = 10d0
  DOUBLE PRECISION, PARAMETER :: fmin = 0d0
  DOUBLE PRECISION, PARAMETER :: beta = 0.1
  INTEGER, PARAMETER          :: N_p  = 20000000
  INTEGER, PARAMETER          :: n_phase_space = 100
  INTEGER, PARAMETER          :: n_z = 16000
  DOUBLE PRECISION, PARAMETER :: d_z = 0.5d0  
  DOUBLE PRECISION, PARAMETER :: d_t = 1.0d0
  DOUBLE PRECISION, PARAMETER :: L_z = DBLE(n_z) * d_z 
  DOUBLE PRECISION, PARAMETER :: d_a = 1d0
  DOUBLE PRECISION, PARAMETER :: d_E = 1d0
  
  !-------------------------------------
  ! variables
  !-------------------------------------
  INTEGER :: clock
  INTEGER :: i_z, i_alpha, i_energy, i_lambda, i_E, i_a, i_v_para, i_v_perp
  INTEGER :: i_tau, i_p
  DOUBLE PRECISION :: a, E
  DOUBLE PRECISION :: z(0 : n_z)
  DOUBLE PRECISION :: B(0 : n_z)
  DOUBLE PRECISION :: rnd, count
  CHARACTER(64)    :: file_distribution
  CHARACTER(64)    :: file_alpha_energy
  CHARACTER(64)    :: file_phase_space
  CHARACTER(64)    :: file_initial_condition
  CHARACTER(64)    :: file_bounce_period

  !------------------------------
  ! for distribution function
  !------------------------------

  DOUBLE PRECISION :: f, fp
  DOUBLE PRECISION :: v
  DOUBLE PRECISION :: p_para, p_perp
  DOUBLE PRECISION :: v_para, v_perp
  DOUBLE PRECISION :: v_0, v_1
  DOUBLE PRECISION :: alpha, alpha_eq, energy ,phi
  DOUBLE PRECISION :: gamma

  DOUBLE PRECISION :: theta
  DOUBLE PRECISION,PARAMETER :: kappa = 2d0
  
  
  !-------------------------------------
  ! for count
  !-------------------------------------
  INTEGER :: count_alpha_energy(0:89, 0:300)
  INTEGER :: weight_count_alpha_energy(0:89, 0:300)
  INTEGER :: count_phase_space(-n_phase_space:0, 0:n_phase_space)
  INTEGER :: total_alpha_energy
  INTEGER :: particle_num_for_plot_a_E
  INTEGER :: phase_space_for_plot
  INTEGER :: particle_num
  DOUBLE PRECISION :: norm_tau_b
  
  !------------------------------
  ! for initial  position
  !------------------------------ 
!  INTEGER,PARAMETER:: myrank = 20
  DOUBLE PRECISION :: lambda, lambda1, lambda2
  DOUBLE PRECISION,PARAMETER :: delta_t = 200
  DOUBLE PRECISION :: r, B1, B2
  DOUBLE PRECISION :: B0, alpha_eq0
  DOUBLE PRECISION :: z_0, z_small, z_large
  DOUBLE PRECISION :: z_p, z_p_1, z_p_2
  DOUBLE PRECISION :: p(0:2)
  DOUBLE PRECISION :: p_1(0:2), p_2(0:2)
  DOUBLE PRECISION :: alpha_1, alpha_2
  DOUBLE PRECISION :: tau_b
  DOUBLE PRECISION :: t_2
  DOUBLE PRECISION :: time
  DOUBLE PRECISION :: B_mirror, z_mirror, n_z_mirror
  DOUBLE PRECISION :: z_p_1_to_z_mirror, z_p_2_to_z_mirror

  !------------------------------
  ! for parallelization
  !------------------------------ 

  !INTEGER, PARAMETER :: n_parallelization = 100

    !--MPI-------------------------------------------------------------------
    include 'mpif.h'
    INTEGER(KIND=4)  ierr,nprocs,myrank
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    !--MPI-------------------------------------------------------------------


  !-------------------------------------------
  !   initial setting
  !-------------------------------------------
  
  do i_z = 0, n_z   
     z(i_z) = DBLE(i_z) * d_z
     call z_to_B(z(i_z), B(i_z))
  end do

  call system_clock(count=clock) 
  clock = 4267529
  call sgrnd(clock)
  write(*,*) clock

  
  !-------------------------------------------
  ! decide distribution function
  !------------------------------------------   
   
   write(file_alpha_energy, '(A33)') 'distribution_alpha_energy_str.dat'
   write(file_phase_space, '(A32)') 'distribution_phase_space_str.dat'
   write(file_bounce_period, '(A17)') 'bounce_period.dat'
   open (11, file = file_alpha_energy)
   open (12, file = file_phase_space)
   open (13, file = file_bounce_period)
   
   count = 0
   count_alpha_energy = 0
   count_phase_space = 0
   
   do while(count <= N_p)
  
      p_para = grnd() * (p_para_max - p_para_min) + p_para_min
      p_perp = grnd() * (p_perp_max - p_perp_min) + p_perp_min
      gamma  = DSQRT(1 + p_para**2 + p_perp**2)
      v_para = p_para / gamma
      v_perp = p_perp / gamma
      energy   = m_e * (gamma - 1.0)
      alpha_eq = pi - DATAN2(p_perp, p_para)
      fp     = grnd() * (fmax - fmin) + fmin

      
!      f = 2d0*pi*v_perp* 1d0 /((2d0*pi)*ave_p_para*ave_p_perp**2) * DEXP(-p_para**2/(2d0*ave_p_para**2)) * &
!           & 1d0/(1d0-beta)*(DEXP(-p_perp**2/(2d0*ave_p_perp**2))-DEXP(-p_perp**2/(2d0*beta*ave_p_perp**2)))

      f = 2d0*pi*v_perp* 1d0 /((2d0*pi)**(3d0/2d0)*ave_p_para*ave_p_perp**2) &
           * DEXP(-p_para**2/(2d0*ave_p_para**2)) * DEXP(-p_perp**2/(2d0*ave_p_perp**2))

      ! f = 2d0*pi*v_perp* 1d0 /( pi*kappa*((ave_p_para/gamma)**2+(ave_p_perp/gamma)**2) )**(3d0/2d0) &
      !      * DGAMMA(kappa+1d0) / DGAMMA(kappa-1d0/2d0) &
      !      * (1d0 + (v_para**2 + v_perp**2)/(kappa*((ave_p_para/gamma)**2+(ave_p_perp/gamma)**2)) )**(-(kappa+1d0))

      
      if( fp <= f )then

            
         count = count + 1
 !        write(10,'(6E15.7)') v_para, v_perp, phi, energy, alpha_eq
            
         do i_a = 0, 89
            a = DBLE(i_a) * d_a
            do i_E = 0, 300
               E = DBLE(i_E) * d_E
               if (a <= alpha_eq*rad2deg .and. alpha_eq*rad2deg <= a+1) then
                  if (E <= energy .and. energy <= E+1) then
                     count_alpha_energy(i_a,i_E) = count_alpha_energy(i_a,i_E) + 1
                    end if
               end if
            end do
         end do
         
         do i_v_para = -n_phase_space, 0
            do i_v_perp = 0, n_phase_space
               v_0 = i_v_para/DBLE(n_phase_space)
               v_1 = i_v_perp/DBLE(n_phase_space)
               if (v_0 <= v_para .and. v_para <= v_0 + 1/DBLE(n_phase_space) ) then
                  if (v_1 <= v_perp .and. v_perp <= v_1 + 1/DBLE(n_phase_space) ) then
                     count_phase_space(i_v_para, i_v_perp) = count_phase_space(i_v_para, i_v_perp) + 1 
                  end if
               end if
            end do
         end do
         
      end if
   end do

   do i_a = 0, 89
      a = DBLE(i_a) * d_a
      do i_E = 0, 300
         E = DBLE(i_E) * d_E
         call energy_to_v(E+d_E/2d0, gamma, v)
         tau_b      = 2d0 * r_eq / DSQRT(1d0 - 1d0/gamma**2) * (1.30d0 - 0.56d0 * DSIN((a+d_a/2d0)*deg2rad))
         norm_tau_b = 2d0 * r_eq / DSQRT(1d0 - 1d0/gamma**2) * (1.30d0 - 0.56d0 * DSIN( 89.5d0*deg2rad))
         weight_count_alpha_energy(i_a,i_E) = count_alpha_energy(i_a,i_E)* tau_b /norm_tau_b !weighting_at_pitch_angle
         particle_num = DBLE(count_alpha_energy(i_a,i_E)) / (tau_b/delta_t)
         particle_num_for_plot_a_E = DBLE(count_alpha_energy(i_a, i_E))/(2d0*pi*DSQRT(1d0-1d0/gamma**2)*DSIN((a+d_a/2d0)*deg2rad))
         write(11,'(2E15.7, 3I10)') a, E, particle_num_for_plot_a_E, count_alpha_energy(i_a,i_E), particle_num
      end do
      write(11,*)''
   end do
   
   
   do i_v_para = -n_phase_space, -1
      do i_v_perp = 1, n_phase_space
         v_0 = i_v_para/DBLE(n_phase_space)
         v_1 = i_v_perp/DBLE(n_phase_space)
         phase_space_for_plot = count_phase_space(i_v_para, i_v_perp) /(2d0*pi*(v_1+0.5d0/DBLE(n_phase_space) ))
         write(12,'(2E15.7, 2I10)') v_0, v_1, phase_space_for_plot, count_phase_space(i_v_para, i_v_perp)
      end do
      write(12,*)''
   end do

   do i_a = 0, 89
      a = DBLE(i_a) * d_a
      do i_E = 0, 300
         E = DBLE(i_E) * d_E
         
         energy   =  DBLE(E+d_E/2d0) 
         gamma    = energy / m_e + 1d0
         tau_b    = 2d0 * r_eq / DSQRT(1d0 - 1d0/gamma**2) * (1.30d0 - 0.56d0 * DSIN((a+d_a/2d0)*deg2rad))
         
         write(13,'(3E15.7)') a, E, tau_b
      end do
      write(13,*)''
   end do
   
   
   
   !close(10)
   close(11)
   close(12)
   close(13)
   !----------------------------------------------------------
   !   decide_z_distribution
   !----------------------------------------------------------

   myrank = myrank + 10
   write(file_initial_condition, '(A17, I3.3, A4)') 'initial_condition', myrank, '.dat'
   open (14, file = file_initial_condition)
   
   particle_num = 0
   tau_b        = 0
   
   do i_a = 5, 89
      a = DBLE(i_a) * d_a

      i_E = myrank  !do i_E = 0, 100
         E = DBLE(i_E) * d_E
         
         if (weight_count_alpha_energy(i_a,i_E) >= 1) then 
            alpha_eq = (DBLE(a) + d_a/2d0)* deg2rad
            energy   =  DBLE(E) + d_E/2d0
            gamma    = energy / m_e + 1d0
            v        = DSQRT(1d0 - 1d0/gamma**2)
            v_para   = v * DCOS(alpha_eq)
            v_perp   = v * DSIN(alpha_eq)
            p(0)     = gamma * v_para
            p(1)     = gamma * v_perp
            p(2)     = 2 * pi * grnd()
            tau_b    = 2d0 * r_eq / DSQRT(1d0 - 1d0/gamma**2) * (1.30d0 - 0.56d0 * DSIN(alpha_eq))
            
            B_mirror = 1.0 / DSIN(alpha_eq)**2
            do i_z = 1, n_z
               z_mirror = z(i_z - 1)
               n_z_mirror = i_z - 1
               if (B(i_z) > B_mirror) exit
            end do
            
            particle_num   = DBLE(weight_count_alpha_energy(i_a,i_E)) / (tau_b/delta_t) 

            if (particle_num >= 1) then

               time = 0d0
               z_p = 0d0
               
               do i_tau = 1, 1000000
                  
                  t_2 = delta_t * i_tau
                  
                  if ( t_2 <= tau_b ) then
                     
                     z_p_1 = z_p
                     p_1   = p
                     
                     do while(time <= t_2)
                        call runge_kutta(z, z_p, p(0:2))
                        z_p_2 = z_p
                        p_2   = p
                        time = time + 1d0 * d_t
                     end do
                     
                     if ( z_p_1 < z_p_2 ) then 
                        z_small = z_p_1
                        z_large = z_p_2
                     else
                        z_small = z_p_2
                        z_large = z_p_1                        
                     end if
                     
                     call z_to_B(z_small, B1)
                     call z_to_B(z_large, B2)                                  
                        
                     z_p_1_to_z_mirror = z_mirror - z_p_1
                     z_p_2_to_z_mirror = z_mirror - z_p_2
                     
                     do i_p = 1, particle_num
                        
                        if (p_1(0) > 0 .and. p_2(0) > 0) then
                           alpha_1  = DASIN(DSQRT(B1 * DSIN(alpha_eq)**2) )
                           alpha_2  = DASIN(DSQRT(B2 * DSIN(alpha_eq)**2) )
                           alpha  = grnd() * (alpha_2 - alpha_1) + alpha_1                             
                           z_0    = grnd() * (z_p_2 - z_p_1) + z_p_1
                        end if
                        
                        if (p_1(0) > 0 .and. p_2(0) < 0) then
                           
                           if (grnd() < (z_p_1_to_z_mirror)/(z_p_1_to_z_mirror + z_p_2_to_z_mirror)) then
                              alpha_1  = DASIN(DSQRT(B1 * DSIN(alpha_eq)**2) )
                              alpha_2  = pi / 2d0
                              alpha  = grnd() * (alpha_2 - alpha_1) + alpha_1
                              z_0    = grnd() * (z_mirror - z_p_1) + z_p_1
                           else
                              alpha_1  = pi / 2d0
                              alpha_2  = pi - DASIN(DSQRT(B2 * DSIN(alpha_eq)**2) )
                              alpha  = grnd() * (alpha_2 - alpha_1) + alpha_1
                              z_0    = grnd() * (z_mirror - z_p_2) + z_p_2
                           end if
                           
                        end if
                        
                        if (p_1(0) < 0 .and. p_2(0) < 0) then                                 
                           alpha_1  = pi - DASIN(DSQRT(B1 * DSIN(alpha_eq)**2) )
                           alpha_2  = pi - DASIN(DSQRT(B2 * DSIN(alpha_eq)**2) )
                           alpha  = grnd() * (alpha_1 - alpha_2) +  alpha_2
                           z_0    = grnd() * (z_p_1 - z_p_2) + z_p_2                              
                        end if

                        call z_to_B(z_0, B0)
                        alpha_eq0 = DASIN(DSQRT(1d0 /B0 * DSIN(alpha)**2) )
                        if (z_0 <= 8000d0 .and. alpha_eq0*rad2deg > 5.6d0) then
                           energy   = grnd() * (DBLE(E) + 1d0 - DBLE(E)) +DBLE(E)
                           gamma    = energy / m_e + 1d0
                           v        = DSQRT(1d0 - 1d0/gamma**2)                         
                           v_para = v * DCOS(alpha)
                           v_perp = v * DSIN(alpha)
                           p_para = gamma * v_para
                           p_perp = gamma * v_perp

                           if (alpha == alpha) then
                              write(14,'(7E15.7)') energy, alpha*rad2deg, z_0, alpha_eq0*rad2deg
                           end if
                              
                        end if
                     end do
                  end if
               end do
            end if
         end if
            
         
         if (mod(i_a,5) == 0) then            
            write(*,*) 'alpha', a, 'energy', myrank
         end if

      end do
   
      write(*,*) "end"

      close(14)

!-MPI-------------------------------------------------------------------
   call MPI_FINALIZE(ierr)      
!-MPI-------------------------------------------------------------------
   
!---------------------------------------------------------------------------------
contains

subroutine runge_kutta(z_f, z_p, u_p)

  implicit none

  DOUBLE PRECISION, INTENT(IN)    :: z_f(0:n_z)
  DOUBLE PRECISION, INTENT(INOUT) :: z_p, u_p(0:2)
  DOUBLE PRECISION :: l1(0:2), l2(0:2), l3(0:2), l4(0:2), u_p_s(0:2)
  DOUBLE PRECISION :: k1, k2, k3, k4
  

  u_p_s(:) = u_p(:)
  
  call dz_dt(u_p_s, k1)
  call force(z_f, z_p, u_p_s, l1)

  call dz_dt(u_p_s + l1 / 2d0, k2)
  call force(z_f, z_p + k1 / 2d0, u_p_s + l1 / 2d0, l2)
  
  call dz_dt(u_p_s + l2 / 2d0, k3)
  call force(z_f, z_p + k2 / 2d0, u_p_s + l2 / 2d0, l3)
  
  call dz_dt(u_p_s + l3, k4)
  call force(z_f, z_p + k3, u_p_s + l3, l4)
 
  u_p(:) = u_p(:) + (l1 + 2d0 * l2 + 2d0 * l3 + l4) * d_t / 6d0
  z_p    = z_p    + (k1 + 2d0 * k2 + 2d0 * k3 + k4) * d_t / 6d0
  
   if (z_p <= 0.d0 ) then
      z_p    = - z_p
      u_p(0) =   DABS(u_p(0))     
      
   else if (z_p >= L_z) then
      z_p    = L_z - (z_p - L_z)
      u_p(0) = - DABS(u_p(0))
     
   end if
  
end subroutine runge_kutta
!--------------------------------------------------------------------------------------------
subroutine force(z, z_p, p, f) 

  implicit none
	
  DOUBLE PRECISION, PARAMETER   :: pi = 4d0*DATAN(1d0)
  DOUBLE PRECISION, INTENT(IN)  :: z(0 : n_z), z_p
  DOUBLE PRECISION, INTENT(IN)  :: p(0:2)
  DOUBLE PRECISION, INTENT(OUT) :: f(0:2)
  
  DOUBLE PRECISION :: gamma, B_p, dB_dz
  DOUBLE PRECISION :: ratio
  INTEGER :: i_z_left, i_z_right
	
  call p_to_gamma(p(0:2), gamma)
  call z_p_to_position(z_p, z, i_z_left, i_z_right, ratio)
 
  call z_to_B(z_p, B_p)
  call z_to_dB_dz(z_p, dB_dz)
  
	
  f(0) = - p(1)**2     * dB_dz / (2d0 * gamma * B_p)
  f(1) = + p(0) * p(1) * dB_dz / (2d0 * gamma * B_p) 
  f(2) = + B_p / gamma                               
  
end subroutine force

!-------------------------------------------------------------------------------------------------
subroutine z_to_dB_dz(z, dB_dz)
  use lshell_setting, only: r_eq
  implicit none
  DOUBLE PRECISION, INTENT(IN)  :: z
  DOUBLE PRECISION, INTENT(OUT) :: dB_dz
  DOUBLE PRECISION :: r, lambda

  call z_to_r_lambda(z, r, lambda)

  dB_dz = 3d0 * DSIN(lambda) / DCOS(lambda)**8 * (3d0 + 5d0 * DSIN(lambda)**2) / (1d0 + 3d0 * DSIN(lambda)**2) / r_eq
end subroutine


!----------------------------------------------------------------------------------------------
subroutine dz_dt(p, v)

  implicit none
  DOUBLE PRECISION, INTENT(IN)  :: p(0:2)
  DOUBLE PRECISION, INTENT(OUT) :: v
  DOUBLE PRECISION :: gamma
  
  call p_to_gamma(p(0:2), gamma)
  v = p(0) / gamma

end subroutine dz_dt

!-----------------------------------------------------------------------------------------------

subroutine p_to_gamma(p, gamma)

  implicit none

  DOUBLE PRECISION, INTENT(IN)  :: p(0:2)
  DOUBLE PRECISION, INTENT(OUT) :: gamma

  gamma = SQRT(1 + p(0)**2 + p(1)**2)

end subroutine p_to_gamma

!-------------------------------------------------------------------------------------------
subroutine energy_to_v(energy, gamma, v)

  implicit none

  DOUBLE PRECISION, INTENT(IN)  :: energy
  DOUBLE PRECISION, INTENT(OUT) :: gamma, v

  gamma = energy / m_e + 1d0
  v = DSQRT(1d0 - 1d0/gamma**2)
  

end subroutine energy_to_v
  
!-----------------------------------------------------------------------------------------------

subroutine z_p_to_position(z_p, z, i_z_left, i_z_right, ratio)
  ! z_p = (1d0 - ratio) * z(i_z_left) + ratio * z(i_z_right)
  implicit none	
  DOUBLE PRECISION, INTENT(IN)  :: z(0 : n_z), z_p
  INTEGER, INTENT(OUT)          :: i_z_left, i_z_right
  DOUBLE PRECISION, INTENT(OUT) :: ratio
  DOUBLE PRECISION :: diff(0 : n_z), d
  INTEGER :: i_min(1)
  
  diff  = ABS(z - z_p) 
  i_min = MINLOC(diff) !- ( n_z + 1 ) !array_z_0~2*n_z+1 
  
 
!    if (i_min(1) > n_z) then
!      write(*,*) '!-------------------------'
!      write(*,*) 'i_min = ', i_min(1)
!      write(*,*) '!-------------------------'
!      i_min(1) = n_z
!    end if

 if (i_min(1) >= n_z) then
    i_z_left  = n_z - 1
    i_z_right = n_z
    ratio = 1d0
    
 else if( i_min(1) <= -n_z) then
    i_z_left  = n_z - 1
    i_z_right = n_z
    ratio = 1d0

          
 else 
    d = z(i_min(1)) - z_p
    if (d > 0) then
       i_z_left  = i_min(1) - 1
       i_z_right = i_min(1)
    else if (d <= 0) then
       i_z_left  = i_min(1)
       i_z_right = i_min(1) + 1
    end if
    
    ratio     = (z_p - z(i_z_left)) / (z(i_z_right) - z(i_z_left))
 end if
end subroutine

!-----------------------------------------------------------------------------------------------


subroutine z_to_r_lambda(z, r, lambda)

  implicit none

  DOUBLE PRECISION, INTENT(IN)  :: z
  DOUBLE PRECISION, INTENT(OUT) :: r, lambda
  DOUBLE PRECISION :: f, g, lambda0, lambda1
  INTEGER :: i

  lambda0 = 1d0
  do i = 1, 1000000 
     if (i == 1000000) then
        write(*, *) "Error: solution is not found. z = ", z	
     end if
     f =  r_eq * ((1d0 / 2d0) * DSIN(lambda0) * DSQRT(3d0 * DSIN(lambda0)**2 + 1d0) &
          & + (1d0 / (2d0 * DSQRT(3d0))) * DLOG(DSQRT(3d0) * DSIN(lambda0) &
          & + DSQRT(3d0 * DSIN(lambda0)**2 + 1d0))) &
          & - z
     g = r_eq * DCOS(lambda0) * DSQRT(3d0 * DSIN(lambda0)**2 + 1d0)

     lambda1 = lambda0 - f / g
     if (DABS(lambda1 - lambda0) <=  1d-3) exit 
     lambda0 = lambda1
  end do
	
  lambda = lambda1
  r      = r_eq * DCOS(lambda)**2

end subroutine z_to_r_lambda

!-----------------------------------------------------------------------------------------------

subroutine z_to_B(z, B)

  implicit none
	
  DOUBLE PRECISION, INTENT(IN)  :: z
  DOUBLE PRECISION, INTENT(OUT) :: B
  DOUBLE PRECISION :: r, lambda
  
  call z_to_r_lambda(z, r, lambda)

  B = DSQRT(1d0 + 3d0 * DSIN(lambda)**2) / DCOS(lambda)**6	 

end subroutine z_to_B


!----------------------------------------------------------------------

end program distribution_calculater

