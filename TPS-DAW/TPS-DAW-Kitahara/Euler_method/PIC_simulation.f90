program PIC_at_Low_PA
  use constant_parameter

  implicit none

  !-------------------------------------------------------------------
  ! initial setting of simulation system
  !------------------------------------------------------------------
  integer, parameter          :: n_t = 1000000     !time_number
  integer, parameter          :: n_z = 100000      !grid_number
  double precision, parameter :: d_t = 0.1d0    !time_interval
  double precision, parameter :: d_z = 0.1d0    !grid_interval
  integer, parameter          :: n_particle = 4


  !--------------------------------------------------------------------
  ! initial setting of wave
  !--------------------------------------------------------------------
  double precision, parameter :: w_p = 4.d0
  double precision, parameter :: w   = 0.4d0
  double precision, parameter :: Bw  = 5.d-4

  !--------------------------------------------------------------------
  ! initial setting of particle
  !-------------------------------------------------------------------
  double precision, parameter :: z0 = 10.d0
  integer, parameter :: gamma_center = 2.d0 !Initial energy [keV]
  integer, parameter :: alpha_center = 120  !initial pitch angle [deg]
  double precision   :: v0,v_para0,v_perp0,gamma0,alpha0

  !---------------------------------------------------------------------
  ! variables
  !---------------------------------------------------------------------
  integer          :: i_energy,i_alpha,i_t,i_phi,i_z,unit2,unit3
  double precision :: B(0:n_z),z(0:n_z),k(0:n_z),psi_k(0:n_z)
  double precision :: r0(0:2),p0(0:2)
  double precision :: r1(0:2),p1(0:2),x(0:2),f_wave(0:2)
  double precision :: p_abs,alpha,alpha1,gamma1,time,zeta0,zeta1
  double precision :: V_R,B0,alpha_R,dB_dz
  character(64)    :: file_output,file_zeta,file_wave

  !---------------------------------------------------------------------
  ! for Runge-Kutta
  !---------------------------------------------------------------------
  double precision :: l1(0:2),l2(0:2),l3(0:2),l4(0:2)
  double precision :: k1,k2,k3,k4,B1,ratio,k_p
  integer          :: i_z_left,i_z_right


  !---------------------------------------------------------------------
  ! for time measurement
  !---------------------------------------------------------------------
  character(10) :: date1,date2,date3
  integer       :: date_time(10)
  character(10) :: command


  !---------------------------------------------------------------------
  ! program start
  !---------------------------------------------------------------------
  
  !-------------------------------------------------------------------
  ! start time
  !-------------------------------------------------------------------
  i_energy = 0.d0
  i_alpha  =0.d0


  write(*,*) 'w_p = ',w_p
  write(*,*) 'w   = ',w
  write(*,*) 'Bw  = ',Bw

  gamma0  = (dble(gamma_center) + dble(i_energy)) / m_e+1.d0
  alpha0  = (dble(alpha_center) + dble(i_alpha)) * deg2rad
  v0      = sqrt(1.d0 -1.d0/gamma0**2)
  v_para0 = v0*cos(alpha0)
  v_perp0 = v0*sin(alpha0)


write(*,*) 'gamma0  = ',gamma0
write(*,*) 'alpha0  = ',alpha0
write(*,*) 'v0      = ',v0
write(*,*) 'v_para0 = ',v_para0
write(*,*) 'v_perp0 = ',v_perp0


  call date_and_time(date1,date2,date3,date_time)
  write(*,*)'--------------------------------'
  write(*,*)'a = ',alpha_center + i_alpha,'K = ',gamma_center+i_energy
  write(*,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
       & date_time(1),'/',date_time(2),'/',date_time(3),' ',&
       & date_time(5),':',date_time(6),':',date_time(7)

  
  !----------------------------------------
  ! initial setting of field (grid)
  !----------------------------------------

  write(command, '(A7, I3.3, A2, I3.3)') 'mkdir K', gamma_center + i_energy, '_a', alpha_center + i_alpha
  call system(command)

  write(file_wave, '(a1,i3.3,a2,i3.3,a9)')'K',gamma_center + i_energy,' _a',alpha_center + i_alpha, '/wave.dat'

  open(10,file=file_wave)

  call z_to_B(z0,B0)
  
  do i_z = 0,n_z
     z(i_z) = dble(i_z)*d_z

     call z_to_B(z(i_z),B(i_z))
     call resonance_velocity(w_p,w,B(i_z),gamma0,V_R)
     call wave_number(w_p,w,B(i_z),k(i_z))  !kaeta
     alpha_R = acos(V_R/v0)*rad2deg

     if (i_z==0) then
        psi_k(i_z) = 0.d0
     else
        psi_k(i_z) = psi_k(i_z-1)-(k(i_z-1)+k(i_z))/2*d_z
     end if

     alpha = sin(alpha0)**2/B0*B(i_z)
     alpha = sqrt(alpha)
     alpha = sign(alpha, cos(alpha0))
     alpha = 180+asin(alpha)*rad2deg


     write(10,'(14E15.7)') z(i_z),B(i_z),alpha,k(i_z),psi_k(i_z),V_R,v0,alpha_R
  end do

  close(10)

!-------------------------------------------
! trace a particle
!-------------------------------------------
!
!$omp parallel private(r0, r1, p0, p1, file_output, file_zeta, unit2, unit3, &
!$omp & date1, date2, date3, date_time, zeta1, zeta0, dB_dz, k_p, &
!$omp & i_t, time, l1, l2, l3, l4, k1, k2, k3, k4, B1, i_z_left, i_z_right, ratio, alpha1, gamma1)
!$omp do
  do i_phi = 1,n_particle

     unit2 = 10  + i_phi
     unit3 = 100 - i_phi

     r0(0) = z0
     r0(1) = 0.d0
     r0(2) = 0.d0
     p0(0) = gamma0*v_para0
     p0(1) = gamma0*v_perp0
     p0(2) = 2.d0*pi/dble(n_particle)*dble(i_phi-1)
     zeta0 = p0(2)-psi_k(n_z)

     write(file_output,'(a1,i3.3,a2,i3.3,a7,i2.2,a4)') &
          &'K',gamma_center+i_energy,'_a',alpha_center+i_alpha,'/output',i_phi,'.dat'
     write(file_zeta,'(a1,i3.3,a2,i3.3,a5,i2.2,a4)') &
          &'K',gamma_center+i_energy,'_a',alpha_center+i_alpha,'/zeta',i_phi,'.dat'
     open(unit = unit2, file = file_output)
     open(unit = unit3, file = file_zeta)

     do i_t = 1,n_t

        !-------------------------------------------
        ! start : runge_kutta
        !-------------------------------------------

        time = dble(i_t)*d_t

        call dz_dt(p0,k1)
        call force_for_PIC(time,w,r0(0),z,n_z,p0,B,Bw,psi_k,k,l1,zeta1,dB_dz)

        call dz_dt(p0+l1/2.d0,k2)
        call force_for_PIC(time,w,r0(0)+k1/2.d0,z,n_z,p0+l1/2.d0,B,Bw,psi_k,k,l2,zeta1,dB_dz)

        call dz_dt(p0+l2/2.d0,k3)
        call force_for_PIC(time,w,r0(0)+k1/2.d0,z,n_z,p0+l1/2.d0,B,Bw,psi_k,k,l2,zeta1,dB_dz)

        call dz_dt(p0+l3,k4)
        call force_for_PIC(time,w,r0(0)+k3,z,n_z,p0+l3,B,Bw,psi_k,k,l4,zeta1,dB_dz)

        p1    = p0    + (l1 + 2.d0 * l2 + 2.d0 * l3 + l4) * d_t/6.d0
        r1(0) = r0(0) + (k1 + 2.d0 * k2 + 2.d0 * k3 + k4) * d_t/6.d0

        if (r1(0) < 1.d0 .or. r1(0) > z0) then
           write(*,*) ' z (dist.) is = ',r1(0)
           write(*,*) ' time step is = ',i_t
        end if
        if (r1(0) < 1.d0 .or. r1(0) > z0) exit
        
        call z_p_to_position(r1(0), z, n_z, i_z_left, i_z_right, ratio)
        B1    = (1.d0-ratio)*B(i_z_left) + ratio*B(i_z_right)
        k_p   = (1.d0-ratio)*k(i_z_left) + ratio*k(i_z_right)
        r1(1) = p1(1) / B1
        r1(2) = p1(2) - pi/2.d0

        !-------------------------------------
        ! end : runge_kutta
        !-------------------------------------

        alpha1 = atan2(p1(1),p1(0))*rad2deg

        call p_to_gamma(p1,gamma1)

        call force_for_PIC(time, w, r1(0), z, n_z, p1, B, Bw, psi_k, k, l1, zeta1, dB_dz)

        write(unit = unit2, fmt='(14e15.7)') time, r1(0), alpha1, p1(0), p1(1), p1(2), modulo(p1(2),2.d0*pi), &
            &  B1, gamma1, zeta1, dB_dz, k_p

        if( (zeta0<pi/8 .and. zeta1>2.0d0*pi-pi/8.d0).or.(zeta0>2.d0*pi-pi/8.d0 .and. zeta1<pi/8) ) then
        write(unit = unit3, fmt = *) ''
        end if
        write(unit = unit3, fmt='(14e15.7)') time, r1(0), alpha1, p1(0), p1(1), p1(2), modulo(p1(2),2.d0*pi), &
            &  B1, gamma1, zeta1, dB_dz, k_p


        r0 = r1
        p0 = p1
        zeta0 = zeta1
     end do


    close(unit = unit2)
    write(*, '(a15,i2.2,a18)')'Calculation of ',i_phi, 'th particle is end.'
    call date_and_time(date1, date2, date3, date_time)
    write(*, '(i4, a1, i2.2, a1, i2.2, a1, i2.2, a1, i2.2, a1, i2.2)') &
         & date_time(1), '/', date_time(2), '/', date_time(3), ' ',date_time(5), ':', date_time(6), ':', date_time(7)
 end do 
!$omp end do
!$omp end parallel
 
end program PIC_at_Low_PA

    








        
