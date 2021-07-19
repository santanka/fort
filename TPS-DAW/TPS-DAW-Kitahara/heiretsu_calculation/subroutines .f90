
subroutine runge_kutta(z_f, phase_f, w_f, Bw_f, z_p, u_p, zeta_p, equator_flag)

  use constants_in_the_simulations, only: d_t, n_z, n_p, L_z, d_z
  implicit none

  DOUBLE PRECISION, INTENT(IN)    :: z_f(0:n_z), phase_f(0:n_z)
  DOUBLE PRECISION, INTENT(IN)    :: w_f(0:n_z), Bw_f(0:n_z)
  DOUBLE PRECISION, INTENT(OUT)   :: zeta_p
  DOUBLE PRECISION, INTENT(INOUT) :: z_p, u_p(0:2)
  LOGICAL, INTENT(OUT)            :: equator_flag

  DOUBLE PRECISION :: l1(0:2), l2(0:2), l3(0:2), l4(0:2), u_p_s(0:2)
  DOUBLE PRECISION :: k1, k2, k3, k4

    if (z_p >= z_f(n_z - 1) .and. u_p(0) >= 0d0) then ! o-kyu-sho-chi for reflective boundary
      u_p(0) = - DABS(u_p(0))
    end if 
    
    u_p_s(:) = u_p(:)

    call dz_dt(u_p_s, k1)
    call force(z_f, phase_f, w_f, Bw_f, z_p, u_p_s, l1, zeta_p)

    call dz_dt(u_p_s + l1 / 2d0, k2)
    call force(z_f, phase_f, w_f, Bw_f, z_p + k1 / 2d0, u_p_s + l1 / 2d0, l2, zeta_p)

    call dz_dt(u_p_s + l2 / 2d0, k3)
    call force(z_f, phase_f, w_f, Bw_f, z_p + k2 / 2d0, u_p_s + l2 / 2d0, l3, zeta_p)

    call dz_dt(u_p_s + l3, k4)
    call force(z_f, phase_f, w_f, Bw_f, z_p + k3, u_p_s + l3, l4, zeta_p)

    u_p(:) = u_p(:) + (l1 + 2d0 * l2 + 2d0 * l3 + l4) * d_t / 6d0
    z_p    = z_p    + (k1 + 2d0 * k2 + 2d0 * k3 + k4) * d_t / 6d0


    equator_flag = .false.

    if (z_p <= 0.d0) then
      z_p    = - z_p
      u_p(0) = - u_p(0)
      equator_flag = .true.
    else if (z_p >= L_z) then
      z_p    = L_z - (z_p - L_z)
      u_p(0) = - Dabs(u_p(0))
    end if

end subroutine

!-----------------------------------------------------------------------------------------------

subroutine write_time(string_in)
  implicit none
  CHARACTER(10) :: date1, date2, date3
  CHARACTER(20) :: string_out
  CHARACTER(20), OPTIONAL, INTENT(IN) :: string_in

  INTEGER       :: date_time(10)
  
  if (present(string_in)) then
    string_out = string_in
  else
    write(string_out, *) '     ' 
  end if
  
  call date_and_time(date1, date2, date3, date_time)
  write(*, '(I4, A1, I2.2, A1, I2.2, A1, I2.2, A1, I2.2, A1, I2.2, A1, A20)') &
  & date_time(1), '/', date_time(2), '/', date_time(3), ' ',         &
  & date_time(5), ':', date_time(6), ':', date_time(7), ' ',         &
  & string_out
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine eno(f0, mu, n_x, f1)
  implicit none

  INTEGER, INTENT(IN)  :: N_x
  DOUBLE PRECISION, INTENT(IN)  :: f0(-N_x:N_x), mu(-N_x:N_x) ! mu = velocity * dt / dx
  DOUBLE PRECISION, INTENT(OUT) :: f1(-N_x:N_x)

  INTEGER          :: i_t, i_x, time_INT, i, k, g
  DOUBLE PRECISION :: D2_k0, D2_k1, D3_k1, D3_k3
  DOUBLE PRECISION :: x, df, df_1, df_2, df_3, c_2, c_3
  
  ! update f(1, x)
  do i_x = - N_x + 3, N_x - 3
    if (mu(i_x) > 0d0) then
      k = i_x - 1
    else
      k = i_x
    end if

    df_1 = f0(k + 1) - f0(k)

    D2_k0 = (f0(k + 1) - 2.d0 * f0(k)     + f0(k - 1)) / 2.d0
    D2_k1 = (f0(k + 2) - 2.d0 * f0(k + 1) + f0(k))     / 2.d0

    if (DABS(D2_k0) <= DABS(D2_k1)) then
      c_2 = D2_k0
      g   = k - 1
    else
      c_2 = D2_k1
      g   = k
    end if

    df_2 = c_2 * DBLE(2 * (i_x - k) - 1)

    D3_k1 = (f0(g + 2) - 3.d0 * f0(g + 1) + 3.d0 * f0(g)     - f0(g - 1)) / 6.d0
    D3_k3 = (f0(g + 3) - 3.d0 * f0(g + 2) + 3.d0 * f0(g + 1) - f0(g))     / 6.d0

    if (DABS(D3_k1) <= DABS(D3_k3)) then
      c_3 = D3_k1
    else
      c_3 = D3_k3
    end if

    df_3 = c_3 * DBLE(3 * (i_x - g)**2 - 6 * (i_x - g) + 2)

    df = df_1 + df_2 + df_3

    f1(i_x) = f0(i_x) - mu(i_x) * df
  end do

end subroutine

!-----------------------------------------------------------------------------------------------
!subroutine runge_kutta(z_f, B0_f, k_f, phase_f, w_f, Bw_f, z_p, u_p, zeta_p, B0_p)

subroutine force(z, psi, freq, Bw, z_p, p, f, zeta)

  use constants_in_the_simulations, only: n_z, wave_existence, w_p

	implicit none
	
	DOUBLE PRECISION, PARAMETER   :: pi = 4d0*DATAN(1d0)
	DOUBLE PRECISION, INTENT(IN)  :: z(0 : n_z), psi(0 : n_z), z_p
	DOUBLE PRECISION, INTENT(IN)  :: p(0:2), freq(0 : n_z), Bw(0 : n_z)
	DOUBLE PRECISION, INTENT(OUT) :: f(0:2)

	DOUBLE PRECISION :: gamma, B_p, dB_dz, zeta, psi_p, k_p, Bw_p, w, fw(0:2)
	DOUBLE PRECISION :: ratio
	INTEGER :: i_z_left, i_z_right
	
	call p_to_gamma(p(0:2), gamma)
	call z_p_to_position(z_p, z, i_z_left, i_z_right, ratio)
		
	call z_to_B(z_p, B_p)
	call z_to_dB_dz(z_p, dB_dz)

	if (wave_existence .eqv. .true.) then	
!	  k_p   = (1d0 - ratio) * k(i_z_left)    + ratio * k(i_z_right)	
	  psi_p = (1d0 - ratio) * psi(i_z_left)  + ratio * psi(i_z_right)
      w     = (1d0 - ratio) * freq(i_z_left) + ratio * freq(i_z_right)
	  call wave_number(w_p, w, B_p, k_p)

	  zeta = MODULO(p(2) - psi_p, 2d0*pi)
	
	  Bw_p = (1d0 - ratio) * Bw(i_z_left)   + ratio * Bw(i_z_right)
	  
	  fw(0) = p(1)/gamma  * Bw_p * DSIN(zeta)
	  fw(1) = (w/k_p - p(0)/gamma) * Bw_p * DSIN(zeta)
	  fw(2) = 1d0 / p(1) * (w/k_P - p(0)/gamma) * Bw_p * DCOS(zeta)
	else
	  fw(:) = 0d0
	end if
	
	f(0) = - p(1)**2     * dB_dz / (2d0 * gamma * B_p) + fw(0)
	f(1) = + p(0) * p(1) * dB_dz / (2d0 * gamma * B_p) + fw(1)
	f(2) = + B_p / gamma                               + fw(2)
end subroutine

!-----------------------------------------------------------------------------------------------

!subroutine force_from_wave(z, p, f, t)
!
!implicit none
!DOUBLE PRECISION, INTENT(IN)  :: z, p(0:2), t
!DOUBLE PRECISION, INTENT(OUT) :: f(0:2)
!DOUBLE PRECISION :: gamma, zeta
!
!DOUBLE PRECISION, PARAMETER   :: w_p = 2d0    ! plasma frequency
!DOUBLE PRECISION, PARAMETER   :: w   = 0.5d0  ! wave frequency
!DOUBLE PRECISION, PARAMETER   :: k   = DSQRT(w**2 + w * w_p**2 / (1d0 - w)) ! = 2.06d0 ! wave number
!DOUBLE PRECISION, PARAMETER   :: Bw  = 1d-4
!
!call p_to_gamma(p(0:2), gamma)
!call relative_phase(z, p, t, zeta)
!
!f(0) =                     p(1)/gamma  * Bw * DSIN(zeta)
!f(1) =              (w/k - p(0)/gamma) * Bw * DSIN(zeta)
!f(2) = 1d0 / p(1) * (w/k - p(0)/gamma) * Bw * DCOS(zeta)
!
!end subroutine

!-----------------------------------------------------------------------------------------------

subroutine z_to_dB_dz(z, dB_dz)
use lshell_setting, only: r_eq
implicit none
DOUBLE PRECISION, INTENT(IN)  :: z
DOUBLE PRECISION, INTENT(OUT) :: dB_dz
DOUBLE PRECISION :: r, lambda

call z_to_r_lambda(z, r, lambda)

dB_dz = 3d0 * DSIN(lambda) / DCOS(lambda)**8 * (3d0 + 5d0 * DSIN(lambda)**2) / (1d0 + 3d0 * DSIN(lambda)**2) / r_eq
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine dz_dt(p, v)

implicit none
DOUBLE PRECISION, INTENT(IN)  :: p(0:2)
DOUBLE PRECISION, INTENT(OUT) :: v
DOUBLE PRECISION :: gamma

call p_to_gamma(p(0:2), gamma)
v = p(0) / gamma
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine p_to_gamma(p, gamma)
implicit none

DOUBLE PRECISION, INTENT(IN)  :: p(0:2)
DOUBLE PRECISION, INTENT(OUT) :: gamma

gamma = SQRT(1 + p(0)**2 + p(1)**2)
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine p_to_energy(p, gamma, energy)
  use constant_parameter, only : m_e

  implicit none

  DOUBLE PRECISION, INTENT(IN)  :: p(0:2)
  DOUBLE PRECISION, INTENT(OUT) :: energy
  DOUBLE PRECISION :: gamma

  gamma  = DSQRT(1 + p(0)**2 + p(1)**2)
  energy = m_e * (gamma - 1.0)
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine r_to_x(r, x)
implicit none

DOUBLE PRECISION, INTENT(IN)  :: r(0:2)
DOUBLE PRECISION, INTENT(OUT) :: x(0:2)

x(0) = r(1) * dcos(r(2))
x(1) = r(1) * dsin(r(2))
x(2) = r(0)
end subroutine
!-----------------------------------------------------------------------------------------------

subroutine z_to_B(z, B)
	implicit none
	
	DOUBLE PRECISION, INTENT(IN)  :: z
	DOUBLE PRECISION, INTENT(OUT) :: B
	DOUBLE PRECISION :: r, lambda

	call z_to_r_lambda(z, r, lambda)
	B = DSQRT(1d0 + 3d0 * DSIN(lambda)**2) / DCOS(lambda)**6	 
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine z_p_to_position(z_p, z, i_z_left, i_z_right, ratio)
  use constants_in_the_simulations, only: n_z
	! z_p = (1d0 - ratio) * z(i_z_left) + ratio * z(i_z_right)
	implicit none	
	DOUBLE PRECISION, INTENT(IN)  :: z(0 : n_z), z_p 
	INTEGER, INTENT(OUT)          :: i_z_left, i_z_right
	DOUBLE PRECISION, INTENT(OUT) :: ratio 
	DOUBLE PRECISION :: diff(0 : n_z), d
	INTEGER :: i_min(1)
	
	diff  = ABS(z - z_p)
	i_min = MINLOC(diff)

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
  use lshell_setting

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
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine resonance_velocity(w_p, w, Omega0, gamma, V_R)
  implicit none

  DOUBLE PRECISION, INTENT(IN)  :: w_p, w, Omega0, gamma
  DOUBLE PRECISION, INTENT(OUT) :: V_R
  DOUBLE PRECISION :: k!, v, a
	
  call wave_number(w_p, w, Omega0, k)
  V_R = (w - Omega0 / gamma) / k

!  v   = DSQRT(1d0 - 1d0 / gamma**2)  
!  if (DABS(V_R) > v) then
!    a = -10.0
!    V_R = ASIN(a)
!  end if
end subroutine

!-----------------------------------------------------------------------------------------------

subroutine vpara_and_vabs_to_alpha(v_para, v_abs, alpha)
  implicit none

  REAL(8), INTENT(IN)  :: v_para, v_abs
  REAL(8), INTENT(OUT) :: alpha
	
  alpha = DACOS(v_para / v_abs)
end subroutine
!
!!-----------------------------------------------------------------------------------------------
!
subroutine gamma_to_vabs(gamma, v_abs)
  implicit none

  DOUBLE PRECISION, INTENT(IN)  :: gamma
  DOUBLE PRECISION, INTENT(OUT) :: v_abs
	
  v_abs = DSQRT(1d0 - 1d0 / gamma**2)
end subroutine
!
!!-----------------------------------------------------------------------------------------------

subroutine wave_number(w_p, w, B, k)
	implicit none

	DOUBLE PRECISION, INTENT(IN)  :: w_p, w, B
	DOUBLE PRECISION, INTENT(OUT) :: k
	
	k = DSQRT(w**2 + w * w_p**2 / (B - w))

end subroutine


!!-----------------------------------------------------------------------------------------------

subroutine group_velocity(w, B, w_p, V_g)

  implicit none
  DOUBLE PRECISION, INTENT(IN)  :: w, B, w_p
  DOUBLE PRECISION, INTENT(OUT) :: V_g

!  V_g = 2.0 * k / (2.0 * w + B * w_p**2 / (B - w)**2)
  V_g = 2.0 * DSQRT(w**2 + w * w_p**2 / (B - w)) / (2.0 * w + B * w_p**2 / (B - w)**2)
end subroutine


!!-----------------------------------------------------------------------------------------------
