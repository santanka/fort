 subroutine force_for_PIC(t,w,z_p,z,n_z,p,B,Bw,psi,k,f,zeta,dB_dz)
   implicit none

   double precision, parameter  :: pi =4.d0*datan(1.d0)
   integer, intent(in)          :: n_z
   double precision, intent(in) :: z(0:n_z),B(0:n_z),psi(0:n_z),k(0:n_z),z_p
   double precision, intent(in) :: p(0:2),t,w,Bw
   double precision, intent(out):: f(0:2),zeta,dB_dz
   
   double precision :: gamma, B_p, psi_p, k_p, psi_p1
   double precision :: ratio, r_dB, dB_dz_left, dB_dz_right
   integer :: i_z_left, i_z_right

   call p_to_gamma(p(0:2), gamma)

   call z_p_to_position(z_p,z,n_z,i_z_left,i_z_right,ratio)

   B_p = (1.d0 -ratio)*B(i_z_left) + ratio*B(i_z_right)

   if (i_z_left == 0) then
      write(*,*) "Hajidesu! z_p =",z_p
      dB_dz_left  = (B(i_z_right) - B(i_z_left)) / (z(i_z_right) - z(i_z_left))
      dB_dz_right = (B(i_z_right) - B(i_z_left)) / (z(i_z_right) - z(i_z_left))
      r_dB = 1.d0    
    else if (ratio <= 0.5d0) then
    	dB_dz_left  = (B(i_z_right - 1) - B(i_z_left - 1)) / (z(i_z_right - 1) - z(i_z_left - 1))
    	dB_dz_right = (B(i_z_right    ) - B(i_z_left))     / (z(i_z_right    ) - z(i_z_left    ))
    	r_dB = ratio + 0.5d0
    else if (ratio > 0.5d0) then 
    	dB_dz_left  = (B(i_z_right    ) - B(i_z_left))     / (z(i_z_right    ) - z(i_z_left    ))
    	dB_dz_right = (B(i_z_right + 1) - B(i_z_left + 1)) / (z(i_z_right + 1) - z(i_z_left + 1))
    	r_dB = ratio - 0.5d0
    end if
    dB_dz = (1d0 - r_dB) * dB_dz_left + r_dB * dB_dz_right

    k_p   = (1.d0-ratio)*k(i_z_left)+ratio*k(i_z_right)
    psi_p = (1.d0-ratio)*psi(i_z_left)+ratio*psi(i_z_right)
    psi_p = w*t+psi_p
    zeta  = modulo(p(2)-psi_p, 2.d0*pi)

    f(0) = -p(1)**2  *dB_dz/(2.d0*gamma*B_p)+       p(1)/gamma *Bw*sin(zeta)
    f(1) =  p(0)*p(1)*dB_dz/(2.d0*gamma*B_p)+(w/k_p-p(0)/gamma)*Bw*sin(zeta)
    f(2) =  B_p/gamma             +1.d0/p(1)*(w/k_p-p(0)/gamma)*Bw*cos(zeta)

  end subroutine force_for_PIC
!-----------------------------------------------------------------------------  
subroutine dz_dt(p,v)
    implicit none
    double precision, intent(in) :: p(0:2)
    double precision, intent(out) :: v
    double precision ::gamma
    
    call p_to_gamma(p(0:2),gamma)
    v = p(0)/gamma
end subroutine dz_dt
!-----------------------------------------------------------------------------  
 subroutine p_to_gamma(p, gamma)
    implicit none
    double precision, intent(in)  :: p(0:2)
    double precision, intent(out) :: gamma

    gamma = sqrt(1.d0 + p(0)**2 + p(1)**2)
 end subroutine p_to_gamma
!----------------------------------------------------------------------------
 subroutine r_to_x(r,x)
   implicit none
   double precision, intent(in)  :: r(0:2)
   double precision, intent(out) :: x(0:2)

   x(0) = r(1)*dcos(r(2))
   x(1) = r(1)*dsin(r(2))
   x(2) = r(0)
 end subroutine r_to_x
!----------------------------------------------------------------------------
 subroutine z_to_B(z,B)
   implicit none
   double precision, intent(in) :: z
   double precision, intent(out) :: B
   double precision :: a,r,lambda

   call z_to_r_lambda(z, r, lambda)
   B = DSQRT(1d0 + 3d0 * DSIN(lambda)**2) / DCOS(lambda)**6
 end subroutine z_to_B
 
!----------------------------------------------------------------------------
 subroutine z_p_to_position(z_p,z,n_z,i_z_left,i_z_right,ratio)    
  implicit none
  double precision, intent(in)  :: z(0:n_z),z_p
  integer, intent(in)           :: n_z
  integer, intent(out)          :: i_z_left, i_z_right
  double precision, intent(out) :: ratio
  double precision              :: diff(0:n_z),d
  integer                       :: i_min(1)

  diff = abs(z-z_p)
  i_min = minloc(diff)

  d=z(i_min(1))-z_p
  if(d>0) then
     i_z_left  = i_min(1)-1
     i_z_right = i_min(1)
  else if (d<=0) then
     i_z_left  = i_min(1)
     i_z_right = i_min(1)+1
  end if
  ratio = (z_p - z(i_z_left)) / (z(i_z_right) - z(i_z_left))
end subroutine z_p_to_position    
!-------------------------------------------------------------------------
subroutine z_to_r_lambda(z,r,lambda)
  use lshell_setting
  implicit none
  double precision, intent(in) :: z
  double precision, intent(out) :: r,lambda
  double precision :: f,g,h,w,lambda0,lambda1
  integer :: i

  lambda = 1.d0
  do i=1,100000
     if(i==100000) then
        write(*,*) "Error: solution is not found. z =", z
     end if

     f = r_eq*(sin(lambda0)/2.d0*sqrt(1.d0+3.d0*sin(lambda0)**2) &
          & + 1.d0/2.d0/sqrt(3.d0) &
          & *log(sqrt(3.d0)*sin(lambda0)+sqrt(3.d0*sin(lambda0)**2)+1.d0))- z
     g = r_eq*cos(lambda0)*sqrt(3.d0*sin(lambda0)**2 + 1.d0)

     lambda1 = lambda0 - f/g
          if (abs(lambda1-lambda0)<= 1.d-3) exit
          lambda0 = lambda1
       end do

       lambda = lambda1
       r = r_eq*cos(lambda)**2
     end subroutine z_to_r_lambda
     
!-------------------------------------------------------------------------
subroutine resonance_velocity(w_p, w, Omega0, gamma, V_R)
  implicit none
  double precision, intent(in) :: w_p,w,Omega0,gamma
  double precision, intent(out) :: V_R
  double precision :: k,v,a

  call wave_number(w_p,w,Omega0,k)
  V_R = (w-Omega0/gamma)/k
  v = sqrt(1.d0-1.d0/gamma**2)
  if (abs(V_R) > v) then
     a = -10.0
     V_R = asin(a)               !???
  end if
end subroutine  
  
!-------------------------------------------------------------------------
subroutine wave_number(w_p,w,Omega0,k)
  implicit none
  double precision, intent(in) :: w_p,w,Omega0
  double precision, intent(out) :: k

  k = sqrt(w**2 + w*w_p**2 / (Omega0 - w))  !kaizen_dekiru

end subroutine wave_number
!--------------------------------------------------------------------------
  

