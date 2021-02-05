program SVC_Temperature

implicit none

!�S�̂�ʂ��Ă̕���
integer :: i, s, j   !do�\���p
integer, parameter :: N = 219 !grid��
integer, parameter :: Z = 180 !���x���grid��
integer, parameter :: kind = 6 !���q�퐔
integer, parameter :: spece = 3 !�^�E���q�퐔
real*8, parameter :: pi = 4.d0*atan(1.d0) !�~����
character(len=80) :: dummy !�g�p���Ȃ������p

!mag_FA�֘A
real*8, dimension(N) :: lambda, BB, dteq

!SVC_parameter�֘A
real*8 :: ee, cc, mu0, GG, alpha, ME, MassE, RE, omE, Req, MassI, ep0

!SVC_IC�֘A
real*8, dimension(N) :: ss, sR
real*8, dimension(N-1) :: ds !grid��

!SVC_BC�֘A
character(len=1), dimension(kind) :: so !�N��
real*8, dimension(kind) :: Tp, Ta, ch, mass
integer, dimension(kind) :: kk !���q��
integer, dimension(kind) :: sigmab !�N����grid

!SVC_min�֘A
real*8, dimension(N) :: pe !�Ód�|�e���V����
real*8, dimension(kind) :: sigma !���E�ł̐����x

!subroutine�A���p
real*8, dimension(kind, Z) :: mu !�f�M�s�ϗ�
real*8, dimension(N, kind) :: UU !�ʒu�G�l���M�[
real*8, dimension(N, kind, Z, 3) :: a3 !accessibility(1:min, 2:lim, 3:max)
real*8 :: ff !���z�֐�(function)
real*8, dimension(N, kind, Z) :: theta1 !a�ϕ�(�����x)
real*8, dimension(N, kind) :: num !�����x
real*8, dimension(N) :: rhov !�d�ז��x
real*8, dimension(N) :: rhop !Poisson�������̓d�ז��x
real*8, dimension(N) :: cvg !�����l
real*8, dimension(N, 3) :: VA !Alfven���x
real*8, dimension(N, kind, Z) :: theta2 !a�ϕ�(���ϗ���)
real*8, dimension(N, kind) :: aave !���ϗ���(a)
real*8, dimension(N, kind, Z, 2) :: theta3 !a�ϕ�(���x)
real*8, dimension(N, kind) :: Tperp, Tpara !���x[eV]
real*8, dimension(N, spece) :: taave !�^�E���ϗ���
real*8, dimension(N, kind, Z, 2) :: theta4 !a�ϕ�(�^�E���x)
real*8, dimension(N, spece) :: tTperp, tTpara, tT !�^�E���x
real*8, dimension(N, spece) :: bs !�e���l
real*8, dimension(N) :: beta !���l
real*8, dimension(N) :: tTall !���x



!mag_FA�̒��o
open(50, file="mag_FA_E_L=4_NS.csv", action="read", status="old")
do i = 1, N
 read(50, *) lambda(i), BB(i), dteq(i)
enddo !i
close(50)


!SVC_parameter�̒��o
open(51, file="SVC_parameter_E_L=4.csv", action="read", status="old")
read(51, *) dummy, ee
read(51, *) dummy, cc
read(51, *) 
read(51, *) dummy, mu0
read(51, *) dummy, GG
read(51, *) dummy, alpha
read(51, *) dummy, ME
read(51, *) dummy, MassE
read(51, *) dummy, RE
read(51, *) dummy, omE
read(51, *) dummy, Req
read(51, *) dummy, MassI
close(51)

ep0 = 1.d0/mu0/cc**2.d0
alpha = alpha*3.d0


!SVC_IC�̒��o
open(52, file="SVC_IC_E_L=4_NS.csv", action="read", status="old")
read(52, *) !1�s�Ƃ΂�
do i = 1, N
 read(52, *) dummy, dummy, ss(i), sR(i)
enddo !i
close(52)

do i = 1, N-1
 ds(i) = ss(i+1) - ss(i) !grid�ԋ���
enddo !i


!SVC_BC�̒��o
open(53, file="SVC_BC_E_L=4_NS.csv", action="read", status="old")
read(53, *) !1�s�Ƃ΂�
do i = 1, kind
 read(53, *) dummy, dummy, Tp(i), Ta(i), so(i), ch(i), mass(i), kk(i)
enddo !i
close(53)

Tp = Tp*ee
ch = ch*ee

do i = 1, kind
 if(so(i) == "I") sigmab(i) = N
 if(so(i) == "M") sigmab(i) = 1
enddo !i


!SVC_min�̒��o
open(54, file="SVC_min_E_L=4_NS.csv", action="read", status="old")
read(54, *) !1�s�Ƃ΂�
do i = 1, N
 read(54, *) pe(i)
enddo !i

do i = 1, kind
 read(54, *) sigma(i)
enddo !i
close(54)


!�f�M�s�ϗ�
call AI(N, Z, kind, alpha, Tp, BB, sigmab, mu)

!�ʒu�G�l���M�[
call EP(N, kind, GG, MassE, Req, omE, MassI, pe, lambda, dteq, ch, mass, UU)

!accessibility(1:min, 2:lim, 3:max)
call access(N, Z, kind, alpha, Tp, Ta, BB, UU, mu, sigmab, a3)

!a�ɂ��Đϕ�(�����x)
call TH1(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta1)

!mu�ɂ��Đϕ�(�����x)
call NM1(N, kind, Z, BB, mass, mu, theta1, num)

!�d�ז��x
call RR(N, kind, ch, num, rhov)

!Poisson������
call Poisson(N, ep0, pe, ds, rhop)

!�����`�F�b�N
call CV(N, kind, rhov, rhop, ch, num, cvg)

!Alfven���x
call Alf(N, kind, mu0, cc, BB, ch, mass, num, VA)

!a�ɂ��Đϕ�(���ϗ���)
call TH2(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta2)

!mu�ɂ��Đϕ�(���ϗ���)
call NM2(N, kind, Z, BB, mass, mu, theta2, num, aave)

!a�ɂ��Đϕ�(���x)
call TH3(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, aave, a3, theta3)

!mu�ɂ��Đϕ�(���x)
call NM3(N, kind, Z, ee, BB, mass, mu, theta3, num, Tperp, Tpara)

!���q�킲�Ƃ̗���
call ryusoku(N, kind, spece, aave, num, UU, kk, sigmab, taave)

!���q�킲�Ƃ̉��x(a�ϕ�)
call TH4(N, kind, Z, spece, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, kk, taave, a3, theta4)

!���q�킲�Ƃ̉��x(mu�ϕ�)
call NM4(N, kind, Z, spece, ee, mu0, BB, mass, mu, theta4, num, UU, kk, sigmab, tTperp, tTpara, tT, bs, beta, tTall)


!�t�@�C����
open(60, file="SVC_E_L=4_NS_all.csv")
do i = 1, N
 write(60, 70) ss(i), sR(i), lambda(i), dteq(i), BB(i), pe(i), UU(i, :), num(i, :), rhov(i), rhop(i), cvg(i), &
  & VA(i, :), aave(i, :), Tperp(i, :), Tpara(i, :), taave(i, :), tTperp(i, :), tTpara(i, :), tT(i, :), tTall(i), &
  & bs(i, :), beta(i)
enddo !i
close(60)

open(61, file="SVC_E_L=4_NS_a3_all.csv")
do s = 1, kind
 do j = 1, Z
  write(61, 71) "s = ", s, "j = ", j
  do i = 1, N
   write(61, 72) "i = ", i, a3(i, s, j, :)
  enddo !i
  write(61, *) "  "
  write(61, *) "  "
 enddo !j
enddo !s
close(61)

70 format(1PE25.15E3, 58(',', 1PE25.15E3)) !5*kind+5*spece+14
71 format(A4, I2, ',', 2x, A4, I2)
72 format(A4, I3, ',' 1PE25.15E3, 2(',', 1PE25.15E3))

end program SVC_Temperature





!subroutine & function

!�f�M�s�ϗ�
subroutine AI(N, Z, kind, alpha, Tp, BB, sigmab, mu)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp
 real*8, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(kind, Z), intent(out) :: mu
 
 integer :: j, s
 
 
 do s = 1, kind
  do j = 1, Z
   mu(s, j) = 1.d-20*(alpha*Tp(s)/BB(sigmab(s))/1.d-20)**(dble(j-1)/dble(Z-1))
  enddo !j
 enddo !s
 
 return
 
end subroutine AI


!�ʒu�G�l���M�[
subroutine EP(N, kind, GG, MassE, Req, omE, MassI, pe, lambda, dteq, ch, mass, UU)
 implicit none
 integer, intent(in) :: N, kind
 real*8, intent(in) :: GG, MassE, Req, omE, MassI
 real*8, dimension(N), intent(in) :: pe, lambda, dteq
 real*8, dimension(kind), intent(in) :: ch, mass
 
 real*8, dimension(N, kind), intent(out) :: UU
 
 integer :: i, s
 
 
 do i = 1, N
  do s = 1, kind
   !�Ód�|�e���V����
   UU(i, s) = ch(s)*pe(i)
   !�f���d��
   UU(i, s) = UU(i, s) - GG*MassE*mass(s)/Req/cos(lambda(i))**2.d0
   !�f�����S��
   UU(i, s) = UU(i, s) - 5.d-1*mass(s)*omE**2.d0*(Req*cos(lambda(i))**3.d0)**2.d0
   !�q���d��
   if(MassI /= 0.d0) then
     UU(i, s) = UU(i, s) - GG*MassI*mass(s)/dteq(i)
   endif
  enddo !s
 enddo !i
 
 return
 
end subroutine EP


!accessibility
subroutine access(N, Z, kind, alpha, Tp, Ta, BB, UU, mu, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp, Ta
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(N, kind, Z, 3), intent(out) :: a3
 
 real*8, dimension(N, kind, Z) :: UB
 integer, dimension(N, kind, Z) :: lima
 integer :: i, s, j, k
 integer :: t, r !�t���[
 real*8 :: f !�t���[
 
 
 !UB�쐬
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    UB(i, s, j) = UU(i, s) + mu(s, j)*BB(i)
   enddo !j
  enddo !s
 enddo !i
 
 
 !amin
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"�N��
    if(sigmab(s) == 1) then
      t = 1
      do k = 1, i
       if(UB(k, s, j) > UB(t, s, j)) t = k
      enddo !k
      a3(i, s, j, 1) = sqrt(UB(t, s, j) - UB(i, s, j))
    endif
    
    !"I"�N��
    if(sigmab(s) == N) then
      t = N
      do k = i, N
       if(UB(k, s, j) > UB(t, s, j)) t = k
      enddo !k
      a3(i, s, j, 1) = sqrt(UB(t, s, j) - UB(i, s, j))
    endif
    
    if(a3(i, s, j, 1) == 0.d0) a3(i, s, j, 1) = 1.d-17
    
   enddo !j
  enddo !s
 enddo !i
 
 
 !alim
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"�N��
    if(sigmab(s) == 1) then
      if(i == 1) a3(i, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
      if(i /= 1) then
        f = UB(1, s, j) + alpha*Tp(s)*Ta(s) - UB(i, s, j)
        if(f <= 0.d0) a3(i, s, j, 2) = 0.d0
        if(f > 0.d0) a3(i, s, j, 2) = sqrt(f)
        if(a3(i-1, s, j, 2) < a3(i-1, s, j, 1)) a3(i, s, j, 2) = 0.d0
      endif
    endif
    
    !"I"�N��
    if(sigmab(s) == N) then
      k = N+1-i
      if(k == N) a3(k, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
      if(k /= N) then
        f = UB(N, s, j) + alpha*Tp(s)*Ta(s) - UB(k, s, j)
        if(f <= 0.d0) a3(k, s, j, 2) = 0.d0
        if(f > 0.d0) a3(k, s, j, 2) = sqrt(f)
        if(a3(k+1, s, j, 2) < a3(k+1, s, j, 1)) a3(k, s, j, 2) = 0.d0
      endif
    endif
    
   enddo !j
  enddo !s
 enddo !i
 
 
 !accessibility�̒���
 do s = 1, kind
  do j = 1, Z
   
   !"M"�N��
   if(sigmab(s) == 1) then
     do i = 1, N
      if(i /= 1 .and. lima(i-1, s, j) == 1) lima(i, s, j) = 1
      if(a3(i, s, j, 2) <= a3(i, s, j, 1)) lima(i, s, j) = 1
     enddo !i
   endif
   
   !"I"�N��
   if(sigmab(s) == N) then
     do i = 1, N
      k = N+1-i
      if(k /= N .and. lima(k+1, s, j) == 1) lima(k, s, j) = 1
      if(a3(k, s, j, 2) <= a3(k, s, j, 1)) lima(k, s, j) = 1
     enddo !i
   endif
   
  enddo !j
 enddo !s
 
 
 !amax
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"�N��
    if(sigmab(s) == 1) then
      if(i == N .or. lima(i, s, j) == 1) a3(i, s, j, 3) = 0.d0
      if(i /= N .and. lima(i, s, j) == 0) then
        t = N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        r = N
        do k = i+1, N
         if(lima(k-1, s, j) == 0 .and. lima(k, s, j) == 1) r = k
        enddo !k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k = i+1, r
         if(UB(k, s, j) > UB(t, s, j)) t = k
        enddo !k
        a3(i, s, j, 3) = UB(t, s, j) - UB(i, s, j)
        if(a3(i, s, j, 3) < 0.d0) a3(i, s, j, 3) = 0.d0
        if(a3(i, s, j, 3) >= 0.d0) a3(i, s, j, 3) = sqrt(a3(i, s, j, 3))
      endif
    endif
    
    !"I"�N��
    if(sigmab(s) == N) then
      if(i == 1 .or. lima(i, s, j) == 1) a3(i, s, j, 3) = 0.d0
      if(i /= 1 .and. lima(i, s, j) == 0) then
        t = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        r = 1
        do k = 1, i-1
         if(lima(k, s, j) == 1 .and. lima(k+1, s, j) == 0) r = k
        enddo !k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k = r, i-1
         if(UB(k, s, j) > UB(t, s, j)) t = k
        enddo !k
        a3(i, s, j, 3) = UB(t, s, j) - UB(i, s, j)
        if(a3(i, s, j, 3) < 0.d0) a3(i, s, j, 3) = 0.d0
        if(a3(i, s, j, 3) >= 0.d0) a3(i, s, j, 3) = sqrt(a3(i, s, j, 3))
      endif
    endif
    
   enddo !j
  enddo !s
 enddo !i
 
 return
 
end subroutine access


!���z�֐�
real*8 function ff(sigma, BB, B1, Tp, Ta, UU, U1, mu, aa, mass)
 implicit none
 real*8, intent(in) :: sigma, BB, B1, Tp, Ta, UU, U1, mu, aa, mass
 
 real*8, parameter :: pi = 4.0*atan(1.d0)
 
 
 ff = sigma*sqrt((mass/2.d0/pi/Tp)**3.d0/Ta)*exp(-B1*mu/Tp)*exp(-(UU+BB*mu+aa**2.d0-(U1+B1*mu))/Ta/Tp)
 return
 
end function ff



!a�ɂ��Đϕ�(�����x)
subroutine TH1(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta1)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z), intent(out) :: theta1
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !���z�֐�(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) >= a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       thetaL = thetaL + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
        & / 2.d0 * abs(aL(p+1) - aL(p))
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) >= a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       thetaM = thetaM + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
        & / 2.d0 * abs(aM(p+1) - aM(p))
      enddo !p
    endif
    
    theta1(i, s, j) = thetaL + thetaM
   enddo !j
  enddo !s
 enddo !i 
 
 return
 
end subroutine TH1


!mu�ɂ��Đϕ�(�����x)
subroutine NM1(N, kind, Z, BB, mass, mu, theta1, num)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z), intent(in) :: theta1
 
 real*8, dimension(N, kind), intent(out) :: num
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 do i = 1, N
  do s = 1, kind
   
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta1(i, s, j) + theta1(i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   num(i, s) = sqrt((2.d0/mass(s))**3.d0)*BB(i)*pi*nnn
   
  enddo !s
 enddo !i
 
 return
 
end subroutine NM1


!�d�ז��x
subroutine RR(N, kind, ch, num, rhov)
 implicit none
 integer, intent(in) :: N, kind
 real*8, dimension(kind) :: ch
 real*8, dimension(N, kind) :: num
 
 real*8, dimension(N), intent(out) :: rhov
 
 integer :: i, s
 
 
 rhov = 0.d0
 
 do i = 1, N
  do s = 1, kind
   rhov(i) = rhov(i) + ch(s) * num(i, s)
  enddo !s
 enddo !i
 
 return
 
end subroutine RR


!Poisson������
subroutine Poisson(N, ep0, pe, ds, rhop)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: ep0
 real*8, dimension(N), intent(in) :: pe
 real*8, dimension(N-1), intent(in) :: ds
 
 real*8, dimension(N), intent(out) :: rhop
 
 integer :: i
 real*8 :: p
 
 
 rhop = 0.d0
 
 do i = 2, N-1
  p = 2.d0/ds(i)/(ds(i)+ds(i-1))*pe(i+1)
  p = p + 2.d0/ds(i-1)/(ds(i)+ds(i-1))*pe(i-1)
  p = p - 2.d0/ds(i)/ds(i-1)*pe(i)
  rhop(i) = -ep0 * p
 enddo !i
 
 return
 
end subroutine Poisson


!�����`�F�b�N
subroutine CV(N, kind, rhov, rhop, ch, num, cvg)
 implicit none
 integer, intent(in) :: N, kind
 real*8, dimension(N), intent(in) :: rhov, rhop
 real*8, dimension(kind), intent(in) :: ch
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N), intent(out) :: cvg
 
 integer :: i, s
 real*8, dimension(N) :: nume, numi
 
 
 cvg = 0.d0
 nume = 0.d0
 numi = 0.d0
 
 do i = 1, N
  do s = 1, kind
   if(ch(s) < 0) nume(i) = nume(i) + abs(ch(s)*num(i, s))
   if(ch(s) > 0) numi(i) = numi(i) + abs(ch(s)*num(i, s))
  enddo !s
 enddo !i
 
 do i = 1, N
  cvg(i) = (rhov(i)-rhop(i))**2.d0/nume(i)/numi(i)
 enddo !i
 
 return
 
end subroutine CV


!Alfven���x
subroutine Alf(N, kind, mu0, cc, BB, ch, mass, num, VA)
 implicit none
 integer, intent(in) :: N, kind
 real*8, intent(in) :: mu0, cc
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: ch, mass
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N, 3), intent(out) :: VA
 
 integer :: i, s
 real*8 :: w
 
 
 do i = 1, N
  w = 0.d0
  do s = 1, kind
   if(ch(s) > 0) w = w + mass(s)*num(i, s)
  enddo !s
  VA(i, 1) = BB(i)/sqrt(mu0 * w)
  VA(i, 2) = VA(i, 1)/sqrt(1.d0 + (VA(i, 1)/cc)**2.d0)
  VA(i, 3) = VA(i, 2)/cc
 enddo !i
 
 return
 
end subroutine Alf


!a�ɂ��Đϕ�(���ϗ���)
subroutine TH2(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta2)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z), intent(out) :: theta2
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !���z�֐�(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !Vpara�쐬
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0 > 0.d0) then
         thetaL = thetaL + &
& ( aL(p)/abs(aL(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0) * &
& ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
& aL(p+1)/abs(aL(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p+1)**2.d0) * &
&ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
& / 2.d0 * abs(aL(p+1) - aL(p))
       endif
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0 > 0.d0) then
         thetaM = thetaM + &
& ( aM(p)/abs(aM(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0) * &
& ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
& aM(p+1)/abs(aM(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p+1)**2.d0) * &
& ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
& / 2.d0 * abs(aM(p+1) - aM(p))
       endif
      enddo !p
    endif
    
    theta2(i, s, j) = thetaL + thetaM
    
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH2


!mu�ɂ��Đϕ�(���ϗ���)
subroutine NM2(N, kind, Z, BB, mass, mu, theta2, num, aave)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z), intent(in) :: theta2
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N, kind), intent(out) :: aave
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 
 do i = 1, N
  do s = 1, kind
   
   !Vpara�쐬
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta2(i, s, j) + theta2(i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   aave(i, s) = (2.d0/mass(s))**(3.d0/2.d0)*pi*BB(i)/num(i, s)*nnn
   
  enddo !s
 enddo !i
 
 
 return
 
end subroutine NM2


!!!�R�����g�ANaN�͏o�����
!!!���_�����ϗ����A���x�͗��q�킲�Ƃɏo���K�v�����肻���A�܂�������



!a�ɂ��Đϕ�(���x)
subroutine TH3(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, aave, a3, theta3)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(N, kind), intent(in) :: aave
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z, 2), intent(out) :: theta3
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !���z�֐�(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !Tperp�쐬
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       thetaL = thetaL + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
        & / 2.d0 * abs(aL(p+1) - aL(p))
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       thetaM = thetaM + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
        & / 2.d0 * abs(aM(p+1) - aM(p))
      enddo !p
    endif
    
    theta3(i, s, j, 1) = mu(s, j)*(thetaL + thetaM)
    
    
    !Tpara�쐬
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0 > 0.d0) then
         thetaL = thetaL + &
& ( (aL(p)/abs(aL(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
& (aL(p+1)/abs(aL(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p+1)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
& / 2.d0 * abs(aL(p+1) - aL(p))
       endif
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0 > 0.d0) then
         thetaM = thetaM + &
& ( (aM(p)/abs(aM(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
& (aM(p+1)/abs(aM(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p+1)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
& / 2.d0 * abs(aM(p+1) - aM(p))
       endif
      enddo !p
    endif
    
    theta3(i, s, j, 2) = thetaL + thetaM
    
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH3


!mu�ɂ��Đϕ�(���x)
subroutine NM3(N, kind, Z, ee, BB, mass, mu, theta3, num, Tperp, Tpara)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, intent(in) :: ee
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z, 2), intent(in) :: theta3
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N, kind), intent(out) :: Tperp, Tpara
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 
 do i = 1, N
  do s = 1, kind
   
   !Tperp�쐬
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta3(i, s, j, 1) + theta3(i, s, j+1, 1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   
   Tperp(i, s) = 2.d0*sqrt(2.d0/mass(s)**3.d0)*(BB(i)**2.d0)*pi/num(i, s)*nnn /ee
   
   
   !Vpara�쐬
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta3(i, s, j, 2) + theta3(i, s, j+1, 2))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Tpara(i, s) = 4.d0*sqrt(2.d0/mass(s)**3.d0)*BB(i)*pi/num(i, s)*nnn /ee
   
  enddo !s
 enddo !i
 
 
 return
 
end subroutine NM3


subroutine ryusoku(N, kind, spece, aave, num, UU, kk, sigmab, taave)
 implicit none
 integer, intent(in) :: N, kind, spece
 real*8, dimension(N, kind), intent(in) :: aave, num, UU
 integer, dimension(kind), intent(in) :: kk, sigmab
 
 real*8, dimension(N, spece), intent(out) :: taave
 
 integer :: i, s
 real*8, dimension(N, spece) :: nnn
 
 
 taave = 0.d0
 nnn = 0.d0
 
 
 do i = 1, N
  do s = 1, kind
   nnn(i, kk(s)) = nnn(i, kk(s)) + num(i, s)
   if(aave(i, s) == aave(i, s)) then
     !M->I�����𐳂Ƃ���
     if(sigmab(s) == N) taave(i, kk(s)) = taave(i, kk(s)) - &
      & num(i, s)*aave(i, s)
     if(sigmab(s) == 1) taave(i, kk(s)) = taave(i, kk(s)) + &
      & num(i, s)*aave(i, s)
   endif
  enddo !s
 enddo !i
 
 taave = taave/nnn
 
 
 
 return
 
end subroutine ryusoku


!a�ɂ��Đϕ�(�^�E���x)
subroutine TH4(N, kind, Z, spece, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, kk, taave, a3, theta4)
 implicit none
 integer, intent(in) :: N, kind, Z, spece
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab, kk
 real*8, dimension(N, spece), intent(in) :: taave
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z, 2), intent(out) :: theta4
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !���z�֐�(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !Tperp�쐬
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       thetaL = thetaL + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
        & / 2.d0 * abs(aL(p+1) - aL(p))
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       thetaM = thetaM + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
        & / 2.d0 * abs(aM(p+1) - aM(p))
      enddo !p
    endif
    
    theta4(i, s, j, 1) = mu(s, j)*(thetaL + thetaM)
    
    
    !Tpara�쐬
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       if(sigmab(s) == 1) aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
       if(sigmab(s) == N) aL(p) = -a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0 > 0.d0) then
         thetaL = thetaL + &
& ( (aL(p)/abs(aL(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
& (aL(p+1)/abs(aL(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p+1)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
& / 2.d0 * abs(aL(p+1) - aL(p))
       endif
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       if(sigmab(s) == 1) aM(p) = -a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
       if(sigmab(s) == N) aM(p) = a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0 > 0.d0) then
         thetaM = thetaM + &
& ( (aM(p)/abs(aM(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
& (aM(p+1)/abs(aM(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p+1)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
& / 2.d0 * abs(aM(p+1) - aM(p))
       endif
      enddo !p
    endif
    
    theta4(i, s, j, 2) = thetaL + thetaM
    
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH4


!mu�ɂ��Đϕ�(���x)
subroutine NM4(N, kind, Z, spece, ee, mu0, BB, mass, mu, theta4, num, UU, kk, sigmab, tTperp, tTpara, tT, bs, beta, tTall)
 implicit none
 integer, intent(in) :: N, kind, Z, spece
 real*8, intent(in) :: ee, mu0
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z, 2), intent(in) :: theta4
 real*8, dimension(N, kind), intent(in) :: num, UU
 integer, dimension(kind), intent(in) :: kk, sigmab
 
 real*8, dimension(N, spece), intent(out) :: tTperp, tTpara, tT, bs
 real*8, dimension(N), intent(out) :: beta, tTall
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8, dimension(N, kind) :: Tperp, Tpara
 real*8, dimension(N, spece) :: nm
 real*8, dimension(N) :: nn
 
 
 do i = 1, N
  do s = 1, kind
   
   !Tperp�쐬
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta4(i, s, j, 1) + theta4(i, s, j+1, 1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   
   Tperp(i, s) = 2.d0*sqrt(2.d0/mass(s)**3.d0)*(BB(i)**2.d0)*pi/num(i, s)*nnn /ee
   
   
   !Tpara�쐬
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta4(i, s, j, 2) + theta4(i, s, j+1, 2))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Tpara(i, s) = 4.d0*sqrt(2.d0/mass(s)**3.d0)*BB(i)*pi/num(i, s)*nnn /ee
   
  enddo !s
 enddo !i
 
 
 tTperp = 0.d0
 tTpara = 0.d0
 nm = 0.d0
 
 
 do i = 1, N
  do s = 1, kind
   nm(i, kk(s)) = nm(i, kk(s)) + num(i, s)
   if(Tperp(i, s) == Tperp(i, s)) then
     tTperp(i, kk(s)) = tTperp(i, kk(s)) + num(i, s)*Tperp(i, s)
   endif
   if(Tpara(i, s) == Tpara(i, s)) then
     tTpara(i, kk(s)) = tTpara(i, kk(s)) + num(i, s)*Tpara(i, s)
   endif
  enddo !s
 enddo !i
 
 tTperp = tTperp/nm
 tTpara = tTpara/nm
 
 tT = (tTperp*2.d0+tTpara)/3.d0
 
 bs = 2.d0*mu0*nm*tT*ee
 beta = 0.d0
 do i = 1, N
  bs(i, :) = bs(i, :)/BB(i)**2.d0
  do s = 1, spece
   if(bs(i, s) == bs(i, s)) beta(i) = beta(i) + bs(i, s)
  enddo !s
 enddo !i
 
 nn = 0.d0
 do s = 1, spece
  nn = nn + nm(:, s)
 enddo !s
 
 tTall = beta*BB**2.d0/2.d0/mu0/nn /ee
 
 
 return
 
end subroutine NM4

