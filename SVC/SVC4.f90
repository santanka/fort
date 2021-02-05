program SVC4

implicit none

!�S�̂�ʂ��Ă̕���
integer :: i, j !do�p
integer, parameter :: N = 79 !grid��
integer, parameter :: Z = 60 !���x���grid��
real*8, parameter :: pi = 4.d0*atan(1.d0) !�~����
character(len=80) :: dummy !�g�p���Ȃ���p

!�
real*8 :: B, e, c, m !��������

!mag_FA�֘A
real*8, dimension(N) :: LMD, BB, dIo

!SVC1_parameter�֘A
real*8, dimension(12) :: parameter
real*8 :: mu0, GG, alpha, MJ, MassJ, RJ, omJ, Req, MIo, ep0

!SVC_IC�֘A
real*8, dimension(N) :: pe0, dk, dR
real*8, dimension(N-1) :: dx

!SVC_BC�֘A
character(len=1), dimension(8) :: so
real*8, dimension(8) :: Nd, Tp, Ta, ch, mass

!sigma�֘A
real*8, dimension(8) :: sigma
integer, dimension(8) :: sigmab

!iteration�֘A
real*8, dimension(N) :: npe
integer :: itn
real*8, dimension(3, N) :: pep, rhov, rhop, cvg
real*8, dimension(3, N, 8) :: UU, num
real*8, dimension(3, N, 8, Z) :: UB
real*8, dimension(8, Z) :: mu
real*8, dimension(3, N, 8, Z, 3) :: a3
real*8, dimension(3, N, 8, Z) :: theta
real*8 :: cvn


!mag_FA�̒��o
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i) !LMD:���C�ܓx,BB:�������x,dIo:Io�Ƃ̋���
enddo
close(40)

B = BB(1)
BB = BB

!SVC1_parameter�̒��o
open(30, file="SVC1_parameter.csv", action="read", status="old")
do i = 1, 12
 read(30, *) dummy, parameter(i)
end do
close(30)

e = parameter(1) !�d�C�f��
c = parameter(2) !����
m = parameter(3) !�d�q�̎���
mu0 = parameter(4) !�^��̓�����
GG = parameter(5) !���L���͒萔
alpha = parameter(6) !alpha
MJ = parameter(7) !Jupiter�̎��C�o�Ɏq���[�����g
MassJ = parameter(8) !Jupiter�̎���
RJ = parameter(9) !Jupiter���a
omJ = parameter(10) !Jupiter���]�p���g��
Req = parameter(11) !Jupiter��req
MIo = parameter(12) !Io�̎���

ep0 = 1.d0/mu0/c/c !�^��̗U�d��

!SVC_IC�̒��o
open(10, file="SVC_IC.csv", action="read", status="old")
read(10, *) !1�s�ǂݔ�΂�
do i = 1, N
 read(10, *) dummy, pe0(i), dk(i), dR(i)
end do
close(10)

do i = 1, N-1
 dx = dk(i+1) - dk(i)
enddo

pe0 = pe0 !�����Ód�|�e���V����
dk = dk !i=1����̋���
dx = dx !grid�ԋ���

!�Â�����͂��߂遥
open(15, file="SVC_penum.csv", action="read", status="old")
do i = 1, N
 !if(i >= 45 .and. i <= 55) read(15, *) dummy
 !if(i < 45 .or. 55 < i) read(15, *) pe0(i)
 read(15, *) pe0(i)
enddo !i
close(15)


!SVC_BC�̒��o
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1�s�ǂݔ�΂�
do i = 1, 8
 read(20, *) dummy, Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do
close(20)

!Ta:Tpara/Tperp,so:�N��,ch:�d��
Nd = Nd !�����x
Tp = Tp*e !perp�������x
mass = mass !����
ch = ch*e

!sigma�쐬
sigma = Nd
do i = 1, 8
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo


!iteration �X�^�[�g

pep(1, :) = pe0
itn = 0

call AI(Z, alpha, N, Tp, BB, sigmab, mu)

do !itn
 itn = itn + 1
 
 !�|�e���V����
 call EP(N, ch, pep(1, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(1, :, :))
 
 !�Ód�|�e���V��������
 call pepm(pep, e, m, c, N, UU(1, :, :), ch)
 
 !�|�e���V��������
 call EP(N, ch, pep(2, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(2, :, :))
 call EP(N, ch, pep(3, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(3, :, :))
 
 !�G�l���M�[(�ʒu+Kperp)
 call UUBB(N, Z, UU, BB, mu, UB)
 
 !accessibility(1:min, 2:lim, 3:max)
 call aaa(N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 
 !a�ɂ��Đϕ�(theta)
 call TH(N, Z, a3, sigma, Tp, Ta, sigmab, UU, BB, mu, theta)
 
 !mu�ɂ��Đϕ�(num)
 call NN(N, Z, mu, theta, num)
 
 !�d�ז��x
 call RR(N, ch, num, rhov)
 
 !Poisson������
 call Poisson(N, pep, ep0, dx, rhop)
 
 !�����`�F�b�N
 call CV(N, ch, rhov, num, cvg, cvn)
 print *, itn, cvn
 print *, "    "
 
 !�t�@�C��������
 open(50, file = "SVC_penum.csv")
 do i = 1, N
  write(50, 72) pep(1, i), num(1, i, :), rhov(1, i), &
                & rhop(1, i), cvg(1, i)
 enddo !i
 close(50)
 
 open(60, file="SVC_pote.csv")
 do i = 1, N
  write(60, 92) pep(1, i), UU(1, i, :)
 enddo
 close(60)
 
 open(70, file="SVC_a3.csv")
 do i = 1, N
  write(70, 62) a3(1, i, 1, 1, :), alpha*Tp(1)/BB(sigmab(1))
 enddo !i
 close(70)
 
 
 !�����`�F�b�N
 if(cvn < 1.d-6) then
   print *, "finish"
   exit
 endif
 
 !Newton�@(npe)
 call Newtonphi(N, pep, cvg, npe, c, m, e)
 
 !NaN�`�F�b�N
 do i = 2, N-1
  if(isnan(npe(i)) .or. isnan(cvg(1, i))) then
    print *, "NaN����"
    stop
  endif
  if(npe(i) == npe(i)-1.d0 .or. cvg(1, i) == cvg(1, i)-1.d0) then
    print *, "infinity����"
    stop
  endif
 enddo !i
 
 
 !�X�V
 pep(1, :) = npe
 
enddo !itn


62 format(E25.15E3, 3(',', E25.15E3))
72 format(E25.15E3, 11(',', E25.15E3))
82 format(E25.15E3, 21(',', 1x, E25.15E3))
92 format(E25.15E3, 8(',', E25.15E3))


end program SVC4





!subroutine, function

!�f�M�s�ϗ�
subroutine AI(Z, alpha, N, Tp, BB, sigmab, mu)
 implicit none
 integer, intent(in) :: N, Z
 real*8, intent(in) :: alpha
 real*8, dimension(8), intent(in) :: Tp
 real*8, dimension(N), intent(in) :: BB
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(8, Z), intent(out) :: mu
 integer :: i, j
 
 do i = 1, 8
  do j = 1, Z
   mu(i, j) = 1.d-15*(alpha*Tp(i)/BB(sigmab(i))/1.d-15)**(dble(j-1)/dble(Z-1))
  enddo !j
 enddo !i
 
 return
 
end subroutine AI


!�|�e���V����
subroutine EP(N, ch, pe, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch, mass
 real*8, dimension(N), intent(in) :: pe, LMD, dIo
 real*8, intent(in) :: GG, MassJ, Req, omJ, MIo
 real*8, dimension(N, 8), intent(out) :: UU
 integer ::i, s
 
 do i = 1, N
  do s = 1, 8
   !�Ód
   UU(i, s) = ch(s)*pe(i)
   !�ؐ��d��
   UU(i, s) = UU(i, s) - GG*MassJ*mass(s)/Req/(cos(LMD(i))**2.d0)
   !�ؐ����S��
   UU(i, s) = UU(i,s) - mass(s)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0
   !�C�I�d��
   UU(i, s) = UU(i, s) - GG*MIo*mass(s)/dIo(i)
  enddo !s
 enddo !i
 
 do s = 1, 8
  UU(:, s) = UU(:, s) - UU(1, s)
 enddo !s
 
 return
 
end subroutine EP


!�Ód�|�e���V��������
subroutine pepm(pep, e, m, c, N, UU, ch)
 implicit none
 real*8, intent(in) :: e, m, c
 integer, intent(in) :: N
 real*8, dimension(N, 8), intent(in) :: UU
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N), intent(inout) :: pep
 real*8, dimension(N) :: dp
 integer :: i, s
 real*8 :: CC
 
 do i = 1, N
  if(pep(1, i) == 0.d0) dp(i) = 1.d-3
  if(pep(1, i) >1.0d4) dp(i) = pep(1, i)*1.d-10
  if(pep(1, i) <= 1.0d4) dp(i) = pep(1, i)*1.d-9
 enddo !i
 
 do i = 1, N
  pep(2, i) = pep(1, i) + dp(i)
  pep(3, i) = pep(1, i) - dp(i)
  !if(pep(1, i) /= 0.d0) then
  !  pep(2, i) = pep(1, i)*1.000001d0
  !  pep(3, i) = pep(1, i)*0.999999d0
  !  
  ! else if(pep(1, i) == 0.d0) then
  !  CC = abs(UU(i, 1)/ch(1))
  !  do s = 2, 8
  !   if(CC < abs(UU(i, s)/ch(s))) CC = abs(UU(i, s)/ch(s))
  !  enddo !s
  !  
  !  pep(2, i) = CC**1.d-2
  !  pep(3, i) = -CC**1.d-2
  !endif
 enddo !i
 
 return
 
end subroutine pepm


!�G�l���M�[(�ʒu+Kperp)
subroutine UUBB(N, Z, UU, BB, mu, UB)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, 8, Z), intent(out) :: UB
 integer :: ud, i, s, j
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     UB(ud, i, s, j) = UU(ud, i, s) + BB(i)*mu(s, j)
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine UUBB


subroutine aaa(N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z
 real*8, intent(in) :: alpha
 real*8, dimension(8), intent(in) :: Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8, Z) :: UB
 real*8, dimension(3, N, 8, Z, 3) :: a3
 integer, dimension(3, N, 8, Z) :: lima
 integer :: ud, i, s, j, k
 real*8 :: CC
 
 lima = 0
 
 !amin
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
    
    !�C�I�N��
    if(sigmab(s) == 1) then
      if(i == 1) a3(ud, i, s, j, 1) = 0.d0
      if(i == 2) then
        if(UB(ud, 2, s, j) >= UB(1, 1, s, j)) a3(ud, i, s, j, 1) = 0.d0
        if(UB(ud, 2, s, j) < UB(1, 1, s, j)) then
          a3(ud, i, s, j, 1) = sqrt(UB(1, 1, s, j) - UB(ud, i, s, j))
        endif
      endif
      if(i >= 3) then
        CC = UB(ud, i, s, j)
        do k = 1, i-1
         if(UB(1, k, s, j) > CC) CC = UB(1, k, s, j)
        enddo !k
        a3(ud, i, s, j, 1) = sqrt(CC - UB(ud, i, s, j))
      endif
    endif
    
    !�ؐ��N��
    if(sigmab(s) == N) then
      if(i == N) a3(ud, i, s, j, 1) = 0.d0
      if(i == N-1) then
        if(UB(ud, N-1, s, j) >= UB(1, N, s, j)) a3(ud, i, s, j, 1) = 0.d0
        if(UB(ud, N-1, s, j) < UB(1, N, s, j)) then
          a3(ud, N-1, s, j, 1) = sqrt(UB(1, N, s, j) - UB(ud, N-1, s, j))
        endif
      endif
      if(i <= N-2) then
        CC = UB(ud, i, s, j)
        do k = i+1, N
         if(UB(1, k, s, j) > CC) CC = UB(1, k, s, j)
        enddo !k
        a3(ud, i, s, j, 1) = sqrt(CC - UB(ud, i, s, j))
      endif
    endif
    
    if(a3(ud, i, s, j, 1) == 0.d0) a3(ud, i, s, j, 1) = 1.d-15
    
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 !alim
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     
     CC = 0.d0
     
     !�C�I�N��
     if(sigmab(s) == 1) then
       if(i == 1) a3(ud, i, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
       if(i /= 1) then
         CC = UB(1, 1, s, j) + alpha*Tp(s)*Ta(s) - UB(ud, i, s, j)
         if(CC <= 0.d0) a3(ud, i, s, j, 2) = 0.d0
         if(CC > 0.d0) a3(ud, i, s, j, 2) = sqrt(CC)
         
         if(a3(1, i-1, s, j, 2) < a3(1, i-1, s, j, 1)) a3(ud, i, s, j, 2) = 0.d0
       endif
     endif
     
     !�ؐ��N��
     if(sigmab(s) == N) then
       k = N+1-i
       if(k == N) a3(ud, k, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
       if(k /= N) then
         CC = UB(1, N, s, j) + alpha*Tp(s)*Ta(s) - UB(ud, k, s, j)
         if(CC <= 0.d0) a3(ud, k, s, j, 2) = 0.d0
         if(CC > 0.d0) a3(ud, k, s, j, 2) = sqrt(CC)
         
         if(a3(1, k+1, s, j, 2) < a3(1, k+1, s, j, 1)) a3(ud, k, s, j, 2) = 0.d0
       endif
     endif
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 !accessibility�̒���
 do ud = 1, 3
  do s = 1, 8
   do j = 1, Z
    
    !�C�I�N��
    if(sigmab(s) == 1) then
      do i = 1, N
       if(i /= 1 .and. lima(ud, i-1, s, j) == 1) lima(ud, i, s, j) = 1
       if(a3(ud, i, s, j, 2) <= a3(ud, i, s, j, 1)) lima(ud, i, s, j) = 1
      enddo !i
    endif
    
    !�ؐ��N��
    if(sigmab(s) == N) then
      do i = 1, N
       if(i /= 1 .and. lima(ud, N+2-i, s, j) == 1) lima(ud, N+1-i, s, j) = 1
       if(a3(ud, N+1-i, s, j, 2) <= a3(ud, N+1-i, s, j, 1)) lima(ud, N+1-i, s, j) = 1
      enddo !i
    endif
    
   enddo !j
  enddo !s
 enddo !ud
 
 
 !amax
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     
     CC = 0.d0
     
     !�C�I�N��
     if(sigmab(s) == 1) then
       if(i == N .or. lima(ud, i, s, j) == 1) a3(ud, i, s, j, 3) = 0.d0
       if(i+1 == N .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, N, s, j) - UB(ud, N-1, s, j)
         if(CC <= 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC > 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i+1 < N .and. lima(ud, i, s, j) == 0) then
         do k = i+1, N
          if(CC < UB(1, k, s, j)) CC = UB(1, k, s, j)
         enddo !k
         CC  = CC - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
     endif
     
     !�ؐ��N��
     if(sigmab(s) == N) then
       if(i == 1 .or. lima(ud, i, s, j) == 1) a3(ud, i, s, j, 3) = 0.d0
       if(i == 2 .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, 1, s, j) - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i > 2 .and. lima(ud, i, s, j) == 0) then
         do k = 1, i-1
          if(CC < UB(1, k, s, j)) CC = UB(1, k, s, j)
         enddo !k
         CC = CC - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
     endif
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine aaa


!���z�֐�
real*8 function ff(sigma, BB, B1, Ta, Tp, UU, U1, mu, aa)
 implicit none
 real*8, intent(in) :: sigma, BB, B1, Ta, Tp, UU, U1, mu, aa
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 ff = sigma*BB/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-(UU+BB*mu+aa**2.d0-(U1+B1*mu))/Ta/Tp)
 !ff = sigma*B1/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-aa/(U1-UU+Tp*(1.d0-BB/B1)+Tp*Ta))
 return

end function ff


!a�ɂ��Đϕ�(theta)
subroutine TH(N, Z, a3, sigma, Tp, Ta, sigmab, UU, BB, mu, theta)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(3, N, 8, Z, 3), intent(in) :: a3
 real*8, dimension(8), intent(in) :: sigma, Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, 8, Z), intent(out) :: theta
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 integer :: ud, i, s, j, p
 real*8 :: ff !function
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     thetaL = 0.d0
     thetaM = 0.d0
     
     !alim
     if(a3(ud, i, s, j, 2) > a3(ud, i, s, j, 1)) then
       do p = 1, Z/2
        aL(p) = a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 2)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
       enddo !p
       !i���N���ł͂Ȃ�
       if(sigmab(s) /= i) then
         do p = 1, Z/2-1
          thetaL = thetaL + (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aL(p)) &
           & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aL(p+1))) &
           & /2.d0*(aL(p+1)-aL(p))
         enddo !p
       endif
       !i���N��
       if(sigmab(s) == i) then
         do p = 1, Z/2-1
          thetaL = thetaL + (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aL(p)) &
           & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aL(p+1))) &
           & /2.d0*(aL(p+1)-aL(p))
         enddo !p
       endif
     endif
     
     !amax
     if(a3(ud, i, s, j, 3) > a3(ud, i, s, j, 1)) then
       do p = 1, Z/2
        aM(p) = -(a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 3)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
       enddo !p
       !i���N���ł͂Ȃ�
       if(sigmab(s) /= i) then
         do p = 1, Z/2-1
          thetaM = thetaM + (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aM(p)) &
           & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aM(p+1))) &
           & /2.d0*(aM(p)-aM(p+1))
         enddo !p
       endif
       !i���N��
       if(sigmab(s) == i) then
         do p = 1, Z/2-1
          thetaM = thetaM + (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aM(p)) &
           & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aM(p+1))) &
           & /2.d0*(aM(p)-aM(p+1))
         enddo !p
       endif
     endif
     
     theta(ud, i, s, j) = thetaL + thetaM
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine TH


!mu�ɂ��Đϕ�(num)
subroutine NN(N, Z, mu, theta, num)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, 8, Z), intent(in) :: theta
 real*8, dimension(3, N, 8), intent(out) :: num
 integer :: ud, i, s, j
 real*8 :: nnn
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    
    nnn = 0.d0
    
    do j = 1, Z-1
     nnn = nnn + (theta(ud, i, s, j) + theta(ud, i, s, j+1))/2.d0*(mu(s, j+1) - mu(s, j))
    enddo !j
    
    num(ud, i, s) = nnn
    
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine NN


!�d�ז��x
subroutine RR(N, ch, num, rhov)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N, 8), intent(in) :: num
 real*8, dimension(3, N), intent(out) :: rhov
 integer :: ud, i, s
 
 rhov = 0.d0
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    rhov(ud, i) = rhov(ud, i) + ch(s)*num(ud, i, s)
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine RR


!�|�A�\��������
subroutine Poisson(N, pep, ep0, dx, rhop)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(3, N), intent(in) ::pep
 real*8, dimension(N-1), intent(in) :: dx
 real*8, intent(in) :: ep0
 real*8, dimension(3, N), intent(out) :: rhop
 integer :: ud, i
 real*8 :: pp
 
 rhop = 0.d0
 
 do ud = 1, 3
  do i = 2, N-1, 1
   pp = 2.d0/dx(i)/(dx(i)+dx(i-1))*pep(ud, i+1)
   pp = pp + 2.d0/dx(i-1)/(dx(i)+dx(i-1))*pep(ud, i-1)
   pp = pp + 2.d0/dx(i)/dx(i-1)*pep(ud, i)
   rhop(ud, i) = -ep0*pp
  enddo !i
 enddo !ud
 
 return
 
end subroutine Poisson


!�����`�F�b�N
subroutine CV(N, ch, rhov, num, cvg, cvn)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N), intent(in) :: rhov
 real*8, dimension(3, N, 8), intent(in) :: num
 real*8, dimension(3, N), intent(out) :: cvg
 real*8, intent(out) :: cvn
 integer :: ud, i
 real*8, dimension(3, N) :: nume, numi
 
 cvg = 0.d0
 cvn = 0.d0
 
 do ud = 1, 3
  do i = 1, N
   nume(ud, i) = abs(ch(2)*num(ud, i, 2) + ch(7)*num(ud, i, 7) + ch(8)*num(ud, i, 8))
   numi(ud, i) = ch(1)*num(ud, i, 1) + ch(3)*num(ud, i, 3) + ch(4)*num(ud, i, 4) &
    & + ch(5)*num(ud, i, 5) + ch(6)*num(ud, i, 6)
  enddo !i
 enddo !ud
 
 do ud = 1, 3
  do i = 1, N
   cvg(ud, i) = rhov(ud, i)**2.d0/nume(ud, i)/numi(ud, i)
   if(ud == 1) cvn = cvn + cvg(ud, i)
  enddo !i
 enddo !ud
 
 cvn = sqrt(cvn/dble(N))
 
 return
 
end subroutine CV


!Newton�@(npe)
subroutine Newtonphi(N, pep, cvg, npe, c, m, e)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(3, N), intent(in) :: pep, cvg
 real*8, dimension(N), intent(out) :: npe
 real*8, intent(in) :: c, m, e
 integer :: i
 real*8 :: CC
 
 CC = 0.d0
 npe = pep(1, :)
 
 do i = 2, N-1
  CC = CC + ((pep(2, i)-pep(3, i))/(cvg(2, i)-cvg(3, i))*cvg(1, i))**2.d0
 enddo !i
 CC = sqrt(CC)
 
 if(CC >= 1.d0) then
   do i = 2, N-1
    npe(i) = npe(i) - (pep(2, i)-pep(3, i))/(cvg(2, i)-cvg(3, i))*cvg(1, i)/1.d2!/CC
   enddo !i
 endif
 if(CC < 1.d0) then
   do i = 2, N-1
    npe(i) = npe(i) - (pep(2, i)-pep(3, i))/(cvg(2, i)-cvg(3, i))*cvg(1, i)/1.d2!*CC
   enddo !i
 endif
 
 open(70, file="SVC_npe.csv")
 do i = 1, N
  write(70, 52) (pep(2, i)-pep(3, i))*m*c*c/e, (cvg(2, i)-cvg(3, i)), &
         & cvg(1, i), (pep(2, i)-pep(3, i))/(cvg(2, i)-cvg(3, i))*cvg(1, i)*m*c*c/e, npe(i)*m*c*c/e
 enddo !i
 close(70)
 
 return
 
 52 format(E25.15E3, 4(',', E25.15E3))
 
end subroutine Newtonphi


