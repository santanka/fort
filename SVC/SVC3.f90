program SVC3

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
real*8, dimension(N) :: pe, npe
integer :: itn
real*8, dimension(3, N) :: pep
real*8, dimension(3, N, N) :: rhov, cvg
real*8, dimension(3, N, 8) :: UU
real*8, dimension(3, N, N, 8) :: num
real*8, dimension(8, Z) :: mu
real*8, dimension(3, N, N, 8, Z, 3) :: a3
real*8, dimension(3, N, N, 8, Z) :: theta
real*8 :: cvn
real*8, dimension(8) :: nsig
real*8, dimension(N-2, N-2) :: Jacobian
real*8, dimension(N-2, N-1) :: JA
real*8, dimension(N-2) :: Jcvg
real*8, dimension(3, 2) :: sigx



!mag_FA�̒��o
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i)
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

ep0 = 1.d0/mu0/c/c

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

pe0 = pe0
dk = dk
dx = dx

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

Nd = Nd
Tp = Tp*e
mass = mass
ch = ch*e

!sigma�쐬
sigma = Nd
do i = 1, 8
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo


!iteration���s��
!�f�M�s�ϗ�
call AI(Z, alpha, N, Tp, BB, sigmab, mu)


!iteration�@�X�^�[�g
pe = pe0
pep(1, :) = pe
itn = 0

do !itn
 itn = itn + 1
 
 !�|�e���V����
 call EP(N, ch, pep(1, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(1, :, :))
 
 !�Ód�|�e���V��������
 call pepm(pep, e, m, c, N, UU(1, :, :), ch)
 
 !�|�e���V����
 call EP(N, ch, pep(2, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(2, :, :))
 call EP(N, ch, pep(3, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(3, :, :))
 
 !accessibility�ݒ�(1:min1,2:lim,3:max)
 call aaa(Z, N, alpha, mu, Tp, Ta, UU, BB, sigmab, a3)
 
 !sigma����
 call sigud(sigma, sigx)
 
 !a�ɂ��Đϕ�(theta)
 call TH(Z, N, a3, sigma, sigx, Tp, Ta, sigmab, mu, UU, BB, theta)
 
 !mu�ɂ��Đϕ�(num)
 call NN(Z, N, mu, theta, num)
 
 !�d�ז��x
 call RR(N, ch, num, rhov)
 
 !�����`�F�b�N(�v�Z)
 call CV(N, ch, rhov, num, cvg, cvn)
 print *, itn, cvn
 print *, "      "
 
 !Alfven���x
 call Alf(N, mu0, BB, mass, num(1, 1, :, :), c)
 
 
 
 !�t�@�C��������(�t���[�Y�΍�)
 open(50, file="SVC_penum.csv")
 do i = 1, N
  write(50, 72) pe(i), num(1, 1, i, :), rhov(1, 1, i), &
           & cvg(1, 1, i)
 enddo
 close(50)
 
 open(60, file="SVC_pote.csv")
 do i = 1, N
  write(60, 92) pe(i)*m*c*c/e, UU(1, i, :)*m*c*c
 enddo
 close(60)
 
 
 !----------������ԃ`�F�b�N----------
 if(itn == 1) then
   open(90, file="SVC_itn1.csv")
   do i = 1, N
    write(90, 82) pe(i), num(1, 1, i, :), rhov(1, 1, i), &
                & (UU(1, i, :)+BB(i)*mu(:, 1)), a3(1, 1, i, 7, 1, :), &
                & theta(1, 1, i, 7, 1)
   enddo
   close(90)
   
   !������ԃ`�F�b�N�݂̂̏ꍇ
   !exit
 
 endif
 !------------------------------------
 
 
 !�����`�F�b�N
 if(cvn < 1.d-6) then
   print *, "finish"
   exit
 endif
 
 !Newton�@
 
 !Jacobian����
 call Jacob(N, cvg, pep, Jacobian)
 
 !Gauss-Newton�@
 call GN(N, Jacobian, cvg, Jcvg)
 open(30, file="SVC_Jacob.csv")
 do i = 2, N-1
  write(30, 62) Jacobian(i-1, :)
 enddo !i
 close(30)
 
 !Gauss�̏����@
 call Gauss(N, Jcvg, Jacobian, pep(1, :), npe, JA)
 open(43, file="SVC_JA.csv")
 do i = 1, N-2
  write(43, 52) JA(i, :N-2), JA(i, N-1)
 enddo !i
 close(43)
 
 call Newtonsig(N, sigx, sigma, cvg, nsig)
 
 
 !NaN�`�F�b�N
 do i = 1, N
  if(isnan(npe(i))) then
    print *, "NaN����"
    stop
  endif
  if(npe(i) == npe(i)-1.d0) then
    print *, "infinity����"
    stop
  endif
 enddo !i

 
 
 !�X�V
 pep(1, :) = npe
 sigma = nsig
 
enddo !itn


52 format(E25.15E3, 75(',', E25.15E3))
62 format(E25.15E3, 76(',', E25.15E3))
72 format(E25.15E3, 10(',', E25.15E3))
82 format(E25.15E3, 21(',', 1x, E25.15E3))
92 format(E25.15E3, 8(',', E25.15E3))



end program SVC3










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
   mu(i, j) = 1.d-20*(alpha*Tp(i)/BB(sigmab(i))/1.d-20)**(dble(j-1)/dble(Z-1))
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


!accessibility�ݒ�
subroutine aaa(Z, N, alpha, mu, Tp, Ta, UU, BB, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z
 real*8, intent(in) :: alpha
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(8), intent(in) :: Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(3, N, N, 8, Z, 3), intent(out) :: a3
 real*8, dimension(3, N, 8, Z) :: UB
 real*8, dimension(3, N, N, 8, Z) :: UB2
 integer :: ud, w, i, s, j, k
 integer, dimension(3, N, N, 8, Z) :: lima
 real*8 :: CC
 
 lima = 0
 
 !�G�l���M�[(Bpara����)
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     UB(ud, i, s, j) = UU(ud, i, s) + BB(i)*mu(s, j)
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 do w = 1, N
  do i = 1, N
   if(i /= w) then
     do ud = 1, 3
      do s = 1, 8
       do j = 1, Z
        UB2(ud, w, i, s, j) = UB(1, i, s, j)
       enddo !j
      enddo !s
     enddo !ud
   endif
   if(i == w) then
     do ud = 1, 3
      do s = 1, 8
       do j = 1, Z
        UB2(ud, w, i, s, j) = UB(ud, i, s, j)
       enddo !j
      enddo !s
     enddo !ud
   endif
  enddo !i
 enddo !w
 
 
 !amin
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    do s = 1, 8
     do j = 1, Z
      
      if(sigmab(s) == 1) then !�C�I�N��
        if(i == 1) a3(ud, w, i, s, j, 1) = 0.d0
        if(i /= 1) then
          CC = 0.d0
          do k = 1, i
           if(UB2(ud, w, k, s, j) > CC) CC = UB2(ud, w, k, s, j)
          enddo
          a3(ud, w, i, s, j, 1) = sqrt(CC - UB2(ud, w, i, s, j))
        endif
      endif
      
      if(sigmab(s) == N) then !�ؐ��N��
        if(i == N) a3(ud, w, i, s, j, 1) = 0.d0
        if(i /= N) then
          CC = 0.d0
          do k = i, N
           if(UB2(ud, w, k, s, j) > CC) CC = UB2(ud, w, k, s, j)
          enddo
          a3(ud, w, i, s, j, 1) = sqrt(CC - UB2(ud, w, i, s, j))
        endif
      endif
      
     enddo !j
    enddo !s
   enddo !i
  enddo !w
 enddo !ud
 
 !alim
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    do s = 1, 8
     do j = 1, Z
      
      CC = 0.d0
      if(sigmab(s) == 1) then !�C�I�N��
       if(i == 1) CC = alpha*Tp(s)*Ta(s)
       if(i /= 1) CC = UB2(1, w, 1, s, j) + alpha*Tp(s)*Ta(s) - UB2(ud, w, i, s, j)
       if(CC < 0.d0) a3(ud, w, i, s, j, 2) = 0.d0
       if(CC >= 0.d0) a3(ud, w, i, s, j, 2) = sqrt(CC)
       if(i /= 1) then
         if(a3(ud, w, i-1, s, j, 2)-a3(ud, w, i-1, s, j, 1) < 0.d0) a3(ud, w, i, s, j, 2) = 0.d0
       endif
        
       else if(sigmab(s) == N) then !�ؐ��N��
       k = N+1-i
       if(k == N) CC = alpha*Tp(s)*Ta(s)
       if(k /= N) CC = UB2(1, w, N, s, j) + alpha*Tp(s)*Ta(s) - UB2(ud, w, k, s, j)
       if(CC < 0.d0) a3(ud, w, k, s, j, 2) = 0.d0
       if(CC >= 0.d0) a3(ud, w, k, s, j, 2) = sqrt(CC)
       if(k /= N) then
         if(a3(ud, w, k+1, s, j, 2)-a3(ud, w, k+1, s, j, 1) < 0.d0) a3(ud, w, k, s, j, 2) = 0.d0
       endif
      endif
      
      
     enddo !j
    enddo !s
   enddo !i
  enddo !w
 enddo !ud
 
 !accessibility�̒���
 do ud = 1, 3
  do w = 1, N
   do s = 1, 8
    do j = 1, Z
     if(sigmab(s) == 1) then !�C�I�N��
       do i = 1, N
        if(i /= 1) then
          if(lima(ud, w, i-1, s, j) == 1) lima(ud, w, i, s, j) = 1
        endif
        if(a3(ud, w, i, s, j, 2)-a3(ud, w, i, s, j, 1) <= 0.d0) then
          lima(ud, w, i, s, j) = 1
        endif
       enddo !i
     endif
     
     if(sigmab(s) == N) then !�ؐ��N��
       do i = 1, N
        k = N+1-i
        if(k /= 1) then
          if(lima(ud, w, k-1, s, j) == 1) lima(ud, w, k, s, j) = 1
        endif
        if(a3(ud, w, k, s, j, 2)-a3(ud, w, k, s, j, 1) <= 0.d0) then
          lima(ud, w, k, s, j) = 1
        endif
       enddo !i
     endif
    enddo !j
   enddo !s
  enddo !w
 enddo !ud
 
 !amax
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    do s = 1, 8
     do j = 1, Z
      CC = 0.d0
      if(sigmab(s) == 1) then !�C�I�N��
        if(i+1 == N .and. lima(ud, w, i, s, j) == 0) then
          CC = UB2(ud, w, N, s, j) - UB2(ud, w, i, s, j)
          if(CC < 0.d0) a3(ud, w, i, s, j, 3) = 0.d0
          if(CC >= 0.d0) a3(ud, w, i, s, j, 3) = sqrt(CC)
        endif
        if(i+1 < N .and. lima(ud, w, i, s, j) == 0) then
          CC = UB2(ud, w, i+1, s, j)
          do k = i+1, N
           if(CC < UB2(ud, w, k, s, j)) CC = UB2(ud, w, k, s, j)
          enddo
          CC = CC - UB2(ud, w, i, s, j)
          if(CC < 0.d0) a3(ud, w, i, s, j, 3) = 0.d0
          if(CC >= 0.d0) a3(ud, w, i, s, j, 3) = sqrt(CC)
        endif
        if(i == N .and. lima(ud, w, i, s, j) == 0) a3(ud, w, i, s, j, 3) = 0.d0
        if(lima(ud, w, i, s, j) == 1) a3(ud, w, i, s, j, 3) = 0.d0
        
       else if(sigmab(s) == N) then !�ؐ��N��
        if(i == 2 .and. lima(ud, w, i, s, j) == 0) then
          CC = UB2(ud, w, 1, s, j) - UB2(ud, w, i, s, j)
          if(CC < 0.d0) a3(ud, w, i, s, j, 3) = 0.d0
          if(CC >= 0.d0) a3(ud, w, i, s, j, 3) = sqrt(CC)
        endif
        if(i > 2 .and. lima(ud, w, i, s, j) == 0) then
          CC = UB2(ud, w, 1, s, j)
          do k = 1, i-1
           if(CC < UB2(ud, w, k, s, j)) CC = UB2(ud, w, k, s, j)
          enddo !k
          CC = CC - UB2(ud, w, i, s, j)
          if(CC < 0.d0) a3(ud, w, i, s, j, 3) = 0.d0
          if(CC >= 0.d0) a3(ud, w, i, s, j, 3) = sqrt(CC)
        endif
        if(i == 1 .and. lima(ud, w, i, s, j) == 0) a3(ud, w, i, s, j, 3) = 0.d0
        if(lima(ud, w, i, s, j) == 1) a3(ud, w, i, s, j, 3) = 0.d0
      endif
      
      
     enddo !j
    enddo !s
   enddo !i
  enddo !w
 enddo !ud
 
 return
 
end subroutine aaa


!���z�֐�
real*8 function ff(sigma, BB, B1, Ta, Tp, UU, U1, mu, aa)
 implicit none
 real*8, intent(in) :: sigma, BB, B1, Ta, Tp, UU, U1, mu, aa
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 ff = sigma*BB/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-(UU+BB*mu+aa**2.d0-(U1+B1*mu))/Ta/Tp)
 
 return

end function ff


!a�ɂ��Đϕ�(theta)
subroutine TH(Z, N, a3, sigma, sigx, Tp, Ta, sigmab, UU, BB, mu, theta)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(3, N, N, 8, Z, 3), intent(in) :: a3
 real*8, dimension(8), intent(in) :: sigma, Tp, Ta
 real*8, dimension(3, 2), intent(in) :: sigx
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(3, N, N, 8) :: UU2
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, N, 8, Z), intent(out) :: theta
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 integer :: ud, i, w, s, j, p, numx
 real*8 :: ff !function
 
 do w = 1, N
  do i = 1, N
   if(i /= w) then
     do ud = 1, 3
      do j = 1, 8
       UU2(ud, w, i, j) = UU(1, i, j)
      enddo !j
     enddo !ud
   endif
   if(i == w) then
     do ud = 1, 3
      do j = 1, 8
       UU2(ud, w, i, j) = UU(ud, i, j)
      enddo !j
     enddo !ud
   endif
  enddo !i
 enddo !w
 
 do ud = 1, 3
  do w = 1, N
  do i = 1, N
   do s = 1, 8
    do j = 1, Z       
     thetaL = 0.d0
     thetaM = 0.d0
     
     if((i == N .and. s == 2) .or. (i == 1 .and. s == 3)) then
       if(s == 2) numx = 1
       if(s == 3) numx = 2
       
       !alim
       if(a3(1, w, i, s, j, 2) > a3(1, w, i, s, j, 1)) then
         do p = 1, Z/2
          aL(p) = a3(1, w, i, s, j, 1)*(a3(1, w, i, s, j, 2)/a3(1, w, i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
         enddo !p
         do p = 1, Z/2-1
          thetaL = thetaL + &
           &(ff(sigx(ud, numx), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(1, w, i, s), UU2(1, w, i, s), mu(s, j), aL(p)) &
           & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(1, w, i, s), UU2(1, w, i, s), mu(s, j), aL(p+1))) &
           & /2.d0*abs(aL(p+1)-aL(p))
         enddo !p
       endif
       !amax
       if(a3(ud, w, i, s, j, 3) > a3(ud, w, i, s, j, 1)) then
         do p = 1, Z/2
          aM(p) = -(a3(1, w, i, s, j, 1)*(a3(1, w, i, s, j, 3)/a3(1, w, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
         enddo !p
         do p = 1, Z/2-1
          thetaM = thetaM + &
           &(ff(sigx(ud, numx), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(1, w, i, s), UU2(1, w, i, s), mu(s, j), aM(p)) &
           & + ff(sigx(ud, numx), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(1, w, i, s), UU2(1, w, i, s), mu(s, j), aM(p+1))) &
           & /2.d0*abs(aM(p)-aM(p+1))
         enddo !p
       endif
       theta(ud, w, i, s, j) = thetaL + thetaM
       
      else
       !alim
       if(a3(ud, w, i, s, j, 2) > a3(ud, w, i, s, j, 1)) then
         do p = 1, Z/2
          aL(p) = a3(ud, w, i, s, j, 1)*(a3(ud, w, i, s, j, 2)/a3(ud, w, i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
         enddo !p
         !i���N���ł͂Ȃ�
         if(sigmab(s) /= i) then
           do p = 1, Z/2-1
            thetaL = thetaL + &
             & (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(1, w, sigmab(s), s), mu(s, j), aL(p)) + &
             & ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(1, w, sigmab(s), s), mu(s, j), aL(p+1))) &
             & /2.d0*abs(aL(p+1)-aL(p))
           enddo !p
         endif
         !i���N��
         if(sigmab(s) == i) then
           do p = 1, Z/2-1
            thetaL = thetaL + &
             &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(ud, w, i, s), mu(s, j), aL(p)) + &
             &ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(ud, w, i, s), mu(s, j), aL(p+1))) &
             & /2.d0*abs(aL(p+1)-aL(p))
           enddo !p
         endif
       endif
       
       !amax
       if(a3(ud, w, i, s, j, 3) > a3(ud, w, i, s, j, 1)) then
         do p = 1, Z/2
          aM(p) = -(a3(ud, w, i, s, j, 1)*(a3(ud, w, i, s, j, 3)/a3(ud, w, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
         enddo !p
         !i���N���ł͂Ȃ�
         if(sigmab(s) /= i) then
           do p = 1, Z/2-1
            thetaM = thetaM + &
             &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(1, w, sigmab(s), s), mu(s, j), aM(p)) + &
             &ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(1, w, sigmab(s), s), mu(s, j), aM(p+1))) &
             & /2.d0*abs(aM(p)-aM(p+1))
           enddo !p
         endif
         !i���N��
         if(sigmab(s) == i) then
           do p = 1, Z/2-1
            thetaM = thetaM + &
             &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(ud, w, i, s), mu(s, j), aM(p)) &
             & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU2(ud, w, i, s), UU2(ud, w, i, s), mu(s, j), aM(p+1))) &
             & /2.d0*abs(aM(p)-aM(p+1))
           enddo !p
         endif
       endif
       
       theta(ud, w, i, s, j) = thetaL + thetaM
     endif
    enddo !j
   enddo !s
  enddo !i
  enddo !w
 enddo !ud
 
 return
 
end subroutine TH


!mu�ɂ��Đϕ�(num)
subroutine NN(Z, N, mu, theta, num)
 implicit none
 integer, intent(in) :: Z, N
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, N, 8, Z), intent(in) :: theta
 real*8, dimension(3, N, N, 8), intent(out) :: num
 integer :: ud, w, i, s, j
 real*8 :: nnn
 
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    do s = 1, 8
     
     nnn = 0.d0
     
     do j = 1, Z-1
      nnn = nnn + (theta(ud, w, i, s, j) + theta(ud, w, i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
     enddo !j
     
     num(ud, w, i, s) = nnn
     
    enddo !s
   enddo !i
  enddo !w
 enddo !ud
 
 return
 
end subroutine NN


!�d�ז��x
subroutine RR(N, ch, num, rhov)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N, N, 8), intent(in) :: num
 real*8, dimension(3, N, N), intent(out) :: rhov
 integer :: ud, w, i, s
 
 rhov = 0.d0
 
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    do s = 1, 8
     rhov(ud, w, i) = rhov(ud, w, i) + ch(s)*num(ud, w, i, s)
    enddo !s
   enddo !i
  enddo !w
 enddo !ud
 
 return
 
end subroutine RR


!�����`�F�b�N
subroutine CV(N, ch, rhov, num, cvg, cvn)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N, N), intent(in) :: rhov
 real*8, dimension(3, N, N, 8), intent(in) :: num
 real*8, dimension(3, N, N), intent(out) :: cvg
 real*8, intent(out) :: cvn
 real*8, dimension(3, N, N) :: nume, numi
 integer :: ud, w, i
 
 cvg = 0.d0
 cvn = 0.d0
 
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    nume(ud, w, i) = abs(ch(2)*num(ud, w, i, 2) + ch(7)*num(ud, w, i, 7) + ch(8)*num(ud, w, i, 8))
    numi(ud, w, i) = ch(1)*num(ud, w, i, 1) + ch(3)*num(ud, w, i, 3) + ch(4)*num(ud, w, i, 4) + &
                   & ch(5)*num(ud, w, i, 5) + ch(6)*num(ud, w, i, 6)
   enddo !i
  enddo !w
 enddo !ud
 
 do ud = 1, 3
  do w = 1, N
   do i = 1, N
    cvg(ud, w, i) = ((rhov(ud, w, i)**2.d0)/numi(ud, w, i)/numi(ud, w, i))
   enddo !i
  enddo !w
 enddo !ud
 
 do i = 2, N-1
  cvn = cvn + cvg(1, i, i)
 enddo !i
 cvn = sqrt(cvn/dble(N-2))
 
 return
 
end subroutine CV


subroutine Jacob(N, cvg, pep, Jacobian)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(3, N, N), intent(in) :: cvg
 real*8, dimension(3, N), intent(in) :: pep
 real*8, dimension(N-2, N-2), intent(out) :: Jacobian
 integer :: w, i
 
 do w = 2, N-1
  do i = 2, N-1
   Jacobian(i-1, w-1) = (cvg(2, w, i)-cvg(3, w, i))/(pep(2, i)-pep(3, i))
  enddo !i
 enddo !w
 
 
 return
 
end subroutine Jacob


subroutine GN(N, Jacobian, cvg, Jcvg)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(N-2, N-2), intent(inout) :: Jacobian
 real*8, dimension(3, N, N), intent(in) :: cvg
 real*8, dimension(N-2), intent(out) :: Jcvg
 real*8, dimension(N-2, N-2) :: JJ
 real*8 :: CC
 integer :: i, j, k
 
 do i = 1, N-2
  do j = 1, N-2
   CC = 0.d0
   do k = 1, N-2
    CC = CC + Jacobian(k, i)*Jacobian(k, j)
   enddo !k
   JJ(i, j) = CC
  enddo !j
 enddo !i
 
 Jacobian = JJ
 
 do i = 1, N-2
  CC = 0.d0
  do k = 1, N-2
   CC = CC + Jacobian(k, i)*cvg(1, 1, k+1)
  enddo !k
  Jcvg(i) = CC
 enddo !i
 
 return
 
end subroutine GN


subroutine Gauss(N, Jcvg, Jacobian, pe, npe, JA)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(N-2), intent(in) :: Jcvg
 real*8, dimension(N-2, N-2), intent(in) :: Jacobian
 real*8, dimension(N), intent(in) :: pe
 real*8, dimension(N), intent(out) :: npe
 real*8, dimension(N-2, N-1), intent(out) :: JA
 real*8 :: CC
 integer :: i, j, k, ii
 
 do i = 1, N-2
  do j = 1, N-1
   if(j /= N-1) JA(i, j) = Jacobian(i, j)
   if(j == N-1) JA(i, j) = -Jcvg(i)
  enddo !j
 enddo !i
 
 !�O�i����
 do i = 1, N-2
  JA(i, :) = JA(i, :)/JA(i, i)
  if(i < N-3) then
    do k = i+1, N-2
     JA(k, :) = JA(k, :) - JA(k, i)*JA(i, :)
    enddo !k
   else if(i == N-3) then
    JA(i+1, :) = JA(i+1, :) - JA(i+1, i)*JA(i, :)
  endif
 enddo !i
 
 !��ޑ��
 do i = 1, N-4
  do k = 1, N-2-i
   JA(k, :) = JA(k, :) - JA(k, N-1-i)*JA(N-1-i, :)
  enddo !k
 enddo !i
 JA(1, :) = JA(1, :) - JA(1, 2)*JA(2, :)
 
 !npe�쐬
 npe = pe
 
 CC = 0.d0
 do i = 1, N-2
  if(CC < abs(JA(i, N-1))) CC = abs(JA(i, N-1))
 enddo !i
 
 do i = 1, N-2
  if(CC > 1.d0) then
    if(i+1 <= 55) npe(i+1) = npe(i+1) + JA(i, N-1)/CC/1.d3
    if(i+1 > 55) npe(i+1) = npe(i+1) + JA(i, N-1)/CC
   else
    if(i+1 <= 55) npe(i+1) = npe(i+1) + JA(i, N-1)/1.d3
    if(i+1 > 55) npe(i+1) = npe(i+1) + JA(i, N-1)
  endif
 enddo !i
 
 return
 
end subroutine Gauss


!sigma�����쐬
subroutine sigud(sigma, sigx)
 implicit none
 real*8, dimension(8), intent(in) :: sigma
 real*8, dimension(3, 2), intent(out) :: sigx
 integer :: s
 
 sigx(1, 1) = sigma(2)
 sigx(1, 2) = sigma(3)
 do s = 1, 2
  sigx(2, s) = sigx(1, s)*1.00001d0
  sigx(3, s) = sigx(1, s)*0.99999d0
 enddo !s
 
 return
 
end subroutine sigud

!Newton�@(nsig)
subroutine Newtonsig(N, sigx, sigma, cvg, nsig)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(3, 2), intent(in) :: sigx
 real*8, dimension(8), intent(in) :: sigma
 real*8, dimension(3, N), intent(in) :: cvg
 real*8, dimension(8), intent(out) :: nsig
 integer :: s
 
 nsig = sigma
 
 nsig(2) = sigma(2) - (sigx(2, 1) - sigx(3, 1))/(cvg(2, N) - cvg(3, N))*cvg(1, N)
 nsig(3) = sigma(3) - (sigx(2, 2) - sigx(3, 2))/(cvg(2, 1) - cvg(3, 1))*cvg(1, 1)
 
 open(44, file="SVC_nsig.csv")
 do s = 1, 8
  write(44, 43) nsig(s), nsig(s)-sigma(s)
 enddo !s
 close(44)
 
 return
 
 43 format(E25.15E3, ',', E25.15E3)
 
end subroutine Newtonsig


!Alfven���x
subroutine Alf(N, mu0, BB, mass, num, c)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: mu0, c
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8), intent(in) :: mass
 real*8, dimension(3, N, 8), intent(in) :: num
 real*8, dimension(N) :: Va, Vr
 integer :: i, s
 real*8 :: CC
 
 do i = 1, N
  CC =  0.d0
  do s = 1, 8
   CC = CC + mass(s)*num(1, i, s)
  enddo !s
  Va(i) = BB(i)/sqrt(mu0*CC)
  Vr(i) = Va(i)/sqrt(1.d0 + (Va(i)/c)**2.d0)
 enddo !i
 
 open(99, file="SVC_AW.csv")
 do i = 1, N
  write(99, 98) Va(i), Vr(i), Vr(i)/c
 enddo !i
 close(99)
 
 return
 
 98 format(E25.15E3, 2(',', E25.15E3))
 
end subroutine Alf

