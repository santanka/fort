program PDS_J

implicit none


!///�g�p�O�̒�������///
!�S�̂�ʂ��Ă̕���
!�Q�ƃt�@�C��(�ǉ��ŁA�Â�����n�߂� or min�Q�ƕ����̒���)
!�t�H�[�}�b�g


!�S�̂�ʂ��Ă̕���
integer :: i, j, s !do�p
integer, parameter :: N = 94 !grid��
integer, parameter :: Z = 100 !���x���grid��
real*8, parameter :: pi = 4.d0*atan(1.d0) !�~����
character(len=80) :: dummy !�g�p���Ȃ���p
integer :: MV = 28 !�ɏ��l(not parameter)
integer :: HATL = 61 !���i(not parameter)
integer :: LATL = 75 !��i(not parameter)�AHATL <= LATL
integer, parameter :: kind = 8 !���q�퐔
integer, parameter :: iodo = 2 !�ؐ����̐����x����
integer, parameter :: mado = 3 !�C�I���̐����x����

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
character(len=1), dimension(kind) :: so
real*8, dimension(kind) :: Nd, Tp, Ta, ch, mass

!sigma�֘A
real*8, dimension(kind) :: sigma
integer, dimension(kind) :: sigmab

!iteration�֘A
real*8, dimension(N) :: npe
integer :: itn, check
integer :: XX = 0
real*8, dimension(3, N) :: pep, rhov, rhop, cvg
real*8, dimension(3, N, kind) :: UU, num
real*8, dimension(3, N, kind, Z) :: UB
real*8, dimension(kind, Z) :: mu
real*8, dimension(3, N, kind, Z, 3) :: a3
real*8, dimension(3, N, kind, Z) :: theta
real*8, dimension(3, 2) :: sigx
real*8, dimension(kind) :: nsig
real*8 :: cvn, DDD, CC
real*8 :: DD = 1.d0 !�w�K��


!mag_FA�̒��o
open(40, file="mag_FA_J_kai.csv", action="read", status="old")
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
MJ = parameter(7) !�f���̎��C�o�Ɏq���[�����g
MassJ = parameter(8) !�f���̎���
RJ = parameter(9) !�f�����a
omJ = parameter(10) !�f�����]�p���g��
Req = parameter(11) !�f����req
MIo = parameter(12) !�q���̎���

ep0 = 1.d0/mu0/c/c !�^��̗U�d��

!SVC_IC�̒��o
open(10, file="PDS_IC_J_30kV_kai.csv", action="read", status="old")
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

!SVC_BC�̒��o
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1�s�ǂݔ�΂�
do i = 1, kind
 read(20, *) dummy, Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do
close(20)

!Ta:Tpara/Tperp,so:�N��,ch:�d��
Nd = Nd !�����x
Tp = Tp*e !perp�������x
mass = mass !����
ch = ch*e
sigma = Nd


!sigmab�쐬
do i = 1, kind
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo

!iteration �X�^�[�g

DDD = 0.d0

open(33, file = "PDS_J_30kV_min_kai.csv")
read(33, *) DDD
!------------------------------
do i = 1, N
 read(33, *) pe0(i)
enddo !i
do s = 1, kind
 read(33, *) sigma(s)
enddo !s
!read(33, *) MV
!read(33, *) HATL
!read(33, *) LATL
!------------------------------
close(33)
print *, DDD


!�Â�����͂��߂遥
open(15, file="PDS_J_penum.csv", action="read", status="old")
do i = 1, N
 read(15, *) CC
 pe0(i) = CC
enddo !i
close(15)


open(18, file="PDS_J_nsig.csv", action="read", status="old")
do s = 1, kind
 read(18, *) sigma(s)
enddo !s
close(18)

pep(1, :) = pe0
itn = 0


call AI(kind, Z, alpha, N, Tp, BB, sigmab, mu)

do !itn
 itn = itn + 1
 
 !�|�e���V����
 call EP(kind, N, ch, pep(1, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(1, :, :), itn)
 
 !�Ód�|�e���V��������
 call pepm(kind, LATL, HATL, pep, e, m, c, N, UU(1, :, :), ch, DD)
 
 !�|�e���V��������
 call EP(kind, N, ch, pep(2, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(2, :, :), 0)
 call EP(kind, N, ch, pep(3, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(3, :, :), 0)
 
 !�G�l���M�[(�ʒu+Kperp)
 call UUBB(kind, N, Z, UU, BB, mu, UB)
 
 !accessibility(1:min, 2:lim, 3:max)
 call aaa(kind, N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 
 !sigma����
 call sigud(kind, iodo, mado, sigma, sigx)
 
 !a�ɂ��Đϕ�(theta)
 call TH(kind, iodo, mado, N, Z, a3, sigma, sigx, Tp, Ta, sigmab, UU, BB, mu, theta)
 
 !mu�ɂ��Đϕ�(num)
 call NN(kind, N, Z, mu, theta, num)
 
 !�d�ז��x
 call RR(kind, N, ch, num, rhov)
 
 !Poisson������
 call Poisson(N, pep, ep0, dx, rhop)
 
 !�����`�F�b�N
 call CV(kind, N, ch, rhov, num, cvg, cvn, rhop)
 print *, itn, cvn
 
 !Alfven���x
 call Alf(kind, N, mu0, BB, mass, ch, num, c)
 
 !�t�@�C��������
 open(50, file = "PDS_J_penum.csv")
 do i = 1, N
  write(50, 72) pep(1, i), num(1, i, :), rhov(1, i), &
                & rhop(1, i), cvg(1, i)
 enddo !i
 close(50)
 
 open(60, file="PDS_J_pote.csv")
 do i = 1, N
  write(60, 92) pep(1, i), UU(1, i, :)
 enddo
 close(60)
 
 
 !�Œ�����l
 if(DDD > cvn .or. DDD == 0.d0) then
   DDD = cvn
   print *, DDD
   open(33, file = "PDS_J_30kV_min_kai.csv")
   write(33, '(1PE25.15E3)') cvn
   do i = 1, N
    write(33, '(1PE25.15E3)') pep(1, i)
   enddo !i
   do s = 1, kind
    write(33, '(1PE25.15E3)') sigma(s)
   enddo !s
   write(33, '(I3)') MV
   write(33, '(I3)') HATL
   write(33, '(I3)') LATL
   close(33)
 endif
 
 !�����`�F�b�N
 if((DD < 1.d-5 .or. cvn < 1.d-6) .and. check == 0) then
   print *, "finish"
   exit
 endif
 
 !Newton�@
 call Newtonphi(MV, LATL, HATL, N, pep, cvg, npe, c, m, e, cvn, itn, check, DD, XX)
 call Newtonsig(iodo, mado, kind, N, sigx, sigma, cvg, nsig, cvn, check, DD)
 
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
 sigma = nsig
 
 
enddo !itn


72 format(1PE25.15E3, 11(',', 1PE25.15E3)) !kind+4
92 format(1PE25.15E3, 8(',', 1PE25.15E3)) !kind+1


end program PDS_J





!subroutine, function

!�f�M�s�ϗ�
subroutine AI(kind, Z, alpha, N, Tp, BB, sigmab, mu)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp
 real*8, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(kind, Z), intent(out) :: mu
 integer :: s, j
 
 do s = 1, kind
  do j = 1, Z
   mu(s, j) = 1.d-20*(alpha*Tp(s)/BB(sigmab(s))/1.d-20)**(dble(j-1)/dble(Z-1))
  enddo !j
 enddo !s
 
 return
 
end subroutine AI


!�|�e���V����
subroutine EP(kind, N, ch, pe, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU, itn)
 implicit none
 integer, intent(in) :: kind, N, itn
 real*8, dimension(kind), intent(in) :: ch, mass
 real*8, dimension(N), intent(in) :: pe, LMD, dIo
 real*8, intent(in) :: GG, MassJ, Req, omJ, MIo
 real*8, dimension(N, kind), intent(out) :: UU
 integer ::i, s
 
 do i = 1, N
  do s = 1, kind
   !�ؐ��d��
   UU(i, s) = - GG*MassJ*mass(s)/Req/(cos(LMD(i))**2.d0)
   !�ؐ����S��
   UU(i, s) = UU(i,s) - mass(s)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0
   !�C�I�d��
   if(MIo /= 0.d0) then
     UU(i, s) = UU(i, s) - GG*MIo*mass(s)/dIo(i)
   endif
  enddo !s
 enddo !i
 
 do s = 1, kind
  UU(:, s) = UU(:, s) - UU(N, s)
 enddo !s
 
 if(itn == 1) then
   open(44, file="PDS_J_phizeropote.csv")
   do i = 1, N
    write(44, 86) UU(i, :)
   enddo !s
   close(44)
 endif
 
 do i = 1, N
  do s = 1, kind
   !�Ód
   UU(i, s) = UU(i, s)+ch(s)*pe(i)
  enddo !s
 enddo !i
 
 return
 
 86 format(1PE25.15E3, 7(',', 1PE25.15E3)) !kind
end subroutine EP


!�Ód�|�e���V��������
subroutine pepm(kind, LATL, HATL, pep, e, m, c, N, UU, ch, DD)
 implicit none
 real*8, intent(in) :: e, m, c, DD
 integer, intent(in) :: N, kind, LATL, HATL
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: ch
 real*8, dimension(3, N), intent(inout) :: pep
 real*8, dimension(N) :: dp
 integer :: i, s
 real*8 :: CC
 
 do i = 1, N
  if(HATL-1 >= i)dp(i) = 1.d-6
  if(HATL <= i .and. i <= LATL-1) dp(i) = 1.d-2
  if(LATL <= i) dp(i) = 1.d-7
 enddo !i
 
 do i = 1, N
  pep(2, i) = pep(1, i) + dp(i)
  pep(3, i) = pep(1, i) - dp(i)
 enddo !i
 
 return
 
end subroutine pepm


!�G�l���M�[(�ʒu+Kperp)
subroutine UUBB(kind, N, Z, UU, BB, mu, UB)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, dimension(3, N, kind), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(3, N, kind, Z), intent(out) :: UB
 integer :: ud, i, s, j
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
    do j = 1, Z
     UB(ud, i, s, j) = UU(ud, i, s) + BB(i)*mu(s, j)
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine UUBB


!accessability
subroutine aaa(kind, N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 implicit none
 integer, intent(in) :: kind, N, Z
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp, Ta
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(3, N, kind, Z) :: UB
 real*8, dimension(3, N, kind, Z, 3) :: a3
 integer :: ud, i, s, j, k, t
 
 
 !amin
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
    do j = 1, Z
    
    !�C�I�N��
    if(sigmab(s) == 1) then
      if(i == 1) a3(ud, i, s, j, 1) = 0.d0
      if(i /= 1) then
        t = 1
        do k = 1, i-1
         if(UB(1, k, s, j) > UB(1, t, s, j)) t = k
        enddo !k
        if(UB(ud, i, s, j) > UB(1, t, s, j)) then
          a3(ud, i, s, j, 1) = sqrt(UB(ud, i, s, j) - UB(1, sigmab(s), s, j))
         else if(UB(ud, i, s, j) <= UB(1, t, s, j)) then
          a3(ud, i, s, j, 1) = sqrt(UB(1, t, s, j) - UB(1, sigmab(s), s, j))
        endif
      endif
    endif
    
    !�ؐ��N��
    if(sigmab(s) == N) then
      if(i == N) a3(ud, i, s, j, 1) = 0.d0
      if(i /= N) then
        t = N
        do k = i+1, N
         if(UB(1, k, s, j) > UB(1, t, s, j)) t = k
        enddo !k
        if(UB(ud, i, s, j) > UB(1, t, s, j)) then
          a3(ud, i, s, j, 1) = sqrt(UB(ud, i, s, j) - UB(1, sigmab(s), s, j))
         else if(UB(ud, i, s, j) <= UB(1, t, s, j)) then
          a3(ud, i, s, j, 1) = sqrt(UB(1, t, s, j) - UB(1, sigmab(s), s, j))
        endif
      endif
    endif
    if(a3(ud, i, s, j, 1) == 0.d0) a3(ud, i, s, j, 1) = 1.d-17
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 !alim
 do s = 1, kind
  a3(:, :, s, :, 2) = sqrt(alpha*Tp(s)*Ta(s))
 enddo !s
 
 
 !amax
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
    do j = 1, Z
     
     !�C�I�N��
     if(sigmab(s) == 1) then
       if(i == N) a3(ud, i, s, j, 3) = 0.d0
       if(i /= N) then
         t = i+1
         do k = i+1, N
          if(UB(1, k, s, j) > UB(1, t, s, j)) t = k
         enddo !k
         a3(ud, i, s, j, 3)  = UB(1, t, s, j) - UB(1, sigmab(s), s, j)
         if(a3(ud, i, s, j, 3) < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(a3(ud, i, s, j, 3) >= 0.d0) a3(ud, i, s, j, 3) = sqrt(a3(ud, i, s, j, 3))
       endif
     endif
     
     !�ؐ��N��
     if(sigmab(s) == N) then
       if(i == 1) a3(ud, i, s, j, 3) = 0.d0
       if(i /= 1) then
         t = 1
         do k = 1, i-1
          if(UB(1, k, s, j) > UB(1, t, s, j)) t = k
         enddo !k
         a3(ud, i, s, j, 3)  = UB(1, t, s, j) - UB(1, sigmab(s), s, j)
         if(a3(ud, i, s, j, 3) < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(a3(ud, i, s, j, 3) >= 0.d0) a3(ud, i, s, j, 3) = sqrt(a3(ud, i, s, j, 3))
       endif
     endif
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 open(70, file="PDS_J_a3.csv")
 do i = 1, Z
  write(70, 62) a3(1, 30, i, 30, :)
 enddo !i
 close(70)
 
 return
 
62 format(1PE25.15E3, 2(',', 1PE25.15E3)) !3�Œ�
 
end subroutine aaa


!���z�֐�
real*8 function ff(sigma, B1, Ta, Tp, mu, aa)
 implicit none
 real*8, intent(in) :: sigma, B1, Ta, Tp, mu, aa
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 ff = sigma*B1/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-aa**2.d0/Ta/Tp)
 return

end function ff


!a�ɂ��Đϕ�(theta)
subroutine TH(kind, iodo, mado, N, Z, a3, sigma, sigx, Tp, Ta, sigmab, UU, BB, mu, theta)
 implicit none
 integer, intent(in) :: kind, iodo, mado, N, Z
 real*8, dimension(3, N, kind, Z, 3), intent(in) :: a3
 real*8, dimension(kind), intent(in) :: sigma, Tp, Ta
 real*8, dimension(3, 2), intent(in) :: sigx
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(3, N, kind), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(3, N, kind, Z), intent(out) :: theta
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 integer :: ud, i, s, j, p, numx
 real*8 :: ff !function
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
    do j = 1, Z
     thetaL = 0.d0
     thetaM = 0.d0
     
     if((i == N .and. s == iodo) .or. (i == 1 .and. s == mado)) then
       if(s == iodo) numx = 1
       if(s == mado) numx = 2
       
       !alim
       if(a3(ud, i, s, j, 2) > a3(ud, i, s, j, 1)) then
         do p = 1, Z/2
          aL(p) = a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 2)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
         enddo !p
         do p = 1, Z/2-1
          thetaL = thetaL + (ff(sigx(ud, numx), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aL(p)) &
           & + ff(sigx(ud, numx), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aL(p+1))) &
           & /2.d0*abs(aL(p+1)-aL(p))
         enddo !p
       endif
       !amax
       if(a3(ud, i, s, j, 3) > a3(ud, i, s, j, 1)) then
         do p = 1, Z/2
          aM(p) = -(a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 3)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
         enddo !p
         do p = 1, Z/2-1
          thetaM = thetaM + (ff(sigx(ud, numx), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aM(p)) &
           & + ff(sigx(ud, numx), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aM(p+1))) &
           & /2.d0*abs(aM(p)-aM(p+1))
         enddo !p
       endif
       theta(ud, i, s, j) = thetaL + thetaM
       
      else
       !alim
       if(a3(ud, i, s, j, 2) > a3(ud, i, s, j, 1)) then
         do p = 1, Z/2
          aL(p) = a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 2)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
         enddo !p
         do p = 1, Z/2-1
          thetaL = thetaL + &
           & (ff(sigma(s), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aL(p)) + &
           & ff(sigma(s), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aL(p+1))) &
           & /2.d0*abs(aL(p+1)-aL(p))
         enddo !p
       endif
       
       !amax
       if(a3(ud, i, s, j, 3) > a3(ud, i, s, j, 1)) then
         do p = 1, Z/2
          aM(p) = -(a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 3)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
         enddo !p
         do p = 1, Z/2-1
          thetaM = thetaM + (ff(sigma(s), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aM(p)) &
           & + ff(sigma(s), BB(sigmab(s)), Ta(s), Tp(s), mu(s, j), aM(p+1))) &
           & /2.d0*abs(aM(p)-aM(p+1))
         enddo !p
       endif
       
       theta(ud, i, s, j) = thetaL + thetaM
     endif
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine TH


!mu�ɂ��Đϕ�(num)
subroutine NN(kind, N, Z, mu, theta, num)
 implicit none
 integer, intent(in) :: kind, N, Z
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(3, N, kind, Z), intent(in) :: theta
 real*8, dimension(3, N, kind), intent(out) :: num
 integer :: ud, i, s, j
 real*8 :: nnn
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
    
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
subroutine RR(kind, N, ch, num, rhov)
 implicit none
 integer, intent(in) :: N, kind
 real*8, dimension(kind), intent(in) :: ch
 real*8, dimension(3, N, kind), intent(in) :: num
 real*8, dimension(3, N), intent(out) :: rhov
 integer :: ud, i, s
 
 rhov = 0.d0
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
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
   pp = pp - 2.d0/dx(i)/dx(i-1)*pep(ud, i)
   rhop(ud, i) = -ep0*pp
  enddo !i
 enddo !ud
 
 return
 
end subroutine Poisson


!�����`�F�b�N
subroutine CV(kind, N, ch, rhov, num, cvg, cvn, rhop)
 implicit none
 integer, intent(in) :: N, kind
 real*8, dimension(kind), intent(in) :: ch
 real*8, dimension(3, N), intent(in) :: rhov, rhop
 real*8, dimension(3, N, kind), intent(in) :: num
 real*8, dimension(3, N), intent(out) :: cvg
 real*8, intent(out) :: cvn
 integer :: ud, i, s
 real*8, dimension(3, N) :: nume, numi
 
 cvg = 0.d0
 cvn = 0.d0
 
 nume = 0.d0
 numi = 0.d0
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
    if(ch(s) < 0) nume(ud, i) = nume(ud, i) + abs(ch(s)*num(ud, i, s))
    if(ch(s) > 0) numi(ud, i) = numi(ud, i) + abs(ch(s)*num(ud, i, s))
   enddo !s
  enddo !i
 enddo !ud
 
 do ud = 1, 3
  do i = 1, N
   cvg(ud, i) = (rhov(ud, i)-rhop(ud, i))**2.d0/nume(ud, i)/numi(ud, i)
   if(ud == 1) cvn = cvn + cvg(ud, i)
   cvg(ud, i) = cvg(ud, i)**(1.d0/2.d0)
  enddo !i
 enddo !ud
 
 cvn = sqrt(cvn/dble(N))
 
 return
 
end subroutine CV


!Newton�@(npe)
subroutine Newtonphi(MV, LATL, HATL, N, pep, cvg, npe, c, m, e, cvn, itn, check, DD, XX)
 implicit none
 integer, intent(in) :: N, itn
 real*8, dimension(3, N), intent(in) :: pep, cvg
 real*8, dimension(N), intent(out) :: npe
 real*8, intent(in) :: c, m, e
 integer, intent(inout) :: MV, LATL, HATL, XX
 integer, dimension(LATL) :: kyoku
 real*8, intent(out) :: DD
 real*8 :: cvn, A
 real*8, save :: CL, DL, LL
 
 integer, intent(inout) :: check
 integer :: i, j, mxn, kyokuchi, search, che
 real*8 :: CC, CC1, CC2, CC3, mx, ser
 
 !�w�K���ݒ�
 if(itn == 1) then
   CL = cvn
   LL = 4.d1
   check = 0
   che = 0
 endif
 !if(LL > 8.d1) LL = 1.d0
 if(CL > cvn .and. itn /= 1 .and. ((cvn < 1.d0 .and. XX == 1) .or. check /= 0)) then
   CL = cvn !�ŏ���cvn���L��
   LL = LL - 1.d0
   che = 0
   check = 0
  else if(CL <= cvn .and. itn /= 1 .and. ((cvn < 1.d0 .and. XX == 1) .or. check /= 0)) then
   che = che + 1
   if(DL > cvn) LL = LL - 1.d-1
   if(DL <= cvn) then
     LL = LL + 5.d0
     check = check + 1
   endif
 endif
 
 !LL�̉���
 if(che > 50 .or. check > 5 .or. (CL*3.d0 < cvn .and. XX ==1)) then
   DL = cvn !�ЂƂO��cvn���L��
   open(33, file = "PDS_J_30kV_min_kai.csv")
    read(33, *)
    do i = 1, N
      read(33, *) npe(i)
    enddo !i
   close(33)
   LL = LL + 1.d1
   check = 6
   print *, "jump"
   goto 100
 endif
 
 DL = cvn !�ЂƂO��cvn���L��
 
 DD = 1.1**(-LL+1.d0) !�w�K��
 XX = 0 !�w�K���N�����Z�b�g
 
 CC1 = 0.d0
 CC2 = 0.d0
 npe = pep(1, :)
 
 if(cvn > 1.d0) cvn = 1.d0
 
 do i = 2, HATL-1
  if(((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. cvn < cvg(1, i) .and. abs(cvg(2, i)-cvg(3, i)) > 1.d-10) then           !
    if(abs(cvg(2, i)-cvg(3, i)) > 1.d-12) then
      if(CC1 < abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))) then
        CC1 = abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))
      endif
    endif
  endif
 enddo !i
 if(LATL > HATL) then
  do i = HATL, LATL-1
   if(((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. cvn < cvg(1, i) .and. abs(cvg(2, i)-cvg(3, i)) > 1.d-10) then                !
     if(CC2 < abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))) then
       if(abs(cvg(2, i)-cvg(3, i)) > 1.d-12) then
         CC2 = abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))
       endif
     endif
   endif
  enddo !i
 endif
 do i = LATL, N-1
  if(((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. cvn < cvg(1, i) .and. abs(cvg(2, i)-cvg(3, i)) > 1.d-10) then              !
    if(abs(cvg(2, i)-cvg(3, i)) > 1.d-12) then
      if(CC3 < abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))) then
        CC3 = abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))
      endif
    endif
  endif
 enddo !i
 
print *, CC1
print *, CC2
print *, CC3
 
 do i = 2, N-1
 if(i <= LATL-1 .or. LATL <= i) then !i <= LATL-1 .or. LATL <= i
  
  !�ω��̒���
  if((cvg(1, i) > cvg(3, i) .or. cvg(1, i) > cvg(2, i)) .and. cvn*0.25 < cvg(1, i)) then                       ! 
  if(abs(cvg(2, i)-cvg(3, i)) > 1.d-13) then
  
  CC = ((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(1, i)-pep(3, i))/(cvg(1, i)-cvg(3, i)))/2.d0*cvg(1, i)
  
  if(CC*DD < abs(npe(i))*5.d0) then
    npe(i) = npe(i) - CC*DD
  
!!CC1�̒���
!  if(i < HATL .and. CC > 1.d0 .and. cvn > 1.d0 .and. CC/CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC/CC*DD*1.d-1
!   else if(i < HATL .and. CC > 1.d0 .and. cvn <= 1.d0 .and. CC/CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC/CC*DD*1.d-1
!   else if(i < HATL .and. CC <= 1.d0 .and. cvn > 1.d-3 .and. CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC*DD
!   else if(i < HATL .and. CC <= 1.d0 .and. cvn <= 1.d-3 .and. CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC*DD
!  
!!CC2�̒���
!   else if((LATL > HATL) .and. (i >= HATL .and. i < LATL) .and. CC > 1.d0 .and. CC/CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC/CC*DD
!   else if((LATL > HATL) .and. (i >= HATL .and. i < LATL) .and. CC <= 1.d0 .and. CC*DD < abs(npe(i))*5.d0) then
!    XX = 0 !�w�K���N��
!    npe(i) = npe(i) - CC*DD
!  
!!CC3�̒���
!   else if(i >= LATL .and. CC3 > 1.d0 .and. cvn > 1.d0 .and. CC/CC*DD*1.d-2 < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC/CC*DD*1.d-3
!   else if(i >= LATL .and. CC3 > 1.d0 .and. cvn <= 1.d0 .and. CC/CC*DD*1.d-2 < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC/CC*DD*1.d-3
!   else if(i >= LATL .and. CC3 <= 1.d0 .and. cvn > 1.d-3 .and. CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC*DD
!   else if(i >= LATL .and. CC3 <= 1.d0 .and. cvn <= 1.d-3 .and. CC*DD < abs(npe(i))*5.d0) then
!    npe(i) = npe(i) - CC*DD
!   
  endif
  
!  if(abs(cvg(2, i)-cvg(3, i)) <= 1.d-13 .and. npe(i) < 29000.) npe(i) = npe(i)*0.99
  
  if(npe(i) < -3.d1 .or. npe(i) > 3.1d4) npe(i) = (pep(1, i+1)+pep(1, i-1))/2.d0
  
  endif
  
  if(abs(cvg(2, i)-cvg(3, i)) <= 1.d-12) then
    A = ((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(1, i)-pep(3, i))/(cvg(1, i)-cvg(3, i)))/2.d0
    if(isnan(npe(i) - pep(1, i)*1.d-2*A/abs(A))) then
      npe(i) = pep(1, i)
     else
      npe(i) = npe(i) - pep(1, i)*1.d-2*A/abs(A)
    endif
  endif
  
   else if((LATL > HATL) .and. (i >= HATL .and. i < LATL)) then
     XX = 0 !�w�K���N��
  endif
 
 endif !�����I�ȃC�^���[�V����
 enddo !i
 
 
 !MV�܂ŁA��ɐÓd�|�e���V�����͓d�����Ɍ������Č�������Ɖ���
 
 search = 1
 do while(search == 1)
  search = 0
  do i = 1, MV-1
   if(npe(i) < npe(i+1)) then
     if(search == 0) search = 1
     if(i == 1) then
       npe(i+1) = 2.d0*npe(i)-npe(i+1)
      else
       ser = npe(i)
       npe(i) = npe(i+1)
       npe(i+1) = ser
     endif
   endif
  enddo !i
 enddo
 
 
 !HATL�ȍ~�A��ɐÓd�|�e���V�����͓d�����Ɍ������đ�������Ɖ���
 search = 1
 do while(search == 1)
  search = 0
  do i = MV+1, N-1
   if(npe(i) > npe(i+1)) then
     if(search == 0) search = 1
     if(i == N-1) then
       npe(i) = 2.d0*npe(i+1) - npe(i)
      else
       ser = npe(i)
       npe(i) = npe(i+1)
       npe(i+1) = ser
     endif
   endif
  enddo !i
 enddo
 
 !if(npe(MV) > npe(MV+1) .or. npe(MV) > npe(MV-1)) then
 !  if(npe(MV+1) > npe(MV-1)) then
 !    ser = npe(MV-1)
 !    npe(MV-1) = npe(MV)
 !    npe(MV) = ser
 !    MV = MV-1
 !   else if(npe(MV-1) > npe(MV+1)) then
 !    ser = npe(MV+1)
 !    npe(MV+1) = npe(MV)
 !    npe(MV) =ser
 !    MV = MV+1
 !  endif
 !  print *, "MV is changed to ", MV
 !endif
 
 100 print *, "  "
 
 
 !HATL,LATL�̕␳
 !if(HATL < LATL) then
 ! if(npe(HATL) < ((npe(HATL+1)+npe(LATL-1))/2.d0+npe(HATL-1))/2.d0 .and. npe(HATL) < 3.d0*npe(HATL-1)-2.d0*npe(HATL-2)) then
 !   HATL = HATL + 1
 !   print *, "HATL is changed to ", HATL
 ! endif
 ! if(npe(HATL-1) > ((npe(HATL)+npe(LATL-1))/2.d0+npe(HATL-2))/2.d0 .and. npe(HATL-1) > 3.d0*npe(HATL-2)-2.d0*npe(HATL-3)) then
 !   HATL = HATL - 1
 !   print *, "HATL is changed to ", HATL
 ! endif
 ! if(npe(LATL) < ((npe(LATL+1)+npe(N))/2.d0+npe(LATL-1))/2.d0 .and. npe(LATL) < 4.d0*npe(LATL+1)-3.d0*npe(LATL+2)) then
 !   LATL = LATL + 1
 !   print *, "LATL is changed to ", LATL
 ! endif
 ! if(npe(LATL-1) > ((npe(LATL)+npe(N))/2.d0+npe(LATL-2))/2.d0 .and. npe(LATL-1) > 4.d0*npe(LATL)-3.d0*npe(LATL+1)) then
 !   LATL = LATL - 1
 !   print *, "LATL is changed to ", LATL
 ! endif
 !endif
 
 if(HATL == LATL) then
   if(npe(LATL) < 3.d0*npe(HATL+1)-2.d0*npe(HATL+2)) then
     HATL = HATL + 1
   endif
   if(npe(LATL-1) > 3.d0*npe(LATL-2)-2.d0*(LATL-3)) then
     LATL = LATL - 1
   endif
 endif
 
 
 open(70, file="PDS_npe.csv")
 do i = 2, N-1
  write(70, 52) (pep(2, i)-pep(3, i)), (cvg(2, i)-cvg(3, i)), &
 &((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(1, i)-pep(3, i))/(cvg(1, i)-cvg(3, i)))/2.d0*cvg(1, i), &
 &cvg(1, i), cvg(2, i), cvg(3, i), npe(i)-pep(1, i), npe(i)
 enddo !i
 close(70)
 
 
 print *, npe(HATL-1)-pep(1, HATL-1), npe(HATL)-pep(1, HATL), LL, DD, check
 print *, "  "
 return
 
 52 format(1PE25.15E3, 7(',', 1PE25.15E3)) !8�Œ�
 
end subroutine Newtonphi


!sigma�����쐬
subroutine sigud(kind, iodo, mado, sigma, sigx)
 implicit none
 integer, intent(in) :: kind, iodo, mado
 real*8, dimension(kind), intent(in) :: sigma
 real*8, dimension(3, 2), intent(out) :: sigx
 integer :: s
 
 sigx(1, 1) = sigma(iodo)
 sigx(1, 2) = sigma(mado)
 do s = 1, 2
  sigx(2, s) = sigx(1, s)*1.00000001d0
  sigx(3, s) = sigx(1, s)*0.99999999d0
 enddo !s
 
 return
 
end subroutine sigud

!Newton�@(nsig)
subroutine Newtonsig(iodo, mado, kind, N, sigx, sigma, cvg, nsig, cvn, check, DD)
 implicit none
 integer, intent(in) :: N, iodo, mado, kind
 real*8, dimension(3, 2), intent(in) :: sigx
 real*8, dimension(kind), intent(in) :: sigma
 real*8, dimension(3, N), intent(in) :: cvg
 real*8, intent(in) :: DD
 real*8, dimension(kind), intent(out) :: nsig
 integer, intent(inout) :: check
 real*8 :: cvn, CC
 integer :: s, i
 
 !jump
 if(check > 5) then
   open(33, file = "PDS_J_30kV_min_kai.csv")
    do i = 1, N+1
      read(33, *)
    enddo !i
    do s = 1, kind
      read(33, *) nsig(s)
    enddo !s
   close(33)
   check = 0
   goto 200
 endif
 
 
 nsig = sigma
 
 if(cvn > 1.d0) cvn = 1.d0
 
 if((cvg(1, N) > cvg(3, N)) .or. (cvg(1, N) > cvg(2, N))) then
   CC = ((sigx(2, 1) - sigx(1, 1))/(cvg(2, N) - cvg(1, N))+(sigx(3, 1) - sigx(1, 1))/(cvg(3, N) - cvg(1, N)))/2.d0*cvg(1, N)&
&*DD
   if(abs(CC) < sigma(iodo)*1.d-1) nsig(iodo) = sigma(iodo) - CC*1.d1
   if(abs(CC) >= sigma(iodo)*1.d-1) nsig(iodo) = sigma(iodo) - CC/abs(CC)*sigma(iodo)
 endif
 if((cvg(1, 1) > cvg(3, 1)) .or. (cvg(1, 1) > cvg(2, 1))) then
   CC = ((sigx(2, 2) - sigx(1, 2))/(cvg(2, 1) - cvg(1, 1))+(sigx(3, 2) - sigx(1, 2))/(cvg(3, 1) - cvg(1, 1)))/2.d0*cvg(1, 1)&
&*DD
   if(abs(CC) < sigma(mado)*1.d-1) nsig(mado) = sigma(mado) - CC*1.d1
   if(abs(CC) >= sigma(mado)*1.d-1) nsig(mado) = sigma(mado) - CC/abs(CC)*sigma(mado)
 endif
 
 open(44, file="PDS_J_nsig.csv")
 do s = 1, kind
  write(44, 43) nsig(s), nsig(s)-sigma(s)
 enddo !s
 close(44)
 
 200 print *, "  "
 
 return
 
 43 format(1PE25.15E3, ',', 1PE25.15E3) !2�Œ�
 
end subroutine Newtonsig


!Alfven���x
subroutine Alf(kind, N, mu0, BB, mass, ch, num, c)
 implicit none
 integer, intent(in) :: N, kind
 real*8, intent(in) :: mu0, c
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass, ch
 real*8, dimension(3, N, kind), intent(in) :: num
 real*8, dimension(N) :: Va, Vr
 integer :: i, s
 real*8 :: CC
 
 do i = 1, N
  CC =  0.d0
  do s = 1, kind
   CC = CC + mass(s)*num(1, i, s)
  enddo !s
  Va(i) = BB(i)/sqrt(mu0*CC)
  Vr(i) = Va(i)/sqrt(1.d0 + (Va(i)/c)**2.d0)
 enddo !i
 
 open(99, file="PDS_J_AW.csv")
 do i = 1, N
  write(99, 98) Va(i), Vr(i), Vr(i)/c
 enddo !i
 close(99)
 
 return
 
 98 format(1PE25.15E3, 2(',', 1PE25.15E3)) !3�Œ�
 
end subroutine Alf
