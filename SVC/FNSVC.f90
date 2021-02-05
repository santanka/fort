program FNSVC

implicit none !���@�@���@��

!�S�̂�ʂ��Ă̕���
integer :: i !�����^�萔
integer, parameter :: N = 79 !�����^�萔
real*8, parameter :: pi = 4.*atan(1.)
character(len=80) :: dummy

!�
real*8 :: B, e, c, m

!mag_FA�֘A�̕���
real*8, dimension(N) :: LMD, BB, dIo

!SVC1_parameter�֘A�̕���
real*8, dimension(12) :: parameter
real*8 :: mu0, GG, alpha, MJ, MassJ, RJ, omJ, Req, MIo

!SVC_IC�֘A�̕���
real*8, dimension(N) :: pe0, dk, dR !�����^�s��
integer, dimension(N) :: gr !�����^�s��

!SVC_BC�֘A�̕���
character(len=15), dimension(8) :: sp
character(len=1), dimension(8) :: so
real*8, dimension(8) :: Nd, Tp, Ta, ch, mass

!iteration�ȍ~�̕���
real*8, dimension(8) :: sigma, sigmaN
integer, dimension(8) :: sigmab
real*8, dimension(N) :: pe, pem, pep, phi, UU, peN
integer :: itn, ii, ss, jj, pp, ud, iii
real*8, dimension(60) :: mu, THT
real*8 :: amin, amax, alim, nis, rhn, GF, theta
real*8, dimension(3, N, 8) :: number
real*8, dimension(3, N) :: RHO

!function
real*8 :: ep, AI, yy

real*8 :: ff
integer, parameter :: z = 250
real*8, dimension(z) :: linsmu
real*8, dimension(2*z) :: linsa
real*8, dimension(N, 8, 2*z) :: laa
real*8, dimension(N, 8, z) :: lmu
real*8, dimension(N, 8, z, 2*z) :: FUNCTION
real*8 :: sum
real*8, dimension(N, 8) :: SUMMATION

!mag_FA�̒��o
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i)
enddo

B = BB(1)
BB = BB/B

!SVC1_parameter�̒��o
open(30, file="SVC1_parameter.csv", action="read", status="old")
do i = 1, 12
 read(30, *) dummy, parameter(i)
end do
e = parameter(1) !�d�C�f��
c = parameter(2) !����
m = parameter(3) !�d�q�̎���
mu0 = parameter(4)*e**3.d0*B/m**2.d0/c !�^��̓�����
GG = parameter(5)*e*B/c**3.d0 !���L���͒萔
alpha = parameter(6) !alpha
MJ = parameter(7)*B/m/c**2.d0 !Jupiter�̎��C�o�Ɏq���[�����g
MassJ = parameter(8)/m !Jupiter�̎���
RJ = parameter(9)*e*B/m/c !Jupiter���a
omJ = parameter(10)*m/e/B !Jupiter���]�p���g��
Req = parameter(11)*e*B/m/c !Jupiter��req
MIo = parameter(12)/m !Io�̎���

!SVC_IC�̒��o
open(10, file="SVC_IC.csv", action="read", status="old")
read(10, *) !1�s�ǂݔ�΂�
do i = 1, N
 read(10, *) gr(i), pe0(i), dk(i), dR(i)
end do

pe0 = 0.d0!pe0*e/m/c/c
dk = dk*e*B/m/c

!SVC_BC�̒��o
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1�s�ǂݔ�΂�
do i = 1, 8
 read(20, *) sp(i), Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do

Nd = Nd*(m*c/e/B)**3.d0
Tp = Tp*e/m/c/c
print *, Tp
mass = mass/m

!sigma�쐬
sigma = Nd
do i = 1, 8
 if(so(i) == "I") sigmab(i) = N
 if(so(i) == "M") sigmab(i) = 1
 print *, i, sigma(i), sigmab(i), Tp(i)
enddo

!���z�֐��f�[�^�쐬

pe = pe0

do ii = 1, N
 
 do ss = 1, 8
  do iii = 1, N
   UU(iii) = ep(ch(ss), pe(iii), GG, MassJ, mass(ss), Req, LMD(iii), omJ, dIo(iii), MIo)
  enddo
  call linspace(1.d-5, 1.d0, 2*z, linsa) !sqrt(alpha*Ta(ss)*Tp(ss))
  call linspace(1.d-10, alpha*Tp(ss)/BB(ii), z, linsmu) !
  laa(ii, ss, :) = linsa
  lmu(ii, ss, :) = linsmu
  
  
  do jj = 1, z
   do pp = 1, 2*z
    call FUNC(ss, N, ii, sigma(ss), BB, Ta(ss), Tp(ss), UU, sigmab(ss), mass(ss), ff, sum, linsmu(jj), linsa(pp))
    FUNCTION(ii, ss, jj, pp) = ff
    if(jj == 1 .and. pp == 1) SUMMATION(ii, ss) = sum
   enddo
  enddo
 enddo
! print *, ii
! print *, "   "
enddo

do i = 1, N
 !print *, i, pe(i)
 print *, i, ch(5)*pe(i)*m*c*c/e
 print *, i, -GG*MassJ*mass(5)/Req/(cos(LMD(i))**2.d0)*m*c*c/e
 print *, i, -mass(5)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0*m*c*c/e
 print *, i, -GG*MIo*mass(5)/dIo(i)*m*c*c/e
 print *, i, "U=", ep(ch(5), pe(i), GG, MassJ, mass(5), Req, LMD(i), omJ, dIo(i), MIo)*m*c*c/e
 !print *, i, "B=", BB(i)*B
 print *, "    "
 !print *, "    "
enddo

!�t�@�C��������
open(60, file="FNSVC.csv", status="unknown")
do i = 1, 2*z
 write(60, 82) FUNCTION(70, 7, :, i)
enddo
open(70, file="FNSVC_SUM.csv", status="unknown")
do i = 1, N
 write(70, 92) SUMMATION(i, :)
enddo
open(80, file="FNSVC_a.csv", status="unknown")
do i = 1, 2*z
 write(80, 72) laa(70, 7, i)*sqrt(m*c*c)
enddo
open(90, file="FNSVC_mu.csv")
do i = 1, z
 write(90, 72) lmu(70, 7, i)*(m*c*c/B)
enddo

!�I������
close(10)
close(20)
close(30)
close(40)
close(60)
close(70)
close(80)
close(90)

!�t�H�[�}�b�g
72 format(E25.15E3)
82 format(E25.15E3, 249(',', 1x, E25.15E3))
92 format(E25.15E3, 7(',', 1x, E25.15E3))

end program FNSVC



!subroutine
subroutine linspace(a, b, N, lins)
 implicit none
 real*8, intent(in) :: a, b
 integer, intent(in) :: N
 real*8, dimension(N), intent(out) :: lins
 integer :: i
 do i = 1, N
  lins(i) = a + (b-a)/dble(N-1)*dble(i-1)
 enddo
 return
end subroutine linspace

subroutine logspace(a, b, N, lins)
 implicit none
 real*8, intent(in) :: a, b
 integer, intent(in) :: N
 real*8, dimension(N), intent(out) :: lins
 integer :: i
 do i = 1, N
  lins(i) = a*(b/a)**(dble(i-1)/dble(N-1))
 enddo
 return
end subroutine logspace



!function
!�|�e���V�����쐬
real*8 function ep(qq, phi, GG, MassJ, mass, Req, lam, omJ, dIo, MIo)
 real*8, intent(in) :: qq, phi, GG, MassJ, mass, Req, lam, omJ, dIo, MIo
 ep = qq*phi - GG*MassJ*mass/Req/(cos(lam)**2.d0) - mass*(omJ**2.d0)*(Req**2.d0)*(cos(lam)**6.d0)/2.d0 - GG*MIo*mass/dIo
 return
end function ep

!���z�֐�
real*8 function fff(sigma, B, Bb, Ta, Tp, U, Ub, mu, a)
 real*8, intent(in) :: sigma, B, Bb, Ta, Tp, U, Ub
 real*8, intent(in) :: mu, a
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 fff = sigma*B/sqrt(pi*Ta*Tp)/Tp*exp(-Bb*mu/Tp)*exp(-(U+B*mu+a*a-(Ub+Bb*mu))/Ta/Tp)
 return
end function fff


subroutine FUNC(ss, N, ii, sigma, B, Ta, Tp, U, sigmab, mass, ff, sum, mu, aa)
 implicit none
 integer, intent(in) :: ss, N, ii, sigmab
 real*8, intent(in) :: sigma, Ta, Tp, mass, mu, aa
 real*8, dimension(N), intent(in) :: B, U
 real*8, intent(out) :: ff, sum
 real*8 :: fff !function
 
 sum = sigma*B(ii)*Ta/(B(sigmab)*(Ta-1.d0)+B(ii))*exp((U(sigmab)-U(ii))/Ta/Tp)
 
 ff = fff(sigma, B(ii), B(sigmab), Ta, Tp, U(ii), U(sigmab), mu, aa) / sum
 
 return
 
end subroutine FUNC
