program mag_FA

implicit none

!SVC_IC�֘A�̕���
integer, parameter :: N = 219 !�����^�萔
real*8, dimension(N) :: pe0, dk, dR !�����^�s��
integer, dimension(N) :: gr !�����^�s��

!�S�̒ʂ��Ă̕���
real*8 :: gf
integer :: i, j !�����^�萔
real*8 :: Rj, mu0, MJ, Req, RI
real*8 :: lambda, s, nl, dl, BE, L
real*8, dimension(N) :: LMD, BB, dIo
real*8, parameter :: pi = 4.*atan(1.)

!SVC_parameter�֘A�̕���
real*8, dimension(13) :: parameter
character(len=80) :: dummy

!SVC_parameter�̒��o
open(30, file="PDS_parameter_E_L=4_NS.csv", action="read", status="old")
do i = 1, 13
 read(30, *) dummy, parameter(i)
end do
mu0 = parameter(4)
Rj = parameter(9)
MJ = parameter(7)
RI = parameter(13)
Req = parameter(11)

print *, mu0, MJ, pi, Rj

BE = mu0*MJ/(4.d0*pi*Rj**3.d0)
L = Req/Rj
print *, BE, L
print *, BE/L**3.d0*sqrt(1.d0+3.d0*sin(0.d0)**2.d0)/cos(0.d0)**6.d0

!SVC_IC�̒��o
open(10, file="SVC_IC_E_L=4_NS.csv", action="read", status="old")
read(10, *) !1�s�ǂݔ�΂�

!�s��Ƀf�[�^����
do i = 1, N
 read(10, *) gr(i), pe0(i), dk(i), dR(i)
end do


!Newton method
nl = -acos(sqrt(Rj/Req)) !grid1�ł�lambda(�n�_)

do i = 1, N
 s = dk(i)
 j = 0
 do
  lambda = nl
  call fg(lambda, nl, gf, s/Req)
  j = j + 1
  if(abs(gf) < 1.d-8) exit
  if(mod(j,1000000) == 0) print *, i, j, gf, nl
 enddo
 LMD(i) = nl
 dIo(i) = sqrt((Req*sin(nl)*cos(nl)**2.d0)**2.d0+(Req-Req*cos(nl)**3.d0)**2.d0)
 print *, i, nl
enddo

!�������x�̓��o
BE = mu0*MJ/4.d0/pi/Rj**3.d0
L = Req/Rj
do i = 1, N
 BB(i) = BE/L**3.d0*sqrt(1.d0+3.d0*sin(LMD(i))**2.d0)/cos(LMD(i))**6.d0
enddo

!�t�@�C��������
open(40, file="mag_FA_E_L=4_NS.csv")
do i = 1, N
 write(40, 57) LMD(i), BB(i), dIo(i)
enddo


!�I������
close(10)
close(30)
close(40)

57 format(1PE25.15, ',', 1x, 1PE25.15, ',', 1x, 1PE25.15)

end program mag_FA

!�֐�
real*8 function hh(x,s)
 real*8, intent(in) :: x,s
 hh = sin(x)*sqrt(1.d0+3.d0*sin(x)**2.d0)/2.d0+asinh(sqrt(3.d0)*sin(x))/2.d0/sqrt(3.d0)-s
 return
end function

real*8 function hhh(x)
 real*8, intent(in) :: x
 hhh = cos(x)*sqrt(1.d0+3.d0*sin(x)**2.d0)
 return
end function

!Newton method�̒��g
subroutine fg(x,nx,gf,s)
 implicit none
 real*8, intent(in) :: x,s
 real*8, intent(out) :: nx,gf
 real*8 :: hh,hhh
 nx = x - hh(x, s)/hhh(x)
 gf = hh(nx,s)
 return
end subroutine fg



