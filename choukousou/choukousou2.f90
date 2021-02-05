program choukousou2

implicit none

integer, parameter :: N = 1560 !1890-2019�N
real*8, dimension(13) :: xfilt !�����l�̐���l
real*8, dimension(13, 13) :: pfilt !�����l�̕��U
real*8, dimension(13, 13) :: QQ !�V�X�e���m�C�Y�̕��U
real*8 :: RR !�ϑ��m�C�Y�̕��U
integer :: i, j, k, l, m
character(len=80) :: dummy1, dummy2 !�g�p���Ȃ������p
real*8, dimension(N) :: yy !�ϑ��f�[�^
real*8, dimension(13) :: xpred !�\���l�̍s��
real*8, dimension(13, 13) :: ppred !�\���l�̍s��
real*8 :: HH(13) = reshape( (/1., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./), shape(HH) )
real*8, dimension(13, 13) :: FF
real*8, dimension(13) :: kmat !�J���}���Q�C��
real*8, dimension(13) :: BB
real*8, dimension(13, 13) ::  DD, EE
real*8, dimension(13, 13) :: II !�P�ʍs��
real*8 :: CC

II = 0.d0
do i = 1, 13
 II(i, i) = 1.d0
enddo !i

xfilt = 0.d0

FF = 0.d0
FF(1, 1) = 2.d0
FF(1, 2) = -1.d0
FF(2, 1) = 1.d0
do i = 3, 13
 FF(3, i) = -1.d0
 if(i /= 13) then
   FF(i+1, i) = 1.d0
 endif
enddo !i

pfilt = 0.d0
QQ = 0.d0
RR = 1.d0**2.d0
do i = 1, 13
 pfilt(i, i) = 1.d1**2.d0
 QQ(i, i) = 1.d1**2.d0
enddo !i


open(50, file="yamagata-tem-ave-month.csv", action="read", status="old")
do i = 1, N
 read(50, *) yy(i), dummy1, dummy2
enddo !i
close(50)

open(60, file="yamagata-tem-ts.csv")
do i = 1, N
 !�����\��
 xpred = matmul(FF, xfilt)
 ppred = matmul(FF, pfilt)
 ppred = matmul(ppred, transpose(FF)) + QQ
 
 !�t�B���^�����O
 kmat = matmul(ppred, HH)
 BB = matmul(HH, ppred)
 CC = dot_product(BB, HH) + RR
 
 kmat = kmat/CC
 xfilt = xpred + (yy(i)-dot_product(HH, xpred))*kmat
 
 do j = 1, 13
  do k = 1, 13
   EE(j, k) = kmat(j)*HH(k)
  enddo !k
 enddo !j
 
 pfilt = matmul(II-EE, ppred)
 
 write(60, 61) 1.889d3+dble(i-1)/1.2d1, yy(i), xfilt(1)+xfilt(3), xfilt(1), xfilt(3), QQ(1, 1), RR
 
enddo !i
close(60)

61 format(1PE25.15E3, 6(',', 1PE25.15E3))


end program choukousou2


