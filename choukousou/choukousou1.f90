program choukousou1

implicit none

integer, parameter :: N = 130 !1890-2019�N
real*8 :: xfilt !����l
real*8 :: pfilt = 1.d1**2.d0 !�����l�̕��U
real*8 :: QQ = 1.d0**2.d0 !�V�X�e���m�C�Y�̕��U
real*8 :: RR = 1.d0 !�ϑ��m�C�Y�̕��U
integer :: i, j
character(len=80) :: dummy1, dummy2, dummy3 !�g�p���Ȃ������p
real*8, dimension(N) :: yy !�ϑ��f�[�^(�N����)
real*8, dimension(12) :: month !�ϑ��f�[�^(����)
real*8 :: sum
real*8 :: xpred, ppred !�\���l
real*8 :: kmat !�J���}���Q�C��


open(50, file="yamagata-tem-ave-month.csv", action="read", status="old")
do i = 1, N
 sum = 0.d0
 do j = 1, 12
  read(50, *) month(j), dummy2, dummy3
  sum = sum + month(j)
 enddo !j
 yy(i) = sum/1.2d1
enddo !i
close(50)

xfilt = yy(1)

open(60, file="yamagata-tem-ave-year.csv")
do i = 1, N
 !�����\��
 xpred = xfilt
 ppred = pfilt + QQ
 
 !�t�B���^�����O
 kmat = ppred/(ppred+RR)
 xfilt = xpred + kmat * (yy(i) - xpred)
 pfilt = (1.d0 - kmat) * ppred
 
 write(60, 61) 1889+i, yy(i), xfilt, xfilt-2.d0*sqrt(pfilt), xfilt+2.d0*sqrt(pfilt), QQ, RR

enddo !i
close(60)

61 format(I4, 6(',', 1PE25.15E3))

end program choukousou1


