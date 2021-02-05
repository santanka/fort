program choukousou1

implicit none

integer, parameter :: N = 130 !1890-2019年
real*8 :: xfilt !推定値
real*8 :: pfilt = 1.d1**2.d0 !初期値の分散
real*8 :: QQ = 1.d0**2.d0 !システムノイズの分散
real*8 :: RR = 1.d0 !観測ノイズの分散
integer :: i, j
character(len=80) :: dummy1, dummy2, dummy3 !使用しない文字用
real*8, dimension(N) :: yy !観測データ(年平均)
real*8, dimension(12) :: month !観測データ(月別)
real*8 :: sum
real*8 :: xpred, ppred !予測値
real*8 :: kmat !カルマンゲイン


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
 !一期先予測
 xpred = xfilt
 ppred = pfilt + QQ
 
 !フィルタリング
 kmat = ppred/(ppred+RR)
 xfilt = xpred + kmat * (yy(i) - xpred)
 pfilt = (1.d0 - kmat) * ppred
 
 write(60, 61) 1889+i, yy(i), xfilt, xfilt-2.d0*sqrt(pfilt), xfilt+2.d0*sqrt(pfilt), QQ, RR

enddo !i
close(60)

61 format(I4, 6(',', 1PE25.15E3))

end program choukousou1


