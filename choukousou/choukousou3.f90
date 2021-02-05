program choukousou3

implicit none

integer, parameter :: N = 130 !1890-2019年
integer, parameter :: M = 10000 !アンサンブル数
real*8, dimension(M) :: xfilt !推定値
real*8 :: QQ = 1.d-1**2.d0 !システムノイズの分散
real*8 :: RR = 1.d-1 !観測ノイズの分散
integer :: i, j
character(len=80) :: dummy1, dummy2, dummy3 !使用しない文字用
real*8, dimension(N) :: yy !観測データ(年平均)
real*8, dimension(12) :: month !観測データ(月別)
real*8 :: sum
real*8, dimension(M) :: xpred !推定値の予測値
real*8, dimension(M) :: vv, ww!システムノイズ群
real*8 :: ppred !分散の予測値平均
real*8 :: pfilt, xx
real*8 :: kmat !カルマンゲイン
real*8 :: xbar


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

open(60, file="yamagata-tem-ave-year-an.csv")
do i = 1, N
 !予測
 call boxmuller(M, 0.d0, sqrt(QQ), vv)
 xpred = xfilt + vv
 
 !フィルタリング
 xbar = 0.d0
 do j = 1, M
  xbar = xbar + xpred(j)
 enddo !j
 xbar = xbar/dble(M)
 ppred = 0.d0
 do j = 1, M
  ppred = ppred + (xpred(j)-xbar)**2.d0
 enddo !j
 ppred = ppred/dble(M-1)
 kmat = ppred/(ppred+RR)
 
 call boxmuller(M, 0.d0, sqrt(RR), ww)
 xfilt = xpred + kmat*(yy(i)-xpred+ww)
 
 xx = 0.d0
 do j = 1, M
  xx = xx + xfilt(j)
 enddo !j
 xx = xx/dble(M)
 pfilt = (1.d0-kmat)*ppred
 
 write(60, 61) 1889+i, yy(i), xx, xx-2.d0*sqrt(pfilt), xx+2.d0*sqrt(pfilt), QQ, RR


enddo !i
close(60)

61 format(I4, 6(',', 1PE25.15E3))

end program choukousou3

subroutine boxmuller(M, mu, sigma, ran)
 implicit none
 integer, intent(in) :: M !乱数数
 real*8, intent(in) :: mu !平均
 real*8, intent(in) :: sigma !標準偏差
 real*8, dimension(M) :: ran !乱数
 
 integer :: i
 real*8, dimension(M) :: al, be
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 
 call random_number(al)
 call random_number(be)
 
 ran = sigma*sqrt(-log(al**2.d0))*sin(2.d0*pi*be) + mu
 
 return
 
end subroutine boxmuller

