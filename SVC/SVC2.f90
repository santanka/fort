program SVC2

implicit none

!全体を通しての文字
integer :: i, j !do用
integer, parameter :: N = 79 !grid数
real*8, parameter :: pi = 4.d0*atan(1.d0) !円周率
character(len=80) :: dummy !使用しない列用

!基準
real*8 :: B, e, c, m !無次元化

!mag_FA関連
real*8, dimension(N) :: LMD, BB, dIo

!SVC1_parameter関連
real*8, dimension(12) :: parameter
real*8 :: mu0, GG, alpha, MJ, MassJ, RJ, omJ, Req, MIo, ep0

!SVC_IC関連
real*8, dimension(N) :: pe0, dk, dR
real*8, dimension(N-1) :: dx

!SVC_BC関連
character(len=1), dimension(8) :: so
real*8, dimension(8) :: Nd, Tp, Ta, ch, mass

!sigma関連
real*8, dimension(8) :: sigma
integer, dimension(8) :: sigmab

!iteration関連
real*8, dimension(N) :: pe, npe
integer :: itn
real*8, dimension(3, N) :: pep, rhov, rhop, cvg
real*8, dimension(3, N, 8) :: UU, num
real*8, dimension(8, 120) :: mu
real*8, dimension(3, N, 8, 120, 3) :: a3
real*8, dimension(3, N, 8, 120) :: theta
real*8 :: cvn
real*8, dimension(8) :: nsig


!mag_FAの抽出
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i) !LMD:磁気緯度,BB:磁束密度,dIo:Ioとの距離
enddo
close(40)

B = BB(1)
BB = BB/B

!SVC1_parameterの抽出
open(30, file="SVC1_parameter.csv", action="read", status="old")
do i = 1, 12
 read(30, *) dummy, parameter(i)
end do
close(30)

e = parameter(1) !電気素量
c = parameter(2) !光速
m = parameter(3) !電子の質量
mu0 = parameter(4)*e**3.d0*B/m**2.d0/c !真空の透磁率
GG = parameter(5)*e*B/c**3.d0 !万有引力定数
alpha = parameter(6) !alpha
MJ = parameter(7)*B/m/c**2.d0 !Jupiterの磁気双極子モーメント
MassJ = parameter(8)/m !Jupiterの質量
RJ = parameter(9)*e*B/m/c !Jupiter半径
omJ = parameter(10)*m/e/B !Jupiter自転角周波数
Req = parameter(11)*e*B/m/c !Jupiterのreq
MIo = parameter(12)/m !Ioの質量

ep0 = 1.d0/mu0/c/c !真空の誘電率

!SVC_ICの抽出
open(10, file="SVC_IC.csv", action="read", status="old")
read(10, *) !1行読み飛ばし
do i = 1, N
 read(10, *) dummy, pe0(i), dk(i), dR(i)
end do
close(10)

do i = 1, N-1
 dx = dk(i+1) - dk(i)
enddo

pe0 = pe0*e/m/c**2.d0 !初期静電ポテンシャル
dk = dk*e*B/m/c !i=1からの距離
dx = dx*e*B/m/c !grid間距離

!SVC_BCの抽出
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1行読み飛ばし
do i = 1, 8
 read(20, *) dummy, Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do
close(20)

!Ta:Tpara/Tperp,so:起源,ch:電荷
Nd = Nd*(m*c/e/B)**3.d0 !数密度
Tp = Tp*e/m/c/c !perp方向温度
mass = mass/m !質量

!sigma作成
sigma = Nd
do i = 1, 8
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo


!iteration　スタート
pe = pe0
pep(1, :) = pe
itn = 0

do !itn
 itn = itn + 1
 
 !ポテンシャル
 call EP(N, ch, pe, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(1, :, :))
 
 !静電ポテンシャル差分
 call pepm(pe, e, m, c, N, UU(1, :, :), ch, pep)
 
 !ポテンシャル
 call EP(N, ch, pep(2, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(2, :, :))
 call EP(N, ch, pep(3, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(3, :, :))
 
 !断熱不変量
 call AI(alpha, N, Tp, BB, sigmab, mu)
 
 !accessibility設定(1:min1,2:lim,3:max)
 call aaa(N, alpha, mu, Tp, Ta, UU, BB, sigmab, a3)
 
 !aについて積分(theta)
 call TH(N, a3, sigma, Tp, Ta, sigmab, mu, UU, BB, theta)
 
 !muについて積分(num)
 call NN(N, mu, theta, num)
 
 !電荷密度
 call RR(N, ch, num, rhov)
 
 !Poisson方程式
 call Poisson(N, pep, ep0, dx, rhop)
 
 
 !----------初期状態チェック----------
 if(itn == 1) then
   open(90, file="SVC_itn1.csv")
   do i = 1, N
    write(90, 82) pe(i)*m*c*c/e, num(1, i, :)*(e*B/m/c)**3.d0, rhov(1, i)*e*(e*B/m/c)**3.d0, &
                & (UU(1, i, :)+BB(i)*mu(:, 1))*m*c*c/e, a3(1, i, 1, 1, :)*sqrt(m)*c, &
                & theta(1, i, 1, 1)*(e*B/m/c)**3.d0*B/m/c/c
   enddo
   close(90)
   
   !初期状態チェックのみの場合
   !exit
 
 endif
 !------------------------------------
 
 
 !収束チェック(計算)
 call CV(N, ch, rhov, rhop, num, cvg, cvn)
 print *, itn, cvn
 print *, "      "
 
 
 !ファイル化処理(フリーズ対策)
 open(50, file="SVC_penum.csv")
 do i = 1, N
  write(50, 72) pe(i)*m*c*c/e, num(1, i, :)*(e*B/m/c)**3.d0, rhov(1, i)*e*(e*B/m/c)**3.d0, &
           & rhop(1, i)*e*(e*B/m/c)**3.d0, cvg(1, i)
 enddo
 close(50)
 
 open(60, file="SVC_pote.csv")
 do i = 1, N
  write(60, 92) pe(i)*m*c*c/e, UU(1, i, :)*m*c*c
 enddo
 close(60)
 
 
 !収束チェック
 if(cvn < 1.d-6) then
   print *, "finish"
   exit
 endif
 
 !Newton法(npe)
 call Newtonphi(N, pep, cvg, npe)
 
 !Newton法(nsig)
 call Newtonsig(N, ch, sigma, sigmab, num, cvg, nsig)
 
 
 !NaNチェック
 do i = 1, N
  if(isnan(npe(i)) .or. isnan(cvn)) then
    print *, "NaN発生"
    stop
  endif
  if(npe(i) == npe(i)-1.d0 .or. cvn == cvn-1.d0) then
    print *, "infinity 発生"
    stop
  endif
 enddo
 
 
 !更新
 pe = npe
 !sigma = nsig
 
enddo !itn


72 format(E25.15E3, 11(',', E25.15E3))
82 format(E25.15E3, 21(',', 1x, E25.15E3))
92 format(E25.15E3, 8(',', E25.15E3))


end program SVC2




!subroutine, function

!静電ポテンシャル差分
subroutine pepm(pe, e, m, c, N, UU, ch, pep)
 implicit none
 real*8, parameter :: dp = 1.d-1
 real*8, intent(in) :: e, m, c
 integer, intent(in) :: N
 real*8, dimension(N), intent(in) :: pe
 real*8, dimension(N, 8), intent(in) :: UU
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N), intent(out) :: pep
 integer :: i, s
 real*8 :: CC
 
 do i = 1, N
  pep(1, i) = pe(i)
  CC = abs(UU(i, 1)/ch(1)-pe(i))
  do s = 2, 8
   if(CC < abs(UU(i, s)/ch(s)-pe(i))) CC = abs(UU(i, s)/ch(s)-pe(i))
  enddo
  if(pe(i) == 0.d0) then
    pep(2, i) = pe(i) + 1.d0*CC*1.d-1
    pep(3, i) = pe(i) - 1.d0*CC*1.d-1
   else
    pep(2, i) = pe(i)*1.01d0
    pep(3, i) = pe(i)*0.99d0
  endif
 enddo !i
 return
end subroutine pepm


!ポテンシャル
subroutine EP(N, ch, pep, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch, mass
 real*8, dimension(N), intent(in) :: LMD, dIo
 real*8, dimension(N), intent(in) :: pep
 real*8, intent(in) :: GG, MassJ, Req, omJ, MIo
 real*8, dimension(N, 8), intent(out) :: UU
 integer :: j, k
 
 do j = 1, N
  do k = 1, 8
   !静電
   UU(j, k) = ch(k)*pep(j)
   !木星重力
   UU(j, k) = UU(j, k) - GG*MassJ*mass(k)/Req/(cos(LMD(j))**2.d0)
   !木星遠心力
   UU(j, k) = UU(j, k) - mass(k)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(j))**6.d0)/2.d0
   !イオ重力
   UU(j, k) = UU(j, k) - GG*MIo*mass(k)/dIo(j)
  enddo !k
 enddo !j
 
 do k = 1, 8
  UU(:, k) = UU(:, k) - UU(1, k)
 enddo !k
 
 return
 
end subroutine EP


!断熱不変量
subroutine AI(alpha, N, Tp, BB, sigmab, mu)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: alpha
 real*8, dimension(8), intent(in) :: Tp
 real*8, dimension(N), intent(in) :: BB
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(8, 120), intent(out) :: mu
 integer :: i, j
 
 do i = 1, 8
  do j = 1, 120
   mu(i, j) = alpha*Tp(i)/BB(sigmab(i))*dble(j-1)/dble(120-1)
  enddo !j
 enddo !i
 
 return
 
end subroutine AI


!accessibility設定
subroutine aaa(N, alpha, mu, Tp, Ta, UU, BB, sigmab, a3)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: alpha
 real*8, dimension(8, 120), intent(in) :: mu
 real*8, dimension(8), intent(in) :: Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(3, N, 8, 120, 3), intent(out) :: a3
 real*8, dimension(3, N, 8, 120) :: UB
 integer :: ud, i, s, j, k
 integer, dimension(3, N, 8, 120) :: lima
 real*8 :: CC
 
 lima = 0
 
 !エネルギー(Bpara除く)
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, 120
     UB(ud, i, s, j) = UU(ud, i, s) + BB(i)*mu(s, j)
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 !amin
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, 120
     
     if(sigmab(s) == 1) then !イオ起源
       if(i == 1) a3(ud, i, s, j, 1) = 0.d0
       if(i == 2) then
         if(UB(ud, 2, s, j) >= UB(1, 1, s, j)) a3(ud, i, s, j, 1) = 0.d0
         if(UB(ud, 2, s, j) < UB(1, 1, s, j)) then
           a3(ud, i, s, j, 1) = sqrt(UB(1, 1, s, j)-UB(ud, i, s, j))
         endif
       endif
       if(i >= 3) then
         CC = UB(ud, i, s, j)
         do k = 1, i-1
          if(UB(1, k, s, j) > CC) CC = UB(1, k, s, j)
         enddo
         a3(ud, i, s, j, 1) = sqrt(CC - UB(ud, i, s, j))
       endif
     endif
     
     if(sigmab(s) == N) then !木星起源
       if(i == N) a3(ud, i, s, j, 1) = 0.d0
       if(i == N-1) then
         if(UB(ud, N-1, s, j) >= UB(1, N, s, j)) a3(ud, i, s, j, 1) = 0.d0
         if(UB(ud, N-1, s, j) < UB(1, N, s, j)) then
           a3(ud, i, s, j, 1) = sqrt(UB(1, N, s, j)-UB(ud, N-1, s, j))
         endif
       endif
       if(i <= N-2) then
         CC = UB(ud, i, s, j)
         do k = i+1, N
          if(UB(1, k, s, j) > CC) CC = UB(1, k, s, j)
         enddo
         a3(ud, i, s, j, 1) = sqrt(CC - UB(ud, i, s, j))
       endif
     endif
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud

!alim
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, 120
     
     CC = 0.d0
     if(sigmab(s) == 1) then !イオ起源
       if(i == 1) CC = alpha*Tp(s)*Ta(s)
       if(i /= 1) CC = UB(1, 1, s, j) + alpha*Tp(s)*Ta(s) - UB(ud, i, s, j)
       if(CC < 0.d0) a3(ud, i, s, j, 2) = 0.d0
       if(CC >= 0.d0) a3(ud, i, s, j, 2) = sqrt(CC)
       if(i /= 1) then
         if(a3(ud, i-1, s, j, 2)-a3(ud, i-1, s, j, 1) < 0.d0) a3(ud, i, s, j, 2) = 0.d0
       endif
       
      else if(sigmab(s) == N) then !木星起源
       k = N+1-i
       if(k == N) CC = alpha*Tp(s)*Ta(s)
       if(k /= N) CC = UB(1, N, s, j) + alpha*Tp(s)*Ta(s) - UB(ud, k, s, j)
       if(CC < 0.d0) a3(ud, k, s, j, 2) = 0.d0
       if(CC >= 0.d0) a3(ud, k, s, j, 2) = sqrt(CC)
       if(k /= N) then
         if(a3(ud, k+1, s, j, 2)-a3(ud, k+1, s, j, 1) < 0.d0) a3(ud, k, s, j, 2) = 0.d0
       endif
     endif
     
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 !accessibilityの調整
 do ud = 1, 3
  do s = 1, 8
   do j = 1, 120
    if(sigmab(s) == 1) then !イオ起源
      do i = 1, N
       if(i /= 1) then
         if(lima(ud, i-1, s, j) == 1) lima(ud, i, s, j) = 1
       endif
       if(a3(ud, i, s, j, 2)-a3(ud, i, s, j, 1) <= 0.d0) then
         lima(ud, i, s, j) = 1
       endif
      enddo !i
    endif
    
    if(sigmab(s) == N) then !木星起源
      do i = 1, N
       k = N+1-i
       if(k /= 1) then
         if(lima(ud, k-1, s, j) == 1) lima(ud, k, s, j) = 1
       endif
       if(a3(ud, k, s, j, 2)-a3(ud, k, s, j, 1) <= 0.d0) then
         lima(ud, k, s, j) = 1
       endif
      enddo !i
    endif
   enddo !j
  enddo !s
 enddo !ud
 
!amax
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, 120
     CC = 0.d0
     if(sigmab(s) == 1) then !イオ起源
       if(i+1 == N .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, N, s, j) - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i+1 < N .and. lima(ud, i, s, j) == 0) then
         CC = 0.d0
         do k = i+1, N
          if(CC < UB(1, k, s, j)) CC = UB(1, k, s, j)
         enddo
         CC = CC - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i == N .and. lima(ud, i, s, j) == 0) a3(ud, i, s, j, 3) = 0.d0
       if(lima(ud, i, s, j) == 1) a3(ud, i, s, j, 3) = 0.d0

      else if(sigmab(s) == N) then !木星起源
       if(i == 2 .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, 1, s, j) - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i > 2 .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, 1, s, j)
         do k = 1, i-1
          if(CC < UB(1, k, s, j)) CC = UB(1, k, s, j)
         enddo
         CC = CC - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i == 1 .and. lima(ud, i, s, j) == 0) a3(ud, i, s, j, 3) = 0.d0
       if(lima(ud, i, s, j) == 1) a3(ud, i, s, j, 3) = 0.d0
     endif
     
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine aaa


!分布関数
real*8 function ff(sigma, BB, B1, Ta, Tp, UU, U1, mu, aa)
 implicit none
 real*8, intent(in) :: sigma, BB, B1, Ta, Tp, UU, U1, mu, aa
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 ff = sigma*BB/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-(UU+BB*mu+aa**2.d0-(U1+B1*mu))/Ta/Tp)
 !ff = sigma*B1/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-aa/(U1-UU+Tp*(1.d0-BB/B1)+Tp*Ta))
 return

end function ff


!aについて積分(theta)
subroutine TH(N, a3, sigma, Tp, Ta, sigmab, mu, UU, BB, theta)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(3, N, 8, 120, 3), intent(in) :: a3
 real*8, dimension(8), intent(in) :: sigma, Tp, Ta
 integer, dimension(8) :: sigmab
 real*8, dimension(8, 120), intent(in) :: mu
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(3, N, 8, 120), intent(out) :: theta
 real*8 :: thetaL, thetaM
 real*8, dimension(60) :: aL, aM
 integer :: ud, i, s, j, p
 real*8 :: ff !function
 
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, 120
     thetaL = 0.d0
     thetaM = 0.d0
     
     !alim
     if(a3(ud, i, s, j, 2) > a3(ud, i, s, j, 1)) then
       do p = 1, 60
        aL(p) = a3(ud, i, s, j, 1) + (a3(ud, i, s, j, 2)-a3(ud, i, s, j, 1))*dble(p-1)/dble(60-1)
       enddo
       if(sigmab(s) /= i) then
         do p = 1, 60-1
          thetaL = thetaL + &
           &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aL(p)) + &
           & ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aL(p+1))) &
           & /2.d0 *abs(aL(p+1)-aL(p))
         enddo
       endif
       if(sigmab(s) == i) then
         do p = 1, 60-1
          thetaL = thetaL + &
           &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, sigmab(s), s), mu(s, j), aL(p)) + &
           & ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, sigmab(s), s), mu(s, j), aL(p+1))) &
           & /2.d0 *abs(aL(p+1)-aL(p))
         enddo
       endif
     endif
     
     !amax
     if(a3(ud, i, s, j, 3) > a3(ud, i, s, j, 1)) then
       do p = 1, 60
        aM(p) = -(a3(ud, i, s, j, 1) + (a3(ud, i, s, j, 3)-a3(ud, i, s, j, 1))*dble(p-1)/dble(60-1))
       enddo
       if(sigmab(s) /= i) then
         do p = 1, 60-1
          thetaM = thetaM + &
           &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aM(p)) + &
           & ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aM(p+1))) &
           & /2.d0 *abs(aM(p)-aM(p+1))
         enddo
       endif
       if(sigmab(s) == i) then
         do p = 1, 60-1
          thetaM = thetaM + &
           &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, sigmab(s), s), mu(s, j), aM(p)) + &
           & ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, sigmab(s), s), mu(s, j), aM(p+1))) &
           & /2.d0 *abs(aM(p)-aM(p+1))
         enddo
       endif
     endif
     
     theta(ud, i, s, j) = thetaL + thetaM
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine TH


!muについて積分(num)
subroutine NN(N, mu, theta, num)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8, 120), intent(in) :: mu
 real*8, dimension(3, N, 8, 120), intent(in) :: theta
 real*8, dimension(3, N, 8), intent(out) :: num
 integer :: ud, i, s, j
 real*8 :: nnn
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    
    nnn = 0.d0
    
    do j = 1, 120-1
     nnn = nnn + (theta(ud, i, s, j) + theta(ud, i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
    enddo !j
    
    num(ud, i, s) = nnn
    
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine NN


!電荷密度
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


!収束チェック
subroutine CV(N, ch, rhov, rhop, num, cvg, cvn)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N), intent(in) :: rhov
 real*8, dimension(3, N), intent(in) :: rhop
 real*8, dimension(3, N, 8), intent(in) :: num
 real*8, dimension(3, N), intent(out) :: cvg
 real*8, dimension(3, N) :: nume, numi
 real*8, intent(out) :: cvn
 integer :: ud, i
 
 cvg = 0.d0
 cvn = 0.d0
 
 do ud = 1, 3
  do i = 1, N
   nume(ud, i) = abs(ch(2)*num(ud, i, 2) + ch(7)*num(ud, i, 7) + ch(8)*num(ud, i, 8))
   numi(ud, i) = ch(1)*num(ud, i, 1) + ch(3)*num(ud, i, 3) + ch(4)*num(ud, i, 4) + &
               & ch(5)*num(ud, i, 5) + ch(6)*num(ud, i, 6)
  enddo !i
 enddo !ud
 
 do ud = 1, 3
  do i = 2, N-1
   cvg(ud, i) = rhov(ud, i)**2.d0/(nume(ud, i)*nume(ud, i))
   if(ud == 1) cvn = cvn + abs(cvg(ud, i))
  enddo !i
  cvg(ud, 1) = rhov(ud, 1)**2.d0/(nume(ud, 1)*nume(ud, 1))
  cvg(ud, N) = rhov(ud, N)**2.d0/(nume(ud, N)*nume(ud, N))
  !if(ud == 1) cvn = cvn + abs(cvg(ud, 1)) + abs(cvg(ud, N))
 enddo !ud
 
 cvn = (cvn/dble(N-2))**(1.d0/2.d0)
 
 return
 
end subroutine CV


!Newton法(npe)
subroutine Newtonphi(N, pep, cvg, npe)
 integer, intent(in) :: N
 real*8, dimension(3, N), intent(in) :: pep, cvg
 real*8, dimension(N), intent(out) :: npe
 real*8 :: CC
 integer :: i
 
 npe = pep(1, :)
 
 do i = 2, N-1
  CC = (pep(2, i)-pep(3, i))/(cvg(2, i)-cvg(3, i))*cvg(1, i)
  if(abs(CC) <= abs(pep(1, i)*1.01d0) .and. pep(1, i) /= 0.d0) then
    npe(i) = pep(1, i) - CC
  endif
  if(pep(1, i) == 0.d0 .and. abs(CC) < 1.d-5) then
    npe(i) = pep(1, i) - CC
  endif
  if(pep(1, i) == 0.d0 .and. abs(CC) >= 1.d-5) then
    npe(i) = pep(1, i) - CC/abs(CC)*1.d-5
  endif
  if(abs(CC) > abs(pep(1, i)*1.01d0) .and. pep(1, i) /= 0.d0) then
    if(abs(pep(1, i)) < 1.d-6) npe(i) = pep(1, i) - CC/abs(CC)*abs(pep(1, i))
    if(abs(pep(1, i)) >= 1.d-6) npe(i) = pep(1, i) - CC/abs(CC)*abs(pep(1, i)/2.d0)
  endif
  
 enddo
 
 return
 
end subroutine Newtonphi


!Newton法(nsig)
subroutine Newtonsig(N, ch, sigma, sigmab, num, cvg, nsig)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch, sigma
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8), intent(in) :: num
 real*8, dimension(3, N), intent(in) :: cvg
 real*8, dimension(8), intent(out) :: nsig
 integer :: s, siMn, siIn
 real*8 :: siM, siI
 
 nsig = sigma
 siM = 0.d0
 siI = 0.d0
 siMn = 0
 siIn = 0
 
 do s = 1, 8
  
  !イオ起源
  if(sigmab(s) == 1 .and. num(1, N, s) == 0.d0 .and. siM < abs(ch(s)*num(1, 1, s))) then
    siM = abs(ch(s)*num(1, 1, s))
    siMn = s
  
  !木星起源
   else if(sigmab(s) == N .and. num(1, 1, s) == 0.d0 .and. siI < abs(ch(s)*num(1, N, s))) then
    siI = abs(ch(s)*num(1, N, s))
    siIn = s
  endif
 enddo !s
 
 if(siMn /= 0) then
   nsig(siMn) = sigma(siMn) - (num(2, 1, siMn)-num(3, 1, siMn))/(cvg(2, 1)-cvg(3, 1))*cvg(1, 1)
   if(nsig(siMn) <= 0.d0 .or. isnan(nsig(siMn))) print *, "nsigM-ERROR"
   if(cvg(2, 1)-cvg(3, 1) == 0.d0) print *, "ERROR"
 endif
 if(siIn /= 0) then
   nsig(siIn) = sigma(siIn) - (num(2, N, siIn)-num(3, N, siIn))/(cvg(2, N)-cvg(3, N))*cvg(1, N)
   if(nsig(siIn) <= 0.d0 .or. isnan(nsig(siIn))) print *, "nsigI-ERROR"
 endif
 
 return
 
end subroutine Newtonsig

