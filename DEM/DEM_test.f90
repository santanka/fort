program DEM_test

implicit none

!grid数
integer, parameter :: N = 1000


!doで使う数
integer :: i, k


!DEM定数枠
real*8, parameter :: GG = 6.6743015d-11 !万有引力定数
real*8, parameter :: kB = 1.38064852d-23 !Boltzmann constant


!DEM可変定数枠
real*8, parameter :: C = 1.d0 !電子とイオンの温度比
real*8, parameter :: D = 1.d7 !温度調整用可変距離


!DEM入力定数枠
!!星関連
real*8, parameter :: ME = 5.9722d24 !星(地球)の質量
real*8, parameter :: RE = 6.371d6 !星(地球)の半径
real*8, parameter :: RAV = 7.292115024135739d-5 !星(地球)の自転角速度
real*8, parameter :: R0 = 6.871d6 !半径+電離圏高度(500km仮定)
real*8, parameter :: L = 4.d0 !L値
real*8, parameter :: Req = RE*L


!!プラズマ関連
integer, parameter :: kind = 3 !プラズマ種数
real*8, dimension(kind, 5) :: pl !プラズマデータまとめ
real*8, parameter :: e = 1.602176634d-19 !電気素量


!DEM変数・関数枠
real*8, dimension(kind) :: HH !電離圏でのスケールハイト
real*8 :: sM !sの上限
real*8 :: pass !function、sを求める
real*8, dimension(N) :: ss !pass
real*8 :: Ts !function、プラズマの温度関数
real*8, dimension(N) :: z !temperature-modified geopotential height
real*8, dimension(kind, N) :: nn !プラズマ数密度profile
real*8 :: lam0 !footprintでの磁気緯度
real*8, dimension(N) :: lambda




!プラズマデータまとめ
pl(1, 1) = 9.1093897d-31 !電子の質量
pl(1, 2) = 2.001d11 !1.00646*2.d11 !電子の電離圏数密度
pl(1, 3) = -1.d0 !電子の価数
pl(1, 4) = 0.1 * e !1.5d3*kB !電子の電離圏温度(J)
pl(1, 5) = 1.d-1 * e !1.5d3*kB !電子の磁気赤道温度(J)

print *, pl(1, 4)

pl(2, 1) = 2.677950266103d-26 !酸素イオンの質量
pl(2, 2) = 2.d11 !酸素イオンの電離圏数密度
pl(2, 3) = 1.d0 !酸素イオンの価数
pl(2, 4) = pl(1, 4) * C !酸素イオンの電離圏温度(J)
pl(2, 5) = pl(1, 5) * C !酸素イオンの磁気赤道温度(J)

pl(3, 1) = 1.67262192369d-27 !水素イオンの質量
pl(3, 2) = 1.d8 !1.6d-4*2.d11 !水素イオンの電離圏数密度
pl(3, 3) = 1.d0 !水素イオンの価数
pl(3, 4) = pl(1, 4) * C !水素イオンの電離圏温度(J)
pl(3, 5) = pl(1, 5) * C !水素イオンの磁気赤道温度(J)

!pl(4, 1) = 3.34524384738d-27 !ヘリウムイオンの質量
!pl(4, 2) = 6.3d-3*2.d11 !0.d0 !ヘリウムイオンの電離圏数密度
!pl(4, 3) = 1.d0 !ヘリウムイオンの価数
!pl(4, 4) = pl(1, 4) * C !ヘリウムイオンの電離圏温度(J)
!pl(4, 5) = pl(1, 5) * C !ヘリウムイオンの磁気赤道温度(J)


!スケールハイト
call scaleheight(kind, pl, GG, ME, R0, HH)

print *, HH

!lam0決定
lam0 = asin(sqrt(1.d0-R0/Req))

!lambda
do i = 1, N
 lambda(i) = lam0*dble(i-1)/dble(N-1)
enddo !i

!ss, sMの決定
do i = 1, N
 ss(i) = pass(Req, R0, lambda(i))
enddo !i
sM = ss(1)

!z
do i = 1, N
 call geo(R0, Req, lambda(i), RAV, GG, ME, pl(1, 4), sM, pl(1, 5), D, z(i))
enddo !i

!プラズマ数密度
do i = 1, N
 call ne(kind, C, lambda(i), ss(i), sM, D, z(i), pl, HH, nn(1, i))
 do k = 2, kind
  call ni(k, kind, C, lambda(i), ss(i), sM, D, z(i), pl, HH, nn(k, i))
 enddo !k
enddo !i

do i = 1, N
 print *, lambda(i), ss(i), z(i), Ts(ss(i), sM, pl(1, 5), pl(1, 4), D), nn(1, i), nn(2, i), nn(3, i), &
& (nn(2, i)+nn(3, i))/nn(1, i)
enddo !i

open(11, file = "DEM_E_L=4.csv")
do i = 1, N
 write(11, 21) lambda(i), ss(i), z(i), Ts(ss(i), sM, pl(1, 5), pl(1, 4), D)/e, nn(:, i)
enddo !i

call Alfven(kind, N, pl, nn, lambda, Req, RE)


21 format(1PE25.15E3, 7(',', 1PE25.15E3))


end program DEM_test





!function
!!sを求める
real*8 function pass(Req, R0, lambda)
 real*8, intent(in) :: Req, R0, lambda
 real*8 :: lam0
 
 lam0 = asin(sqrt(1.d0-R0/Req))
 
 pass = Req/2.d0/sqrt(3.d0)*(sqrt(3.d0)*sin(lam0)*sqrt(3.d0*sin(lam0)**2.d0+1.d0) + asinh(sqrt(3.d0)*sin(lam0)))
 pass = pass - Req/2.d0/sqrt(3.d0)*(sqrt(3.d0)*sin(lambda)*sqrt(3.d0*sin(lambda)**2.d0+1.d0) + asinh(sqrt(3.d0)*sin(lambda)))
 
 return
end function


!!プラズマの温度関数
real*8 function Ts(s, sM, TsT, Ts0, D)
 implicit none
 real*8, intent(in) :: s, sM, TsT, Ts0, D
 
 Ts = Ts0 + (TsT-Ts0)*(1.d0-exp(-s/D))/(1.d0-exp(-sM/D))
 !Ts = Ts0
 
 return
end function


!subroutine
!!スケールハイト
subroutine scaleheight(kind, pl, GG, ME, R0, HH)
 implicit none
 integer, intent(in) :: kind
 real*8, dimension(kind, 5), intent(in) :: pl
 real*8, intent(in) :: GG, ME, R0
 real*8, dimension(kind), intent(out) :: HH
 
 integer :: k
 
 do k = 1, kind
  HH(k) = R0**2.d0*pl(k, 4)/GG/ME/pl(k, 1)
 enddo !k
 
 return
end subroutine scaleheight


!!温度補正ジオポテンシャル高度
subroutine geo(R0, Req, lambda, RAV, GG, ME, Ts0, sM, TsT, D, z)
 implicit none
 real*8, intent(in) :: R0, Req, lambda, RAV, GG, ME, Ts0, sM, TsT, D
 real*8, intent(out) :: z
 real*8 :: Ts, pass !function
 real*8 :: ss1, ss2, lam0, dl, nl1, nl2, XX, YY
 integer :: NX = 1000
 integer :: i
 
 XX = 0
 YY = 0
 
 lam0 = asin(sqrt(1.d0-R0/Req))
 dl = (lam0-lambda)/dble(NX-1)
 
 do i = 1, NX-1
  nl1 = lambda + dl*dble(i-1)
  nl2 = lambda + dl*dble(i)
  ss1 = pass(Req, R0, nl1)
  ss2 = pass(Req, R0, nl2)
  XX = XX + (sin(nl1)/cos(nl1)**3.d0*Ts0/Ts(ss1, sM, TsT, Ts0, D) &
      & +sin(nl2)/cos(nl2)**3.d0*Ts0/Ts(ss2, sM, TsT, Ts0, D))/2.d0*dl
  YY = YY + (sin(nl1)*cos(nl1)**5.d0*Ts0/Ts(ss1, sM, TsT, Ts0, D) &
      & +sin(nl2)*cos(nl2)**5.d0*Ts0/Ts(ss2, sM, TsT, Ts0, D))/2.d0*dl
 enddo !i
 
 z = 2.d0*R0**2.d0/Req*XX
 z = z - 3.d0*(RAV*R0*Req)**2.d0/GG/ME*YY
 
 return
end subroutine geo


!!電子数密度
subroutine ne(kind, C, lambda, ss, sM, D, z, pl, HH, nn)
 implicit none
 integer, intent(in) :: kind
 real*8, intent(in) :: C, lambda, ss, sM, D, z
 real*8, dimension(kind, 5), intent(in) :: pl
 real*8, dimension(kind) :: HH
 real*8, intent(out) :: nn
 real*8 :: Ts !function
 real*8 :: EE
 integer :: k
 
 EE = 0.d0
 do k = 2, kind
  EE = EE + pl(k, 3)*pl(k, 2)/pl(1, 2)*exp(-z/HH(k))
 enddo !k
 
 nn = pl(1, 2)*pl(1, 4)/Ts(ss, sM, pl(1, 5), pl(1, 4), D)
 nn = nn * exp(-C/(C+1.d0)*z/HH(1))
 nn = nn * EE**(1.d0/(C+1.d0))
 
 return
end subroutine ne


!!イオン数密度
subroutine ni(kk, kind, C, lambda, ss, sM, D, z, pl, HH, nn)
 implicit none
 integer, intent(in) :: kk, kind
 real*8, intent(in) :: C, lambda, ss, sM, D, z
 real*8, dimension(kind, 5), intent(in) :: pl
 real*8, dimension(kind) :: HH
 real*8, intent(out) :: nn
 real*8 :: Ts !function
 real*8 :: EE
 integer :: k
 
 EE = 0.d0
 do k = 2, kind
  EE = EE + pl(k, 3)*pl(k, 2)/pl(1, 2)*exp(-z/HH(k))
 enddo !k
 
 nn = pl(kk, 2)*(pl(1, 4)/Ts(ss, sM, pl(1, 5), pl(1, 4), D))**C
 nn = nn * exp(-C/(C+1.d0)*z/HH(1))
 nn = nn * exp(-z/HH(kk))
 nn = nn / EE**(C/(1.d0+C))
 
 return
end subroutine ni


!!Alfven速度
subroutine Alfven(kind, N, pl, nn, lambda, Req, RE)
 implicit none
 integer, intent(in) :: kind, N
 real*8, dimension(kind, 5) :: pl
 real*8, dimension(kind, N), intent(in) :: nn
 real*8, dimension(N), intent(out) :: lambda
 real*8, intent(in) :: Req, RE
 real*8, dimension(N) :: BB, VA, Vr, isum
 integer :: i, s
 real*8 :: L, BE
 real*8, parameter :: mu0 = 1.25663706143592d-6 !真空の透磁率
 real*8, parameter :: c = 2.99792458d8 !光速
 real*8, parameter :: ME = 7.8d22 !地球の磁気双極子モーメント
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 L = Req/RE
 BE = mu0*ME/4.d0/pi/RE**3.d0
 isum = 0.d0
 
 do i = 1, N
  BB(i) = BE/L**3.d0*sqrt(1.d0+3.d0*sin(lambda(i))**2.d0)/cos(lambda(i))**6.d0
  do s = 1, kind
   if(pl(s, 3) > 0.d0) isum(i) = isum(i) + nn(s, i)*pl(s, 1)
  enddo !s
 enddo !i
 
 do i = 1, N
  VA(i) = BB(i)/sqrt(mu0*isum(i))
  Vr(i) = VA(i)/sqrt(1.d0 + (VA(i)/c)**2.d0)
 enddo !i
 
 open(53, file = "DEM_AW_E_L=4.csv")
 do i = 1, N
  write(53, 60) VA(i), Vr(i), Vr(i)/c
 enddo !i
 close(53)
 
 return
 
 60 format(1PE25.15E3, 2(',', 1PE25.15E3))
 
end subroutine Alfven
