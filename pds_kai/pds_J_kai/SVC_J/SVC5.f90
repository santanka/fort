program SVC5

implicit none

!全体を通しての文字
integer :: i, j, s !do用
integer, parameter :: N = 79 !grid数
integer, parameter :: Z = 60 !速度空間grid数
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
real*8, dimension(N) :: npe
integer :: itn
real*8, dimension(3, N) :: pep, rhov, rhop, cvg
real*8, dimension(3, N, 8) :: UU, num
real*8, dimension(3, N, 8, Z) :: UB
real*8, dimension(8, Z) :: mu
real*8, dimension(3, N, 8, Z, 3) :: a3
real*8, dimension(3, N, 8, Z) :: theta
real*8, dimension(3, 2) :: sigx
real*8, dimension(8) :: nsig
real*8 :: cvn, DDD, CC


!mag_FAの抽出
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i) !LMD:磁気緯度,BB:磁束密度,dIo:Ioとの距離
enddo
close(40)

B = BB(1)
BB = BB

!SVC1_parameterの抽出
open(30, file="SVC1_parameter.csv", action="read", status="old")
do i = 1, 12
 read(30, *) dummy, parameter(i)
end do
close(30)

e = parameter(1) !電気素量
c = parameter(2) !光速
m = parameter(3) !電子の質量
mu0 = parameter(4) !真空の透磁率
GG = parameter(5) !万有引力定数
alpha = parameter(6) !alpha
MJ = parameter(7) !Jupiterの磁気双極子モーメント
MassJ = parameter(8) !Jupiterの質量
RJ = parameter(9) !Jupiter半径
omJ = parameter(10) !Jupiter自転角周波数
Req = parameter(11) !Jupiterのreq
MIo = parameter(12) !Ioの質量

ep0 = 1.d0/mu0/c/c !真空の誘電率

!SVC_ICの抽出
open(10, file="SVC_IC_J_30kV_75_50.csv", action="read", status="old")
read(10, *) !1行読み飛ばし
do i = 1, N
 read(10, *) dummy, pe0(i), dk(i), dR(i)
end do
close(10)

do i = 1, N-1
 dx = dk(i+1) - dk(i)
enddo

pe0 = pe0 !初期静電ポテンシャル
dk = dk !i=1からの距離
dx = dx !grid間距離

!SVC_BCの抽出
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1行読み飛ばし
do i = 1, 8
 read(20, *) dummy, Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do
close(20)

!Ta:Tpara/Tperp,so:起源,ch:電荷
Nd = Nd !数密度
Tp = Tp*e !perp方向温度
mass = mass !質量
ch = ch*e
sigma = Nd

!つづきからはじめる▼
!open(15, file="SVC_penum.csv", action="read", status="old")
!do i = 1, N
 !if(i >= 45 .and. i <= 55) read(15, *) dummy
 !if(i < 45 .or. 55 < i) read(15, *) pe0(i)
! read(15, *) CC
 !if(i <= 62) pe0(i) = CC/2
 !if(i > 62) pe0(i) = CC
! pe0(i) = CC
!enddo !i
!close(15)

!open(18, file="SVC_nsig.csv", action="read", status="old")
!do s = 1, 8
! read(18, *) sigma(s)
!enddo !s
!close(18)



!sigmab作成
do i = 1, 8
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo


!iteration スタート

pep(1, :) = pe0
itn = 0
DDD = 1.d0

!open(33, file = "SVC_min_30kV.csv")
!read(33, *) DDD
!------------------------------
!do i = 1, N
! read(33, *) pep(1, i)
!enddo !i
!do s = 1, 8
! read(33, *) sigma(s)
!enddo !s
!------------------------------
!close(33)

call AI(Z, alpha, N, Tp, BB, sigmab, mu)

do !itn
 itn = itn + 1
 
 !ポテンシャル
 call EP(N, ch, pep(1, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(1, :, :), itn)
 
 !静電ポテンシャル差分
 call pepm(pep, e, m, c, N, UU(1, :, :), ch)
 
 !ポテンシャル差分
 call EP(N, ch, pep(2, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(2, :, :), 0)
 call EP(N, ch, pep(3, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(3, :, :), 0)
 
 !エネルギー(位置+Kperp)
 call UUBB(N, Z, UU, BB, mu, UB)
 
 !accessibility(1:min, 2:lim, 3:max)
 call aaa(N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 
 !sigma差分
 call sigud(sigma, sigx)
 
 !aについて積分(theta)
 call TH(N, Z, a3, sigma, sigx, Tp, Ta, sigmab, UU, BB, mu, theta)
 
 !muについて積分(num)
 call NN(N, Z, mu, theta, num)
 
 !電荷密度
 call RR(N, ch, num, rhov)
 
 !Poisson方程式
 call Poisson(N, pep, ep0, dx, rhop)
 
 !収束チェック
 call CV(N, ch, rhov, num, cvg, cvn)
 if(mod(itn, 50) == 1) then
  print *, itn, cvn
  print *, "    "
 endif
  
 !Alfven速度
 call Alf(N, mu0, BB, mass, num, c)
 
 !ファイル化処理
 open(60, file="SVC_pote.csv")
 do i = 1, N
  write(60, 92) pep(1, i), UU(1, i, :)
 enddo
 close(60)
 
 open(70, file="SVC_a3.csv")
 do i = 1, N
  write(70, 62) a3(1, i, 1, 1, :), alpha*Tp(1)/BB(sigmab(1))
 enddo !i
 close(70)
 
 !最低収束値
  if(DDD > cvn) then
    DDD = cvn
    open(33, file = "SVC_min.csv")
    write(33, '(1PE25.15E3)') cvn
    do i = 1, N
      write(33, '(1PE25.15E3)') pep(1, i)
    enddo !i
    do s = 1, 8
      write(33, '(1PE25.15E3)') sigma(s)
    enddo !s
    close(33)
    open(50, file = "SVC_penum.csv")
    do i = 1, N
      write(50, 72) pep(1, i), num(1, i, :), rhov(1, i), rhop(1, i), cvg(1, i)
    enddo !i
    close(50)
  endif
 
 !収束チェック
 if(cvn < 1.d-6) then
   print *, "finish"
   exit
 endif
 
 !Newton法
 call Newtonphi(N, pep, cvg, npe, c, m, e, cvn, itn)
 call Newtonsig(N, sigx, sigma, cvg, nsig, cvn)
 
 !NaNチェック
 do i = 2, N-1
  if(isnan(npe(i)) .or. isnan(cvg(1, i))) then
    print *, "NaN発生"
    stop
  endif
  if(npe(i) == npe(i)-1.d0 .or. cvg(1, i) == cvg(1, i)-1.d0) then
    print *, "infinity発生"
    stop
  endif
 enddo !i
 
 
 !更新
 pep(1, :) = npe
 sigma = nsig
 
 
enddo !itn


62 format(1PE25.15E3, 3(',', 1PE25.15E3))
72 format(1PE25.15E3, 11(',', 1PE25.15E3))
82 format(1PE25.15E3, 21(',', 1x, 1PE25.15E3))
92 format(1PE25.15E3, 8(',', 1PE25.15E3))


end program SVC5





!subroutine, function

!断熱不変量
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


!ポテンシャル
subroutine EP(N, ch, pe, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU, itn)
 implicit none
 integer, intent(in) :: N, itn
 real*8, dimension(8), intent(in) :: ch, mass
 real*8, dimension(N), intent(in) :: pe, LMD, dIo
 real*8, intent(in) :: GG, MassJ, Req, omJ, MIo
 real*8, dimension(N, 8), intent(out) :: UU
 integer ::i, s
 
 do i = 1, N
  do s = 1, 8
   !木星重力
   UU(i, s) = - GG*MassJ*mass(s)/Req/(cos(LMD(i))**2.d0)
   !木星遠心力
   UU(i, s) = UU(i,s) - mass(s)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0
   !イオ重力
   UU(i, s) = UU(i, s) - GG*MIo*mass(s)/dIo(i)
  enddo !s
 enddo !i
 
 do s = 1, 8
  UU(:, s) = UU(:, s) - UU(1, s)
 enddo !s
 
 if(itn == 1) then
   open(44, file="SVC_phizeropote.csv")
   do i = 1, N
    write(44, 86) UU(i, :)
   enddo !s
   close(44)
 endif
 
 do i = 1, N
  do s = 1, 8
   !静電
   UU(i, s) = UU(i, s)+ch(s)*pe(i)
  enddo !s
 enddo !i
 
 return
 
 86 format(E25.15E3, 7(',', E25.15E3))
end subroutine EP


!静電ポテンシャル差分
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
  !if(pep(1, i) == 0.d0) dp(i) = 1.d-3
  !if(pep(1, i) >2.5d4) dp(i) = pep(1, i)*1.d-10
  !if(pep(1, i) <= 2.5d4 .and. pep(1, i) > 1.d4) dp(i) = pep(1, i)*1.d-10
  !if(pep(1, i) <= 1.d4) dp(i) = pep(1, i)*1.d-4
  dp(i) = 1d-7!pep(1, i)*1.d-4
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


!エネルギー(位置+Kperp)
subroutine UUBB(N, Z, UU, BB, mu, UB)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, 8, Z), intent(out) :: UB
 integer :: ud, i, s, j
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     UB(ud, i, s, j) = UU(ud, i, s) + BB(i)*mu(s, j)
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine UUBB


subroutine aaa(N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z
 real*8, intent(in) :: alpha
 real*8, dimension(8), intent(in) :: Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8, Z) :: UB
 real*8, dimension(3, N, 8, Z, 3) :: a3
 integer, dimension(3, N, 8, Z) :: lima
 integer :: ud, i, s, j, k
 real*8 :: CC
 
 lima = 0
 
 !amin
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
    
    !イオ起源
    if(sigmab(s) == 1) then
      if(i == 1) a3(ud, i, s, j, 1) = 0.d0
      if(i == 2) then
        if(UB(ud, 2, s, j) >= UB(1, 1, s, j)) a3(ud, i, s, j, 1) = 0.d0
        if(UB(ud, 2, s, j) < UB(1, 1, s, j)) then
          a3(ud, i, s, j, 1) = sqrt(UB(1, 1, s, j) - UB(ud, i, s, j))
        endif
      endif
      if(i >= 3) then
        CC = UB(ud, i, s, j)
        do k = 1, i-1
         if(UB(1, k, s, j) > CC) CC = UB(1, k, s, j)
        enddo !k
        a3(ud, i, s, j, 1) = sqrt(CC - UB(ud, i, s, j))
      endif
    endif
    
    !木星起源
    if(sigmab(s) == N) then
      if(i == N) a3(ud, i, s, j, 1) = 0.d0
      if(i == N-1) then
        if(UB(ud, N-1, s, j) >= UB(1, N, s, j)) a3(ud, i, s, j, 1) = 0.d0
        if(UB(ud, N-1, s, j) < UB(1, N, s, j)) then
          a3(ud, N-1, s, j, 1) = sqrt(UB(1, N, s, j) - UB(ud, N-1, s, j))
        endif
      endif
      if(i <= N-2) then
        CC = UB(ud, i, s, j)
        do k = i+1, N
         if(UB(1, k, s, j) > CC) CC = UB(1, k, s, j)
        enddo !k
        a3(ud, i, s, j, 1) = sqrt(CC - UB(ud, i, s, j))
      endif
    endif
    
    if(a3(ud, i, s, j, 1) == 0.d0) a3(ud, i, s, j, 1) = 1.d-15
    
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 !alim
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     
     CC = 0.d0
     
     !イオ起源
     if(sigmab(s) == 1) then
       if(i == 1) a3(ud, i, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
       if(i /= 1) then
         CC = UB(1, 1, s, j) + alpha*Tp(s)*Ta(s) - UB(ud, i, s, j)
         if(CC <= 0.d0) a3(ud, i, s, j, 2) = 0.d0
         if(CC > 0.d0) a3(ud, i, s, j, 2) = sqrt(CC)
         
         if(a3(1, i-1, s, j, 2) < a3(1, i-1, s, j, 1)) a3(ud, i, s, j, 2) = 0.d0
       endif
     endif
     
     !木星起源
     if(sigmab(s) == N) then
       k = N+1-i
       if(k == N) a3(ud, k, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
       if(k /= N) then
         CC = UB(1, N, s, j) + alpha*Tp(s)*Ta(s) - UB(ud, k, s, j)
         if(CC <= 0.d0) a3(ud, k, s, j, 2) = 0.d0
         if(CC > 0.d0) a3(ud, k, s, j, 2) = sqrt(CC)
         
         if(a3(1, k+1, s, j, 2) < a3(1, k+1, s, j, 1)) a3(ud, k, s, j, 2) = 0.d0
       endif
     endif
     
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 
 !accessibilityの調整
 do ud = 1, 3
  do s = 1, 8
   do j = 1, Z
    
    !イオ起源
    if(sigmab(s) == 1) then
      do i = 1, N
       if(i /= 1 .and. lima(ud, i-1, s, j) == 1) lima(ud, i, s, j) = 1
       if(a3(ud, i, s, j, 2) <= a3(ud, i, s, j, 1)) lima(ud, i, s, j) = 1
      enddo !i
    endif
    
    !木星起源
    if(sigmab(s) == N) then
      do i = 1, N
       if(i /= 1 .and. lima(ud, N+2-i, s, j) == 1) lima(ud, N+1-i, s, j) = 1
       if(a3(ud, N+1-i, s, j, 2) <= a3(ud, N+1-i, s, j, 1)) lima(ud, N+1-i, s, j) = 1
      enddo !i
    endif
    
   enddo !j
  enddo !s
 enddo !ud
 
 
 !amax
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z
     
     CC = 0.d0
     
     !イオ起源
     if(sigmab(s) == 1) then
       if(i == N .or. lima(ud, i, s, j) == 1) a3(ud, i, s, j, 3) = 0.d0
       if(i+1 == N .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, N, s, j) - UB(ud, N-1, s, j)
         if(CC <= 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC > 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i+1 < N .and. lima(ud, i, s, j) == 0) then
         do k = i+1, N
          if(CC < UB(1, k, s, j)) CC = UB(1, k, s, j)
         enddo !k
         CC  = CC - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
     endif
     
     !木星起源
     if(sigmab(s) == N) then
       if(i == 1 .or. lima(ud, i, s, j) == 1) a3(ud, i, s, j, 3) = 0.d0
       if(i == 2 .and. lima(ud, i, s, j) == 0) then
         CC = UB(1, 1, s, j) - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
       if(i > 2 .and. lima(ud, i, s, j) == 0) then
         do k = 1, i-1
          if(CC < UB(1, k, s, j)) CC = UB(1, k, s, j)
         enddo !k
         CC = CC - UB(ud, i, s, j)
         if(CC < 0.d0) a3(ud, i, s, j, 3) = 0.d0
         if(CC >= 0.d0) a3(ud, i, s, j, 3) = sqrt(CC)
       endif
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
 return

end function ff


!aについて積分(theta)
subroutine TH(N, Z, a3, sigma, sigx, Tp, Ta, sigmab, UU, BB, mu, theta)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(3, N, 8, Z, 3), intent(in) :: a3
 real*8, dimension(8), intent(in) :: sigma, Tp, Ta
 real*8, dimension(3, 2), intent(in) :: sigx
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, 8, Z), intent(out) :: theta
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 integer :: ud, i, s, j, p, numx
 real*8 :: ff !function
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    do j = 1, Z       
     thetaL = 0.d0
     thetaM = 0.d0
     
     if((i == N .and. s == 2) .or. (i == 1 .and. s == 3)) then
       if(s == 2) numx = 1
       if(s == 3) numx = 2
       
       !alim
       if(a3(1, i, s, j, 2) > a3(1, i, s, j, 1)) then
         do p = 1, Z/2
          aL(p) = a3(1, i, s, j, 1)*(a3(1, i, s, j, 2)/a3(1, i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
         enddo !p
         do p = 1, Z/2-1
          thetaL = thetaL + (ff(sigx(ud, numx), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(1, i, s), UU(1, i, s), mu(s, j), aL(p)) &
           & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(1, i, s), UU(1, i, s), mu(s, j), aL(p+1))) &
           & /2.d0*abs(aL(p+1)-aL(p))
         enddo !p
       endif
       !amax
       if(a3(ud, i, s, j, 3) > a3(ud, i, s, j, 1)) then
         do p = 1, Z/2
          aM(p) = -(a3(1, i, s, j, 1)*(a3(1, i, s, j, 3)/a3(1, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
         enddo !p
         do p = 1, Z/2-1
          thetaM = thetaM + (ff(sigx(ud, numx), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(1, i, s), UU(1, i, s), mu(s, j), aM(p)) &
           & + ff(sigx(ud, numx), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(1, i, s), UU(1, i, s), mu(s, j), aM(p+1))) &
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
         !iが起源ではない
         if(sigmab(s) /= i) then
           do p = 1, Z/2-1
            thetaL = thetaL + &
             & (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aL(p)) + &
             & ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aL(p+1))) &
             & /2.d0*abs(aL(p+1)-aL(p))
           enddo !p
         endif
         !iが起源
         if(sigmab(s) == i) then
           do p = 1, Z/2-1
            thetaL = thetaL + &
             &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aL(p)) + &
             &ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aL(p+1))) &
             & /2.d0*abs(aL(p+1)-aL(p))
           enddo !p
         endif
       endif
       
       !amax
       if(a3(ud, i, s, j, 3) > a3(ud, i, s, j, 1)) then
         do p = 1, Z/2
          aM(p) = -(a3(ud, i, s, j, 1)*(a3(ud, i, s, j, 3)/a3(ud, i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
         enddo !p
         !iが起源ではない
         if(sigmab(s) /= i) then
           do p = 1, Z/2-1
            thetaM = thetaM + &
             &(ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aM(p)) + &
             &ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(1, sigmab(s), s), mu(s, j), aM(p+1))) &
             & /2.d0*abs(aM(p)-aM(p+1))
           enddo !p
         endif
         !iが起源
         if(sigmab(s) == i) then
           do p = 1, Z/2-1
            thetaM = thetaM + (ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aM(p)) &
             & + ff(sigma(s), BB(i), BB(sigmab(s)), Ta(s), Tp(s), UU(ud, i, s), UU(ud, i, s), mu(s, j), aM(p+1))) &
             & /2.d0*abs(aM(p)-aM(p+1))
           enddo !p
         endif
       endif
       
       theta(ud, i, s, j) = thetaL + thetaM
     endif
    enddo !j
   enddo !s
  enddo !i
 enddo !ud
 
 return
 
end subroutine TH


!muについて積分(num)
subroutine NN(N, Z, mu, theta, num)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(3, N, 8, Z), intent(in) :: theta
 real*8, dimension(3, N, 8), intent(out) :: num
 integer :: ud, i, s, j
 real*8 :: nnn
 
 do ud = 1, 3
  do i = 1, N
   do s = 1, 8
    
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


!ポアソン方程式
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


!収束チェック
subroutine CV(N, ch, rhov, num, cvg, cvn)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch
 real*8, dimension(3, N), intent(in) :: rhov
 real*8, dimension(3, N, 8), intent(in) :: num
 real*8, dimension(3, N), intent(out) :: cvg
 real*8, intent(out) :: cvn
 integer :: ud, i
 real*8, dimension(3, N) :: nume, numi
 
 cvg = 0.d0
 cvn = 0.d0
 
 do ud = 1, 3
  do i = 1, N
   nume(ud, i) = abs(ch(2)*num(ud, i, 2) + ch(7)*num(ud, i, 7) + ch(8)*num(ud, i, 8))
   numi(ud, i) = ch(1)*num(ud, i, 1) + ch(3)*num(ud, i, 3) + ch(4)*num(ud, i, 4) &
    & + ch(5)*num(ud, i, 5) + ch(6)*num(ud, i, 6)
  enddo !i
 enddo !ud
 
 do ud = 1, 3
  do i = 1, N
   cvg(ud, i) = rhov(ud, i)**2.d0/nume(ud, i)/numi(ud, i)
   if(ud == 1 .and. i /= 1 .and. i /= N) cvn = cvn + cvg(ud, i)
   cvg(ud, i) = sqrt(cvg(ud, i))
  enddo !i
 enddo !ud
 
 cvn = sqrt(cvn/dble(N-2))
 
 return
 
end subroutine CV


!Newton法(npe)
subroutine Newtonphi(N, Phih, cvg, nPhi, c, m, e, cvn, itn)
 implicit none
 integer, intent(in) :: N, itn
 real*8, dimension(3, N), intent(in) :: Phih, cvg
 real*8, dimension(N), intent(out) :: nPhi
 real*8, intent(in) :: c, m, e, cvn
 
 integer :: i, MV, k, CP, search, MV_min, MV_max
 real*8 :: CC, ser
 
 CP = 10000 !固定値なし
 MV = 39
 MV_min = 1
 MV_max = N

 do i = 2, N-1
  CC = 0.d0
  if(i /= CP) then
    if(cvg(2, i) < cvg(1, i) .or. cvg(3, i) < cvg(1, i)) then
      if(cvg(2, i) == cvg(1, i)) then
        CC = (Phih(3, i)-Phih(1, i))/(cvg(3, i)-cvg(1, i))*cvg(1, i)
       else if(cvg(3, i) == cvg(1, i)) then
        CC = (Phih(2, i)-Phih(1, i))/(cvg(2, i)-cvg(1, i))*cvg(1, i)
       else
        CC = ((Phih(3, i)-Phih(1, i))/(cvg(3, i)-cvg(1, i)) + (Phih(2, i)-Phih(1, i))/(cvg(2, i)-cvg(1, i)))/2.d0*cvg(1, i)
      endif
    endif
    if(abs(CC) < cvg(1, i) .and. cvg(1, i) <= 1.d0 .and. CC == CC) then
      nPhi(i) = Phih(1, i) - CC
     else if(abs(CC) >= cvg(1, i) .and. cvg(1, i) <= 1.d0) then
      if(Phih(1, i) > 1.d2 .and. Phih(1, i) < 2.9995d4) then
        nPhi(i) = Phih(1, i) - CC/abs(CC)*1.d1*sqrt(cvg(1, i))
       else
        if((abs(CC) < cvg(1, i) .and. cvg(1, i) <= 1.d-4) .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC
         else if(abs(CC) >= cvg(1, i) .and. cvg(1, i) <= 1.d-4 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC/abs(CC)*cvg(1, i)
         else if(abs(CC) < 1.d-4 .and. cvg(1, i) <= 1.d-3 .and. cvg(1, i) > 1.d-4 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC
         else if(abs(CC) >= 1.d-4 .and. cvg(1, i) <= 1.d-3 .and. cvg(1, i) > 1.d-4 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC/abs(CC)*1d-3
         else if(abs(CC) < 1.d-3 .and. cvg(1, i) <= 1.d-2 .and. cvg(1, i) > 1.d-3 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC
         else if(abs(CC) >= 1.d-3 .and. cvg(1, i) <= 1.d-2 .and. cvg(1, i) > 1.d-3 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC/abs(CC)*1d-3
         else if(abs(CC) < 1.d-2 .and. cvg(1, i) > 1.d-2 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC
         else if(abs(CC) >= 1.d-2 .and. cvg(1, i) > 1.d-2 .and. i /= CP) then
          nPhi(i) = Phih(1, i) - CC/abs(CC)*1d-2
        endif
      endif
     else if(abs(CC) < 1.d0 .and. cvg(1, i) > 1.d0) then
      if(Phih(1, i) > 1.d2 .and. Phih(1, i) < 2.9995d4) then
        nPhi(i) = Phih(1, i) - CC/abs(CC)*1.d1*sqrt(cvg(1, i))
      else
        nPhi(i) = Phih(1, i) - CC
      endif
     else if(abs(CC) >= 1.d0 .and. cvg(1, i) > 1.d0) then
      if(Phih(1, i) > 1.d2 .and. Phih(1, i) < 2.9995d4) then
        nPhi(i) = Phih(1, i) - CC/abs(CC)*1.d1*sqrt(cvg(1, i))
      else
        nPhi(i) = Phih(1, i) - CC/abs(CC)*1d0
      endif
    endif
  endif
  if(nPhi(i) /= nPhi(i)) then
    nPhi(i) = Phih(1, i)
  endif
 enddo !i
   
  search = 1
    do while(search == 1)
     search = 0
     do k = 1, 5
      do i = MV_min, MV-1
       if(nPhi(i) < nPhi(i+1)) then
         if(search == 0) search = 1
         if(i == MV_min) then
           nPhi(i+1) = 2.d0*nPhi(i)-nPhi(i+1)
          else !if(i /= CP .and. i+1 /= CP) then
           ser = nPhi(i)
           nPhi(i) = nPhi(i+1)
           nPhi(i+1) = ser
         endif
       endif
      enddo !i
     enddo
    enddo
  
  search = 1
    do while(search == 1)
     search = 0
     do k = 1, 5
      do i = MV, MV_max-2
       if(nPhi(i+1) < nPhi(i)) then
         if(search == 0) search = 1
         if(i == MV_max-1) then
           nPhi(i) = 2.d0*nPhi(i+1)-nPhi(i)
          else !if(i /= CP .and. i+1 /= CP) then
           ser = nPhi(i)
           nPhi(i) = nPhi(i+1)
           nPhi(i+1) = ser
         endif
       endif
      enddo !i
     enddo !k
    enddo
  
  nPhi(N) = Phih(1, N)


 open(70, file="SVC_npe.csv")
 do i = 2, N-1
  write(70, 52) (Phih(2, i)-Phih(3, i)), (cvg(2, i)-cvg(3, i)), cvg(1, i), cvg(2, i), cvg(3, i),&
&nPhi(i)-Phih(1, i), nPhi(i)
 enddo !i
 close(70)
 
 !do i = 2, N-1
 ! if(mod(itn, 77)+1 /= i .and. i /= mxn) npe(i) = pep(1, i)
 !enddo
 if (mod(itn, 50) == 1) then
  print *, mod(itn, 77)+1, nPhi(55)-Phih(1, 55)!, mxn, npe(mxn)-pep(1, mxn)
  print *, "  "
 endif

 return
 
 52 format(1PE25.15E3, 7(',', 1PE25.15E3))
 
end subroutine Newtonphi


!sigma差分作成
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

!Newton法(nsig)
subroutine Newtonsig(N, sigx, sigma, cvg, nsig, cvn)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(3, 2), intent(in) :: sigx
 real*8, dimension(8), intent(in) :: sigma
 real*8, dimension(3, N), intent(in) :: cvg
 real*8, dimension(8), intent(out) :: nsig
 real*8, intent(in) :: cvn
 integer :: s
 
 nsig = sigma
 
 nsig(2) = sigma(2) - &
&((sigx(2, 1) - sigx(1, 1))/(cvg(2, N) - cvg(1, N))+(sigx(3, 1) - sigx(1, 1))/(cvg(3, N) - cvg(1, N)))/2.d0*cvg(1, N)!*cvn
 nsig(3) = sigma(3) - &
&((sigx(2, 2) - sigx(1, 2))/(cvg(2, 1) - cvg(1, 1))+(sigx(3, 2) - sigx(1, 2))/(cvg(3, 1) - cvg(1, 1)))/2.d0*cvg(1, 1)!*cvn
 
 open(44, file="SVC_nsig.csv")
 do s = 1, 8
  write(44, 43) nsig(s), nsig(s)-sigma(s)
 enddo !s
 close(44)
 
 return
 
 43 format(E25.15E3, ',', E25.15E3)
 
end subroutine Newtonsig


!Alfven速度
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

