program SVC_E_L

implicit none

!///使用前の調整項目///
!全体を通しての文字
!参照ファイル(追加で、つづきから始める or min参照部分の調整)
!フォーマット


!全体を通しての文字
integer :: i, j, s !do用
integer, parameter :: N = 219 !grid数
integer, parameter :: Z = 60 !速度空間grid数
real*8, parameter :: pi = 4.d0*atan(1.d0) !円周率
character(len=80) :: dummy !使用しない列用
integer :: HATL = 208 !下段(not parameter)
integer :: LATL = 208 !上段(not parameter)、HATL <= LATL
integer, parameter :: kind = 5 !粒子種数
integer, parameter :: iodo = 3 !電離圏側の数密度調整
integer, parameter :: mado = 4 !磁気圏側の数密度調整

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
character(len=1), dimension(kind) :: so
real*8, dimension(kind) :: Nd, Tp, Ta, ch, mass

!sigma関連
real*8, dimension(kind) :: sigma
integer, dimension(kind) :: sigmab

!iteration関連
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
real*8 :: DD = 1.d0 !学習率


!mag_FAの抽出
open(40, file="mag_FA_E_L=4.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i) !LMD:磁気緯度,BB:磁束密度,dIo:Ioとの距離
enddo
close(40)

B = BB(1)
BB = BB

!SVC1_parameterの抽出
open(30, file="SVC_parameter_E_L=4.csv", action="read", status="old")
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
MJ = parameter(7) !Earthの磁気双極子モーメント
MassJ = parameter(8) !Earthの質量
RJ = parameter(9) !Earth半径
omJ = parameter(10) !Earth自転角周波数
Req = parameter(11) !Earthのreq
MIo = parameter(12) !衛星(今回なし)の質量

ep0 = 1.d0/mu0/c/c !真空の誘電率

!SVC_ICの抽出
open(10, file="SVC_IC_E_L=4.csv", action="read", status="old")
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
open(20, file="SVC_BC_E_L=4.csv", action="read", status="old")
read(20, *) !1行読み飛ばし
do i = 1, kind
 read(20, *) dummy, Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do
close(20)

!Ta:Tpara/Tperp,so:起源,ch:電荷
Nd = Nd !数密度
Tp = Tp*e !perp方向温度
mass = mass !質量
ch = ch*e
sigma = Nd


!sigmab作成
do i = 1, kind
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo

!iteration スタート

DDD = 0.d0

open(33, file = "SVC_min_E_L=4.csv")
read(33, *) DDD
!------------------------------
do i = 1, N
 read(33, *) pe0(i)
enddo !i
do s = 1, kind
 read(33, *) sigma(s)
enddo !s
read(33, *) HATL
read(33, *) LATL
!------------------------------
close(33)



!つづきからはじめる▼
!open(15, file="SVC_penum_E_L=4.csv", action="read", status="old")
!do i = 1, N
! read(15, *) CC
! pe0(i) = CC
!enddo !i
!close(15)


!open(18, file="SVC_nsig_E_L=4.csv", action="read", status="old")
!do s = 1, kind
! read(18, *) sigma(s)
!enddo !s
!close(18)

pep(1, :) = pe0
itn = 0


call AI(kind, Z, alpha, N, Tp, BB, sigmab, mu)

do !itn
 itn = itn + 1
 
 !ポテンシャル
 call EP(kind, N, ch, pep(1, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(1, :, :), itn)
 
 !静電ポテンシャル差分
 call pepm(kind, LATL, HATL, pep, e, m, c, N, UU(1, :, :), ch, DD)
 
 !ポテンシャル差分
 call EP(kind, N, ch, pep(2, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(2, :, :), 0)
 call EP(kind, N, ch, pep(3, :), GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU(3, :, :), 0)
 
 !エネルギー(位置+Kperp)
 call UUBB(kind, N, Z, UU, BB, mu, UB)
 
 !accessibility(1:min, 2:lim, 3:max)
 call aaa(kind, N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 
 !sigma差分
 call sigud(kind, iodo, mado, sigma, sigx)
 
 !aについて積分(theta)
 call TH(kind, iodo, mado, N, Z, a3, sigma, sigx, Tp, Ta, sigmab, UU, BB, mu, theta)
 
 !muについて積分(num)
 call NN(kind, N, Z, mu, theta, num)
 
 !電荷密度
 call RR(kind, N, ch, num, rhov)
 
 !Poisson方程式
 call Poisson(N, pep, ep0, dx, rhop)
 
 !収束チェック
 call CV(kind, N, ch, rhov, num, cvg, cvn, rhop)
 print *, itn, cvn
 
 !Alfven速度
 call Alf(kind, N, mu0, BB, mass, ch, num, c)
 
 !ファイル化処理
 open(50, file = "SVC_penum_E_L=4.csv")
 do i = 1, N
  write(50, 72) pep(1, i), num(1, i, :), rhov(1, i), &
                & rhop(1, i), cvg(1, i)
 enddo !i
 close(50)
 
 open(60, file="SVC_pote_E_L=4.csv")
 do i = 1, N
  write(60, 92) pep(1, i), UU(1, i, :)
 enddo
 close(60)
 
 
 !最低収束値
 if(DDD > cvn .or. DDD == 0.d0) then
   DDD = cvn
   open(33, file = "SVC_min_E_L=4.csv")
   write(33, '(1PE25.15E3)') cvn
   do i = 1, N
    write(33, '(1PE25.15E3)') pep(1, i)
   enddo !i
   do s = 1, kind
    write(33, '(1PE25.15E3)') sigma(s)
   enddo !s
   write(33, '(I3)') HATL
   write(33, '(I3)') LATL
   close(33)
 endif
 
 !収束チェック
 if((DD < 1.d-5 .or. cvn < 1.d-6) .and. check == 0) then
   print *, "finish"
   exit
 endif
 
 !Newton法
 call Newtonphi(LATL, HATL, N, pep, cvg, npe, c, m, e, cvn, itn, check, DD, XX)
 call Newtonsig(iodo, mado, kind, N, sigx, sigma, cvg, nsig, cvn, check, DD)
 
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


72 format(1PE25.15E3, 8(',', 1PE25.15E3)) !penum
92 format(1PE25.15E3, 6(',', 1PE25.15E3)) !pote


end program SVC_E_L





!subroutine, function

!断熱不変量
subroutine AI(kind, Z, alpha, N, Tp, BB, sigmab, mu)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp
 real*8, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(kind, Z), intent(out) :: mu
 integer :: i, j
 
 do i = 1, kind
  do j = 1, Z
   mu(i, j) = 1.d-20*(alpha*Tp(i)/BB(sigmab(i))/1.d-20)**(dble(j-1)/dble(Z-1))
  enddo !j
 enddo !i
 
 return
 
end subroutine AI


!ポテンシャル
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
   !木星重力
   UU(i, s) = - GG*MassJ*mass(s)/Req/(cos(LMD(i))**2.d0)
   !木星遠心力
   UU(i, s) = UU(i,s) - mass(s)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0
   !イオ重力
   if(MIo /= 0.d0) then
     UU(i, s) = UU(i, s) - GG*MIo*mass(s)/dIo(i)
   endif
  enddo !s
 enddo !i
 
 do s = 1, kind
  UU(:, s) = UU(:, s) - UU(N, s)
 enddo !s
 
 if(itn == 1) then
   open(44, file="SVC_phizeropote.csv")
   do i = 1, N
    write(44, 86) UU(i, :)
   enddo !s
   close(44)
 endif
 
 do i = 1, N
  do s = 1, kind
   !静電
   UU(i, s) = UU(i, s)+ch(s)*pe(i)
  enddo !s
 enddo !i
 
 return
 
 86 format(1PE25.15E3, 4(',', 1PE25.15E3))
end subroutine EP


!静電ポテンシャル差分
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
  if(abs(pep(1, i)) >= 1.d1) dp(i) = 1.d-3*DD
  if(abs(pep(1, i)) < 1.d1) dp(i) = 1.d-4*pep(1, i)*DD
!  if(i >= LATL) dp(i) = 1.d-3*DD
!  if(i == HATL) dp(i) = 5.d1
 enddo !i
 
 do i = 1, N
  pep(2, i) = pep(1, i) + dp(i)
  pep(3, i) = pep(1, i) - dp(i)
 enddo !i
 
 return
 
end subroutine pepm


!エネルギー(位置+Kperp)
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
 integer, dimension(3, N, kind, Z) :: lima
 integer :: ud, i, s, j, k
 real*8 :: CC
 
 lima = 0
 
 !amin
 do ud = 1, 3
  do i = 1, N
   do s = 1, kind
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
   do s = 1, kind
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
  do s = 1, kind
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
   do s = 1, kind
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
 
 open(70, file="SVC_a3.csv")
 do i = 1, N
  write(70, 62) a3(1, i, 1, 1, :)
 enddo !i
 close(70)
 
 return
 
62 format(1PE25.15E3, 2(',', 1PE25.15E3))
 
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


!電荷密度
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
   cvg(ud, i) = sqrt(cvg(ud, i))
  enddo !i
 enddo !ud
 
 cvn = sqrt(cvn/dble(N))
 
 return
 
end subroutine CV


!Newton法(npe)
subroutine Newtonphi(LATL, HATL, N, pep, cvg, npe, c, m, e, cvn, itn, check, DD, XX)
 implicit none
 integer, intent(in) :: N, itn
 real*8, dimension(3, N), intent(in) :: pep, cvg
 real*8, dimension(N), intent(out) :: npe
 real*8, intent(in) :: c, m, e
 integer, intent(inout) :: LATL, HATL, XX
 integer, dimension(LATL) :: kyoku
 real*8, intent(out) :: DD
 real*8 :: cvn, A
 real*8, save :: CL, DL, LL
 
 integer, intent(inout) :: check
 integer :: i, j, mxn, kyokuchi, serch, che
 real*8 :: CC1, CC2, CC3, mx, ser
 
 !学習率設定
 if(itn == 1) then
   CL = cvn
   LL = 4.d1
   check = 0
   che = 0
 endif
 !if(LL > 8.d1) LL = 1.d0
 if(CL > cvn .and. itn /= 1 .and. ((cvn < 1.d0 .and. XX == 1) .or. check /= 0)) then
   CL = cvn !最小のcvnを記憶
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
 
 !LLの下限
 if(che > 50 .or. check > 5 .or. (CL*3.d0 < cvn .and. XX ==1)) then
   DL = cvn !ひとつ前のcvnを記憶
   open(33, file = "SVC_min_E_L=4.csv")
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
 
 DL = cvn !ひとつ前のcvnを記憶
 
 DD = 1.1**(-LL+1.d0) !学習率
 XX = 0 !学習率起動リセット
 
 CC1 = 0.d0
 CC2 = 0.d0
 npe = pep(1, :)
 
 if(cvn > 1.d0) cvn = 1.d0
 
 do i = 2, HATL-1
  if(cvn < cvg(1, i)) then           !((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. 
    if(abs(cvg(2, i)-cvg(3, i)) > 1.d-12) then
      if(CC1 < abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))) then
        CC1 = abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))
      endif
    endif
  endif
 enddo !i
 if(LATL > HATL) then
  do i = HATL, LATL-1
   if(cvn < cvg(1, i)) then                !((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. 
     if(CC2 < abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))) then
       if(abs(cvg(2, i)-cvg(3, i)) > 1.d-12) then
         CC2 = abs(((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i))
       endif
     endif
   endif
  enddo !i
 endif
 do i = LATL, N-1
  if(cvn < cvg(1, i)) then              !((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. 
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
 !do i = HATL, LATL-1
  
  !変化の調整
  if(((cvg(1, i) > cvg(3, i)) .or. (cvg(1, i) > cvg(2, i))) .or. (cvn < cvg(1, i))) then                       !
  if(abs(cvg(2, i)-cvg(3, i)) > 1.d-12) then
  
!CC1の調整
  if(i < HATL .and. CC1 > 1.d0 .and. sqrt(cvn) > 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)/CC1*1.d0*DD
  if(i < HATL .and. CC1 > 1.d0 .and. sqrt(cvn) <= 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)/CC1*sqrt(cvn)*DD
  if(i < HATL .and. CC1 <= 1.d0 .and. sqrt(cvn) > 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)*1.d0*DD
  if(i < HATL .and. CC1 <= 1.d0 .and. sqrt(cvn) <= 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)*sqrt(cvn)*DD
  
!CC2の調整
  if((LATL > HATL) .and. (i >= HATL .and. i < LATL) .and. CC2 > 1.d0) then
    npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)/CC2*DD*sqrt(CC2)
  endif
  if((LATL > HATL) .and. (i >= HATL .and. i < LATL) .and. CC2 <= 1.d0) then
    XX = 1 !学習率起動
    npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)*DD!*sqrt(CC2)
  endif
  
!CC3の調整
  if(i >= LATL .and. CC3 > 1.d0 .and. sqrt(cvn) > 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)/CC3*1.d0*DD
  if(i >= LATL .and. CC3 > 1.d0 .and. sqrt(cvn) <= 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)/CC3*sqrt(cvn)*DD
  if(i >= LATL .and. CC3 <= 1.d0 .and. sqrt(cvn) > 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)*1.d0*DD
  if(i >= LATL .and. CC3 <= 1.d0 .and. sqrt(cvn) <= 1.d0) npe(i) = npe(i) - &
&((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0*cvg(1, i)*sqrt(cvn)*DD
  
  endif
  
  if(abs(cvg(2, i)-cvg(3, i)) <= 1.d-12) then
    A = ((pep(2, i)-pep(1, i))/(cvg(2, i)-cvg(1, i))+(pep(3, i)-pep(1, i))/(cvg(3, i)-cvg(1, i)))/2.d0
    if(isnan(npe(i) - pep(1, i)*1.d-2*A/abs(A))) then
      npe(i) = pep(1, i)
     else
      npe(i) = npe(i) - pep(1, i)*1.d-2*A/abs(A)
    endif
  endif
  
   else if((LATL > HATL) .and. (i >= HATL .and. i < LATL)) then
     XX = 1 !学習率起動
  endif
 enddo !i
 
 
 !HATLまで、極小値は一つしかないと仮定
 !kyoku = 0
 !kyokuchi = 0
 !do i = 2, HATL-1
 ! if((npe(i-1)-npe(i))*(npe(i)-npe(i+1)) < 0.d0) then
 !   kyoku(i) = 1
 !   if(kyokuchi == 0) kyokuchi = i
 ! endif
 !enddo !i
! 
! do i = 2, HATL-1
!  if(i /= HATL-1 .and. kyoku(i) == 1 .and. i /= kyokuchi) then
!    npe(i) = (npe(i-1)+npe(i+1))/2.d0
!  endif
!  if(i == HATL-1 .and. kyoku(i) == 1 .and. i /= kyokuchi) then
!    npe(i) = - npe(kyokuchi) - 2.d0*npe(i-1)
!  endif
! enddo !i
! 
 !HATLまで、常に静電ポテンシャルは電離圏に向かって減少すると仮定
 do j = 1, 10000
  serch = 0
  do i = 1, HATL-1
   if(npe(i-1) < npe(i)) then
     serch = 1
     if(i == HATL-1) then
       npe(i) = 2.d0*npe(i-2)-npe(i-1)
      else if(i == 2) then
       npe(i) = - npe(i)
      else
       ser = npe(i-1)
       npe(i-1) = npe(i)
       npe(i) = ser
     endif
   endif
  enddo !i
  if(serch == 0) then
    goto 100
  endif
 enddo
 
 
 
 
! 
! !HATL以降、常に静電ポテンシャルは電離圏に向かって増加すると仮定
 do j = 1, 10000
  serch = 0
  do i = HATL, N-1
   if(npe(i-1) > npe(i)) then
     serch = 1
     if(i == HATL .or. i == LATL) then
       npe(i) = 2.d0*npe(i+1)-npe(i+2)
      else if(i == HATL-1 .or. i== LATL-1) then
       npe(i) = 2.d0*npe(i-1)-npe(i-2)
      else if(i == 2) then
       npe(i) = npe(i)
      else
       ser = npe(i-1)
       npe(i-1) = npe(i)
       npe(i) = ser
     endif
   endif
  enddo !i
  if(serch == 0) then
    goto 100
  endif
 enddo
 
 100 print *, " "
 
 !HATL,LATLの補正
 if(HATL < LATL) then
  if(npe(HATL) < ((npe(HATL+1)+npe(LATL-1))/2.d0+npe(HATL-1))/2.d0 .and. npe(HATL) < 2.d0*npe(HATL-1)-1.d0*npe(HATL-2)) then
    HATL = HATL + 1
    print *, "HATL is changed to ", HATL
  endif
  if(npe(HATL-1) > ((npe(HATL)+npe(LATL-1))/2.d0+npe(HATL-2))/2.d0 .and. npe(HATL-1) > 2.d0*npe(HATL-2)-1.d0*npe(HATL-3)) then
    HATL = HATL - 1
    print *, "HATL is changed to ", HATL
  endif
  if(npe(LATL) < ((npe(LATL+1)+npe(N))/2.d0+npe(LATL-1))/2.d0 .and. npe(LATL) < 2.d0*npe(LATL+1)-1.d0*npe(LATL+2)) then
    LATL = LATL + 1
    print *, "LATL is changed to ", LATL
  endif
  if(npe(LATL-1) > ((npe(LATL)+npe(N))/2.d0+npe(LATL-2))/2.d0 .and. npe(LATL-1) > 2.d0*npe(LATL)-1.d0*npe(LATL+1)) then
    LATL = LATL - 1
    print *, "LATL is changed to ", LATL
  endif
 endif
 
 if(HATL == LATL) then
   if(npe(LATL) < 3.d0*npe(HATL+1)-2.d0*npe(HATL+2)) then
     HATL = HATL + 1
   endif
   if(npe(LATL-1) > 3.d0*npe(LATL-2)-2.d0*(LATL-3)) then
     LATL = LATL - 1
   endif
 endif
 
 
 open(70, file="SVC_npe.csv")
 do i = 2, N-1
  write(70, 52) (pep(2, i)-pep(3, i)), (cvg(2, i)-cvg(3, i)), cvg(1, i), cvg(2, i), cvg(3, i),&
&npe(i)-pep(1, i), npe(i)
 enddo !i
 close(70)
 
 
 print *, npe(HATL-1)-pep(1, HATL-1), npe(HATL)-pep(1, HATL), LL, DD, check
 print *, "  "
 return
 
 52 format(1PE25.15E3, 6(',', 1PE25.15E3))
 
end subroutine Newtonphi


!sigma差分作成
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

!Newton法(nsig)
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
   open(33, file = "SVC_min_E_L=4.csv")
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
&*sqrt(cvn)*DD
   if(abs(CC) < sigma(iodo)*1.d-1) nsig(iodo) = sigma(iodo) - CC
   if(abs(CC) >= sigma(iodo)*1.d-1) nsig(iodo) = sigma(iodo) - CC/abs(CC)*sigma(iodo)*1.d-1
 endif
 if((cvg(1, 1) > cvg(3, 1)) .or. (cvg(1, 1) > cvg(2, 1))) then
   CC = ((sigx(2, 2) - sigx(1, 2))/(cvg(2, 1) - cvg(1, 1))+(sigx(3, 2) - sigx(1, 2))/(cvg(3, 1) - cvg(1, 1)))/2.d0*cvg(1, 1)&
&*sqrt(cvn)*DD
   if(abs(CC) < sigma(mado)*1.d-1) nsig(mado) = sigma(mado) - CC
   if(abs(CC) >= sigma(mado)*1.d-1) nsig(mado) = sigma(mado) - CC/abs(CC)*sigma(mado)*1.d-1
 endif
 
 open(44, file="SVC_nsig_E_L=4.csv")
 do s = 1, kind
  write(44, 43) nsig(s), nsig(s)-sigma(s)
 enddo !s
 close(44)
 
 200 print *, "  "
 
 return
 
 43 format(1PE25.15E3, ',', 1PE25.15E3)
 
end subroutine Newtonsig


!Alfven速度
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
   if(ch(s) > 0) CC = CC + mass(s)*num(1, i, s)
  enddo !s
  Va(i) = BB(i)/sqrt(mu0*CC)
  Vr(i) = Va(i)/sqrt(1.d0 + (Va(i)/c)**2.d0)
 enddo !i
 
 open(99, file="SVC_AW_E_L=4.csv")
 do i = 1, N
  write(99, 98) Va(i), Vr(i), Vr(i)/c
 enddo !i
 close(99)
 
 return
 
 98 format(1PE25.15E3, 2(',', 1PE25.15E3))
 
end subroutine Alf

