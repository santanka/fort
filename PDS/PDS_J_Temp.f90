program SVC_Temperature

implicit none

!全体を通しての文字
integer :: i, s, j   !do構文用
integer, parameter :: N = 94 !grid数
integer, parameter :: Z = 300 !速度空間grid数
integer, parameter :: kind = 8 !粒子種数
integer, parameter :: spece = 5 !真・粒子種数
real*8, parameter :: pi = 4.d0*atan(1.d0) !円周率
character(len=80) :: dummy !使用しない文字用

!mag_FA関連
real*8, dimension(N) :: lambda, BB, dteq

!SVC_parameter関連
real*8 :: ee, cc, mu0, GG, alpha, ME, MassE, RE, omE, Req, MassI, ep0

!SVC_IC関連
real*8, dimension(N) :: ss, sR
real*8, dimension(N-1) :: ds !grid間

!SVC_BC関連
character(len=1), dimension(kind) :: so !起源
real*8, dimension(kind) :: Tp, Ta, ch, mass
integer, dimension(kind) :: kk !粒子種
integer, dimension(kind) :: sigmab !起源のgrid

!SVC_min関連
real*8, dimension(N) :: pe !静電ポテンシャル
real*8, dimension(kind) :: sigma !境界での数密度

!subroutine連絡用
real*8, dimension(kind, Z) :: mu !断熱不変量
real*8, dimension(N, kind) :: UU !位置エネルギー
real*8, dimension(N, kind, Z, 3) :: a3 !accessibility(1:min, 2:lim, 3:max)
real*8 :: ff !分布関数(function)
real*8, dimension(N, kind, Z) :: theta1 !a積分(数密度)
real*8, dimension(N, kind) :: num !数密度
real*8, dimension(N) :: rhov !電荷密度
real*8, dimension(N) :: rhop !Poisson方程式の電荷密度
real*8, dimension(N) :: cvg !収束値
real*8, dimension(N, 3) :: VA !Alfven速度
real*8, dimension(N, kind, Z) :: theta2 !a積分(平均流速)
real*8, dimension(N, kind) :: aave !平均流速(a)
real*8, dimension(N, kind, Z, 2) :: theta3 !a積分(温度)
real*8, dimension(N, kind) :: Tperp, Tpara, TTT !温度[eV]
real*8, dimension(N, kind) :: Pperp, Ppara, PP !圧力
real*8, dimension(N, spece) :: taave !真・平均流速
real*8, dimension(N, kind, Z, 2) :: theta4 !a積分(真・温度)
real*8, dimension(N, spece) :: tTperp, tTpara, tT !真・温度
real*8, dimension(N, kind) :: bs !各β値
real*8, dimension(N, spece) :: tbs !各β値
real*8, dimension(N) :: beta, tbeta !β値
real*8, dimension(N) :: Tall, tTall !温度
real*8, dimension(N, spece) :: tPperp, tPpara, tPP !真・圧力
real*8, dimension(N) :: Pall, tPall !総プラズマ圧



!mag_FAの抽出
open(50, file="mag_FA_J_kai.csv", action="read", status="old")
do i = 1, N
 read(50, *) lambda(i), BB(i), dteq(i)
enddo !i
close(50)


!SVC_parameterの抽出
open(51, file="SVC1_parameter.csv", action="read", status="old")
read(51, *) dummy, ee
read(51, *) dummy, cc
read(51, *) 
read(51, *) dummy, mu0
read(51, *) dummy, GG
read(51, *) dummy, alpha
read(51, *) dummy, ME
read(51, *) dummy, MassE
read(51, *) dummy, RE
read(51, *) dummy, omE
read(51, *) dummy, Req
read(51, *) dummy, MassI
close(51)

ep0 = 1.d0/mu0/cc**2.d0


!SVC_ICの抽出
open(52, file="PDS_IC_J_30kV_kai.csv", action="read", status="old")
read(52, *) !1行とばし
do i = 1, N
 read(52, *) dummy, dummy, ss(i), sR(i)
enddo !i
close(52)

do i = 1, N-1
 ds(i) = ss(i+1) - ss(i) !grid間距離
enddo !i


!SVC_BCの抽出
open(53, file="SVC_BC.csv", action="read", status="old")
read(53, *) !1行とばし
do i = 1, kind
 read(53, *) dummy, dummy, Tp(i), Ta(i), so(i), ch(i), mass(i), kk(i)
enddo !i
close(53)

Tp = Tp*ee
ch = ch*ee

do i = 1, kind
 if(so(i) == "I") sigmab(i) = N
 if(so(i) == "M") sigmab(i) = 1
enddo !i


!SVC_minの抽出
open(54, file="PDS_J_30kV_min_kai.csv", action="read", status="old")
read(54, *) !1行とばし
do i = 1, N
 read(54, *) pe(i)
enddo !i

do i = 1, kind
 read(54, *) sigma(i)
enddo !i
close(54)


!断熱不変量
call AI(N, Z, kind, alpha, Tp, BB, sigmab, mu)
print *, "Finish AI"

!位置エネルギー
call EP(N, kind, GG, MassE, Req, omE, MassI, pe, lambda, dteq, ch, mass, UU)
print *, "Finish EP"

!accessibility(1:min, 2:lim, 3:max)
call access(N, Z, kind, alpha, Tp, Ta, mass, BB, UU, mu, sigmab, a3)
print *, "Finish access"

!aについて積分(数密度)
call TH1(N, kind, Z, Tp, Ta, a3, theta1)
print *, "Finish TH1"

!muについて積分(数密度)
call NM1(kind, N, Z, Tp, sigma, BB, sigmab, mu, theta1, num)
print *, "Finish NM1"

!電荷密度
call RR(N, kind, ch, num, rhov)
print *, "Finish RR"

!Poisson方程式
call Poisson(N, ep0, pe, ds, rhop)
print *, "Finish Poisson"

!収束チェック
call CV(N, kind, rhov, rhop, ch, num, cvg)
print *, "Finish CV"

!Alfven速度
call Alf(N, kind, mu0, cc, BB, ch, mass, num, VA)
print *, "Finish Alf"

!aについて積分(平均流速)
call TH2(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta2)
print *, "Finish TH2"

!muについて積分(平均流速)
call NM2(N, kind, Z, BB, mass, mu, theta2, num, sigmab, aave)
print *, "Finish NM2"

!aについて積分(温度)
call TH3(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, aave, a3, theta3)
print *, "Finish TH3"

!muについて積分(温度)
call NM3(N, kind, Z, ee, mu0, BB, mass, mu, theta3, num, sigmab, Tperp, Tpara, TTT, Tall, Pperp, Ppara, PP, Pall, bs, beta)
print *, "Finish NM3"

!粒子種ごとの流速
call ryusoku(N, kind, spece, aave, num, UU, kk, sigmab, taave)
print *, "Finish ryusoku"

!粒子種ごとの温度(a積分)
call TH4(N, kind, Z, spece, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, kk, taave, a3, theta4)
print *, "Finish TH4"

!粒子種ごとの温度(mu積分)
call NM4(N, kind, Z, spece, ee, mu0, BB, mass, mu, theta4, num, UU, kk, sigmab, tPperp, &
& tPpara, tPP, tPall, tTperp, tTpara, tT, tbs, tbeta, tTall)
print *, "Finish NM4"


!ファイル化
open(60, file="PDS_J_30kV_all_kai.csv")
do i = 1, N
 write(60, 70) ss(i), sR(i), lambda(i), dteq(i), BB(i), pe(i), UU(i, :), num(i, :), rhov(i), rhop(i), cvg(i), &
  & VA(i, :), aave(i, :), Pperp(i, :), Ppara(i, :), PP(i, :), Pall(i), Tperp(i, :), Tpara(i, :), TTT(i, :), Tall(i), &
  & bs(i, :), beta(i), taave(i, :), tPperp(i, :), tPpara(i, :), tPP(i, :), tPall(i), tTperp(i, :), tTpara(i, :), &
  & tT(i, :), tTall(i), tbs(i, :), tbeta(i)
enddo !i
close(60)
print *, "Finish PDS_J_30kV_all_kai.csv"

open(61, file="PDS_J_30kV_a3_all_kai.csv")
do s = 1, kind
 do j = 1, Z
  do i = 1, N
   write(61, 72) s, j, i, a3(i, s, j, :), mu(s, j)
  enddo !i
 enddo !j
enddo !s
close(61)
print *, "Finish PDS_J_30kV_a3_all_kai.csv"


!call disfun(N, kind, Z, sigmab, ee, BB, Tp, Ta, mass, sigma, mu, UU, a3)

call disfun2(N, kind, Z, sigmab, ee, BB, Tp, Ta, mass, sigma, mu, UU, a3)


70 format(1PE25.15E3, 137(',', 1PE25.15E3)) !10*kind+8*spece+18
72 format(3(I3, ','), 1PE25.15E3, 3(',', 1PE25.15E3))

end program SVC_Temperature





!subroutine & function

!断熱不変量
subroutine AI(N, Z, kind, alpha, Tp, BB, sigmab, mu)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp
 real*8, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(kind, Z), intent(out) :: mu
 
 integer :: j, s
 
 
 do s = 1, kind
  do j = 1, Z
   mu(s, j) = 1.d-30*(alpha*Tp(s)/BB(sigmab(s))/1.d-30)**(dble(j-1)/dble(Z-1))
  enddo !j
 enddo !s
 
 return
 
end subroutine AI


!位置エネルギー
subroutine EP(N, kind, GG, MassE, Req, omE, MassI, pe, lambda, dteq, ch, mass, UU)
 implicit none
 integer, intent(in) :: N, kind
 real*8, intent(in) :: GG, MassE, Req, omE, MassI
 real*8, dimension(N), intent(in) :: pe, lambda, dteq
 real*8, dimension(kind), intent(in) :: ch, mass
 
 real*8, dimension(N, kind), intent(out) :: UU
 
 integer :: i, s
 
 
 do i = 1, N
  do s = 1, kind
   !静電ポテンシャル
   UU(i, s) = ch(s)*pe(i)
   !惑星重力
   UU(i, s) = UU(i, s) - GG*MassE*mass(s)/Req/cos(lambda(i))**2.d0
   !惑星遠心力
   UU(i, s) = UU(i, s) - 5.d-1*mass(s)*omE**2.d0*(Req*cos(lambda(i))**3.d0)**2.d0
   !衛星重力
   if(MassI /= 0.d0) then
     UU(i, s) = UU(i, s) - GG*MassI*mass(s)/dteq(i)
   endif
  enddo !s
 enddo !i
 
 return
 
end subroutine EP


!accessibility
subroutine access(N, Z, kind, alpha, Tp, Ta, mass, BB, UU, mu, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp, Ta, mass
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(N, kind, Z, 3), intent(out) :: a3
 
 real*8, dimension(N, kind, Z) :: UB
 integer :: i, s, j, k
 integer :: t !フリー
 real*8 :: f !フリー
 
 
 !UB作成
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    UB(i, s, j) = UU(i, s) + mu(s, j)*BB(i)
   enddo !j
  enddo !s
 enddo !i
 
 
 !amin
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"起源
    if(sigmab(s) == 1) then
      t = 1
      do k = 1, i
       if(UB(k, s, j) > UB(t, s, j)) t = k
      enddo !k
      a3(i, s, j, 1) = sqrt(UB(t, s, j) - UB(sigmab(s), s, j))
    endif
    
    !"I"起源
    if(sigmab(s) == N) then
      t = N
      do k = i, N
       if(UB(k, s, j) > UB(t, s, j)) t = k
      enddo !k
      a3(i, s, j, 1) = sqrt(UB(t, s, j) - UB(sigmab(s), s, j))
    endif
    
    !aminの値調整
    if(a3(i, s, j, 1) == 0.d0) a3(i, s, j, 1) = sqrt(mass(s)/2.d0)*1.d-3
    
   enddo !j
  enddo !s
 enddo !i
 
 !alim
 do s = 1, kind
  a3(:, s, :, 2) = sqrt(alpha*Tp(s)*Ta(s))
 enddo !s
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
!    if(a3(i, s, j, 2) < a3(i, s, j, 1)*1.d2) a3(i, s, j, 2) = a3(i, s, j, 1)*1.d2
   enddo !j
  enddo !s
 enddo !i
 
 !amax
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"起源
    if(sigmab(s) == 1) then
      if(i == N) a3(i, s, j, 3) = 0.d0
      if(i /= N) then
        t = i+1
        do k = i+1, N
         if(UB(k, s, j) > UB(t, s, j)) t = k
        enddo !k
        a3(i, s, j, 3) = UB(t, s, j) - UB(sigmab(s), s, j)
        if(a3(i, s, j, 3) < 0.d0) a3(i, s, j, 3) = 0.d0
        if(a3(i, s, j, 3) >= 0.d0) a3(i, s, j, 3) = sqrt(a3(i, s, j, 3))
      endif
    endif
    
    !"I"起源
    if(sigmab(s) == N) then
      if(i == 1) a3(i, s, j, 3) = 0.d0
      if(i /= 1) then
        t = 1
        do k = 1, i-1
         if(UB(k, s, j) > UB(t, s, j)) t = k
        enddo !k
        a3(i, s, j, 3) = UB(t, s, j) - UB(sigmab(s), s, j)
        if(a3(i, s, j, 3) < 0.d0) a3(i, s, j, 3) = 0.d0
        if(a3(i, s, j, 3) >= 0.d0) a3(i, s, j, 3) = sqrt(a3(i, s, j, 3))
      endif
    endif
    
   enddo !j
  enddo !s
 enddo !i
 
 return
 
end subroutine access


!分布関数
real*8 function ff(sigma, BB, B1, Tp, Ta, UU, U1, mu, aa, mass)
 implicit none
 real*8, intent(in) :: sigma, BB, B1, Tp, Ta, UU, U1, mu, aa, mass
 
 real*8, parameter :: pi = 4.0*atan(1.d0)
 
 
 ff = sigma*sqrt((mass/2.d0/pi/Tp)**3.d0/Ta)*exp(-B1*mu/Tp)*exp(-aa**2.d0/Ta/Tp)
 return
 
end function ff



!aについて積分(数密度)
subroutine TH1(N, kind, Z, Tp, Ta, a3, theta1)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(kind), intent(in) :: Tp, Ta
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z), intent(out) :: theta1
 
 integer :: i, s, j
 real*8 :: thetaL, thetaM
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    thetaL = 0.d0
    thetaM = 0.d0
    
    !emit
    thetaL = 1-erf(a3(i, s, j, 1)/sqrt(Tp(s)*Ta(s)))
    
    !bounce
    if(a3(i, s, j, 3) > a3(i, s, j, 1)) then
      thetaM = erf(a3(i, s, j, 3)/sqrt(Tp(s)*Ta(s))) - erf(a3(i, s, j, 1)/sqrt(Tp(s)*Ta(s)))
    endif
    
    theta1(i, s, j) = thetaL + thetaM
   enddo !j
  enddo !s
 enddo !i 
 
 return
 
end subroutine TH1


!muについて積分(数密度)
subroutine NM1(kind, N, Z, Tp, sigma, BB, sigmab, mu, theta1, num)
 implicit none
 integer, intent(in) :: kind, N, Z
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(kind), intent(in) :: Tp, sigma
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z), intent(in) :: theta1
 real*8, dimension(N, kind), intent(out) :: num
 real*8, dimension(N, kind, Z) :: th
 integer :: i, s, j, numx
 
 do s = 1, kind
  do j = 1, Z
   th(:, s, j) = theta1(:, s, j)*exp(-mu(s, j)*BB(sigmab(s))/Tp(s))
  enddo !j
 enddo !s
 
 num = 0.d0
 
 do s = 1, kind
  do j = 1, Z-1
   num(:, s) = num(:, s) + (th(:, s, j)+th(:, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
  enddo !j
  
  num(:, s) = sigma(s)*num(:, s)*BB(sigmab(s))/2.d0/Tp(s)
  
 enddo !s
 
 return
 
end subroutine NM1


!電荷密度
subroutine RR(N, kind, ch, num, rhov)
 implicit none
 integer, intent(in) :: N, kind
 real*8, dimension(kind) :: ch
 real*8, dimension(N, kind) :: num
 
 real*8, dimension(N), intent(out) :: rhov
 
 integer :: i, s
 
 
 rhov = 0.d0
 
 do i = 1, N
  do s = 1, kind
   rhov(i) = rhov(i) + ch(s) * num(i, s)
  enddo !s
 enddo !i
 
 return
 
end subroutine RR


!Poisson方程式
subroutine Poisson(N, ep0, pe, ds, rhop)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: ep0
 real*8, dimension(N), intent(in) :: pe
 real*8, dimension(N-1), intent(in) :: ds
 
 real*8, dimension(N), intent(out) :: rhop
 
 integer :: i
 real*8 :: p
 
 
 rhop = 0.d0
 
 do i = 2, N-1
  p = 2.d0/ds(i)/(ds(i)+ds(i-1))*pe(i+1)
  p = p + 2.d0/ds(i-1)/(ds(i)+ds(i-1))*pe(i-1)
  p = p - 2.d0/ds(i)/ds(i-1)*pe(i)
  rhop(i) = -ep0 * p
 enddo !i
 
 return
 
end subroutine Poisson


!収束チェック
subroutine CV(N, kind, rhov, rhop, ch, num, cvg)
 implicit none
 integer, intent(in) :: N, kind
 real*8, dimension(N), intent(in) :: rhov, rhop
 real*8, dimension(kind), intent(in) :: ch
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N), intent(out) :: cvg
 
 integer :: i, s
 real*8, dimension(N) :: nume, numi
 
 
 cvg = 0.d0
 nume = 0.d0
 numi = 0.d0
 
 do i = 1, N
  do s = 1, kind
   if(ch(s) < 0) nume(i) = nume(i) + abs(ch(s)*num(i, s))
   if(ch(s) > 0) numi(i) = numi(i) + abs(ch(s)*num(i, s))
  enddo !s
 enddo !i
 
 do i = 1, N
  cvg(i) = (rhov(i)-rhop(i))**2.d0/nume(i)/numi(i)
 enddo !i
 
 return
 
end subroutine CV


!Alfven速度
subroutine Alf(N, kind, mu0, cc, BB, ch, mass, num, VA)
 implicit none
 integer, intent(in) :: N, kind
 real*8, intent(in) :: mu0, cc
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: ch, mass
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N, 3), intent(out) :: VA
 
 integer :: i, s
 real*8 :: w
 
 
 do i = 1, N
  w = 0.d0
  do s = 1, kind
   if(ch(s) > 0) w = w + mass(s)*num(i, s)
  enddo !s
  VA(i, 1) = BB(i)/sqrt(mu0 * w)
  VA(i, 2) = VA(i, 1)/sqrt(1.d0 + (VA(i, 1)/cc)**2.d0)
  VA(i, 3) = VA(i, 2)/cc
 enddo !i
 
 return
 
end subroutine Alf


!aについて積分(平均流速)
subroutine TH2(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta2)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(N, kind, Z, 3), intent(inout) :: a3
 
 real*8, dimension(N, kind, Z), intent(out) :: theta2
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !分布関数(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !aminの値調整
    if(a3(i, s, j, 1) == 0.d0) a3(i, s, j, 1) = sqrt(mass(s)/2.d0)*1.d-3
    
    !Vpara作成
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0 > 0.d0) then
         thetaL = thetaL + &
& ( aL(p)/abs(aL(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0) * &
& ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
& aL(p+1)/abs(aL(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p+1)**2.d0) * &
&ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
& / 2.d0 * abs(aL(p+1) - aL(p))
       endif
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0 > 0.d0) then
         thetaM = thetaM + &
& ( aM(p)/abs(aM(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0) * &
& ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
& aM(p+1)/abs(aM(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p+1)**2.d0) * &
& ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
& / 2.d0 * abs(aM(p+1) - aM(p))
       endif
      enddo !p
    endif
    
    theta2(i, s, j) = thetaL + thetaM
    
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH2


!muについて積分(平均流速)
subroutine NM2(N, kind, Z, BB, mass, mu, theta2, num, sigmab, aave)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z), intent(in) :: theta2
 real*8, dimension(N, kind), intent(in) :: num
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(N, kind), intent(out) :: aave
 
 integer :: i, s, j
 real*8 :: nnn, a
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 
 a = -1.d0 !NaNの前処理
 
 do i = 1, N
  do s = 1, kind
   
   !Vpara作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta2(i, s, j) + theta2(i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   if(num(i, s) /= 0.d0) aave(i, s) = (2.d0/mass(s))**(3.d0/2.d0)*pi*BB(sigmab(s))/num(i, s)*nnn
   if(num(i, s) == 0.d0) aave(i, s) =  sqrt(a) !NaNを代入
   
  enddo !s
 enddo !i
 
 
 return
 
end subroutine NM2


!!!コメント、NaNは出るもの
!!!問題点→平均流速、温度は粒子種ごとに出す必要がありそう、まだ未実装



!aについて積分(温度)
subroutine TH3(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, aave, a3, theta3)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(N, kind), intent(in) :: aave
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z, 2), intent(out) :: theta3
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM, a
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !分布関数(function)
 
 
 a = -1.d0 !NaN前処理
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !Tperp作成
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       thetaL = thetaL + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
        & / 2.d0 * abs(aL(p+1) - aL(p))
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       thetaM = thetaM + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
        & / 2.d0 * abs(aM(p+1) - aM(p))
      enddo !p
    endif
    
    theta3(i, s, j, 1) = mu(s, j)*(thetaL + thetaM)
    
    
    !Tpara作成
    thetaL = 0.d0
    thetaM = 0.d0
    
    if(aave(i, s) == aave(i, s)) then
      
      !alim
      if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
        do p = 1, Z/2
         aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
        enddo !p
        
        do p = 1, Z/2-1
         if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0 > 0.d0) then
           thetaL = thetaL + &
& ( (aL(p)/abs(aL(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
& (aL(p+1)/abs(aL(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p+1)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
& / 2.d0 * abs(aL(p+1) - aL(p))
         endif
        enddo !p
      endif
      
      !amax
      if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
        do p = 1, Z/2
         aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
        enddo !p
        
        do p = 1, Z/2-1
         if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0 > 0.d0) then
           thetaM = thetaM + &
& ( (aM(p)/abs(aM(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
& (aM(p+1)/abs(aM(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p+1)**2.d0)-aave(i, s))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
& / 2.d0 * abs(aM(p+1) - aM(p))
         endif
        enddo !p
      endif
      
      theta3(i, s, j, 2) = thetaL + thetaM
      
     else if(aave(i, s) == aave(i, s)) then
      theta3(i, s, j, 2) = sqrt(a)
    endif
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH3


!muについて積分(温度)
subroutine NM3(N, kind, Z, ee, mu0, BB, mass, mu, theta3, num, sigmab, Tperp, Tpara, TTT, Tall, Pperp, Ppara, PP, Pall, bs, beta)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, intent(in) :: ee, mu0
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z, 2), intent(in) :: theta3
 real*8, dimension(N, kind), intent(in) :: num
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(N, kind), intent(out) :: Tperp, Tpara, Pperp, Ppara, TTT, PP, bs
 real*8, dimension(N), intent(out) :: Pall, beta, Tall
 
 integer :: i, s, j
 real*8 :: nnn, a
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8, dimension(N) :: nm
 
 
 a = -1.d0 !NaN前処理
 
 do i = 1, N
  do s = 1, kind
   
   !Tperp作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta3(i, s, j, 1) + theta3(i, s, j+1, 1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Pperp(i, s) = 2.d0*sqrt(2.d0/mass(s)**3.d0)*(BB(sigmab(s))**2.d0)*pi*nnn
   
   if(num(i, s) == num(i, s)) then
     Tperp(i, s) = Pperp(i, s)/num(i, s) /ee
    else if(num(i, s) /= num(i, s)) then
     Tperp(i, s) = sqrt(a)
   endif
   
   
   !Tpara作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta3(i, s, j, 2) + theta3(i, s, j+1, 2))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Ppara(i, s) = 4.d0*sqrt(2.d0/mass(s)**3.d0)*BB(sigmab(s))*pi*nnn
   
   if(num(i, s) == num(i, s)) then
     Tpara(i, s) = Ppara(i, s)/num(i, s) /ee
    else if(num(i, s) /= num(i, s)) then
     Tpara(i, s) = sqrt(a)
   endif
   
   if(Ppara(i, s) == Ppara(i, s)) PP(i, s) = (2.d0*Pperp(i, s) + Ppara(i, s))/3.d0
   if(Ppara(i, s) /= Ppara(i, s)) PP(i, s) = sqrt(a)
   
   if(num(i, s) == num(i, s) .and. PP(i, s) == PP(i, s)) then
     TTT(i, s) = PP(i, s)/num(i, s) /ee
    else
     TTT(i, s) = sqrt(a)
   endif
   
   if(PP(i, s) == PP(i, s)) Pall(i) = Pall(i) + PP(i, s)
   
   bs(i, s) = 2.d0*mu0*PP(i, s)/BB(i)**2.d0
   if(PP(i, s) == PP(i, s)) then
     bs(i, s) = 2.d0*mu0*PP(i, s)/BB(i)**2.d0
    else
     bs(i, s) = sqrt(a)
   endif
   if(bs(i, s) == bs(i, s)) beta(i) = beta(i) + bs(i, s)
   
   if(num(i, s) == num(i, s)) nm(i) = nm(i) + num(i, s)
   
  enddo !s
  
  Tall(i) = Pall(i)/nm(i) /ee
  
 enddo !i
 
 
 return
 
end subroutine NM3


subroutine ryusoku(N, kind, spece, aave, num, UU, kk, sigmab, taave)
 implicit none
 integer, intent(in) :: N, kind, spece
 real*8, dimension(N, kind), intent(in) :: aave, num, UU
 integer, dimension(kind), intent(in) :: kk, sigmab
 
 real*8, dimension(N, spece), intent(out) :: taave
 
 integer :: i, s
 real*8, dimension(N, spece) :: nnn
 real*8 :: a
 
 
 a = -1.d0 !NaNの前処理
 taave = 0.d0
 nnn = 0.d0
 
 
 do i = 1, N
  do s = 1, kind
   nnn(i, kk(s)) = nnn(i, kk(s)) + num(i, s)
   if(aave(i, s) == aave(i, s)) then
     !M->I向きを正とする
     if(sigmab(s) == N) taave(i, kk(s)) = taave(i, kk(s)) - &
      & num(i, s)*aave(i, s)
     if(sigmab(s) == 1) taave(i, kk(s)) = taave(i, kk(s)) + &
      & num(i, s)*aave(i, s)
   endif
  enddo !s
 enddo !i
 
 do i = 1, N
  do s = 1, spece
   if(nnn(i, s) /= 0.d0) taave(i, s) = taave(i, s)/nnn(i, s)
   if(nnn(i, s) == 0.d0) taave(i, s) = sqrt(a) !NaNを代入
  enddo !s
 enddo !i
 
 
 return
 
end subroutine ryusoku


!aについて積分(真・温度)
subroutine TH4(N, kind, Z, spece, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, kk, taave, a3, theta4)
 implicit none
 integer, intent(in) :: N, kind, Z, spece
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab, kk
 real*8, dimension(N, spece), intent(in) :: taave
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z, 2), intent(out) :: theta4
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !分布関数(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !Tperp作成
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       thetaL = thetaL + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
        & / 2.d0 * abs(aL(p+1) - aL(p))
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       aM(p) = -(a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      
      do p = 1, Z/2-1
       thetaM = thetaM + &
        & ( ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
        & ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
        & / 2.d0 * abs(aM(p+1) - aM(p))
      enddo !p
    endif
    
    theta4(i, s, j, 1) = mu(s, j)*(thetaL + thetaM)
    
    
    !Tpara作成
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) > a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
      do p = 1, Z/2
       if(sigmab(s) == 1) aL(p) = a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
       if(sigmab(s) == N) aL(p) = -a3(i, s, j, 1)*(a3(i, s, j, 2)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0 > 0.d0) then
         thetaL = thetaL + &
& ( (aL(p)/abs(aL(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p), mass(s)) + &
& (aL(p+1)/abs(aL(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aL(p+1)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aL(p+1), mass(s)) ) &
& / 2.d0 * abs(aL(p+1) - aL(p))
       endif
      enddo !p
    endif
    
    !amax
    if(a3(i, s, j, 3) > a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
      do p = 1, Z/2
       if(sigmab(s) == 1) aM(p) = -a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
       if(sigmab(s) == N) aM(p) = a3(i, s, j, 1)*(a3(i, s, j, 3)/a3(i, s, j, 1))**(dble(p-1)/dble(Z/2-1))
      enddo !p
      
      do p = 1, Z/2-1
       if(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0 > 0.d0) then
         thetaM = thetaM + &
& ( (aM(p)/abs(aM(p))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p), mass(s)) + &
& (aM(p+1)/abs(aM(p+1))*sqrt(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+aM(p+1)**2.d0)-taave(i, kk(s)))**2.d0 &
& *ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), aM(p+1), mass(s)) ) &
& / 2.d0 * abs(aM(p+1) - aM(p))
       endif
      enddo !p
    endif
    
    theta4(i, s, j, 2) = thetaL + thetaM
    
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH4


!muについて積分(温度)
subroutine NM4(N, kind, Z, spece, ee, mu0, BB, mass, mu, theta4, num, UU, kk, &
 & sigmab, tPperp, tPpara, tPP, tPall, tTperp, tTpara, tT, tbs, tbeta, tTall)
 implicit none
 integer, intent(in) :: N, kind, Z, spece
 real*8, intent(in) :: ee, mu0
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z, 2), intent(in) :: theta4
 real*8, dimension(N, kind), intent(in) :: num, UU
 integer, dimension(kind), intent(in) :: kk, sigmab
 
 real*8, dimension(N, spece), intent(out) :: tTperp, tTpara, tT, tbs, tPperp, tPpara, tPP
 real*8, dimension(N), intent(out) :: tbeta, tTall, tPall
 
 integer :: i, s, j
 real*8 :: nnn, a
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8, dimension(N, kind) :: Tperp, Tpara
 real*8, dimension(N, spece) :: nm
 real*8, dimension(N) :: nn
 
 
 a = -1.d0 !NaNの前処理
 
 do i = 1, N
  do s = 1, kind
   
   !Tperp作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta4(i, s, j, 1) + theta4(i, s, j+1, 1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   
   if(num(i, s) /= 0.d0) Tperp(i, s) = 2.d0*sqrt(2.d0/mass(s)**3.d0)*(BB(sigmab(s))**2.d0)*pi/num(i, s)*nnn
   if(num(i, s) == 0.d0) Tperp(i, s) = sqrt(a) !NaNを代入
   
   
   !Tpara作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta4(i, s, j, 2) + theta4(i, s, j+1, 2))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   if(num(i, s) /= 0.d0) Tpara(i, s) = 4.d0*sqrt(2.d0/mass(s)**3.d0)*BB(sigmab(s))*pi/num(i, s)*nnn
   if(num(i, s) == 0.d0) Tpara(i, s) = sqrt(a) !NaNを代入
   
  enddo !s
 enddo !i
 
 
 tPperp = 0.d0
 tPpara = 0.d0
 nm = 0.d0
 
 
 do i = 1, N
  do s = 1, kind
   nm(i, kk(s)) = nm(i, kk(s)) + num(i, s)
   if(Tperp(i, s) == Tperp(i, s)) then
     tPperp(i, kk(s)) = tPperp(i, kk(s)) + num(i, s)*Tperp(i, s)
   endif
   if(Tpara(i, s) == Tpara(i, s)) then
     tPpara(i, kk(s)) = tPpara(i, kk(s)) + num(i, s)*Tpara(i, s)
   endif
  enddo !s
 enddo !i
 
 tPP = (tPperp*2.d0+tPpara)/3.d0
 
 do i = 1, N
  do s = 1, spece
   if(nm(i, s) /= 0.d0) then
     tTperp(i, s) = tPperp(i, s)/nm(i, s) /ee
     tTpara(i, s) = tPpara(i, s)/nm(i, s) /ee
     tT(i, s) = tPP(i, s)/nm(i, s) /ee
    else if(nm(i, s) == 0.d0) then
     tTperp(i, s) = sqrt(a)
     tTpara(i, s) = sqrt(a)
     tT(i, s) = sqrt(a)
   endif
  enddo !s
 enddo !i
 
 do i = 1, N
  do s = 1, spece
   if(tPP(i, s) == tPP(i, s)) tPall(i) = tPall(i) + tPP(i, s)
  enddo !s
 enddo !i
 
 do s = 1, spece
  tbs(:, s) = 2.d0*mu0*tPP(:, s)/BB**2.d0
 enddo !s
 
 tbeta = 2.d0*mu0*tPall/BB**2.d0
 
 nn = 0.d0
 do s = 1, spece
  nn = nn + nm(:, s)
 enddo !s
 
 tTall = tPall/nn /ee
 
 
 return
 
end subroutine NM4


subroutine disfun(N, kind, Z, sigmab, ee, BB, Tp, Ta, mass, sigma, mu, UU, a3)
 implicit none
 integer, intent(in) :: N, kind, Z
 integer, dimension(kind), intent(in) :: sigmab
 real*8, intent(in) :: ee
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: Tp, Ta, mass
 real*8, dimension(kind), intent(in) ::sigma
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 character(len=128) :: filename, filename1, filename2
 integer :: i, j, o, p, q, s
 real*8 :: AA, vperpb, vparab, vperp, vpara, kperp, kpara, kall, CC
 
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8 :: ff !分布関数
 
 
 do s = 1, kind
  do i = 1, N
   if(i == 29) then !##
   !ファイル名作成
   write(filename1, *) s
   write(filename2, *) i
   filename = 'PDS_J_30kV_disfun_s='//trim(adjustl(filename1))//'_i='//trim(adjustl(filename2))//'.csv'
   print *, filename
   
   open(62, file=filename)
   
   do j = 1, Z
    if(mod(j, 3) == 1) then
      do p = 2, 3
       if(p == 2 .or. (p == 3 .and. a3(i, s, j, 3) > a3(i, s, j, 1))) then
         do q = 1, Z/2
          if(mod(q, 3) == 1) then
            if((p == 2 .and. sigmab(s) == 1) .or. (p == 3 .and. sigmab(s) == N)) then
              AA = a3(i, s, j, 1)*(a3(i, s, j, p)/a3(i, s, j, 1))**(dble(q-1)/dble(Z/2-1))
             else if((p == 2 .and. sigmab(s) == N) .or. (p == 3 .and. sigmab(s) == 1)) then
              AA = - a3(i, s, j, 1)*(a3(i, s, j, p)/a3(i, s, j, 1))**(dble(q-1)/dble(Z/2-1))
            endif
            
            CC = ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), AA, mass(s))
            
            vpara = sqrt(2.d0/mass(s))*AA/abs(AA)*sqrt(abs(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+AA**2.d0))
            vperp = sqrt(2.d0*mu(s, j)*BB(i)/mass(s))
            vparab = sqrt(2.d0/mass(s))*AA
            vperpb = sqrt(2.d0*mu(s, j)*BB(sigmab(s))/mass(s))
            kpara = mass(s)*vpara**2.d0 /2.d0 /ee
            kperp = mass(s)*vperp**2.d0 /2.d0 /ee
            kall = (kpara+2.d0*kperp)/3.d0
            
            write(62, 74) i, vperpb, vparab, vperp, vpara, kperp, kpara, kall, CC, 2*pi*vperp*CC
          endif
         enddo !q
       endif
      enddo !p
    endif
   enddo !j
   
   close(62)
   endif !##
  enddo !i
 enddo !s
 
 
 return
 
 74 format(I3, 9(',', 1PE25.15E3))
 
end subroutine disfun


subroutine disfun2(N, kind, Z, sigmab, ee, BB, Tp, Ta, mass, sigma, mu, UU, a3)
 implicit none
 integer, intent(in) :: N, kind, Z
 integer, dimension(kind), intent(in) :: sigmab
 real*8, intent(in) :: ee
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: Tp, Ta, mass
 real*8, dimension(kind), intent(in) ::sigma
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 character(len=128) :: filename, filename1, filename2
 integer :: i, j, o, p, q, s
 real*8 :: AA, vperpb, vparab, vperp, vpara, kperp, kpara, kall, CC
 
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8 :: ff !分布関数
 
 
 do s = 1, kind
  do j = 1, Z
   if(j == 1) then !##
   !ファイル名作成
   write(filename1, *) s
   write(filename2, '(1PE25.15E3)') mu(s, j)
   filename = 'PDS_J_30kV_disfun_s='//trim(adjustl(filename1))//'_mu='//trim(adjustl(filename2))//'.csv'
   print *, filename
   
   open(62, file=filename)
   
   do i = 1, N
    do p = 2, 3
     if(p == 2 .or. (p == 3 .and. a3(i, s, j, 3) > a3(i, s, j, 1))) then
       do q = 1, Z/2
        if(mod(q, 3) == 1) then
          if((p == 2 .and. sigmab(s) == 1) .or. (p == 3 .and. sigmab(s) == N)) then
            AA = a3(i, s, j, 1)*(a3(i, s, j, p)/a3(i, s, j, 1))**(dble(q-1)/dble(Z/2-1))
           else if((p == 2 .and. sigmab(s) == N) .or. (p == 3 .and. sigmab(s) == 1)) then
            AA = - a3(i, s, j, 1)*(a3(i, s, j, p)/a3(i, s, j, 1))**(dble(q-1)/dble(Z/2-1))
          endif
          
          CC = ff(sigma(s), BB(i), BB(sigmab(s)), Tp(s), Ta(s), UU(i, s), UU(sigmab(s), s), mu(s, j), AA, mass(s))
          
          vpara = sqrt(2.d0/mass(s))*AA/abs(AA)*sqrt(abs(UU(sigmab(s), s)-UU(i, s)+mu(s, j)*(BB(sigmab(s))-BB(i))+AA**2.d0))
          vperp = sqrt(2.d0*mu(s, j)*BB(i)/mass(s))
          vparab = sqrt(2.d0/mass(s))*AA
          vperpb = sqrt(2.d0*mu(s, j)*BB(sigmab(s))/mass(s))
          kpara = mass(s)*vpara**2.d0 /2.d0 /ee
          kperp = mass(s)*vperp**2.d0 /2.d0 /ee
          kall = (kpara+2.d0*kperp)/3.d0
          
          write(62, 74) i, vperpb, vparab, vperp, vpara, kperp, kpara, kall, CC, 2*pi*vperp*CC
        endif
       enddo !q
     endif
    enddo !p
   enddo !i
   
   close(62)
   endif !##
  enddo !j
 enddo !s
 
 
 return
 
 74 format(I3, 9(',', 1PE25.15E3))
 
end subroutine disfun2








