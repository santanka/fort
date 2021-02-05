program SVC_Temperature

implicit none

!全体を通しての文字
integer :: i, s, j   !do構文用
integer, parameter :: N = 219 !grid数
integer, parameter :: Z = 180 !速度空間grid数
integer, parameter :: kind = 6 !粒子種数
integer, parameter :: spece = 3 !真・粒子種数
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
real*8, dimension(N, kind) :: Tperp, Tpara !温度[eV]
real*8, dimension(N, spece) :: taave !真・平均流速
real*8, dimension(N, kind, Z, 2) :: theta4 !a積分(真・温度)
real*8, dimension(N, spece) :: tTperp, tTpara, tT !真・温度
real*8, dimension(N, spece) :: bs !各β値
real*8, dimension(N) :: beta !β値
real*8, dimension(N) :: tTall !温度



!mag_FAの抽出
open(50, file="mag_FA_E_L=4_NS.csv", action="read", status="old")
do i = 1, N
 read(50, *) lambda(i), BB(i), dteq(i)
enddo !i
close(50)


!SVC_parameterの抽出
open(51, file="SVC_parameter_E_L=4.csv", action="read", status="old")
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
alpha = alpha*3.d0


!SVC_ICの抽出
open(52, file="SVC_IC_E_L=4_NS.csv", action="read", status="old")
read(52, *) !1行とばし
do i = 1, N
 read(52, *) dummy, dummy, ss(i), sR(i)
enddo !i
close(52)

do i = 1, N-1
 ds(i) = ss(i+1) - ss(i) !grid間距離
enddo !i


!SVC_BCの抽出
open(53, file="SVC_BC_E_L=4_NS.csv", action="read", status="old")
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
open(54, file="SVC_min_E_L=4_NS.csv", action="read", status="old")
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

!位置エネルギー
call EP(N, kind, GG, MassE, Req, omE, MassI, pe, lambda, dteq, ch, mass, UU)

!accessibility(1:min, 2:lim, 3:max)
call access(N, Z, kind, alpha, Tp, Ta, BB, UU, mu, sigmab, a3)

!aについて積分(数密度)
call TH1(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta1)

!muについて積分(数密度)
call NM1(N, kind, Z, BB, mass, mu, theta1, num)

!電荷密度
call RR(N, kind, ch, num, rhov)

!Poisson方程式
call Poisson(N, ep0, pe, ds, rhop)

!収束チェック
call CV(N, kind, rhov, rhop, ch, num, cvg)

!Alfven速度
call Alf(N, kind, mu0, cc, BB, ch, mass, num, VA)

!aについて積分(平均流速)
call TH2(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta2)

!muについて積分(平均流速)
call NM2(N, kind, Z, BB, mass, mu, theta2, num, aave)

!aについて積分(温度)
call TH3(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, aave, a3, theta3)

!muについて積分(温度)
call NM3(N, kind, Z, ee, BB, mass, mu, theta3, num, Tperp, Tpara)

!粒子種ごとの流速
call ryusoku(N, kind, spece, aave, num, UU, kk, sigmab, taave)

!粒子種ごとの温度(a積分)
call TH4(N, kind, Z, spece, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, kk, taave, a3, theta4)

!粒子種ごとの温度(mu積分)
call NM4(N, kind, Z, spece, ee, mu0, BB, mass, mu, theta4, num, UU, kk, sigmab, tTperp, tTpara, tT, bs, beta, tTall)


!ファイル化
open(60, file="SVC_E_L=4_NS_all.csv")
do i = 1, N
 write(60, 70) ss(i), sR(i), lambda(i), dteq(i), BB(i), pe(i), UU(i, :), num(i, :), rhov(i), rhop(i), cvg(i), &
  & VA(i, :), aave(i, :), Tperp(i, :), Tpara(i, :), taave(i, :), tTperp(i, :), tTpara(i, :), tT(i, :), tTall(i), &
  & bs(i, :), beta(i)
enddo !i
close(60)

open(61, file="SVC_E_L=4_NS_a3_all.csv")
do s = 1, kind
 do j = 1, Z
  write(61, 71) "s = ", s, "j = ", j
  do i = 1, N
   write(61, 72) "i = ", i, a3(i, s, j, :)
  enddo !i
  write(61, *) "  "
  write(61, *) "  "
 enddo !j
enddo !s
close(61)

70 format(1PE25.15E3, 58(',', 1PE25.15E3)) !5*kind+5*spece+14
71 format(A4, I2, ',', 2x, A4, I2)
72 format(A4, I3, ',' 1PE25.15E3, 2(',', 1PE25.15E3))

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
   mu(s, j) = 1.d-20*(alpha*Tp(s)/BB(sigmab(s))/1.d-20)**(dble(j-1)/dble(Z-1))
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
subroutine access(N, Z, kind, alpha, Tp, Ta, BB, UU, mu, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z, kind
 real*8, intent(in) :: alpha
 real*8, dimension(kind), intent(in) :: Tp, Ta
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 
 real*8, dimension(N, kind, Z, 3), intent(out) :: a3
 
 real*8, dimension(N, kind, Z) :: UB
 integer, dimension(N, kind, Z) :: lima
 integer :: i, s, j, k
 integer :: t, r !フリー
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
      a3(i, s, j, 1) = sqrt(UB(t, s, j) - UB(i, s, j))
    endif
    
    !"I"起源
    if(sigmab(s) == N) then
      t = N
      do k = i, N
       if(UB(k, s, j) > UB(t, s, j)) t = k
      enddo !k
      a3(i, s, j, 1) = sqrt(UB(t, s, j) - UB(i, s, j))
    endif
    
    if(a3(i, s, j, 1) == 0.d0) a3(i, s, j, 1) = 1.d-17
    
   enddo !j
  enddo !s
 enddo !i
 
 
 !alim
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"起源
    if(sigmab(s) == 1) then
      if(i == 1) a3(i, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
      if(i /= 1) then
        f = UB(1, s, j) + alpha*Tp(s)*Ta(s) - UB(i, s, j)
        if(f <= 0.d0) a3(i, s, j, 2) = 0.d0
        if(f > 0.d0) a3(i, s, j, 2) = sqrt(f)
        if(a3(i-1, s, j, 2) < a3(i-1, s, j, 1)) a3(i, s, j, 2) = 0.d0
      endif
    endif
    
    !"I"起源
    if(sigmab(s) == N) then
      k = N+1-i
      if(k == N) a3(k, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
      if(k /= N) then
        f = UB(N, s, j) + alpha*Tp(s)*Ta(s) - UB(k, s, j)
        if(f <= 0.d0) a3(k, s, j, 2) = 0.d0
        if(f > 0.d0) a3(k, s, j, 2) = sqrt(f)
        if(a3(k+1, s, j, 2) < a3(k+1, s, j, 1)) a3(k, s, j, 2) = 0.d0
      endif
    endif
    
   enddo !j
  enddo !s
 enddo !i
 
 
 !accessibilityの調整
 do s = 1, kind
  do j = 1, Z
   
   !"M"起源
   if(sigmab(s) == 1) then
     do i = 1, N
      if(i /= 1 .and. lima(i-1, s, j) == 1) lima(i, s, j) = 1
      if(a3(i, s, j, 2) <= a3(i, s, j, 1)) lima(i, s, j) = 1
     enddo !i
   endif
   
   !"I"起源
   if(sigmab(s) == N) then
     do i = 1, N
      k = N+1-i
      if(k /= N .and. lima(k+1, s, j) == 1) lima(k, s, j) = 1
      if(a3(k, s, j, 2) <= a3(k, s, j, 1)) lima(k, s, j) = 1
     enddo !i
   endif
   
  enddo !j
 enddo !s
 
 
 !amax
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
    !"M"起源
    if(sigmab(s) == 1) then
      if(i == N .or. lima(i, s, j) == 1) a3(i, s, j, 3) = 0.d0
      if(i /= N .and. lima(i, s, j) == 0) then
        t = N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        r = N
        do k = i+1, N
         if(lima(k-1, s, j) == 0 .and. lima(k, s, j) == 1) r = k
        enddo !k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k = i+1, r
         if(UB(k, s, j) > UB(t, s, j)) t = k
        enddo !k
        a3(i, s, j, 3) = UB(t, s, j) - UB(i, s, j)
        if(a3(i, s, j, 3) < 0.d0) a3(i, s, j, 3) = 0.d0
        if(a3(i, s, j, 3) >= 0.d0) a3(i, s, j, 3) = sqrt(a3(i, s, j, 3))
      endif
    endif
    
    !"I"起源
    if(sigmab(s) == N) then
      if(i == 1 .or. lima(i, s, j) == 1) a3(i, s, j, 3) = 0.d0
      if(i /= 1 .and. lima(i, s, j) == 0) then
        t = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        r = 1
        do k = 1, i-1
         if(lima(k, s, j) == 1 .and. lima(k+1, s, j) == 0) r = k
        enddo !k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k = r, i-1
         if(UB(k, s, j) > UB(t, s, j)) t = k
        enddo !k
        a3(i, s, j, 3) = UB(t, s, j) - UB(i, s, j)
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
 
 
 ff = sigma*sqrt((mass/2.d0/pi/Tp)**3.d0/Ta)*exp(-B1*mu/Tp)*exp(-(UU+BB*mu+aa**2.d0-(U1+B1*mu))/Ta/Tp)
 return
 
end function ff



!aについて積分(数密度)
subroutine TH1(N, kind, Z, BB, UU, sigma, mass, Tp, Ta, mu, sigmab, a3, theta1)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(N, kind), intent(in) :: UU
 real*8, dimension(kind), intent(in) :: sigma, mass, Tp, Ta
 real*8, dimension(kind, Z), intent(in) :: mu
 integer, dimension(kind), intent(in) :: sigmab
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z), intent(out) :: theta1
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !分布関数(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    thetaL = 0.d0
    thetaM = 0.d0
    
    !alim
    if(a3(i, s, j, 2) >= a3(i, s, j, 1) .and. a3(i, s, j, 2) /= 0.d0) then
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
    if(a3(i, s, j, 3) >= a3(i, s, j, 1) .and. a3(i, s, j, 3) /= 0.d0) then
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
    
    theta1(i, s, j) = thetaL + thetaM
   enddo !j
  enddo !s
 enddo !i 
 
 return
 
end subroutine TH1


!muについて積分(数密度)
subroutine NM1(N, kind, Z, BB, mass, mu, theta1, num)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z), intent(in) :: theta1
 
 real*8, dimension(N, kind), intent(out) :: num
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 do i = 1, N
  do s = 1, kind
   
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta1(i, s, j) + theta1(i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   num(i, s) = sqrt((2.d0/mass(s))**3.d0)*BB(i)*pi*nnn
   
  enddo !s
 enddo !i
 
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
 real*8, dimension(N, kind, Z, 3), intent(in) :: a3
 
 real*8, dimension(N, kind, Z), intent(out) :: theta2
 
 integer :: i, s, j, p
 real*8 :: thetaL, thetaM
 real*8, dimension(Z/2) :: aL, aM
 real*8 :: ff !分布関数(function)
 
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    
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
subroutine NM2(N, kind, Z, BB, mass, mu, theta2, num, aave)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z), intent(in) :: theta2
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N, kind), intent(out) :: aave
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 
 do i = 1, N
  do s = 1, kind
   
   !Vpara作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta2(i, s, j) + theta2(i, s, j+1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   aave(i, s) = (2.d0/mass(s))**(3.d0/2.d0)*pi*BB(i)/num(i, s)*nnn
   
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
    
    theta3(i, s, j, 1) = mu(s, j)*(thetaL + thetaM)
    
    
    !Tpara作成
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
    
    
   enddo !j
  enddo !s
 enddo !i 
 
 
 return
 
end subroutine TH3


!muについて積分(温度)
subroutine NM3(N, kind, Z, ee, BB, mass, mu, theta3, num, Tperp, Tpara)
 implicit none
 integer, intent(in) :: N, kind, Z
 real*8, intent(in) :: ee
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z, 2), intent(in) :: theta3
 real*8, dimension(N, kind), intent(in) :: num
 
 real*8, dimension(N, kind), intent(out) :: Tperp, Tpara
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 
 
 do i = 1, N
  do s = 1, kind
   
   !Tperp作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta3(i, s, j, 1) + theta3(i, s, j+1, 1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   
   Tperp(i, s) = 2.d0*sqrt(2.d0/mass(s)**3.d0)*(BB(i)**2.d0)*pi/num(i, s)*nnn /ee
   
   
   !Vpara作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta3(i, s, j, 2) + theta3(i, s, j+1, 2))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Tpara(i, s) = 4.d0*sqrt(2.d0/mass(s)**3.d0)*BB(i)*pi/num(i, s)*nnn /ee
   
  enddo !s
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
 
 taave = taave/nnn
 
 
 
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
subroutine NM4(N, kind, Z, spece, ee, mu0, BB, mass, mu, theta4, num, UU, kk, sigmab, tTperp, tTpara, tT, bs, beta, tTall)
 implicit none
 integer, intent(in) :: N, kind, Z, spece
 real*8, intent(in) :: ee, mu0
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(kind), intent(in) :: mass
 real*8, dimension(kind, Z), intent(in) :: mu
 real*8, dimension(N, kind, Z, 2), intent(in) :: theta4
 real*8, dimension(N, kind), intent(in) :: num, UU
 integer, dimension(kind), intent(in) :: kk, sigmab
 
 real*8, dimension(N, spece), intent(out) :: tTperp, tTpara, tT, bs
 real*8, dimension(N), intent(out) :: beta, tTall
 
 integer :: i, s, j
 real*8 :: nnn
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8, dimension(N, kind) :: Tperp, Tpara
 real*8, dimension(N, spece) :: nm
 real*8, dimension(N) :: nn
 
 
 do i = 1, N
  do s = 1, kind
   
   !Tperp作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta4(i, s, j, 1) + theta4(i, s, j+1, 1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   
   Tperp(i, s) = 2.d0*sqrt(2.d0/mass(s)**3.d0)*(BB(i)**2.d0)*pi/num(i, s)*nnn /ee
   
   
   !Tpara作成
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (theta4(i, s, j, 2) + theta4(i, s, j+1, 2))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Tpara(i, s) = 4.d0*sqrt(2.d0/mass(s)**3.d0)*BB(i)*pi/num(i, s)*nnn /ee
   
  enddo !s
 enddo !i
 
 
 tTperp = 0.d0
 tTpara = 0.d0
 nm = 0.d0
 
 
 do i = 1, N
  do s = 1, kind
   nm(i, kk(s)) = nm(i, kk(s)) + num(i, s)
   if(Tperp(i, s) == Tperp(i, s)) then
     tTperp(i, kk(s)) = tTperp(i, kk(s)) + num(i, s)*Tperp(i, s)
   endif
   if(Tpara(i, s) == Tpara(i, s)) then
     tTpara(i, kk(s)) = tTpara(i, kk(s)) + num(i, s)*Tpara(i, s)
   endif
  enddo !s
 enddo !i
 
 tTperp = tTperp/nm
 tTpara = tTpara/nm
 
 tT = (tTperp*2.d0+tTpara)/3.d0
 
 bs = 2.d0*mu0*nm*tT*ee
 beta = 0.d0
 do i = 1, N
  bs(i, :) = bs(i, :)/BB(i)**2.d0
  do s = 1, spece
   if(bs(i, s) == bs(i, s)) beta(i) = beta(i) + bs(i, s)
  enddo !s
 enddo !i
 
 nn = 0.d0
 do s = 1, spece
  nn = nn + nm(:, s)
 enddo !s
 
 tTall = beta*BB**2.d0/2.d0/mu0/nn /ee
 
 
 return
 
end subroutine NM4

