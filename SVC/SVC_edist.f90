program SVC4

implicit none

!全体を通しての文字
integer :: i, j, s !do用
integer, parameter :: N = 79 !grid数
integer, parameter :: Z = 1000 !速度空間grid数
integer, parameter :: AC = 57 !分布チェックするgrid
integer, parameter :: AS = 7 !分布チェックするplasma
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
real*8, dimension(N) :: pep, rhov, rhop, cvg
real*8, dimension(N, 8) :: UU, num
real*8, dimension(N, 8, Z) :: UB
real*8, dimension(8, Z) :: mu
real*8, dimension(N, 8, Z, 3) :: a3
real*8, dimension(N, 8, Z) :: theta
real*8, dimension(2) :: sigx
real*8, dimension(8) :: nsig
real*8 :: cvn, DDD


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
alpha = 5.d2
MJ = parameter(7) !Jupiterの磁気双極子モーメント
MassJ = parameter(8) !Jupiterの質量
RJ = parameter(9) !Jupiter半径
omJ = parameter(10) !Jupiter自転角周波数
Req = parameter(11) !Jupiterのreq
MIo = parameter(12) !Ioの質量

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
open(15, file="SVC_penum.csv", action="read", status="old")
do i = 1, N
 !if(i >= 45 .and. i <= 55) read(15, *) dummy
 !if(i < 45 .or. 55 < i) read(15, *) pe0(i)
 read(15, *) pe0(i)
enddo !i
close(15)

open(18, file="SVC_nsig.csv", action="read", status="old")
do s = 1, 8
 read(18, *) sigma(s)
enddo !s
close(18)



!sigmab作成
do i = 1, 8
 if(so(i) == "M") sigmab(i) = 1
 if(so(i) == "I") sigmab(i) = N
enddo


!iteration スタート

pep = pe0
itn = 0

open(33, file = "SVC_min.csv")
read(33, *) DDD
!------------------------------
do i = 1, N
 read(33, *) pep(i)
enddo !i
do s = 1, 8
 read(33, *) sigma(s)
enddo !s
!------------------------------
close(33)

call AI(Z, alpha, N, Tp, BB, sigmab, mu)

!ポテンシャル
call EP(N, ch, pep, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU)

!エネルギー(位置+Kperp)
call UUBB(N, Z, UU, BB, mu, UB)

!accessibility(1:min, 2:lim, 3:max)
call aaa(N, Z, alpha, Tp, Ta, UB, sigmab, a3)

call dis(N, Z, AC, AS, a3, sigma, Tp, Ta, sigmab, mass, UU, BB, mu)


62 format(E25.15E3, 3(',', E25.15E3))
72 format(E25.15E3, 11(',', E25.15E3))
82 format(E25.15E3, 21(',', 1x, E25.15E3))
92 format(E25.15E3, 8(',', E25.15E3))


end program SVC4





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
   mu(i, j) = 1.d-30*(alpha*Tp(i)/BB(sigmab(i))/1.d-30)**(dble(j-1)/dble(Z-1))
   !mu(i, j) = alpha*Tp(i)/BB(sigmab(i))*dble(j-1)/dble(Z-1)
  enddo !j
 enddo !i
 
 return
 
end subroutine AI


!ポテンシャル
subroutine EP(N, ch, pe, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch, mass
 real*8, dimension(N), intent(in) :: pe, LMD, dIo
 real*8, intent(in) :: GG, MassJ, Req, omJ, MIo
 real*8, dimension(N, 8), intent(out) :: UU
 integer ::i, s
 
 do i = 1, N
  do s = 1, 8
   !静電
   UU(i, s) = ch(s)*pe(i)
   !木星重力
   UU(i, s) = UU(i, s) - GG*MassJ*mass(s)/Req/(cos(LMD(i))**2.d0)
   !木星遠心力
   UU(i, s) = UU(i,s) - mass(s)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0
   !イオ重力
   UU(i, s) = UU(i, s) - GG*MIo*mass(s)/dIo(i)
  enddo !s
 enddo !i
 
 do s = 1, 8
  UU(:, s) = UU(:, s) - UU(1, s)
 enddo !s
 
 return
 
end subroutine EP


!エネルギー(位置+Kperp)
subroutine UUBB(N, Z, UU, BB, mu, UB)
 implicit none
 integer, intent(in) :: N, Z
 real*8, dimension(N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 real*8, dimension(N, 8, Z), intent(out) :: UB
 integer :: i, s, j
 
 do i = 1, N
  do s = 1, 8
   do j = 1, Z
    UB(i, s, j) = UU(i, s) + BB(i)*mu(s, j)
   enddo !j
  enddo !s
 enddo !i
 
 return
 
end subroutine UUBB


subroutine aaa(N, Z, alpha, Tp, Ta, UB, sigmab, a3)
 implicit none
 integer, intent(in) :: N, Z
 real*8, intent(in) :: alpha
 real*8, dimension(8), intent(in) :: Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(N, 8, Z) :: UB
 real*8, dimension(N, 8, Z, 3) :: a3
 integer, dimension(N, 8, Z) :: lima
 integer :: ud, i, s, j, k
 real*8 :: CC
 
 lima = 0
 
 !amin
 do i = 1, N
  do s = 1, 8
   do j = 1, Z
   
   !イオ起源
   if(sigmab(s) == 1) then
     if(i == 1) a3(i, s, j, 1) = 0.d0
     if(i == 2) then
       if(UB(2, s, j) >= UB(1, s, j)) a3(i, s, j, 1) = 0.d0
       if(UB(2, s, j) < UB(1, s, j)) then
         a3(i, s, j, 1) = sqrt(UB(1, s, j) - UB(i, s, j))
       endif
     endif
     if(i >= 3) then
       CC = UB(i, s, j)
       do k = 1, i-1
        if(UB(k, s, j) > CC) CC = UB(k, s, j)
       enddo !k
       a3(i, s, j, 1) = sqrt(CC - UB(i, s, j))
     endif
   endif
   
   !木星起源
   if(sigmab(s) == N) then
     if(i == N) a3(i, s, j, 1) = 0.d0
     if(i == N-1) then
       if(UB(N-1, s, j) >= UB(N, s, j)) a3(i, s, j, 1) = 0.d0
       if(UB(N-1, s, j) < UB(N, s, j)) then
         a3(N-1, s, j, 1) = sqrt(UB(N, s, j) - UB(N-1, s, j))
       endif
     endif
     if(i <= N-2) then
       CC = UB(i, s, j)
       do k = i+1, N
        if(UB(k, s, j) > CC) CC = UB(k, s, j)
       enddo !k
       a3(i, s, j, 1) = sqrt(CC - UB(i, s, j))
     endif
   endif
   
   if(a3(i, s, j, 1) == 0.d0) a3(i, s, j, 1) = 1.d-15
   
   enddo !j
  enddo !s
 enddo !i
 
 !alim
 do i = 1, N
  do s = 1, 8
   do j = 1, Z
    
    CC = 0.d0
    
    !イオ起源
    if(sigmab(s) == 1) then
      if(i == 1) a3(i, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
      if(i /= 1) then
        CC = UB(1, s, j) + alpha*Tp(s)*Ta(s) - UB(i, s, j)
        if(CC <= 0.d0) a3(i, s, j, 2) = 0.d0
        if(CC > 0.d0) a3(i, s, j, 2) = sqrt(CC)
        
        if(a3(i-1, s, j, 2) < a3(i-1, s, j, 1)) a3(i, s, j, 2) = 0.d0
      endif
    endif
    
    !木星起源
    if(sigmab(s) == N) then
      k = N+1-i
      if(k == N) a3(k, s, j, 2) = sqrt(alpha*Tp(s)*Ta(s))
      if(k /= N) then
        CC = UB(N, s, j) + alpha*Tp(s)*Ta(s) - UB(k, s, j)
        if(CC <= 0.d0) a3(k, s, j, 2) = 0.d0
        if(CC > 0.d0) a3(k, s, j, 2) = sqrt(CC)
        
        if(a3(k+1, s, j, 2) < a3(k+1, s, j, 1)) a3(k, s, j, 2) = 0.d0
      endif
    endif
    
   enddo !j
  enddo !s
 enddo !i
 
 !accessibilityの調整
 do s = 1, 8
  do j = 1, Z
   
   !イオ起源
   if(sigmab(s) == 1) then
     do i = 1, N
      if(i /= 1 .and. lima(i-1, s, j) == 1) lima(i, s, j) = 1
      if(a3(i, s, j, 2) <= a3(i, s, j, 1)) lima(i, s, j) = 1
     enddo !i
   endif
   
   !木星起源
   if(sigmab(s) == N) then
     do i = 1, N
      if(i /= 1 .and. lima(N+2-i, s, j) == 1) lima(N+1-i, s, j) = 1
      if(a3(N+1-i, s, j, 2) <= a3(N+1-i, s, j, 1)) lima(N+1-i, s, j) = 1
     enddo !i
   endif
   
  enddo !j
 enddo !s
 
 
 !amax
 do i = 1, N
  do s = 1, 8
   do j = 1, Z
    
    CC = 0.d0
     
    !イオ起源
    if(sigmab(s) == 1) then
      if(i == N .or. lima(i, s, j) == 1) a3(i, s, j, 3) = 0.d0
      if(i+1 == N .and. lima(i, s, j) == 0) then
        CC = UB(N, s, j) - UB(N-1, s, j)
        if(CC <= 0.d0) a3(i, s, j, 3) = 0.d0
         if(CC > 0.d0) a3(i, s, j, 3) = sqrt(CC)
      endif
      if(i+1 < N .and. lima(i, s, j) == 0) then
        do k = i+1, N
         if(CC < UB(k, s, j)) CC = UB(k, s, j)
        enddo !k
        CC  = CC - UB(i, s, j)
        if(CC < 0.d0) a3(i, s, j, 3) = 0.d0
        if(CC >= 0.d0) a3(i, s, j, 3) = sqrt(CC)
      endif
    endif
    
    !木星起源
    if(sigmab(s) == N) then
      if(i == 1 .or. lima(i, s, j) == 1) a3(i, s, j, 3) = 0.d0
      if(i == 2 .and. lima(i, s, j) == 0) then
        CC = UB(1, s, j) - UB(i, s, j)
        if(CC < 0.d0) a3(i, s, j, 3) = 0.d0
        if(CC >= 0.d0) a3(i, s, j, 3) = sqrt(CC)
      endif
      if(i > 2 .and. lima(i, s, j) == 0) then
        do k = 1, i-1
         if(CC < UB(k, s, j)) CC = UB(k, s, j)
        enddo !k
        CC = CC - UB(i, s, j)
        if(CC < 0.d0) a3(i, s, j, 3) = 0.d0
        if(CC >= 0.d0) a3(i, s, j, 3) = sqrt(CC)
      endif
    endif
    
   enddo !j
  enddo !s
 enddo !i
 
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


!電子分布
subroutine dis(N, Z, AC, AS, a3, sigma, Tp, Ta, sigmab, mass, UU, BB, mu)
 implicit none
 integer, intent(in) :: N, Z, AC, AS
 real*8, dimension(N, 8, Z, 3), intent(in) :: a3
 real*8, dimension(8), intent(in) :: sigma, Tp, Ta
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(8), intent(in) :: mass
 real*8, dimension(N, 8), intent(in) :: UU
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(8, Z), intent(in) :: mu
 integer :: j, p
 real*8, dimension(Z, Z/2) :: al, am
 real*8 :: ff !function
 real*8, parameter :: pi = 4.d0*atan(1.d0) !円周率
 real*8 :: LL, MM
 
 !alim
 do j = 1, Z
  if(a3(AC, AS, j, 2) > a3(AC, AS, j, 1)) then
    do p = 1, Z/2
     al(j, p) = a3(AC, AS, j, 1)*(a3(AC, AS, j, 2)/a3(AC, AS, j, 1))**(dble(p-1)/dble(Z/2-1))
     !al(j, p) = a3(AC, AS, j, 1)+(a3(AC, AS, j, 2)-a3(AC, AS, j, 1))*(dble(p-1)/dble(Z/2-1))
    enddo !p
   else
    al(j, :) = 0.d0
  endif
 enddo !j
 
 !amax
 do j = 1, Z
  if(a3(AC, AS, j, 3) > a3(AC, AS, j, 1)) then
    do p = 1, Z/2
     am(j, p) = -a3(AC, AS, j, 1)*(a3(AC, AS, j, 3)/a3(AC, AS, j, 1))**(dble(p-1)/dble(Z/2-1))
     !am(j, p) = -(a3(AC, AS, j, 1)+(a3(AC, AS, j, 3)-a3(AC, AS, j, 1))*(dble(p-1)/dble(Z/2-1)))
    enddo !p
   else
    am(j, :) = 0.d0
  endif
 enddo !j
 
 !ファイル化
 open(74, file="SVC_edist.csv")
 do j = 1, Z
  do p = 1, Z/2
   if(a3(AC, AS, j, 2) > a3(AC, AS, j, 1)) then
     LL = ff(sigma(AS), BB(AC), BB(sigmab(AS)), Ta(AS), Tp(AS), UU(AC, AS), UU(sigmab(AS), AS), mu(AS, j), al(j, p))&
       &/(2.d0/mass(AS))**(3.d0/2.d0)/pi/BB(AC)*sqrt(2.d0*mu(AS, j)*BB(AC)/mass(AS))*2.*pi
     if(LL < 1.d-20) LL = 0.d0
    else
     LL = 0.d0
   endif
   write(74, 72) al(j, p)/sqrt(mass(AS)/2.d0), sqrt(2.d0*mu(AS, j)*BB(AC)/mass(AS)), &
                &LL
  enddo !p
 enddo !j
 
 do j = 1, Z
  do p = 1, Z/2
   if(a3(AC, AS, j, 3) > a3(AC, AS, j, 1)) then
     MM = ff(sigma(AS), BB(AC), BB(sigmab(AS)), Ta(AS), Tp(AS), UU(AC, AS), UU(sigmab(AS), AS), mu(AS, j), am(j, p))&
       &/(2.d0/mass(AS))**(3.d0/2.d0)/pi/BB(AC)*sqrt(2.d0*mu(AS, j)*BB(AC)/mass(AS))*2.*pi
     if(MM < 1.d-20) MM = 0.d0
    else
     MM = 0.d0
   endif
   write(74, 72) am(j, p)/sqrt(mass(AS)/2.d0), sqrt(2.d0*mu(AS, j)*BB(AC)/mass(AS)), &
                &MM
  enddo !p
 enddo !j
 
 close(74)
 
 return
 
 72 format(E25.15E3, 2(',', E25.15E3))
 
end subroutine dis


