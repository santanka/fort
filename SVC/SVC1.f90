program SVC1

implicit none !い　つ　も　の

!全体を通しての文字
integer :: i, j !整数型定数
integer, parameter :: N = 79 !整数型定数
real*8, parameter :: pi = 4.*atan(1.)
character(len=80) :: dummy

!基準
real*8 :: B, e, c, m

!mag_FA関連の文字
real*8, dimension(N) :: LMD, BB, dIo

!SVC1_parameter関連の文字
real*8, dimension(12) :: parameter
real*8 :: mu0, GG, alpha, MJ, MassJ, RJ, omJ, Req, MIo

!SVC_IC関連の文字
real*8, dimension(N) :: pe0, dk, dR !実数型行列
integer, dimension(N) :: gr !整数型行列

!SVC_BC関連の文字
character(len=15), dimension(8) :: sp
character(len=1), dimension(8) :: so
real*8, dimension(8) :: Nd, Tp, Ta, ch, mass

!iteration以降の文字
real*8, dimension(8) :: sigma, sigmaN
integer, dimension(8) :: sigmab, sigmac
real*8, dimension(N) :: pe, phi, peN
integer :: itn, ii, ss, jj, pp, ud, iii
real*8, dimension(60) :: mu
real*8 :: amin, amax, alim, nis, rhn, GF, theta
real*8, dimension(3, N, 8) :: number, UU, nu
real*8, dimension(3, N) :: RHO, pep
real*8, dimension(3, N, 8, 60, 3) :: a3
real*8, dimension(3, N, 8, 60) :: THT
real*8, dimension(3, N) :: iter


!function
real*8 :: AI

!mag_FAの抽出
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i)
enddo

B = BB(1)
BB = BB/B

!SVC1_parameterの抽出
open(30, file="SVC1_parameter.csv", action="read", status="old")
do i = 1, 12
 read(30, *) dummy, parameter(i)
end do
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

!SVC_ICの抽出
open(10, file="SVC_IC.csv", action="read", status="old")
read(10, *) !1行読み飛ばし
do i = 1, N
 read(10, *) gr(i), pe0(i), dk(i), dR(i)
end do

pe0 = 0.d0!pe0*e/m/c**2.d0
dk = dk*e*B/m/c

!SVC_BCの抽出
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1行読み飛ばし
do i = 1, 8
 read(20, *) sp(i), Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do

Nd = Nd*(m*c/e/B)**3.d0
Tp = Tp*e/m/c/c
mass = mass/m

!sigma作成
sigma = Nd
do i = 1, 8
 if(ch(i) > 0) sigmab(i) = 1
 if(ch(i) < 0) sigmab(i) = N
 if(so(i) == "M") sigmac(i) = 1
 if(so(i) == "I") sigmac(i) = N
enddo


!iteration スタート
pe = pe0
itn = 0
do
 itn = itn + 1
 call pepm(pe, pep, e, m, c) !静電ポテンシャル
 call ep(N, ch, pep, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU) !ポテンシャル
 
 do i = 1, 3
  do j = 1, 8
   UU(i, :, j) = UU(i, :, j) - UU(i, 1, j)
  enddo
 enddo
 
 !UU = 0.d0
 call SUMMATION(N, sigma, BB, UU, Ta, Tp, nu, sigmac) !規格化因子
 
 do ud = 1, 3 !静電ポテンシャルの選択
  do ii = 1, N !磁力線の各gridで考える
   do ss = 1, 8 !各粒子種で考える
    do jj = 1, 60 !各Vperp(断熱不変量)で考える
     
     mu(jj) = AI(alpha, Tp(ss), BB(sigmac(ss)), jj) !断熱不変量作成
     call aaa(ii, N, alpha, mu(jj), Tp(ss), Ta(ss), UU(ud, :, ss), BB, amin, amax, alim, sigmac(ss)) !a関連作成
     a3(ud, ii, ss, jj, 1) = amin
     a3(ud, ii, ss, jj, 2) = amax
     a3(ud, ii, ss, jj, 3) = alim
     call TH(N, ii, amin, amax, alim, sigma(ss), Ta(ss), Tp(ss), mu(jj), UU(ud, :, ss), BB, theta, nu(ud, ii, ss), sigmac(ss)) !theta作成
     
     THT(ud, ii, ss, jj) = theta
    enddo !jj
    call NN(alpha, Tp(ss), BB(sigmac(ss)), THT(ud, ii, ss, :), nis) !数密度作成
    number(ud, ii, ss) = nis!*nu(ud, ii, ss)
   enddo !ss
   call RR(ch, number(ud, ii, :), rhn) !電荷密度作成
   RHO(ud, ii) = rhn
  enddo !ii
  
 enddo !ud
 
do i = 1, N
 !print *, i, pe(i)
 print *, i, ch(5)*pe(i)*m*c*c/e
 print *, i, -GG*MassJ*mass(5)/Req/(cos(LMD(i))**2.d0)*m*c*c/e
 print *, i, -mass(5)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(i))**6.d0)/2.d0*m*c*c/e
 print *, i, -GG*MIo*mass(5)/dIo(i)*m*c*c/e
 print *, i, BB(i)*AI(alpha, Tp(5), BB(sigmac(5)), 40)*m*c*c/e
 print *, "    "
 !print *, "    "
enddo



 do i =1, N 
!  print *, i, " s ", nu(:, i, 2), "          ", nu(:, i, 7), " s ", i
!  print *, i, " p ", pe(i)
!  print *, i, " U ", UU(:, i, 2), "          ", UU(:, i, 2), " U ", i
  print *, i, " a ", a3(1, i, 3, 40, :), "          ", a3(1, i, 5, 40, :), " a ", i
!  print *, i, " t ", THT(:, i, 2, 1), "          ", THT(:, i, 7, 1), " t ", i
!  print *, i, " n ", number(:, i, 2), "          ", number(:, i, 7), " n ", i
!  print *, i, " r ", RHO(:, i)
 enddo
 print *, "    "
 
 if(itn == 1) then
   open(90, file="SVC_Uzero.csv")
   do i = 1, N
    write(90, 82) number(1, i, :)*(e*B/m/c)**3.d0, RHO(1, i)*e*(e*B/m/c)**3.d0, UU(1, i, :)*m*c*c/e
   enddo
   exit
 endif
 
 call yy(N, e, RHO, number, GF, iter) !収束値
 
 
 call NewtonPHI(N, pep(1,:), pep(2,:), pep(3,:), RHO(1, :), RHO(2, :), RHO(3, :), iter, peN) !ニュートン法(静電ポテンシャル)
 
 call NewtonSIGMA(N, sigma, RHO, so, number, iter, sigmaN) !ニュートン法(sigma)
 
 pe = peN
 sigma = sigmaN
 
! print *, "grid", "     ", "ionosperic e-", "     ", "magnetospheric cold e-"
! do i = 1, N
!  print *, i, number(1, i, 2)*(e*B/m/c)**3.d0, number(1, i, 7)*(e*B/m/c)**3.d0
! enddo
! print *, "   "
 
 
 
 
!ファイル化処理(フリーズ対策)
 open(50, file="SVC_penum.csv")
 do i = 1, N
  write(50, 72) pe(i)*m*c*c/e, number(1, i, :)*(e*B/m/c)**3.d0, RHO(1, i)*e*(e*B/m/c)**3.d0, iter(1, i)
 enddo
 close(50)
 
 open(60, file="SVC_pote.csv")
 do i = 1, N
  write(60, 82) pe(i)*m*c*c/e, UU(1, i, :)*m*c*c
 enddo
 close(60)
 
 
!収束確認
 print *, itn, GF
 print *, "   "
 if (GF <= 1.d-6) go to 200 !ニュートン法の収束
! if (itn == 50) go to 300 !一応
 do i = 1, N
  if (isnan(pe(i))) then
    go to 100
  endif
 enddo
 
enddo

300 print *, "itn == 50"
100 print *, "NaN発生"
200 print *, "finish"



!終了処理
close(10)
close(20)
close(30)
close(40)

!フォーマット
72 format(E25.15E3, 10(',', 1x, E25.15E3))
82 format(E25.15E3, 16(',', 1x, E25.15E3))
92 format(E25.15E3, 7(',', 1x, E25.15E3))

end program SVC1






!function

!断熱不変量データ作成
real*8 function AI(alpha, T, B, j)
 implicit none
 real*8, intent(in) :: alpha, T, B
 integer, intent(in) :: j
 AI = (1.d-30)*(alpha*T/B/1.d-30)**(dble(j-1)/dble(60-1))
 return
end function AI


!分布関数
real*8 function ff(sigma, B, B1, Ta, Tp, U, U1, mu, aa, nu)
 implicit none
 real*8, intent(in) :: sigma, B, B1, Ta, Tp, U, U1, mu, aa, nu
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 ff = sigma*B/sqrt(pi*Ta*Tp**3.d0)*exp(-B1*mu/Tp)*exp(-(U+B*mu+aa**2.d0-(U1+B1*mu))/Ta/Tp)!/nu
 return
end function ff


!iteration条件
subroutine yy(N, ee, rr, Num, GF, iter)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: ee
 real*8, dimension(3, N) :: rr !可変
 real*8, dimension(3, N, 8), intent(in) :: Num
 real*8, dimension(3, N) :: NE, NP
 real*8, dimension(3, N), intent(out) :: iter
 real*8, intent(out) :: GF
 integer :: i, j
 do j = 1, 3
  do i = 1, N
   NP(j, i) = Num(j, i, 1) + Num(j, i, 3) + Num(j, i, 4) + Num(j, i, 5)*2.d0 + Num(j, i, 6)
   NE(j, i) = Num(j, i, 2) + Num(j, i, 7) + Num(j, i, 8)
!   if(NP(j, i) == 0.d0) then
!     print *, j, i, "NP = 0"
!     if(j == 1) NP(j, i) = 1.d-50
!     if(j == 2) NP(j, i) = 1.d-40
!     if(j == 3) NP(j, i) = 1.d-60
!   endif
!   if(NE(j, i) == 0.d0) then
!     print *, j, i, "NE = 0"
!     if(j == 1) NE(j, i) = 1.d-50
!     if(j == 2) NE(j, i) = 1.d-60
!     if(j == 3) NE(j, i) = 1.d-40
!   endif
!   if(rr(j, i) == 0.d0) then
!     print *, j, i, "rho = 0"
!     if(j == 1) rr(j, i) = 1.d-30
!     if(j == 2) rr(j, i) = 1.d-20
!     if(j == 3) rr(j, i) = 1.d-40
!   endif
  enddo
 enddo
 
 GF = 0.d0
 do j = 1, 3
  do i = 1, N
   iter(j, i) = rr(j, i)**2.d0/NE(j, i)**2.d0
   if(j == 1) GF = GF + iter(j, i)
  enddo
 enddo
 
 GF = GF/dble(N)
 GF = sqrt(GF)
 
 return
end subroutine yy



!subroutine

!ポテンシャル作成
subroutine ep(N, ch, pep, GG, MassJ, mass, Req, LMD, omJ, dIo, MIo, UU)
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: ch, mass
 real*8, dimension(N), intent(in) :: LMD, dIo
 real*8, dimension(3, N), intent(in) :: pep
 real*8, intent(in) :: GG, MassJ, Req, omJ, MIo
 real*8, dimension(3, N, 8), intent(out) :: UU
 integer :: i, j, k
 
 do i = 1, 3
  do j = 1, N
   do k = 1, 8
    UU(i, j, k) = ch(k)*pep(i, j)
    UU(i, j, k) = UU(i, j, k) - GG*MassJ*mass(k)/Req/(cos(LMD(j))**2.d0)
    UU(i, j, k) = UU(i, j, k) - mass(k)*(omJ**2.d0)*(Req**2.d0)*(cos(LMD(j))**6.d0)/2.d0
    UU(i, j, k) = UU(i, j, k) - GG*MIo*mass(k)/dIo(j)
   enddo !k
  enddo !j
 enddo !i 
 return
 
end subroutine ep


!電位ポテンシャルの差分作成
subroutine pepm(p, pp, e, m, c)
 implicit none
 real*8, parameter :: dp = 1.d0
 real*8, intent(in) :: e, m, c
 integer, parameter :: N = 79
 real*8, dimension(N), intent(in) :: p
 real*8, dimension(3, N), intent(out) :: pp
 integer :: i
 do i = 1, N
  pp(1, i) = p(i)
  if(p(i) <= dp*e/m/c**2.d0) then
    pp(2, i) = p(i)+dp*e/m/c**2.d0
    pp(3, i) = p(i)-dp*e/m/c**2.d0
   else
    pp(2, i) = p(i)*1.5d0
    pp(3, i) = p(i)*0.5d0
  endif
 enddo
 return 
end subroutine pepm

!分布関数規格化因子
subroutine SUMMATION(N, sigma, BB, UU, Ta, Tp, nu, sigmab)
 implicit none
 integer, intent(in) :: N
 integer, dimension(8), intent(in) :: sigmab
 real*8, dimension(8), intent(in) :: sigma, Ta, Tp
 real*8, dimension(N), intent(in) :: BB
 real*8, dimension(3, N, 8), intent(in) :: UU
 real*8, dimension(3, N, 8), intent(out) :: nu
 integer :: i, j, k
 
 do k = 1, 3
  do i = 1, N
   do j = 1, 8
    nu(k, i, j) = sigma(j)*BB(i)*Ta(j)/(BB(sigmab(j))*(Ta(j)-1.d0)+BB(i))*exp((UU(k, sigmab(j), j)-UU(k, i, j))/Ta(j)/Tp(j))
   enddo !j
  enddo !i
 enddo !k
 
 return
end subroutine SUMMATION

!amin,amax,alim作成
subroutine aaa(ii, N, alpha, mu, T, Ta, UU, BB, amin, amax, alim, sigmab)
 implicit none
 integer, intent(in) :: N, ii, sigmab
 real*8, intent(in) :: alpha, mu, T, Ta
 real*8, dimension(N), intent(in) :: UU, BB
 real*8, intent(out) :: amin, amax, alim
 real*8, dimension(N) :: Bmu, UB
 real*8 :: AA, CC
 integer :: i
 
 do i = 1, N
  Bmu(i) = BB(i)*mu
 enddo
 
 UB = UU + Bmu
 
 if(sigmab == 1) then
   if (ii >= 2) then
     CC = UB(1)
     do i = 1, ii
      if (CC < UB(i)) then
       CC = UB(i)
      end if
     enddo
     amin = sqrt(CC - UB(ii))
    else
     amin = 0.d0
   end if
   
   if (N > ii+1) then
     CC = UB(ii+1)
     do i =ii+1, N
      if (CC < UB(i)) CC = UB(i)
     enddo
     amax = CC - UB(ii)
     if(amax < 0.d0) amax = 0.d0
     if(amax >= 0.d0) amax = sqrt(amax)
    else if (N == ii+1) then
     amax = UB(ii+1)-UB(ii)
     if(amax < 0.d0) amax = 0.d0
     if(amax >= 0.d0) amax = sqrt(amax)
    else
     amax = 0.d0
   end if
   
   alim = UB(1)+alpha*T*Ta-UB(ii)
   if(alim < 0.d0) alim = 0.d0
   if(alim >= 0.d0) alim = sqrt(alim)



    !要検証
  elseif(sigmab == N) then
   if (ii /= N) then
     CC = UB(ii)
     do i = ii, N
      if (CC < UB(i)) then
       CC = UB(i)
      end if
     enddo
     amin = sqrt(CC - UB(ii))
    else
     amin = 0.d0
   end if
   
   if (ii >= 3) then
     CC = UB(1)
     do i = 1, ii-1
      if (CC < UB(i)) CC = UB(i)
     enddo
     amax = CC - UB(ii)
     if(amax < 0.d0) amax = 0.d0
     if(amax >= 0.d0) amax = sqrt(amax)
    else if (ii == 2) then
     amax = UB(1)-UB(ii)
     if(amax < 0.d0) amax = 0.d0
     if(amax >= 0.d0) amax = sqrt(amax)
    else
     amax = 0.d0
   end if
   
   alim = UB(N)+alpha*T*Ta-UB(ii)
   if(alim < 0.d0) alim = 0.d0
   if(alim >= 0.d0) alim = sqrt(alim)
   
 endif
 return
end subroutine aaa


!theta作成
subroutine TH(N, ii, amin, amax, alim, sigma, Ta, Tp, mu, UU, BB, theta, nu, sigmab)
 implicit none
 integer, intent(in) :: N, ii, sigmab
 real*8, dimension(N), intent(in) :: UU, BB
 real*8, intent(in) :: amax, alim, sigma, Ta, Tp, mu, nu
 real*8 :: amin !可変
 real*8, intent(out) :: theta
 real*8 :: thetaL, thetaM, the
 real*8, dimension(30) :: aL, aM
 integer :: pp
 real*8 :: ff !function
 
 thetaL = 0.d0
 thetaM = 0.d0
 
 if (alim > amin) then
  if(amin == 0.d0) amin = 1.d-30
  do pp = 1, 30-1
   aL(pp) = amin*((alim/amin)**(dble(pp-1)/dble(30-1)))
  enddo
  do pp = 1, 30-1
   the = ff(sigma, BB(ii), BB(sigmab), Ta, Tp, UU(ii), UU(sigmab), mu, aL(pp), nu)! + &
               ! &ff(sigma, BB(ii), BB(sigmab), Ta, Tp, UU(ii), UU(sigmab), mu, aL(pp+1), nu))/2.d0
   thetaL = thetaL + the*amin*((alim/amin)**(dble(pp)/dble(30-1))-(alim/amin)**(dble(pp-1)/dble(30-1)))
  enddo
 end if
 
 if (amax > amin .and. alim /= 0.d0) then
  if(amin == 0.d0) amin = 1.d-30
  do pp = 1, 30-1
   aM(pp) = amin*(amax/amin)**(dble(pp-1)/dble(30-1))
  enddo
  do pp = 1, 30-1
   the = ff(sigma, BB(ii), BB(sigmab), Ta, Tp, UU(ii), UU(sigmab), mu, -aM(pp), nu) !+ &
                !&ff(sigma, BB(ii), BB(sigmab), Ta, Tp, UU(ii), UU(sigmab), mu, aM(pp+1), nu))/2.d0
   thetaM = thetaM + the*amin*((amax/amin)**(dble(pp)/dble(30-1))-(amax/amin)**(dble(pp-1)/dble(30-1)))
  enddo
 end if
 
 theta = thetaL + thetaM
 
 return
end subroutine TH


!数密度作成
subroutine NN(alpha, Tp, B, THT, nis)
 implicit none
 real*8, intent(in) :: alpha, Tp, B
 real*8, dimension(60), intent(in) :: THT
 real*8, intent(out) :: nis
 integer :: p
 
 nis = 0.d0
 
 do p = 1, 60-1
  nis = nis + THT(p)*(1.d-30)*((alpha*Tp/B/1.d-30)**(dble(p)/dble(60-1)) &
            & - (alpha*Tp/B/1.d-30)**(dble(p-1)/dble(60-1)))  
                    !nis = nis + THT(p)*alpha*Tp/B/6.d1
 enddo
 
 return
end subroutine NN


!電荷密度作成
subroutine RR(qq, nn, rho)
 implicit none
 real*8, dimension(8), intent(in) :: qq, nn
 real*8, intent(out) :: rho
 integer :: i
 
 rho = 0.d0
 do i = 1, 8
  rho = rho + qq(i)*nn(i)
 enddo
 
 return
end subroutine RR


!ニュートン法
subroutine NewtonPHI(N, HH, Hu, Hd, II, Iu, Id, iter, HNew)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(N), intent(in) :: HH, Hu, Hd, II, Id
 real*8, dimension(3, N) :: iter
 real*8, dimension(N), intent(out) :: HNew
 real*8, dimension(N) :: Iu !可変
 integer :: i
 
 HNew(N) = HH(N)
 
 do i = 2, N-1
  if(iter(2, i) == iter(3, i)) then
    iter(2, i) = iter(2, i) + iter(2, i)*1.d-15
    print *,  i, "caution PHI"
  endif
  HNew(i) = HH(i) - (Hu(i) - Hd(i))/(iter(2, i) - iter(3, i))*iter(1, i)/1.d15
 enddo
 
 return
end subroutine NewtonPHI

subroutine NewtonSIGMA(N, SS, RR, so, Num, iter, sigmaN)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: SS
 real*8, dimension(3, N), intent(in) :: iter
 real*8, dimension(3, N) :: RR !可変
 character(len=1), dimension(8), intent(in) :: so 
 real*8, dimension(3, N, 8), intent(in) :: Num
 real*8, dimension(8), intent(out) :: sigmaN
 real*8 :: siM, siI
 integer :: s, siMs, siIs
 
 sigmaN = SS
 siM = 0.d0
 siI = 0.d0
 siMs = 0
 siIs = 0
 do s = 1, 8
  if(so(s) == "I" .and. Num(1, 1, s) == 0.d0 .and. siI < Num(1, N, s)) then
    siI = Num(1, N, s)
    siIs = s
   else if(so(s) == "M" .and. Num(1, N, s) == 0.d0 .and. siM < Num(1, 1, s)) then
    siM = Num(1, 1, s)
    siMs = s
  endif
 enddo
 
 if(RR(2, N) == RR(3, N)) then
   RR(2, N) = RR(2, N) + RR(2, N)*1.d-15
   print *, siIs, "caution SIGMA"
 endif
 if(RR(2, 1) == RR(3, 1)) then
   RR(2, 1) = RR(2, 1) + RR(2, 1)*1.d-15
   print *, siMs, "caution SIGMA"
 endif
 
 if(siIs /= 0) then
   sigmaN(siIs) = SS(siIs) - (Num(2, N, siIs) - Num(3, N, siIs))/&
      &(iter(2, N) - iter(3, N))*iter(1, N)/1.d10
 endif
 if(siMs /= 0) then
   sigmaN(siMs) = SS(siMs) - (Num(2, 1, siMs) - Num(3, 1, siMs))/&
      &(iter(2, 1) - iter(3, 1))*iter(1, 1)/1.d10
 endif
 return
end subroutine NewtonSIGMA

