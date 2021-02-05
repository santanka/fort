program SVC

implicit none !い　つ　も　の

!全体を通しての文字
integer :: i !整数型定数
real*8, parameter :: pi = 4.*atan(1.)
character(len=80) :: dummy

!SVC_parameter関連の文字
real*8, dimension(11) :: parameter
real*8 :: RJ, mu0, ME, MJ, Req, alpha, GG, omJ, MassJ, MIo, ee

!SVC_IC関連の文字
integer, parameter :: N = 79 !整数型定数
real*8, dimension(N) :: pe0, dk, dR !実数型行列
integer, dimension(N) :: gr !整数型行列

!SVC_BC関連の文字
character(len=15), dimension(8) :: sp
character(len=1), dimension(8) :: so
real*8, dimension(8) :: Nd, Tp, Ta, ch, mass

!mag_FA関連の文字
real*8, dimension(N) :: LMD, BB, dIo

!iteration以降の文字
real*8, dimension(8) :: sigma, sigmaN
real*8, dimension(N) :: pe, pem, pep, phi, UU, peN
integer :: itn, ii, ss, jj, pp, ud, iii
real*8, dimension(60) :: mu, THT
real*8 :: amin, amax, alim, nis, rhn, GF, theta
real*8, dimension(3, N, 8) :: number
real*8, dimension(3, N) :: RHO

!function
real*8 :: ep, AI, yy


!SVC_parameterの抽出
open(30, file="SVC_parameter.csv", action="read", status="old")
do i = 1, 11
 read(30, *) dummy, parameter(i)
end do
mu0 = parameter(1)
RJ = parameter(2)
ME = parameter(3)
MJ = parameter(4)
Req = parameter(5)
alpha = parameter(6)
GG = parameter(7)
omJ = parameter(8)
MIo = parameter(9)
MassJ = parameter(10)
ee = parameter(11)


!SVC_ICの抽出
open(10, file="SVC_IC.csv", action="read", status="old")
read(10, *) !1行読み飛ばし

!行列にデータを代入
do i = 1, N
 read(10, *) gr(i), pe0(i), dk(i), dR(i)
end do


!SVC_BCの抽出
open(20, file="SVC_BC.csv", action="read", status="old")
read(20, *) !1行読み飛ばし

do i = 1, 8
 read(20, *) sp(i), Nd(i), Tp(i), Ta(i), so(i), ch(i), mass(i)
end do

ch = ch*ee ![C]に変換
Tp = Tp*ee ![eV]を[J]に変換


!mag_FAの抽出
open(40, file="mag_FA.csv", action="read", status="old")
do i = 1, N
 read(40, *) LMD(i), BB(i), dIo(i)
enddo


!初期sigma作成
do i = 1, 8
 if (so(i) == "I") sigma(i) = Nd(2)
 if (so(i) == "M") sigma(i) = Nd(3)
enddo


!iteration スタート
pe = pe0
itn = 0
do
 itn = itn + 1
 call pepm(pe, pep, pem)
 do ud = 1, 3 !静電ポテンシャルの選択
  if (ud == 1) phi = pe
  if (ud == 2) phi = pep
  if (ud == 3) phi = pem
  do ii = 1, N !磁力線の各gridで考える
   
   do ss = 1, 8 !各粒子種で考える
    do iii = 1, N
     UU(iii) = ep(ch(ss), phi(iii), GG, MassJ, mass(ss), Req, LMD(iii), omJ, dIo(iii), MIo) !ポテンシャル作成
    enddo !iii
    do jj = 1, 60 !各Vperp(断熱不変量)で考える
     mu(jj) = AI(alpha, Tp(ss), BB(1), jj) !断熱不変量作成
     call aaa(ii, N, alpha, mu(jj), Tp(ss), Ta(ss), UU, BB, amin, amax, alim) !a関連作成
     
     call TH(N, ii, amin, amax, alim, sigma(ss), Ta(ss), Tp(ss), mu(jj), UU, BB, theta) !theta作成
     
     THT(jj) = theta
    enddo !jj
    call NN(alpha, Tp(ss), BB(1), THT, nis) !数密度作成
    number(ud, ii, ss) = nis
   enddo !ss
   call RR(ch, number(ud, ii, :), rhn) !電荷密度作成
   RHO(ud, ii) = rhn
  enddo !ii
  
 enddo !ud
 
 call NewtonPHI(N, pe, pep, pem, RHO(1, :), RHO(2, :), RHO(3, :), peN) !ニュートン法(静電ポテンシャル)
 
 call NewtonSIGMA(N, sigma, RHO, so, number, sigmaN) !ニュートン法(sigma)
 
 do i = 1, N
 print *, i, number(1, i, 2), number(1, i, 3)
 enddo
 print *, "     "
 
 pe = peN
 sigma = sigmaN
 
 GF = yy(N, ee, RHO(1, :), number(1, :, :))
 if (GF <= 1.d-7) exit !ニュートン法の収束
 if (itn == 1000) exit !一応
 print *, itn, GF
 print *, "   "
 if (isnan(GF)) then
  write(*,*) "NaN発生"
  exit
 endif

enddo


!ファイル化処理
open(50, file="SVC_penum.csv")
do i = 1, N
 write(50, 72) pe(i), number(1, i, :)
enddo


!終了処理
close(10)
close(20)
close(30)
close(40)
close(50)

!フォーマット
72 format(E15.7, 8(',', 1x, E15.7))

end program SVC






!function
!ポテンシャル作成
real*8 function ep(qq, phi, GG, MassJ, mass, Req, lam, omJ, dIo, MIo)
 real*8, intent(in) :: qq, phi, GG, MassJ, mass, Req, lam, omJ, dIo, MIo
 ep = qq*phi - GG*MassJ*mass/Req/(cos(lam)**2.) + mass*(omJ**2.)*(Req**2.)*(cos(lam)**6.)/2. !- GG*MIo*mass/dIo 
 return
end function ep


!断熱不変量データ作成
real*8 function AI(alpha, T, B, j)
 implicit none
 real*8, intent(in) :: alpha, T, B
 integer, intent(in) :: j
 real*8 :: jj
 jj = dble(j) !倍精度に変換
 AI = alpha*T/B*(jj-1.)/60.
 return
end function AI


!分布関数
real*8 function ff(sigma, B, B1, Ta, Tp, U, U1, mu, aa)
 implicit none
 real*8, intent(in) :: sigma, B, B1, Ta, Tp, U, U1, mu, aa
 real*8, parameter :: pi = 4.*atan(1.)
 ff = sigma*B/sqrt(pi*Ta*Tp**3)*exp(-B1*mu/Tp)*exp(-(U+B*mu+aa**2-(U1+B1*mu))/Ta/Tp)
 return
end function ff


!iteration条件
real*8 function yy(N, ee, rr, Num)
 implicit none
 integer, intent(in) :: N
 real*8, intent(in) :: ee
 real*8, dimension(N), intent(in) :: rr
 real*8, dimension(N, 8), intent(in) :: Num
 real*8, dimension(N) :: NE
 integer :: i
 
 do i = 1, N
  NE(i) = Num(i, 2) + Num(i, 7) + Num(i, 8)
 enddo
 
 yy = 0.
 do i = 1, N
  yy = yy + (rr(i)/ee/NE(i))**2.
 enddo
 
 yy = yy/dble(N)
 yy = sqrt(yy)
 
 return
end function yy






!subroutine
!電位ポテンシャルの差分作成
subroutine pepm(p, pp, pm)
 implicit none
 real*8, parameter :: dp = 100
 integer, parameter :: N = 79
 real*8, dimension(N), intent(in) :: p
 real*8, dimension(N), intent(out) :: pp, pm
 integer :: i
 pp = p + dp
 pm = p - dp
 
 return 
end subroutine pepm


!amin,amax,alim作成
subroutine aaa(ii, N, alpha, mu, T, Ta, UU, BB, amin, amax, alim)
 implicit none
 integer, intent(in) :: N, ii
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
 
 if (ii >= 2) then
   CC = UB(1)
   do i = 1, ii
    if (CC < UB(i)) then
     CC = UB(i)
    end if
   enddo
   amin = sqrt(CC - UB(ii))
  else
   amin = 0.
 end if
 
 if (N > ii+1) then
   CC = UB(ii+1)
   do i =ii+1, N
    if (CC < UB(i)) CC = UB(i)
   enddo
   amax = sqrt(CC - UB(ii))
  else if (N == ii+1) then
   amax = sqrt(UB(ii+1)-UB(ii))
  else
   amax = 0.
 end if
 
 alim = sqrt(UB(1)+alpha*T*Ta-UB(ii))
 
 return
end subroutine aaa


!theta作成
subroutine TH(N, ii, amin, amax, alim, sigma, Ta, Tp, mu, UU, BB, theta)
 implicit none
 integer, intent(in) :: N, ii
 real*8, dimension(N), intent(in) :: UU, BB
 real*8, intent(in) :: amin, amax, alim, sigma, Ta, Tp, mu
 real*8, intent(out) :: theta
 real*8 :: thetaL, thetaM
 real*8, dimension(30) :: aL, aM
 integer :: pp
 real*8 :: ff !function
 
 thetaL = 0
 thetaM = 0
 
 if (alim > amin) then
  do pp = 1, 29
   aL(pp) = amin + (alim - amin)*dble(pp-1)/30.
   thetaL = thetaL + ff(sigma, BB(ii), BB(1), Ta, Tp, UU(ii), UU(1), mu, aL(pp))*(alim - amin)/30.
  enddo
 end if
 
 if (amax > amin) then
  do pp = 1, 29
   aM(pp) = amin + (amax - amin)*dble(pp-1)/30.
   thetaM = thetaM + ff(sigma, BB(ii), BB(1), Ta, Tp, UU(ii), UU(1), mu, aM(pp))*(amax - amin)/30.
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
 
 nis = 0.
 
 do p = 1, 59
  nis = nis + THT(p)*alpha*Tp/60./B
 enddo
 
 return
end subroutine NN


!電荷密度作成
subroutine RR(qq, nn, rho)
 implicit none
 real*8, dimension(8), intent(in) :: qq, nn
 real*8, intent(out) :: rho
 integer :: i
 
 rho = 0.
 do i = 1, 8
  rho = rho + qq(i)*nn(i)
 enddo
 
 return
end subroutine RR


!ニュートン法
subroutine NewtonPHI(N, HH, Hu, Hd, II, Iu, Id, HNew)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(N), intent(in) :: HH, Hu, Hd, II, Iu, Id
 real*8, dimension(N), intent(out) :: HNew
 integer :: i
 
! print *, Id
! print *, "  "
! print *, Iu
 
 do i = 1, N
  HNew(i) = HH(i) - (Hu(i) - Hd(i))/(Iu(i)-Id(i))*II(i)
 enddo
 
 return
end subroutine NewtonPHI

subroutine NewtonSIGMA(N, SS, RR, so, Num, sigmaN)
 implicit none
 integer, intent(in) :: N
 real*8, dimension(8), intent(in) :: SS
 real*8, dimension(3, N), intent(in) :: RR
 character(len=1), dimension(8), intent(in) :: so 
 real*8, dimension(3, N, 8), intent(in) :: Num
 real*8, dimension(8), intent(out) :: sigmaN
 integer :: s, b
 
 do s = 1, 8
  if(so(s) == "I") b = N
  if(so(s) == "M") b = 1
  sigmaN(s) = SS(s) - (Num(2, b, s) - Num(3, b, s))/(RR(2, b) - RR(3, b))*RR(1, b)
 enddo
 
 return
end subroutine NewtonSIGMA


