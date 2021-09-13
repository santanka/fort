program pds_J_kai_Temp

implicit none

!mainの文字
integer :: i, s, j !do文用
integer, parameter :: N = 999 !grid数
integer, parameter :: Z = 500 !速度空間grid数
double precision, parameter :: pi = 4.d0*atan(1.d0) !円周率
character(len=80) :: dummy !使用しない文字用

!parameterの導入
double precision, parameter :: ee = 1.602176634d-19 !電気素量
double precision, parameter :: cc = 2.99792458d+08 !光速
double precision, parameter :: mu0 = 1.25663706143592d-06 !真空の透磁率
double precision, parameter :: ep0 = 8.8541878128d-12 !真空の誘電率
double precision, parameter :: GG = 6.6743015d-11 !万有引力定数
double precision, parameter :: alpha = 12.d0 !vperpのlimit
double precision, parameter :: cc_rate = 1.d0 !光速の割合

!惑星のデータ
double precision, parameter :: Mdp = 1.16d+27 !惑星の磁気双極子モーメント
double precision, parameter :: Mp = 1.8982d+27 !惑星の質量
double precision, parameter :: Rp = 7.1492d+07 !惑星の半径
double precision, parameter :: omg = 1.75851813802955d-04 !惑星の自転角周波数
double precision, parameter :: Lp = 5.85d0 !対象の磁力線のL値
double precision, parameter :: lamN = acos(sqrt((1.d0+5.d5/Rp)/Lp)) !boundary:NのMLT
double precision, parameter :: lamS = -acos(sqrt((1.d0+5.d5/Rp)/Lp)) !boundary:SのMLT

!衛星のデータ(設定しない場合は値を0に)
double precision, parameter :: Ms = 8.931938d22 !衛星の質量
double precision, parameter :: Ls = 5.9d0 !衛星軌道のL値

!minファイルの導入
integer, parameter :: kind = 10
character(len=128) :: filemin = 'pds_J_kai_3_min.csv'
double precision :: cvn !収束値
double precision, dimension(N) :: lam !MLT
double precision, dimension(N) :: ss !磁力線上の座標
double precision, dimension(N) :: BB !磁束密度
double precision, dimension(N) :: Phi !静電ポテンシャル
double precision, dimension(N, kind) :: nd !数密度
double precision, dimension(N) :: rhod !電荷密度(分布関数)
double precision, dimension(N) :: rhop !電荷密度(Poisson方程式)
double precision, dimension(N) :: cvg !収束値
double precision, dimension(kind) :: sig !数密度
double precision, dimension(kind) :: Tperp !boundaryでのperp温度[eV]
double precision, dimension(kind) :: Tpara !boundaryでのpara温度[eV]
double precision, dimension(kind) :: chr !電荷[/ee]
double precision, dimension(kind) :: mass !質量
integer, dimension(kind) :: ijn !injection number

!計算
double precision, dimension(N) :: rr !惑星中心との距離
double precision, dimension(N) :: rs !衛星中心との距離
double precision, dimension(N) :: Pperpalli !イオン圧(perpendicular)
double precision, dimension(N) :: Pparaalli !イオン圧(parallel)
double precision, dimension(N) :: PPalli !イオン圧ecision, dimension(N) :: rs !衛星中心との距離
double precision, dimension(N, 3) :: VA !Alfven速度
double precision, dimension(kind, Z) :: mu !第1断熱不変量 
double precision, dimension(N, kind) :: UU !位置エネルギー
double precision, dimension(N, kind, Z) :: amin !amin
double precision, dimension(N, kind, Z) :: alim !alim
double precision, dimension(N, kind, Z) :: amax !amax
double precision, dimension(N, kind) :: Vpara !平均流速(parallel)
double precision, dimension(N, kind) :: Pperp !圧力(perpendicular)
double precision, dimension(N, kind) :: Ppara !圧力(parallel)
double precision, dimension(N, kind) :: PP !圧力
double precision, dimension(N) :: Pperpall !全圧(perpendicular)
double precision, dimension(N) :: Pparaall !全圧(parallel)
double precision, dimension(N) :: PPall !全圧
double precision, dimension(N) :: Pperpalle !電子圧(perpendicular)
double precision, dimension(N) :: Pparaalle !電子圧(parallel)
double precision, dimension(N) :: PPalle !電子圧
double precision, dimension(N) :: betaperp !β値(perpendicular)
double precision, dimension(N) :: betapara !β値(parallel)
double precision, dimension(N) :: beta !β値
double precision, dimension(N) :: betaperpi !イオンβ値(perpendicular)
double precision, dimension(N) :: betaparai !イオンβ値(parallel)
double precision, dimension(N) :: betai !イオンβ値
double precision, dimension(N) :: betaperpe !電子β値(perpendicular)
double precision, dimension(N) :: betaparae !電子β値(parallel)
double precision, dimension(N) :: betae !電子β値
double precision, dimension(N) :: border !IAWとKAWの境界値
double precision, dimension(N) :: lari !イオンLarmor半径
double precision, dimension(N) :: lare !電子Larmor半径
double precision, dimension(N) :: ski !イオン慣性長
double precision, dimension(N) :: ske !電子慣性長


!保存ファイル
character(len=128) :: fileall = 'pds_J_kai_3_all.csv'
!!lam(1), ss(2), BB(3), Phi(4), nd(5:kind+4), rhod(kind+5), rhop(kind+6), cvg(kind+7), VA(kind+8:kind+10),
!!Vpara(kind+11:2*kind+10), Pperp(2*kind+11:3*kind+10), Ppara(3*kind+11:4*kind+10), PP(4*kind+11:5*kind+10),
!!Pperpall(5*kind+11), Pparaall(5*kind+12), PPall(5*kind+13), Pperpalli(5*kind+14), Pparaalli(5*kind+15), PPalli(5*kind+16),
!!Pperpalle(5*kind+17), Pparaalle(5*kind+18), PPalle(5*kind+19), betaperp(5*kind+20), betapara(5*kind+21), beta(5*kind+22),
!!betaperpi(5*kind+23), betaparai(5*kind+24), betai(5*kind+25), betaperpe(5*kind+26), betaparae(5*kind+27), betae(5*kind+28),
!!border(5*kind+29), lari(5*kind+30), lare(5*kind+31), ski(5*kind+32), ske(5*kind+33)
92 format(1PE25.15E3, 82(',', 1PE25.15E3)) !5*kind+33

integer, parameter :: channel1 = 0 !座標固定, vperp vs vparaの分布関数(0:off, 1:on)
integer, parameter :: N1 = 129
integer, parameter :: channel2 = 0 !mu固定, MLT vs vparaの分布関数(0:off, 1:on)
integer, parameter :: Z2 = 1
character(len=128) :: fileoption = 'pds_E_kai_3_disfun_'


!/////データ抽出/////
open(20, file=filemin, action="read", status="old")
read(20, *) cvn
do i = 1, N
 read(20, *) lam(i), ss(i), BB(i), Phi(i), nd(i, :), rhod(i), rhop(i), cvg(i)
enddo !i
do s = 1, kind
 read(20, *) sig(s), Tperp(s), Tpara(s), chr(s), mass(s)
enddo !s
do s = 1, kind
 read(20, *) ijn(s)
enddo !s
close(20)

chr = chr*ee !次元を[C]に
Tperp = Tperp*ee !次元を[J]に
Tpara = Tpara*ee !次元を[J]に


!/////計算/////
!惑星中心との距離
rr = Lp*Rp*cos(lam)**2.d0

!衛星中心との距離
if(Ms /= 0) then
  rs = sqrt(rr**2.d0+(Ls*Rp)**2.d0-2.d0*rr*Ls*Rp*cos(lam))
endif

!ALfven速度
call Alf(N, kind, mu0, cc, BB, mass, nd, VA)

!断熱不変量
call AI(N, Z, kind, cc, cc_rate, mass, BB, ijn, mu)

!位置エネルギー
call EP(N, kind, GG, Mp, mass, rr, omg, lam, Ms, rs, chr, Phi, UU)

!Accessibility
call access(N, Z, kind, cc_rate, cc, mass, UU, mu, BB, ijn, amin, alim, amax)

!平均流速
call ryusoku(N, Z, kind, UU, nd, mu, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, Vpara)

!圧力
!perpendicular
call pressureperp(N, Z, kind, mu, BB, Tperp, Tpara, sig, amin, amax, ijn, Pperp)

!parallel
call pressurepara(N, Z, kind, UU, Vpara, mu, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, Ppara)

!全圧他
call pressure(N, kind, chr, Pperp, Ppara, PP, Pperpall, Pparaall, PPall, Pperpalli, Pparaalli, PPalli, &
& Pperpalle, Pparaalle, PPalle)

!β値
betaperp = 2.d0*mu0*Pperpall/BB**2.d0
betapara = 2.d0*mu0*Pparaall/BB**2.d0
beta = 2.d0*mu0*PPall/BB**2.d0
betaperpi = 2.d0*mu0*Pperpalli/BB**2.d0
betaparai = 2.d0*mu0*Pparaalli/BB**2.d0
betai = 2.d0*mu0*PPalli/BB**2.d0
betaperpe = 2.d0*mu0*Pperpalle/BB**2.d0
betaparae = 2.d0*mu0*Pparaalle/BB**2.d0
betae = 2.d0*mu0*PPalle/BB**2.d0

!IAWとKAWの境界値&Larmor半径&慣性長
call bor(N, kind, ep0, cc, BB, Pperpalli, Pperpalle, nd, mass, chr, border, lari, lare, ski, ske)


!ファイル化処理
open(50, file = fileall)
do i = 1, N
 write(50, 92) lam(i), ss(i), BB(i), Phi(i), nd(i, :), rhod(i), rhop(i), cvg(i), VA(i, :), Vpara(i, :), Pperp(i, :), Ppara(i, :), &
& PP(i, :), Pperpall(i), Pparaall(i), PPall(i), Pperpalli(i), Pparaalli(i), PPalli(i), Pperpalle(i), Pparaalle(i), PPalle(i), &
& betaperp(i), betapara(i), beta(i), betaperpi(i), betaparai(i), betai(i), betaperpe(i), betaparae(i), betae(i), border(i), &
& lari(i), lare(i), ski(i), ske(i)
enddo !i
close(50)



!座標固定, vperp vs vparaの分布関数
if(channel1 == 1) then
  call disfun1(N, Z, kind, N1, ee, UU, nd, mu, lam, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, fileoption)
endif

!mu固定, MLT vs vparaの分布関数
if(channel2 == 1) then
  call disfun2(N, Z, kind, Z2, ee, UU, nd, mu, lam, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, fileoption)
endif


end program pds_J_kai_Temp






!subroutine & function

!Alfven速度
subroutine Alf(N, kind, mu0, cc, BB, mass, nd, VA)
 implicit none
 integer, intent(in) :: N, kind
 double precision, intent(in) :: mu0, cc
 double precision, dimension(N), intent(in) :: BB
 double precision, dimension(kind), intent(in) :: mass
 double precision, dimension(N, kind), intent(in) :: nd
 double precision, dimension(N, 3), intent(out) :: VA
 
 integer :: s
 double precision, dimension(N) :: ww
 
 ww = 0.d0
 
 do s = 1, kind
  ww = ww + mass(s)*nd(:, s)
 enddo !s
 
 VA(:, 1) = BB/sqrt(mu0*ww)
 VA(:, 2) = VA(:, 1)/sqrt(1.d0 + (VA(:, 1)/cc)**2.d0)
 VA(:, 3) = VA(:, 2)/cc
 
 return
 
end subroutine Alf


!断熱不変量
subroutine AI(N, Z, kind, cc, cc_rate, mass, BB, ijn, mu)
 implicit none
 integer, intent(in) :: N, Z, kind
 double precision, intent(in) :: cc, cc_rate
 double precision, dimension(kind), intent(in) :: mass
 double precision, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(kind, Z), intent(out) :: mu
 
 integer :: j
 
 do j = 1, Z
  if(j == 1) then
    mu(:, j) = 0.d0
   else if(j /= 1) then
    mu(:, j) = 1.d-30*(mass*(cc_rate*cc)**2.d0/2.d0/BB(ijn)/1.d-30)**(dble(j-2)/dble(Z-2))
  endif
 enddo !j
 
 return
 
end subroutine AI


!位置エネルギー
subroutine EP(N, kind, GG, Mp, mass, rr, omg, lam, Ms, rs, chr, Phi, UU)
 implicit none
 integer, intent(in) :: N, kind
 double precision, intent(in) :: GG, Mp, omg, Ms
 double precision, dimension(kind), intent(in) :: mass, chr
 double precision, dimension(N), intent(in) :: rr, lam, rs
 double precision, dimension(N), intent(in) :: Phi
 double precision, dimension(N, kind), intent(out) :: UU
 
 integer :: i
 
 do i = 1, N
  UU(i, :) = -GG*Mp*mass/rr(i) !惑星重力
  UU(i, :) = UU(i, :) - mass*(omg*rr(i)*cos(lam(i)))**2.d0/2.d0 !惑星遠心力
  if(Ms /= 0) UU(i, :) = UU(i, :) - GG*Ms*mass/rs(i) !衛星重力
  UU(i, :) = UU(i, :) + chr*Phi(i) !クーロン力
 enddo !i
 
 return
 
end subroutine EP


!Accessibility
subroutine access(N, Z, kind, cc_rate, cc, mass, UU, mu, BB, ijn, amin, alim, amax)
 implicit none
 integer, intent(in) :: N, Z, kind
 double precision, intent(in) :: cc_rate, cc
 double precision, dimension(kind), intent(in) :: mass
 double precision, dimension(N, kind), intent(in) :: UU
 double precision, dimension(kind, Z), intent(in) :: mu
 double precision, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(N, kind, Z), intent(out) :: amin, alim, amax
 
 integer :: i, s, j, k, t
 double precision, dimension(N, kind, Z) :: EE !UU+mu*BB
 double precision :: WW
 
 
 !エネルギーの和(UU+mu*BB)
 do i = 1, N
  do s = 1, kind
   EE(i, s, :) = UU(i, s) + mu(s, :)*BB(i)
  enddo !s
 enddo !i
 
 !amin
 do s = 1, kind
  do i = 1, N
   if(i < ijn(s)) then
     do j = 1, Z
      t = i+1
      do k = i+1, ijn(s)
       if(EE(k, s, j) > EE(t, s, j)) t = k
      enddo !k
      if(EE(i, s, j) > EE(t, s, j)) then
        amin(i, s, j) = sqrt(EE(i, s, j) - EE(ijn(s), s, j))
       else if(EE(i, s, j) <= EE(t, s, j)) then
        amin(i, s, j) = sqrt(EE(t, s, j) - EE(ijn(s), s, j))
      endif
     enddo !j
     
    else if(i == ijn(s)) then
     amin(ijn(s), s, :) = 0.d0
     
    else if(i > ijn(s)) then
     do j = 1, Z
      t = ijn(s)
      do k = ijn(s), i-1
       if(EE(k, s, j) > EE(t, s, j)) t = k
      enddo !k
      if(EE(i, s, j) > EE(t, s, j)) then
        amin(i, s, j) = sqrt(EE(i, s, j) - EE(ijn(s), s, j))
       else if(EE(i, s, j) <= EE(t, s, j)) then
        amin(i, s, j) = sqrt(EE(t, s, j) - EE(ijn(s), s, j))
      endif
     enddo !j
   endif
  enddo !i
 enddo !s
 
 !alim
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    if(EE(i, s, j) < EE(ijn(s), s, j)) then
      alim(i, s, j) = sqrt(EE(i, s, j)-EE(ijn(s), s, j)+mass(s)/2.d0*(cc_rate*cc)**2.d0)
     else if(EE(i, s, j) >= EE(ijn(s), s, j)) then
      alim(i, s, j) = sqrt(mass(s)/2.d0*(cc_rate*cc)**2.d0)
    endif
   enddo !j
  enddo !s
 enddo !i
 
 !amax
 do s = 1, kind
  do i = 1, N
   if((i <= ijn(s) .and. ijn(s) == N) .or. (i < ijn(s) .and. ijn(s) /= N)) then
     if(i == 1) then
       amax(1, s, :) = 0.d0
      else if(i /= 1) then
       do j = 1, Z
        t = 1
        do k = 1, i-1
         if(EE(k, s, j) > EE(t, s, j)) t = k
        enddo !k
        WW = EE(t, s, j) - EE(ijn(s), s, j)
        if(WW <= 0.d0) then
          amax(i, s, j) = 0.d0
         else if(WW > 0.d0 .and. WW < mass(s)/2.d0*(cc_rate*cc)**2.d0) then
          amax(i, s, j) = sqrt(WW)
         else if(WW >= mass(s)/2.d0*(cc_rate*cc)**2.d0) then
          amax(i, s, j) = sqrt(mass(s)/2.d0)*cc_rate*cc
        endif
       enddo !j
     endif
     
    else if((i >= ijn(s) .and. ijn(s) == 1) .or. (i > ijn(s) .and. ijn(s) /= 1)) then
     if(i == N) then
       amax(N, s, :) = 0.d0
      else if(i /= N) then
       do j = 1, Z
        t = i+1
        do k = i+1, N
         if(EE(k, s, j) > EE(t, s, j)) t = k
        enddo !k
        WW = EE(t, s, j) - EE(ijn(s), s, j)
        if(WW <= 0.d0) then
          amax(i, s, j) = 0.d0
         else if(WW > 0.d0 .and. WW < mass(s)/2.d0*(cc_rate*cc)**2.d0) then
          amax(i, s, j) = sqrt(WW)
         else if(WW >= mass(s)/2.d0*(cc_rate*cc)**2.d0) then
          amax(i, s, j) = sqrt(mass(s)/2.d0)*cc_rate*cc
        endif
       enddo !j
     endif
     
    else if(i == ijn(s) .and. ijn(s) /= 1 .and. ijn(s) /= N) then !例外処理
     amax(i, s, :) = alim(i, s, :)
     
   endif
  enddo !i
 enddo !s
 
 return
 
end subroutine access


!平均流速
subroutine ryusoku(N, Z, kind, UU, nd, mu, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, Vpara)
 implicit none
 integer, intent(in) :: N, Z, kind
 double precision, dimension(N, kind), intent(in) :: UU, nd
 double precision, dimension(kind, Z), intent(in) :: mu
 double precision, dimension(N), intent(in) :: BB
 double precision, dimension(kind), intent(in) :: mass, Tpara, Tperp, sig
 double precision, dimension(N, kind, Z), intent(in) :: amin, alim, amax
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(N, kind), intent(out) :: Vpara
 
 double precision, dimension(N, kind, Z) :: EE
 double precision, dimension(Z/2) :: aL, aM
 double precision, dimension(Z) :: TH
 double precision :: a, ww, thetaL, thetaM, nnn, LL1, LL2, LL3, LL4, MM1, MM2, MM3, MM4
 integer :: i, s, j, p, check
 double precision, parameter :: pi = 4.d0*atan(1.d0)
 
 
 !エネルギーの和(UU+mu*BB)
 do i = 1, N
  do s = 1, kind
   EE(i, s, :) = UU(i, s) + mu(s, :)*BB(i)
  enddo !s
 enddo !i
 
 !NaNの前処理
 a = -1.d0
 
 
 do i = 1, N
  do s = 1, kind
   if(nd(i, s) == 0.d0) then
     Vpara(i, s) = sqrt(a) !NaNを代入
    else if(nd(i, s) /= 0.d0) then
     TH = 0.d0
     if(i < ijn(s)) then
       ww = -1.d0
      else if(i >= ijn(s)) then
       ww = 1.d0
     endif
     do j = 1, Z
      check = 0 !amaxの積分をするかのチェック
      if(amin(i, s, j) == 0.d0) then
        aL(1) = 0.d0
        do p = 2, Z/2
         aL(p) = ww * alim(i, s, j)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
        enddo !p
        if(amax(i, s, j) /= 0.d0) then
          check = 1 !amaxの積分を実行
          aM(1) = 0.d0
          do p = 2, Z/2
           aM(p) = -ww * amax(i, s, j)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
          enddo !p
         else if(amax(i, s, j) == 0.d0) then
          aM = 0.d0
        endif
        
       else if(amin(i, s, j) /= 0.d0) then
        do p = 1, Z/2
         aL(p) = ww * amin(i, s, j)*((alim(i, s, j)/amin(i, s, j))**(dble(p-1)/dble(Z/2-1)))
        enddo !p
        if(amax(i, s, j) > amin(i, s, j)) then
          check = 1 !amaxの積分を実行
          do p = 1, Z/2
           aM(p) = -ww * amin(i, s, j)*((amax(i, s, j)/amin(i, s, j))**(dble(p-1)/dble(Z/2-1)))
          enddo !p
         else if(amax(i, s, j) <= amin(i, s, j)) then
          aM = 0.d0
        endif
      endif
      
      
      thetaL = 0.d0
      thetaM = 0.d0
      
      do p = 1, Z/2-1
       LL1 = sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aL(p)**2.d0))
       LL2 = sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aL(p+1)**2.d0))
       LL3 = exp(-(aL(p)**2.d0)/Tpara(s))
       LL4 = exp(-(aL(p+1)**2.d0)/Tpara(s))
       if(LL1 > 0.d0 .and. LL2 > 0.d0 .and. LL3 >= 0.d0 .and. LL4 >= 0.d0) then
         thetaL = thetaL + aL(p+1)/abs(aL(p+1))*(LL1*LL3 + LL2*LL4) / 2.d0 * abs(aL(p+1)-aL(p))
       endif
       MM1 = sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aM(p)**2.d0))
       MM2 = sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aM(p+1)**2.d0))
       MM3 = exp(-(aM(p)**2.d0)/Tpara(s))
       MM4 = exp(-(aM(p+1)**2.d0)/Tpara(s))
       if(check == 1 .and. MM1 > 0.d0 .and. MM2 > 0.d0 .and. MM3 >= 0.d0 .and. MM4 >= 0.d0) then
         thetaM = thetaM + aM(p+1)/abs(aM(p+1))*(MM1*MM3 + MM2*MM4) / 2.d0 * abs(aM(p+1)-aM(p))
       endif
      enddo !p
      
      TH(j) = (thetaL + thetaM)*exp(-mu(s, j)*BB(ijn(s))/Tperp(s))
      
     enddo !j
     
     nnn = 0.d0
     
     do j = 1, Z-1
      nnn = nnn + (TH(j) + TH(j+1))/2.d0*(mu(s, j+1)-mu(s, j))
     enddo !j
     
     Vpara(i, s) = sqrt(2.d0/pi/mass(s)/Tpara(s))*sig(s)*BB(ijn(s))/Tperp(s)/nd(i, s) * nnn
   endif
  enddo !s
 enddo !i
 
 return
 
end subroutine ryusoku


!perpendiicular圧力
subroutine pressureperp(N, Z, kind, mu, BB, Tperp, Tpara, sig, amin, amax, ijn, Pperp)
 implicit none
 integer, intent(in) :: N, Z, kind
 double precision, dimension(kind, Z), intent(in) :: mu
 double precision, dimension(N), intent(in) :: BB
 double precision, dimension(kind), intent(in) :: Tperp, Tpara, sig
 double precision, dimension(N, kind, Z), intent(in) :: amin, amax
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(N, kind), intent(out) :: Pperp
 
 integer :: i, s, j
 double precision, dimension(N, kind, Z) :: TL, TM, TT
 
 
 do s = 1, kind
  TL(:, s, :) = 1.d0 - erf(amin(:, s, :)/sqrt(Tpara(s)))
 enddo !s
 
 do i = 1, N
  do s = 1, kind
   do j = 1, Z
    if(amax(i, s, j) > amin(i, s, j)) then
      TM(i, s, j) = erf(amax(i, s, j)/sqrt(Tpara(s))) - erf(amin(i, s, j)/sqrt(Tpara(s)))
     else if(amax(i, s, j) <= amin(i, s, j)) then
      TM(i, s, j) = 0.d0
    endif
   enddo !j
  enddo !s
 enddo !i
 
 do s = 1, kind
  do j = 1, Z
   TT(:, s, j) = (TL(:, s, j) + TM(:, s, j)) * mu(s, j) * exp(-mu(s, j)*BB(ijn(s))/Tperp(s))
  enddo !j
 enddo !s
 
 Pperp = 0.d0
 
 do s = 1, kind
  do j = 1, Z-1
   Pperp(:, s) = Pperp(:, s) + (TT(:, s, j)+TT(:, s, j+1))/2.d0 * (mu(s, j+1)-mu(s, j))
  enddo !j
  Pperp(:, s) = sig(s)*BB(ijn(s))*BB/Tperp(s)*Pperp(:, s)
 enddo !s
 
 return
 
end subroutine pressureperp


!parallel圧力
subroutine pressurepara(N, Z, kind, UU, Vpara, mu, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, Ppara)
 implicit none
 integer, intent(in) :: N, Z, kind
 double precision, dimension(N, kind), intent(in) :: UU, Vpara
 double precision, dimension(kind, Z), intent(in) :: mu
 double precision, dimension(N), intent(in) :: BB
 double precision, dimension(kind), intent(in) :: mass, Tpara, Tperp, sig
 double precision, dimension(N, kind, Z), intent(in) :: amin, alim, amax
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(N, kind), intent(out) :: Ppara
 
 double precision, dimension(N, kind, Z) :: EE
 double precision, dimension(Z/2) :: aL, aM
 double precision, dimension(Z) :: TH
 double precision :: a, ww, thetaL, thetaM, nnn, LL1, LL2, LL3, LL4, MM1, MM2, MM3, MM4
 integer :: i, s, j, p, check
 double precision, parameter :: pi = 4.d0*atan(1.d0)
 
 
 !エネルギーの和(UU+mu*BB)
 do i = 1, N
  do s = 1, kind
   EE(i, s, :) = UU(i, s) + mu(s, :)*BB(i)
  enddo !s
 enddo !i
 
 
 do i = 1, N
  do s = 1, kind
   TH = 0.d0
   if(i < ijn(s)) then
     ww = -1.d0
    else if(i >= ijn(s)) then
     ww = 1.d0
   endif
   do j = 1, Z
    check = 0 !amaxの積分をするかのチェック
    if(amin(i, s, j) == 0.d0) then
      aL(1) = 0.d0
      do p = 2, Z/2
       aL(p) = ww * alim(i, s, j)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
      enddo !p
      if(amax(i, s, j) /= 0.d0) then
        check = 1 !amaxの積分を実行
        aM(1) = 0.d0
        do p = 2, Z/2
         aM(p) = -ww * amax(i, s, j)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
        enddo !p
       else if(amax(i, s, j) == 0.d0) then
        aM = 0.d0
      endif
      
     else if(amin(i, s, j) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = ww * amin(i, s, j)*((alim(i, s, j)/amin(i, s, j))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      if(amax(i, s, j) > amin(i, s, j)) then
        check = 1 !amaxの積分を実行
        do p = 1, Z/2
         aM(p) = -ww * amin(i, s, j)*((amax(i, s, j)/amin(i, s, j))**(dble(p-1)/dble(Z/2-1)))
        enddo !p
       else if(amax(i, s, j) <= amin(i, s, j)) then
        aM = 0.d0
      endif
    endif
    
    
    thetaL = 0.d0
    thetaM = 0.d0
    
    do p = 1, Z/2-1
     LL1 = (sqrt(2.d0/mass(s))*aL(p+1)/abs(aL(p+1))*sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aL(p)**2.d0)) - Vpara(i, s))**2.d0
     LL2 = (sqrt(2.d0/mass(s))*aL(p+1)/abs(aL(p+1))*sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aL(p+1)**2.d0)) - Vpara(i, s))**2.d0
     LL3 = exp(-(aL(p)**2.d0)/Tpara(s))
     LL4 = exp(-(aL(p+1)**2.d0)/Tpara(s))
     if(LL1 > 0.d0 .and. LL2 > 0.d0 .and. LL3 >= 0.d0 .and. LL4 >= 0.d0) then
       thetaL = thetaL + (LL1*LL3 + LL2*LL4) / 2.d0 * abs(aL(p+1)-aL(p))
     endif
     MM1 = (sqrt(2.d0/mass(s))*aM(p+1)/abs(aM(p+1))*sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aM(p)**2.d0)) - Vpara(i, s))**2.d0
     MM2 = (sqrt(2.d0/mass(s))*aM(p+1)/abs(aM(p+1))*sqrt(EE(ijn(s), s, j)-EE(i, s, j)+(aM(p+1)**2.d0)) - Vpara(i, s))**2.d0
     MM3 = exp(-(aM(p)**2.d0)/Tpara(s))
     MM4 = exp(-(aM(p+1)**2.d0)/Tpara(s))
     if(check == 1 .and. MM1 > 0.d0 .and. MM2 > 0.d0 .and. MM3 >= 0.d0 .and. MM4 >= 0.d0) then
       thetaM = thetaM + (MM1*MM3 + MM2*MM4) / 2.d0 * abs(aM(p+1)-aM(p))
     endif
    enddo !p
    
    TH(j) = (thetaL + thetaM)*exp(-mu(s, j)*BB(ijn(s))/Tperp(s))
    
   enddo !j
   
   nnn = 0.d0
   
   do j = 1, Z-1
    nnn = nnn + (TH(j) + TH(j+1))/2.d0*(mu(s, j+1)-mu(s, j))
   enddo !j
   
   Ppara(i, s) = sig(s)*BB(ijn(s))*mass(s)/Tperp(s)/sqrt(pi*Tpara(s)) * nnn
  enddo !s
 enddo !i
 
 return
 
end subroutine pressurepara


!圧力
subroutine pressure(N, kind, chr, Pperp, Ppara, PP, Pperpall, Pparaall, PPall, Pperpalli, Pparaalli, PPalli, &
& Pperpalle, Pparaalle, PPalle)
 implicit none
 integer, intent(in) :: N, kind
 double precision, dimension(kind), intent(in) :: chr
 double precision, dimension(N, kind), intent(in) :: Pperp, Ppara
 double precision, dimension(N, kind), intent(out) :: PP
 double precision, dimension(N), intent(out) :: Pperpall, Pparaall, PPall, Pperpalli, Pparaalli, PPalli
 double precision, dimension(N), intent(out) :: Pperpalle, Pparaalle, PPalle
 
 integer :: s
 
 
 PP = (2.d0*Pperp + Ppara)/3.d0 !圧力
 
 Pperpall = 0.d0 !全圧(perp)
 Pparaall = 0.d0 !全圧(para)
 PPall = 0.d0 !全圧
 Pperpalli = 0.d0 !イオン圧(perp)
 Pparaalli = 0.d0 !イオン圧(para)
 PPalli = 0.d0 !イオン圧
 Pperpalle = 0.d0 !電子圧(perp)
 Pparaalle = 0.d0 !電子圧(para)
 PPalle = 0.d0 !電子圧
 
 
 do s = 1, kind 
  Pperpall = Pperpall + Pperp(:, s)
  Pparaall = Pparaall + Ppara(:, s)
  PPall = PPall + PP(:, s)
  if(chr(s) > 0.d0) then
    Pperpalli = Pperpalli + Pperp(:, s)
    Pparaalli = Pparaalli + Ppara(:, s)
    PPalli = PPalli + PP(:, s)
   else if(chr(s) < 0.d0) then
    Pperpalle = Pperpalle + Pperp(:, s)
    Pparaalle = Pparaalle + Ppara(:, s)
    PPalle = PPalle + PP(:, s)
  endif
 enddo !s
 
 return
 
end subroutine pressure


!IAWとKAWの境界値&Larmor半径&慣性長
subroutine bor(N, kind, ep0, cc, BB, Pperpalli, Pperpalle, nd, mass, chr, border, lari, lare, ski, ske)
 implicit none
 integer, intent(in) :: N, kind
 double precision, intent(in) :: ep0, cc
 double precision, dimension(N), intent(in) :: BB, Pperpalli, Pperpalle
 double precision, dimension(N, kind), intent(in) :: nd
 double precision, dimension(kind), intent(in) :: mass, chr
 double precision, dimension(N), intent(out) :: border, lari, lare, ski, ske
 
 double precision, dimension(N) :: mp, mn, np, nn, chp, chn
 integer :: s
 
 border = 0.d0
 mp = 0.d0
 mn = 0.d0
 np = 0.d0
 nn = 0.d0
 chp = 0.d0
 chn = 0.d0
 
 do s = 1, kind
  if(chr(s) > 0.d0) then
    mp = mp + mass(s)*nd(:, s)
    np = np + nd(:, s)
    chp = chp + chr(s)*nd(:, s)
   else if(chr(s) < 0.d0) then
    mn = mn + mass(s)*nd(:, s)
    nn = nn + nd(:, s)
    chn = chn + chr(s)*nd(:, s)
  endif
 enddo !s
 
 border = mn/nn / (mp/np)
 
 lari = sqrt(2.d0 * mp/np * Pperpalli / np) / (chp/np) / BB
 lare = sqrt(2.d0 * mn/nn * Pperpalle / nn) / abs(chn/nn) / BB
 
 ski = cc / (chp/np) * sqrt(ep0 * mp/np / np)
 ske = cc / abs(chn/nn) * sqrt(ep0 * mn/nn / nn)
 
 return
 
end subroutine bor


!座標固定, vperp vs vparaの分布関数
subroutine disfun1(N, Z, kind, N1, ee, UU, nd, mu, lam, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, fileoption)
  implicit none
  integer, intent(in) :: N, Z, kind, N1
  double precision, intent(in) :: ee
  double precision, dimension(N, kind), intent(in) :: UU, nd
  double precision, dimension(kind, Z), intent(in) :: mu
  double precision, dimension(N), intent(in) :: lam, BB
  double precision, dimension(kind), intent(in) :: mass, Tpara, Tperp, sig
  double precision, dimension(N, kind, Z), intent(in) :: amin, alim, amax
  integer, dimension(kind), intent(in) :: ijn
  character(len=128), intent(in) :: fileoption
  
  double precision, dimension(Z/2) :: aL, aM
  double precision :: a, ww, vpara, vperp, vparab, vperpb, kpara, kperp, kall, ffb, ffi
  integer :: i, s, j, p, check
  character(len=128) :: filename, filename1, filename2
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  74 format(1PE25.15E3, 9(',', 1PE25.15E3))
  
  
  do s = 1, kind
   !ファイル名作成
   write(filename1, *) s
   write(filename2, *) lam(N1) /pi*180.
   filename = trim(adjustl(fileoption))//'s='//trim(adjustl(filename1))//'_MLAT='//trim(adjustl(filename2))//'.csv'
   print *, filename
   
   open(62, file=filename)
   if(N1 < ijn(s)) then
     ww = -1.d0
    else if(N1 >= ijn(s)) then
     ww = 1.d0
   endif
   do j = 1, Z
    if(mod(j, 3) == 1) then
      check = 0 !amaxチェック
      if(amin(N1, s, j) == 0.d0) then
        aL(1) = 0.d0
        do p = 2, Z/2
         aL(p) = ww * alim(N1, s, j)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
        enddo !p
        if(amax(N1, s, j) /= 0.d0) then
          check = 1 !amax実行
          aM(1) = 0.d0
          do p = 2, Z/2
           aM(p) = -ww * amax(N1, s, j)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
          enddo !p
         else if(amax(N1, s, j) == 0.d0) then
          aM = 0.d0
        endif
        
       else if(amin(N1, s, j) /= 0.d0) then
        do p = 1, Z/2
         aL(p) = ww * amin(N1, s, j)*((alim(N1, s, j)/amin(N1, s, j))**(dble(p-1)/dble(Z/2-1)))
        enddo !p
        if(amax(N1, s, j) > amin(N1, s, j)) then
          check = 1 !amaxの積分を実行
          do p = 1, Z/2
           aM(p) = -ww * amin(N1, s, j)*((amax(N1, s, j)/amin(N1, s, j))**(dble(p-1)/dble(Z/2-1)))
          enddo !p
         else if(amax(N1, s, j) <= amin(N1, s, j)) then
          aM = 0.d0
        endif
      endif
      
      do p = 1, Z/2
       if(mod(p, 3) == 2 .and. nd(N1, s) /= 0.d0) then
         ffb = sig(s)/nd(N1, s)*(mass(s)/2.d0/pi)**(3.d0/2.d0)/Tperp(s)/sqrt(Tpara(s)) &
 & *exp(-mu(s, j)*BB(ijn(s))/Tperp(s))*exp(-aL(p)**2.d0/Tpara(s))
         vpara = sqrt(2.d0/mass(s))*sqrt(UU(ijn(s), s)-UU(N1, s)+mu(s, j)*(BB(ijn(s))-BB(N1))+aL(p)**2.d0)*aL(p)/abs(aL(p))
         vperp = sqrt(2.d0*BB(N1)*mu(s, j)/mass(s))
         vparab = sqrt(2.d0/mass(s))*aL(p)
         vperpb = sqrt(2.d0*BB(ijn(s))*mu(s, j)/mass(s))
         kpara = mass(s)*vpara**2.d0 /2.d0 /ee
         kperp = mass(s)*vperp**2.d0 /2.d0 /ee
         kall = (kpara+kperp*2.d0)/3.d0
         
         if(ffb /= 0.d0) then
           if(vparab /= 0.d0) then
             ffi = ffb * BB(ijn(s))/BB(N1) * abs(vpara/vparab)
           else
             ffi = 0d0
           end if
 
           write(62, 74) lam(N1), vpara, vperp, vparab, vperpb, kpara, kperp, kall, ffb, ffi
         endif
 
       endif
      enddo !p
      
      if(check == 1) then
        do p = 1, Z/2
         if(mod(p, 3) == 2 .and. nd(N1, s) /= 0.d0) then
           ffb = sig(s)/nd(N1, s)*(mass(s)/2.d0/pi)**(3.d0/2.d0)/Tperp(s)/sqrt(Tpara(s)) &
 & *exp(-mu(s, j)*BB(ijn(s))/Tperp(s))*exp(-aM(p)**2.d0/Tpara(s))
           vpara = sqrt(2.d0/mass(s))*sqrt(UU(ijn(s), s)-UU(N1, s)+mu(s, j)*(BB(ijn(s))-BB(N1))+aM(p)**2.d0)*aM(p)/abs(aM(p))
           vperp = sqrt(2.d0*BB(N1)*mu(s, j)/mass(s))
           vparab = sqrt(2.d0/mass(s))*aM(p)
           vperpb = sqrt(2.d0*BB(ijn(s))*mu(s, j)/mass(s))
           kpara = mass(s)*vpara**2.d0 /2.d0 /ee
           kperp = mass(s)*vperp**2.d0 /2.d0 /ee
           kall = (kpara+kperp*2.d0)/3.d0
           
           if(ffb /= 0.d0) then
             if(vparab /= 0.d0) then
               ffi = ffb * BB(ijn(s))/BB(N1) * abs(vpara/vparab)
             else
               ffi = 0d0
             end if
   
             write(62, 74) lam(N1), vpara, vperp, vparab, vperpb, kpara, kperp, kall, ffb, ffi
           endif
           
         endif
        enddo !p
      endif
    endif
   enddo !j
   
   close(62)
   
  enddo !s
  
  return
  
 end subroutine disfun1
 
 
 !mu固定, MLT vs vparaの分布関数
 subroutine disfun2(N, Z, kind, Z2, ee, UU, nd, mu, lam, BB, mass, Tpara, Tperp, sig, amin, alim, amax, ijn, fileoption)
  implicit none
  integer, intent(in) :: N, Z, kind, Z2
  double precision, intent(in) :: ee
  double precision, dimension(N, kind), intent(in) :: UU, nd
  double precision, dimension(kind, Z), intent(in) :: mu
  double precision, dimension(N), intent(in) :: lam, BB
  double precision, dimension(kind), intent(in) :: mass, Tpara, Tperp, sig
  double precision, dimension(N, kind, Z), intent(in) :: amin, alim, amax
  integer, dimension(kind), intent(in) :: ijn
  character(len=128), intent(in) :: fileoption
  
  double precision, dimension(Z/2) :: aL, aM
  double precision :: a, ww, vpara, vperp, vparab, vperpb, kpara, kperp, kall, ffb, ffi
  integer :: i, s, j, p, check
  character(len=128) :: filename, filename1, filename2
  double precision, parameter :: pi = 4.d0*atan(1.d0)
  74 format(1PE25.15E3, 9(',', 1PE25.15E3))
  
  
  do s = 1, kind
   !ファイル名作成
   write(filename1, *) s
   write(filename2, *) mu(s, Z2)
   filename = trim(adjustl(fileoption))//'s='//trim(adjustl(filename1))//'_mu='//trim(adjustl(filename2))//'.csv'
   print *, filename
   
   open(62, file=filename)
   
   do i = 1, N
    if(i < ijn(s)) then
      ww = -1.d0
     else if(i >= ijn(s)) then
      ww = 1.d0
    endif
    
    check = 0 !amaxチェック
    if(amin(i, s, Z2) == 0.d0) then
      aL(1) = 0.d0
      do p = 2, Z/2
       aL(p) = ww * alim(i, s, Z2)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
      enddo !p
      if(amax(i, s, Z2) /= 0.d0) then
        check = 1 !amax実行
        aM(1) = 0.d0
        do p = 2, Z/2
         aM(p) = -ww * amax(i, s, Z2)*1.d-8*(1.d8**(dble(p-2)/dble(Z/2-2)))
        enddo !p
       else if(amax(i, s, Z2) == 0.d0) then
        aM = 0.d0
      endif
      
     else if(amin(i, s, Z2) /= 0.d0) then
      do p = 1, Z/2
       aL(p) = ww * amin(i, s, Z2)*((alim(i, s, Z2)/amin(i, s, Z2))**(dble(p-1)/dble(Z/2-1)))
      enddo !p
      if(amax(i, s, Z2) > amin(i, s, Z2)) then
        check = 1 !amaxの積分を実行
        do p = 1, Z/2
         aM(p) = -ww * amin(i, s, Z2)*((amax(i, s, Z2)/amin(i, s, Z2))**(dble(p-1)/dble(Z/2-1)))
        enddo !p
       else if(amax(i, s, Z2) <= amin(i, s, Z2)) then
        aM = 0.d0
      endif
    endif
    
    do p = 1, Z/2
     if(mod(p, 3) == 2 .and. nd(i, s) /= 0.d0) then
       ffb = sig(s)/nd(i, s)*(mass(s)/2.d0/pi)**(3.d0/2.d0)/Tperp(s)/sqrt(Tpara(s)) &
 & *exp(-mu(s, Z2)*BB(ijn(s))/Tperp(s))*exp(-aL(p)**2.d0/Tpara(s))
       vpara = sqrt(2.d0/mass(s))*sqrt(UU(ijn(s), s)-UU(i, s)+mu(s, Z2)*(BB(ijn(s))-BB(i))+aL(p)**2.d0)*aL(p)/abs(aL(p))
       vperp = sqrt(2.d0*BB(i)*mu(s, Z2)/mass(s))
       vparab = sqrt(2.d0/mass(s))*aL(p)
       vperpb = sqrt(2.d0*BB(ijn(s))*mu(s, Z2)/mass(s))
       kpara = mass(s)*vpara**2.d0 /2.d0 /ee
       kperp = mass(s)*vperp**2.d0 /2.d0 /ee
       kall = (kpara+kperp*2.d0)/3.d0
       
       if(ffb /= 0.d0) then
         if(vparab /= 0.d0) then
           ffi = ffb * BB(ijn(s))/BB(i) * abs(vpara/vparab)
         else
           ffi = 0d0
         end if
 
         write(62, 74) lam(i), vpara, vperp, vparab, vperpb, kpara, kperp, kall, ffb, ffi
       endif
 
     endif
    enddo !p
    
    if(check == 1) then
      do p = 1, Z/2
       if(mod(p, 3) == 2 .and. nd(i, s) /= 0.d0) then
         ffb = sig(s)/nd(i, s)*(mass(s)/2.d0/pi)**(3.d0/2.d0)/Tperp(s)/sqrt(Tpara(s)) &
 & *exp(-mu(s, Z2)*BB(ijn(s))/Tperp(s))*exp(-aM(p)**2.d0/Tpara(s))
         vpara = sqrt(2.d0/mass(s))*sqrt(UU(ijn(s), s)-UU(i, s)+mu(s, Z2)*(BB(ijn(s))-BB(i))+aM(p)**2.d0)*aM(p)/abs(aM(p))
         vperp = sqrt(2.d0*BB(i)*mu(s, Z2)/mass(s))
         vparab = sqrt(2.d0/mass(s))*aM(p)
         vperpb = sqrt(2.d0*BB(ijn(s))*mu(s, Z2)/mass(s))
         kpara = mass(s)*vpara**2.d0 /2.d0 /ee
         kperp = mass(s)*vperp**2.d0 /2.d0 /ee
         kall = (kpara+kperp*2.d0)/3.d0
         
         if(ffb /= 0.d0) then
           if(vparab /= 0.d0) then
             ffi = ffb * BB(ijn(s))/BB(i) * abs(vpara/vparab)
           else
             ffi = 0d0
           end if
 
           write(62, 74) lam(i), vpara, vperp, vparab, vperpb, kpara, kperp, kall, ffb, ffi
         endif
         
       endif
      enddo !p
    endif
   enddo !i
   
   close(62)
   
  enddo !s
  
  return
  
 end subroutine disfun2
 
 
 