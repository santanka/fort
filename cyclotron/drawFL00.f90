implicit none !暗黙の型宣言無効

real :: rr,lam,xx,zz,req !実数型変数
integer :: ii !整数型変数

real, parameter :: pi = 4.0*atan(1.0) !実数型定数の円周率
real, parameter :: factp = 180.0/pi !radとdegreeの変換用

!双極子磁場の磁力線の連続関数
open(21,file="dataFL00a.csv",status="unknown") !csvファイル作成,コード21

!初期位置
lam = -45.0/factp !磁気緯度
req = 1.0/cos(lam)**2 !磁気赤道面での長径距離
xx = cos(lam)
zz = sin(lam)
write(21,42) xx,zz

do ii = 1,1800 !繰り返し処理
  lam = lam + 1.0/factp
  rr = req*cos(lam)**2
  xx = rr*cos(lam)
  zz = rr*sin(lam)
  write(21,42) xx,zz

  if(rr <= 1.0) exit
enddo

close(21) !編集中のcsvファイル21を閉じる

!地球表面の描写
open(12,file="dataFL00aGS.csv",status="unknown") !csvファイル作成,コード12

do ii=0,360,2 !0-360を間隔2で繰り返す
  lam = dble(ii)/factp !dbleで倍精度実数型に変換
  xx = cos(lam)
  zz = sin(lam)
  write(12,42) xx,zz
enddo

close(12) !編集中のcsvファイル12を閉じる

42 format(E12.4, ',', 1x, E12.4) !コード42形式フォーマットの定義

end
