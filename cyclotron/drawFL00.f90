implicit none !�Öق̌^�錾����

real :: rr,lam,xx,zz,req !�����^�ϐ�
integer :: ii !�����^�ϐ�

real, parameter :: pi = 4.0*atan(1.0) !�����^�萔�̉~����
real, parameter :: factp = 180.0/pi !rad��degree�̕ϊ��p

!�o�Ɏq����̎��͐��̘A���֐�
open(21,file="dataFL00a.csv",status="unknown") !csv�t�@�C���쐬,�R�[�h21

!�����ʒu
lam = -45.0/factp !���C�ܓx
req = 1.0/cos(lam)**2 !���C�ԓ��ʂł̒��a����
xx = cos(lam)
zz = sin(lam)
write(21,42) xx,zz

do ii = 1,1800 !�J��Ԃ�����
  lam = lam + 1.0/factp
  rr = req*cos(lam)**2
  xx = rr*cos(lam)
  zz = rr*sin(lam)
  write(21,42) xx,zz

  if(rr <= 1.0) exit
enddo

close(21) !�ҏW����csv�t�@�C��21�����

!�n���\�ʂ̕`��
open(12,file="dataFL00aGS.csv",status="unknown") !csv�t�@�C���쐬,�R�[�h12

do ii=0,360,2 !0-360���Ԋu2�ŌJ��Ԃ�
  lam = dble(ii)/factp !dble�Ŕ{���x�����^�ɕϊ�
  xx = cos(lam)
  zz = sin(lam)
  write(12,42) xx,zz
enddo

close(12) !�ҏW����csv�t�@�C��12�����

42 format(E12.4, ',', 1x, E12.4) !�R�[�h42�`���t�H�[�}�b�g�̒�`

end
