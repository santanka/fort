program pds_E_kai !�n��

implicit none

!main�̕����Ƃ�
integer :: h, i, j, s, itn !do���p
integer, parameter :: N = 501 !grid��
integer, parameter :: Z = 500 !perp���xgrid��
integer, parameter :: MV = 251 !�ɏ��l�ƂȂ���W
integer, parameter :: CP = 0 !boundary�ȊO�̐Ód�|�e���V�������ω��_
double precision, parameter :: pi = 4.d0*atan(1.d0) !�~����
character(len=80) :: dummy !�ǂ܂Ȃ������p

!parameter�̓���
double precision, parameter :: ee = 1.602176634d-19 !�d�C�f��
double precision, parameter :: cc = 2.99792458d+08 !����
double precision, parameter :: mu0 = 1.25663706143592d-06 !�^��̓�����
double precision, parameter :: ep0 = 8.8541878128d-12 !�^��̗U�d��
double precision, parameter :: GG = 6.6743015d-11 !���L���͒萔
double precision, parameter :: cc_rate = 1.d0 !�����̊���

!�f���̃f�[�^
double precision, parameter :: Mdp = 7.8d+22 !�f���̎��C�o�Ɏq���[�����g
double precision, parameter :: Mp = 5.9722d+24 !�f���̎���
double precision, parameter :: Rp = 6.371d+06 !�f���̔��a
double precision, parameter :: omg = 7.292115024135739d-05 !�f���̎��]�p���g��
double precision, parameter :: Lp = 1.5d1 !�Ώۂ̎��͐���L�l
double precision, parameter :: lamN = acos(sqrt((1.d0+5.d5/Rp)/Lp)) !boundary:N��MLT
double precision, parameter :: lamS = -acos(sqrt((1.d0+5.d5/Rp)/Lp)) !boundary:S��MLT

!�q���̃f�[�^(�ݒ肵�Ȃ��ꍇ�͒l��0��)
double precision, parameter :: Ms = 0 !�q���̎���
double precision, parameter :: Ls = 0 !�q���O����L�l

!Boundary Condition(���q�̃f�[�^)
character(len=128) :: fileBC = 'pds_BC_E_kai_L=15_4.csv'
integer, parameter :: kind = 12 !���q�퐔
integer, parameter :: nsc = 4 !boundary:N�ł̐����x���������闱�q
integer, parameter :: ssc = 9 !boundary:S�ł̐����x���������闱�q
double precision, dimension(kind) :: sig !�����x
double precision, dimension(kind) :: Tperp !boundary�ł�perp���x[eV]
double precision, dimension(kind) :: Tpara !boundary�ł�para���x[eV]
double precision, dimension(kind) :: chr !�d��[/ee]
double precision, dimension(kind) :: mass !����
integer, dimension(kind) :: ijn !injection number

!Initial Condition(�����Ód�|�e���V�������z)
character(len=128) :: fileIC = 'pds_IC_E_kai_L=15_4.csv'
double precision, dimension(N) :: Phi !�Ód�|�e���V�������z

!��̐ݒ�
double precision, dimension(N) :: lam !MLT
double precision, dimension(N) :: rr !�f�����S�Ƃ̋���
double precision, dimension(N) :: rs !�q�����S�Ƃ̋���
double precision, dimension(N) :: ss !���͐���̍��W
double precision, dimension(N-1) :: ds !ss���W�Ԋu
double precision, dimension(N) :: BB !�������x

!�f�M�s�ϗ�
double precision, dimension(kind, Z) :: mu !��1�f�M�s�ϗ� 

!iteration
double precision, dimension(3, N) :: Phih !�������Ód�|�e���V����
double precision, dimension(3, kind) :: sigg !�����������x
double precision, dimension(3, N, kind) :: UU !�ʒu�G�l���M�[
double precision, dimension(3, N, kind, Z) :: amin !amin
double precision, dimension(3, N, kind, Z) :: amax !amax
double precision, dimension(3, N, kind) :: nd !�����x
double precision, dimension(3, N) :: rhod !�d�ז��x
double precision, dimension(3, N) :: rhodp !�d�ז��x(+)
double precision, dimension(3, N) :: rhodm !�d�ז��x(-)
double precision, dimension(3, N) :: rhop !�d�ז��x(Poisson eq.)
double precision, dimension(3, N) :: cvg !�����l
double precision :: cvn !�����l
double precision :: DDD !�Œ�����l
double precision, dimension(N) :: nPhi !�Ód�|�e���V�����X�V�l
double precision, dimension(kind) :: nsig !�����x�X�V�l


!�ۑ��t�@�C����
integer, parameter :: channel = 1 !�͂��߂���(1)or�Â�����(2)
character(len=128) :: fileresult = 'pds_E_kai_L=15_4_result.csv'
72 format(1PE25.15E3, 18(',', 1PE25.15E3)) !kind+7
character(len=128) :: filepote = 'pds_E_kai_L=15_4_potential.csv'
82 format(1PE25.15E3, 15(',', 1PE25.15E3)) !kind+4
character(len=128) :: filemin = 'pds_E_kai_L=15_4_min.csv'
92 format(1PE25.15E3) !1
character(len=128) :: filecheck = 'pds_E_kai_L=15_4_check.csv'
52 format(1PE25.15E3, 4(',', 1PE25.15E3)) !5(double precision)
42 format(I5)
32 format(1PE25.15E3, 7(',', 1PE25.15E3)) !8(double precision)


!/////�f�[�^���o/////
!BC�̒��o
open(30, file=fileBC, action="read", status="old")
read(30, *) !1�s�ǂݔ�΂�
do s = 1, kind
 read(30, *) dummy, sig(s), Tperp(s), Tpara(s), chr(s), mass(s), ijn(s)
enddo !s
close(30)
chr = chr*ee !������[C]��
Tperp = Tperp*ee !������[J]��
Tpara = Tpara*ee !������[J]��

if(channel == 1) then
  !IC�̒��o
  open(40, file=fileIC, action="read", status="old")
  do i = 1, N
   read(40, *) Phi(i)
  enddo !i
  close(40)
  DDD = 0.d0
  
 else if(channel == 2) then
  !�Œ�l�X�V
  open(40, file=filemin)
  read(40, *) DDD
  do i = 1, N
   read(40, *) lam(i), ss(i), BB(i), Phi(i), nd(1, i, :), rhod(1, i), rhop(1, i), cvg(1, i)
  enddo !i
  do s = 1, kind
   read(40, *) sig(s)
  enddo !s
endif

!/////��̐ݒ�/////

!MLT�̕���
do i = 1, N
 lam(i) = lamS + (lamN-lamS)*dble(i-1)/dble(N-1)
enddo !i

!�f�����S�Ƃ̋���
rr = Lp*Rp*cos(lam)**2.d0

!�q�����S�Ƃ̋���
if(Ms /= 0) then
  rs = sqrt(rr**2.d0+(Ls*Rp)**2.d0-2.d0*rr*Ls*Rp*cos(lam))
endif

!���͐���̍��W
ss = Lp*Rp*(sin(lam)*sqrt(1.d0+3.d0*sin(lam)**2.d0)/2.d0 + asinh(sqrt(3.d0)*sin(lam))/2.d0/sqrt(3.d0))

!ss���W�Ԋu
do i = 1, N-1
 ds(i) = ss(i+1) - ss(i)
enddo !i

!�������x
BB = mu0*Mdp/4.d0/pi/(Lp*Rp)**3.d0 * sqrt(1.d0+3.d0*sin(lam)**2.d0) / cos(lam)**6.d0


!�f�M�s�ϗ�
call AI(N, Z, kind, cc_rate, cc, mass, BB, ijn, mu)


!/////�ȉ�iteration/////
itn = 0

do !itn
 itn = itn + 1
 
 !�Ód�|�e���V��������
 call pepm(N, Phi, Phih)
 
 !�����x����
 call sigx(kind, nsc, ssc, sig, sigg)
 
 !�ʒu�G�l���M�[
 call EP(N, kind, GG, Mp, mass, rr, omg, lam, Ms, rs, chr, Phih, UU)
 
 !Accessibility
 call access(N, Z, kind, cc, cc_rate, mass, UU, mu, BB, ijn, amin, amax)
 
 !�ϕ�
 call dense(N, Z, kind, Tperp, Tpara, sigg, mu, BB, ijn, amin, amax, nd)
 
 !�d�ז��x(���z�֐�)
 call chde(N, kind, chr, nd, rhod, rhodp, rhodm)
 
 !�d�ז��x(Poisson������)
 call chpo(N, ep0, ds, Phih, rhop)
 
 !�����l
 call CV(N, rhod, rhodp, rhodm, rhop, cvg, cvn)
 print *, itn, cvn
 
 !�t�@�C��������
 open(50, file = fileresult)
 do i = 1, N
  write(50, 72) lam(i), ss(i), BB(i), Phi(i), nd(1, i, :), rhod(1, i), rhop(1, i), cvg(1, i)
 enddo !i
 close(50)
 
 open(60, file = filepote)
 do i = 1, N
  write(60, 82) lam(i), ss(i), BB(i), Phi(i), UU(1, i, :)
 enddo !i
 close(60)
 
 !�Œ�����l
 if(DDD >= cvn .or. DDD == 0.d0) then
   DDD = cvn
   print *, "cvn min", DDD
   open(33, file = filemin)
   write(33, 92) cvn
   do i = 1, N
    write(33, 72) lam(i), ss(i), BB(i), Phi(i), nd(1, i, :), rhod(1, i), rhop(1, i), cvg(1, i)
   enddo !i
   do s = 1, kind
    write(33, 52) sig(s), Tperp(s)/ee, Tpara(s)/ee, chr(s)/ee, mass(s)
   enddo !s
   do s = 1, kind
    write(33, 42) ijn(s)
   enddo !s
   close(33)
 endif
 
 !�����`�F�b�N
 !if(cvn < 1.d-5) then
 !  print *, "finish"
 !  stop
 !endif
 
 !Newton�@(�Ód�|�e���V����)
 call NewtonPhi(N, MV, CP, Phih, cvg, nPhi)
 
 !Newton�@(�����x)
 call Newtonsig(N, kind, ijn, sigg, cvg, nsig)
 
 open(70, file = filecheck)
 write(70, 32) sigg(:, ssc), nsig(ssc), nsig(ssc)-sigg(1, ssc), cvg(:, 1)
 do i = 2, N-1
  write(70, 32) Phih(:, i), nPhi(i), nPhi(i)-Phih(1, i), cvg(:, i)
 enddo !i
 write(70, 32) sigg(:, nsc), nsig(nsc), nsig(nsc)-sigg(1, nsc), cvg(:, N)
 close(70)
 
 !NaN, infinity�`�F�b�N
 do i = 1, N
  if(nPhi(i) /= nPhi(i) .or. nPhi(i) == nPhi(i)-1.d0) then
    if(nPhi(i) /= nPhi(i)) then
      print *, "NaN����"
     else if(nPhi(i) == nPhi(i)-1.d0) then
      print *, "infinity����"
    endif
    stop
  endif
  if(cvg(1, i) /= cvg(1, i) .or. cvg(1, i) == cvg(1, i)-1.d0) then
    if(cvg(1, i) /= cvg(1, i)) then
      print *, "NaN����"
     else if(cvg(1, i) == cvg(1, i)-1.d0) then
      print *, "infinity����"
    endif
    stop
  endif
 enddo !i
 
 !�X�V
 Phi = nPhi
 sig = nsig
 
 
enddo !itn


end program pds_E_kai





!subroutine & function

subroutine AI(N, Z, kind, cc_rate, cc, mass, BB, ijn, mu)
  implicit none
  integer, intent(in) :: N, Z, kind
  double precision, intent(in) :: cc_rate, cc
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
 

!�Ód�|�e���V��������
subroutine pepm(N, Phi, Phih)
 implicit none
 integer, intent(in) :: N
 double precision, dimension(N), intent(in) :: Phi
 double precision, dimension(3, N), intent(out) :: Phih
 
 Phih(1, :) = Phi
 Phih(2, :) = Phi + 1.d-7
 Phih(3, :) = Phi - 1.d-7
 
 return
 
end subroutine pepm


!�����x����
subroutine sigx(kind, nsc, ssc, sig, sigg)
 implicit none
 integer, intent(in) :: kind, nsc, ssc
 double precision, dimension(kind), intent(in) :: sig
 double precision, dimension(3, kind), intent(out) :: sigg
 
 sigg(1, :) = sig
 sigg(2, :) = sig
 sigg(3, :) = sig
 sigg(2, nsc) = sig(nsc)*(1.d0+1.d-8)
 sigg(3, nsc) = sig(nsc)*(1.d0-1.d-8)
 sigg(2, ssc) = sig(ssc)*(1.d0+1.d-8)
 sigg(3, ssc) = sig(ssc)*(1.d0-1.d-8)
 
 return
 
end subroutine sigx


!�ʒu�G�l���M�[
subroutine EP(N, kind, GG, Mp, mass, rr, omg, lam, Ms, rs, chr, Phih, UU)
 implicit none
 integer, intent(in) :: N, kind
 double precision, intent(in) :: GG, Mp, omg, Ms
 double precision, dimension(kind), intent(in) :: mass, chr
 double precision, dimension(N), intent(in) :: rr, lam, rs
 double precision, dimension(3, N), intent(in) :: Phih
 double precision, dimension(3, N, kind), intent(out) :: UU
 
 integer :: h, i
 
 do h = 1, 3
  do i = 1, N
   UU(h, i, :) = -GG*Mp*mass/rr(i) !�f���d��
   UU(h, i, :) = UU(h, i, :) - mass*(omg*rr(i)*cos(lam(i)))**2.d0/2.d0 !�f�����S��
   if(Ms /= 0) UU(h, i, :) = UU(h, i, :) - GG*Ms*mass/rs(i) !�q���d��
   UU(h, i, :) = UU(h, i, :) + chr*Phih(h, i) !�N�[������
  enddo !i
 enddo !h
 
 return
 
end subroutine EP


!Accessibility
subroutine access(N, Z, kind, cc, cc_rate, mass, UU, mu, BB, ijn, amin, amax)
  implicit none
  integer, intent(in) :: N, Z, kind
  double precision, intent(in) :: cc, cc_rate
  double precision, dimension(kind), intent(in) :: mass
  double precision, dimension(3, N, kind), intent(in) :: UU
  double precision, dimension(kind, Z), intent(in) :: mu
  double precision, dimension(N), intent(in) :: BB
  integer, dimension(kind), intent(in) :: ijn
  double precision, dimension(3, N, kind, Z), intent(out) :: amin, amax
  
  integer :: h, i, s, j, k, t
  double precision, dimension(3, N, kind, Z) :: EE !UU+mu*BB
  double precision :: WW
  
  
  !�G�l���M�[�̘a(UU+mu*BB)
  do h = 1, 3
   do i = 1, N
    do s = 1, kind
     EE(h, i, s, :) = UU(h, i, s) + mu(s, :)*BB(i)
    enddo !s
   enddo !i
  enddo !h
  
  !amin
  do s = 1, kind
   do i = 1, N
    if(i < ijn(s)) then
      do h = 1, 3
       do j = 1, Z
        t = i+1
        do k = i+1, ijn(s)
         if(EE(1, k, s, j) > EE(1, t, s, j)) t = k
        enddo !k
        if(EE(h, i, s, j) > EE(1, t, s, j)) then
          amin(h, i, s, j) = sqrt(EE(h, i, s, j) - EE(1, ijn(s), s, j))
         else if(EE(h, i, s, j) <= EE(1, t, s, j)) then
          amin(h, i, s, j) = sqrt(EE(1, t, s, j) - EE(1, ijn(s), s, j))
        endif
       enddo !j
      enddo !h
      
     else if(i == ijn(s)) then
      amin(:, ijn(s), s, :) = 0.d0
      
     else if(i > ijn(s)) then
      do h = 1, 3
       do j = 1, Z
        t = ijn(s)
        do k = ijn(s), i-1
         if(EE(1, k, s, j) > EE(1, t, s, j)) t = k
        enddo !k
        if(EE(h, i, s, j) > EE(1, t, s, j)) then
          amin(h, i, s, j) = sqrt(EE(h, i, s, j) - EE(1, ijn(s), s, j))
         else if(EE(h, i, s, j) <= EE(1, t, s, j)) then
          amin(h, i, s, j) = sqrt(EE(1, t, s, j) - EE(1, ijn(s), s, j))
        endif
       enddo !j
      enddo !h
    endif
   enddo !i
  enddo !s
  
  !amax
  do s = 1, kind
   do i = 1, N
    if((i <= ijn(s) .and. ijn(s) == N) .or. (i < ijn(s) .and. ijn(s) /= N)) then
      if(i == 1) then
        amax(:, 1, s, :) = 0.d0
       else if(i /= 1) then
        do j = 1, Z
         t = 1
         do k = 1, i-1
          if(EE(1, k, s, j) > EE(1, t, s, j)) t = k
         enddo !k
         if(i == ijn(s) .and. ijn(s) == N) then
           do h = 1, 3
            WW = EE(1, t, s, j) - EE(h, ijn(s), s, j)
            if(WW <= 0.d0) then
              amax(h, i, s, j) = 0.d0
             else if(WW > 0.d0 .and. WW < mass(s)/2.d0*(cc_rate*cc)**2.d0) then
              amax(h, i, s, j) = sqrt(WW)
             else if(WW >= mass(s)/2.d0*(cc_rate*cc)**2.d0) then
              amax(h, i, s, j) = sqrt(mass(s)/2.d0)*cc_rate*cc
            endif
           enddo !h
          else
           WW = EE(1, t, s, j) - EE(1, ijn(s), s, j)
           if(WW <= 0.d0) then
             amax(:, i, s, j) = 0.d0
            else if(WW > 0.d0 .and. WW < mass(s)/2.d0*(cc_rate*cc)**2.d0) then
             amax(:, i, s, j) = sqrt(WW)
            else if(WW >= mass(s)/2.d0*(cc_rate*cc)**2.d0) then
             amax(:, i, s, j) = sqrt(mass(s)/2.d0)*cc_rate*cc
           endif
         endif
        enddo !j
      endif
      
     else if((i >= ijn(s) .and. ijn(s) == 1) .or. (i > ijn(s) .and. ijn(s) /= 1)) then
      if(i == N) then
        amax(:, N, s, :) = 0.d0
       else if(i /= N) then
        do j = 1, Z
         t = i+1
         do k = i+1, N
          if(EE(1, k, s, j) > EE(1, t, s, j)) t = k
         enddo !k
         if(i == ijn(s) .and. ijn(s) == 1) then
           do h = 1, 3
            WW = EE(1, t, s, j) - EE(h, ijn(s), s, j)
            if(WW <= 0.d0) then
              amax(h, i, s, j) = 0.d0
             else if(WW > 0.d0 .and. WW < mass(s)/2.d0*(cc_rate*cc)**2.d0) then
              amax(h, i, s, j) = sqrt(WW)
             else if(WW >= mass(s)/2.d0*(cc_rate*cc)**2.d0) then
              amax(h, i, s, j) = sqrt(mass(s)/2.d0)*cc_rate*cc
            endif
           enddo !h
          else
           WW = EE(1, t, s, j) - EE(1, ijn(s), s, j)
           if(WW <= 0.d0) then
             amax(:, i, s, j) = 0.d0
            else if(WW > 0.d0 .and. WW < mass(s)/2.d0*(cc_rate*cc)**2.d0) then
             amax(:, i, s, j) = sqrt(WW)
            else if(WW >= mass(s)/2.d0*(cc_rate*cc)**2.d0) then
             amax(:, i, s, j) = sqrt(mass(s)/2.d0)*cc_rate*cc
           endif
         endif
        enddo !j
      endif
      
     else if(i == ijn(s) .and. ijn(s) /= 1 .and. ijn(s) /= N) then !��O����
      amax(:, i, s, :) = 1.d100
      
    endif
   enddo !i
  enddo !s
  
  return
  
 end subroutine access


subroutine dense(N, Z, kind, Tperp, Tpara, sigg, mu, BB, ijn, amin, amax, nd)
 implicit none
 integer, intent(in) :: N, Z, kind
 double precision, dimension(kind), intent(in) :: Tperp, Tpara
 double precision, dimension(3, kind), intent(in) :: sigg
 double precision, dimension(kind, Z), intent(in) :: mu
 double precision, dimension(N), intent(in) :: BB
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(3, N, kind, Z), intent(in) :: amin, amax
 double precision, dimension(3, N, kind), intent(out) :: nd
 
 integer :: h, i, s, j
 double precision, dimension(3, N, kind, Z) :: TL, TM, TT
 
 do s = 1, kind
  TL(:, :, s, :) = 1.d0 - erf(amin(:, :, s, :)/sqrt(Tpara(s)))
 enddo !s
 
 do h = 1, 3
  do i = 1, N
   do s = 1, kind
    do j = 1, Z
     if(amax(h, i, s, j) > amin(h, i, s, j)) then
       TM(h, i, s, j) = erf(amax(h, i, s, j)/sqrt(Tpara(s))) - erf(amin(h, i, s, j)/sqrt(Tpara(s)))
      else if(amax(h, i, s, j) <= amin(h, i, s, j)) then
       TM(h, i, s, j) = 0.d0
     endif
    enddo !j
   enddo !s
  enddo !i
 enddo !h
 
 
 do s = 1, kind
  do j = 1, Z
   TT(:, :, s, j) = (TL(:, :, s, j) + TM(:, :, s, j)) * exp(-mu(s, j)*BB(ijn(s))/Tperp(s))
  enddo !j
 enddo !s
 
 nd = 0.d0
 
 do h = 1, 3
  do i = 1, N
   if(i /= 1 .and. i /= N) then
     do j = 1, Z-1
      nd(h, i, :) = nd(h, i, :) + (TT(h, i, :, j)+TT(h, i, :, j+1))/2.d0*(mu(:, j+1)-mu(:, j))
     enddo !j
     nd(h, i, :) = sigg(1, :)*BB(ijn)/2.d0/Tperp*nd(h, i, :)
    else if(i == 1 .or. i == N) then
     do j = 1, Z-1
      nd(h, i, :) = nd(h, i, :) + (TT(1, i, :, j)+TT(1, i, :, j+1))/2.d0*(mu(:, j+1)-mu(:, j))
     enddo !j
     nd(h, i, :) = sigg(h, :)*BB(ijn)/2.d0/Tperp*nd(h, i, :)
   endif
  enddo !i
 enddo !h
 
 return
 
end subroutine dense


subroutine chde(N, kind, chr, nd, rhod, rhodp, rhodm)
 implicit none
 integer, intent(in) :: N, kind
 double precision, dimension(kind), intent(in) :: chr
 double precision, dimension(3, N, kind), intent(in) :: nd
 double precision, dimension(3, N), intent(out) :: rhod, rhodp, rhodm
 
 integer :: s
 
 rhod = 0.d0
 rhodp = 0.d0
 rhodm = 0.d0
 
 do s = 1, kind
  rhod = rhod + chr(s)*nd(:, :, s)
  if(chr(s) > 0.d0) then
    rhodp = rhodp + chr(s)*nd(:, :, s)
   else if(chr(s) < 0.d0) then
    rhodm = rhodm + chr(s)*nd(:, :, s)
  endif
 enddo !s
 
 return
 
end subroutine chde


subroutine chpo(N, ep0, ds, Phih, rhop)
 implicit none
 integer, intent(in) :: N
 double precision, intent(in) :: ep0
 double precision, dimension(N-1), intent(in) :: ds
 double precision, dimension(3, N), intent(in) :: Phih
 double precision, dimension(3, N), intent(out) :: rhop
 
 integer :: i
 
 do i = 1, N
  if(i == 1 .or. i == N) then
    rhop(:, i) = 0.d0
   else if(i /= 1 .and. i /= N) then
    rhop(:, i) = Phih(:, i-1)/ds(i-1)/(ds(i-1)+ds(i))
    rhop(:, i) = rhop(:, i) + Phih(:, i+1)/ds(i)/(ds(i-1)+ds(i))
    rhop(:, i) = rhop(:, i) - Phih(:, i)/ds(i-1)/ds(i)
  endif
 enddo !i
 
 rhop = -2.d0*ep0*rhop
 
 return
 
end subroutine chpo


subroutine CV(N, rhod, rhodp, rhodm, rhop, cvg, cvn)
 implicit none
 integer, intent(in) :: N
 double precision, dimension(3, N), intent(in) :: rhod, rhodp, rhodm, rhop
 double precision, dimension(3, N), intent(out) :: cvg
 double precision, intent(out) :: cvn
 
 integer :: i
 
 cvg = sqrt((rhod-rhop)**2.d0/rhodp/(-rhodm))
 
 cvn = 0.d0
 
 do i = 1, N
  cvn = cvn + cvg(1, i)**2.d0
 enddo !i
 
 cvn = sqrt(cvn/dble(N))
 
 return
 
end subroutine CV


!Newton�@(�Ód�|�e���V����)
subroutine NewtonPhi(N, MV, CP, Phih, cvg, nPhi)
 implicit none
 integer, intent(in) :: N, MV, CP
 double precision, dimension(3, N), intent(in) :: Phih, cvg
 double precision, dimension(N), intent(out) :: nPhi
 
 integer :: i, k, search
 double precision :: CC, ser
 
 
 nPhi = Phih(1, :)
 
 do i = 2, N-1
  CC = 0.d0
  if((cvg(2, i) < cvg(1, i) .or. cvg(3, i) < cvg(1, i)) .and. i /= CP) then
    if(cvg(2, i) == cvg(1, i)) then
      CC = (Phih(3, i)-Phih(1, i))/(cvg(3, i)-cvg(1, i))*cvg(1, i)
     else if(cvg(3, i) == cvg(1, i)) then
      CC = (Phih(2, i)-Phih(1, i))/(cvg(2, i)-cvg(1, i))*cvg(1, i)
     else
      CC = ((Phih(3, i)-Phih(1, i))/(cvg(3, i)-cvg(1, i)) + (Phih(2, i)-Phih(1, i))/(cvg(2, i)-cvg(1, i)))/2.d0*cvg(1, i)
    endif
  endif
  if((abs(CC) < cvg(1, i) .and. cvg(1, i) <= 1.d-2) .and. i /= CP) then
    nPhi(i) = Phih(1, i) - CC
   else if(abs(CC) >= cvg(1, i) .and. cvg(1, i) <= 1.d-2) then
    nPhi(i) = Phih(1, i) - CC/abs(CC)*cvg(1, i)
   else if(abs(CC) < 1.d-2 .and. cvg(1, i) > 1.d-2) then
    nPhi(i) = Phih(1, i) - CC
   else if(abs(CC) >= 1.d-2 .and. cvg(1, i) > 1.d-2) then
    nPhi(i) = Phih(1, i) - CC/abs(CC)*1.d-2
  endif
 enddo !i
 
 search = 1
 if(MV > 1 .and. MV <= N) then
   do while(search == 1)
    search = 0
    do k = 1, 5
     do i = 1, MV-1
      if(nPhi(i) < nPhi(i+1) .and. (i /= CP .and. i+1 /= CP)) then
        if(search == 0) search = 1
        if(i == 1) then
          nPhi(i+1) = 2.d0*nPhi(i)-nPhi(i+1)
         else
          ser = nPhi(i)
          nPhi(i) = nPhi(i+1)
          nPhi(i+1) = ser
        endif
      endif
     enddo !i
    enddo
   enddo
 endif
 
 search = 1
 if(MV < N .and. MV >= 1) then
   do while(search == 1)
    search = 0
    do k = 1, 5
     do i = MV, N-1
      if(nPhi(i+1) < nPhi(i) .and. (i /= CP .and. i+1 /= CP)) then
        if(search == 0) search = 1
        if(i == N-1) then
          nPhi(i) = 2.d0*nPhi(i+1)-nPhi(i)
         else
          ser = nPhi(i)
          nPhi(i) = nPhi(i+1)
          nPhi(i+1) = ser
        endif
      endif
     enddo !i
    enddo !k
   enddo
 endif
 
 return
 
end subroutine NewtonPhi


!Newton�@(�����x)
subroutine Newtonsig(N, kind, ijn, sigg, cvg, nsig)
 implicit none
 integer, intent(in) :: N, kind
 integer, dimension(kind), intent(in) :: ijn
 double precision, dimension(3, kind), intent(in) :: sigg
 double precision, dimension(3, N), intent(in) :: cvg
 double precision, dimension(kind), intent(out) :: nsig
 
 integer :: s
 double precision :: CC
 
 do s = 1, kind
  if(sigg(1, s) /= sigg(2, s)) then
    CC = 0.d0
    if(cvg(2, ijn(s)) < cvg(1, ijn(s)) .or. cvg(3, ijn(s)) < cvg(1, ijn(s))) then
      if(cvg(2, ijn(s)) == cvg(1, ijn(s))) then
        CC = (sigg(3, s)-sigg(1, s))/(cvg(3, ijn(s))-cvg(1, ijn(s)))*cvg(1, ijn(s))
       else if(cvg(3, ijn(s)) == cvg(1, ijn(s))) then
        CC = (sigg(2, s)-sigg(1, s))/(cvg(2, ijn(s))-cvg(1, ijn(s)))*cvg(1, ijn(s))
       else
        CC = ((sigg(3, s)-sigg(1, s))/(cvg(3, ijn(s))-cvg(1, ijn(s)))+(sigg(2, s)-sigg(1, s))/(cvg(2, ijn(s))-cvg(1, ijn(s))))&
        &/2.d0*cvg(1, ijn(s))
      endif
    endif
    
    if(abs(CC) < sigg(1, s)*1.d-1) then
      nsig(s) = sigg(1, s) - CC
     else if(abs(CC) >= sigg(1, s)*1.d-1) then
      nsig(s) = sigg(1, s) - CC/abs(CC)*sigg(1, s)*1.d-1
    endif
    
   else if(sigg(1, s) == sigg(2, s)) then
    nsig(s) = sigg(1, s)
    
  endif
 enddo !s
 
 return
 
end subroutine Newtonsig

