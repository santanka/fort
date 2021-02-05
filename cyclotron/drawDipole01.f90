implicit none

real :: Bx,Bz,Br,Bl,rr,lam,xx,zz,B,dd
real, parameter :: pi = 4.*atan(1.), factp = 180./pi, L = 6.371e6
real, parameter :: mu0 = 1.256e-6, Bm = 8.05e22, Be = mu0*Bm/4./pi/L**3
integer :: ii,jj
integer, parameter :: istep = 10000
character :: filename*30,chj*3,chd*6

print *, "Input magnetic latitude number of start point (degree)" 
read *, jj !‰Šú¥‹CˆÜ“x“ü—Í

write(chj,80) jj !jj‚ğ•¶š—ñ‚É•ÏŠ·
chj = adjustl(chj) 
!adjustl():string‚ğ¶‚ÉŠñ‚¹‚é

print *, "Input distance between spots (*Re)" 
read *, dd !delta“ü—Í

write(chd,70) dd !dd‚ğ•¶š—ñ‚É•ÏŠ·
chd = adjustl(chd) 
!adjustl():string‚ğ¶‚ÉŠñ‚¹‚é

filename = "dataDP01a("//chj//","//chd//").csv"
filename = adjustl(filename)

print *, filename
open(10,file=filename,status="unknown")

lam = dble(jj)/factp
rr = 1.
xx = rr*cos(lam)
zz = rr*sin(lam)
write(10,90) xx,zz

do ii = 1,istep
  Br = -2.*sin(lam)*Be/rr**3
  Bl = cos(lam)*Be/rr**3
  Bx = Br*cos(lam) - Bl*sin(lam)
  Bz = Br*sin(lam) + Bl*cos(lam)
  B = sqrt(Bx**2 + Bz**2)

  xx = xx + Bx/B*dd
  zz = zz + Bz/B*dd

  write(10,90) xx,zz

  rr = sqrt(xx**2 + zz**2)
  lam = asin(zz/rr)

  if(ii /= 1 .and. rr < 1.) exit

enddo

close(10)

70 format(ES6.0)
80 format(I3)
90 format(E12.4, ',', 1x, E12.4)

end