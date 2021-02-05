program dipoleRK

implicit none

real :: Bx,Bz,Br,Bl,lam,xx,zz,dd,rr
real :: Bx0,Bz0,Bx1,Bz1
real :: kx1,kz1
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
chd = adjustl(chd) !adjustl():string‚ğ¶‚ÉŠñ‚¹‚é

filename = "dataDPrk2a("//chj//","//chd//").csv"
filename = adjustl(filename)
print *, filename

open(10,file=filename,status="unknown")

lam = dble(jj)/factp
rr = 1.
xx = rr*cos(lam)
zz = rr*sin(lam)
write(10,90) xx,zz
do ii = 1,istep
  call calcDP(xx,zz,Bx0,Bz0)

  kx1 = dd*Bx0/2.
  kz1 = dd*Bz0/2.

  call calcDP(xx+kx1,zz+kz1,Bx1,Bz1)

  xx = xx + dd*Bx1
  zz = zz + dd*Bz1

  write(10,90) xx,zz

  rr = sqrt(xx**2 + zz**2)
  print *,rr

  if(ii /= 1 .and. rr <= 1.) exit

enddo

close(10)

70 format(ES6.0)
80 format(I3)
90 format(E12.4, ',', 1x, E12.4)

end program dipoleRK

subroutine calcDP(xx,zz,Bx,Bz)
  real :: rr,lam
  real, parameter :: pi = 4.*atan(1.), factp = 180./pi, L = 6.371e6
  real, parameter :: mu0 = 1.256e-6, Bm = 8.05e22, Be = mu0*Bm/4./pi/L**3

  rr = sqrt(xx**2 + zz**2)
  lam = asin(zz/rr)

  Br = -2.*sin(lam)*Be/rr**3
  Bl = cos(lam)*Be/rr**3
  Bx = (Br*cos(lam) - Bl*sin(lam))/sqrt(Br**2 + Bl**2)
  Bz = (Br*sin(lam) + Bl*cos(lam))/sqrt(Br**2 + Bl**2)

  return

  end subroutine calcDP


