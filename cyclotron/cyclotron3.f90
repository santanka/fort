program cyclotron3

implicit none

real, parameter :: c = 2.99792458e8, we = 2.8e4, q0 = 1.6021766208e-19
real, parameter :: m0 = 9.10938356e-31, pi = 4.*atan(1.) !,v0 = 5.9e5
real :: dt,x,y,z,vx0,vy0,vz0,vx1,vy1,vz1,r,lam,Bx,By,Bz,v0
real :: tx,ty,tz,sx,sy,sz,t,W
real, parameter :: mu0 = pi*4.e-7, Me = 8.05e22, RE = 6378.1e3, factp = 180./pi
character :: chdt*6,filename*19
integer :: ii,imx

print *, "Input time interval"
read *, dt !ŽžŠÔŠÔŠu“ü—Í
write(chdt,70) dt
chdt = adjustl(chdt)

print *, "Input iteration number"
read *, imx

print *, "Input Energy [eV]"
read *, W

print *, "Input equator pitch angle [deg]"
read *, lam

filename = "datacy3("//chdt//").csv"
filename = adjustl(filename)
print *, filename

open(10,file=filename,status="unknown")
open(20,file="datacyvel3.csv",status="unknown")

v0 = sqrt(2.*(W*q0)/m0)

x = 3.*RE
y = 0.
z = 0.
r = x
vx0 = v0*sin(lam/factp)*cos(45./factp)
vy0 = v0*sin(lam/factp)*sin(45./factp)
vz0 = v0*cos(lam/factp)

write(10,90) x/RE,y/RE,z/RE

do ii = 1, imx
  x = x + vx0*dt/we
  y = y + vy0*dt/we
  z = z + vz0*dt/we

  write(10,90) x/RE,y/RE,z/RE
  write(20,90) vx0,vy0,vz0

  r = sqrt(x**2+y**2+z**2)
  Bx = -3./4.*mu0*Me/pi/(r**5)*x*z
  By = -3./4.*mu0*Me/pi/(r**5)*y*z
  Bz = mu0*Me/4./pi/(r**3)*(1.-3.*z**2/(r**2))

  tx = q0*dt/we/2./m0*Bx
  ty = q0*dt/we/2./m0*By
  tz = q0*dt/we/2./m0*Bz
  t = sqrt(tx**2+ty**2+tz**2)
  sx = 2*tx/(1.+t**2)
  sy = 2*ty/(1.+t**2)
  sz = 2*tz/(1.+t**2)

  vx1 = vx0 + vy0*tz - vz0*ty
  vy1 = vy0 + vz0*tx - vx0*tz
  vz1 = vz0 + vx0*ty - vy0*tx

  vx0 = vx0 + vy1*sz - vz1*sy
  vy0 = vy0 + vz1*sx - vx1*sz
  vz0 = vz0 + vx1*sy - vy1*sx

  if(r < RE) exit

enddo

close(10)
close(20)

print *, "finish"

70 format(ES6.0)
90 format(ES18.9, 2(',', 1x,ES18.9))

end program cyclotron3
