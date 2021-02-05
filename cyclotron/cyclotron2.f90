program cyclotron2 !Buneman-Boris's law

implicit none

real, parameter :: c = 2.99792458e8, we = 2.8e4, B0 = 1.e-6, q0 = 1.6021766208e-19
real, parameter :: m0 = 9.10938356e-31, pi = 4.*atan(1.),v0 = 5.9e5
real :: dt,x,y,vx0,vy0,vx1,vy1,t,s
character :: chdt*6,filename*19
integer :: ii

print *, "Input time interval"
read *, dt !ŽžŠÔŠÔŠu“ü—Í
write(chdt,70) dt
chdt = adjustl(chdt)

filename = "datacy2("//chdt//").csv"
filename = adjustl(filename)
print *, filename

open(10,file=filename,status="unknown")
open(20,file="datacy2velocity.csv",status="unknown")

print *, q0*B0/m0/we

x = 0.
y = 0.

vx0 = v0/c
vy0 = 0.

write(10,90) x,y
write(20,90) vx0,vy0

do ii = 1,nint(1./dt)
  x = x + vx0*dt
  y = y + vy0*dt

  write(10,90) x,y

  t = dt/2.*(q0*B0/m0)/we

  vx1 = vx0 + vy0*t
  vy1 = vy0 - vx0*t

  s = 2.*t/(1.+t**2)

  vx0 = vx0 + vy1*s
  vy0 = vy0 - vx1*s

  write(20,90) vx0,vy0

enddo

close(10)
close(20)

70 format(ES6.0)
80 format(I3)
90 format(ES18.9, ',', 1x, ES18.9)

end program cyclotron2
