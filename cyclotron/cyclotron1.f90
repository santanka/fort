program cyclotron1 !leap-frog law

implicit none

real, parameter :: c = 2.99792458e8, we = 2.8e4, B0 = 1.e-6, q0 = 1.6021766208e-19
real, parameter :: m0 = 9.10938356e-31, pi = 4.*atan(1.),v0 = 5.9e5
real :: dt,x,y,vx0,vy0,vx1,vy1,vx2,vy2,v
character :: chdt*6,filename*19
integer :: ii

print *, "Input time interval"
read *, dt !ŽžŠÔŠÔŠu“ü—Í
write(chdt,70) dt
chdt = adjustl(chdt)

filename = "datacy1("//chdt//").csv"
filename = adjustl(filename)
print *, filename

open(10,file=filename,status="unknown")

v = v0/c
vx0 = v*sin(-2.*pi/dt/we)
vy0 = v*cos(-2.*pi/dt/we)

vx1 = 0.
vy1 = 1.

x = 0.
y = 0.

write(10,90) x,y

do ii = 1,6*nint(1/dt)
  vx2 = vx0 - vy1*dt
  vy2 = vy0 + vx1*dt
  x = x + vx2*dt
  y = y + vy2*dt

  write(10,90) x,y

  vx0 = vx1
  vy0 = vy1
  vx1 = vx2
  vy1 = vy2

  vx2 = vx0 - vy1*dt
  vy2 = vy0 + vx1*dt

  vx0 = vx1
  vy0 = vy1
  vx1 = vx2
  vy1 = vy2

enddo

close(10)

70 format(ES6.0)
80 format(I3)
90 format(ES18.9, ',', 1x, ES18.9)

end program cyclotron1
