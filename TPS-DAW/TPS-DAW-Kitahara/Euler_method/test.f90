  integer, parameter :: i=1000
  double precision :: time
  character(10) :: date1,date2,date3
  integer       :: date_time(10)
  double precision :: a0(0:i), a1(0:i)
  call date_and_time(date1,date2,date3,date_time)
  do j=0,i
     a1(i)=j*1.d-1  
     write(*,*) time!date_time(1),date_time(2),date_time(3),date_time(4),date_time(5),date_time(6),date_time(7)
  end do
   end program
   
