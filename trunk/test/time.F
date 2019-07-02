      program test_cpu_time
      integer :: k
      real :: start, finish
      k = -3
      call cpu_time(start)
      print '("k=",i3," sign", i3)',k, isgn(k)
      call cpu_time(finish)
      print '("Time = ",f6.3," seconds.")',finish-start
      end program test_cpu_time
