      program test_jacobi
      use error
      use kinds
      implicit none
      real*8, external ::  jacobi_w
      real*8 :: a,b,x
      real*8, dimension(0:10,10,10) :: eval
      integer :: i, nord, nfile, w
      write(*,*) 'set a,b,x, nord'
      read (*,*) a,b,x, nord
      write(*,*) 'a,b,x, nord, jacobi',
     .  a,b,x, nord, jacobi_w(a,b,x,nord)

      stop 0
      
      do w=0,10
        do nord=1,10
          do i=1, 10
            x = -1.d0 + 2.d0 * i / 10.d0
            eval(w,nord, i) = jacobi_w(real(w), 1.d0, x, nord)
          enddo
        enddo
      enddo
      
      do w=0,10
        do nord=1,10
          write(*,*) 'weight, nord = ', w,nord
          write(*,7001) eval(w,nord,1:10)
        enddo
      enddo
 7001 format(10(E12.5,2x))

      end program test_jacobi
      
