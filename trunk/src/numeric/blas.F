c-----------------------------------------------------------------------
c> @date Feb 2023
c-----------------------------------------------------------------------
      subroutine decomp(n,na,a,ip,Iflag)
      implicit none
c
      integer :: n,na
      real(8) :: a(na,n)
      integer :: ip(n)
      integer :: Iflag
c
      integer :: j,k,imax,inc
      real(8) :: c,pivmin
c
      real(8), parameter :: eps = 1.0d-15
c
      integer, external :: idamax
c
#if HP3D_DEBUG
      integer :: i
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,7001)
 7001   format('decomp: n,na = ',2i4)
        do i=1,n
          write(*,7002) i,a(i,1:n)
 7002     format('i = ',i4,2x,10e13.5,10(/,10x,10e13.5))
        enddo
        call pause
      endif
#endif
c
      pivmin = 1.d30
c
      Iflag=0
      inc=1
c
      do k=1,n
        imax = idamax(n-k+1,a(k,k),inc)
        ip(k) = k+imax-1
        if (abs(a(ip(k),k)).le.1.d-10) then
          write(*,*)'k,ip(k) = ',k,ip(k)
          write(*,*)'K,PIVOT=', k,abs(a(ip(k),k))
          call pause
        endif
        pivmin = dmin1(pivmin,abs(a(ip(k),k)) )
        if (dabs(a(ip(k),k)) > eps) go to 20
        Iflag=1
        return
c
 20     continue
        c = a(ip(k),k)
        a(ip(k),k) = a(k,k)
        a(k,k) = c
        c=1.d0/c
c
        call dscal(n-k,c,a(k+1,k),inc)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7003) k
 7003     format('decomp: k = ',i4)
          do i=1,n
            write(*,7002) i,a(i,1:n)
          enddo
          call pause
        endif
#endif
c
        do j=k+1,n
          c = a(ip(k),j)
          a(ip(k),j) = a(k,j)
          a(k,j) = c
          call daxpy(n-k ,-a(k,j),a(k+1,k),inc,a(k+1,j),inc)
        enddo
c
      enddo
c
      end subroutine decomp
c
c-----------------------------------------------------------------------
c> @date Feb 2023
c-----------------------------------------------------------------------
      subroutine gauss2(n,na,a,ip,b,x)
c
      implicit none
c
      integer :: n,na
      real(8) :: a(na,n),b(n),x(n)
      integer :: ip(n)
c
      integer :: inc,j,k
      real(8) :: c
c
#if HP3D_DEBUG
      integer :: i
      integer :: iprint
      iprint = 0
c
      if (iprint.eq.1) then
        write(*,7001)
 7001   format('gauss2: n,na = ',2i4)
        do i=1,n
          write(*,7002) i,a(i,1:n)
 7002     format('i = ',i4,2x,10e12.5)
        enddo
        write(*,7003) ip(1:n)
 7003   format('ip = ',30i4)
        write(*,7004)
 7004   format('b = ')
        write(*,7005) b(1:n)
 7005   format(10e12.5)
        call pause
      endif
#endif
c
      inc=1
c
      do k=1,n-1
        c=b(ip(k))
        b(ip(k))=b(k)
        b(k)=c
        call daxpy(n-k,-b(k),a(k+1,k),inc,b(k+1),inc)
      enddo

      do j=n,1,-1
        x(j)=b(j)/a(j,j)
        call daxpy(j-1,-x(j),a(1,j),inc,b,inc)
      enddo
c
      end subroutine gauss2
c
c-----------------------------------------------------------------------
c
c   routine name       - gausse
c
c-----------------------------------------------------------------------
c
c   purpose            - solution of a full matrix linear system
c                        of equations using gauss elimination
c                        without pivoting
c
c   arguments
!        in:      gk   - n by n matrix
!                 igk  - row dimension of matrix gk exactly as
!                        specified in the dimension statement in the
!                        calling program
!                 gf   - vector of length n (the right-hand side)
!                 n    - number of equations
!       out:      u    - vector of length n containing the solution
c
c   required routines  - rhsub,tri
c
c-----------------------------------------------------------------------
c> @date Feb 2023
c-----------------------------------------------------------------------
      subroutine gausse(gk,igk,gf,u,n)
c
      implicit none
c
      integer :: igk,n
      real(8) :: gk(igk,*),gf(*),u(*)
c
      call tri(gk,igk,n)
      call rhsub(gk,u,gf,igk,n)
c
      end subroutine gausse
c
c-----------------------------------------------------------------------
c> @date Feb 2023
c-----------------------------------------------------------------------
      subroutine matinv(n,na,a,ip,Work)
c
      implicit none
c
      integer :: n,na
      real(8) :: a(na,n)
      integer :: ip(n)
      real(8) :: Work(n)
c
      integer :: i,k
c
      do i=n,2,-1
        a(i,i)=1/a(i,i)
        call dscal (i-1, -a(i,i), a(1,i), 1)
        a(i-1,i)=a(i-1,i)/a(i-1,i-1)

        do k=i-1,2,-1
          call daxpy (k-1, -a(k,i), a(1,k), 1, a(1,i), 1)
          a(k-1,i)=a(k-1,i)/a(k-1,k-1)
        enddo
      enddo
c
      a(1,1)=1/a(1,1)
c
      do i=n-1,1,-1
        call dcopy(   i, a(1,i), 1, Work     , 1)
        call dcopy( n-i, 0.d0  , 0, Work(i+1), 1)
c
        do k=i+1,n
          call daxpy( n, -a(k,i), a(1,k), 1, Work, 1)
        enddo
c
        call dcopy( n, Work, 1, a(1,i), 1)
        call dswap( n, a(1,i), 1, a(1,ip(i)), 1)
      enddo
c
      end subroutine matinv
c
c-----------------------------------------------------------------------
c> @date Feb 2023
c-----------------------------------------------------------------------
      function sdot1(n,sx,Incx,sy,Incy)
c
      implicit none
c
      real(8) :: sdot1
      integer :: n
      real(8) :: sx(*),sy(*)
      integer :: Incx,Incy
c
      integer :: jx,jy,i
c
      sdot1=0.d0
      if(n.le.0) return
c
      jx=1
      if(Incx.lt.0) jx=1+(n-1)*(-Incx)
      jy=1
      if(Incy.lt.0) jy=1+(n-1)*(-Incy)
c
      do i=1,n
        sdot1 = sdot1+sx(jx)*sy(jy)
        jx=jx+Incx
        jy=jy+Incy
      enddo
c
      end function sdot1
c
c-----------------------------------------------------------------------
c> @date Feb 2023
c-----------------------------------------------------------------------
      function snrm21(n,sx,Incx)
c
      implicit none
c
      real(8) :: snrm21
      integer :: n
      real(8) :: sx(*)
      integer :: Incx
c
      real(8) :: sn
      integer :: i,jx
c
      sn=0.d0
      jx=1
      if(Incx.lt.0) jx=1+(n-1)*(-Incx)

      do i=1,n
        sn = sn + sx(jx)*sx(jx)
        jx=jx+Incx
      enddo

      snrm21=sqrt(sn)

      end function snrm21
