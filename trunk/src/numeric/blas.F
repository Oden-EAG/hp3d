
      subroutine decomp(n,na,a,ip,Iflag)
#include "syscom.blk"

      dimension a(na,n),ip(n)
c
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
c
      pivmin = 1.d30
c
      Iflag=0
      inc=1
c
      do k=1,n
        imax = idamax(n-k+1,a(k,k),inc)
        ip(k) = k+imax-1
        if (abs(a(ip(k),k)).le.1.e-10) then
          write(*,*)'k,ip(k) = ',k,ip(k)
          write(*,*)'K,PIVOT=', k,abs(a(ip(k),k))
          call pause
        endif
        pivmin = dmin1(pivmin,abs(a(ip(k),k)) )
        if (abs(a(ip(k),k)).ne.0.d0) go to 20
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
        if (iprint.eq.1) then
          write(*,7003) k
 7003     format('decomp: k = ',i4)
          do i=1,n
            write(*,7002) i,a(i,1:n)
          enddo
          call pause
        endif
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
c
      end subroutine




      subroutine gauss2(n,na,a,ip,b,x)
#include "syscom.blk"
      dimension a(na,n),b(n),x(n),ip(n)
c
      iprint=0
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
c
      inc=1
c
      do 10 k=1,n-1
      c=b(ip(k))
      b(ip(k))=b(k)
      b(k)=c
 10   call daxpy(n-k,-b(k),a(k+1,k),inc,b(k+1),inc)

      do 20 j=n,1,-1
      x(j)=b(j)/a(j,j)
 20   call daxpy(j-1,-x(j),a(1,j),inc,b,inc)
c
c
      end subroutine



c-----------------------------------------------------------------------
c
c   routine name       - gausse
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - april 2, 1988
c
c   purpose            - solution of a full matrix linear system
c                        of equations using gauss elimination
c                        without pivoting
c
c   usage              - call gausse(gk,igk,gf,u,n)
c
c   arguments in: gk   - n by n matrix
c                 igk  - row dimension of matrix gk exactly as
c                        specified in the dimension statement in the
c                        calling program
c                 gf   - vector of length n (the right-hand side)
c                 n    - number of equations
c             out:u    - vector of length n containing the solution
c
c   required routines  - rhsub,tri
c
c-----------------------------------------------------------------------
c
      subroutine gausse(gk,igk,gf,u,n)
c
#include "syscom.blk"
c
c
      dimension gk(igk,*),gf(*),u(*)
c
      call tri(gk,igk,n)
      call rhsub(gk,u,gf,igk,n)
c
      return
      end


      subroutine matinv(n,na,a,ip,Work)
#include "syscom.blk"
      dimension a(na,n),ip(n),Work(n)

      do 10 i=n,2,-1
      a(i,i)=1/a(i,i)
      call sscal (i-1, -a(i,i), a(1,i), 1)
      a(i-1,i)=a(i-1,i)/a(i-1,i-1)

      do 30 k=i-1,2,-1
      call saxpy (k-1, -a(k,i), a(1,k), 1, a(1,i), 1)
 30   a(k-1,i)=a(k-1,i)/a(k-1,k-1)

 10   continue

      a(1,1)=1/a(1,1)


      do 50 i=n-1,1,-1
      call scopy(   i, a(1,i), 1, Work     , 1)
      call scopy( n-i, 0.    , 0, Work(i+1), 1)

      do 70 k=i+1,n
 70   call saxpy( n, -a(k,i), a(1,k), 1, Work, 1)

      call scopy( n, Work, 1, a(1,i), 1)
      call sswap( n, a(1,i), 1, a(1,ip(i)), 1)

 50   continue

      return
      end


      function sdot1(n,sx,Incx,sy,Incy)
#include "syscom.blk"
      dimension sx(*),sy(*)

      sdot1=0.
      if(n.le.0) return

      jx=1
      if(Incx.lt.0) jx=1+(n-1)*(-Incx)
      jy=1
      if(Incy.lt.0) jy=1+(n-1)*(-Incy)

      do 10 i=1,n
      sdot1 = sdot1+sx(jx)*sy(jy)
      jx=jx+Incx
      jy=jy+Incy
 10   continue

      return
      end



      function snrm21(n,sx,Incx)
#include "syscom.blk"
      dimension sx(*)

      sn=0.
      jx=1
      if(Incx.lt.0) jx=1+(n-1)*(-Incx)

      do 10 i=1,n
      sn   =sn   +sx(jx)*sx(jx)
      jx=jx+Incx
 10   continue

      snrm21=sqrt(sn   )

      return
      end

