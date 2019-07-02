c----------------------------------------------------------------------
c
c     routine name      - shape_edgeE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of edge H(curl) shape
c                         functions in edge affine coordinates
c
c     arguments:
c
c     in:
c             Vshap     - edge vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             N         - space dimension
c             Nord      - polynomial order (for H1 element)
c
c     out:
c             NrdofE    - number of shape functions
c             ShapE     - values of the shape functions at the point
c             CurlE     - the corresponding curl
c
c-----------------------------------------------------------------------
c
      subroutine shape_edgeE(Vshap,Vgrad,N,Nord, NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:1),Vgrad(N,0:1),
     .          ShapE(N,MAXP),CurlE(2*N-3,MAXP)
c
c  ...work space for routine Legendre
      dimension P(0:MAXP)
c
c  ...Whitney shape function, its curl (half)
      dimension vecW(N),crossW(2*N-3)
c
      iprint=0
c
c  ...dimension of curl
      m = 2*N-3
c
c  ...Whitney shape function
      vecW(1:N) = Vshap(0)*Vgrad(1:N,1) - Vshap(1)*Vgrad(1:N,0)
c
c  ...related curl
      call cross(N,Vgrad(1:N,0),Vgrad(1:N,1), crossW)
      if (iprint.eq.1) then
        write(*,*) 'Vgrad(1:N,0) = ',Vgrad(1:N,0)
        write(*,*) 'Vgrad(1:N,1) = ',Vgrad(1:N,1)
        write(*,*) 'crossW = ',crossW(1:m)
      endif
c
c  ...compute shifted scaled Legendre polynomials
      x = Vshap(1); t = Vshap(0) + Vshap(1)
      call LegendreNEW(x,t,Nord-1, P)
c
      k=0
      do np=0,Nord-1
        k=k+1
        ShapE(1:N,k) = P(np)*vecW(1:N)
        CurlE(1:m,k) = (np+2)*P(np)*crossW(1:m)
      enddo
      NrdofE = k
c
c
      if (iprint.eq.1) then
        write(*,7001) Nord, Vshap(0:1)
 7001   format('shape_edgeE: Nord = ',i2,' Vshap = ',2F8.3)
        do k=1,NrdofE
          write(*,7002) k,ShapE(1:N,k),CurlE(1:m,k)
 7002     format(' k = ',i2,' ShapE,CurlE = ',6e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape_edgeE
c
c----------------------------------------------------------------------
c
c     routine name      - shape_edgeE_orient
c
c----------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of edge H(curl) shape
c                         functions in edge affine coordinates
c                         accounting for orientation
c
c     arguments:
c
c     in:
c             Vshap     - edge vertex shape functions (affine
c                         coordinates)
c             Vgrad     - the corresponding gradients
c             N         - space dimension
c             Nord      - polynomial order (for H1 element)
c             Norient   - edge orientation
c
c     out:
c             NrdofE    - number of shape functions
c             ShapE     - values of the shape functions at the point
c             CurlE     - the corresponding curl
c
c-----------------------------------------------------------------------
c
      subroutine shape_edgeE_orient(Vshap,Vgrad,N,Nord,Norient,
     .                              NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Vshap(0:1),Vgrad(N,0:1),
     .          ShapE(N,MAXP),CurlE(2*N-3,MAXP)
c
c     edge affine coordinates
      dimension vshap_loc(0:1),vgrad_loc(N,0:1)
c
      iprint=0
c
      select case(Norient)
      case(0)
        n0=0; n1=1
      case(1)
        n0=1; n1=0
      end select
      vshap_loc(0) = Vshap(n0); vgrad_loc(1:N,0) = Vgrad(1:N,n0)
      vshap_loc(1) = Vshap(n1); vgrad_loc(1:N,1) = Vgrad(1:N,n1)
c
      call shape_edgeE(vshap_loc,vgrad_loc,N,Nord, NrdofE,ShapE,CurlE)
c
      if (iprint.eq.1) then
        write(*,7001) Nord, Norient, Vshap(0:1)
 7001   format('shape_edgeE_orient: Nord,Norient = ',2i2,
     .         ' Vshap = ',2F8.3)
        m = 2*N-3
        do k=1,NrdofE
          write(*,7002) k,ShapE(1:N,k),CurlE(1:m,k)
 7002     format(' k = ',i2,' ShapE,CurlE = ',6e12.5)
        enddo
        call pause
      endif
c
c
c
      end subroutine shape_edgeE_orient


