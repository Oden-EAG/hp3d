c--------------------------------------------------------------------
c
c     routine name      - shape3EE
c
c--------------------------------------------------------------------
c
c     latest revision:  - May 14
c
c     purpose:          - routine returns values of 3D DISCONTINUOUS
c                         Hcurl shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Type         - element type
c          Xi           - master element coordinates
c          Nord         - polynomial order for the middle node
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3EE(Type,Xi,Nord, NrdofE,ShapE,CurlE)
c
      use parametersDPG
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(3), ShapE(3,MAXbrickEE),CurlE(3,MAXbrickEE)
c

      select case(Type)
      case('bric','mdlb')
        call shape3EE_bric(Xi,Nord, NrdofE,ShapE,CurlE)
      case('tetr','mdln')
        call shape3EE_tetr(Xi,Nord, NrdofE,ShapE,CurlE)
      case('pris','mdlp')
        call shape3EE_pris(Xi,Nord, NrdofE,ShapE,CurlE)
      case('pyra','mdld')
        call shape3EE_pyra(Xi,Nord, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape3EE: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3EE
c
c--------------------------------------------------------------------
c
c     routine name      - shape3EE_bric
c
c--------------------------------------------------------------------
c
c     latest revision:  - May 14
c
c     purpose:          - routine returns values of 3D brick element
c                         DISCONTINOUS H1 shape functions and their
c                         derivatives
c
c     arguments:
c
c     in:
c          Xi           - master brick coordinates from (0,1)^3
c          Nord         - middle node order
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3EE_bric(Xi,Nord, NrdofE,ShapE,CurlE)
c
      use parametersDPG
#include "syscom.blk"
c
      dimension Xi(3), ShapE(3,MAXbrickEE),CurlE(3,MaxbrickEE)
c
c  ...1D shape functions
      dimension shapH1(MAXPP+1),dshapH1(MAXPP+1),
     .          shapH2(MAXPP+1),dshapH2(MAXPP+1),
     .          shapH3(MAXPP+1),dshapH3(MAXPP+1),
     .          shapQ1(MAXP),shapQ2(MAXP),shapQ3(MAXP)
c
c
      iprint=0
c
c  ...initiate counter for shape functions
      k=0
c
c  ...shape functions are tensor products of 1D shape functions
      call decode(Nord, norda,nord3)
      call decode(norda, nord1,nord2)
c
      call shape1H(Xi(1),nord1, nrdofH1,shapH1,dshapH1)
      call shape1H(Xi(2),nord2, nrdofH2,shapH2,dshapH2)
      call shape1H(Xi(3),nord3, nrdofH3,shapH3,dshapH3)
c
      call shape1Q(Xi(1),nord1, nrdofQ1,shapQ1)
      call shape1Q(Xi(2),nord2, nrdofQ2,shapQ2)
      call shape1Q(Xi(3),nord3, nrdofQ3,shapQ3)
c
c
c  ...shape functions with values along the x-axis
      do l=1,nrdofH3
        do j=1,nrdofH2
          do i=1,nrdofQ1
            k=k+1
            ShapE(1,k) = shapQ1(i)*shapH2(j)*shapH3(l)
            ShapE(2,k) = 0.d0
            ShapE(3,k) = 0.d0
            CurlE(1,k)   = 0.d0
            CurlE(2,k)   =  shapQ1(i)* shapH2(j)*dshapH3(l)
            CurlE(3,k)   = -shapQ1(i)*dshapH2(j)* shapH3(l)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the y-axis
      do l=1,nrdofH3
        do j=1,nrdofQ2
          do i=1,nrdofH1
            k=k+1
            ShapE(1,k) = 0.d0
            ShapE(2,k) = shapH1(i)*shapQ2(j)*shapH3(l)
            ShapE(3,k) = 0.d0
            CurlE(1,k)   = -shapH1(i)* shapQ2(j)*dshapH3(l)
            CurlE(2,k)   =  0.d0
            CurlE(3,k)   = dshapH1(i)* shapQ2(j)* shapH3(l)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the z-axis
      do l=1,nrdofQ3
        do j=1,nrdofH2
          do i=1,nrdofH1
            k=k+1
            ShapE(1,k) = 0.d0
            ShapE(2,k) = 0.d0
            ShapE(3,k) = shapH1(i)*shapH2(j)*shapQ3(l)
            CurlE(1,k)   =   shapH1(i)*dshapH2(j)*shapQ3(l)
            CurlE(2,k)   = -dshapH1(i)* shapH2(j)*shapQ3(l)
            CurlE(3,k)   = 0.d0
          enddo
        enddo
      enddo
c
      NrdofE = k
c
c
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3EE: Xi = ',3f8.3)
        do k=1,NrdofEE
          write(*,7002) k,ShapE(1:3,k),CurlE(1:3,k)
 7002     format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape3EE_bric
c
c--------------------------------------------------------------------
c
c     routine name      - shape3EE_tetr
c
c--------------------------------------------------------------------
c
c     latest revision:  -
c
c     purpose:          - routine returns values of 3D tetra element
c                         DISCONTINUOUS H1 shape functions and their
c                         derivatives
c
c     arguments:
c
c     in:
c          Xi           - master tetra coordinates
c          Nord         - middle node order
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3EE_tetr(Xi,Nord, NrdofE,ShapE,CurlE)
c
      use parametersDPG
      use element_data
#include "syscom.blk"
c
      dimension Xi(3), ShapE(3,MAXtetraEE),CurlE(3,MAXtetraEE)
c
c
c
      end subroutine shape3EE_tetr
c
c--------------------------------------------------------------------
c
c     routine name      - shape3EE_pris
c
c--------------------------------------------------------------------
c
c     latest revision:  -
c
c     purpose:          - routine returns values of 3D prism element
c                         DISCONTINUOUS H1 shape functions and their
c                         derivatives
c
c     arguments:
c
c     in:
c          Xi           - master prisk coordinates from (0,1)^3
c          Nord         - middle node order
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE       - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3EE_pris(Xi,Nord, NrdofE,ShapE,CurlE)
c
      use parametersDPG
      use element_data
#include "syscom.blk"
c
      dimension Xi(3), ShapE(3,MAXprismEE),CurlE(3,MAXprismEE)
c
c
      end subroutine shape3EE_pris
c
c--------------------------------------------------------------------
c
c     routine name      - shape3EE_pyra
c
c--------------------------------------------------------------------
c
c     latest revision:  -
c
c     purpose:          - routine returns values of 3D pyramid element
c                         DISCONTINUOUS H1 shape functions and their
c                         derivatives
c
c     arguments:
c
c     in:
c          Xip          - master pyramid coordinates
c          Nord         - middle node order
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3EE_pyra(Xip,Nord, NrdofE,ShapE,CurlE)
c
      use parametersDPG
      use element_data
#include "syscom.blk"
c
      dimension Xip(3), ShapE(3,MAXpyramEE),CurlE(3,MAXpyramEE)
c
c  ...common factor, modified coordinates
      dimension psi(4),gradpsi(3,4),alpha(2,4),gradalpha(3,2,4)

c
c
      end subroutine shape3EE_pyra
c





