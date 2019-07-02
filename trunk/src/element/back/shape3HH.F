c--------------------------------------------------------------------
c
c     routine name      - shape3HH
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 14
c
c     purpose:          - routine returns values of 3D DISCONTINUOUS
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Type         - element type
c          Xi           - master element coordinates
c          Nord         - polynomial order for the middle node
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3HH(Type,Xi,Nord, NrdofH,ShapH,GradH)
c
      use parametersDPG
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(3), ShapH(MAXbrickHH),GradH(3,MAXbrickHH)
c

      select case(Type)
      case('bric','mdlb')
        call shape3HH_bric(Xi,Nord, NrdofH,ShapH,GradH)
      case('tetr','mdln')
        call shape3HH_tetr(Xi,Nord, NrdofH,ShapH,GradH)
      case('pris','mdlp')
        call shape3HH_pris(Xi,Nord, NrdofH,ShapH,GradH)
      case('pyra','mdld')
        call shape3HH_pyra(Xi,Nord, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape3HH: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3HH
c
c--------------------------------------------------------------------
c
c     routine name      - shape3HH_bric
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
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
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3HH_bric(Xi,Nord, NrdofH,ShapH,GradH)
c
      use parametersDPG
#include "syscom.blk"
c
      dimension Xi(3), ShapH(MAXbrickHH),GradH(3,MaxbrickHH)
c
c  ...1D shape functions
      dimension shapH1(MAXPP+1),dshapH1(MAXPP+1),
     .          shapH2(MAXPP+1),dshapH2(MAXPP+1),
     .          shapH3(MAXPP+1),dshapH3(MAXPP+1)
c
      iprint=0
c
c  ...initiate counter for shape functions
      k=0
c
c  ...shape functions are tensor products of 1D shape functions
      call decode(Nord, norda,nord3)
      call decode(norda, nord1,nord2)
      call shape1H(Xi(1),nord1, nrdofH1,shapH1,dshapH1)
      call shape1H(Xi(2),nord2, nrdofH2,shapH2,dshapH2)
      call shape1H(Xi(3),nord3, nrdofH3,shapH3,dshapH3)
c
ccc      shapH1(1) = 1.d0;       dshapH1(1) = 0.d0
ccc      shapH1(2) = Xi(1)-.5d0; dshapH1(2) = 1.d0
ccc      call shape1Hb(Xi(1),nord1,nv,shapH1(3:MAXPP+1),dshapH1(3:MAXPP+1))
ccc      nrdofH1 = nord1+1
c
ccc      shapH2(1) = 1.d0;       dshapH2(1) = 0.d0
ccc      shapH2(2) = Xi(2)-.5d0; dshapH2(2) = 1.d0
ccc      call shape1Hb(Xi(2),nord2,nv,shapH2(3:MAXPP+1),dshapH2(3:MAXPP+1))
ccc      nrdofH2 = nord2+1
c
ccc      shapH3(1) = 1.d0;       dshapH3(1) = 0.d0
ccc      shapH3(2) = Xi(3)-.5d0; dshapH3(2) = 1.d0
ccc      call shape1Hb(Xi(3),nord3,nv,shapH3(3:MAXPP+1),dshapH3(3:MAXPP+1))
ccc      nrdofH3 = nord3+1
c
c
c
      do l=1,nrdofH3
        do j=1,nrdofH2
          do i=1,nrdofH1
            k=k+1
            ShapH(k)   =  shapH1(i)* shapH2(j)* shapH3(l)
            GradH(1,k) = dshapH1(i)* shapH2(j)* shapH3(l)
            GradH(2,k) =  shapH1(i)*dshapH2(j)* shapH3(l)
            GradH(3,k) =  shapH1(i)* shapH2(j)*dshapH3(l)
          enddo
        enddo
      enddo
c
      NrdofH = k
c
c
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3)
 7001   format('shape3HH: Xi = ',3f8.3)
        do k=1,NrdofHH
          write(*,7002) k,ShapH(k),GradH(1:3,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shape3HH_bric
c
c--------------------------------------------------------------------
c
c     routine name      - shape3HH_tetr
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
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3HH_tetr(Xi,Nord, NrdofH,ShapH,GradH)
c
      use parametersDPG
      use element_data
#include "syscom.blk"
c
      dimension Xi(3), ShapH(MAXtetraHH),GradH(3,MAXtetraHH)
c
c
c
      end subroutine shape3HH_tetr
c
c--------------------------------------------------------------------
c
c     routine name      - shape3HH_pris
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
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH       - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3HH_pris(Xi,Nord, NrdofH,ShapH,GradH)
c
      use parametersDPG
      use element_data
#include "syscom.blk"
c
      dimension Xi(3), ShapH(MAXprismHH),GradH(3,MAXprismHH)
c
c
      end subroutine shape3HH_pris
c
c--------------------------------------------------------------------
c
c     routine name      - shape3HH_pyra
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
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3HH_pyra(Xip,Nord, NrdofH,ShapH,GradH)
c
      use parametersDPG
      use element_data
#include "syscom.blk"
c
      dimension Xip(3), ShapH(MAXpyramHH),GradH(3,MAXpyramHH)
c
c  ...common factor, modified coordinates
      dimension psi(4),gradpsi(3,4),alpha(2,4),gradalpha(3,2,4)

c
c
      end subroutine shape3HH_pyra
c





