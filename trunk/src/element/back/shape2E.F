c--------------------------------------------------------------------
c
c     routine name      - shape2E_1
c
c--------------------------------------------------------------------
c
c     latest revision:  - Dec 13
c
c     purpose:          - routine returns values of a 2D element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c             Type      - element type
c             Xi        - master quad coordinates
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2E_1(Type,Xi,Norder,Norient, NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),Norder(5),Norient(4),
     .          ShapE(2,2*(MAXP+1)*MAXP),CurlE(2*(MAXP+1)*MAXP)
c
      select case(Type)
      case('tria','mdlt')
        call shape2E_tria(Xi,Norder,Norient, NrdofE,ShapE,CurlE)
      case('quad','mdlq')
        call shape2E_quad(Xi,Norder,Norient, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape2E: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2E_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape2E_quad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 2D quad element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c             Xi        - master quad coordinates from (0,1)x(0,1)
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2E_quad(Xi,Norder,Norient, NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(2),Norder(5),Norient(4),
     .          ShapE(2,2*(MAXP+1)*MAXP),CurlE(2*(MAXP+1)*MAXP)
c
c  ...blending functions and quad2edge projections
      dimension psi(4),gradpsi(2,4),te(4),gradte(2,4)
c
c  ...work space for shape_edgeQ_orient
      dimension shap_edgeQ(MAXP)
cc
c  ...work space for routines shape1Hb, shape1Q
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1),
     .          shapQ1(MAXP),shapQ2(MAXP)

      iprint=1
c
      k=0
c
c  ...edge shape functions
      call blend_project_quad(Xi, psi,gradpsi,te,gradte)
      do ie=1,4
        call shape1Q_orient(te(ie),Norder(ie),Norient(ie),
     .                      nrdof_edgeQ,shap_edgeQ)
        curl = gradpsi(1,ie)*gradte(2,ie)-gradpsi(2,ie)*gradte(1,ie)
        do ke=1,nrdof_edgeQ
          k=k+1
          ShapE(1:2,k) = psi(ie)*shap_edgeQ(ke)*gradte(1:2,ie)
          CurlE(k) = curl*shap_edgeQ(ke)
        enddo
      enddo
c
c  ...bubbles
      call decode(Nord, nord1,nord2)
      call shape1Hb(Xi(1),nord1, nrdofH1,shapH1,dshapH1)
      call shape1Hb(Xi(2),nord2, nrdofH2,shapH2,dshapH2)
      call shape1Q(Xi(1),nord1, nrdofQ1,shapQ1)
      call shape1Q(Xi(2),nord2, nrdofQ2,shapQ2)
c
c  ...horizontal shape functions
      do j=1,nrdofH2
        do i=1,nrdofQ1
          k=k+1
          ShapE(1,k) = shapQ1(i)*shapH2(j)
          ShapE(2,k) = 0.d0
          CurlE(k)   = - shapQ1(i)*dshapH2(j)
        enddo
      enddo
c
c  ...vertical shape functions
      do j=1,nrdofQ2
        do i=1,nrdofH1
          k=k+1
          ShapE(1,k) = 0.d0
          ShapE(2,k) = shapH1(i)*shapQ2(j)
          CurlE(k)   = dshapH1(i)*shapQ2(j)
        enddo
      enddo
c
      NrdofE=k
c
c
      if (iprint.eq.1) then
        write(*,7000) Xi(1:2),Norder(1:5),Norient(1:4)
 7000   format('shape2E_quad: Xi = ',2f8.3,/,
     .         'Norder  = ',4i2,i3,/,
     .         'Norient = ',4i2)
        k=0
        do ie=1,4
          ndof = Norder(ie)
          if (ndof.gt.0) then
            write(*,7001) ie
 7001       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapE(1:2,k),CurlE(k)
 7002         format('k = ',i3,' ShapE, CurlE = ',2e12.5,3x,e12.5)
            enddo
          endif
        enddo
        call decode(Norder(5), nord1,nord2)
        ndof = nord1*(nord2-1) + (nord1-1)*nord2
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapE(1:2,k),CurlE(k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape2E_quad
c
c--------------------------------------------------------------------
c
c     routine name      - shape2E_tria
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 2D triangular
c                         element H(curl) shape functions and their
c                         curl (Nedelec triangle fo the first type)
c
c     arguments:
c
c     in:
c             Xi        - master tria coordinates from (0,1)x(0,1)
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2E_tria(Xi,Norder,Norient, NrdofE,ShapE,CurlE)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(2),Norder(4),Norient(3),
     .          ShapE(2,MAXtriaE),CurlE(MAXtriaE)
c
c  ...barycentric coordinates and their derivatives
      dimension vshap(3),vgrad(1:2,3)
c
      iprint=0
c
c  ...barycentric coordinates = vertex shape functions
      vshap(1) = 1.d0-Xi(1)-Xi(2); vgrad(1:2,1)=-1.d0
      vshap(2) = Xi(1); vgrad(1,2) = 1.d0; vgrad(2,2) = 0.d0
      vshap(3) = Xi(2); vgrad(1,3) = 0.d0; vgrad(2,3) = 1.d0
c
      k=0
c
c  ...edge shape functions
      do ie=1,3
        nv1 = TRIAN_EDGE_TO_VERT(1,ie)
        nv2 = TRIAN_EDGE_TO_VERT(2,ie)
        call shape_edgeE_orient((/vshap(nv1),vshap(nv2)/),
     .                          (/vgrad(1:2,nv1),vgrad(1:2,nv2)/),
     .                          2,Norder(ie),Norient(ie),
     .                          nrdof_edgeE,ShapE(1:2,k+1),CurlE(k+1))
        k=k+nrdof_edgeE
      enddo
c
c  ...bubbles
      kk=k
      call shape_triaE(vshap,vgrad,2,Norder(4),
     .                 NrdofE,ShapE(1:2,k+1),CurlE(k+1))
      NrdofE = kk+NrdofE
c
c
      if (iprint.eq.1) then
        write(*,7000) Xi(1:2),Norder(1:4),Norient(1:3)
 7000   format('shape2E_tria: Xi = ',2f8.3,/,
     .         'Norder  = ',3i2,i3,/,
     .         'Norient = ',3i2)
        k=0
        do ie=1,3
          ndof = Norder(ie)
          if (ndof.gt.0) then
            write(*,7001) ie
 7001       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapE(1:2,k),CurlE(k)
 7002         format('k = ',i3,' ShapE, CurlE = ',2e12.5,3x,e12.5)
            enddo
          endif
        enddo
        nord = Norder(4)
        ndof = nord*(nord-1)
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapE(1:2,k),CurlE(k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape2E_tria
c
c--------------------------------------------------------------------
c
c     routine name      - shape2Eb_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 14
c
c     purpose:          - routine returns values of a 2D element
c                         Hcurl bubble shape functions and their curl
c                         accounting for orientation (neede for PB
c                         interpolation of Dirichlet data)
c
c     arguments:
c
c     in:
c             Type      - element type
c             Xi        - master element coordinates
c             Nord      - polynomial order (of the middle node)
c             Norient   - element orientation
c
c     out:
c             NrdofE    - number of dof
c             ShapE     - values of the shape functions at the point
c             CurlE     - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2Eb_orient(Type,Xi,Nord,Norient,
     .                           NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),
     .          ShapE(2,MAXmdlqE),CurlE(MAXmdlqE)
c
c  ...affine coordinates for the master triangle
      dimension vshap(0:2), vgrad(2,0:2)
c
c  ...work space for shape1Hb_orient
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1),
     .          shapQ1(MAXP),shapQ2(MAXP)
c
c  ...translating quad orientation into product of 1D orientation
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
      select case(Type)
      case('tria','mdlt')
        vgrad = 0.d0
        vshap(0) = 1.d0-Xi(1)-Xi(2); vgrad(1:2,0) = -1.d0
        vshap(1) = Xi(1); vgrad(1,1) = 1.d0
        vshap(2) = Xi(2); vgrad(2,2) = 1.d0
        call shape_triaE_orient(vshap,vgrad,2,Nord,Norient,
     .                          NrdofE,ShapE,CurlE)
      case('quad','mdlq','rect')
        call decode(Nord, nord1,nord2)
        NrdofE = (nord1-1)*nord2 + nord1*(nord2-1)
        if (NrdofE.eq.0) return
c
        iswap = norient221(1,Norient)
        norient1 = norient221(2,Norient)
        norient2 = norient221(3,Norient)
        select case(iswap)
c
c  .....unswapped coordinates
        case(0)
          xil1 = Xi(1); xil2 = Xi(2)
c
c  .....swapped coordinates
        case(1)
          xil1 = Xi(2); xil2 = Xi(1)
        end select
        call shape1Hb_orient(xil1,nord1,norient1,
     .                       nrdofH1,shapH1,dshapH1)
        call shape1Hb_orient(xil2,nord2,norient2,
     .                       nrdofH2,shapH2,dshapH2)
c
        call shape1Q_orient(xil1,nord1,norient1, nrdofQ1,shapQ1)
        call shape1Q_orient(xil2,nord2,norient2, nrdofQ2,shapQ2)
c
        k=0
c
c  .....horizontal shape functions
        do j=1,nrdofH2
          do i=1,nrdofQ1
            k=k+1
            select case(iswap)
            case(0)
              shapE(1,k) = shapQ1(i)*shapH2(j)
              shapE(2,k) = 0.d0
              curlE(k)   = - shapQ1(i)*dshapH2(j)
            case(1)
              shapE(2,k) = shapQ1(i)*shapH2(j)
              shapE(1,k) = 0.d0
              curlE(k)   =   shapQ1(i)*dshapH2(j)
            end select
          enddo
        enddo
c
c  .....vertical shape functions
        do j=1,nrdofQ2
          do i=1,nrdofH1
            k=k+1
            select case(iswap)
            case(0)
              ShapE(1,k) = 0.d0
              ShapE(2,k) = shapH1(i)*shapQ2(j)
              CurlE(k)   = dshapH1(i)*shapQ2(j)
            case(1)
              ShapE(2,k) = 0.d0
              ShapE(1,k) = shapH1(i)*shapQ2(j)
              CurlE(k)   = -dshapH1(i)*shapQ2(j)
            end select
          enddo
        enddo
        NrdofE = k
      end select
c
      if (iprint.eq.1) then
        write(*,7000) Xi(1:2),Nord,Norient
 7000   format('shape2Eb_orient: Xi = ',2f8.3,/,
     .         'Nord  = ',i3,' Norient = ',i2)
        do k=1,NrdofE
          write(*,7002) k,ShapE(1:2,k),CurlE(k)
 7002     format('k = ',i3,' ShapE = ',2e12.5,' CurlE = ', e12.5)
        enddo
        call pause
      endif
c
      end subroutine shape2Eb_orient
c







