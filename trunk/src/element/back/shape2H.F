c--------------------------------------------------------------------
c
c     routine name      - shape2H_1
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of a 2D element
c                         shape functions and their derivatives
c
c     arguments:
c
c     in:
c             Type      - element type
c             Xi        - master quad coordinates from (0,1)x(0,1)
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH     - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2H_1(Type,Xi,Norder,Norient, NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),Norder(5),Norient(4),
     .          ShapH((MAXP+1)**2),GradH(2,(MAXP+1)**2)
c
      select case(Type)
      case('tria','mdlt')
        call shape2H_tria(Xi,Norder,Norient, NrdofH,ShapH,GradH)
      case('quad','mdlq')
        call shape2H_quad(Xi,Norder,Norient, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape2H: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2H_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape2H_quad
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 2D quad element
c                         shape functions and their derivatives
c
c     arguments:
c
c     in:
c             Xi        - master quad coordinates from (0,1)x(0,1)
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH    - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2H_quad(Xi,Norder,Norient, NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(2),Norder(5),Norient(4),
     .          ShapH((MAXP+1)**2),GradH(2,(MAXP+1)**2)
c
c  ...1D vertex shape functions
      dimension shap1(2),dshap1(2),
     .          shap2(2),dshap2(2)
c
c  ...1D bubbles
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1)
c
c  ...blending functions and quad2edge projections
      dimension psi(4),gradpsi(2,4),te(4),gradte(2,4)
c
c  ...work space for shape_edgeH_orient
      dimension shap_edgeH(MAXP-1),dshap_edgeH(MAXP-1)
c
c  ...definition of 2D vertex shape function in terms of 1D shape
c     functions
      dimension ivert(2,4)
      data ivert /1,1, 2,1, 2,2, 1,2/
c
      iprint=0
c
      call shape1H(Xi(1),1, nvoid,shap1,dshap1)
      call shape1H(Xi(2),1, nvoid,shap2,dshap2)
c
      k=0
c
c  ...vertex shape functions
      do iv=1,4
        i = ivert(1,iv); j = ivert(2,iv)
        k=k+1
        ShapH(k)   =  shap1(i)* shap2(j)
        GradH(1,k) = dshap1(i)* shap2(j)
        GradH(2,k) =  shap1(i)*dshap2(j)
      enddo
c
c  ...edge shape functions
      call blend_project_quad(Xi, psi,gradpsi,te,gradte)
      do ie=1,4
        call shape1Hb_orient(te(ie),Norder(ie),Norient(ie),
     .                       nrdof_edgeH,shap_edgeH,dshap_edgeH)
        do ke=1,nrdof_edgeH
          k=k+1
          ShapH(k) = psi(ie)*shap_edgeH(ke)
          GradH(1:2,k) = gradpsi(1:2,ie)*shap_edgeH(ke)
     .                 + psi(ie)*dshap_edgeH(ke)*gradte(1:2,ie)
        enddo
      enddo
c
c  ...bubbles = tensor products of 1D bubbles
      call decode(Norder(5), nord1,nord2)
      call shape1Hb(Xi(1),nord1, nrdofH1,shapH1,dshapH1)
      call shape1Hb(Xi(2),nord2, nrdofH2,shapH2,dshapH2)
c
      do j=1,nrdofH2
        do i=1,nrdofH1
          k=k+1
          ShapH(k)   =  shapH1(i)* shapH2(j)
          GradH(1,k) = dshapH1(i)* shapH2(j)
          GradH(2,k) =  shapH1(i)*dshapH2(j)
        enddo
      enddo
      NrdofH=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:2),Norder(1:5),Norient(1:4)
 7001   format('shape2H_quad: Xi = ',2f8.3,/,
     .         'Norder  = ',4i2,i3,/,
     .         'Norient = ',4i2)
        write(*,7002)
 7002   format('VERTEX SHAPE FUNCTIONS = ')
        do iv=1,4
          k=iv
          write(*,7003) k,ShapH(k),GradH(1:2,k)
 7003     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,2e12.5)
        enddo
        do ie=1,4
          ndof = Norder(ie)-1
          if (ndof.gt.0) then
            write(*,7004) ie
 7004       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapH(k),GradH(1:2,k)
            enddo
          endif
        enddo
        call decode(Norder(5), nord1,nord2)
        ndof = (nord1-1)*(nord2-1)
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapH(k),GradH(1:2,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape2H_quad
c
c--------------------------------------------------------------------
c
c     routine name      - shape2H_tria
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 2D tria element
c                         shape functions and their derivatives
c
c     arguments:
c
c     in:
c             Xi        - master tria coordinates from (0,1)x(0,1)
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH    - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2H_tria(Xi,Norder,Norient, NrdofH,ShapH,GradH)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(2),Norder(4),Norient(3),
     .          ShapH((MAXP+1)*(MAXP+2)/2),GradH(2,(MAXP+1)*(MAXP+2)/2)
c
      iprint=0
c
c  ...barycentric coordinates = vertex shape functions
      ShapH(1) = 1.d0-Xi(1)-Xi(2); GradH(1:2,1)=-1.d0
      ShapH(2) = Xi(1); GradH(1,2) = 1.d0; GradH(2,2) = 0.d0
      ShapH(3) = Xi(2); GradH(1,3) = 0.d0; GradH(2,3) = 1.d0
      if (Norder(4).eq.1) go to 999
c
      k=3
c
c  ...edge shape functions
      do ie=1,3
        nv1 = TRIAN_EDGE_TO_VERT(1,ie)
        nv2 = TRIAN_EDGE_TO_VERT(2,ie)
        call shape_edgeH_orient((/ShapH(nv1),ShapH(nv2)/),
     .                          (/GradH(1:2,nv1),GradH(1:2,nv2)/),
     .                          2,Norder(ie),Norient(ie),
     .                          nrdof_edgeH,ShapH(k+1),GradH(1:2,k+1))
        k=k+nrdof_edgeH
      enddo
c
c  ...bubbles
      kk=k
      call shape_triaH(ShapH(1:3),GradH(1:2,1:3),2,Norder(4),
     .                 NrdofH,ShapH(k+1),GradH(1:2,k+1))
      NrdofH = kk+NrdofH
c
c
 999  continue
      if (iprint.eq.1) then
        write(*,7001) Xi(1:2),Norder(1:4),Norient(1:3)
 7001   format('shape2H_tria: Xi = ',2f8.3,/,
     .         'Norder  = ',3i2,2xi2,/,
     .         'Norient = ',3i2)
        write(*,7002)
 7002   format('VERTEX SHAPE FUNCTIONS = ')
        do iv=1,3
          k=iv
          write(*,7003) k,ShapH(k),GradH(1:2,k)
 7003     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,2e12.5)
        enddo
        do ie=1,3
          ndof = Norder(ie)-1
          if (ndof.gt.0) then
            write(*,7004) ie
 7004       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapH(k),GradH(1:2,k)
            enddo
          endif
        enddo
        nord = Norder(4)
        ndof = (nord-1)*(nord-2)/2
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapH(k),GradH(1:2,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape2H_tria
c
c--------------------------------------------------------------------
c
c     routine name      - shape2Hb_orient
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 14
c
c     purpose:          - routine returns values of a 2D element
c                         bubble shape functions and their derivatives
c                         accounting for orientation
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
c             NrdofH    - number of dof
c             ShapH     - values of the shape functions at the point
c             GradH     - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2Hb_orient(Type,Xi,Nord,Norient,
     .                           NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),
     .          ShapH(MAXmdlqH),GradH(2,MAXmdlqH)
c
c  ...affine coordinates for the master triangle
      dimension vshap(0:2), vgrad(2,0:2)
c
c  ...work space for shape1Hb_orient
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1)
c
c  ...translating quad orientation into product of 1D orientation
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      select case(Type)
      case('tria','mdlt')
        vgrad = 0.d0
        vshap(0) = 1.d0-Xi(1)-Xi(2); vgrad(1:2,0) = -1.d0
        vshap(1) = Xi(1); vgrad(1,1) = 1.d0
        vshap(2) = Xi(2); vgrad(2,2) = 1.d0
        call shape_triaH_orient(vshap,vgrad,2,Nord,Norient,
     .                          NrdofH,ShapH,GradH)
      case('quad','mdlq','rect')
        call decode(Nord, nord1,nord2)
        NrdofH = (nord1-1)*(nord2-1)
        if (NrdofH.eq.0) return
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
        k=0
        do j=1,nrdofH2
          do i=1,nrdofH1
            k=k+1
            ShapH(k) = shapH1(i)*shapH2(j)
            select case(iswap)
            case(0)
              GradH(1,k) = dshapH1(i)* shapH2(j)
              GradH(2,k) =  shapH1(i)*dshapH2(j)
            case(1)
              GradH(1,k) =  shapH1(i)*dshapH2(j)
              GradH(2,k) = dshapH1(i)* shapH2(j)
            end select
          enddo
        enddo
      case default
        write(*,*)'shape2Hb_orient: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2Hb_orient
c




