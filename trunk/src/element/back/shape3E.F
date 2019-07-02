c--------------------------------------------------------------------
c
c     routine name      - shape3E
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of 3D element
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
      subroutine shape3E_1(Type,Xi,Norder,Norient_edge,Norient_face,
     .                     NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(3),Norder(19),Norient_edge(12),Norient_face(6),
     .          ShapE(3,MAXbrickE),CurlE(3,MAXbrickE)
c
      select case(Type)
      case('bric','mdlb')
        call shape3E_bric(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofE,ShapE,CurlE)
      case('tetr','mdln')
        call shape3E_tetr(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofE,ShapE,CurlE)
      case('pris','mdlp')
        call shape3E_pris(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofE,ShapE,CurlE)
      case('pyra','mdld')
        call shape3E_pyra(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape3E: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3E_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape3E_bric
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D brick element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c          Xi           - master bric coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3E_bric(Xi,Norder,Norient_edge,Norient_face,
     .                        NrdofE,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),Norder(19),Norient_edge(12),Norient_face(6),
     .          ShapE(3,MaxbrickE),CurlE(3,MaxbrickE)
c
c  ...blending functions and projections
      dimension psie(12),gradpsie(3,12),te(12),gradte(3,12),
     .          psif(6),gradpsif(3,6),tf(2,6),gradtf(2,3,6)
c
c  ...work space for shape_edgeQ_orient
      dimension shap_edgeQ(MAXP),vec(3)
c
c  ...work space for face bubbles
      dimension shapEf(2)
c
c  ...work space for faces
      dimension vect(3,2)
c
c  ...work space for routines shape1Hb, shape1Q
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1),
     .          shapH3(MAXP-1),dshapH3(MAXP-1),
     .          shapQ1(MAXP),shapQ2(MAXP),shapQ3(MAXP)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
      k=0
c
c  ...determine edge blending functions and projections
      call blend_project_edges_bric(Xi, psie,gradpsie,te,gradte)
c
c  ...edge shape functions
      do ie=1,12
        call shape1Q_orient(te(ie),Norder(ie),Norient_edge(ie),
     .                      nrdof_edgeQ,shap_edgeQ)
        call cross_product(gradpsie(1:3,ie),gradte(1:3,ie), vec)
        do ke=1,nrdof_edgeQ
          k=k+1
          ShapE(1:3,k) = psie(ie)*shap_edgeQ(ke)*gradte(1:3,ie)
          CurlE(1:3,k) = vec(1:3)*shap_edgeQ(ke)
        enddo
      enddo
c
c  ...determine face blending functions and projections
      call blend_project_faces_bric(Xi, psif,gradpsif,tf,gradtf)
c
c  ...face shape functions
      do if=1,6
c
c  .....order is returned in the element system of coordinates
        call decode(Norder(12+if), nord1,nord2)
        ndof = nord1*(nord2-1)+(nord1-1)*nord2
        if (ndof.eq.0) cycle
c
        norient1 = norient221(2,Norient_face(if))
        norient2 = norient221(3,Norient_face(if))
        iswap = norient221(1,Norient_face(if))
        select case(iswap)
c
c  .....unswapped coordinates
        case(0)
          xif1 = tf(1,if); xif2 = tf(2,if)
c
c  .....swapped coordinates
        case(1)
          xif1 = tf(2,if); xif2 = tf(1,if)
        end select
c
        call shape1Hb_orient(xif1,nord1,norient1,
     .                       nrdofH1,shapH1,dshapH1)
        call shape1Hb_orient(xif2,nord2,norient2,
     .                       nrdofH2,shapH2,dshapH2)
        call shape1Q_orient(xif1,nord1,norient1, nrdofQ1,shapQ1)
        call shape1Q_orient(xif2,nord2,norient2, nrdofQ2,shapQ2)
c
        call cross_product(gradpsif(1:3,if),gradtf(1,1:3,if),
     .                     vect(1:3,1))
        call cross_product(gradpsif(1:3,if),gradtf(2,1:3,if),
     .                     vect(1:3,2))
        call cross_product(gradtf(1,1:3,if),gradtf(2,1:3,if), vec)
c
c  .....horizontal (in face owned coordinates) shape functions
        do j=1,nrdofH2
          do i=1,nrdofQ1
            k=k+1
            select case(iswap)
            case(0)
              shapEf(1) = shapQ1(i)*shapH2(j)
              shapEf(2) = 0.d0
              curlEf    = - shapQ1(i)*dshapH2(j)
            case(1)
              shapEf(2) = shapQ1(i)*shapH2(j)
              shapEf(1) = 0.d0
              curlEf    =   shapQ1(i)*dshapH2(j)
            end select
            ShapE(1:3,k) = psif(if)*(shapEf(1)*gradtf(1,1:3,if)
     .                              +shapEf(2)*gradtf(2,1:3,if))

            CurlE(1:3,k) = shapEf(1)*vect(1:3,1)
     .                   + shapEf(2)*vect(1:3,2)
     .                   + psif(if)*curlEf*vec(1:3)
          enddo
        enddo
c
c  .....vertical (in face owned coordinates) shape functions
        do j=1,nrdofQ2
          do i=1,nrdofH1
            k=k+1
            select case(iswap)
            case(0)
              ShapEf(1) = 0.d0
              ShapEf(2) = shapH1(i)*shapQ2(j)
              CurlEf    = dshapH1(i)*shapQ2(j)
            case(1)
              ShapEf(2) = 0.d0
              ShapEf(1) = shapH1(i)*shapQ2(j)
              CurlEf    = -dshapH1(i)*shapQ2(j)
            end select
            ShapE(1:3,k) = psif(if)*(shapEf(1)*gradtf(1,1:3,if)
     .                              +shapEf(2)*gradtf(2,1:3,if))

            CurlE(1:3,k) = shapEf(1)*vect(1:3,1)
     .                   + shapEf(2)*vect(1:3,2)
     .                   + psif(if)*curlEf*vec(1:3)
          enddo
        enddo
      enddo
c
c  ...bubbles...........................................................
      call decode(Norder(19), norda,nord3)
      call decode(norda, nord1,nord2)
      call shape1Hb(Xi(1),nord1, nrdofH1,shapH1,dshapH1)
      call shape1Hb(Xi(2),nord2, nrdofH2,shapH2,dshapH2)
      call shape1Hb(Xi(3),nord3, nrdofH3,shapH3,dshapH3)
      call shape1Q(Xi(1),nord1, nrdofQ1,shapQ1)
      call shape1Q(Xi(2),nord2, nrdofQ2,shapQ2)
      call shape1Q(Xi(3),nord3, nrdofQ3,shapQ3)
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
      NrdofE=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(1:19),
     .                Norient_edge(1:12),Norient_face(1:6)
 7001   format('shape3E_bric: Xi = ',3f8.3,/,
     .         'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/,
     .         'orient = ',3(4i2,2x),2i3,2x,4i3)
        k=0
        do ie=1,12
          ndof = Norder(ie)
          if (ndof.gt.0) then
            write(*,7002) ie
 7002       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
 7003         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        do if=1,6
          call decode(Norder(12+if), nordh,nordv)
          ndof = (nordh-1)*nordv + nordh*(nordv-1)
          if (ndof.gt.0) then
            write(*,7004) if
 7004       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
            enddo
          endif
        enddo
        call decode(Norder(19), nord,nord3)
        call decode(nord, nord1,nord2)
        ndof =  nord1   *(nord2-1)*(nord3-1)
     .       + (nord1-1)* nord2   *(nord3-1)
     .       + (nord1-1)*(nord2-1)* nord3
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3E_bric
c
c--------------------------------------------------------------------
c
c     routine name      - shape3E_tetr
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D tetra element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c          Xi           - master tetr coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3E_tetr(Xi,Norder,Norient_edge,Norient_face,
     .                        NrdofE,ShapE,CurlE)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(3),Norder(11),Norient_edge(6),Norient_face(4),
     .          ShapE(3,MAXtetraE),CurlE(3,MAXtetraE)
c
c  ...barycentric coordinates and their derivatives
      dimension vshap(4),vgrad(1:3,4)
c
      iprint=0
c
c  ...barycentric coordinates = vertex shape functions
      vshap(1) = 1.d0-Xi(1)-Xi(2)-Xi(3); vgrad(1:3,1)=-1.d0
      vgrad(1:3,2:4) = 0.d0
      vshap(2) = Xi(1); vgrad(1,2) = 1.d0
      vshap(3) = Xi(2); vgrad(2,3) = 1.d0
      vshap(4) = Xi(3); vgrad(3,4) = 1.d0
c
      k=0
c
c  ...edge shape functions
      do ie=1,6
        nv1 = TETRA_EDGE_TO_VERT(1,ie)
        nv2 = TETRA_EDGE_TO_VERT(2,ie)
        call shape_edgeE_orient((/vshap(nv1),vshap(nv2)/),
     .                          (/vgrad(1:3,nv1),vgrad(1:3,nv2)/),
     .                          3,Norder(ie),Norient_edge(ie),
     .                          nrdof_edgeE,ShapE(1,k+1),CurlE(1,k+1))
        k=k+nrdof_edgeE
      enddo
c
c  ...face shape functions
      do if=1,4
        nv1 = TETRA_FACE_TO_VERT(1,if)
        nv2 = TETRA_FACE_TO_VERT(2,if)
        nv3 = TETRA_FACE_TO_VERT(3,if)
        call shape_triaE_orient((/vshap(nv1),vshap(nv2),vshap(nv3)/),
     .              (/vgrad(1:3,nv1),vgrad(1:3,nv2),vgrad(1:3,nv3)/),
     .                               3,Norder(6+if),Norient_face(if),
     .                        nrdof_faceE,ShapE(1,k+1),CurlE(1,k+1))
        k=k+nrdof_faceE
      enddo
c
c  ...bubbles
      kk=k
      call shape_tetrE(vshap,vgrad,Norder(11),
     .                 NrdofE,ShapE(1,k+1),CurlE(1,k+1))
      NrdofE = kk+NrdofE
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(1:11),
     .                Norient_edge(1:6),Norient_face(1:4)
 7001   format('shape3E_tetr: Xi = ',3f8.3,/,
     .         'Norder = ',6(i2,2x),3x,4(i2,2x),2x,i2/,
     .         'orient = ',6(i2,2x),3x,4(i2,2x))
        k=0
        do ie=1,6
          ndof = Norder(ie)
          if (ndof.gt.0) then
            write(*,7002) ie
 7002       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
 7003         format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
            enddo
          endif
        enddo
        do if=1,4
          nord = Norder(6+if)
          ndof = nord*(nord-1)
          if (ndof.gt.0) then
            write(*,7004) if
 7004       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
            enddo
          endif
        enddo
        nord = Norder(11)
        ndof = (nord-2)*(nord-1)*nord/2
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3E_tetr
c
c--------------------------------------------------------------------
c
c     routine name      - shape3E_pris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D prism element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c          Xi           - master pris coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3E_pris(Xi,Norder,Norient_edge,Norient_face,
     .                        NrdofE,ShapE,CurlE)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(3),Norder(15),Norient_edge(9),Norient_face(5),
     .          ShapE(3,MAXprismE),CurlE(3,MAXprismE)
c
c  ...triangle vertex shape functions, 1D shape functions
      dimension shapHt(3),gradHt(2,3), shapH1(2),dshapH1(2)
c
c  ...1D H1 bubble and L2 shape functions
      dimension shapH1b(MAXP-1),dshapH1b(MAXP-1),shapQ1(MAXP)
c
c  ...work space for shape_egdeE, shape_triaE
      dimension shapEt(2,MAXmdltE),curlEt(MAXmdltE)
c
c  ...work space for shape_triaH_orient, shape_triaH
      dimension shapHtb(MAXmdltH),gradHtb(2,MAXmdltH)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
c  ...triangle vertex shape functions (affine coordinates)
      shapHt(1) = 1.d0-Xi(1)-Xi(2); gradHt(1:2,1)=-1.d0
      shapHt(2) = Xi(1); gradHt(1,2) = 1.d0; gradHt(2,2) = 0.d0
      shapHt(3) = Xi(2); gradHt(1,3) = 0.d0; gradHt(2,3) = 1.d0
c
c  ...1D vertex shape functions (in z)
      shapH1(1) = 1.d0-Xi(3); dshapH1(1) = -1.d0
      shapH1(2) =      Xi(3); dshapH1(2) =  1.d0
c
c
      k=0
c
c  ...horizontal edges
      ie=0
      do j=1,2
        do i=1,3
          ie=ie+1
          nv1 = TRIAN_EDGE_TO_VERT(1,i)
          nv2 = TRIAN_EDGE_TO_VERT(2,i)
          call shape_edgeE_orient((/shapHt(nv1),shapHt(nv2)/),
     .                          (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                          2,Norder(ie),Norient_edge(ie),
     .                          nrdof_edgeE,shapEt,curlEt)
          do ke=1,nrdof_edgeE
            k=k+1
            ShapE(1:2,k) = shapEt(1:2,ke)* shapH1(j)
            ShapE(3,k) = 0.d0
            CurlE(1,k) = - shapEt(2,ke)  *dshapH1(j)
            CurlE(2,k) =   shapEt(1,ke)  *dshapH1(j)
            CurlE(3,k) =   curlEt(ke)    * shapH1(j)
          enddo
        enddo
      enddo
c
c  ...vertical edges
      do i=1,3
        ie=ie+1
        call shape1Q_orient(Xi(3),Norder(ie),Norient_edge(ie),
     .                      nrdof_edgeE,shapQ1)
        do ke=1,nrdof_edgeE
          k=k+1
          ShapE(1:2,k) =   0.d0
          ShapE(3,k)   =   shapHt(i)  *shapQ1(ke)
          CurlE(1,k)   =   gradHt(2,i)*shapQ1(ke)
          CurlE(2,k)   = - gradHt(1,i)*shapQ1(ke)
          CurlE(3,k)   =   0.d0
        enddo
      enddo
c
c  ...horizontal faces
      if=0
      do j=1,2
        if=if+1
        call shape_triaE_orient(shapHt,gradHt,2,Norder(9+if),
     .         Norient_face(if),  nrdof_faceE,shapEt,curlEt)
        do ke=1,nrdof_faceE
          k=k+1
          ShapE(1:2,k) = shapEt(1:2,ke)* shapH1(j)
          ShapE(3,k) = 0.d0
          CurlE(1,k) = - shapEt(2,ke)  *dshapH1(j)
          CurlE(2,k) =   shapEt(1,ke)  *dshapH1(j)
          CurlE(3,k) =   curlEt(ke)    * shapH1(j)
        enddo
      enddo
c
c  ...vertical faces
      do ie=1,3
        if=if+1
        call decode(Norder(9+if), nord1,nord2)
        ndof = nord1*(nord2-1) + (nord1-1)*nord2
        if (ndof.eq.0) cycle
c
        nv1 = TRIAN_EDGE_TO_VERT(1,ie)
        nv2 = TRIAN_EDGE_TO_VERT(2,ie)
c
        iswap    = norient221(1,Norient_face(if))
        norient1 = norient221(2,Norient_face(if))
        norient2 = norient221(3,Norient_face(if))
c
        select case(iswap)
        case(0)
          call shape_edgeH_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                             2,nord1,norient1,
     .                             nrdofHtb,shapHtb,gradHtb)
          call shape_edgeE_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                             2,nord1,norient1,
     .                             nrdofEt,shapEt,curlEt)
          call shape1Hb_orient(Xi(3),nord2,norient2,
     .                         nrdofH1b,shapH1b,dshapH1b)
          call shape1Q_orient(Xi(3),nord2,norient2,
     .                        nrdofE1,shapQ1)
c
c  .......first family
          do j=1,nrdofH1b
            do i=1,nrdofEt
              k=k+1
              ShapE(1:2,k) = shapEt(1:2,i)* shapH1b(j)
              ShapE(3,k) = 0.d0
              CurlE(1,k) = - shapEt(2,i)  *dshapH1b(j)
              CurlE(2,k) =   shapEt(1,i)  *dshapH1b(j)
              CurlE(3,k) =   curlEt(i)    * shapH1b(j)
            enddo
          enddo
c
c  .......second family
          do j=1,nrdofE1
            do i=1,nrdofHtb
              k=k+1
              ShapE(1:2,k) = 0.d0
              ShapE(3,k) =   shapHtb(i)  *shapQ1(j)
              CurlE(1,k) =   gradHtb(2,i)*shapQ1(j)
              CurlE(2,k) = - gradHtb(1,i)*shapQ1(j)
              CurlE(3,k) = 0.d0
            enddo
          enddo
c
c  .....the ordering of face shape functions changes when the face coordinates
c       are swapped
        case(1)
          call shape_edgeH_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                             2,nord2,norient2,
     .                             nrdofHtb,shapHtb,gradHtb)
          call shape_edgeE_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                             2,nord2,norient2,
     .                             nrdofEt,shapEt,curlEt)
          call shape1Hb_orient(Xi(3),nord1,norient1,
     .                         nrdofH1b,shapH1b,dshapH1b)
          call shape1Q_orient(Xi(3),nord1,norient1,
     .                        nrdofE1,shapQ1)
c
c  .......first family
          do i=1,nrdofHtb
            do j=1,nrdofE1
              k=k+1
              ShapE(1:2,k) =  0.d0
              ShapE(3,k) =   shapHtb(i)  *shapQ1(j)
              CurlE(1,k) =   gradHtb(2,i)*shapQ1(j)
              CurlE(2,k) = - gradHtb(1,i)*shapQ1(j)
              CurlE(3,k) = 0.d0
            enddo
          enddo
c
c  .......second family
          do i=1,nrdofEt
            do j=1,nrdofH1b
              k=k+1
              ShapE(1:2,k) = shapEt(1:2,i)* shapH1b(j)
              ShapE(3,k) = 0.d0
              CurlE(1,k) = - shapEt(2,i)  *dshapH1b(j)
              CurlE(2,k) =   shapEt(1,i)  *dshapH1b(j)
              CurlE(3,k) =   curlEt(i)    * shapH1b(j)
            enddo
          enddo
        end select
c
c  ...end of loop through vertical faces
      enddo
c
c  ...bubbles
      call decode(Norder(15), nord1,nord2)
      call shape_triaH(shapHt,gradHt,2,nord1, nrdofHtb,shapHtb,gradHtb)
      call shape_triaE(shapHt,gradHt,2,nord1, nrdofEt,shapEt,curlEt)
      call shape1Hb(Xi(3),nord2, nrdofH1b,shapH1b,dshapH1b)
      call shape1Q(Xi(3),nord2, nrdofE1,shapQ1)
c
c  ...first family
      do j=1,nrdofH1b
        do i=1,nrdofEt
          k=k+1
          ShapE(1:2,k) = shapEt(1:2,i)* shapH1b(j)
          ShapE(3,k) = 0.d0
          CurlE(1,k) = - shapEt(2,i)  *dshapH1b(j)
          CurlE(2,k) =   shapEt(1,i)  *dshapH1b(j)
          CurlE(3,k) =   curlEt(i)    * shapH1b(j)
        enddo
      enddo
c
c  ...second family
      do j=1,nrdofE1
        do i=1,nrdofHtb
          k=k+1
          ShapE(1:2,k) =  0.d0
          ShapE(3,k) =   shapHtb(i)  *shapQ1(j)
          CurlE(1,k) =   gradHtb(2,i)*shapQ1(j)
          CurlE(2,k) = - gradHtb(1,i)*shapQ1(j)
          CurlE(3,k) = 0.d0
        enddo
      enddo
c
      NrdofE=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(1:15),
     .                Norient_edge(1:9),Norient_face(1:5)
 7001   format('shape3E_pris: Xi = ',3f8.3,/,
     .         'Norder = ',3(3i2,2x),2x,2i2,2x,3i3,3x,i3,/,
     .         'orient = ',3(3i2,2x),2x,2i2,2x,3i3)
        k=0
        do ie=1,9
          ndof = Norder(ie)
          if (ndof.eq.0) cycle
          write(*,7002) ie
 7002     format('SHAPE FUNCTIONS FOR EDGE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
 7003       format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
          enddo
        enddo
        do if=1,2
          nord = Norder(9+if)
          ndof = nord*(nord-1)
          if (ndof.eq.0) cycle
          write(*,7004) if
 7004     format('SHAPE FUNCTIONS FOR HORIZONTAL FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        enddo
        do if=3,5
          call decode(Norder(9+if), nordh,nordv)
          ndof = (nordh-1)*nordv + nordh*(nordv-1)
          if (ndof.eq.0) cycle
          write(*,7005) if
 7005     format('SHAPE FUNCTIONS FOR FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        enddo
        call decode(Norder(15), nord1,nord2)
        ndof =  nord1*(nord1-1)*(nord2-1)+(nord1-2)*(nord1-1)/2*nord2
        if (ndof.gt.0) then
          write(*,7006)
 7006     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3E_pris
c
c--------------------------------------------------------------------
c
c     routine name      - shape3E_pyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D pyramid element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c          Xip          - master pyra coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofE       - number of dof
c          ShapE        - values of the shape functions at the point
c          CurlE        - curl of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3E_pyra(Xip,Norder,Norient_edge,Norient_face,
     .                        NrdofE,ShapE,CurlE)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xip(3),Norder(14),Norient_edge(8),Norient_face(5),
     .          ShapE(3,MAXpyramE),CurlE(3,MAXpyramE)
c
c  ...common factor, modified coordinates, vertex shape functions
      dimension psi(4),gradpsi(3,4),alpha(2,4),gradalpha(3,2,4),
     .          vshap(5),gradv(3,5),vec(3),vec1(3)
c
c  ...work space for shape_edgeH_orient
      dimension shapH1b(MAXP-1),gradH1b(MAXP-1),
     .          shapH2b(MAXP-1),gradH2b(MAXP-1)
c
c  ...work space for shape1Q_orient
      dimension shapQ1(MAXP),shapQ2(MAXP)
c
      double precision zeta,z
c
c  ...gradients of infinite element coordinates
      dimension gradx(3),grady(3),gradz(3)
c
c  ...work space for ILegendre and Legendre
      dimension rL1(2:MAXP),dLdx1(2:MAXP),P1(0:MAXP),
     .          rL2(2:MAXP),dLdx2(2:MAXP),P2(0:MAXP),
     .          rL3(1:MAXP+1,1:MAXP),dLdx3(1:MAXP+1,1:MAXP),
ccc P.Gatto, Nov 24, 2014
ccc     .          rL3(1:MAXP+1,1:MAXP-1),dLdx3(1:MAXP+1,1:MAXP-1),
     .          P3(0:MAXP)
c
c
c  ...work space for the pullback
      dimension Piola_curl(3,3),Evec(3),
     .          Piola_div(3,3), Vvec(3)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
      xi = Xip(1); eta = Xip(2); zeta = Xip(3)
      ShapH = 0.d0; GradH = 0.d0
      call blend_pyram(Xip,x,gradx,y,grady,z,gradz,
     .                 psi,gradpsi,alpha,gradalpha,
     .                 vshap,gradv)
c
c  ...pull  back from x,y,zeta to xi,eta,zeta
      Piola_curl(1,1:3) = (/1.d0, 0.d0, 0.d0/)
      Piola_curl(2,1:3) = (/0.d0, 1.d0, 0.d0/)
      Piola_curl(3,1:3) = (/x,y,1.d0-zeta/)
      Piola_curl = Piola_curl/(1.d0-zeta)
c
      Piola_div(1,1:3) = (/1.d0-zeta, 0.d0, -x/)
      Piola_div(2,1:3) = (/0.d0, 1.d0-zeta, -y/)
      Piola_div(3,1:3) = (/0.d0,0.d0,1.d0/)
      Piola_div = Piola_div/(1.d0-zeta)**2
      k=0
c
c  ...horizontal edges
      do ie=1,4
        call shape_edgeE_orient(alpha(1:2,ie),
     .                          gradalpha(1:3,1:2,ie),
     .                          3,Norder(ie),Norient_edge(ie),
     .                          nrdof_edgeE,ShapE(1,k+1),CurlE(1,k+1))
        do j=1,nrdof_edgeE
          k=k+1
          call cross_product(gradpsi(1:3,ie),ShapE(1:3,k), vec)
          CurlE(1:3,k) = psi(ie)*CurlE(1:3,k) + vec(1:3)
          ShapE(1:3,k) = psi(ie)*ShapE(1:3,k)
        enddo
      enddo
c
c  ...vertical edges
      do ie=5,8
        nv1 = ie-4; nv2 = 5
        call shape_edgeE_orient((/vshap(nv1),vshap(nv2)/),
     .                          (/gradv(1:3,nv1),gradv(1:3,nv2)/),
     .                          3,Norder(ie),Norient_edge(ie),
     .                          nrdof_edgeE,ShapE(1,k+1),CurlE(1,k+1))
        k=k+nrdof_edgeE
      enddo
c
c  ...base face...................................................
      if=1
      call decode(Norder(9), nord1,nord2)
      iswap = norient221(1,Norient_face(if))
      norient1 = norient221(2,Norient_face(if))
      norient2 = norient221(3,Norient_face(if))
      select case(iswap)
c
c  ...unswapped coordinates
      case(0)
        xif1 = x; xif2 = y
c
c  ...swapped coordinates
      case(1)
        xif1 = y; xif2 = x
      end select
c
      call shape1Hb_orient(xif1,nord1,norient1,nrdofH1b,shapH1b,gradH1b)
      call shape1Hb_orient(xif2,nord2,norient2,nrdofH2b,shapH2b,gradH2b)
      call shape1Q_orient (xif1,nord1,norient1, nrdofQ1,shapQ1)
      call shape1Q_orient (xif2,nord2,norient2, nrdofQ2,shapQ2)
c
      ksave = k
c
c  ...horizontal (in face owned coordinates) shape functions
      do j=1,nrdofH2b
        do i=1,nrdofQ1
          n = max(i,j+1)+1
          fact = (1.d0-zeta)**n
          dfact = -n*(1.d0-zeta)**(n-1) ! watch for the minus sign
          k=k+1
          select case(iswap)
          case(0)
            ShapE(1,k) = shapQ1(i)*shapH2b(j)*fact
            ShapE(2:3,k) = 0.d0
            CurlE(1,k) = 0.d0
            CurlE(2,k) =  shapQ1(i)*shapH2b(j)*dfact
            CurlE(3,k) = -shapQ1(i)*gradH2b(j)* fact
          case(1)
            ShapE(1,k) = 0.d0
            ShapE(2,k) = shapQ1(i)*shapH2b(j)* fact
            ShapE(3,k) = 0.d0
            CurlE(1,k) = -shapQ1(i)*shapH2b(j)*dfact
            CurlE(2,k) =  0.d0
            CurlE(3,k) =  shapQ1(i)*gradH2b(j)* fact
          end select
        enddo
      enddo
c
c  ...vertical (in face owned coordinates) shape functions
      do j=1,nrdofQ2
        do i=1,nrdofH1b
          n = max(i+1,j)+1
          fact = (1.d0-zeta)**n
          dfact = -n*(1.d0-zeta)**(n-1)
          k=k+1
          select case(iswap)
          case(0)
            ShapE(1,k) = 0.d0
            ShapE(2,k) = shapH1b(i)*shapQ2(j)*fact
            ShapE(3,k) = 0.d0
            CurlE(1,k) = -shapH1b(i)*shapQ2(j)*dfact
            CurlE(2,k) =  0.d0
            CurlE(3,k) =  gradH1b(i)*shapQ2(j)* fact
          case(1)
            ShapE(1,k) = shapH1b(i)*shapQ2(j)*fact
            ShapE(2,k) = 0.d0
            ShapE(3,k) = 0.d0
            CurlE(1,k) =  0.d0
            CurlE(2,k) =  shapH1b(i)*shapQ2(j)*dfact
            CurlE(3,k) = -gradH1b(i)*shapQ2(j)* fact
          end select
        enddo
      enddo
      NrdofE = k
c
c  ...pull back base shape functions from x,y,zeta to xi,eta,zeta
      do k1=ksave+1,k
        Evec(1:3) = ShapE(1:3,k1)
        ShapE(1:3,k1) = 0.d0
        do l=1,3
          ShapE(1:3,k1) = ShapE(1:3,k1) + Piola_curl(1:3,l)*Evec(l)
        enddo
        Vvec(1:3) = CurlE(1:3,k1)
        CurlE(1:3,k1) = 0.d0
        do l=1,3
          CurlE(1:3,k1) = CurlE(1:3,k1) + Piola_div(1:3,l)*Vvec(l)
        enddo
      enddo
c
c  ...vertical faces
      do ie=1,4
        if=ie+1
        call shape_triaE_orient((/alpha(1,ie),alpha(2,ie),vshap(5)/),
     .      (/gradalpha(1:3,1,ie),gradalpha(1:3,2,ie),gradv(1:3,5)/),
     .                               3,Norder(8+if),Norient_face(if),
     .                         nrdof_faceE,ShapE(1,k+1),CurlE(1,k+1))
c
c  .....multiply by the common factor blending function
        do j=1,nrdof_faceE
          k=k+1
          call cross_product(gradpsi(1:3,ie),ShapE(1:3,k), vec)
          CurlE(1:3,k) = CurlE(1:3,k)*psi(ie) + vec(1:3)
          ShapE(1:3,k) = ShapE(1:3,k)*psi(ie)
        enddo
      enddo
c
c  ...bubbles...........................................................
      ksave = k
      nord = Norder(14)
      call ILegendre(   x,1.d0,nord,0, rL1,dLdx1,void)
      call ILegendre(   y,1.d0,nord,0, rL2,dLdx2,void)
      call IJacobi(zeta,1.d0,nord+1,nord,0, rL3,dLdx3,void)
c
      call LegendreNEW(   x,1.d0,nord, P1)
      call LegendreNEW(   y,1.d0,nord, P2)
      call LegendreNEW(zeta,1.d0,nord, P3)
c
      if (iprint.eq.2) then
        rL1(2)=x*(1.d0-x); dLdx1(2)=1.d0-2.d0*x
        rL2(2)=y*(1.d0-y); dLdx2(2)=1.d0-2.d0*y
        rL3(1:3,1)=zeta; dLdx3(1:3,1) = 1.d0
        P1(0) = 1.d0; P1(1) = 1.d0-2.d0*x
        P2(0) = 1.d0; P2(1) = 1.d0-2.d0*y
        P3(0) = 1.d0; P3(1) = zeta-0.5d0
      endif
c
c  ...first family
      if (iprint.eq.2) write(*,*) 'Family I'
      do m=1,nord-1
        do j=2,nord
          do i=0,nord-1
            ij1 = max(i+1,j)+1
            n = max(i+1,j,m+1)+1
            k=k+1
            if (iprint.eq.2) write(*,*) 'i,j,m,n-m = ',i,j,m,n-m
            ShapE(1,k) = P1(i)*rL2(j)*rL3(ij1,m)*(1-zeta)**(n-m)
            ShapE(2:3,k) = 0.d0
            CurlE(1,k) = 0.d0
            CurlE(2,k) = P1(i)*rL2(j)*
     .                  (dLdx3(ij1,m)*(1-zeta)**(n-m)
     .                  -(n-m)*rL3(ij1,m)*(1-zeta)**(n-m-1))
            CurlE(3,k) = -P1(i)*dLdx2(j)*rL3(ij1,m)*(1-zeta)**(n-m)
          enddo
        enddo
      enddo
      if (iprint.eq.2) then
        write(*,*) 'ShapE(1,k) = ',ShapE(1,k),
     .             (1.d0-2.d0*x)*y*(1.d0-y)*zeta*(1.d0-zeta)**2

        write(*,*) 'CurlE(1:3,k) = ',CurlE(1:3,k),
     .             0.d0,
     .           (1.d0-2.d0*x)*y*(1.d0-y)*(1.d0-4.d0*zeta+3.d0*zeta**2),
     .           -(1.d0-2.d0*x)*(1.d0-2.d0*y)*zeta*(1.d0-zeta)**2
      endif
c
c  ...second family
      if (iprint.eq.2) write(*,*) 'Family II'
      do m=1,nord-1
        do j=0,nord-1
          do i=2,nord
            ij1 = max(i,j+1)+1
            n = max(i,j+1,m+1)+1
            k=k+1
            if (iprint.eq.2) write(*,*) 'i,j,m,n-m = ',i,j,m,n-m
            ShapE(1,k) = 0.d0
            ShapE(2,k) = rL1(i)*P2(j)*rL3(ij1,m)*(1-zeta)**(n-m)
            ShapE(3,k) = 0.d0
            CurlE(1,k) = - rL1(i)*P2(j)*
     .                  ( dLdx3(ij1,m)*(1-zeta)**(n-m)
     .                   -(n-m)*rL3(ij1,m)*(1-zeta)**(n-m-1))
            CurlE(2,k) = 0.d0
            CurlE(3,k) = dLdx1(i)*P2(j)*rL3(ij1,m)*(1-zeta)**(n-m)
          enddo
        enddo
      enddo
      if (iprint.eq.2) then
        write(*,*) 'ShapE(2,k) = ',ShapE(2,k),
     .             x*(1.d0-x)*(1.d0-2.d0*y)*zeta*(1.d0-zeta)**2
      endif
c
c  ...third family
      if (iprint.eq.2) write(*,*) 'Family III'
      do m=0,nord-2
        do j=2,nord
          do i=2,nord
            n = max(i,j,m+2)-1
            k=k+1
            if (iprint.eq.2) write(*,*) 'i,j,m,n-m = ',i,j,m,n-m
            ShapE(1:2,k) = 0.d0
            ShapE(3,k) = rL1(i)*rL2(j)*P3(m)*(1-zeta)**(n-m)
c
            CurlE(1,k) = rL1(i)*dLdx2(j)*P3(m)*(1-zeta)**(n-m)
            CurlE(2,k) = -dLdx1(i)*rL2(j)*P3(m)*(1-zeta)**(n-m)
            CurlE(3,k) = 0.d0
          enddo
        enddo
      enddo
      if (iprint.eq.2) then
        write(*,*) 'ShapE(3,k) = ',ShapE(3,k),
     .             x*(1.d0-x)*y*(1.d0-y)*(1.d0-zeta)
        write(*,*) 'CurlE(1:3,k) = ',CurlE(1:3,k),
     .             x*(1.d0-x)*(1.d0-2.d0*y)*(1.d0-zeta),
     .             -(1.d0-2.d0*x)*y*(1.d0-y)*(1.d0-zeta), 0.d0
      endif
c
c  ...fourth family
      if (iprint.eq.2) write(*,*) 'Family IV'
      do j=2,nord
        do i=2,nord
          n = max(i,j)-1
          k=k+1
          if (iprint.eq.2) write(*,*) 'i,j,n = ',i,j,n
          aux = zeta*(1.d0-zeta)
          ShapE(1,k) = dLdx1(i)*rL2(j)*aux
          ShapE(2,k) = rL1(i)*dLdx2(j)*aux
          ShapE(3,k) = -rL1(i)*rL2(j)
c
c  .......curl so far
          CurlE(1,k) = -rL1(i)*dLdx2(j)*2.d0*(1.d0-zeta)
          CurlE(2,k) =  dLdx1(i)*rL2(j)*2.d0*(1.d0-zeta)
          CurlE(3,k) = 0.d0
c
c  .......multiply by factor zeta**n
          fact = zeta**n
          vec(1:3) = (/0.d0,0.d0,n*zeta**(n-1)/)
          call cross_product(vec,ShapE(1:3,k), vec1)
          ShapE(1:3,k) = fact*ShapE(1:3,k)
          CurlE(1:3,k) = fact*CurlE(1:3,k) + vec1(1:3)
        enddo
      enddo
      if (iprint.eq.2) then
        r = x*(1.d0-x)*y*(1.d0-y)
        drdx = (1.d0-2.d0*x)*y*(1.d0-y)
        drdy = x*(1.d0-x)*(1.d0-2.d0*y)
        write(*,*) 'ShapE(1:3,k) = ',ShapE(1:3,k)
        write(*,*) drdx*zeta**2*(1.d0-zeta),
     .             drdy*zeta**2*(1.d0-zeta), - r*zeta
        write(*,*) 'CurlE(1:3,k) = ',CurlE(1:3,k),
     .             -drdy*(zeta+2.d0*zeta-3.d0*zeta**2),
     .              drdx*(zeta+2.d0*zeta-3.d0*zeta**2),0.d0
      endif

c
      NrdofE = k
c
c  ...pull bubbles back
      do k=ksave+1,NrdofE
        Evec(1:3) = ShapE(1:3,k)
        ShapE(1:3,k) = 0.d0
        do l=1,3
          ShapE(1:3,k) = ShapE(1:3,k) + Piola_curl(1:3,l)*Evec(l)
        enddo
        Vvec(1:3) = CurlE(1:3,k)
        CurlE(1:3,k) = 0.d0
        do l=1,3
          CurlE(1:3,k) = CurlE(1:3,k) + Piola_div(1:3,l)*Vvec(l)
        enddo
      enddo
c
c
      if (iprint.ge.1) then
        write(*,7001) Xip(1:3),Norder(1:14),
     .                Norient_edge(1:8),Norient_face(1:5)
 7001   format('shape3E_pyra: Xip = ',3f8.3,/,
     .         'Norder = ',8(i2,2x),3x,i3,2x,4(i2,2x),2x,i2/,
     .         'orient = ',8(i2,2x),3x,i3,2x,4(i2,2x))
        k=0
        do ie=1,8
          ndof = Norder(ie)
          if (ndof.eq.0) cycle
          write(*,7002) ie
 7002     format('SHAPE FUNCTIONS FOR EDGE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
 7003       format('k = ',i3,' ShapE, CurlE = ',3e12.5,3x,3e12.5)
          enddo
        enddo
        call decode(Norder(9), nord1,nord2)
        ndof = nord1*(nord2-1) + (nord1-1)*nord2
        if (ndof.gt.0) then
          write(*,7004)
 7004     format('SHAPE FUNCTIONS FOR BASE FACE')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        endif
        do if=2,5
          nord = Norder(8+if)
          ndof = nord*(nord-1)
          if (ndof.eq.0) cycle
          write(*,7005) if
 7005     format('SHAPE FUNCTIONS FOR FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        enddo
        nord = Norder(14)
        ndof = 3*(nord-1)**2*nord
        if (ndof.gt.0) then
          write(*,7006)
 7006     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapE(1:3,k),CurlE(1:3,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3E_pyra




