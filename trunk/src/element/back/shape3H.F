c--------------------------------------------------------------------
c
c     routine name      - shape3H
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of 3D element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Type         - element type
c          Xi           - master element coordinates
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3H_1(Type,Xi,Norder,Norient_edge,Norient_face,
     .                     NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(3),Norder(19),Norient_edge(12),Norient_face(6),
     .          ShapH(MAXbrickH),GradH(3,MAXbrickH)
c

      select case(Type)
      case('bric','mdlb')
        call shape3H_bric(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofH,ShapH,GradH)
      case('tetr','mdln')
        call shape3H_tetr(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofH,ShapH,GradH)
      case('pris','mdlp')
        call shape3H_pris(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofH,ShapH,GradH)
      case('pyra','mdld')
        call shape3H_pyra(Xi,Norder,Norient_edge,Norient_face,
     .                    NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape3H: Type = ', Type
        stop
      end select
c
      end subroutine shape3H_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape3H_bric
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 3D brick element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Xi           - master brick coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH       - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3H_bric(Xi,Norder,Norient_edge,Norient_face,
     .                        NrdofH,ShapH,GradH)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),Norder(19),Norient_edge(12),Norient_face(6),
     .          ShapH(MAXbrickH),GradH(3,MaxbrickH)
c
c  ...blending functions and projections
      dimension psie(12),gradpsie(3,12),te(12),gradte(3,12),
     .          psif(6),gradpsif(3,6),tf(2,6),gradtf(2,3,6)
c
c  ...work space for face bubbles
      dimension dshapHf(2)
c
c  ...work space for shape1H, shape1Hb_orient, shape1Hb
      dimension shapH1(MAXP-1),dshapH1(MAXP-1),
     .          shapH2(MAXP-1),dshapH2(MAXP-1),
     .          shapH3(MAXP-1),dshapH3(MAXP-1)
c
c  ...translating quad orientation into product of 1D orientation
      dimension norient221(1:3,0:7)
c
c  ...definition of 3D vertex shape function in terms of 1D shape
c     functions
      dimension ivert(3,8)
      data ivert /1,1,1, 2,1,1, 2,2,1, 1,2,1,
     .            1,1,2, 2,1,2, 2,2,2, 1,2,2/
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
c  ...initiate counter for shape functions
      k=0
c
c  ...vertex shape functions
      call shape1H(Xi(1),1, nrdofH1,shapH1,dshapH1)
      call shape1H(Xi(2),1, nrdofH2,shapH2,dshapH2)
      call shape1H(Xi(3),1, nrdofH3,shapH3,dshapH3)
      do iv=1,8
        i = ivert(1,iv); j = ivert(2,iv); l = ivert(3,iv)
        k=k+1
        ShapH(k) = shapH1(i)*shapH2(j)*shapH3(l)
        GradH(1,k) = dshapH1(i)* shapH2(j)* shapH3(l)
        GradH(2,k) =  shapH1(i)*dshapH2(j)* shapH3(l)
        GradH(3,k) =  shapH1(i)* shapH2(j)*dshapH3(l)
      enddo
c
c  ...determine edge blending functions and projections
      call blend_project_edges_bric(Xi, psie,gradpsie,te,gradte)
c
c  ...edge shape functions
      do ie=1,12
        if (Norder(ie).gt.1) then
          call shape1Hb_orient(te(ie),Norder(ie),Norient_edge(ie),
     .                         nrdof_edgeH,shapH1,dshapH1)
          do ke=1,nrdof_edgeH
            k=k+1
            ShapH(k) = psie(ie)*shapH1(ke)
            GradH(1:3,k) = gradpsie(1:3,ie)*shapH1(ke)
     .                    + psie(ie)*dshapH1(ke)*gradte(1:3,ie)
          enddo
        endif
      enddo
c
c  ...determine face blending functions and projections
      call blend_project_faces_bric(Xi, psif,gradpsif,tf,gradtf)
c
c  ...face shape functions
      do if=1,6
        call decode(Norder(12+if), nord1,nord2)
        ndof = (nord1-1)*(nord2-1)
        if (ndof.eq.0) cycle
c
        iswap = norient221(1,Norient_face(if))
        norient1 = norient221(2,Norient_face(if))
        norient2 = norient221(3,Norient_face(if))
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
        call shape1Hb_orient(xif1,nord1,norient1,
     .                       nrdofH1,shapH1,dshapH1)
        call shape1Hb_orient(xif2,nord2,norient2,
     .                       nrdofH2,shapH2,dshapH2)
c
        do j=1,nrdofH2
          do i=1,nrdofH1
            k=k+1
            ShapH(k) = psif(if)*shapH1(i)*shapH2(j)
            select case(iswap)
            case(0)
              dshapHf(1) = dshapH1(i)* shapH2(j)
              dshapHf(2) =  shapH1(i)*dshapH2(j)
            case(1)
              dshapHf(1) =  shapH1(i)*dshapH2(j)
              dshapHf(2) = dshapH1(i)* shapH2(j)
            end select
            GradH(1:3,k) = gradpsif(1:3,if)*shapH1(i)*shapH2(j)
     .                    + psif(if)*(dshapHf(1)*gradtf(1,1:3,if)
     .                               +dshapHf(2)*gradtf(2,1:3,if))
          enddo
        enddo
      enddo
c
c  ...bubbles = tensor products of 1D bubles
      call decode(Norder(19), norda,nord3)
      call decode(norda, nord1,nord2)
      call shape1Hb(Xi(1),nord1, nrdofH1,shapH1,dshapH1)
      call shape1Hb(Xi(2),nord2, nrdofH2,shapH2,dshapH2)
      call shape1Hb(Xi(3),nord3, nrdofH3,shapH3,dshapH3)
c
      do l=1,nrdofH3
        do j=1,nrdofH2
          do i=1,nrdofH1
            k=k+1
            ShapH(k) = shapH1(i)*shapH2(j)*shapH3(l)
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
        write(*,7001) Xi(1:3),Norder(1:19),
     .                Norient_edge(1:12),Norient_face(1:6)
 7001   format('shape2H_quad: Xi = ',3f8.3,/,
     .         'Norder = ',3(4i2,2x),2i3,2x,4i3,3x,i4,/,
     .         'orient = ',3(4i2,2x),2i3,2x,4i3)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do iv=1,8
          k=iv
          write(*,7002) k,ShapH(k),GradH(1:3,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do ie=1,12
          ndof = Norder(ie)-1
          if (ndof.gt.0) then
            write(*,7011) ie
 7011       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapH(k),GradH(1:3,k)
            enddo
          endif
        enddo
        do if=1,6
          call decode(Norder(12+if), nordh,nordv)
          ndof = (nordh-1)*(nordv-1)
          if (ndof.gt.0) then
            write(*,7012) if
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapH(k),GradH(1:3,k)
            enddo
          endif
        enddo
        call decode(Norder(19), nord,nord3)
        call decode(nord, nord1,nord2)
        ndof = (nord1-1)*(nord2-1)*(nord3-1)
        if (ndof.gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3H_bric
c
c--------------------------------------------------------------------
c
c     routine name      - shape3H_tetr
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D tetra element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Xi           - master tetrk coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3H_tetr(Xi,Norder,Norient_edge,Norient_face,
     .                        NrdofH,ShapH,GradH)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(3),Norder(11),Norient_edge(6),Norient_face(4),
     .          ShapH(MAXtetraH),GradH(3,MAXtetraH)
c
      dimension dt(3)
c
c  ...work space for shape_edgeH_orient
      dimension shap_edgeH(MAXP-1),dshap_edgeH(MAXP-1),
     .                             dshap_edgeHdt(MAXP-1)
c
c  ...work space for shape_triaH_orient
      dimension tf(2),shap_faceH((MAXP-2)*(MAXP-1)/2),
     .          dshap_faceH(2,(MAXP-2)*(MAXP-1)/2),
     .          dshap_faceHdt((MAXP-2)*(MAXP-1)/2)
c
c
      iprint=0
c
c  ...vertex shape functions
      ShapH(1) = 1.d0 - Xi(1) - Xi(2) - Xi(3)
      GradH(1:3,1) = -1.d0
      GradH(1:3,2:4) = 0.d0
      ShapH(2) = Xi(1); GradH(1,2) = 1.d0
      ShapH(3) = Xi(2); GradH(2,3) = 1.d0
      ShapH(4) = Xi(3); GradH(3,4) = 1.d0
      k=4
c
c  ...edge shape functions
      do ie=1,6
        if (Norder(ie).le.1) cycle
        nv1 = TETRA_EDGE_TO_VERT(1,ie)
        nv2 = TETRA_EDGE_TO_VERT(2,ie)
        call shape_edgeH_orient((/ShapH(nv1),ShapH(nv2)/),
     .                          (/GradH(1:3,nv1),GradH(1:3,nv2)/),
     .                           3,Norder(ie),Norient_edge(ie),
     .                           nrdof_edgeH,ShapH(k+1),GradH(1:3,k+1))
       k=k+nrdof_edgeH
      enddo
c
c  ...face shape functions
      do if=1,4
        if (Norder(6+if).le.2) cycle
        nv1 = TETRA_FACE_TO_VERT(1,if)
        nv2 = TETRA_FACE_TO_VERT(2,if)
        nv3 = TETRA_FACE_TO_VERT(3,if)
        call shape_triaH_orient((/ShapH(nv1),ShapH(nv2),ShapH(nv3)/),
     .              (/GradH(1:3,nv1),GradH(1:3,nv2),GradH(1:3,nv3)/),
     .                               3,Norder(6+if),Norient_face(if),
     .                         nrdof_faceH,ShapH(k+1),GradH(1:3,k+1))
        k=k+nrdof_faceH
      enddo
c
c  ...bubbles
      call shape_tetrH(Xi,Norder(11), NrdofH,ShapH(k+1),GradH(1,k+1))
      if (iprint.eq.1) write(*,*) 'shape3H_tetr: k,NrdofH =',k,NrdofH
      NrdofH = k+NrdofH
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(1:11),
     .                Norient_edge(1:6),Norient_face(1:4)
 7001   format('shape3H_tetr: Xi = ',3f8.3,/,
     .         'Norder = ',6i2,2x,4i2,2x,i2,/,
     .         'orient = ',6i2,2x,4i2)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do iv=1,4
          k=iv
          write(*,7002) k,ShapH(k),GradH(1:3,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do ie=1,6
          ndof = Norder(ie)-1
          if (ndof.gt.0) then
            write(*,7011) ie
 7011       format('SHAPE FUNCTIONS FOR EDGE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapH(k),GradH(1:3,k)
            enddo
          endif
        enddo
        do if=1,4
          nord = Norder(6+if)
          ndof = (nord-1)*(nord-2)/2
          if (ndof.gt.0) then
            write(*,7012) if
 7012       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapH(k),GradH(1:3,k)
            enddo
          endif
        enddo
        nord = Norder(11)
        ndof = (nord-1)*(nord-2)*(nord-3)/6
        if (ndof.gt.0) then
          write(*,7013)
 7013     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        endif
ccc        call pause
      endif
c
c
      end subroutine shape3H_tetr
c
c--------------------------------------------------------------------
c
c     routine name      - shape3H_pris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 3D prism element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Xi           - master prisk coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH       - derivatives of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3H_pris(Xi,Norder,Norient_edge,Norient_face,
     .                        NrdofH,ShapH,GradH)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(3),Norder(15),Norient_edge(9),Norient_face(5),
     .          ShapH(MAXprismH),GradH(3,MAXprismH)
c
c  ...triangle vertex shape functions, 1D shape functions
      dimension shapHt(3),gradHt(2,3),shapH1(2),dshapH1(2)
c
c  ...work space for shape_triaH, shape1Hb_orient,shape1Hb
      dimension shapHtb(MAXmdltH),gradHtb(2,MAXmdltH),
     .          shapH1b(MAXP-1),dshapH1b(MAXP-1)
c
c  ...translating quad orientation into product of 1D orientation
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
c  ...initiate counter for shape functions
      k=0
c
c  ...vertex shape functions = tensor products of triangle shape
c     functions and 1D shape functions in z direction
      shapHt(1) = 1.d0-Xi(1)-Xi(2); gradHt(1:2,1) = -1.d0
      shapHt(2) = Xi(1); gradHt(1,2) = 1.d0; gradHt(2,2) = 0.d0
      shapHt(3) = Xi(2); gradHt(1,3) = 0.d0; gradHt(2,3) = 1.d0
      shapH1(1) = 1.d0-Xi(3); dshapH1(1) = -1.d0
      shapH1(2) =      Xi(3); dshapH1(2) =  1.d0
      do j=1,2
        do i=1,3
          k=k+1
          ShapH(k) = shapHt(i)*shapH1(j)
          GradH(1:2,k) = gradHt(1:2,i)*shapH1(j)
          GradH(3,k) = shapHt(i)*dshapH1(j)
        enddo
      enddo
c
c  ...horizontal edges
      ie=0
      do j=1,2
        do i=1,3
          ie=ie+1
          if (Norder(ie).eq.1) cycle
c
          nv1 = TRIAN_EDGE_TO_VERT(1,i)
          nv2 = TRIAN_EDGE_TO_VERT(2,i)
          call shape_edgeH_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                            2,Norder(ie),Norient_edge(ie),
     .                            nrdof_edgeH,shapHtb,gradHtb)
          do ke=1,nrdof_edgeH
            k=k+1
            ShapH(k) = shapHtb(ke)*shapH1(j)
            GradH(1:2,k) = gradHtb(1:2,ke)*shapH1(j)
            GradH(3,k) = shapHtb(ke)*dshapH1(j)
          enddo
        enddo
      enddo
c
c  ...vertical edges
      do i=1,3
        ie=ie+1
        if (Norder(ie).eq.1) cycle
c
        call shape1Hb_orient(Xi(3),Norder(ie),Norient_edge(ie),
     .                       nrdof_edgeH,shapH1b,dshapH1b)
        do ke=1,nrdof_edgeH
          k=k+1
          ShapH(k) = shapHt(i)*shapH1b(ke)
          GradH(1:2,k) = gradHt(1:2,i)*shapH1b(ke)
          GradH(3,k) = shapHt(i)*dshapH1b(ke)
        enddo
      enddo
c
c  ...horizontal faces
      if=0
      do j=1,2
        if=if+1
        if (Norder(9+if).lt.3) cycle
        call shape_triaH_orient(shapHt,gradHt,2,Norder(9+if),
     .                          Norient_face(if),
     .                          nrdof_faceH,shapHtb,gradHtb)
        do kf=1,nrdof_faceH
          k=k+1
          ShapH(k) = shapHtb(kf)*shapH1(j)
          GradH(1:2,k) = gradHtb(1:2,kf)*shapH1(j)
          GradH(3,k) = shapHtb(kf)*dshapH1(j)
        enddo
      enddo
c
c  ...vertical faces
      do ie=1,3
        if=if+1
c
c  .....nord1,nord2 refer to horizontal and vertical order in element
c       coordinates...
        call decode(Norder(9+if), nord1,nord2)
        ndof = (nord1-1)*(nord2-1)
        if (ndof.eq.0) cycle
c
        nv1 = TRIAN_EDGE_TO_VERT(1,ie)
        nv2 = TRIAN_EDGE_TO_VERT(2,ie)
c
        iswap = norient221(1,Norient_face(if))
        norient1 = norient221(2,Norient_face(if))
        norient2 = norient221(3,Norient_face(if))
        call shape_edgeH_orient((/shapHt(nv1),shapHt(nv2)/),
     .                          (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                           2,nord1,norient1,
     .                           nrdofHtb,shapHtb,gradHtb)
        call shape1Hb_orient(Xi(3),nord2,norient2,
     .                       nrdofH1b,shapH1b,dshapH1b)
c
c  .....the ordering of face shape functions changes when the face coordinates
c       are swapped
        select case(iswap)
c
c  .....unswapped coordinates
        case(0)
          do j=1,nrdofH1b
            do i=1,nrdofHtb
              k=k+1
              ShapH(k) = shapHtb(i)*shapH1b(j)
              GradH(1:2,k) = gradHtb(1:2,i)*shapH1b(j)
              GradH(3,k) = shapHtb(i)*dshapH1b(j)
            enddo
          enddo
        case(1)
          do i=1,nrdofHtb
            do j=1,nrdofH1b
              k=k+1
              ShapH(k) = shapHtb(i)*shapH1b(j)
              GradH(1:2,k) = gradHtb(1:2,i)*shapH1b(j)
              GradH(3,k) = shapHtb(i)*dshapH1b(j)
            enddo
          enddo
        end select
      enddo
c
c  ...bubbles = tensor products of triangle and 1D bubbles
      call decode(Norder(15), nord1,nord2)
      call shape_triaH(shapHt,gradHt,2,nord1, nrdofHtb,shapHtb,gradHtb)
      call shape1Hb(Xi(3),nord2, nrdofH1b,shapH1b,dshapH1b)
      do j=1,nrdofH1b
        do i=1,nrdofHtb
          k=k+1
          ShapH(k) = shapHtb(i)*shapH1b(j)
          GradH(1:2,k) = gradHtb(1:2,i)*shapH1b(j)
          GradH(3,k) = shapHtb(i)*dshapH1b(j)
        enddo
      enddo
c
      NrdofH = k
c
c
      if (iprint.ge.1) then
        write(*,7001) Xi(1:3),Norder(1:15),
     .                Norient_edge(1:9),Norient_face(1:5)
 7001   format('shape3H_pris: Xi = ',3f8.3,/,
     .         'Norder = ',3(3i2,2x),2x,2i2,2x,3i3,3x,i3,/,
     .         'orient = ',3(3i2,2x),2x,2i2,2x,3i3)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do iv=1,6
          k=iv
          write(*,7002) k,ShapH(k),GradH(1:3,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do ie=1,9
          ndof = Norder(ie)-1
          if (ndof.eq.0) cycle
          write(*,7011) ie
 7011     format('SHAPE FUNCTIONS FOR EDGE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        enddo
        do if=1,2
          nord = Norder(9+if)
          ndof = (nord-2)*(nord-1)/2
          if (ndof.eq.0) cycle
          write(*,7012) if
 7012     format('SHAPE FUNCTIONS FOR HORIZONTAL FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        enddo
        do if=3,5
          call decode(Norder(9+if), nordh,nordv)
          ndof = (nordh-1)*(nordv-1)
          if (ndof.eq.0) cycle
          write(*,7013) if-2
 7013     format('SHAPE FUNCTIONS FOR VERTICAL FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        enddo
        call decode(Norder(15), nord1,nord2)
        ndof = (nord1-2)*(nord1-1)/2*(nord2-1)
        if (ndof.gt.0) then
          write(*,7014)
 7014     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3H_pris
c
c--------------------------------------------------------------------
c
c     routine name      - shape3H_pyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D pyramid element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c          Xip          - master pyramid coordinates
c          Norder       - polynomial order
c          Norient_edge - edge orientation
c          Norient_face - face orientation
c
c     out:
c          NrdofH       - number of dof
c          ShapH        - values of the shape functions at the point
c          GradH        - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3H_pyra(Xip,Norder,Norient_edge,Norient_face,
     .                        NrdofH,ShapH,GradH)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xip(3),Norder(14),Norient_edge(8),Norient_face(5),
     .          ShapH(MAXpyramH),GradH(3,MAXpyramH)
c
c  ...common factor, modified coordinates
      dimension psi(4),gradpsi(3,4),alpha(2,4),gradalpha(3,2,4)
c
c  ...work space for shape1Hb_orient
      dimension shapH1b(MAXP-1),gradH1b(MAXP-1),
     .          shapH2b(MAXP-1),gradH2b(MAXP-1)
c
c  ...work space for ILegendre and IJacobi
      dimension rL1(2:MAXP),dLdx1(2:MAXP),
     .          rL2(2:MAXP),dLdx2(2:MAXP),
     .          rL3(1:MAXP+1,1:MAXP),dLdx3(1:MAXP+1,1:MAXP)
ccc P.Gatto, Nov 24, 2014
ccc     .          rL3(1:MAXP+1,1:MAXP-1),dLdx3(1:MAXP+1,1:MAXP-1)
c
c  ...work space for the pullback
      dimension Piola_curl(3,3),grad(3)
c
      double precision zeta,z
c
c  ...gradients of infinite element coordinates
      dimension gradx(3),grady(3),gradz(3)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      real*8,parameter :: eps = 1.d-12
c
      iprint=0
c
      xi = Xip(1); eta = Xip(2); zeta = Xip(3)
c
c  ...avoid singularity at zeta = 1
      if (abs(zeta - 1.d0) < eps)  zeta=1.d0-eps
c
      call blend_pyram(Xip,x,gradx,y,grady,z,gradz,
     .                 psi,gradpsi,alpha,gradalpha,
     .                 ShapH,GradH)
c
c  ...pullback from x,y,zeta to xi,eta,zeta for gradients
      Piola_curl(1,1:3) = (/1.d0, 0.d0, 0.d0/)
      Piola_curl(2,1:3) = (/0.d0, 1.d0, 0.d0/)
      Piola_curl(3,1:3) = (/x,y,1.d0-zeta/)
      Piola_curl = Piola_curl/(1.d0-zeta)
c
      k=5
c
c  ...horizontal edges shape functions
      do ie=1,4
        if (Norder(ie).le.1) cycle
        call shape_edgeH_orient(alpha(1:2,ie),
     .                          gradalpha(1:3,1:2,ie),
     .                          3,Norder(ie),Norient_edge(ie),
     .                          nrdof_edgeH,ShapH(k+1),GradH(1:3,k+1))
c
c  .....multiply by the common factor blending function
        do j=1,nrdof_edgeH
          k=k+1
          GradH(1:3,k) = GradH(1:3,k)*psi(ie) + ShapH(k)*gradpsi(1:3,ie)
          ShapH(k) = ShapH(k)*psi(ie)
        enddo
      enddo
c
c  ...vertical edges shape functions
      do ie=5,8
        if (Norder(ie).le.1) cycle
        nv1 = PYRAM_EDGE_TO_VERT(1,ie)
        nv2 = PYRAM_EDGE_TO_VERT(2,ie)
        call shape_edgeH_orient((/ShapH(nv1),ShapH(nv2)/),
     .                          (/GradH(1:3,nv1),GradH(1:3,nv2)/),
     .                           3,Norder(ie),Norient_edge(ie),
     .                           nrdof_edgeH,ShapH(k+1),GradH(1:3,k+1))
        k=k+nrdof_edgeH
      enddo
c
c  ...horizontal face, use x,y,zeta frame
      call decode(Norder(9), nord1,nord2)
      iswap = norient221(1,Norient_face(1))
      norient1 = norient221(2,Norient_face(1))
      norient2 = norient221(3,Norient_face(1))
      call shape1Hb_orient(x,nord1,norient1, nrdofH1b,shapH1b,gradH1b)
      call shape1Hb_orient(y,nord2,norient2, nrdofH2b,shapH2b,gradH2b)
      ksave = k
c
c  ...the ordering of face shape functions changes when the face coordinates
c     are swapped
      select case(iswap)
c
c  ...unswapped coordinates
      case(0)
        do j=1,nrdofH2b
          do i=1,nrdofH1b
            n = max(i,j)+1 ! polynomial order in x,y is i+1,j+1
            fact = (1.d0-zeta)**n
            dfact = -n*(1.d0-zeta)**(n-1) ! watch for the minus sign
            k=k+1
            ShapH(k)   = shapH1b(i)*shapH2b(j)* fact
            GradH(1,k) = gradH1b(i)*shapH2b(j)* fact
            GradH(2,k) = shapH1b(i)*gradH2b(j)* fact
            GradH(3,k) = shapH1b(i)*shapH2b(j)*dfact
          enddo
        enddo
      case(1)
        do i=1,nrdofH1b
          do j=1,nrdofH2b
            n = max(i,j)+1
            fact = (1.d0-zeta)**n
            dfact = -n*(1.d0-zeta)**(n-1)
            k=k+1
            ShapH(k)   = shapH1b(i)*shapH2b(j)* fact
            GradH(1,k) = gradH1b(i)*shapH2b(j)* fact
            GradH(2,k) = shapH1b(i)*gradH2b(j)* fact
            GradH(3,k) = shapH1b(i)*shapH2b(j)*dfact
          enddo
        enddo
      end select
c
c  ...pullback the base shape functions gradients from x,y,zeta to xi,eta,zeta
      do k1=ksave+1,k
        grad(1:3) = GradH(1:3,k1)
        GradH(1:3,k1) = 0.d0
        do l=1,3
          GradH(1:3,k1) = GradH(1:3,k1) + Piola_curl(1:3,l)*grad(l)
        enddo
      enddo
c
c  ...vertical faces
      do ie=1,4
        if=ie+1
        if (Norder(8+if).le.2) cycle
        call shape_triaH_orient((/alpha(1,ie),alpha(2,ie),ShapH(5)/),
     .      (/gradalpha(1:3,1,ie),gradalpha(1:3,2,ie),GradH(1:3,5)/),
     .                               3,Norder(8+if),Norient_face(if),
     .                         nrdof_faceH,ShapH(k+1),GradH(1:3,k+1))
c
c  .....multiply by the common factor blending function
        do j=1,nrdof_faceH
          k=k+1
          GradH(1:3,k) = GradH(1:3,k)*psi(ie) + ShapH(k)*gradpsi(1:3,ie)
          ShapH(k) = ShapH(k)*psi(ie)
        enddo
      enddo
c
c  ...bubbles
      ksave = k
      nord = Norder(14)
      call ILegendre(   x,1.d0,nord,0, rL1,dLdx1,void)
      call ILegendre(   y,1.d0,nord,0, rL2,dLdx2,void)
      call IJacobi(zeta,1.d0,nord+1,nord,0, rL3,dLdx3,void)


      if (iprint.eq.2) then
        rL1(2)=x*(1.d0-x); dLdx1(2)=1.d0-2.d0*x
        rL2(2)=y*(1.d0-y); dLdx2(2)=1.d0-2.d0*y
        rL3(1:3,1)=zeta; dLdx3(1:3,1) = 1.d0
      endif



c
      do m=1,nord-1
        do j=2,nord
          do i=2,nord
            ij1=max(i,j)+1
            n = max(i,j,m+1)
            fact  = (1.d0 -zeta)**(n-m)
            dfact = - (n-m)*(1.d0 -zeta)**(n-m-1)
            k=k+1
            ShapH(k)   = rL1(i)  *rL2(j)  *rL3(ij1,m)  *fact
            GradH(1,k) = dLdx1(i)*rL2(j)  *rL3(ij1,m)  *fact
            GradH(2,k) = rL1(i)  *dLdx2(j)*rL3(ij1,m)  *fact
            GradH(3,k) = rL1(i)*rL2(j)
     .                 *(dLdx3(ij1,m)*fact + rL3(ij1,m)*dfact)
          enddo
        enddo
      enddo

cc      if (Nord.ge.2) then
cc      write(*,*) 'ShapH(k) = ',ShapH(k),
cc     .            x*(1.d0-x)*y*(1.d0-y)*zeta*(1.d0-zeta)
cc      write(*,*) 'GradH(1,k) = ',GradH(1,k),
cc     .            (1.d0-2.d0*x)*y*(1.d0-y)*zeta*(1.d0-zeta)
cc      write(*,*) 'GradH(2,k) = ',GradH(2,k),
cc     .            x*(1.d0-x)*(1.d0-2.d0*y)*zeta*(1.d0-zeta)
cc      write(*,*) 'GradH(3,k) = ',GradH(3,k),
cc     .            x*(1.d0-x)*y*(1.d0-y)*(1.d0-2.d0*zeta)
cc      call pause
cc      endif
c
c  ...transform the gradients to the xi,eta,zeta frame
      NrdofH = k
      do k=ksave+1,NrdofH
        grad(1:3) = GradH(1:3,k)
        GradH(1:3,k) = 0.d0
        do l=1,3
          GradH(1:3,k) = GradH(1:3,k) + Piola_curl(1:3,l)*grad(l)
        enddo
      enddo
c
c
      if (iprint.ge.1) then
        write(*,7001) Xip(1:3),Norder(1:14),
     .                Norient_edge(1:8),Norient_face(1:5)
 7001   format('shape3H_pyra: Xip = ',3f8.3,/,
     .         'Norder = ',8i2,2x,i3,2x,4i2,2x,i2,/,
     .         'orient = ',8i2,2x,i3,2x,4i2)
        write(*,7010)
 7010   format('VERTEX SHAPE FUNCTIONS = ')
        do iv=1,5
          k=iv
          write(*,7002) k,ShapH(k),GradH(1:3,k)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,3e12.5)
        enddo
        do ie=1,8
          nord = Norder(ie)
          if (nord.eq.0) cycle
          ndof = nord-1
          if (ndof.eq.0) cycle
          write(*,7011) ie
 7011     format('SHAPE FUNCTIONS FOR EDGE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        enddo
        if (Norder(9).gt.0) then
          call decode(Norder(9), nord1,nord2)
          ndof = (nord1-1)*(nord2-1)
          if (ndof.gt.0) then
            write(*,7012)
 7012       format('SHAPE FUNCTIONS FOR BASE FACE ')
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapH(k),GradH(1:3,k)
            enddo
          endif
        endif
        do if=2,5
          nord = Norder(8+if)
          if (nord.eq.0) cycle
          ndof = (nord-1)*(nord-2)/2
          if (ndof.eq.0) cycle
          write(*,7013) if
 7013     format('SHAPE FUNCTIONS FOR VERTICAL FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7002) k,ShapH(k),GradH(1:3,k)
          enddo
        enddo
        nord = Norder(14)
        if (nord.gt.0) then
          ndof = (nord-1)**3
          if (ndof.gt.0) then
            write(*,7014)
 7014       format('BUBBLES = ')
            do j=1,ndof
              k=k+1
              write(*,7002) k,ShapH(k),GradH(1:3,k)
            enddo
          endif
        endif
        call pause
      endif
c
c
      end subroutine shape3H_pyra
c





