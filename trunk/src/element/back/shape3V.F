c--------------------------------------------------------------------
c
c     routine name      - shape3V
c
c--------------------------------------------------------------------
c
c     latest revision:  - Oct 14
c
c     purpose:          - routine returns values of 3D elements
c                         H(div) shape functions and their div
c
c     arguments:
c
c     in:
c          Type         - element type
c          Xi           - master bric coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_face - face orientation
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - div of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3V_1(Type,Xi,Norder,Norient_face,
     .                     NrdofV,ShapV,DivV)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(3),Norder(19),Norient_face(6),
     .          ShapV(3,MAXbrickV),DivV(MAXbrickV)
c
c
      select case(Type)
      case('bric','mdlb')
        call shape3V_bric(Xi,Norder,Norient_face, NrdofV,ShapV,DivV)
      case('tetr','mdln')
        call shape3V_tetr(Xi,Norder,Norient_face, NrdofV,ShapV,DivV)
      case('pris','mdlp')
        call shape3V_pris(Xi,Norder,Norient_face, NrdofV,ShapV,DivV)
      case('pyra','mdld')
        call shape3V_pyra(Xi,Norder,Norient_face, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape3V: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3V_1
c
c--------------------------------------------------------------------
c
c     routine name      - shape3V_bric
c
c--------------------------------------------------------------------
c
c     latest revision:  - Beb 14
c
c     purpose:          - routine returns values of 3D brick element
c                         H(div) shape functions and their div
c
c     arguments:
c
c     in:
c          Xi           - master bric coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_face - face orientation
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - div of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3V_bric(Xi,Norder,Norient_face,
     .                        NrdofV,ShapV,DivV)
c
      use parameters
#include "syscom.blk"
c
      dimension Xi(3),Norder(19),Norient_edge(12),Norient_face(6),
     .          ShapV(3,MaxbrickV),DivV(MAXbrickV)
c
c  ...face blending functions and projections
      dimension psif(6),gradpsif(3,6),tf(2,6),gradtf(2,3,6)
c
c  ...work space for faces
      dimension vec(3)
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
c  ...determine face blending functions and projections
      call blend_project_faces_bric(Xi, psif,gradpsif,tf,gradtf)
c
c  ...face shape functions
      do if=1,6
c
        call decode(Norder(12+if), nord1,nord2) ! order acoounts for orientation
        ndof = nord1*nord2
        if (ndof.eq.0) cycle
c
        call cross_product(gradtf(1,1:3,if),gradtf(2,1:3,if), vec)
        call dot_product(gradpsif(1:3,if),vec, prod)
c
        iswap = norient221(1,Norient_face(if))
        norient1 = norient221(2,Norient_face(if))
        norient2 = norient221(3,Norient_face(if))
        select case(iswap)
        case(0); xif1 = tf(1,if); xif2 = tf(2,if)
        case(1); xif1 = tf(2,if); xif2 = tf(1,if)
        end select
c
c  .....needed for Piola transform
        select case(Norient_face(if))
        case(0,1,2,3); det =  1.d0
        case(4,5,6,7); det = -1.d0
        end select
c
        call shape1Q_orient(xif1,nord1,norient1, nrdofQ1,shapQ1)
        call shape1Q_orient(xif2,nord2,norient2, nrdofQ2,shapQ2)
        do j=1,nrdofQ2
          do i=1,nrdofQ1
            shapQf = shapQ1(i)*shapQ2(j)*det
            k=k+1
            ShapV(1:3,k) = psif(if)*shapQf*vec(1:3)
            DivV(k) = shapQf*prod
          enddo
        enddo
      enddo
c
c  ...bubbles
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
      do l=1,nrdofQ3
        do j=1,nrdofQ2
          do i=1,nrdofH1
            k=k+1
            ShapV(1,k) = shapH1(i)*shapQ2(j)*shapQ3(l)
            ShapV(2,k) = 0.d0
            ShapV(3,k) = 0.d0
            DivV(k) = dshapH1(i)*shapQ2(j)*shapQ3(l)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the y-axis
      do l=1,nrdofQ3
        do j=1,nrdofH2
          do i=1,nrdofQ1
            k=k+1
            ShapV(1,k) = 0.d0
            ShapV(2,k) = shapQ1(i)*shapH2(j)*shapQ3(l)
            ShapV(3,k) = 0.d0
            DivV(k) = shapQ1(i)*dshapH2(j)*shapQ3(l)
          enddo
        enddo
      enddo
c
c  ...shape functions with values along the z-axis
      do l=1,nrdofH3
        do j=1,nrdofQ2
          do i=1,nrdofQ1
            k=k+1
            ShapV(1,k) = 0.d0
            ShapV(2,k) = 0.d0
            ShapV(3,k) = shapQ1(i)*shapQ2(j)*shapH3(l)
            DivV(k) = shapQ1(i)*shapQ2(j)*dshapH3(l)
          enddo
        enddo
      enddo
c
      NrdofV=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(13:19),Norient_face(1:6)
 7001   format('shape3V_bric: Xi = ',3f8.3,/,
     .         'Norder       = ',2i3,2x,4i3,3x,i4,/,
     .         'Norient_face = ',2i3,2x,4i3)
        k=0
        do if=1,6
          call decode(Norder(12+if), nordh,nordv)
          ndof = nordh*nordv
          if (ndof.gt.0) then
            write(*,7002) if
 7002       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapV(1:3,k),DivV(k)
 7003         format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
            enddo
          endif
        enddo
        call decode(Norder(19), nord,nord3)
        call decode(nord, nord1,nord2)
        ndof = (nord1-1)* nord2   * nord3
     .       +  nord1   *(nord2-1)* nord3
     .       +  nord1   * nord2   *(nord3-1)

        if (ndof.gt.0) then
          write(*,7004)
 7004     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3V_bric
c
c--------------------------------------------------------------------
c
c     routine name      - shape3V_tetr
c
c--------------------------------------------------------------------
c
c     latest revision:  - Nov 13
c
c     purpose:          - routine returns values of 3D tetra element
c                         H(div) shape functions and their div
c
c     arguments:
c
c     in:
c          Xi           - master tetr coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_face - face orientation
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - div of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3V_tetr(Xi,Norder,Norient_face,
     .                        NrdofV,ShapV,DivV)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(3),Norder(11),Norient_face(6),
     .          ShapV(3,MAXtetraV),DivV(MAXtetraV)
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
c  ...face shape functions
      do if=1,4
        nv1 = TETRA_FACE_TO_VERT(1,if)
        nv2 = TETRA_FACE_TO_VERT(2,if)
        nv3 = TETRA_FACE_TO_VERT(3,if)
        call shape_triaV_orient((/vshap(nv1),vshap(nv2),vshap(nv3)/),
     .              (/vgrad(1:3,nv1),vgrad(1:3,nv2),vgrad(1:3,nv3)/),
     .                                 Norder(6+if),Norient_face(if),
     .                           nrdofV_face,ShapV(1,k+1),DivV(k+1))
        k=k+nrdofV_face
      enddo
c
c  ...bubbles
      call shape_tetrV(vshap,vgrad,Norder(11),
     .                 NrdofV,ShapV(1,k+1),DivV(k+1))
      NrdofV = k+NrdofV
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(7:11),Norient_face(1:4)
 7001   format('shape3V_tetr: Xi = ',3f8.3,/,
     .         'Norder       = ',4i3,3x,i3,/,
     .         'Norient_face = ',4i3)
        k=0
        do if=1,4
          nord=Norder(6+if)
          ndof = nord*(nord+1)/2
          if (ndof.gt.0) then
            write(*,7002) if
 7002       format('SHAPE FUNCTIONS FOR FACE = ',i2)
            do j=1,ndof
              k=k+1
              write(*,7003) k,ShapV(1:3,k),DivV(k)
 7003         format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
            enddo
          endif
        enddo
        nord = Norder(11)
        ndof = (nord-1)*nord*(nord+1)/2
        if (ndof.gt.0) then
          write(*,7004)
 7004     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3V_tetr
c
c--------------------------------------------------------------------
c
c     routine name      - shape3V_pris
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 3D prism element
c                         H(div) shape functions and their div
c
c     arguments:
c
c     in:
c          Xi           - master pris coordinates from (0,1)^3
c          Norder       - polynomial order
c          Norient_face - face orientation
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - div of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3V_pris(Xi,Norder,Norient_face,
     .                        NrdofV,ShapV,DivV)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xi(3),Norder(15),Norient_edge(9),Norient_face(5),
     .          ShapV(3,MAXprismV),DivV(MAXprismV)
c
c  ...triangle vertex shape functions, 1D shape functions
      dimension shapHt(3),gradHt(2,3),shapH1(2),dshapH1(2)
c
c  ...work space for shape2Q_tria_orient, shape1Q_orient
      dimension shapQt(MAXmdltQ), shapQ2(MAXP)
c
c  ...work space for shape_edgeE_orient, shape_triaE
      dimension shapEt(2,MAXmdltE),curlEt(MAXmdltE)
c
c  ...work space for shape1Hb
      dimension shapH1b(MAXP-1),dshapH1b(MAXP-1)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
c  ...triangle and 1D vertex shape functions
      shapHt(1) = 1.d0-Xi(1)-Xi(2); gradHt(1:2,1)=-1.d0
      shapHt(2) = Xi(1); gradHt(1,2) = 1.d0; gradHt(2,2) = 0.d0
      shapHt(3) = Xi(2); gradHt(1,3) = 0.d0; gradHt(2,3) = 1.d0
      shapH1(1) = 1.d0-Xi(3); dshapH1(1) = -1.d0
      shapH1(2) =      Xi(3); dshapH1(2) =  1.d0
c
      k=0
c
c  ...horizontal faces
      if=0
      do j=1,2
        if=if+1
        call shape2Q_tria_orient(Xi,Norder(9+if),Norient_face(if),
     .                           nrdofQt,shapQt)
        do kf=1,nrdofQt
          k=k+1
          ShapV(1:2,k) = 0.d0
          ShapV(3,k) = shapQt(kf)*shapH1(j)
          DivV(k) = shapQt(kf)*dshapH1(j)
        enddo
      enddo
c
c  ...vertical faces
      do ie=1,3
        if=if+1
        call decode(Norder(9+if), nord1,nord2)
c
        nv1 = TRIAN_EDGE_TO_VERT(1,ie)
        nv2 = TRIAN_EDGE_TO_VERT(2,ie)
c
        iswap = norient221(1,Norient_face(if))
        norient1 = norient221(2,Norient_face(if))
        norient2 = norient221(3,Norient_face(if))
c
c  .....needed for Piola transform
        select case(Norient_face(if))
        case(0,1,2,3); det =  1.d0
        case(4,5,6,7); det = -1.d0
        end select
c

        select case(iswap)
        case(0)
          call shape_edgeE_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                              2,nord1,norient1,
     .                              nrdofEt,shapEt,curlEt)
          call shape1Q_orient(Xi(3),nord2,norient2, nrdofQ2,shapQ2)
          do j=1,nrdofQ2
            do i=1,nrdofEt
              k=k+1
c
c  ...........rotate the 2D H(curl) functions
              ShapV(1,k) =   shapEt(2,i)*shapQ2(j)*det
              ShapV(2,k) = - shapEt(1,i)*shapQ2(j)*det
              ShapV(3,k) = 0.d0
              DivV(k) = curlEt(i)*shapQ2(j)*det
            enddo
          enddo
c
c  .....swapped face axes
        case(1)
          call shape_edgeE_orient((/shapHt(nv1),shapHt(nv2)/),
     .                            (/gradHt(1:2,nv1),gradHt(1:2,nv2)/),
     .                              2,nord2,norient2,
     .                              nrdofEt,shapEt,curlEt)
          call shape1Q_orient(Xi(3),nord1,norient1, nrdofQ2,shapQ2)
c
c  .......change ordering
          do i=1,nrdofEt
            do j=1,nrdofQ2
              k=k+1
c
c  ...........rotate the 2D H(curl) functions
              ShapV(1,k) =   shapEt(2,i)*shapQ2(j)*det
              ShapV(2,k) = - shapEt(1,i)*shapQ2(j)*det
              ShapV(3,k) = 0.d0
              DivV(k) = curlEt(i)*shapQ2(j)*det
            enddo
          enddo
        end select
c
c  ...end of loop through vertical faces
      enddo
c
c  ...bubbles...........................................................
      call decode(Norder(15), nord1,nord2)
c
c  ...determine triangle H(curl) bubbles and L2 shape functions
      call shape_triaE(shapHt,gradHt,2,nord1, nrdofEt,shapEt,curlEt)
      call shape2Q_tria(Xi,nord1, nrdofQt,shapQt)
c
c  ...determine 1D H1 bubbles and L2 shape functions
      call shape1Hb(Xi(3),nord2, nrdofH1b,shapH1b,dshapH1b)
      call shape1Q(Xi(3),nord2, nrdofQ2,shapQ2)
c
c  ...first family
      do j=1,nrdofQ2
        do i=1,nrdofEt
          k=k+1
c
c  .......rotate 2D H(curl) functions
          ShapV(1,k) =   shapEt(2,i)*shapQ2(j)
          ShapV(2,k) = - shapEt(1,i)*shapQ2(j)
          ShapV(3,k) = 0.d0
          DivV(k) = curlEt(i)*shapQ2(j)
        enddo
      enddo
c
c  ...second family
      do j=1,nrdofH1b
        do i=1,nrdofQt
          k=k+1
          ShapV(1:2,k) = 0.d0
          ShapV(3,k) = shapQt(i)*shapH1b(j)
          DivV(k) = shapQt(i)*dshapH1b(j)
        enddo
      enddo
c
      NrdofV=k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3),Norder(10:15),Norient_face(1:5)
 7001   format('shape3V_pris: Xi = ',3f8.3,/,
     .         'Norder       = ',2i3,2x,3i3,3x,i3,/,
     .         'Norient_face = ',2i3,2x,3i3)
        k=0
        do if=1,2
          nord = Norder(9+if)
          ndof = nord*(nord+1)/2
          write(*,7002) if
 7002     format('SHAPE FUNCTIONS FOR BASE FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
 7003       format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
          enddo
        enddo
        do if=3,5
          call decode(Norder(9+if), nordh,nordv)
          ndof = nordh*nordv
          write(*,7004) if-2
 7004     format('SHAPE FUNCTIONS FOR VERTICAL FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
          enddo
        enddo
        call decode(Norder(15), nord1,nord2)
        ndof = nord1*(nord1-1)*nord2 + nord1*(nord1+1)/2*(nord2-1)
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3V_pris
c
c--------------------------------------------------------------------
c
c     routine name      - shape3V_pyra
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of 3D pyramid element
c                         H(div) shape functions and their div
c
c     arguments:
c
c     in:
c          Xip          - master pyramid coordinates
c          Norder       - polynomial order
c          Norient_face - face orientation
c
c     out:
c          NrdofV       - number of dof
c          ShapV        - values of the shape functions at the point
c          DivV         - div of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape3V_pyra(Xip,Norder,Norient_face,
     .                        NrdofV,ShapV,DivV)
c
      use parameters
      use element_data
#include "syscom.blk"
c
      dimension Xip(3),Norder(14),Norient_face(5),
     .          ShapV(3,MAXpyramV),DivV(MAXpyramV)
c
c  ...common factor, modified coordinates, vertex shape functions
      dimension psi(4),gradpsi(3,4),alpha(2,4),gradalpha(3,2,4),
     .          vshap(5),gradv(3,5),vec(3),vec1(3)
c
c  ...work space for shape1Q_orient
      dimension shapQ1(MAXP), shapQ2(MAXP)
c
      double precision zeta,z
      integer a
c
c  ...gradients of infinite element coordinates
      dimension gradx(3),grady(3),gradz(3)
c
c  ...work space for ILegendre and Legendre
      dimension rL1(2:MAXP),dLdx1(2:MAXP),p1(0:MAXP),
     .          rL2(2:MAXP),dLdx2(2:MAXP),p2(0:MAXP)
c
c  ...work space for Jacobi and IJacobi
      dimension rL3(1:MAXP+1,1:MAXP),dLdx3(1:MAXP+1,1:MAXP),
ccc P. Gatto, Nov 14, 2014
ccc      dimension rL3(1:MAXP+1,1:MAXP-1),dLdx3(1:MAXP+1,1:MAXP-1),
     .          p3(0:MAXP+1,0:MAXP)
c
c  ...work space for the pullback
      dimension Piola_div(3,3), Vvec(3)
c
c  ...translating quad orientation into product of 1D orientations
      dimension norient221(1:3,0:7)
      data norient221/0,0,0, 1,1,0, 0,1,1, 1,0,1,
     .                1,0,0, 0,1,0, 1,1,1, 0,0,1/
c
      iprint=0
c
      xi = Xip(1); eta = Xip(2); zeta = Xip(3)
      call blend_pyram(Xip,x,gradx,y,Grady,z,Gradz,
     .                 psi,gradpsi,alpha,gradAlpha,
     .                 vshap,gradv)
c
c  ...pull  back from x,y,zeta to xi,eta,zeta
      Piola_div(1,1:3) = (/1.d0-zeta, 0.d0, -x/)
      Piola_div(2,1:3) = (/0.d0, 1.d0-zeta, -y/)
      Piola_div(3,1:3) = (/0.d0,0.d0,1.d0/)
      rjac = (1.d0-zeta)**2
      Piola_div = Piola_div/rjac
c
      k=0
c
c  ...base face
      if=1
      call decode(Norder(9), nord1,nord2)
      iswap = norient221(1,Norient_face(if))
      norient1 = norient221(2,Norient_face(if))
      norient2 = norient221(3,Norient_face(if))
      select case(iswap)
      case(0); xif1 = x; xif2 = y
      case(1); xif1 = y; xif2 = x
      end select
c
c  ...needed for Piola transform
      select case(Norient_face(if))
      case(0,1,2,3); det =  1.d0
      case(4,5,6,7); det = -1.d0
      end select
c
      call shape1Q_orient (xif1,nord1,norient1, nrdofQ1,shapQ1)
      call shape1Q_orient (xif2,nord2,norient2, nrdofQ2,shapQ2)
c
      ksave=k
      do j=1,nrdofQ2
        do i=1,nrdofQ1
          n = max(i,j)+2 ! order in x,y is i+1,j+1
          fact = (1.d0-zeta)**n
          dfact = -n*(1.d0-zeta)**(n-1) ! watch for the minus sign
          k=k+1
          ShapV(1:2,k) = 0.d0
          ShapV(3,k) = shapQ1(i)*shapQ2(j)*det* fact
          DivV(k)    = shapQ1(i)*shapQ2(j)*det*dfact
        enddo
      enddo
c
c  ...pull back base shape functions from x,y,zeta to xi,eta,zeta
      do k1=ksave+1,k
        Vvec(1:3) = ShapV(1:3,k1)
        ShapV(1:3,k1) = 0.d0
        do l=1,3
          ShapV(1:3,k1) = ShapV(1:3,k1) + Piola_div(1:3,l)*Vvec(l)
        enddo
        DivV(k1) = DivV(k1)/rjac
      enddo
c
c  ...vertical faces
      do ie=1,4
        if=ie+1
        call shape_triaV_orient((/alpha(1,ie),alpha(2,ie),vshap(5)/),
     .      (/gradalpha(1:3,1,ie),gradalpha(1:3,2,ie),gradv(1:3,5)/),
     .                                 Norder(8+if),Norient_face(if),
     .                          nrdof_faceV,ShapV(1:3,k+1),DivV(k+1))
c
c  .....multiply by the common factor blending function
        do j=1,nrdof_faceV
          k=k+1
          call dot_product(gradpsi(1:3,ie),ShapV(1:3,k), s)
          DivV(k) = DivV(k)*psi(ie) + s
          ShapV(1:3,k) = ShapV(1:3,k)*psi(ie)
        enddo
      enddo
c
c  ...bubbles
      ksave = k
      nord = Norder(14)
      call ILegendre(   x,1.d0,nord,0, rL1,dLdx1,void)
      call ILegendre(   y,1.d0,nord,0, rL2,dLdx2,void)
      call IJacobi(zeta,1.d0,Nord,Nord,0, rL3,dLdx3,void)
      call LegendreNEW(   x,1.d0,nord, p1)
      call LegendreNEW(   y,1.d0,nord, p2)
      maxalpha = 2*Nord1
      call Jacobi(zeta,1.d0,maxalpha,Nord, p3)
c
c  ...first family
      do m=0,nord-1
        do j=0,nord-1
          do i=2,nord
            n = max(i,j+1,m+1) ! corrected
            a = max(i,j+1)-1
c            jalpha = 2*max(i,j+1)-1 ! Jacobi weight, alpha is already used
            k=k+1
            ShapV(1,k) = rL1(i)*P2(j)*P3(a,m)*(1.d0-zeta)**(n-m)
            ShapV(2:3,k) = 0.d0
            Divv(k) =  dLdx1(i)*P2(j)*P3(a,m)*(1.d0-zeta)**(n-m)
          enddo
        enddo
      enddo
c
c  ...second family
      do m=0,nord-1
        do j=2,nord
          do i=0,nord-1
            n = max(i+1,j,m+1) ! corrected
            a = max(i+1,j)-1
c            jalpha = 2*max(i+1,j)-1
            k=k+1
            ShapV(1,k) = 0.d0; ShapV(3,k) = 0.d0
            ShapV(2,k) = P1(i)*rL2(j)*P3(a,m)*(1.d0-zeta)**(n-m)
            Divv(k) =  P1(i)*dLdx2(j)*P3(a,m)*(1.d0-zeta)**(n-m)
          enddo
        enddo
      enddo
c
c  ...third family
      do m=1,nord  ! corrected
        do j=0,nord-1
          do i=0,nord-1
            n = max(i+1,j+1,m)+2
            index = max(i+1,j+1) !determines weight for rL3
            k=k+1
            ShapV(1:2,k) = 0.d0
            ShapV(3,k) = P1(i)*P2(j)*rL3(index,m)*(1.d0-zeta)**(n-m)
            Divv(k) =  P1(i)*P2(j)*
     .                (dLdx3(index,m)*(1.d0-zeta)**(n-m)
     .                -(n-m)*rL3(index,m)*(1.d0-zeta)**(n-m-1))
          enddo
        enddo
      enddo
c
c  ...pullback element bubbles
      do k1=ksave+1,k
        Vvec(1:3) = ShapV(1:3,k1)
        ShapV(1:3,k1) = 0.d0
        do l=1,3
          ShapV(1:3,k1) = ShapV(1:3,k1) + Piola_div(1:3,l)*Vvec(l)
        enddo
        DivV(k1) = DivV(k1)/rjac
      enddo
      NrdofV = k
c
c
      if (iprint.eq.1) then
        write(*,7001) Xip(1:3),Norder(9:14),Norient_face(1:5)
 7001   format('shape3V_pyra: Xip = ',3f8.3,/,
     .         'Norder       = ',i3,2x,4i3,3x,i2,/,
     .         'Norient_face = ',i3,2x,4i3)
        k=0
        call decode(Norder(9), nord1,nord2)
        ndof = nord1*nord2
        write(*,7002)
 7002   format('SHAPE FUNCTIONS FOR BASE FACE ')
        do j=1,ndof
          k=k+1
          write(*,7003) k,ShapV(1:3,k),DivV(k)
 7003     format('k = ',i3,' ShapV, DivV = ',3e12.5,3x,e12.5)
        enddo
        do if=2,5
          nord=Norder(8+if)
          ndof = nord*(nord+1)/2
          if (ndof.eq.0) cycle
          write(*,7004) if
 7004     format('SHAPE FUNCTIONS FOR FACE = ',i2)
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
          enddo
        enddo
        nord = Norder(14)
        ndof = nord*(nord+1)*(nord+2)/6
        if (ndof.gt.0) then
          write(*,7005)
 7005     format('BUBBLES = ')
          do j=1,ndof
            k=k+1
            write(*,7003) k,ShapV(1:3,k),DivV(k)
          enddo
        endif
        call pause
      endif
c
c
      end subroutine shape3V_pyra



