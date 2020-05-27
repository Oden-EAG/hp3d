c----------------------------------------------------------------------
c
c   routine name       - shape3nE
c
c----------------------------------------------------------------------
c
c   latest revision    - Sep 10
c
c   purpose            - 3D Nedelec's tetrahedron of the 1st type,
c                        routine evaluates H(curl) shape functions
c                        and derivatives using Sabine's construction
c                        based on scaled integrated Legendre
c                        polynomials and integrated Legendre
c                        polynomials
c
c
c   arguments :
c     in:
c            Xi        - master element coordinates
c            Norder    - order of approximation for the nodes
c            Ne_orient - edge orientations
c            Nf_orient - face orientations
c     out:
c            NrdofE    - number of element dof
c            ShapE     - values of shape functions
c            CurlE     - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shape3nE(Xi,Norder,Ne_orient,Nf_orient,
     .                    NrdofE,ShapE,CurlE)
c
      use element_data
      use parameters
      use control
#include "syscom.blk"
c
c  ...allow for exceeding the order of approximation
      integer, parameter :: ndelta = 2
c
c  ...input output parameters
      dimension Xi(3),
     .          Norder(19), Ne_orient(12), Nf_orient(6),
     .          ShapE(3,MAXbrickE),CurlE(3,MAXbrickE)
c
c  ...local variables
      dimension jtemp(3)
      dimension dp_dXi(3),dq_dXi(3),ds_dXi(3),dt_dXi(3)
      dimension aux1(3),aux2(3)
c
c  ...block u
      dimension uu(0:MAXP+ndelta),duu(2,0:MAXP+ndelta)
      dimension duu_dXi(3,0:MAXP+ndelta)
c
c  ...block v
      dimension vv(0:MAXP+ndelta),dvv(2,0:MAXP+ndelta)
      dimension dvv_dXi(3,0:MAXP+ndelta)
c
c  ...polynomial blocks u,v,w
      dimension poly_u(1:MAXP+ndelta),dpoly_u(1:3,1:MAXP+ndelta)
      dimension poly_v(1:MAXP+ndelta,1:MAXP+ndelta),
     .          dpoly_v(1:3,1:MAXP+ndelta,1:MAXP+ndelta)
      dimension poly_w(1:MAXP+ndelta,1:MAXP+ndelta,1:MAXP+ndelta),
     .        dpoly_w(1:3,1:MAXP+ndelta,1:MAXP+ndelta,1:MAXP+ndelta)
c
c  ...affine coordinates
      dimension shapH(4),dshapH(1:3,4)
c
c----------------------------------------------------------------------
c
      iprint=0
c
      ShapE(1:3,1:MAXbrickE) = 0.d0
      CurlE(1:3,1:MAXbrickE) = 0.d0
c
c  ...check if the input point is within the master tet
      if (((Xi(1)+Xi(2)+Xi(3)).gt.1.d0+GEOM_TOL).or.
     .    (Xi(1).lt.-GEOM_TOL).or.
     .    (Xi(2).lt.-GEOM_TOL).or.
     .    (Xi(3).lt.-GEOM_TOL)) then
        write(*,7001)
 7001   format(' shape3nE: WARNING !! POINT Xi OUTSIDE OF MASTER TET')
        write(*,7002) Xi(1:3)
 7002   format(' shape3nE: Xi = ',3(e12.5,2x))
        call pause
      endif
c
c  ...calculate the affine coordinates and their gradients
      shapH(1)      =  1.d0 - Xi(1) - Xi(2) - Xi(3)
      dshapH(1:3,1) = -1.d0
      do i=1,3
        shapH(1+i)      = Xi(i)
        dshapH(1:3,1+i) = 0.d0; dshapH(i,1+i) = 1.d0
      enddo
c
c
c======================================================================
c     F  I  R  S  T      T  Y  P  E                                   |
c======================================================================
c
c  ...initiate shape functions counter
      k=0
c
c  ...mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop over edges
      do i=1,6
c
c  .....determine edge # of dof
        ndofE = Norder(i)
c
        if (ndofE.gt.0) then
c
c  .......get vertices specifing local edge orientation
          iv1 = TETRA_EDGE_TO_VERT(1,i); iv2 = TETRA_EDGE_TO_VERT(2,i)
c
c  .......determine endpoints in the global orientation

          if (Ne_orient(i).eq.1) call swap(iv1,iv2)

c
c  .......edge Whitney shape function
          k=k+1
          ShapE(1:3,k) = dshapH(1:3,iv1)*shapH(iv2) -
     .                   dshapH(1:3,iv2)*shapH(iv1)
          call cross_product(dshapH(1:3,iv2),dshapH(1:3,iv1),
     .                       CurlE(1:3,k))
          CurlE(1:3,k) = 2.d0*CurlE(1:3,k)
c
c  .......edge coordinates
          p      = shapH(     iv1) - shapH(     iv2)
          dp_dXi = dshapH(1:3,iv1) - dshapH(1:3,iv2)
          q      = shapH(     iv1) + shapH(     iv2)
          dq_dXi = dshapH(1:3,iv1) + dshapH(1:3,iv2)
c
c  .......compute Scaled Integrated Legendre polynomials:
c         L_i^{S}  ;  i = 0,Norder
          call Legendre_S_I(p,q,Norder(i)+1, uu,duu)
c
c  .......compute shape functions and their curl
          do j=1,ndofE-1
            k=k+1
            ShapE(1:3,k) = duu(1,j+1)*dp_dXi(1:3) +
     .                     duu(2,j+1)*dq_dXi(1:3)
c  -------> CurlE(1:3,k) = 0
          enddo
        endif
c
c  ...end of loop over edges
      enddo
c
c  ...mid-face nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop over faces
      do ii=1,4
c
        if (Norder(6+ii).ge.2) then
c
c  .......get face vertices
          do jj=1,3
            jtemp(jj) = TETRA_FACE_TO_VERT(jj,ii)
          enddo
c
c  .......arrange vertices according to face global orientation
          iv1 = jtemp(TRIAN_L2G(1,Nf_orient(ii)))
          iv2 = jtemp(TRIAN_L2G(2,Nf_orient(ii)))
          iv3 = jtemp(TRIAN_L2G(3,Nf_orient(ii)))
c
c  .......compute Scaled Integrated Legendre polynomials
          p      = shapH(     iv1) - shapH(     iv2)
          dp_dXi = dshapH(1:3,iv1) - dshapH(1:3,iv2)
          q      = shapH(     iv1) + shapH(     iv2)
          dq_dXi = dshapH(1:3,iv1) + dshapH(1:3,iv2)
c         L_i^{S}  ;  i = 0,Norder
          call Legendre_S_I(p,q,Norder(6+ii), uu,duu)
c
c  .......compute Scaled Legendre polynomials
          s      = -shapH(     iv1) - shapH(     iv2) + shapH(     iv3)
          ds_dXi = -dshapH(1:3,iv1) - dshapH(1:3,iv2) + dshapH(1:3,iv3)
          t      =  shapH(     iv1) + shapH(     iv2) + shapH(     iv3)
          dt_dXi =  dshapH(1:3,iv1) + dshapH(1:3,iv2) + dshapH(1:3,iv3)
c         l_i^{S}  ;  i =  0,Norder
          call Legendre_S(s,t,Norder(6+ii), vv,dvv)
c
c  .......compute block vv, \nabla uu, \nabla vv
          do i=0,MAXP+ndelta
            duu_dXi(1:3,i) =duu(1,i)*dp_dXi(1:3) + duu(2,i)*dq_dXi(1:3)
            dvv_dXi(1:3,i) =dshapH(1:3,iv3)*vv(i) + shapH(iv3)*(
     .                      dvv(1,i)*ds_dXi(1:3) + dvv(2,i)*dt_dXi(1:3))
            vv(i)          =shapH(iv3)*vv(i)
          enddo
c
c  .......loop through the total order first (ij+2)
          do ij=0,Norder(6+ii)-2
ccc            write(*,*) 'TOTAL ORDER ij+2 = ',ij+2
            j=ij
c
c  .........COMPLETION with TYPE III....................................
ccc            write(*,*) 'TYPE III: j = ',j
            k=k+1
            ShapE(1:3,k) = (dshapH(1:3,iv1)*shapH(iv2) -
     .                      dshapH(1:3,iv2)*shapH(iv1))*vv(j)
c
c===========================================================================
c  ...curl TYPE III =
c     \nabla v x (\lambda_2 \nabla \lambda_1 - \lambda_1 \nabla \lambda_2) +
c            2 v (\nabla \lambda_2 x \nabla \lambda_1) =
c     \nabla v x aux2 + 2v * aux1
c===========================================================================
c
c           aux2 = (\lambda_2 \nabla \lambda_1 - \lambda_1 \nabla \lambda_2)
            aux2(1:3) = shapH(iv2)*dshapH(1:3,iv1) -
     .                  shapH(iv1)*dshapH(1:3,iv2)
c           aux1 = \nabla \lamda_2 x \nabla \lambda_1
            call cross_product(dshapH(1:3,iv2),dshapH(1:3,iv1),
     .                         aux1(1:3))
c           \nabla v x aux1
            call cross_product(dvv_dXi(1:3,j),aux2(1:3), CurlE(1:3,k))
c           \nabla v x aux2 + 2v * aux1
            CurlE(1:3,k) = CurlE(1:3,k) + 2.d0*vv(j)*aux1(1:3)
c
c
c  .........COMPLETION with TYPE II....................................
            do i=0,ij
              j=ij-i
ccc              write(*,*) 'TYPE II: i,j = ',i,j
              k=k+1
              ShapE(1:3,k) = duu_dXi(1:3,i+2)* vv(        j) -
     .                        uu(        i+2)*dvv_dXi(1:3,j)
c
c  ...........curl TYPE II = 2 \nabla v x \nabla u
              call cross_product(dvv_dXi(1:3,j),
     .                           duu_dXi(1:3,i+2), CurlE(1:3,k))
              CurlE(1:3,k) = 2.d0*CurlE(1:3,k)
            enddo
c
c
c  .........COMPLETION with TYPE I (gradients).........................
            do i=0,ij
              j=ij-i
c
c  ...........(tet of the first type !)
              if (ij.ne.Norder(6+ii)-2) then
ccc                write(*,*) 'TYPE III: i,j = ',i,j
                k=k+1
                ShapE(1:3,k) = duu_dXi(1:3,i+2)* vv(        j) +
     .                          uu(        i+2)*dvv_dXi(1:3,j)
c  -----------> CurlE(1:3,k) = 0
              endif
            enddo
c
c  .......end of loop over total order
          enddo
c
c  .....end if (face_order.ge.2)
        endif
c
c  ...end of loop over faces
      enddo
c
c
c
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nordp = Norder(6+4+1)
      if (nordp.le.2) go to 999
      call build_u(shapH(1:2),dshapH(1:3,1:2),nordp+1, poly_u,dpoly_u)
      call build_v(shapH(3:4),dshapH(1:3,3:4),nordp+1, poly_v,dpoly_v)
      call build_w(shapH(4),  dshapH(1:3,4)  ,nordp+1, poly_w,dpoly_w)
ccc      write(*,7240) poly_u(2),dpoly_u(1:3,2)
ccc      write(*,7240) poly_v(2,1),dpoly_v(1:3,2,1)
ccc      write(*,7240) poly_w(2,1,1),dpoly_w(1:3,2,1,1)
ccc      write(*,*) '-------------------------------'
ccc      write(*,7240) -2.d0*shapH(1)*shapH(2),
ccc     .            -2.d0*(dshapH(1:3,1)*shapH(2)+shapH(1)*dshapH(1:3,2))
ccc      write(*,7240) 2.d0*shapH(3),2.d0*dshapH(1:3,3)
ccc      write(*,7240) 2.d0*shapH(4),2.d0*dshapH(1:3,4)
ccc 7240 format('shape3E: u, grad u = ',e12.5,3x,3e12.5)
ccc      call pause
c
c  ...loop through the total order of polynomials
      do ijkk=4,nordp+1
c
ccc        write(*,*) 'TOTAL ORDER ijkk-1 = ',ijkk-1
c
c  .....loop over index k
        do kk=1,ijkk-3
c  .......TYPE II......................................................
          j=ijkk-kk-2
          k=k+1
ccc          write(*,*) 'TYPE II: k, i,j,kk = ',k,2,j,kk
          ShapE(1:3,k) = (shapH(2)*dshapH(1:3,1) -
     .                    shapH(1)*dshapH(1:3,2))*
     .                    poly_v(1,j)*poly_w(1,j,kk)
          aux1(1:3) = shapH(2)*dshapH(1:3,1) -
     .                shapH(1)*dshapH(1:3,2)
          call cross_product(dpoly_v(1:3,1,j),aux1(1:3),
     .                       CurlE(1:3,k))
          CurlE(1:3,k) = poly_w(1,j,kk)*CurlE(1:3,k)
          call cross_product(dpoly_w(1:3,1,j,kk),aux1(1:3),
     .                           aux2(1:3))
          CurlE(1:3,k) = CurlE(1:3,k) + poly_v(1,j)*aux2(1:3)
          call cross_product(dshapH(1:3,2),dshapH(1:3,1), aux1(1:3))
          CurlE(1:3,k) = CurlE(1:3,k) +
     .                   2.d0*poly_v(1,j)*poly_w(1,j,kk)*aux1(1:3)
c
c  .......loop over index j
          do j=1,ijkk-kk-2
            i=ijkk-kk-j
c
c  .........TYPE III...................................................
            k=k+1
ccc            write(*,*) 'TYPE III: k, i,j,kk = ',k,i,j,kk
            ShapE(1:3,k) = dpoly_u(1:3,i)*poly_v(i,j)*poly_w(i,j,kk)
            call cross_product(dpoly_u(1:3,i),dpoly_v(1:3,i,j),
     .                         aux1(1:3))
            call cross_product(dpoly_u(1:3,i),dpoly_w(1:3,i,j,kk),
     .                         aux2(1:3))
            CurlE(1:3,k) = -poly_w(i,j,kk)*aux1(1:3)
     .                     -poly_v(i,j)*aux2(1:3)
c
c  .........TYPE IV....................................................
            k=k+1
ccc            write(*,*) 'TYPE IV: k, i,j,kk = ',k,i,j,kk
            ShapE(1:3,k) = (dpoly_u(1:3,i)*poly_v(i,j) +
     .                      poly_u(i)*dpoly_v(1:3,i,j))*poly_w(i,j,kk)
            call cross_product(dpoly_u(1:3,i),dpoly_w(1:3,i,j,kk),
     .                         aux1(1:3))
            call cross_product(dpoly_v(1:3,i,j),dpoly_w(1:3,i,j,kk),
     .                         aux2(1:3))
            CurlE(1:3,k) = - poly_v(i,j)*aux1(1:3) - poly_u(i)*aux2(1:3)
c
c  .......end of loop over index j
          enddo
c  .....end of loop over index k
        enddo
c
c
c  .....TYPE I (tet of the first type !!)..............................
        if (ijkk.lt.nordp+1) then
c  .......loop over over index k
          do kk=1,ijkk-3
c  .........loop over index j
            do j=1,ijkk-kk-2
              i=ijkk-kk-j
              k=k+1
ccc              write(*,*) 'TYPE I: k, i,j,kk = ',k,i,j,kk
              ShapE(1:3,k) = dpoly_u(1:3,i)*poly_v(i,j)*poly_w(i,j,kk)+
     .                       poly_u(i)*dpoly_v(1:3,i,j)*poly_w(i,j,kk)+
     .                       poly_u(i)*poly_v(i,j)*dpoly_w(1:3,i,j,kk)
c  ------->   CurlE(1:3,k) = 0
c  .........end of loop over index k
            enddo
c  .......end of loop over index j
          enddo
        endif
c
c  ...end of loop over oder of approximation
      enddo
c
c
c  ...save total number of dof's
 999  NrdofE = k
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7300) Xi(1:3)
 7300   format(' shape3nE:   Xi = ',3(e12.5,2x))
        do n=1,NrdofE
          write(*,7100) n,ShapE(1:3,n),CurlE(1:3,n)
 7100     format(' n = ',i3,'; ShapE,CurlE = ',2(3(e12.5,2x),2x))
        enddo
ccc        call pause
      endif
c
c
      end subroutine
