c----------------------------------------------------------------------
c
c     routine name       - shape3pE
c
c----------------------------------------------------------------------
c
c     latest revision    - Jan 11
c
c     purpose            - 3D Nedelec's prism of the 1st type,
c                          routine evaluates H(curl) shape functions
c                          and derivatives using Sabine's construction
c                          based on scaled integrated Legendre
c                          polynomials and integrated Legendre
c                          polynomials
c
c     arguments :
c     in:
c       Xi        - master element coordinates
c       Norder    - order of approximation for the nodes
c       Ne_orient - edge orientations
c       Nf_orient - face orientations
c     out:
c       NrdofE    - number of element dof
c       ShapE     - values of shape functions
c       CurlE     - values of derivatives of the shape functions
c                   wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shape3pE(Xi,Norder,Ne_orient,Nf_orient,
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
c  ...``affine'' coordinates
      dimension vlambda(6),dlambda(1:3,6)
      dimension vmu(6),dmu(1:3,6), imu(6)
c
c  ...vertices on an edge or face
      dimension iv(4), jtemp(4)
      dimension aux_a(3,2),aux_b(3,2)
      dimension dp_dXi(3), dq_dXi(3), ds_dXi(3), dt_dXi(3)
c
      dimension uu(0:MAXP+ndelta),duu(2,0:MAXP+ndelta)
      dimension vv(0:MAXP+ndelta),dvv(  0:MAXP+ndelta)
      dimension ww(0:MAXP+ndelta),dww(  0:MAXP+ndelta)
c
c  ...chain rule
      dimension duu_dXi(3,0:MAXP+ndelta)
      dimension dvv_dXi(3,0:MAXP+ndelta)
      dimension dww_dXi(3,0:MAXP+ndelta)
c
c  ...work space for shapHbt_SZ
      dimension shap2H(MAXmdltH),dshap2H(3,MAXmdltH)
c
c  ...work space for shapEbt
      dimension shap2E(3,MAXtriaE),curl2E(3,MAXtriaE)
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...initialization
      ShapE(1:3,1:MAXbrickE) = 0.d0
      CurlE(1:3,1:MAXbrickE) = 0.d0
c
c  ...check if the input point is within the master prism
      if (((Xi(1)+Xi(2)).gt.1.d0+GEOM_TOL).or.
     .     (Xi(1).lt.-GEOM_TOL).or.
     .     (Xi(2).lt.-GEOM_TOL).or.
     .     (Xi(3).gt.1.d0+GEOM_TOL).or.
     .     (Xi(3).lt.-GEOM_TOL)) then
        write(*,7001)
 7001   format(' shape3pE: WARNING !! POINT Xi OUTSIDE OF MASTER PRISM')
        write(*,7002) Xi(1:3)
 7002   format(' shape3pE: Xi = ',3(e12.5,2x))
        call pause
      endif
c
c  ...calculate the affine coordinates and their gradients
      vlambda(1)     =  1.d0 - Xi(1) - Xi(2)
      dlambda(1:2,1) = -1.d0;
      dlambda(3  ,1) =  0.d0
c
      do i=1,2
        vlambda(    i+1) = Xi(i)
        dlambda(1:3,i+1) = 0.d0
        dlambda(i,  i+1) = 1.d0
      enddo
c
      vlambda(    4:6) = vlambda(    1:3)
      dlambda(1:3,4:6) = dlambda(1:3,1:3)
c
      imu(1:3) = 0;                 imu(4:6) = 1
      vmu(    1:3) =  1.d0 - Xi(3); vmu(    4:6) =  Xi(3)
      dmu(1:3,1:3) =  0.d0;         dmu(1:3,4:6) =  0.d0
      dmu(3,  1:3) = -1.d0;         dmu(3,  4:6) =  1.d0
c
c======================================================================
c
c  ...initiate shape functions counter
      k=0
c
c  ...horizontal mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  ...loop over horizontal edges
      do i=1,6
c
c  .....edge # of dof
        ndofE = Norder(i)
c
c  .....get vertices specifing local edge orientation
        iv(1:2) = PRISM_EDGE_TO_VERT(1:2,i)
c
c  .....determine endpoints in the global orientation
        if (Ne_orient(i).eq.1) call swap(iv(1),iv(2))
c
c  .....edge Whitney shape function
        k=k+1
        ShapE(1:3,k) = ( dlambda(1:3,iv(1))*vlambda(iv(2)) -
     .                   dlambda(1:3,iv(2))*vlambda(iv(1)) )*vmu(iv(1))
        call cross_product(
     .      dlambda(1:3,iv(2))*vmu(iv(1))+vlambda(iv(2))*dmu(1:3,iv(1)),
     .      dlambda(1:3,iv(1)),
     .      aux_a(1:3,1) )
        call cross_product(
     .      dlambda(1:3,iv(1))*vmu(iv(1))+vlambda(iv(1))*dmu(1:3,iv(1)),
     .      dlambda(1:3,iv(2)),
     .      aux_a(1:3,2) )
        CurlE(1:3,k) = aux_a(1:3,1) - aux_a(1:3,2)
c
c  .....higher order edge shape functions are gradients
        if (ndofE.ge.2) then
c
c  .......edge coordinates
          p      = vlambda(    iv(1)) - vlambda(    iv(2))
          dp_dXi = dlambda(1:3,iv(1)) - dlambda(1:3,iv(2))
          q      = vlambda(    iv(1)) + vlambda(    iv(2))
          dq_dXi = dlambda(1:3,iv(1)) + dlambda(1:3,iv(2))
c
c  .......compute Scaled Integrated Legendre polynomials:
          call Legendre_S_I(p,q,Norder(i), uu,duu)
c
c  .......compute shape functions and their curl
          do j=1,ndofE-1
            k=k+1
            ShapE(1:3,k) = (duu(1,j+1)*dp_dXi(1:3) +
     .        duu(2,j+1)*dq_dXi(1:3))*vmu(iv(1)) +
     .        uu(  j+1)*dmu(1:3,iv(1))
c  -------> CurlE(1:3,k) = 0
          enddo
        endif
c
c  ...end of loop over horizontal edges
      enddo
c
c  ...vertical mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop over edges
      do i=7,9
c
c  .....determine edge # of dof
        ndofE = Norder(i)
c
c  .......get vertices specifing local edge orientation
          iv(1:2) = PRISM_EDGE_TO_VERT(1:2,i)
c
c  .......determine endpoints in the global orientation
          if (Ne_orient(i).eq.1) call swap(iv(1),iv(2))
c
c  .......edge Whitney shape function
          k=k+1
          ShapE(1:3,k) = vlambda(iv(1))*dmu(1:3,iv(1))
          call cross_product(
     .      dlambda(1:3,iv(1)), dmu(1:3,iv(1)),
     .      CurlE(1:3,k))
c
        if (ndofE.ge.2) then
c
c  .......edge coordinates
          p      = vmu(    iv(1)) - vmu(     iv(2))
          dp_dXi = dmu(1:3,iv(1)) - dmu(1:3, iv(2))
c
c  .......compute Integrated Legendre polynomials:
          call Legendre_I(p,Norder(i), uu,duu(1,:))
c
c  .......higher order edge shape functions are gradients
          do j=1,ndofE-1
            k=k+1
            ShapE(1:3,k) = (duu(1,j+1)*dp_dXi(1:3))*vlambda(    iv(1)) +
     .                       uu(  j+1)             *dlambda(1:3,iv(1))
c  -------> CurlE(1:3,k) = 0
          enddo
        endif
c
c  ...end of loop over vertical edges
      enddo
c
c  ...triangle mid-face nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      shap2E(1:3,1:MAXtriaE)=0.d0 ; curl2E(1:3,1:MAXtriaE)=0.d0
c  ...loop over triangular faces
      do ii=1,2
c
c  .....determine 2D H(curl) bubbles
        call shapEbt(Xi(1:2),Norder(9+ii),Nf_orient(ii),
     .                  vlambda(1:3),dlambda(1:2,1:3),
     .                  nrdof2E,shap2E(1:2,1:MAXtriaE),
     .                          curl2E(  3,1:MAXtriaE) )
c
        do i=1,nrdof2E
          k=k+1
          ShapE(1:2,k) = shap2E(1:2,i)*vmu(2+ii)
          call cross_product( dmu(1:3,2+ii), shap2E(1:3,i),
     .                        CurlE(1:3,k) )
          CurlE(3,k)   = CurlE(3,k) + curl2E(3,i)*  vmu(2+ii)
        enddo


      enddo
c
c  ...quadrilateral mid-face nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c  ...loop over faces
      do ii=3,5
c
        call decode(Norder(9+ii), nordh,nordv)
c
c  .....get face vertices
        jtemp(1:4) = PRISM_FACE_TO_VERT(1:4,ii)
c
c  .....arrange vertices according to face global orientation
        iv(1:4) = jtemp(QUADR_L2G(1:4,Nf_orient(ii)))
c
c  .....set iv(1) iv(2) as horizontal edge
        if ( imu(iv(1)).ne.imu(iv(2)) ) then
          call swap(iv(2),iv(4))
          call swap(nordh, nordv)
        endif
c
c  .....compute Scaled Integrated Legendre polynomials
        p      = vlambda(    iv(1)) - vlambda(    iv(2))
        dp_dXi = dlambda(1:3,iv(1)) - dlambda(1:3,iv(2))
        q      = vlambda(    iv(1)) + vlambda(    iv(2))
        dq_dXi = dlambda(1:3,iv(1)) + dlambda(1:3,iv(2))
        call Legendre_S_I(p,q,nordh, uu,duu)
c
c  .....compute Legendre polynomials
        s      = vmu(    iv(1)) - vmu(    iv(4))
        ds_dXi = dmu(1:3,iv(1)) - dmu(1:3,iv(4))
        call Legendre_I(s,nordv, ww,dww)
c
c  .....the first horizontal shape function is the Whitney function,
c       NOT a gradient
        duu_dXi(1:3,1) = dlambda(1:3,iv(1))*vlambda(    iv(2)) -
     .                   vlambda(    iv(1))*dlambda(1:3,iv(2))
        call cross_product(dlambda(1:3,iv(2)),dlambda(1:3,iv(1)),
     .                     aux_a(1:3,1))
c
c  .....higher order guys are gradients
        do i=2,nordh
          duu_dXi(1:3,i) =duu(1,i)*dp_dXi(1:3) + duu(2,i)*dq_dXi(1:3)
        enddo
c
c  .....the first vertical shape function
        dww_dXi(1:3,1) = dmu(1:3,iv(1))
        call cross_product(dlambda(1:3,iv(1)),dmu(1:3,iv(1)),
     .                     aux_b(1:3,1))
c
c  .....higher order guys are gradients
        do j=2,nordv
          dww_dXi(1:3,j) =dww(j)*ds_dXi(1:3)
        enddo


        select case(Nf_orient(ii))
        case(0,2,5,7)

c
c  .....``horizontal'' shape functions first
        do j=2,nordv
          k=k+1
          ShapE(1:3,k) = duu_dXi(1:3,1)*ww(j)
          call cross_product(dww_dXi(1:3,j),duu_dXi(1:3,1),aux_a(1:3,2))
          CurlE(1:3,k) = 2.d0*aux_a(1:3,1)*ww(j)+aux_a(1:3,2)
          do i=2,nordh
            k=k+1
            ShapE(1:3,k) = duu_dXi(1:3,i)*ww(j)
            call cross_product(dww_dXi(1:3,j), duu_dXi(1:3,i),
     .                         CurlE(1:3,k))
          enddo
        enddo
c
c  .....``vertical'' shape functions next
        do j=1,nordv
          do i=2,nordh
            k=k+1
            ShapE(1:3,k) = uu(i)*dww_dXi(1:3,j)
            call cross_product(duu_dXi(1:3,i), dww_dXi(1:3,j),
     .                         CurlE(1:3,k))
          enddo
        enddo

        case(1,3,4,6)
c
c  .....``vertical'' shape functions next
        do j=1,nordv
          do i=2,nordh
            k=k+1
            ShapE(1:3,k) = uu(i)*dww_dXi(1:3,j)
            call cross_product(duu_dXi(1:3,i), dww_dXi(1:3,j),
     .                         CurlE(1:3,k))
          enddo
        enddo
c
c  .....``horizontal'' shape functions first
        do j=2,nordv
          k=k+1
          ShapE(1:3,k) = duu_dXi(1:3,1)*ww(j)
          call cross_product(dww_dXi(1:3,j),duu_dXi(1:3,1),aux_a(1:3,2))
          CurlE(1:3,k) = 2.d0*aux_a(1:3,1)*ww(j)+aux_a(1:3,2)
          do i=2,nordh
            k=k+1
            ShapE(1:3,k) = duu_dXi(1:3,i)*ww(j)
            call cross_product(dww_dXi(1:3,j), duu_dXi(1:3,i),
     .                         CurlE(1:3,k))
          enddo
        enddo

        end select

c
c  ...end of loop over faces
      enddo
c
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call decode(Norder(15), nordh,nordv)
      shap2E=0.d0 ; curl2E=0.d0 ; dshap2H=0.d0;

c
c  ...determine 2D H1 and H(curl) bubble shape functions
      call shapHbt_SZ(Xi(1:2),nordh,0,vlambda(1:3),dlambda(1:2,1:3),
     .                nrdof2H,shap2H(     1:MAXmdltH),
     .                        dshap2H(1:2,1:MAXmdltH) )
      call shapEbt(Xi(1:2),nordh,0,vlambda(1:3),dlambda(1:2,1:3),
     .                nrdof2E,shap2E(1:2,1:MAXtriaE),
     .                        curl2E(  3,1:MAXtriaE) )
c
c  ...determine integrated Legendre polynomials in z
      s      = vmu(    1) - vmu(    4)
      ds_dXi = dmu(1:3,1) - dmu(1:3,4)
      call Legendre_I(s,nordv, ww,dww)
      do i=0,nordv
        dww_dXi(1:3,i) =dww(i)*ds_dXi(1:3)
      enddo
c
c  ...modify for the whitney
ccc      dww_dXi(1:3,1) = dmu(1:3, 1)*vmu(     4) -
ccc     .                 vmu(     1)*dmu(1:3, 4)
c
c  ...``horizontal'' shape functions
      do j=2,nordv
c
c  .....loop through 2D H(curl) bubbles
        do i=1,nrdof2E
          k=k+1
          ShapE(1:2,k) = shap2E(1:2,i)*ww(j)
          call cross_product(dww_dXi(1:3,j),Shap2E(1:3,i),
     .                       CurlE(1:3,k) )
          CurlE(3,k) = CurlE(3,k)+curl2E(3,  i)*ww(j)
        enddo
      enddo
c
c  ...``vertical'' shape functions
      do j=1,nordv
c
c ......loop through 2D H1 bubbles
        do i=1,nrdof2H
          k=k+1
          ShapE(1:3,k) = shap2H(i)*dww_dXi(1:3,j)
          call cross_product(dshap2H(1:3,i), dww_dXi(1:3,j),
     .                       CurlE(1:3,k) )
        enddo
      enddo
c
c  ...save total number of dof's
 999  NrdofE = k
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7300) Xi(1:3)
 7300   format(' shape3pE:   Xi = ',3(e12.5,2x))
        do n=1,NrdofE
          write(*,7100) n,ShapE(1:3,n),CurlE(1:3,n)
 7100     format(' n = ',i3,'; ShapE,CurlE = ',2(3(e12.5,2x),2x))
        enddo
        call pause
      endif
c
      end subroutine
