c----------------------------------------------------------------------
c     
c     routine name       - shapeEt_0
c     
c----------------------------------------------------------------------
c     
c     latest revision    - May 12
c     
c     purpose            - routine evaluates H(curl) shape functions
c                          for a triangle, using construction
c                          of Sabine Zaglmayr
c     
c     arguments :
c     in:
c       Xi        - master element coordinates
c       Norder    - order of approximation for the nodes
c       Norient   - orienttaion of edges
c     out:
c       NrdofE    - number of element dof
c       ShapE     - values of shape functions
c       CurlE     - values of derivatives of the shape functions
c                   wrt to master element coordinates      
c     
c----------------------------------------------------------------------
c
      subroutine shapeEt_0(Xi,Norder,Norient, NrdofE,ShapE,CurlE)
c     
      use element_data
      use parameters
      use cross_product_module
#include "syscom.blk"
c     
c  ...input output parameters
      dimension Xi(2),Norder(4),Norient(3),
     .  ShapE(2,MAXtriaE),CurlE(MAXtriaE)
c
      integer, parameter :: ndelta = 2
      dimension dp_dXi(2),dq_dXi(2)
      dimension vshap(3),dshap(2,3)
      dimension u(0:MAXP+ndelta), du(2,0:MAXP+ndelta)
c  ...cross checking with tetr
      dimension xi_tetr(3),
     .  norder_tetr(19), ne_orient_tetr(12), nf_orient_tetr(6),
     .  shapE_tetr(3,MAXbrickE),curlE_tetr(3,MAXbrickE)
c     
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...compute the affine coordinates       
      vshap(1) = 1.d0 -  (Xi(1)+Xi(2))
      dshap(1,1) = -1.d0
      dshap(2,1) = -1.d0
c
      vshap(2) = Xi(1)
      dshap(1,2) = 1.d0
      dshap(2,2) = 0.d0
c
      vshap(3) = Xi(2)
      dshap(1,3) = 0.d0
      dshap(2,3) = 1.d0
c
c  ...shape function counter     
      k = 0
c
c  ...mid-edge nodes shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do ie=1,3
c  .....get the edge vertices wrt the local edge orientation
        iv1 =  TRIAN_EDGE_TO_VERT(1,ie); iv2 =  TRIAN_EDGE_TO_VERT(2,ie)
c
        nord = Norder(ie)
        if (nord.gt.0) then
          if (Norient(ie).eq.1) call swap(iv1, iv2)
c
c  .......Whitney function          
          k = k + 1
          ShapE(1:2,k) =
     .      dshap(1:2,iv1)*vshap(iv2) - vshap(iv1)*dshap(1:2,iv2)
          call cross_product2D(dshap(1:2,iv2),dshap(1:2,iv1),CurlE(k))
          CurlE(k) = 2.0*CurlE(k)
c     
c  .......High order
          p      = vshap(     iv1) - vshap(     iv2)
          dp_dXi = dshap(1:2, iv1) - dshap(1:2, iv2)
          q      = vshap(     iv1) + vshap(     iv2)
          dq_dXi = dshap(1:2, iv1) + dshap(1:2, iv2)

          call Legendre_S_I(p,q,nord+1, u,du)
          do j=1,(nord-1)
            k = k + 1
            ShapE(1:2,k) = du(1,j+1)*dp_dXi(1:2) + du(2,j+1)*dq_dXi(1:2)
            CurlE(k) = 0.d0
          enddo
        end if
      enddo
      
c     
c  ...middle node shape functions
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call shapEbt(Xi,Norder(4),0,vshap,dshap,
     .  nrdofEt, ShapE(1:2,k+1),CurlE(k+1))
c
      NrdofE = k + nrdofEt
c
c     ...check with tetr
      if (iprint.eq.1) then
        xi_tetr = 0.d0; xi_tetr(1:2) = Xi(1:2)
c        
        norder_tetr = 0;
        norder_tetr(1:3) = Norder(1:3); norder_tetr(7) = Norder(4);
c
        ne_orient_tetr = 0; nf_orient_tetr = 0;
        ne_orient_tetr(1:3) = Norient(1:3)
c        
        nrdofE_tetr = 0; shapE_tetr = 0.d0; curlE_tetr = 0.d0;
c
        call shape3nE(
     .    xi_tetr,norder_tetr,ne_orient_tetr,nf_orient_tetr,
     .    nrdofE_tetr,shapE_tetr,curlE_tetr)
c
        write(*,7001) Xi,Norder(1:4),Norient(1:3)
 7001   format('shapEbe: Xi,Nord,Norient = ',2(f8.3,2x),4i3,2x,3i3)
c
        write(*,*) 'check dofE = ', nrdofE_tetr, NrdofE
        do k=1,nrdofE_tetr
          do i=1,2
            write(*,7002) k,i, shapE_tetr(i,k), ShapE(i,k),
     .        abs(shapE_tetr(i,k) - ShapE(i,k))
          enddo
          write(*,*) 'curl   = ', curlE_tetr(3,i), CurlE(i),
     .      abs(curlE_tetr(3,i) - CurlE(i))
        enddo

 7002   format('k,i = ',2i2,' Shap3nE = ',e12.5,
     .    ' Shape2tE = ',e12.5,
     .    ' Diff = ',e12.5)
        call pause
      endif
c      
      end subroutine
