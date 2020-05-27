c----------------------------------------------------------------------
c
c     routine name       - shapEbt
c
c----------------------------------------------------------------------
c
c     latest revision    - Jan 11
c
c     purpose            - routine evaluates H(curl) bubbles
c                          for a triangle, using construction
c                          of Sabine Zaglmayr
c
c     arguments :
c     in:
c       Xi        - master element coordinates
c       Nord      - order of approximation for 'mdlt'
c       Nort      - orienttaion of 'mdlt'
c       Vshap     - affine coordinates (vertex H1 shape functions)
c       Dvshap    - derivatives of Vshap
c     out:
c       NrdofE    - number of element dof
c       ShapE     - values of shape functions
c       CurlE     - values of derivatives of the shape functions
c                   wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapEbt(Xi,Nord,Nort,Vshap,Dvshap,
     .                      NrdofE,ShapE,CurlE)
c
      use element_data
      use parameters
      use cross_product_module
#include "syscom.blk"
c
c  ...input output parameters
      dimension Xi(2),Vshap(3),Dvshap(2,3),
     .          ShapE(2,MAXtriaE),CurlE(MAXtriaE)
c
c  ...allow for exceeding the order of approximation
      integer, parameter :: ndelta = 2
c
      dimension dp_dXi(2), dq_dXi(2), ds_dXi(2)
      dimension uu(0:MAXP+ndelta),duu(2,0:MAXP+ndelta),
     .          duu_dXi(2,0:MAXP+ndelta)
      dimension vv(0:MAXP+ndelta),dvv(  0:MAXP+ndelta),
     .          dvv_dXi(2,0:MAXP+ndelta)
c
      integer :: jtemp(1:3) = (/1,2,3/), iv(1:3)
      dimension aux2(2)
c
c-----------------------------------------------------------------------
c
      iprint=0
      if (Nord.eq.1) then
        NrdofE=0; return
      endif
c
      iv(1:3) = jtemp(TRIAN_L2G(1:3,Nort))
c
c  ...compute Scaled Integrated Legendre polynomials
      p      = Vshap(     iv(1)) - Vshap(     iv(2))
      dp_dXi = Dvshap(1:2,iv(1)) - Dvshap(1:2,iv(2))
      q      = Vshap(     iv(1)) + Vshap(     iv(2))
      dq_dXi = Dvshap(1:2,iv(1)) + Dvshap(1:2,iv(2))
      call Legendre_S_I(p,q,Nord, uu,duu)
c
c  ...compute Scaled Legendre polynomials
      s      = -Vshap(     iv(1)) -Vshap(     iv(2)) +Vshap(     iv(3))
      ds_dXi = -Dvshap(1:2,iv(1)) -Dvshap(1:2,iv(2)) +Dvshap(1:2,iv(3))
      call Legendre(s,Nord, vv,dvv)
c
c  ...compute block vv, \nabla uu, \nabla vv
      do i=0,Nord
        duu_dXi(1:2,i) =duu(1,i)*dp_dXi(1:2) + duu(2,i)*dq_dXi(1:2)
        dvv_dXi(1:2,i) = Dvshap(1:2,iv(3))*vv(i)
     .                 + Vshap(iv(3))*dvv(i)*ds_dXi(1:2)
        vv(i)          = Vshap(iv(3))*vv(i)
      enddo
c
      k=0
c
c  ...loop through the total order first (ij+2)
      do ij=0,Nord-2
        j=ij
c
c  .....COMPLETION with TYPE III....................................
        k=k+1
        ShapE(1:2,k) = (Dvshap(1:2,iv(1))*Vshap(iv(2)) -
     .                  Dvshap(1:2,iv(2))*Vshap(iv(1)))*vv(j)
c
c===========================================================================
c  ...curl TYPE III =
c     \nabla v x (\lambda_2 \nabla \lambda_1 - \lambda_1 \nabla \lambda_2) +
c            2 v (\nabla \lambda_2 x \nabla \lambda_1) =
c     \nabla v x aux2 + 2v * aux1
c===========================================================================
c
c       aux2 = (\lambda_2 \nabla \lambda_1 - \lambda_1 \nabla \lambda_2)
        aux2(1:2) = Vshap(iv(2))*Dvshap(1:2,iv(1)) -
     .              Vshap(iv(1))*Dvshap(1:2,iv(2))
c       aux1 = \nabla \lamda_2 x \nabla \lambda_1
        call cross_product2D(Dvshap(1:2,iv(2)),Dvshap(1:2,iv(1)), aux1)
c       \nabla v x aux1
        call cross_product2D(dvv_dXi(1:2,j),aux2(1:2), CurlE(k))
c       \nabla v x aux2 + 2v * aux1
        CurlE(k) = CurlE(k) + 2.d0*vv(j)*aux1
c
c  .....COMPLETION with TYPE II....................................
        do i=0,ij
          j=ij-i
          k=k+1
          ShapE(1:2,k) = duu_dXi(1:2,i+2)* vv(        j) -
     .                    uu(        i+2)*dvv_dXi(1:2,j)
c
c  .......curl TYPE II = 2 \nabla v x \nabla u
          call cross_product2D(dvv_dXi(1:2,j),
     .                         duu_dXi(1:2,i+2), CurlE(k))
          CurlE(k) = 2.d0*CurlE(k)
        enddo
c
c  .....COMPLETION with TYPE I (gradients).........................
        do i=0,ij
          j=ij-i
c
c  .......(triangle of the first type !)
          if (ij.ne.Nord-2) then
            k=k+1
            ShapE(1:2,k) = duu_dXi(1:2,i+2)* vv(        j) +
     .                      uu(        i+2)*dvv_dXi(1:2,j)
c  -------> CurlE(k) = 0
          endif
        enddo
c
c  ...end of loop over total order
      enddo
c
c  ...save total number of dof's
 999  NrdofE = k
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7010) Xi(1:2)
 7010   format('shapEbt:   Xi = ',3(e12.5,2x))
        do n=1,NrdofE
          write(*,7020) n,ShapE(1:2,n),CurlE(n)
 7020     format(' n = ',i3,'; ShapE,CurlE = ',2(2(e12.5,2x),2x))
        enddo
        call pause
      endif
c
      end subroutine
