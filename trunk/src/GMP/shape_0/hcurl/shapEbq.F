c----------------------------------------------------------------------
c
c     routine name       - shapEbq
c
c----------------------------------------------------------------------
c
c     latest revision    - Feb 11
c
c     purpose            - routine evaluates H(curl) bubbles
c                          for a quad, using construction
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
      subroutine shapEbq(Xi,Nord,Nort,
     .                   NrdofE,ShapE,CurlE)
c
      use element_data
      use parameters
#include "syscom.blk"
c
c  ...input output parameters
      dimension Xi(2),
     .          ShapE(2,MAXmdlqE),CurlE(MAXmdlqE)
c
c  ...allow for exceeding the order of approximation
      integer, parameter :: ndelta = 2

      integer :: jtemp(1:4) = (/1,2,3,4/), iv(1:4)
c
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...paolo, Mar 11
      ShapE = 0.d0 ; CurlE = 0.d0
      write(*,*)'shapEbq: return 0 for now'
      return
c
ccc      vsig(1) = 1.d0 - Xi(1) + 1.d0 - Xi(2)
ccc      vsig(2) =        Xi(1) + 1.d0 - Xi(2)
ccc      vsig(3) =        Xi(1) +        Xi(2)
ccc      vsig(4) = 1.d0 - Xi(1) +        Xi(2)
ccc
ccc      iv(1:4) = jtemp(QUADR_L2G(1:4,Nort))
ccc      call decode(Nord, nordh,nordv)
cccc
cccc  ...compute Integrated legendre functions
ccc      p      =  vsig(    iv(1)) -  vsig(    iv(2))
ccc      dp_dXi = dvsig(1:2,iv(1)) - dvsig(1:2,iv(2))
ccc      call legendre_I(p, nordh, uu, duu)
ccc      do j=1,nordh
ccc        duu_dXi(1:2,j) = duu(j)*dp_dXi(1:2)
ccc      enddo
cccc
ccc      q      =  vsig(    iv(1)) -  vsig(    iv(4))
ccc      dq_dXi = dvsig(1:2,iv(1)) - dvsig(1:2,iv(4))
ccc      call legendre_I(q, nordv, vv, dvv)
ccc      do j=1,nordv
ccc        dvv_dXi(1:2,j) = dvv(j)*dq_dXi(1:2)
ccc      enddo
cccc
cccc  ...horizontal first
ccc      do j=2,nordv
ccc        k=k+1
ccc        ShapE(1:3,k) = duu_dXi(1:3,1)*ww(j)
ccc        call cross_product(dww_dXi(1:3,j),duu_dXi(1:3,1),aux_a(1:3,2))
ccc        CurlE(1:3,k) = 2.d0*aux_a(1:3,1)*ww(j)+aux_a(1:3,2)
ccc        do i=2,nordh
ccc          k=k+1
ccc          ShapE(1:3,k) = duu_dXi(1:3,i)*ww(j)
ccc          call cross_product(dww_dXi(1:3,j), duu_dXi(1:3,i),
ccc     .                       CurlE(1:3,k))
ccc        enddo
ccc      enddo
cccc  ...vertical next
ccc      do j=1,nordv
ccc        do i=2,nordh
ccc          k=k+1
ccc          ShapE(1:3,k) = uu(i)*dww_dXi(1:3,j)
ccc          call cross_product(duu_dXi(1:3,i), dww_dXi(1:3,j),
ccc     .                       CurlE(1:3,k))
ccc        enddo
ccc      enddo
      end subroutine

