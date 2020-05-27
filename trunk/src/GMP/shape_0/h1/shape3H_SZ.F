c----------------------------------------------------------------------
c     
c     routine name       - shapHbt_SZ
c     
c----------------------------------------------------------------------
c     
c     latest revision    - Jan 11
c     
c     purpose            - evalute triangle H1 bubble shape functions
c                          using Sabine Zaglmayr construction
c     
c     arguments :
c     in:
c       Xi        - master element coordinates
c       Nord      - order of approximation of 'mdlt'
c       Nort      - node orientation (0-5)
c       Vshap     - affine coordinates (linear vertex shape functions)
c       Dvshap    - derivatives of vertex shape functions
c     out:
c       NrdofH    - number of element dof
c       ShapH     - values of shape functions
c       DshapH    - values of derivatives of the shape functions
c                   wrt to master element coordinates      
c     
c----------------------------------------------------------------------
c     
      subroutine shapHbt_SZ(Xi,Nord,Nort,Vshap,Dvshap,
     .                      NrdofH,ShapH,DshapH)
c
      use element_data
      use parameters
#include "syscom.blk"
c
      dimension Xi(2), Vshap(3),Dvshap(2,3),
     .          ShapH(MAXmdltH),DshapH(2,MAXmdltH)
c
c  ...allow for exceeding the order of approximation  
      integer, parameter :: ndelta = 2
c
      integer :: jtemp(1:3) = (/1,2,3/),iv(1:3)
c
      dimension uu(0:MAXP+ndelta),duu(2,0:MAXP+ndelta),
     .          duu_dXi(2,0:MAXP+ndelta)
      dimension vv(0:MAXP+ndelta),dvv(  0:MAXP+ndelta),
     .          dvv_dXi(2,0:MAXP+ndelta)
c
c  ...derivatives of projections
      dimension dp_dXi(2), dq_dXi(2), ds_dXi(2)
c
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...account for orientation
      iv(1:3) = jtemp(TRIAN_L2G(1:3,Nort))
c     
c  ...compute Scaled Integrated Legendre polynomials
      p      =  Vshap(    iv(1)) -  Vshap(    iv(2))
      dp_dXi = Dvshap(1:2,iv(1)) - Dvshap(1:2,iv(2))
      q      =  Vshap(    iv(1)) +  Vshap(    iv(2))
      dq_dXi = Dvshap(1:2,iv(1)) + Dvshap(1:2,iv(2))
      call Legendre_S_I(p,q,Nord, uu,duu)
c     
c  ...compute Legendre polynomials
      s      = Vshap (    iv(3)) - Vshap (    iv(1)) - Vshap(     iv(2))
      ds_dXi = Dvshap(1:2,iv(3)) - Dvshap(1:2,iv(1)) - Dvshap(1:2,iv(2))
      call Legendre(s,Nord, vv,dvv)
c
      do i=0,Nord
        duu_dXi(1:2,i) = duu(1,i)*dp_dXi(1:2) + duu(2,i)*dq_dXi(1:2)
        dvv_dXi(1:2,i) = dvv(i)*ds_dXi(1:2)
      enddo
c
c  ...initiate shape function counter
      k=0
c
c  ...loop through the total order of u*v
      do ij=2,Nord-1
c
c ......loop through order of i
        do i=2,ij-0
          j=ij-i
          k=k+1
          ShapH(k) = uu(i)*vv(j)*Vshap(3)
          DshapH(1:2,k) = duu_dXi(1:2,i)*vv(j)*Vshap(3)
     .                  + uu(i)*dvv_dXi(1:2,j)*Vshap(3)
     .                  + uu(i)*vv(j)*Dvshap(1:2,3)
        enddo
      enddo
c     
c  ...save total number of dof's      
 999  NrdofH = k
c     
c  ...printing
      if (iprint.eq.1) then
        write(*,7010) Xi(1:2)
 7010   format('shapHbt_SZ:   Xi = ',2(e12.5,2x))       
        do n=1,NrdofH
          write(*,7020) n,ShapH(n),DshapH(1:2,n)
 7020     format(' n = ',i3,'; ShapH,DshapH = ',e12.5,3x,2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine
