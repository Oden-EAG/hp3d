      subroutine eval_vtk(X,Xi,Dxdxi,Dxidx,Rjac,
     .                    NrdofH,ShapH,GradH,ZdofH, ZH,ZgradH)
c
      use parameters
      use physics
c
#include "syscom.blk"
c
      dimension X(3),Xi(3),Dxdxi(3,3),Dxidx(3,3)
      dimension ShapH(MAXbrickH),GradH(3,MAXbrickH)
      dimension ZdofH(NRHVAR,MAXbrickH),ZH(NRHVAR),ZgradH(NRHVAR,3)
c
      dimension zgradhH(3)
c
c-----------------------------------------------------------------------
c
c     H^1 variables and gradient
      do i=1,NRHVAR
c
        ZH(i) = ZERO; zgradhH(:) = ZERO
        do k=1,nrdofH
          ZH(i)      = ZH(i)      + ZdofH(i,k)*ShapH(k)
          zgradhH(:) = zgradhH(:) + ZdofH(i,k)*GradH(:,k)
        enddo
c
        do k=1,3
          ZgradH(i,k) = Dxidx(1,k)*zgradhH(1)
     .                + Dxidx(2,k)*zgradhH(2)
     .                + Dxdxi(3,k)*zgradhH(3)
        enddo
c
      enddo
c
      end subroutine eval_vtk
