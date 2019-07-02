c----------------------------------------------------------------------
c
c   routine name       - curve_OneSur
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 06
c
c   purpose            - routine defines parametrization for a curve
c                        (geometry reconstruction of Cynthia, needs
c                         an explanation)
c
c   arguments :
c     in:
c               No     - the curve number
c               Xi     - reference coordinate  (between 0 and 1)
c     out:
c               X      - physical coordinates of the point
c               Dxdxi  - derivatives of the physical coordinates wrt
c                        reference coordinate
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine curve_OneSur(No,Xi, X,Dxdxi)
c
      use GMP
      use control
#include "syscom.blk"
c
      dimension X(3),Dxdxi(3)
c
      dimension xv(3,2), h(4),dh(4),a1(3),a2(3)
      dimension aij(3,3),xprim(3),dxprdxi(3),center(NDIM)
      dimension radv1(3),radv2(3),direc(3),nrep(2)
c
c  ...implicit curve:
c       nsurf     - surfaces defining the curve
c       sfact     - sign factors to adjust the orientation of surfaces
c       xmid      - midpoint
      dimension Xs(3),nsurf(4),aux(3),ntem(5),dd(3),Xif(2,2),
     .          sfact(4),xmid(3)
      dimension ver1(3),ver2(3),x1(3),x2(3),x3(3),rn(3)
      dimension der1(3),der2(3),der3(3),der4(3),aalpha(3,3),bbeta(3)
      dimension xn1(3),xn2(3),b1(3),b2(3),xnc(3)
c
      character*9 Type
      data eps /10e-8/
c
cc      if (no.eq.1) then
cc        iprint=0
cc      else
cc        iprint = 0
cc      endif
      iprint_mnewt = iprint
c
      if ((CURVES(No)%Type.ne.'OneSur').or.NDIM.ne.3)) then
        write(*,7001)
 7001   format('curve_ImpCir: WRONG CALL')
        stop 1
      endif
c
c  ...get  endpoints for the curve
      do iend = 1,2
        np = CURVES(No)%EndPoNo(iend)
        call pointr(np, xv(1:3,iend))
      enddo
c
      nsurf = CURVES(No)%Idata(1)
      do j=1,3
        X(j) = Xi*xv(j,2) + (1.d0-Xi)*xv(j,1)
      enddo
c
c  ...find the surface value
      do i=1,10
c        .....surface 1
              call surf(nsurf,X, Fval,der1)
              bbeta(1) = -Fval
c        .... surface 2
              nf = CURVES(No)%Idata(2)
              do ii = 1, 3
                 nc = RECTANGLES(nf)%EdgeNo(ii)
                 if(nc.gt.0)then
                    ntem(ii) = CURVES(nc)%EndPoNo(1)
                 else
                    ntem(ii) = CURVES(-nc)%EndPoNo(2)
                 endif
              enddo
              X1(1:3) = POINTS(ntem(1))%Rdata(1:3)
              X2(1:3) = POINTS(ntem(2))%Rdata(1:3)
              X3(1:3) = POINTS(ntem(3))%Rdata(1:3)
              call plane2(X,X1,X2,X3, Fval,der2)
c        .....surface 3 and surface 4
              X1(1:3) = xv(1:3,1)
              data X2/1,0,0/
              data X3/0,1,0/
              call plane2(X,X1,X2,X3, Fval3,der3)
              X1(1:3) = xv(1:3,2)
              call plane2(X,X1,X2,X3, Fval4,der4)
c
              bbeta(3) = -(1.d0-Xi)*Fval3 - Xi*Fval4
              do  ii=1,3
                 aalpha(1,ii) = der1(ii)
                 aalpha(2,ii) = der2(ii)
                 aalpha(3,ii) = (1.d0-Xi)*der3(ii) + Xi*der4(ii)
              enddo
c
              errf=0
              do  ii=1 ,3
                 errf=errf+abs(bbeta(ii))
              enddo
              if (errf.le.eps)return
c
              call saruss(aalpha,bbeta, aux)
              bbeta(1) = aux(1)
              bbeta(2) = aux(2)
              bbeta(3) = aux(3)
c
              errx=0
              do ii=1,3
                 errx=errx+abs(bbeta(ii))
                 X(ii)=X(ii)+bbeta(ii)
              enddo
c
           enddo
c
      end


