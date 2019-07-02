c----------------------------------------------------------------------
c
c   routine name       - reconstruct_normal
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 07
c
c   purpose            - routine reconstructs a normal at a surface
c                        point by performing the least squares fit
c                        of a quadratic surface
c
c   arguments
c     in:
c          Np          - a point on the surface
c          List_curves - list of surface curves connected to the point
c          Nrcurves    - number of curves on the list
c     out:
c          Rnorm       - components of the normal
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_normal(Np,List_curves,Nrcurves, Rn)
c
c
      use GMP
      use control
c
#include "syscom.blk"
c
      dimension List_curves(Nrcurves),Rn(3)
c
c  ...the reconstructed hesjan
      dimension acoeff(3,3)
c
c  ...points to fit
      dimension xp(3,20),weight(20)
c
c  ...matrices for a system of linear equations
      dimension aa(6,6),bb(6,3),cc(3,6),dd(3,3),uu(6),aux(6,3)

c  ...work space
      dimension eigen(3),work(37),iwork(18),rnw(3),aGauss(3,3),temp(3)
c
c
c**********************************************************************
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,7006) Nrcurves,Rn
 7006   format(' reconstruct_normal: Nrcurves,Rn = ',i2,2x,3f8.3)
      endif
c
c  ...collect point coordinates in a system with origin shifted to "Np"
      dmin = 1.d30
      do ic=1,Nrcurves
        nc = List_curves(ic)
        do j=1,2
          np1 = CURVES(nc)%EndPoNo(j)
          if (np1.ne.Np) then
            xp(1:3,ic) = POINTS(np1)%Rdata(1:3) - POINTS(Np)%Rdata(1:3)
          endif
        enddo
        call norm(xp(1:3,ic), weight(ic))
        dmin = min(dmin,weight(ic))
      enddo
      do ic=1,Nrcurves
        weight(ic) = (dmin/weight(ic))**.7d0
ccc        weight(ic) = 1.d0
      enddo
      if (iprint.eq.1) then
        write(*,*) 'reconstruct_normal: weight = '
        write(*,7035) weight(1:Nrcurves)
 7035   format(10f8.3)
      endif
c
c  ...set up the matrices for the linear system corresponding to the
c     least squares fit problem
      aa(1:6,1:6) = 0.d0; bb(1:6,1:3) = 0.d0
      jj=0
      do i=1,3
      do j=i,3
        if (i.eq.j) then
          ddd=.5d0
        else
          ddd=1.d0
        endif
        jj=jj+1
c
c  .....loop through the points
        do ic=1,Nrcurves
          do k=1,3
            bb(jj,k) = bb(jj,k) + ddd*xp(i,ic)*xp(j,ic)*xp(k,ic)
     .                *weight(ic)
          enddo
          ii=0
          do k=1,3
          do l=k,3
            ii=ii+1
            if (k.eq.l) then
              ccc=.5d0
            else
              ccc=1.d0
            endif
            aa(jj,ii) = aa(jj,ii)
     .                + ccc*xp(k,ic)*xp(l,ic)*ddd*xp(i,ic)*xp(j,ic)
     .                 *weight(ic)
          enddo
          enddo
        enddo
      enddo
      enddo
c
c  ...compute the maximum diagonal entry
      scale = 0.d0
      do jj=1,6
        scale = max(scale,abs(aa(jj,jj)))
      enddo
c
c  ...add the regularizing term of two magnitudes smaller
      jj=0
      do i=1,3
      do j=i,3
        jj=jj+1
        if (i.eq.j) then
          ddd = .01d0
        else
          ddd = .02d0
        endif
        aa(jj,jj) = aa(jj,jj) + scale*ddd
      enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,7001)
 7001   format(' reconstruct_normal: aa,bb = ')
        do jj=1,6
          write(*,7002) jj,aa(jj,1:6),bb(jj,1:3)
 7002     format(1x,i1,3x,6e12.5,5x,3e12.5)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
      cc(1:3,1:6) = 0.d0; dd(1:3,1:3) = 0.d0
      do j=1,3
c
c  .....loop through the points
        do ic=1,Nrcurves
          do k=1,3
            dd(j,k) = dd(j,k) + xp(j,ic)*xp(k,ic)*weight(ic)
          enddo
          ii=0
          do k=1,3
          do l=k,3
            ii=ii+1
            if (k.eq.l) then
              ccc=.5d0
            else
              ccc=1.d0
            endif
            cc(j,ii) = cc(j,ii)
     .               + ccc*xp(k,ic)*xp(l,ic)*xp(j,ic)*weight(ic)
          enddo
          enddo
        enddo
      enddo
c
c  ...check symmetry
      do i=1,6
      do j=1,3
        diff = abs(bb(i,j)-cc(j,i))
        if (diff.gt.1.d-12) then
          write(*,*)'i,j,diff = ',i,j,diff
          stop 1
        endif
      enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...set up and solve the eigenvalue problem
      call tri(aa,6,6)
      do k=1,3
        call rhsub(aa,aux(1:6,k),bb(1:6,k),6,6)
      enddo
      do k=1,3
        do l=1,3
          s = 0.d0
          do i=1,6
            s = s + cc(k,i)*aux(i,l)
          enddo
          dd(k,l) = dd(k,l) - s
        enddo
      enddo
      if (iprint.eq.1) then
        do k=1,3
          write(*,7008) k,dd(k,1:3)
 7008     format(' k,dd = ',i2,3e12.5)
        enddo
      endif
      call DSYEVD('V','U',3,dd,3,eigen,work,37,iwork,18,info)
      if (iprint.eq.1) then
        write(*,7012) eigen
 7012   format(' reconstruct_normal: eigen = ',3e12.5)
      endif
c
      if (info.ne.0) then
        write(*,*)'reconstruct_normal: info = ',info
        stop 1
      endif
c
c-----------------------------------------------------------------------
c
c  ...pick the eigenvector that minimizes the misfit function
      rnw(1:3) = dd(1:3,1)
      call norm(rnw-Rn, d1)
      call norm(-rnw-Rn, d2)
      if (d1.lt.d2) then
        Rn(1:3) = rnw(1:3)
      else
        Rn(1:3) = -rnw(1:3)
      endif
c
c
      end

