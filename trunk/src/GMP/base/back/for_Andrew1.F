c----------------------------------------------------------------------
c
c   routine name       - reconstruct_curvature
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Oct 07
c
c   purpose            - given a list of points on a surface, routine
c                        reconstructs second derivatives of an implicit
c                        quadratic surface passing through a given
c                        point on the list, along with its normal,
c                        by solving a least squares fit eigenvalue problem
c
c   arguments
c     in:
c          Xp          - a point on the surface
c          Xcoord      - coordinates of points on a surface surrounding
c                        point Xp
c          Nrpoints    - number of points
c          Rn          - normal at the point
c     out:
c          Acoeff      - second derivatives
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_curvature(Xp,Xcoord,Nrpoints,Rn,
     .                                 Acoeff)
c
#include "syscom.blk"
#include "cpmax.blk"
#include "cintr.blk"
c
      dimension Xp(3),Xcoord(3,Nrpoints),Rn(3), Acoeff(3,3)
c
c  ...points to fit
      dimension xp(3,20),weight(20)
c
c  ...matrices for a system of linear equations
      dimension aa(6,6),bb(6,3),cc(3,6),dd(3,3),uu(6),aux(6,3)
c
c  ...work space
      dimension eigen(3),work(37),iwork(18),rnw(3),aacopy(3,3),temp(3)
c
c
c**********************************************************************
c
      iprint=1
      if (Nrpoints.gt.20) then
        write(*,*) 'reconstruct_curvature: INCREASE DIMENSIONS '
        stop 1
      endif
c
      if (iprint.eq.1) then
        write(*,7006) Nrcurves,Rn
 7006   format('reconstruct_curvature: Nrcurves,Rn = ',i2,2x,3f8.3)
      endif
c
c  ...compute the secant vectors
      dmin = 1.d30
      do i=1,Nrpoints
        xp(1:3,i) = Xcoord(1:3,i) - Xp(1:3)
        call norm(xp(1:3,ic), weight(ic))
        dmin = min(dmin,weight(ic))
      enddo
      do ic=1,Nrpoints
        weight(ic) = (dmin/weight(ic))**.7d0
      enddo
      if (iprint.eq.1) then
        write(*,*) 'reconstruct_curvature: weight = '
        write(*,7035) weight(1:Nrpoints)
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
        do ic=1,Nrpoints
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
ccc      do jj=1,6
ccc        scale = max(scale,abs(aa(jj,jj)))
ccc      enddo
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
 7001   format('reconstruct_curvature: aa,bb = ')
        do jj=1,6
          write(*,7002) jj,aa(jj,1:6),bb(jj,1:3)
 7002     format(i1,3x,6e12.5,5x,3e12.5)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
      cc(1:3,1:6) = 0.d0; dd(1:3,1:3) = 0.d0
      do j=1,3
c
c  .....loop through the points
        do ic=1,Nrpoints
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
 7008     format('k,dd = ',i2,3e12.5)
        enddo
      endif
      call DSYEVD( 'V', 'U', 3, dd, 3, eigen, work, 37, iwork,
     $              18, info )
      if (iprint.eq.1) then
        write(*,7012) eigen
 7012   format('reconstruct_curvature: eigen = ',3e12.5)
      endif
c
      if (info.ne.0) then
        write(*,*) 'reconstruct_curvature: info = ',info
        stop 1
      endif
c
c-----------------------------------------------------------------------
c
c  ...choose the eigenvector which minimizes the the misfit function
      fmin = 1.d30
      do iv=1,3
        rnw(1:3) = dd(1:3,iv)
c
c  .....compute the corresponding curvature
        do i=1,6
          s = 0.d0
          do k=1,3
            s = s + aux(i,k)*rnw(k)
          enddo
          uu(i) = s
        enddo
c
c  .....store in the matrix form
        ii=0
        do i=1,3
        do j=i,3
          ii=ii+1
          Acoeff(i,j) = uu(ii)
          Acoeff(j,i) = uu(ii)
        enddo
        enddo
c
        s = 0.d0
        do i=1,3
        do j=1,3
          s = s + rnw(i)*rnw(j)*Acoeff(i,j)
        enddo
        enddo
        do i=1,3
        do j=1,3
          s1 = 0.d0; s2 = 0.d0
          do l=1,3
            s1 = s1 + rnw(l)*Acoeff(l,j)
            s2 = s2 + rnw(l)*Acoeff(l,i)
          enddo
          aacopy(i,j) = Acoeff(i,j) - rnw(i)*s1 - rnw(j)*s2
     .                + rnw(i)*rnw(j)*s
        enddo
        enddo
        if (iprint.eq.1) then
          write(*,*) 'reconstruct_curvature: Acoeff = '
          do i=1,3
            write(*,7030) i,Acoeff(i,1:3)
 7030       format('i = ',i1,2x,6e12.5)
          enddo
          call pause
        endif
c
c  .....compute principal curvatures
        call DSYEVD( 'V', 'U', 3, aacopy, 3, eigen, work, 37, iwork,
     $                18, info )
        if (iprint.eq.1) then
          write(*,7027) eigen(1:3)
 7027     format('reconstruct_curvature: PRINCIPAL CURVATURES = ',
     .            3e12.5)
        endif
c
c  .....compute the misfit function
        do isign=-1,1,2
          s = 0.d0
          do i=1,3
          do j=1,3
            s = s + Acoeff(i,j)**2
          enddo
          enddo
          s = s*scale*0.01d0
c
          do ic=1,Nrpoints
            s1=0.d0
            do k=1,3
              s1 = s1 + isign*rnw(k)*xp(k,ic)
            enddo
            do k=1,3
            do l=1,3
              s1 = s1 + .5d0*Acoeff(k,l)*xp(k,ic)*xp(l,ic)
            enddo
            enddo
            s = s + s1**2
          enddo
          if (iprint.eq.1) then
            write(*,8002) iv,rnw(1:3),s
 8002       format('reconstruct_curvature: iv,rnw,s = ',i2,
     .              3f8.3,2x,e20.12)
          endif
c
c  .......compare
          if (s.lt.fmin) then
            fmin=s; ivmin=iv; isign_min = isign
          endif
        enddo
      enddo
c
c  ...define the normal
      Rn(1:3) = dd(1:3,ivmin)*isign_min
c
c  ...compute the corresponding curvature
      do i=1,6
        s = 0.d0
        do j=1,6
          do k=1,3
            s = s - aa(i,j)*bb(j,k)*Rn(k)
          enddo
        enddo
        uu(i) = s
      enddo
c
c  ...store in the matrix form
      ii=0
      do i=1,3
      do j=i,3
        ii=ii+1
        Acoeff(i,j) = uu(ii)
        Acoeff(j,i) = uu(ii)
      enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,7003) Rn(1:3)
 7003   format('reconstruct_curvature: NEW Rn = ',3f8.3)
        write(*,7004)
 7004   format('reconstruct_curvature: Acoeff = ')
        do i=1,3
          write(*,7005) Acoeff(i,1:3)
 7005     format(3e12.5)
        enddo
        call pause
      endif
c
c
      end
