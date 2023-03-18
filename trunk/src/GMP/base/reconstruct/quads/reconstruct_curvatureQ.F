c----------------------------------------------------------------------
c
c   routine name       - reconstruct_curvature
c
c----------------------------------------------------------------------
c
c   latest revision    - Apr 10
c
c   purpose            - given a patch of elements surrounding a point
c                        with a normal vector, routine reconstructs
c                        second derivatives of an implicit quadratic
c                        surface passing through the point and
c                        orthogonal to the vector, by solving a least
c                        squares fit problem
c
c   arguments
c     in:
c          Np          - a point on the surface
c          List_curves - list of surface curves connected to the point
c          Nrcurves    - number of curves on the list
c          Rn          - suggested normal at the point
c     out:
c          Acoeff      - second derivatives
c          Rn          - ultimate normal at the point
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_curvature(Np,List_curves,Nrcurves,Rn,
     .                                 Acoeff)
c
      use GMP
      use control
c
c  ...removed syscom (declare variables)
      implicit none
c
      dimension List_curves(Nrcurves),Rn(3),Acoeff(3,3)
c
c  ...points to fit
      dimension xp(3,20),weight(20)
c
c  ...matrices for a system of linear equations
      dimension aa(6,6),bb(6,3),cc(3,6),dd(3,3),uu(6),aux(6,3)
c
c  ...work space
      dimension eigen(3),work(37),iwork(18),rnw(3),aGauss(3,3),temp(3)
c
c-----------------------------------------------------------------------
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,7006) Np,Nrcurves,Rn
 7006   format(' reconstruct_curvature: Np = ',i4'; Nrcurves = ',i2,';
     . Rn = ',3(e12.5))
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 0: collect point coordinates in a system with origin shifted to
c          "Np" and compute weights
c
c  ...initialize minimum weight to big number
      dmin = 1.d30
c  ...loop over curves
      do ic = 1,Nrcurves
        nc = List_curves(ic)
        do j = 1,2
          np1 = CURVES(nc)%EndPoNo(j)
          if (np1.ne.Np) then
            xp(1:3,ic) = POINTS(np1)%Rdata(1:3) - POINTS(Np)%Rdata(1:3)
          endif
        enddo
        call norm(xp(1:3,ic), weight(ic))
c  .....updated mininum weight
        dmin = min(dmin,weight(ic))
c  ...end of loop over curves
      enddo
c  ...store weights
      do ic = 1,Nrcurves
        weight(ic) = (dmin/weight(ic))**.7d0
ccc        weight(ic) = 1.d0
      enddo
c  ...printing
      if (iprint.eq.1) then
        write(*,7035) weight(1:Nrcurves)
 7035   format('   weights   = ',20(e12.5,2x))
        write(*,*)''
      endif
c
c---------------------------------------------------------------------
c
c  STEP 1: set up the matrices aa and bb for the linear system
c          corresponding to the least squares fit problem
c
      aa(1:6,1:6) = 0.d0; bb(1:6,1:3) = 0.d0
c  ...set row index for aa & bb
      jj = 0
c  ...loops over Hessian indices
      do i = 1,3
        do j = i,3
c  .......compute coefficient for a_ij
          if (i.eq.j) then
            ddd = 0.5d0
          else
            ddd = 1.d0
          endif
c  .......update row index
          jj = jj + 1
c  .......loop over connected points
          do ic = 1,Nrcurves
            do k = 1,3
              bb(jj,k) = bb(jj,k) +
     .                   ddd*xp(i,ic)*xp(j,ic)*xp(k,ic)*weight(ic)
            enddo
c  .........set column index for aa
            ii = 0
c  .........loops over Hessian indices
            do k = 1,3
              do l = k,3
                ii = ii + 1
                if (k.eq.l) then
                  ccc = 0.5d0
                else
                  ccc = 1.d0
                endif
                aa(jj,ii) = aa(jj,ii)
     .                    + ccc*xp(k,ic)*xp(l,ic)*ddd*xp(i,ic)*xp(j,ic)
     .                     *weight(ic)
c  ...........end of loops over Hessian indices
              enddo
            enddo
c  .......end of loop over connected points
          enddo
c  .....end of loops over Hessian indices
        enddo
      enddo
c
c  ...compute the maximum diagonal entry
      scale = 0.d0
      do jj = 1,6
        scale = max(scale,abs(aa(jj,jj)))
      enddo
c
c  ...add the regularizing term of two magnitudes smaller if needed
c======================================================================
c  REMARK: recall that we are storing the Hessian in the following
c    order: (11,12,13,22,23,33)
c
c           j
c         --------->
c         | X X X
c       i |   X X
c         |     X
c         |
c         V
c
c======================================================================
      if (Nrcurves.lt.6) then
        jj = 0
        do i = 1,3
          do j = i,3
            jj = jj + 1
            if (i.eq.j) then
              ddd = .01d0
            else
              ddd = .02d0
            endif
c  .........add penalty term to diagonal entry of matrix
            aa(jj,jj) = aa(jj,jj) + scale*ddd
          enddo
        enddo
      endif
c  ...printing
      if (iprint.eq.1) then
        do i = 1,6
          write(*,7001)i,aa(i,1:6)
 7001     format('   aa(',i1,',1:6) = ',6(e12.5,2x))
        enddo
        write(*,*)''
        do i = 1,6
          write(*,7002)i,bb(i,1:3)
 7002     format('   bb(',i1,',1:3) = ',3(e12.5,2x))
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 2: set up matrices cc and dd, cc being transpose(bb)
c
      cc(1:3,1:6) = 0.d0;  dd(1:3,1:3) = 0.d0
c
      do j = 1,3
c  .....loop through the points
        do ic = 1,Nrcurves
          do k = 1,3
            dd(j,k) = dd(j,k) + xp(j,ic)*xp(k,ic)*weight(ic)
          enddo
          ii = 0
          do k = 1,3
          do l = k,3
            ii = ii + 1
            if (k.eq.l) then
              ccc = 0.5d0
            else
              ccc = 1.d0
            endif
            cc(j,ii) = cc(j,ii)
     .               + ccc*xp(k,ic)*xp(l,ic)*xp(j,ic)*weight(ic)
          enddo
          enddo
        enddo
      enddo
c
c  ...check symmetry b/w bb and cc
      do i = 1,6
        do j = 1,3
          diff = abs(bb(i,j) - cc(j,i))
          if (diff.gt.1.d-12) then
            write(*,*)'i,j,diff = ',i,j,diff
            stop
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  STEP 3: set up and solve the eigenvalue problem
c
c  ...perform forward Gauss elimination
      call tri(aa,6,6)
c  ...perform backward Gauss substitution
      do k = 1,3
        call rhsub(aa,aux(1:6,k),bb(1:6,k),6,6)
      enddo
c  ...build  D - B'*A^(-1)*B , see eq (4.61), and store in dd(1:3,1:3)
      do k = 1,3
        do l = 1,3
          s = 0.d0
          do i = 1,6
            s = s + cc(k,i)*aux(i,l)
          enddo
          dd(k,l) = dd(k,l) - s
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*)''
        do k = 1,3
          write(*,7008) k,dd(k,1:3)
 7008     format('   dd(',i1,',1:3) = ',3(e12.5,2x))
        enddo
      endif
c  ...compute eigenvalues and eigenvectors
      call DSYEVD('V','U',3,dd,3,eigen,work,37,iwork,18,info)
      if (iprint.eq.1) then
        write(*,*)''
        write(*,7012) eigen
 7012   format('   eigenvalues   = ',3(e12.5,2x))
        do i = 1,3
          write(*,8000)i,dd(1:3,i)
 8000     format('   eigenvector ',i1,' = ',3(e12.5,2x))
        enddo
      endif
c
      if (info.ne.0) then
        write(*,*)'reconstruct_curvature: info = ',info
        stop
      endif
c
c----------------------------------------------------------------------
c
c  STEP 4: choose the eigenvector with the orientation closest to the
c     original normal, i.e. the best candidate
c
      dist = 1.d30
c  ...loop over eigenvectors
      do i = 1,3
        call norm(dd(1:3,i)-Rn, dnew)
        if (dnew.lt.dist) then
          sign =  1.d0;  ibest = i;  dist = dnew
        endif
        call norm(-dd(1:3,i)-Rn, dnew)
        if (dnew.lt.dist) then
          sign = -1.d0;  ibest = i;  dist = dnew
        endif
c  ...end of loop over eigenvectors
      enddo
c
c----------------------------------------------------------------------
c
c  STEP 5: check whether the selected vector indeed minimizes the
c     the misfit function, if Nrcurves >= 6.
c  REMARK: when the penalty term is added, the normal does not
c     minimizes the misfit function any more!
c
      if (Nrcurves.lt.6) then
ccc        ivmin = ibest
        goto 10
      endif

      fmin = 1.d30
c  ...loop over possible choices of normals: original + 3 eigenvectors
ccc      do iv = 0,3
      do iv = 1,3
c  .....pick normal
        select case(iv)
        case(0)
          rnw(1:3) = Rn(1:3)
        case(1,2,3)
          rnw(1:3) = dd(1:3,iv)
        end select
c  .....compute the corresponding 6 entries of the Hessian using
c       eq (4.60):  Hessian = -A^(-1) * B * n
        do i = 1,6
          s = 0.d0
          do k = 1,3
            s = s + aux(i,k)*rnw(k)
          enddo
          uu(i) = - s
        enddo
c  .....store the Hessian in matrix form
        ii = 0
        do i = 1,3
          do j = i,3
            ii = ii+1
            Acoeff(i,j) = uu(ii)
            Acoeff(j,i) = uu(ii)
          enddo
        enddo
c
c
ccc======================================================================
ccc                     G  A  U  S  S '   M  A  P
ccc
ccc    STEP 1: compute s = (A * n) * n
cc        s = 0.d0
cc        do i = 1,3
cc          do j = 1,3
cc            s = s + rnw(i)*rnw(j)*Acoeff(i,j)
cc          enddo
cc        enddo
ccc    STEP 2: compute composition of the projection on the tangent
ccc       plane with the Gauss' map. i.e. eq (2.33)
cc        do i = 1,3
cc          do j = 1,3
cc            s1 = 0.d0; s2 = 0.d0
cc            do l = 1,3
cc              s1 = s1 + rnw(l)*Acoeff(l,j)
cc              s2 = s2 + rnw(l)*Acoeff(l,i)
cc            enddo
cc            aGauss(i,j) = Acoeff(i,j) - rnw(i)*s1 - rnw(j)*s2
cc     .                  + rnw(i)*rnw(j)*s
cc          enddo
cc        enddo
cc        if (iprint.eq.1) then
cc          write(*,8050)iv
cc 8050     format('   case ',i1,':')
cc          do i = 1,3
cc            write(*,7030) i,Acoeff(i,1:3)
cc 7030       format('     Hessian(',i1,',1:3) = ',3(e12.5,2x))
cc          enddo
cc        endif
ccc    STEP 3: compute principal curvatures, i.e. eigenvectors of the
ccc       composition of the projection on the tangent plane with the
ccc       Gauss' map. Recall that one eigenvalue should ALWAYS be 0,
ccc       since there are only 2 principal curvatures!
cc        call DSYEVD('V','U',3,aGauss,3,eigen,work,37,iwork,18,info)
cc        if (iprint.eq.1) then
cc          write(*,7027) eigen(1:3)
cc 7027     format(' reconstruct_curvature: PRINCIPAL CURVATURES = ',
cc     .            3(e12.5,2x))
cc        endif
ccc======================================================================
c
c
c  .....compute scaling for misfit function
        s = 0.d0
        do i = 1,3
          do j = 1,3
            s = s + Acoeff(i,j)**2
          enddo
        enddo
        s = s*scale*0.01d0
c  .....compute misfit function
        do ic = 1,Nrcurves
c  .......compute: s1 = 0.5*(Hess*x)*x + n*x
          s1 = 0.d0
          do k = 1,3
            s1 = s1 + rnw(k)*xp(k,ic)
          enddo
          do k = 1,3
            do l = 1,3
              s1 = s1 + 0.5d0*Acoeff(k,l)*xp(k,ic)*xp(l,ic)
            enddo
          enddo
c  .......compute: s = weight*s1^2
          s = s + weight(ic)*s1**2
        enddo
c
        if (iprint.eq.1) then
          write(*,8002) iv,rnw(1:3),s
 8002     format(' reconstruct_curvature: iv,rnw,s = ',i2,
     .            3f8.3,2x,e20.12)
        endif
c  .....update if necessary
        if (s.lt.fmin) then
          fmin = s;  ivmin = iv
        endif
c  ...end of loop over possible normals
      enddo
c
c  ...check consistency
      if (ibest.ne.ivmin) then
        write(*,8001) ibest,ivmin
 8001   format(' reconstruct_curvature: INCONSISTENCY ',
     .         'ibest,ivmin = ',2i3)
        write(*,9000)Rn(1:3)
 9000   format('   orig normal   = ',3(e12.5,2x))
        do i = 1,3
          write(*,8000)i,dd(1:3,i)
        enddo
ccc        call pause
      endif
c
c----------------------------------------------------------------------
c
c  STEP 6: store definite normal and compute corresponding Hessian

c  ...reset the normal
 10   continue
      Rn(1:3) = sign*dd(1:3,ibest)
cccc  ...determine sign of normal
ccc      call norm(dd(1:3,ivmin)-Rn, dnew)
ccc      sign =  1.d0;  dist = dnew
ccc      call norm(-dd(1:3,ivmin)-Rn, dnew)
ccc      if (dnew.lt.dist) then
ccc       sign = -1.d0
ccc      endif
ccc      Rn(1:3) = sign*dd(1:3,ivmin)
c  ...compute the corresponding Hessian: Hessian = -A^(-1)*B*n
      do i = 1,6
        s = 0.d0
        do j = 1,6
          do k = 1,3
            s = s - aa(i,j)*bb(j,k)*Rn(k)
          enddo
        enddo
        uu(i) = s
      enddo
c  ...store in matrix form
      ii = 0
      do i = 1,3
        do j = i,3
          ii = ii + 1
          Acoeff(i,j) = uu(ii)
          Acoeff(j,i) = uu(ii)
        enddo
      enddo
c  ...printing
      if (iprint.eq.1) then
        write(*,*)'reconstruct_curvature: '
        write(*,7003) Rn(1:3)
 7003   format('   normal         = ',3(e12.5,2x))
        do i = 1,3
          write(*,7005)i,Acoeff(i,1:3)
 7005     format('   Hessian(',i1,',1:3) = ',3(e12.5,2x))
        enddo
        call pause
      endif
c
c
      end
