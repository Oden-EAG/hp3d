c----------------------------------------------------------------------
c
c   routine name       - reconstruct_curve
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine reconstructs a curve minimizing the
c                        L2-norm of the second derivative
c
c   arguments
c     in:
c            Xpoint    - endpoint coordinates
c            Evect     - endpoint unit vectors
c            Rnorm     - endpoint normals scaled with curvatures
c     out:
c            Evect     - endpoint tangent vectors
c            Eprim     - endpoint second derivatives
c            Nerr_flag - error flag
c                        = 0 if no error
c                        = i(=1,2) if i-th tangent vector has been
c                                  reversed
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_curve(Xpoint,Evect,Rnorm,
     .                             Eprim,Nerr_flag)
c
      use control
      use parameters       , only : MAXP
      use gauss_quadrature , only : INITIALIZED,XIGAUS1,WAGAUS1
c
      implicit none
c
      integer :: Nerr_flag
      real(8) :: Xpoint(3,2),Evect(3,2),Eprim(3,2),Rnorm(3,2)
c
c  ...matrices for a system of linear equations
      real(8) :: aa(4,4),bb(4),uu(4),u(4)
c
c  ...shape functions and their derivatives
      real(8) :: vshap(MAXP+1),dvshap(MAXP+1),ddvshap(MAXP+1)
c
c  ...dot(r) vector and its first and second derivatives wrt the
c     optimization parameters
      real(8) :: dr(3),drdl(3,4),d2drdl(3,4,4)
c
c  ...work space
      real(8) :: temp(3),eigen(4),work(9),aacopy(4,4)
c
      integer :: i,info,iter,iwork,j,l,nint
      real(8) :: dij,res,rjac2,rl,rl_old,s,sij,sj,step,wa,xi
c
      integer :: iprint
c**********************************************************************
c
      iprint=0
      Nerr_flag=0
c
c  ...set up quadrature
      if (.NOT. INITIALIZED) call init_gauss_quadrature
c
c----------------------------------------------------------------------
c
c  STEP 0: check input
c
      do i = 1,2
c  .....check the norm of the tangent versors
        call norm(Evect(1:3,i), s)
        if (abs(s-1.d0).gt.GEOM_TOL) then
          write(*,7001)  Evect(1:3,i)
 7001     format(' reconstruct_curve: Evect = ',3f8.3)
          stop 1
        endif
c  .....check orthogonality of normals and tangent versors
        call scalar_product(Evect(1:3,i),Rnorm(1:3,i), s)
        if (s.gt.GEOM_TOL) then
          write(*,7002) i,Evect(1:3,i),Rnorm(1:3,i)
 7002  format(' reconstruct_curve: i,Evect,Rnorm = ',i2,2x,2(3e12.5,2x))
          stop
        endif
      enddo
 10   continue
      if (iprint.eq.1) then
        write(*,7003) Xpoint(1:3,1),Xpoint(1:3,2)
 7003   format(' reconstruct_curve: Xpoint1,2 = ',2(3(e12.5,2x),';',2x))
        write(*,7004) Evect(1:3,1),Evect(1:3,2)
 7004   format(' reconstruct_curve: Evect1,2  = ',2(3(e12.5,2x),';',2x))
        write(*,7005) Rnorm(1:3,1),Rnorm(1:3,2)
 7005   format(' reconstruct_curve: Rnorm1,2  = ',2(3(e12.5,2x),';',2x))
      endif
c
c----------------------------------------------------------------------
c
c  STEP 1: minimize linearized curvature in order to determine
c     u(1) = s'(0), u(2) = s'(1), u(3) = s"(0), u(4) = s"(1)
c
      d2drdl = 0.d0
c
c  ...set up initial values for first and second derivatives
      temp(1:3) = Xpoint(1:3,2) - Xpoint(1:3,1)
      call norm(temp, s)
      u(1) = s; u(2) = s;
      u(3) = 0.d0; u(4) = 0.d0; uu = 0.d0
c
c  ...perform NR iterations
      do iter = 1,10
c  .....update the solution
        step = 1.d0
 30     continue
        u(1:4) = u(1:4) + step*uu(1:4)
c  .....check for reversion of first derivatives and stop NR iterations
c       if this happens
        do i = 1,2
          if (u(i).le.0.d0) then
            Nerr_flag = i
            write(*,7006) i, iter
 7006       format(' reconstruct_curve: REVERSED ',i1,
     .             ' TANGENT VECTOR IN ITERATION ',i2)
ccc            call pause
ccc         iprint=1; go to 10
          endif
        enddo
c  .....if reversion go back to the previous values
        if (Nerr_flag.ne.0) then
          u(1:4) = u(1:4) - step*uu(1:4)
          go to 20
        endif
c
        bb(1:4) = 0.d0; aa(1:4,1:4) = 0.d0; rl = 0.d0
c
c  .....integrate
        nint=9
        do l=1,nint
          xi = XIGAUS1(l,nint)
          wa  = WAGAUS1(l,nint)
          call Gshape1(5,xi, vshap,dvshap,ddvshap)
          dr(1:3) = 0.d0
          do i=1,2
            dr(1:3) = dr(1:3) + Xpoint(1:3,i)*ddvshap(i)
            dr(1:3) = dr(1:3) + u(i)*Evect(1:3,i)*ddvshap(2+i)
            dr(1:3) = dr(1:3) + (u(2+i)*Evect(1:3,i)
     .                        + u(i)**2*Rnorm(1:3,i))*ddvshap(4+i)
            drdl(1:3,i) = Evect(1:3,i)*ddvshap(2+i)
     .                  + 2.d0*u(i)*Rnorm(1:3,i)*ddvshap(4+i)
            drdl(1:3,2+i) = Evect(1:3,i)*ddvshap(4+i)
            d2drdl(1:3,i,i) = 2.d0*Rnorm(1:3,i)*ddvshap(4+i)
          enddo
c
c  .......jacobian
          call scalar_product(dr,dr, rjac2)
          rl = rl + rjac2*wa
c
c  .......loop through equations
          do j=1,4
c
c ..........accumulate for the right-hand side (residual)
            call scalar_product(dr,drdl(1:3,j), sj)
            bb(j) = bb(j) - sj*wa
c
c  .........loop through unknowns
            do i=1,4
c
c  ...........accumulate for the left-hand side
              call scalar_product(dr,d2drdl(1:3,i,j), sij)
              call scalar_product(drdl(1:3,i),drdl(1:3,j), dij)
              aa(j,i) = aa(j,i)
     .                + (dij+sij)*wa
            enddo
          enddo
c
c  .....end of loop through integration points
        enddo
c
c  .....check whether the functional has decreased
        if (iter.gt.1) then
          if (rl.gt.rl_old) then
            write(*,*) 'rl,rl_old,step = ',rl,rl_old,step
            u(1:4) = u(1:4) - step*uu(1:4)
            step = step/2.d0
            go to 30
          endif
        endif
        rl_old = rl
        res = 0.d0
        do j=1,4
          res = res + bb(j)**2
        enddo
        if (iprint.eq.1) then
          write(*,7007) iter,sqrt(res),rl,u(1:4)
 7007     format(' reconstruct_curve: iter,res,rl,u = ',i2,2x,2e12.5,3x,
     .            2(2e12.5,2x))
          call pause
        endif
        if (res.lt.GEOM_TOL) go to 20
c  .....printing
        if (iprint.eq.1) then
          write(*,7008)
 7008     format(' reconstruct_curve: bb,aa = ')
          do i = 1,4
            write(*,7009) i,bb(i),aa(i,1:4)
 7009       format(' i = ',i1,2x,e12.5,3x,4e12.5)
          enddo
          aacopy(1:4,1:4) = aa(1:4,1:4)
          call DSYEVD('N','U',4,aacopy,4,eigen,work,9,iwork,1,info)
          if (info.eq.0) then
            write(*,7010) eigen(1:4)
 7010       format(' reconstruct_curve: EIGENVALUES = ',4e12.5)
          else
            write(*,*) 'info = ',info
            stop
          endif
          call pause
        endif
c
c  .....solve the linearized system
        call gausse(aa,4,bb, uu,4)
        if (iprint.eq.1) then
          write(*,7011) uu(1:4)
 7011     format(' reconstruct_curve: uu = ',4e12.5)
        endif
      enddo
      write(*,*)'reconstruct_curve: HAVE NOT CONVERGED '
ccc      call pause; iprint=1; go to 10
ccc      call pause
 20   continue
c
c----------------------------------------------------------------------
c
c  STEP 2: compute velocities and accelerations at endpoints
c
      do i = 1,2
c  .....using eq (2.39)
        Eprim(1:3,i) = u(2+i)*Evect(1:3,i)+u(i)**2*Rnorm(1:3,i)
c  .....using eq (2.38)
        Evect(1:3,i) = u(i)*Evect(1:3,i)
      enddo
c  ...printing
      if (iprint.eq.1) then
        write(*,7003) Xpoint(1:3,1),Xpoint(1:3,2)
        write(*,7004) Evect(1:3,1),Evect(1:3,2)
        write(*,7012) Eprim(1:3,1),Eprim(1:3,2)
 7012   format(' reconstruct_curve: Eprim = ',2(3f8.3,2x))
        call pause
      endif
c
c
      end subroutine reconstruct_curve
