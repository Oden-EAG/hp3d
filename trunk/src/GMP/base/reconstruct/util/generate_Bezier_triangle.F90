!-----------------------------------------------------------------------------------------
!
      subroutine generate_Bezier_triangle(Nt)
!
!-----------------------------------------------------------------------------------------
!
!     latest revision: Apr 10
!
!     purpose:         routines generate BB control points for a triangle, is a
!                      way such that G^1 continuity b/w adjacent triangles is
!                      satisfied
!
!     arguments:  IN - Nt: a triangle number
!                OUT - BB control points stored in data structure
!
!-----------------------------------------------------------------------------------------
!     MODULES
      use GMP
      use control
      use bezier
!-----------------------------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------------------------
!     DUMMY ARGUMENTS
      integer,intent(in)         :: Nt
!-----------------------------------------------------------------------------------------
!     LOCAL VARIABLES
      real(8),dimension(2,3)     :: alpha
      real(8),dimension(4,3)     :: beta
      real(8),dimension(3,4,3)   :: psi
      real(8),dimension(3,2,3)   :: vel,acc,dac,dda
      real(8),dimension(3,3)     :: nor,dx_deta
      integer                    :: i,j,k,k1,k2,k3,k4,l
      integer                    :: i1,i2,i4,j1,ie,nc,nv,np,np1,np2
      real(8), dimension(3,2)    :: dual
      real(8)                    :: s1,s2,rbeta,rdbeta,rddbeta,dddbeta
      real(8),dimension(3)       :: void,der,ddpsi,x,dx,ddx,dddx,temp,rpsi,rdpsi,dder, &
                                    rddpsi,ddder,dddpsi
      real(8),dimension(3)       :: r_der,r_dder,r_ddder,temp1
      real(8),dimension(4)       :: shapef,dshapef,ddshapef
      integer                    :: iprint,iprint0,iprint1,iprint2,iprint3,iprint4,iprint5,nsub
      real(8)                    :: ddbeta,c,t,s
      integer                    :: it,lab,n,iflag
      real(8),dimension(-3:23)   :: Rdata_save,Rdata_aux
!----------------------------------------------------------------------
!     EXTERNAL PROCEDURES
      integer, external          :: imod
      integer, external          :: bijec
!----------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,2000)Nt
 2000   format(' generate_Bezier_triangle: Nt = ',i4)
      endif
!
!-----------------------------------------------------------------------
!
!  STEP 0: store edge control points and update edge curves
!
      iprint0=0
      if (iprint0.eq.1) then
        write(*,2001)
 2001   format(' -----------------------------------------------------------------')
        write(*,*)'****  STEP 0: update and store edge control points'
      endif
!
!  ...loop over edges
      do ie = 1,3
        nc = TRIANGLES(Nt)%EdgeNo(ie)
!  .....printing
        if (iprint0.eq.1) then
          write(*,1018)ie
 1018     format(' curve ',i1,' control points:')
          do i = -1,6
            write(*,1017)i,CURVES(abs(nc))%Rdata(3*i:3*i+2)
 1017       format('   b',i2,' = ',3(e12.5,2x))
          enddo
          write(*,*)'------------------------------------------------'
        endif
!  .....save edge control points
        Rdata_save = 0.d0; Rdata_aux = 0.d0
        Rdata_save(0:17) = CURVES(abs(nc))%Rdata(0:17)
!  .....degree elevate one time...
        do i = 0,6
          Rdata_aux(3*i:3*i+2) = ((6-i)*Rdata_save(3* i   :3* i   +2) +  &
                                     i *Rdata_save(3*(i-1):3*(i-1)+2))/6.d0
        enddo
!  .....save new control points
        Rdata_save = Rdata_aux
!  .....and one time more
        do i = 0,7
          Rdata_aux(3*i:3*i+2) = ((7-i)*Rdata_save(3* i   :3* i   +2) +  &
                                     i *Rdata_save(3*(i-1):3*(i-1)+2))/7.d0
        enddo
!  .....store new control points
        do i = 0,7
!  .......select edge
          select case(ie)
          case(1)
            if (nc.gt.0) then
              k = bijec(i,0)
            else
              k = bijec(7-i,0)
            endif
          case(2)
            if (nc.gt.0) then
              k = bijec(7-i,i)
            else
              k = bijec(i,7-i)
            endif
          case(3)
            if (nc.gt.0) then
              k = bijec(0,i)
            else
              k = bijec(0,7-i)
            endif
!  .......end of select edge
          endselect
          TRIANGLES(Nt)%Rdata(k:k+2) = Rdata_aux(3*i:3*i+2)
        enddo
!  ...end of loop over edges
      enddo
!  ...printing
      if (iprint0.eq.1) then
        write(*,*)'edge 1 control points'
        do i = 0,7
          k = bijec(i,0)
          write(*,1014)i,TRIANGLES(Nt)%Rdata(k:k+2)
 1014     format('   b',i1,'0 = ',3(e12.5,2x))
        enddo
        write(*,*)'------------------------------------------------'
        write(*,*)'edge 2 control points'
        do i = 0,7
          k = bijec(7-i,i)
          write(*,1015)(7-i),i,TRIANGLES(Nt)%Rdata(k:k+2)
 1015     format('   b',i1,i1,' = ',3(e12.5,2x))
        enddo
        write(*,*)'------------------------------------------------'
        write(*,*)'edge 3 control points'
        do i = 0,7
          k = bijec(0,i)
          write(*,1016)i,TRIANGLES(Nt)%Rdata(k:k+2)
 1016     format('   b0',i1,' = ',3(e12.5,2x))
        enddo
        call pause
      endif
!
!-----------------------------------------------------------------------------------
!
!  STEP 1: loop over triangle edges and store normals, and 1st, 2nd, 3rd and 4th
!     derivatives at endpoints
!
      iprint1=0
      if (iprint1.eq.1) then
        write(*,2001)
        write(*,*)'****  STEP 1: get normals, 1st,2nd,3rd,4th derivatives'
      endif
!
!  ...loop over edges
      do ie = 1,3
        nc = TRIANGLES(Nt)%EdgeNo(ie)
        nv = TRIANGLES(Nt)%VertNo(ie)
!  .....store normal
        select case (POINTS(nv)%Type)
        case ('CoorNrm')
          nor(1:3,ie) = POINTS(nv)%Rdata(4:6)
        case ('SharpPt')
          n = size(POINTS(nv)%Idata)
          iflag = 0
          do i = 1,n
            call decode(POINTS(nv)%Idata(i), it,lab)
            if (it.ne.nt) cycle
            select case(lab)
            case (1)
              nor(1:3,ie) = POINTS(nv)%Rdata(4:6)
              iflag = 1
            case (2)
              nor(1:3,ie) = POINTS(nv)%Rdata(7:9)
              iflag = 1
            case default
              write(*,*)'generate_Bezier_triangle: inconsistent label, lab = ',lab
              stop
            endselect
          enddo
          if (iflag.eq.0) then
            write(*,*)'generate_Bezier_triangle: topology inconsistency!'
            stop
          endif
        case default
          write(*,*)'generate_Bezier_triangle: inconsistent point type!'
          write(*,*)'nv,type = ',nv,POINTS(nv)%Type
          stop
        endselect
!
!  .....1st edge determines polarization of triangle
        if (ie.eq.1) then
          if (nc.gt.0) then
            np = CURVES(nc)%EndPoNo(2)
            call curve_Bezier(nc,0.d0, void,vel(1:3,1,1),acc(1:3,1,1),  &
                                            dac(1:3,1,1),dda(1:3,1,1))
            call curve_Bezier(nc,1.d0, void,vel(1:3,2,1),acc(1:3,2,1),  &
                                            dac(1:3,2,1),dda(1:3,2,1))
          else
            np = CURVES(abs(nc))%EndPoNo(1)
            call curve_Bezier(nc,1.d0, void,vel(1:3,1,1),acc(1:3,1,1),  &
                                            dac(1:3,1,1),dda(1:3,1,1))
            call curve_Bezier(nc,0.d0, void,vel(1:3,2,1),acc(1:3,2,1),  &
                                            dac(1:3,2,1),dda(1:3,2,1))
            vel(1:3,1,1) = -vel(1:3,1,1);  vel(1:3,2,1) = -vel(1:3,2,1)
            dac(1:3,1,1) = -dac(1:3,1,1);  dac(1:3,2,1) = -dac(1:3,2,1)
          endif
!  .....2nd and 3rd edges
        else
!  .......get edge endpoints
          np1 = CURVES(abs(nc))%EndPoNo(1);  np2 = CURVES(abs(nc))%EndPoNo(2)
!  .......compare to endpoints of previous edge
          if (np1.eq.np) then
            call curve_Bezier(nc,0.d0, void,vel(1:3,1,ie),acc(1:3,1,ie),  &
                                            dac(1:3,1,ie),dda(1:3,1,ie))
            call curve_Bezier(nc,1.d0, void,vel(1:3,2,ie),acc(1:3,2,ie),  &
                                            dac(1:3,2,ie),dda(1:3,2,ie))
            np = np2
          elseif (np2.eq.np) then
            call curve_Bezier(nc,1.d0, void,vel(1:3,1,ie),acc(1:3,1,ie),  &
                                            dac(1:3,1,ie),dda(1:3,1,ie))
            call curve_Bezier(nc,0.d0, void,vel(1:3,2,ie),acc(1:3,2,ie),  &
                                            dac(1:3,2,ie),dda(1:3,2,ie))
            vel(1:3,1,ie) = -vel(1:3,1,ie);  vel(1:3,2,ie) = -vel(1:3,2,ie)
            dac(1:3,1,ie) = -dac(1:3,1,ie);  dac(1:3,2,ie) = -dac(1:3,2,ie)
            np = np1
          else
            write(*,*)'trian_G1RecTri: geometry inconsistency!'
            stop
          endif
        endif
!  ...end of loop over edges
      enddo
!
!  CHECK
      do i = 1,3
        call cross_product(vel(1:3,2,imod(i+2,3)),vel(1:3,1,i), temp)
        call normalize(temp)
!!!        write(*,*)'cross prod = ',temp
!!!        write(*,*)'nor        = ',nor(1:3,i)
        call cross_product(temp,nor(1:3,i), temp1)
!!!        write(*,*)'temp       = ',temp1
        call norm(temp1, s)
        if (s.gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency upon input!'
          write(*,112)i,s
 112      format('   i = ',i1,';   s = ',e12.5)
          write(*,113)nor(1:3,i)
 113      format('   nor        = ',3(e12.5,2x))
          call cross_product(vel(1:3,2,imod(i+2,3)),vel(1:3,1,i), temp)
          call normalize(temp)
          write(*,114)temp
 114      format('   cross prod = ',3(e12.5,2x))
          call pause
        endif
      enddo
!
!  ...printing
      if (iprint1.eq.1) then
        do j = 1,3
          nc = TRIANGLES(Nt)%EdgeNo(j)
          write(*,1001)j,nc
 1001     format(' ie = ',i1,'; nc = ',i4)
          do i = 1,2
            write(*,1000)i,vel(1:3,i,j)
 1000       format('   endpoint = ',i1,' --> vel = ',3(e12.5,2x))
            write(*,50)acc(1:3,i,j)
   50       format('                    acc = ',3(e12.5,2x))
            write(*,51)dac(1:3,i,j)
   51       format('                    dac = ',3(e12.5,2x))
            write(*,52)dda(1:3,i,j)
   52       format('                    dda = ',3(e12.5,2x))
          enddo
        enddo
        do j = 1,3
          write(*,1006)j,nor(1:3,j)
 1006     format(' iv = ',i1,'         --> nor = ',3(e12.5,2x))
        enddo
        call pause
      endif
!
!
!-----------------------------------------------------------------------------------
!
!  STEP 2: determine mixed derivative at endpoints
!
      iprint2=0
      if (iprint2.eq.1) then
        write(*,2001)
        write(*,*)'****  STEP 2: determine mixed derivative at vertices'
      endif
!
!  ...loop over edges
      do ie = 1,3
!  .....printing
        if (iprint2.eq.1) then
          nc = TRIANGLES(Nt)%EdgeNo(ie)
          write(*,1001)ie,nc
        endif
!  STEP 2a: determine psi at endpoints, psi = n x tan
        call cross_product(nor(1:3,ie),vel(1:3,1,ie), psi(1:3,1,ie))
        call cross_product(nor(1:3,imod(ie+1,3)),vel(1:3,2,ie), psi(1:3,2,ie))
        call normalize(psi(1:3,1,ie))
        call normalize(psi(1:3,2,ie))
        if (iprint2.eq.1) then
          write(*,1002)psi(1:3,1,ie)
 1002     format('   psi(0)   = ',3(e12.5,2x))
          write(*,999)psi(1:3,2,ie)
 999      format('   psi(1)   = ',3(e12.5,2x))
        endif
!  STEP 2b: determine alpha at endpoints using eq (2.48)
!  .....1st endpoint
        call scalar_product(vel(1:3,1,ie),vel(1:3,2,imod(ie+2,3)), s1)
        call scalar_product(vel(1:3,1,ie),vel(1:3,1,ie), s2)
        alpha(1,ie) = -s1/s2
!  .....2nd endpoint
        call scalar_product(vel(1:3,2,ie),vel(1:3,1,imod(ie+1,3)), s1)
        call scalar_product(vel(1:3,2,ie),vel(1:3,2,ie), s2)
        alpha(2,ie) = s1/s2
        if (iprint2.eq.1) then
          write(*,1003)alpha(1,ie)
 1003     format('   alpha(0) = ',e12.5)
          write(*,998) alpha(2,ie)
 998      format('   alpha(1) = ',e12.5)
        endif
!  STEP 2c: determine beta at endpoints using eq (2.49)
!  .....1st endpoint
        call scalar_product(vel(1:3,2,imod(ie+2,3)),psi(1:3,1,ie), beta(1,ie))
        beta(1,ie) = -beta(1,ie)
!  .....2nd endpoint
        call scalar_product(vel(1:3,1,imod(ie+1,3)),psi(1:3,2,ie), beta(2,ie))
        if (iprint2.eq.1) then
          write(*,1004)beta(1,ie)
 1004     format('   beta(0)  = ',e12.5)
          write(*,997)beta(2,ie)
 997      format('   beta(1)  = ',e12.5)
        endif
!  STEP 2d: compute dual basis (v^1, v^2, n)
!  .....1st dual basis vector for edge ie  : v^1
        call cross_product(nor(1:3,ie),vel(1:3,2,imod(ie+2,3)), dual(1:3,1))
        call scalar_product(dual(1:3,1),vel(1:3,1,ie), s1)
        dual(1:3,1) = dual(1:3,1)/s1
!  .....2nd dual basis vector for edge ie+1: v^2
        call cross_product(nor(1:3,imod(ie+1,3)),vel(1:3,1,imod(ie+1,3)), dual(1:3,2))
        call scalar_product(dual(1:3,2),vel(1:3,2,ie), s1)
        dual(1:3,2) = -dual(1:3,2)/s1
!
!=====================================================================
!  .....check
        call scalar_product(dual(1:3,1),nor(1:3,ie), s1)
        if (abs(s1).gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency!'
          write(*,9001)ie,s1
 9001     format('   ie = ',i1,' :  v^1 * n = ',e12.5)
          call pause
        endif
        call scalar_product(dual(1:3,1),vel(1:3,2,imod(ie+2,3)), s1)
        if (abs(s1).gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency!'
          write(*,9002)ie,s1
 9002     format('   ie = ',i1,' :  v^1 * v_2 = ',e12.5)
          call pause
        endif
        call scalar_product(dual(1:3,1),vel(1:3,1,ie), s1)
        if (abs(s1-1.d0).gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency!'
          write(*,9003)ie,s1
 9003     format('   ie = ',i1,' :  v^1 * v_1 = ',e12.5)
          call pause
        endif
!
        call scalar_product(dual(1:3,2),nor(1:3,imod(ie+1,3)), s1)
        if (abs(s1).gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency!'
          write(*,9001)ie,s1
 9004     format('   ie = ',i1,' :  v^2 * n = ',e12.5)
          call pause
        endif
        call scalar_product(dual(1:3,2),vel(1:3,1,imod(ie+1,3)), s1)
        if (abs(s1).gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency!'
          write(*,9005)ie,s1
 9005     format('   ie = ',i1,' :  v^2 * v_1 = ',e12.5)
          call pause
        endif
        call scalar_product(dual(1:3,2),-vel(1:3,2,ie), s1)
        if (abs(s1-1.d0).gt.GEOM_TOL) then
          write(*,*)'generate_Bezier_triangle: inconsistency!'
          write(*,9006)ie,s1
 9006     format('   ie = ',i1,' :  v^1 * v_1 = ',e12.5)
          call pause
        endif
!=====================================================================
!
!  .....printing
        if (iprint2.eq.1) then
          write(*,1005)dual(1:3,1)
 1005     format('   v^1      = ',3(e12.5,2x))
          write(*,996) dual(1:3,2)
 996      format('   v^2      = ',3(e12.5,2x))
        endif
!  STEP 2e: determine components in the dual basis using eq (2.51)
!  .....component along v^1
        call scalar_product(vel(1:3,1,ie),vel(1:3,1,ie), s1)
        call scalar_product(vel(1:3,1,ie),acc(1:3,1,ie), s2)
        dual(1:3,1) = ((alpha(2,ie)-alpha(1,ie) + 1.d0)*s1 + alpha(1,ie)*s2)*dual(1:3,1)
!  .....component along v^2
        call scalar_product(vel(1:3,2,ie),vel(1:3,2,ie), s1)
        call scalar_product(vel(1:3,2,ie),acc(1:3,2,ie), s2)
        dual(1:3,2) = ((alpha(2,ie)-alpha(1,ie) + 1.d0)*s1 + alpha(2,ie)*s2)*dual(1:3,2)
        if (iprint2.eq.1) then
          write(*,1009)dual(1:3,1)
 1009     format('   v^1 comp = ',3(e12.5,2x))
          write(*,995)dual(1:3,2)
 995      format('   v^2 comp = ',3(e12.5,2x))
        endif
!  STEP 2f: update mixed derivative
        select case(ie)
        case(1)
          i = 1; j = 1; i1 = 5; j1 = 1
        case(2)
          i = 5; j = 1; i1 = 1; j1 = 5
        case(3)
          i = 1; j = 5; i1 = 1; j1 = 1
        endselect
        k = bijec(i,j); l = bijec(i1,j1)
        TRIANGLES(Nt)%Rdata(k:k+2) = TRIANGLES(Nt)%Rdata(k:k+2) + dual(1:3,1)
        TRIANGLES(Nt)%Rdata(l:l+2) = TRIANGLES(Nt)%Rdata(l:l+2) + dual(1:3,2)
!  ...end of loop over edges
      enddo
!  ...printing
      if (iprint2.eq.1) then
        do ie = 1,3
          select case(ie)
          case(1)
            i = 1; j = 1
          case(2)
            i = 5; j = 1
          case(3)
            i = 1; j = 5
          endselect
          k = bijec(i,j)
          write(*,1010)ie,TRIANGLES(Nt)%Rdata(k:k+2)
 1010     format('   ie = ',i1,'; mixed derivative = ',3(e12.5,2x))
        enddo
      endif
!
!======================================================================
!  ...check consistency of radial derivative at endpoints
      do ie = 1,3
        temp = alpha(1,ie)*vel(1:3,1,ie) + beta(1,ie)*psi(1:3,1,ie)
        temp = temp + vel(1:3,2,imod(ie+2,3))
        call norm(temp, s)
        if (s.gt.GEOM_TOL) then
          write(*,*)'  INCONSISTENCY (1), ie = ',ie,'; s = ',s
          call pause
        endif
        temp = alpha(2,ie)*vel(1:3,2,ie) + beta(2,ie)*psi(1:3,2,ie)
        temp = temp - vel(1:3,1,imod(ie+1,3))
        call norm(temp, s)
        if (s.gt.GEOM_TOL) then
          write(*,*)'  INCONSISTENCY (2), ie = ',ie,'; s = ',s
          call pause
        endif
      enddo
!
!  ...loop over edges
      do ie = 1,3
!  .....get curve number and endpoints
        nc = TRIANGLES(Nt)%EdgeNo(ie)
        np1 = CURVES(abs(nc))%EndPoNo(1);  np2 = CURVES(abs(nc))%EndPoNo(2)
!  .....loop over endpoints
        do i = 0,1
          t = i*1.d0
!  .......1st edge determines triangle polarization
          if (ie.eq.1) then
            if (nc.gt.0) then
              call curve_Bezier(nc,t, void,dx)
!  ...........update endpoint if last visit
              if (i.eq.1) np = np2
            else
              call curve_Bezier(abs(nc),1.d0-t, void,dx)
              dx = -dx
!  ...........update endpoint if last visit
              if (i.eq.1) np = np1
            endif
!  .......2nd and 3rd curves
          else
            if (np1.eq.np) then
              call curve_Bezier(abs(nc),t, void,dx)
!  ...........update endpoint if last visit
              if (i.eq.1) np = np2
            else
              call curve_Bezier(abs(nc),1.d0-t, void,dx)
              dx = -dx
!  ...........update endpoint if last visit
              if (i.eq.1) np = np1
            endif
          endif
!  .......compute desired radial derivative
          der = alpha(i+1,ie)*dx + beta(i+1,ie)*psi(1:3,i+1,ie)
!  .......compare
          select case(i)
          case(0)
            r_der = - vel(1:3,2,imod(ie+2,3))
          case(1)
            r_der =   vel(1:3,1,imod(ie+1,3))
          endselect
          temp = der - r_der
          call norm(temp, s1)
          if (s1.gt.GEOM_TOL) then
            write(*,6000)ie,t,der
 6000       format('  INCONSISTENCY (3): ie = ',i1,'; t = ',e12.5,' --> p_der   = ',3(e12.5,2x))
            write(*,6001)t,r_der
 6001       format('                             t = ',e12.5,' --> r_der   = ',3(e12.5,2x))
            call pause
          endif
!  .....end of loop over subdivisions
        enddo
!  ...end of loop over edges
      enddo
!======================================================================
!
!
!----------------------------------------------------------------------------------
!
!  STEP 3: compute beta' and psi' using eq (2.50)
!
      iprint3=0
      if (iprint3.eq.1) then
        write(*,2001)
        write(*,*)'****  STEP 3: determine dbeta and dpsi'
      endif
!
!  ...loop over edges
      do ie = 1,3
!  .....select edge
        select case(ie)
        case(1)
          i = 1; j = 1
        case(2)
          i = 5; j = 1
        case(3)
          i = 1; j = 5
        endselect
        k = bijec(i,j)
!  .....beta'(0) for edge ie
        call scalar_product(TRIANGLES(Nt)%Rdata(k:k+2),psi(1:3,1,ie), s1)
        call scalar_product(acc(1:3,1,ie),psi(1:3,1,ie), s2)
        beta(3,ie) = s1 - alpha(1,ie)*s2
!  .....psi'(0) for edge ie
        psi(1:3,3,ie) = (TRIANGLES(Nt)%Rdata(k:k+2) -                      &
                         (alpha(2,ie)-alpha(1,ie) + 1.d0)*vel(1:3,1,ie) -  &
                         alpha(1,ie)*acc(1:3,1,ie) -                       &
                         beta(3,ie)*psi(1:3,1,ie))/beta(1,ie)
!  .....beta'(1) for edge ie+2
        call scalar_product(TRIANGLES(Nt)%Rdata(k:k+2),psi(1:3,2,imod(ie+2,3)), s1)
        call scalar_product(acc(1:3,2,imod(ie+2,3)),psi(1:3,2,imod(ie+2,3)), s2)
        beta(4,imod(ie+2,3)) = - s1 - alpha(2,imod(ie+2,3))*s2
! ..... psi'(1) for edge ie+2
        psi(1:3,4,imod(ie+2,3)) = (-TRIANGLES(Nt)%Rdata(k:k+2) -                                                     &
                                   (1.d0 + alpha(2,imod(ie+2,3)) - alpha(1,imod(ie+2,3)))*vel(1:3,2,imod(ie+2,3)) -  &
                                   alpha(2,imod(ie+2,3))*acc(1:3,2,imod(ie+2,3)) -                                   &
                                   beta(4,imod(ie+2,3))*psi(1:3,2,imod(ie+2,3)))/beta(2,imod(ie+2,3))
!  ...end of loop over edges
      enddo
!
!  ...check that tangential component of psi' vanishes at endpoints
      do ie = 1,3
        do j = 0,1
          call scalar_product(psi(1:3,3+j,ie),nor(1:3,imod(ie+j,3)), s)
          temp = psi(1:3,3+j,ie) - s*nor(1:3,imod(ie+j,3))
          call norm(temp, s)
          if (s.gt.GEOM_TOL) then
            write(*,8999)j
 8999       format(' generate_Bezier_triangle: tangential component of dpsi at ',i1,' does not vanish')
            write(*,9000)ie,s
 9000       format('   ie = ',i1,'; s = ',e12.5)
            call pause
          endif
        enddo
      enddo
!
!
!!!!!
!!!!!  ...check
!!!!      if (iprint3.eq.1) then
!!!!!  .....loop over edges
!!!!        do ie = 1,3
!!!!!  .......printing
!!!!          nc = TRIANGLES(Nt)%EdgeNo(ie)
!!!!          write(*,1001)ie,nc
!!!!          write(*,1011)beta(3,ie)
!!!! 1011     format('   dbeta(0) = ',e12.5)
!!!!          write(*,994) beta(4,ie)
!!!! 994      format('   dbeta(1) = ',e12.5)
!!!!          write(*,1008)psi(1:3,3,ie)
!!!! 1008     format('   dpsi(0)  = ',3(e12.5,2x))
!!!!          write(*,993) psi(1:3,4,ie)
!!!! 993      format('   dpsi(1)  = ',3(e12.5,2x))
!!!!!  .......check radial derivative
!!!!          nsub = 10
!!!!!  .......loop over subdivisions
!!!!          do i = 0,nsub
!!!!            t = i*1.d0/nsub
!!!!            if (nc.gt.0) then
!!!!              call curve_Bezier(nc,t, x,dx,ddx,dddx)
!!!!            else
!!!!              call curve_Bezier(abs(nc),1.d0-t, x,dx,ddx,dddx)
!!!!              dx = -dx; dddx = -dddx
!!!!            endif
!!!!            rpsi = 0.d0; rbeta = 0.d0
!!!!            rdpsi = 0.d0; rdbeta = 0.d0
!!!!            rddpsi = 0.d0; rddbeta = 0.d0
!!!!            call Hshape1(3,t, shapef,dshapef,ddshapef)
!!!!!  .........loop over dof's
!!!!            do j = 1,4
!!!!              rpsi    = rpsi   + psi(1:3,j,ie)*shapef(j)
!!!!              rdpsi   = rdpsi  + psi(1:3,j,ie)*dshapef(j)
!!!!!              rddpsi  = rddpsi  + psi(1:3,j,ie)*ddshapef(j)
!!!!!
!!!!              rbeta   = rbeta  + beta(j,ie)*shapef(j)
!!!!              rdbeta  = rdbeta + beta(j,ie)*dshapef(j)
!!!!!              rddbeta = rddbeta + beta(j,ie)*ddshapef(j)
!!!!!  .........end of loop over dof's
!!!!            enddo
!!!!            der = (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dx + rbeta*rpsi
!!!!!!!            write(*,5000)t,der
!!!! 5000       format('   t = ',e12.5,' --> p_der = ',3(e12.5,2x))
!!!!!  .......end of loop over subdivision
!!!!          enddo
!!!!!  .....end of loop over edges
!!!!        enddo
!!!!      endif
!
!
!------------------------------------------------------------------------------------------------------
!
!  STEP 4: determine corner control points
!
      iprint4=0
      if (iprint4.eq.1) then
        write(*,2001)
        write(*,*)'****  STEP 4: determine corner control points'
      endif
!  ...loop over edges
      do ie = 1,3
        select case(ie)
        case(1)
          i = 1; j = 1
        case(2)
          i = 5; j = 1
        case(3)
          i = 1; j = 5
        endselect
        k = bijec(i,j)
        if (iprint4.eq.1) then
          write(*,3776)ie, TRIANGLES(Nt)%Rdata(k:k+2)
 3776     format('   iv = ',i1,' --> mixed derivative = ',3(e12.5,2x))
        endif
!  .....select edge
        select case(ie)
        case(1)
          k1 = bijec(0,0);  k2 = bijec(1,0);  k3 = bijec(0,1);  k4 = bijec(1,1)
        case(2)
          k1 = bijec(7,0);  k2 = bijec(6,1);  k3 = bijec(6,0);  k4 = bijec(5,1)
        case(3)
          k1 = bijec(0,7);  k2 = bijec(0,6);  k3 = bijec(1,6);  k4 = bijec(1,5)
        endselect
!  .....accumulate
        TRIANGLES(Nt)%Rdata(k4:k4+2) = TRIANGLES(Nt)%Rdata(k4:k4+2)/42.d0 +  &
                                       TRIANGLES(Nt)%Rdata(k2:k2+2) +        &
                                       TRIANGLES(Nt)%Rdata(k3:k3+2) -        &
                                       TRIANGLES(Nt)%Rdata(k1:k1+2)
!  .....printing
        if (iprint4.eq.1) then
          write(*,1013)ie,TRIANGLES(Nt)%Rdata(k4:k4+2)
 1013     format('   ie = ',i1,';  b = ',3(e12.5,2x))
        endif
!
!  ...end of loop over edges
      enddo
!
!--------------------------------------------------------------------------------------------
!
!  STEP 5: determine b_21, b_41; b_42, b_24; b_14, b_12
!
      iprint5=0
      if (iprint5.eq.1) then
        write(*,2001)
        write(*,*)'****  STEP 5: determine inner control points'
      endif
!
!  ...loop over edges
      do ie = 1,3
!  STEP 5a: compute 2nd derivative of radial derivative at 1st endpoint
!  .....compute beta"(0) and psi"(0)
        call Hshape1(3,0.d0, shapef,dshapef,ddshapef)
        ddbeta = 0.d0; ddpsi = 0.d0
        do i = 1,4
          ddbeta = ddbeta + beta(i,ie)*ddshapef(i)
          ddpsi  = ddpsi  + psi(1:3,i,ie)*ddshapef(i)
        enddo
        if (iprint5.eq.1) then
          write(*,1019)ie,ddbeta,ddpsi
 1019     format('   ie = ',i1,' --> ddbeta(0) = ',e12.5,'; ddpsi(0) = ',3(e12.5,2x))
        endif
!  .....compute derivative
        der = 2.d0*(alpha(2,ie)-alpha(1,ie))*acc(1:3,1,ie) + alpha(1,ie)*dac(1:3,1,ie) +  &
              ddbeta*psi(1:3,1,ie) + 2.d0*beta(3,ie)*psi(1:3,3,ie) + beta(1,ie)*ddpsi
!  .....printing
        if (iprint5.eq.1) then
          write(*,1021)ie,der
 1021     format('   ie = ',i1,' --> p_ddder(0) = ',3(e12.5,2x))
        endif
!  .....store derivative according to edge
        select case(ie)
        case(1)
          TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2) = der
        case(2)
          TRIANGLES(Nt)%Rdata(bijec(4,2):bijec(4,2)+2) = der
        case(3)
          TRIANGLES(Nt)%Rdata(bijec(1,4):bijec(1,4)+2) = der
        endselect
!  STEP 5b: determine 1st control point
!  .....loop over control points involved
        do j = 0,1
          do i = 0,2
            k = 7-i-j
!  .........cycle if on control point
            if ((i.eq.2).and.(j.eq.1)) cycle
!  .........compute coefficient
!!!            call compute_coefficients(0.d0,i,j, c)
            call rad_der_coefficient(0.d0,i,j,2, c)
            select case(ie)
            case(1)
              l = bijec(i,j)
              TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            case(2)
              l = bijec(k,i)
              TRIANGLES(Nt)%Rdata(bijec(4,2):bijec(4,2)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(4,2):bijec(4,2)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            case(3)
              l = bijec(j,k)
              TRIANGLES(Nt)%Rdata(bijec(1,4):bijec(1,4)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(1,4):bijec(1,4)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            endselect
          enddo
!  .....end of loop over control points involved
        enddo
!  .....divide by coefficient
!!!        call compute_coefficients(0.d0,2,1, c)
        call rad_der_coefficient(0.d0,2,1,2, c)
        select case(ie)
        case(1)
          TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2)/c
          if (iprint5.eq.1) then
            write(*,1030)ie,TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2)
 1030       format('   ie = ',i1,' --> b21 = ',3(e12.5,2x))
          endif
        case(2)
          TRIANGLES(Nt)%Rdata(bijec(4,2):bijec(4,2)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(4,2):bijec(4,2)+2)/c
          if (iprint5.eq.1) then
            write(*,1031)ie,TRIANGLES(Nt)%Rdata(bijec(3,2):bijec(3,2)+2)
 1031       format('   ie = ',i1,' --> b32 = ',3(e12.5,2x))
          endif
        case(3)
          TRIANGLES(Nt)%Rdata(bijec(1,4):bijec(1,4)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(1,4):bijec(1,4)+2)/c
          if (iprint5.eq.1) then
            write(*,1032)ie,TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2)
 1032       format('   ie = ',i1,' --> b13 = ',3(e12.5,2x))
          endif
       endselect
!  STEP 5c: compute 2nd deterivative of radial derivative at 2nd endpoint
!  .....compute beta"(1) and psi"(1)
        call Hshape1(3,1.d0, shapef,dshapef,ddshapef)
        ddbeta = 0.d0; ddpsi = 0.d0
        do i = 1,4
          ddbeta = ddbeta + beta(i,ie)*ddshapef(i)
          ddpsi  = ddpsi  + psi(1:3,i,ie)*ddshapef(i)
        enddo
!  .....printing
        if (iprint5.eq.1) then
          write(*,1020)ie,ddbeta,ddpsi
 1020     format('   ie = ',i1,' --> ddbeta(1) = ',e12.5,'; ddpsi(1) = ',3(e12.5,2x))
        endif
!  .....compute derivative
        der = 2.d0*(alpha(2,ie)-alpha(1,ie))*acc(1:3,2,ie) + alpha(2,ie)*dac(1:3,2,ie) +  &
              ddbeta*psi(1:3,2,ie) + 2.d0*beta(4,ie)*psi(1:3,4,ie) + beta(2,ie)*ddpsi
!  .....printing
        if (iprint5.eq.1) then
          write(*,4021)ie,der
 4021     format('   ie = ',i1,' --> p_ddder(1) = ',3(e12.5,2x))
        endif
!  .....store derivative
        select case(ie)
        case(1)
          TRIANGLES(Nt)%Rdata(bijec(4,1):bijec(4,1)+2) = der
        case(2)
          TRIANGLES(Nt)%Rdata(bijec(2,4):bijec(2,4)+2) = der
        case(3)
          TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2) = der
        endselect
!
!  STEP 5d: compute 2nd control point
!  .....loop over involved control points
        do k = 0,2
          do j = 0,1
            i = 7-j-k
!  .........cycle if on control point
            if ((k.eq.2).and.(j.eq.1)) cycle
!  .........compute coefficient
!!!            call compute_coefficients(1.d0,i,j, c)
            call rad_der_coefficient(1.d0,i,j,2, c)
            select case(ie)
            case(1)
              l = bijec(i,j)
              TRIANGLES(Nt)%Rdata(bijec(4,1):bijec(4,1)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(4,1):bijec(4,1)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            case(2)
              l = bijec(k,i)
              TRIANGLES(Nt)%Rdata(bijec(2,4):bijec(2,4)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(2,4):bijec(2,4)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            case(3)
              l = bijec(j,k)
              TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            endselect
          enddo
!  .....end of loop over involved control points
        enddo
!  .....divide by coefficient
!!!        call compute_coefficients(1.d0,4,1, c)
        call rad_der_coefficient(1.d0,4,1,2, c)
        select case(ie)
        case(1)
          TRIANGLES(Nt)%Rdata(bijec(4,1):bijec(4,1)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(4,1):bijec(4,1)+2)/c
          if (iprint5.eq.1) then
            write(*,1033)ie,TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2)
 1033       format('   ie = ',i1,' --> b31 = ',3(e12.5,2x))
          endif
        case(2)
          TRIANGLES(Nt)%Rdata(bijec(2,4):bijec(2,4)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(2,4):bijec(2,4)+2)/c
          if (iprint5.eq.1) then
            write(*,1034)ie,TRIANGLES(Nt)%Rdata(bijec(2,3):bijec(2,3)+2)
 1034       format('   ie = ',i1,' --> b23 = ',3(e12.5,2x))
          endif
        case(3)
          TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2)/c
          if (iprint5.eq.1) then
            write(*,1035)ie,TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2)
 1035       format('   ie = ',i1,' --> b12 = ',3(e12.5,2x))
          endif
       endselect
!
!  ...end of loop over edges
      enddo
!
!------------------------------------------------------------------------------------
!
!  STEP 6: determine b_31; b_33; b_13
!
!  ...loop over edges
      do ie = 1,3
!  STEP 6a: compute 3rd derivative of radial derivative at 1st endpoint
!  .....compute beta"(0) and psi"(0)
        call Hshape1(3,0.d0, shapef,dshapef,ddshapef)
        ddbeta = 0.d0; ddpsi = 0.d0
        do i = 1,4
          ddbeta = ddbeta + beta(i,ie)*ddshapef(i)
          ddpsi  = ddpsi  + psi(1:3,i,ie)*ddshapef(i)
        enddo
        dddbeta = 12.d0*(beta(1,ie) - beta(2,ie)) + 6.d0*(beta(3,ie) + beta(4,ie))
        dddpsi  = 12.d0*(psi(1:3,1,ie) - psi(1:3,2,ie)) + 6.d0*(psi(1:3,3,ie) + psi(1:3,4,ie))
!  .....compute derivative
        der = 3.d0*(alpha(2,ie)-alpha(1,ie))*dac(1:3,1,ie) + alpha(1,ie)*dda(1:3,1,ie) +        &
                    dddbeta*psi(1:3,1,ie) + 3.d0*ddbeta*psi(1:3,3,ie) + 3.d0*beta(3,ie)*ddpsi + &
                    beta(1,ie)*dddpsi
!  .....store derivative according to edge
        select case(ie)
        case(1)
          TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) = der
        case(2)
          TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2) = der
        case(3)
          TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2) = der
        endselect
!  STEP 6b: determine control point
!  .....loop over control points involved
        do j = 0,1
          do i = 0,3
            k = 7-i-j
!  .........cycle if on control point
            if ((i.eq.3).and.(j.eq.1)) cycle
!  .........compute coefficient
            call rad_der_coefficient(0.d0,i,j,3, c)
            select case(ie)
            case(1)
              l = bijec(i,j)
              TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            case(2)
              l = bijec(k,i)
              TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            case(3)
              l = bijec(j,k)
              TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2) =      &
                  TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2) -  &
                  c*TRIANGLES(Nt)%Rdata(l:l+2)
            endselect
          enddo
!  .....end of loop over control points involved
        enddo
!  .....divide by coefficient
        call rad_der_coefficient(0.d0,3,1,3, c)
        select case(ie)
        case(1)
          TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2)/c
        case(2)
          TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2)/c
        case(3)
          TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2) =          &
            TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2)/c
        endselect
!  ...end of loop over edges
      enddo
!
!-----------------------------------------------------------------------------------------
!
!  STEP 7: determine inner control points b_22; b_32; b_23
!
      TRIANGLES(Nt)%Rdata(bijec(2,2):bijec(2,2)+2) =                    &
               (TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) +          &
                TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2)   )/2.d0
!
      TRIANGLES(Nt)%Rdata(bijec(3,2):bijec(3,2)+2) =                    &
               (TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) +          &
                TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2)   )/2.d0
!
      TRIANGLES(Nt)%Rdata(bijec(2,3):bijec(2,3)+2) =                    &
               (TRIANGLES(Nt)%Rdata(bijec(3,3):bijec(3,3)+2) +          &
                TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2)   )/2.d0





!!!!
!!!!  STEP 6: determine inner control point
!!!!
!!!      TRIANGLES(Nt)%Rdata(bijec(2,2):bijec(2,2)+2) =                               &
!!!                              (TRIANGLES(Nt)%Rdata(bijec(2,1):bijec(2,1)+2) +      &
!!!                               TRIANGLES(Nt)%Rdata(bijec(3,1):bijec(3,1)+2) +      &
!!!                               TRIANGLES(Nt)%Rdata(bijec(3,2):bijec(3,2)+2) +      &
!!!                               TRIANGLES(Nt)%Rdata(bijec(2,3):bijec(2,3)+2) +      &
!!!                               TRIANGLES(Nt)%Rdata(bijec(1,3):bijec(1,3)+2) +      &
!!!                               TRIANGLES(Nt)%Rdata(bijec(1,2):bijec(1,2)+2))/6.d0
!!!
!
!  ...print control points
      iprint=0
      if (iprint.eq.1) then
        do j = 0,7
          do i = 0,(7-j)
            l = bijec(i,j)
            write(*,4000)i,j,TRIANGLES(Nt)%Rdata(l:l+2)
 4000       format('   b',i1,i1,' = ',3(e12.5,2x))
          enddo
        enddo
        call pause
      endif

!!!      return








!!!      write(*,*)'*********************** C H E C K I N G  **********************************'
!!!      do ie = 1,3
!!!!       1 ST      E  N  D  P  O  I  N  T
!!!        t = 0.d0
!!!        call compute_radial_derivative_new(Nt,ie,t, r_der,r_dder,r_ddder)
!!!!  .....check radial derivative 1st method
!!!        temp = r_der + vel(1:3,2,imod(ie+2,3))
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,8998)ie,t,s1
!!! 8998     format('   inconsistency (1st meth)! ie = ',i1,', t = ',e6.2,', s = ',e12.5)
!!!          call pause
!!!        endif
!!!        select case(ie)
!!!        case(1)
!!!          k1 = bijec(1,1); k2 = bijec(1,0); k3 = bijec(0,1); k4 = bijec(0,0)
!!!        case(2)
!!!          k1 = bijec(4,1); k2 = bijec(5,1); k3 = bijec(5,0); k4 = bijec(6,0)
!!!        case(3)
!!!          k1 = bijec(1,4); k2 = bijec(0,5); k3 = bijec(1,5); k4 = bijec(0,6)
!!!        endselect
!!!!  .....check radial derivative 2nd method
!!!        temp = 6.d0*(TRIANGLES(Nt)%Rdata(k3:k3+2) - TRIANGLES(Nt)%Rdata(k4:k4+2))
!!!        temp = temp - r_der
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,8995)ie,t,s1
!!! 8995     format('   inconsistency (2nd meth)! ie = ',i1,', t = ',e6.2,', s = ',e12.5)
!!!          call pause
!!!        endif
!!!!  .....check derivative of radial derivative
!!!        temp = 30.d0*(TRIANGLES(Nt)%Rdata(k1:k1+2) - TRIANGLES(Nt)%Rdata(k2:k2+2) -   &
!!!                      TRIANGLES(Nt)%Rdata(k3:k3+2) + TRIANGLES(Nt)%Rdata(k4:k4+2)) -  &
!!!                6.d0*(TRIANGLES(Nt)%Rdata(k2:k2+2) - TRIANGLES(Nt)%Rdata(k4:k4+2))
!!!        temp = temp - r_dder
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,8996)ie,t,s1
!!! 8996     format('   inconsistency with 1st derivative! ie = ',i1,', t = ',e6.2,', s = ',e12.5)
!!!          call pause
!!!        endif
!!!!  .....check agreement with desired radial derivative
!!!        nc = TRIANGLES(Nt)%EdgeNo(ie)
!!!        if (nc.gt.0) then
!!!          call curve_Bezier(nc,t, x,dx,ddx,dddx)
!!!        else
!!!          call curve_Bezier(abs(nc),(1.d0-t), x,dx,ddx,dddx)
!!!          dx = -dx; dddx = -dddx
!!!        endif
!!!!!        rpsi = 0.d0; rbeta = 0.d0; rdpsi = 0.d0; rdbeta = 0.d0
!!!!!        rddpsi = 0.d0; rddbeta = 0.d0
!!!!!        call Hshape1(3,t, shapef,dshapef,ddshapef)
!!!!!        do j = 1,4
!!!!!          rpsi    = rpsi   + psi(1:3,j,ie)*shapef(j)
!!!!!          rdpsi   = rdpsi  + psi(1:3,j,ie)*dshapef(j)
!!!!!          rddpsi  = rddpsi  + psi(1:3,j,ie)*ddshapef(j)
!!!!!
!!!!!          rbeta   = rbeta  + beta(j,ie)*shapef(j)
!!!!!          rdbeta  = rdbeta + beta(j,ie)*dshapef(j)
!!!!!          rddbeta = rddbeta + beta(j,ie)*ddshapef(j)
!!!!!        enddo
!!!!!        der = (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dx + rbeta*rpsi
!!!!!        dder = (alpha(2,ie)-alpha(1,ie))*dx + (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*ddx + &
!!!!!               rdbeta*rpsi + rbeta*rdpsi
!!!!!        ddder = 2.d0*(alpha(2,ie)-alpha(1,ie))*ddx + (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dddx + &
!!!!!                rddbeta*rpsi + 2.d0*rdbeta*rdpsi + rbeta*rddpsi
!!!!
!!!
!!!
!!!!!        temp = vel(1:3,1,ie) - dx
!!!!!        call norm(temp, s)
!!!!!        if (s.gt.GEOM_TOL) then
!!!!!          write(*,*)'inconsistency with velocity (1)! ie = ',ie
!!!!!          write(*,*)'x    = ',x
!!!!!          write(*,*)'dx   = ',dx
!!!!!          write(*,*)'ddx  = ',ddx
!!!!!          write(*,*)'dddx = ',dddx
!!!!!          write(*,*)'vel  = ',vel(1:3,1,ie)
!!!!!          call pause
!!!!!        endif
!!!
!!!
!!!
!!!        der = alpha(1,ie)*dx + beta(1,ie)*psi(1:3,1,ie)
!!!        call compute_radial_derivative_new(Nt,ie,t, r_der,r_dder,r_ddder)
!!!!
!!!        temp = der - r_der
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,*)'ie = ',ie
!!!          write(*,3555)t,der
!!!! 3555       format('   t = ',e12.5,' --> p_der   = ',3(e12.5,2x))
!!!          write(*,992)t,r_der
!!!! 992        format('   t = ',e12.5,' --> r_der   = ',3(e12.5,2x))
!!!          call pause
!!!        endif
!!!!
!!!!       2 ND      E  N  D  P  O  I  N  T
!!!        t = 1.d0
!!!        call compute_radial_derivative_new(Nt,ie,t, r_der,r_dder,r_ddder)
!!!!  .....check radial derivative 1st method
!!!        temp = r_der - vel(1:3,1,imod(ie+1,3))
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,8998)ie,t,s1
!!!          call pause
!!!        endif
!!!        select case(ie)
!!!        case(1)
!!!          k1 = bijec(4,1); k2 = bijec(5,0); k3 = bijec(5,1); k4 = bijec(6,0)
!!!        case(2)
!!!          k1 = bijec(1,4); k2 = bijec(1,5); k3 = bijec(0,5); k4 = bijec(0,6)
!!!        case(3)
!!!          k1 = bijec(1,1); k2 = bijec(0,1); k3 = bijec(1,0); k4 = bijec(0,0)
!!!        endselect
!!!!  .....check radial derivative 2nd method
!!!        temp = 6.d0*(TRIANGLES(Nt)%Rdata(k3:k3+2) - TRIANGLES(Nt)%Rdata(k4:k4+2))
!!!        temp = temp - r_der
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,8995)ie,t,s1
!!!          call pause
!!!        endif
!!!!  .....check derivative of radial derivative
!!!        temp = - 30.d0*(TRIANGLES(Nt)%Rdata(k1:k1+2) - TRIANGLES(Nt)%Rdata(k2:k2+2) -   &
!!!                        TRIANGLES(Nt)%Rdata(k3:k3+2) + TRIANGLES(Nt)%Rdata(k4:k4+2)) +  &
!!!                  6.d0*(TRIANGLES(Nt)%Rdata(k2:k2+2) - TRIANGLES(Nt)%Rdata(k4:k4+2))
!!!        temp = temp - r_dder
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,8996)ie,t,s1
!!!          call pause
!!!        endif
!!!!  .....check agreement with desired radial derivative
!!!        nc = TRIANGLES(Nt)%EdgeNo(ie)
!!!        if (nc.gt.0) then
!!!          call curve_Bezier(nc,t, x,dx,ddx,dddx)
!!!        else
!!!          call curve_Bezier(abs(nc),1.d0-t, x,dx,ddx,dddx)
!!!          dx = -dx; dddx = -dddx
!!!        endif
!!!        rpsi = 0.d0; rbeta = 0.d0; rdpsi = 0.d0; rdbeta = 0.d0
!!!        rddpsi = 0.d0; rddbeta = 0.d0
!!!        call Hshape1(3,t, shapef,dshapef,ddshapef)
!!!        do j = 1,4
!!!          rpsi    = rpsi   + psi(1:3,j,ie)*shapef(j)
!!!          rdpsi   = rdpsi  + psi(1:3,j,ie)*dshapef(j)
!!!          rddpsi  = rddpsi  + psi(1:3,j,ie)*ddshapef(j)
!!!
!!!          rbeta   = rbeta  + beta(j,ie)*shapef(j)
!!!          rdbeta  = rdbeta + beta(j,ie)*dshapef(j)
!!!          rddbeta = rddbeta + beta(j,ie)*ddshapef(j)
!!!        enddo
!!!        der = (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dx + rbeta*rpsi
!!!        dder = (alpha(2,ie)-alpha(1,ie))*dx + (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*ddx + &
!!!               rdbeta*rpsi + rbeta*rdpsi
!!!        ddder = 2.d0*(alpha(2,ie)-alpha(1,ie))*ddx + (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dddx + &
!!!                rddbeta*rpsi + 2.d0*rdbeta*rdpsi + rbeta*rddpsi
!!!!
!!!        call compute_radial_derivative_new(Nt,ie,t, r_der,r_dder,r_ddder)
!!!!
!!!        temp = der - r_der
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,*)'ie = ',ie
!!!          write(*,3555)t,der
!!!! 3555       format('   t = ',e12.5,' --> p_der   = ',3(e12.5,2x))
!!!          write(*,992)t,r_der
!!!! 992        format('   t = ',e12.5,' --> r_der   = ',3(e12.5,2x))
!!!          call pause
!!!        endif
!!!!
!!!        temp = dder - r_dder
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,*)'ie = ',ie
!!!          write(*,3556)t,dder
!!!! 3556       format('   t = ',e12.5,' --> p_dder  = ',3(e12.5,2x))
!!!          write(*,991)t,r_dder
!!!! 991        format('   t = ',e12.5,' --> r_dder  = ',3(e12.5,2x))
!!!          call pause
!!!        endif
!!!!
!!!        temp = ddder - r_ddder
!!!        call norm(temp, s1)
!!!        if (s1.gt.GEOM_TOL) then
!!!          write(*,*)'ie = ',ie
!!!          write(*,3557)t,ddder
!!!! 3557       format('   t = ',e12.5,' --> p_ddder = ',3(e12.5,2x))
!!!          write(*,990)t,r_ddder
!!!! 990        format('   t = ',e12.5,' --> r_ddder = ',3(e12.5,2x))
!!!          call pause
!!!        endif
!!!
!!!
!!!!
!!!      enddo



!
!  ...check matching b/w desired radial derivative and actual radial derivative
!!!      write(*,*)'generate_Bezier_triangle: checking radial derivative...'
      nsub = 10
!  ...loop over edges
      do ie = 1,3
!  .....get curve number and endpoints
        nc = TRIANGLES(Nt)%EdgeNo(ie)
        np1 = CURVES(abs(nc))%EndPoNo(1);  np2 = CURVES(abs(nc))%EndPoNo(2)
!!!        write(*,*)'*************************************************************************'
!!!        write(*,1001)ie,nc
!  .....loop over subdivisions
        do i = 0,nsub
          t = i*1.d0/nsub
!  .......1st edge determines triangle polarization
          if (ie.eq.1) then
            if (nc.gt.0) then
              call curve_Bezier(nc,t, x,dx,ddx,dddx)
!  ...........update endpoint if last visit
              if (i.eq.nsub) np = np2
            else
              call curve_Bezier(abs(nc),1.d0-t, x,dx,ddx,dddx)
              dx = -dx; dddx = -dddx
!  ...........update endpoint if last visit
              if (i.eq.nsub) np = np1
            endif
!  .......2nd and 3rd curves
          else
            if (np1.eq.np) then
              call curve_Bezier(abs(nc),t, x,dx,ddx,dddx)
              if (i.eq.nsub) np = np2
            else
              call curve_Bezier(abs(nc),1.d0-t, x,dx,ddx,dddx)
              dx = -dx;  dddx = -dddx
              if (i.eq.nsub) np = np1
            endif
          endif

          call Hshape1(3,t, shapef,dshapef,ddshapef)
          rpsi = 0.d0; rbeta = 0.d0; rdpsi = 0.d0; rdbeta = 0.d0
          rddpsi = 0.d0; rddbeta = 0.d0
          do j = 1,4
            rpsi    = rpsi    + psi(1:3,j,ie)*shapef(j)
            rdpsi   = rdpsi   + psi(1:3,j,ie)*dshapef(j)
            rddpsi  = rddpsi  + psi(1:3,j,ie)*ddshapef(j)

            rbeta   = rbeta   + beta(j,ie)*shapef(j)
            rdbeta  = rdbeta  + beta(j,ie)*dshapef(j)
            rddbeta = rddbeta + beta(j,ie)*ddshapef(j)
          enddo
          der   = (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dx + rbeta*rpsi
          dder  = (alpha(2,ie)-alpha(1,ie))*dx + (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*ddx + &
                   rdbeta*rpsi + rbeta*rdpsi
          ddder = 2.d0*(alpha(2,ie)-alpha(1,ie))*ddx + (alpha(1,ie)*(1.d0-t)+alpha(2,ie)*t)*dddx + &
                  rddbeta*rpsi + 2.d0*rdbeta*rdpsi + rbeta*rddpsi

!!!          write(*,6051)alpha(1:2,ie)
!!! 6051     format(' alpha   = ',2(e12.5,2x))
!!!          write(*,6050)beta(1:4,ie)
!!! 6050     format(' beta    = ',4(e12.5,2x))
!!!          write(*,6052)psi(1:3,1,ie)
!!! 6052     format(' psi(0)  = ',3(e12.5,2x))
!!!          write(*,6053)psi(1:3,2,ie)
!!! 6053     format(' psi(1)  = ',3(e12.5,2x))
!!!          write(*,6054)psi(1:3,3,ie)
!!! 6054     format(' dpsi(0) = ',3(e12.5,2x))
!!!          write(*,6055)psi(1:3,4,ie)
!!! 6055     format(' dpsi(1) = ',3(e12.5,2x))
!!!          write(*,6056)vel(1:3,1,ie)
!!! 6056     format(' vel(0)  = ',3(e12.5,2x))
!!!          write(*,6057)dx
!!! 6057     format(' dx      = ',3(e12.5,2x))
!!!          write(*,6058)acc(1:3,1,ie)
!!! 6058     format(' acc(0)  = ',3(e12.5,2x))
!!!          write(*,6059)ddx
!!! 6059     format(' ddx     = ',3(e12.5,2x))
!!!          write(*,6060)dac(1:3,1,ie)
!!! 6060     format(' dac(0)  = ',3(e12.5,2x))
!!!          write(*,6061)dddx
!!! 6061     format(' dddx    = ',3(e12.5,2x))
!!!
!!!
!!!
!!!          write(*,6005)rbeta
!!! 6005     format(' rbeta   = ',e12.5)
!!!          write(*,6009)rdbeta
!!! 6009     format(' rdbeta  = ',e12.5)
!!!          write(*,6010)rddbeta
!!! 6010     format(' rddbeta = ',e12.5)
!!!          write(*,6006)rpsi
!!! 6006     format(' rpsi    = ',3(e12.5,2x))
!!!          write(*,6007)rdpsi
!!! 6007     format(' rdpsi   = ',3(e12.5,2x))
!!!          write(*,6008)rddpsi
!!! 6008     format(' rddpsi  = ',3(e12.5,2x))
!!!          call pause
!!!


          call compute_radial_derivative(Nt,ie,t, r_der,r_dder,r_ddder)
!
          temp = der - r_der
          call norm(temp, s1)
          if (s1.gt.GEOM_TOL) then
            write(*,3555)t,der
 3555       format('   t = ',e12.5,' --> p_der   = ',3(e12.5,2x))
            write(*,992)t,r_der
 992        format('   t = ',e12.5,' --> r_der   = ',3(e12.5,2x))
          endif
!
          temp = dder - r_dder
          call norm(temp, s1)
          if (s1.gt.GEOM_TOL) then
            write(*,3556)t,dder
 3556       format('   t = ',e12.5,' --> p_dder  = ',3(e12.5,2x))
            write(*,991)t,r_dder
 991        format('   t = ',e12.5,' --> r_dder  = ',3(e12.5,2x))
          endif
!
          temp = ddder - r_ddder
          call norm(temp, s1)
          if (s1.gt.GEOM_TOL) then
            write(*,3557)t,ddder
 3557       format('   t = ',e12.5,' --> p_ddder = ',3(e12.5,2x))
            write(*,990)t,r_ddder
 990        format('   t = ',e12.5,' --> r_ddder = ',3(e12.5,2x))
          endif
!  .....end of loop over subdivisions
        enddo
!!!        call pause
!  ...end of loop over edges
      enddo

!      call pause
      iflag = 0
      call trian(Nt,(/0.d0,0.d0/), x,dx_deta)
      call cross_product(dx_deta(1:3,1),dx_deta(1:3,2), dx_deta(1:3,3))
      call normalize(dx_deta(1:3,3))
      call cross_product(dx_deta(1:3,3),nor(1:3,1), dx_deta(1:3,1))
      call norm(dx_deta(1:3,1), s)
      if (s.gt.GEOM_TOL) then
        write(*,*)'generate_Bezier_triangle: inconsistency with normal (1)!'
        iflag = 1
      endif
!
      call trian(Nt,(/1.d0,0.d0/), x,dx_deta)
      call cross_product(dx_deta(1:3,1),dx_deta(1:3,2), dx_deta(1:3,3))
      call normalize(dx_deta(1:3,3))
      call cross_product(dx_deta(1:3,3),nor(1:3,2), dx_deta(1:3,1))
      call norm(dx_deta(1:3,1), s)
      if (s.gt.GEOM_TOL) then
        write(*,*)'generate_Bezier_triangle: inconsistency with normal (2)!'
        iflag = 1
      endif
!
      call trian(Nt,(/0.d0,1.d0/), x,dx_deta)
      call cross_product(dx_deta(1:3,1),dx_deta(1:3,2), dx_deta(1:3,3))
      call normalize(dx_deta(1:3,3))
      call cross_product(dx_deta(1:3,3),nor(1:3,3), dx_deta(1:3,1))
      call norm(dx_deta(1:3,1), s)
      if (s.gt.GEOM_TOL) then
        write(*,*)'generate_Bezier_triangle: inconsistency with normal (3)!'
        iflag = 1
      endif
      if (iflag.eq.1) then
        write(*,*)'Nt = ',Nt
        call print_GMP
      endif
!
      end subroutine generate_Bezier_triangle
