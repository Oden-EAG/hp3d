c----------------------------------------------------------------------
c
c   routine name       - reconstruct
c
c----------------------------------------------------------------------
c
c   latest revision    - Apr 10
c
c   purpose            - routine performs G^1 reconstruction of
c                        geometry for a triangulated surface
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine reconstruct
c
      use GMP
      use control
c
#include "syscom.blk"
c
      integer,parameter :: nin = 15
c
      common /crhsub/ iprint_rhsub
      common /ctri/ iprint_tri
      common /creconstruct_curve/ iprint_reconstruct_curve
      logical lexist
c
c  ...data for curve reconstruction, points and normals
      dimension x(3,2)
c
c  ...lists of curves adjacent to a point
      parameter(maxc=15)
      dimension ncurve_list1(maxc), ncurve_list2(maxc)
c
c  ...corresponding list of triangles, vertex number, orientation
      dimension ntria_list(1:3,maxc)
c
c  ...work space
      dimension temp(1:3),rn1(3),rn2(3),dnorm(1:3,2),rnorm(3)
      integer, dimension(200) :: curve_fail
c
c  ...tangent vectors and their derivatives at endpoints of a curve
      dimension e(3,2), eprim(3,2)
c
c  ...for checking continuity of derivatives
      dimension rnorm1(1:3,0:10),rnorm2(1:3,0:10),
     .          eta(3),dxdeta(3,3),void(3)
c
c  ...data for checking continuity of normals
      dimension nr_save(2),nedg_save(2),nblock_save(2)
c
c  ...minimum admissible solid angle
      data pi /3.14159265358979312d0/
      angle_min = 1.d0/9.d0*pi/2.d0
      angle_max = pi - angle_min
c
c----------------------------------------------------------------------
c
      iprint=0
c
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: DEBUGGING...............'
      endif
c
c  ...determine if FINE GRID is available
      inquire(file='./files/fine_grid' , exist=lexist)
      if (.not.lexist) then
        ifine_grid = 0
        write(*,*)'reconstruct: no fine grid detected'
        call pause
      else
        open(unit=nin,file='./files/fine_grid',form='formatted',
     .       access='sequential',status='unknown')
        ifine_grid = 1
        write(*,*)'reconstruct: has opened fine grid file'
      endif
c
c  ...read FINE GRID file
      if (ifine_grid.eq.1) then
c  .....read in # number of points in the fine grid
        read(nin,*)NRPO_FG
        write(*,*)'reconstruct: NRPO_FG = ',NRPO_FG
        ! call pause
        if (NRPO_FG.gt.MAXNP_FG) then
          write(*,*)'reconstruct: increase MAXNP_FG,NRPO_FG = ',NRPO_FG
          stop
        endif
c  .....read in fine grid scale
        read(nin,*)SCALE_FG
        write(*,*)'reconstruct: SCALE_FG = ',SCALE_FG
        ! call pause
c  .....read in points in the fine grid
        do i = 1,NRPO_FG
          allocate(POINTS_FG(i)%Rdata(3), STAT=is)
          if (is.ne.0) then
            write(*,*)'reconstruct: Rdata not allocated for i = ',i
            stop
          endif
          allocate(POINTS_FG(i)%Idata(1), STAT=is)
          if (is.ne.0) then
            write(*,*)'reconstruct: Idata not allocated for i = ',i
            stop
          endif
          read(nin,*)POINTS_FG(i)%Idata(1),POINTS_FG(i)%Rdata(1:3)
        enddo
c  .....close file
        close(nin)
      endif
c
c*********************************************************************
c
c  ...Step 0: Verification of regularity conditions necessary for
c             geometry reconstruction
c
c----------------------------------------------------------------------
c
c  STEP 0a: redefine triangles on the reconstructed surface; a triangle
c     is assumed to lie on the reconstructed surface if it satifies 2
c     conditions: it is on A SURFACE, i.e. it has only one neigboring
c     block, and ALL of its three vertices have been marked as vertex-
c     normal points
c
      if (iprint.ge.1) then
        write(*,*)'reconstruct: redefining triangles...'
      endif

c======================================================================
c  REMARK: this check does not work for sharp edges, i.e. a triangle  |
c     with 3 'CoorNrm' vertices, but inside the domain               |
c======================================================================
cccc  ...loop over triangles
ccc      do nr = 1,NRTRIAN
ccc        ndec = 0
cccc  .....loop over vertices
ccc        do iv = 1,3
ccc          np = TRIANGLES(nr)%VertNo(iv)
ccc          if (POINTS(np)%Type.eq.'CoorNrm') ndec = ndec+1
ccc        enddo
ccc        ilen = 0
cccc  .....loop over adjacent blocks
ccc        do is = 1,2
ccc          if (TRIANGLES(nr)%BlockNo(is).ne.0) ilen = ilen+1
ccc        enddo
cccc  .....select # of adjacent blocks
ccc        select case(ilen)
ccc        case(0)
ccc          write(*,7001)nr
ccc 7001     format('reconstruct: HAVE ENCOUNTERED AN ISOLATED
ccc     . TRIANGLE, WITH NO NEIGHBORING BLOCKS nr = ',i6)
ccc          call pause
ccc        case(1)
cccc  .......incremet ndec, i.e. ndec = 4
ccc          ndec = ndec+1
ccc        case(2)
cccc  .......if all vertces are 'CoorNrm'
ccc          if (ndec.eq.3) then
ccc            write(*,7002) nr
ccc 7002       format('reconstruct: TRIANGLE nr = ',i6,' SPANS
ccc     . ACROSS THE DOMAIN')
ccc          endif
cccc  .....end select # of adjacent blocks
ccc        end select
cccc
cccc  .....select number of 'CoorNrm' vertices + 1
ccc        select case(ndec)
ccc        case(2,3)
cccc  .......use the transfinite interpolation
ccc          TRIANGLES(nr)%Type = 'TransTri'
ccc        case(4)
cccc  .......use the G^1 triangle
ccc          if (iprint.ge.1) then
ccc            write(*,7003) nr
ccc 7003       format(' reconstruct: G1 TRIANGLE nr = ',i6)
ccc          endif
ccc          TRIANGLES(nr)%Type = 'G1RecTri'
ccc        end select
cccc  ...end of loop over triangles
ccc      enddo
c
c  ...printing
      if (iprint.ge.1) then
        write(*,*)'reconstruct: triangles redefined.'
      endif
c
c----------------------------------------------------------------------
c
c  STEP 0b: redefine the curves on the reconstructed surface; a curve
c     is defined to be Hermitian if it has EXACTLY two neighboring
c     G1 triangles, an additional check involves verifying that
c     both endpoints are point-normal points. Curves will be redefined
c     as Bezier curves later.
c

c     NOT REALLY NEEDED....


      if (iprint.ge.1) then
        write(*,*)'reconstruct: redefining curves...'
      endif
c  ...loop over curves
      do nc = 1,NRCURVE
c  .....collect the data for the curve
        nflag = 0
c  .....loop over curve endpoints
        do i = 1,2
          np = CURVES(nc)%EndPoNo(i)
          x(1:3,i) = POINTS(np)%Rdata(1:3)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag = nflag + 1
          endif
c  .....end of loop over endpoints
        enddo
c  .....loop over connected figures
        do i = 1,CURVES(nc)%NrFig
          call decode(iabs(CURVES(nc)%FigNo(i)), nr,lab)
c  .......if attached figure is a triangle
          if (lab.eq.1) then
            if (TRIANGLES(nr)%Type.eq.'G1RecTri') nflag = nflag + 1
          endif
c  .....end of loop over figures
        enddo
c  .....Hermitian curve has 2 'CoorNrm' vertices and 2 'G1RecTri'
        if (nflag.eq.4) then
          if (iprint.ge.2) then
            write(*,7004) nc
 7004       format(' reconstruct: DEFINING HERMITIAN CURVE
     . nc = ',i6)
          endif
          CURVES(nc)%Type='HermCur'
          allocate(CURVES(nc)%Rdata(12), STAT=is)
          if (is.ne.0) then
            write(*,9001)nc
 9001     format(' reconstruct: Rdata not allocated for nc = ',i5)
            stop
          endif
          CURVES(nc)%Rdata(1:3) = x(1:3,2) - x(1:3,1)
          CURVES(nc)%Rdata(4:6) = x(1:3,2) - x(1:3,1)
          CURVES(nc)%Rdata(7:12) = 0.d0
        elseif (nflag.gt.4) then
          write(*,7005) nc,nflag
 7005     format('reconstruct: WRONG TOPOLOGY ! nc,nflag = ',2i6)
          stop
        endif
c  ...end of loop over curves
      enddo
      if (iprint.ge.1) then
        write(*,*)'reconstruct: curves redefined.'
      endif
c
c
c----------------------------------------------------------------------
c
c  STEP 1: Construction of tangent bundles
c
      if (iprint.ge.1) then
        write(*,*)'reconstruct: constructing tangent bundles...'
      endif
c  ...loop over points
      do np = 1,NRPOINT
c  .....skip if not a point with the normal
        if (POINTS(np)%Type.ne.'CoorNrm') cycle
        if (POINTS(np)%NrCurv.eq.0) cycle
        if (iprint.ge.2) then
          write(*,8010)POINTS(np)%Rdata(1:6)
 8010     format(' reconstruct: point, normal = ',3f8.3,2x,3f8.3)
        endif
c  .....collect all Hermite curves connected to the point
        k = 0
        do j = 1,POINTS(np)%NrCurv
          nc = POINTS(np)%CurvNo(j)
          if (CURVES(nc)%Type.eq.'HermCur') then
            k=k+1
            if (k.gt.maxc) then
              write(*,*) 'reconstruct: INCREASE maxc'
              stop
            endif
            ncurve_list1(k) = nc
          endif
        enddo
        nr_connected_curves = k
        if (iprint.ge.2) then
          write(*,8005) np,ncurve_list1(1:nr_connected_curves)
 8005     format(' reconstruct: np,ncurve_list1 = ',
     .            i5,2x,10i5)
          call pause
        endif
c
c
c=======================================================================
c  .....for the point:
c         determine UNIT normal and store it in data structure
c  .....for each curve:
c         determine UNIT tangent vectors and principal normals rescaled
c         by curvature and store them in data structure
        call set_tangents(np,ncurve_list1,nr_connected_curves,1,
     .                    ncurve_list2,ntria_list,
     .                    ifirst,isecond,rn1,rn2)
c=======================================================================
c  REMARK: so far, for each curve, ONLY the UNIT tangent vectors at the
c  endpoints have been determined! The actual velocities at the end-
c  points are still unknown! They will be determined in STEP 2.
c=======================================================================
c
c
c  ...end of loop over points
      enddo
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: tangent bundles reconstructed.'
      endif
c
c----------------------------------------------------------------------
c
c  STEP 2: Reconstruct curves using unit tangent vectors and local
c     curvatures
c


ccc      iprint=1

      if (iprint.ge.1) then
        write(*,*)'reconstruct: reconstructing curves...'
      endif
c  ...loop over all curves
      do nc = 1,NRCURVE
c  .....skip if not a Hermitian curve
        if (CURVES(nc)%Type.ne.'HermCur') cycle
c  .....collect the data for the curve
        do i=1,2
          np = CURVES(nc)%EndPoNo(i)
          x(1:3,i)  = POINTS(np)%Rdata(1:3)
        enddo
c  .....unit tangents
        e(1:3,1) = CURVES(nc)%Rdata(1:3)
        e(1:3,2) = CURVES(nc)%Rdata(4:6)
c  .....principal normals rescaled by the curvature
        dnorm(1:3,1) = CURVES(nc)%Rdata(7:9)
        dnorm(1:3,2) = CURVES(nc)%Rdata(10:12)
c  .....compute velocities (e) and accelerations (eprim)
        call reconstruct_curve(x,e,dnorm,eprim,ii)
c  .....store velocities and accelerations
        CURVES(nc)%Rdata(1:3)   = e(1:3,1)
        CURVES(nc)%Rdata(4:6)   = e(1:3,2)
        CURVES(nc)%Rdata(7:9)   = eprim(1:3,1)
        CURVES(nc)%Rdata(10:12) = eprim(1:3,2)
c  .....check error flag
        if (ii.ne.0) then
          write(*,7039) nc
 7039     format('reconstruct: REVERSED TANGENT FOR nc = ',i6)
          nrfig = CURVES(nc)%NrFig
          write(*,7040) iabs(CURVES(nc)%FigNo(1:nrfig))/10
 7040     format('             ADJACENT TRIANGLES = ',10i6)
          call pause
        endif
c  .....printing
        if (iprint.ge.2) then
          write(*,8001) nc
 8001    format(' reconstruct: HAVE RECONSTRUCTED CURVE nc = ',i5)
          write(*,8002) x(1:3,1),x(1:3,2)
 8002     format(' reconstruct: ENDPOINTS       = ',2(3f8.3,2x))
          write(*,8004) CURVES(nc)%Rdata(1:6)
 8004     format(' reconstruct: TANGENT VECTORS = ',2(3f8.3,2x))
          write(*,8014) CURVES(nc)%Rdata(7:12)
 8014     format(' reconstruct: 2nd DER VECTORS = ',
     .            2(3f8.3,2x))
          call pause
        endif

c
cc        iprint=0
cc  .....collect data for checking consistency b/w Hermite and Bezier curve
c        if (iprint.ge.1) then
c          write(*,9006)nc
c 9006    format(' reconstruct: checking consistency for nc = ',i4)
c          i_sub = 10
c          do i = 0,i_sub
c            r_eta = 1.d0/i_sub
c            call curve_HermCur(nc,r_eta, r_work(i,1:3),r_work(i,4:6))
c          enddo
c        endif
cc
cc  DON'T GENERATE Bezier CURVES JUST YET
c
cc  .....determine Bezier curve control points
c        call generate_Bezier_curve(nc,x,e,eprim)
cc
cc  .....check consistency
c        if (iprint.ge.1) then
c          do i = 0,i_sub
c            r_eta = 1.d0/i_sub
c            call curve_Bezier(nc,r_eta, r_x,r_dx)
c            call norm((r_x(1:3)-r_work(i,1:3)), r_norm)
c            write(*,9007)i,r_work(i,1:3)
c 9007       format(' i = ',i2,'   Hermite = ',3(e12.5,2x))
c            write(*,9008)i,r_x(1:3)
c 9008       format(' i = ',i2,'   Bezier  = ',3(e12.5,2x))
c            write(*,9009)i,r_norm
c 9009       format(' i = ',i2,'   diff    = ',e12.5)
c            if (r_norm.gt.GEOM_TOL) then
c              write(*,9010)i
c 9010         format(' i = ',i2,'   GEOM_TOL exceeded ')
c            endif
c          enddo
c          call pause
c          write(*,*)'--------------------------------------------------'
c        endif
cc
c  ...end of loop through curves
      enddo
c
c  ...printing
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: curves reconstructed.'
        call pause
      endif
c
c
c======================================================================
c  REMARK: right now curves have been fully determined, i.e. the
c  velocities and accelerations at the endpoints have been computed!!!
c
c  For checking the wire frame activate the following line:
ccc      call graphg
c======================================================================
c
c
c----------------------------------------------------------------------
c
c  STEP 3: compute  normal components of mixed derivative at vertices
c     of G^1 triangles and store them in datastructure
c
c  ...allocate Rdata for G1RecTri (36 control points)
      do i = 1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'G1RecTri') cycle
        allocate(TRIANGLES(i)%Rdata(0:107),STAT=is)
        if (is.ne.0) then
          write(*,1111)i
 1111     format(' reconstruct: Rdata not allocated for triangle
     . = ',i4)
          stop
        endif
      enddo
c
ccc      iprint=0
      if (iprint.ge.1) then
        write(*,995)
 995    format(' reconstruct: determining normal component of
     . mixed derivative...')
      endif
c
c  ...loop over points
      do np = 1,NRPOINT
c  .....skip if not a point with the normal
        if (POINTS(np)%Type.ne.'CoorNrm') cycle
        if (POINTS(np)%NrCurv.eq.0) cycle
        if (iprint.eq.3) then
          write(*,8010)POINTS(np)%Rdata(1:6)
        endif
c  .....collect all Hermite curves connected to the point
        k = 0
        do j = 1,POINTS(np)%NrCurv
          nc = POINTS(np)%CurvNo(j)
          if (CURVES(nc)%Type.eq.'HermCur') then
            k = k + 1
            if (k.gt.maxc) then
              write(*,*)'reconstruct: INCREASE maxc'
              stop
            endif
            ncurve_list1(k) = nc
          endif
        enddo
        nr_connected_curves = k
        if (iprint.eq.3) then
          write(*,8005) np,ncurve_list1(1:nr_connected_curves)
          call pause
        endif
c  .....get the patch data
        call set_tangents(np,ncurve_list1,nr_connected_curves,0,
     .                    ncurve_list2,ntria_list,
     .                    ifirst,isecond,rn1,rn2)
c
c  .....if no sharp edge was found
        if (ifirst.eq.0) then
          call reconstruct_patch(np,ncurve_list2,nr_connected_curves,0,
     .                            POINTS(np)%Rdata(4:6),ntria_list)
c
c  .....if a sharp edge was found
        else
c  .......1st    S  U  B  P  A  T  C  H
          call reconstruct_patch(np,ncurve_list2(1:isecond),isecond,1,
     .                            rn1,ntria_list(1:3,1:isecond))
c  .......2nd    S  U  B  P  A  T  C  H
          nr_curv = nr_connected_curves
          call reconstruct_patch(np,ncurve_list2(isecond:nr_curv),
     .                           (nr_curv-isecond+1),1,rn2,
     .                            ntria_list(1:3,isecond:nr_curv))
        endif
c
c  ...end of loop over points
      enddo
c
      if (iprint.ge.1) then
        write(*,996)
 996    format(' reconstruct: normal component of
     . mixed derivative determined')
      endif
c
c-------------------------------------------------------------------------
c
c  STEP 5: generate quintic Bezier curves
c
      if (iprint.ge.1) then
        write(*,*)'reconstruct: generating Bezier curves...'
      endif
c  ...loop over curves
      do nc = 1,NRCURVE
c  .....skip if not a Hermitian curve
        if (CURVES(nc)%Type.ne.'HermCur') cycle
c  .....endpoints
        do i = 1,2
          np = CURVES(nc)%EndPoNo(i)
          x(1:3,i)  = POINTS(np)%Rdata(1:3)
        enddo
c  .....velocities
        e(1:3,1) = CURVES(nc)%Rdata(1:3)
        e(1:3,2) = CURVES(nc)%Rdata(4:6)
c  .....accelerations
        eprim(1:3,1) = CURVES(nc)%Rdata(7:9)
        eprim(1:3,2) = CURVES(nc)%Rdata(10:12)
c  .....determine Bezier curve control points
        call generate_Bezier_curve(nc,x,e,eprim)
c  ...end of loop over curves
      enddo
      if (iprint.ge.1) then
        write(*,*)'reconstruct: have generated Bezier curves'
      endif
c
c-------------------------------------------------------------------------
c
c  STEP 6: generate sextic Bezier triangles
c
      if (iprint.ge.1) then
        write(*,*)'reconstruct: generating Bezier triangles...'
      endif
c
c  ...loop over triangles
      do nt = 1,NRTRIAN
        if (TRIANGLES(nt)%Type.ne.'G1RecTri') cycle
        call generate_Bezier_triangle(nt)
      enddo
c  ...end of loop over triangles
c
      if (iprint.ge.1) then
        write(*,*)'reconstruct: have generated Bezier triangles'
      endif


ccc      return

c
c  ...check compatibility b/w G1Rec triangle and edges parameterizations
      do i = 1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'G1RecTri') cycle
        call check_tri2edgs_compatibility(i)
      enddo
c
c-------------------------------------------------------------------------
c
c  CHECK...
c
      nfail = 0
      ntri = 0
c  ...loop over triangles
      do nt = 1,NRTRIAN
        if (TRIANGLES(nt)%Type.ne.'G1RecTri') cycle
        ntri = ntri + 1
        np1 = TRIANGLES(nt)%VertNo(1)
        np2 = TRIANGLES(nt)%VertNo(2)
        np3 = TRIANGLES(nt)%VertNo(3)
c
c  .....determine point normal
        select case(POINTS(np1)%Type)
        case('CoorNrm')
          rnorm(1:3) = POINTS(np1)%Rdata(4:6)
        case('SharpPt')
c  .......set error flag to 0
          iflag = 0
c  .......determine # of triangles in the sharp patch
          n = size(POINTS(np1)%Idata)
c  .......loop over sharp patch triangles
          do i = 1,n
            call decode(POINTS(np1)%Idata(i), nt1,lab)
            if (nt1.eq.nt) then
              select case(lab)
              case(1); rnorm(1:3) = POINTS(np1)%Rdata(4:6)
              case(2); rnorm(1:3) = POINTS(np1)%Rdata(7:9)
              endselect
              iflag = 1
              exit
            endif
          enddo
c  .......if triangle was not found
          if (iflag.eq.0) then
            write(*,*)'reconstruct: geometry inconsistency!'
            write(*,*)'nt,np = ',nt,np1
            call print_GMP
            stop
          endif
        case default
          write(*,*)'reconstruct: inconsistent point type!'
          write(*,*)'nt,jv,np = ',nt,np1
          call print_GMP
          stop
        endselect
        if (iprint.ge.1) then
          write(*,1000)np1,rnorm
 1000     format(' reconstruct : np1 = ',i8,'; rnorm = ',3(e12.5,2x))
        endif
c
c  .....determine triangle normal
        call trian(nt,(/0.d0,0.d0/), void,dxdeta)
        call cross_product(dxdeta(1:3,1),dxdeta(1:3,2),
     .                                    dxdeta(1:3,3))
        call cross_product(dxdeta(1:3,3),rnorm(1:3), dxdeta(1:3,1))
        call norm(dxdeta(1:3,1), s)
        if (s.gt.GEOM_TOL) then
          nfail = nfail + 1
ccc          write(*,*)'reconstruct: INCONSISTENCY!!'
ccc          write(*,*)'np,type = ',np1,POINTS(np1)%Type
        endif
c  ...end of loop over triangles
      enddo
      if (nfail.gt.0) then
        write(*,*)'======================================='
        write(*,*)'nfail = ',nfail,'; ntri = ',ntri
        write(*,*)'======================================='
      call pause
      endif






c
c  Step 7: check the continuity of normals
c
      icheck=1
      iprint=0
c  ...skip to next check if icheck is off
      if (icheck.eq.0) goto 60
c
      nfail = 0
c  ...loop through all curves
      do nc = 1,NRCURVE

        if (nc.eq.10065) then
          iprint=1
        else
          iprint=0
        endif

c  .....skip if not a quintic Bezier curve
        if (CURVES(nc)%Type.ne.'5Bezier') cycle
c  .....skip if curve is a sharpe edge
        call check_edge(nc, idec); if (idec.eq.2) cycle
ccc        write(*,8062) nc
 8062   format(' reconstruct: checking normal along nc = ',i5)
c  .....get the beginning point
        np = CURVES(nc)%EndPoNo(1)
c  .....loop through connected triangles
        ifound = 0
        do ir = 1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(ir)), nr,lab)
c  .......cycle if not a G1RecTri
          if (TRIANGLES(nr)%Type.ne.'G1RecTri') cycle
c  .......loop through edges of the triangle
          do ie = 1,3
            nc1 = iabs(TRIANGLES(nr)%EdgeNo(ie))
            if (nc1.ne.nc) cycle
            ifound = ifound + 1
            nr_save(ifound) = nr; nedg_save(ifound) = ie
c  .........loop through adjacent blocks
            do is = 1,2
              if (TRIANGLES(nr)%BlockNo(is).ne.0) then
                nblock_save(ifound) = TRIANGLES(nr)%BlockNo(is)
              endif
            enddo
c
c  for each subdivision point collect normals
            nsub = 2
c  .........loop over subdivision points
            do i = 0,nsub
              xi = i*1.d0/nsub
              select case(ie)
              case(1)
                if (TRIANGLES(nr)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = xi; eta(2) = 0.d0
              case(2)
                if (TRIANGLES(nr)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = 1.d0 - xi; eta(2) = xi
              case(3)
                if (TRIANGLES(nr)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = 0.d0; eta(2) = xi
              end select
              call trian(nr,eta, void,dxdeta)
              call cross_product(dxdeta(1:3,1),dxdeta(1:3,2), temp)
              call normalize(temp)
c  ...........use 1st subdivision point to determine triangle orientation
              if (i.eq.0) then
                call scalar_product(temp,POINTS(np)%Rdata(4:6), s)
                factor = 1.d0
c  .............sign function returns 1 if s > 0, -1 otherwise
                rsign = sign(factor,s)
              endif
              if (iprint.ge.1) then
         write(*,2001)ifound,nr_save(ifound),nedg_save(ifound),xi,rsign,
     . temp
 2001    format(' ifound = ',i1,'; nt = ',i5,'; ie = ',i1,'; xi = ',
     . e12.5,' --> sign = ',e12.5,'; nor = ',3(e12.5,2x))
              endif
              select case(ifound)
              case(1)
                rnorm1(1:3,i) = rsign*temp(1:3)
              case(2)
                rnorm2(1:3,i) = rsign*temp(1:3)
              case default
                write(*,*) 'reconstruct: INCONSISTENCY '
                stop
              end select
c  .........end of loop over subdivision points
            enddo
c  .......end of loop over triangle edges
          enddo
c  .....end of loop over connected triangles
        enddo
c  .....printing
       if (iprint.ge.1) then
         do i = 0,nsub
           write(*,2045)rnorm1(1:3,i)
 2045      format('   rnorm1 = ',3(e12.5,2x))
           write(*,2046)rnorm2(1:3,i)
 2046      format('   rnorm2 = ',3(e12.5,2x))
           enddo
           call pause
         endif
c
        iflag = 0
c  for each subdivision point compare normals
c  .....loop over subdivisions
        do i = 0,nsub
          if (iprint.ge.1) then
            write(*,2045)rnorm1(1:3,i)
            write(*,2046)rnorm2(1:3,i)
            call pause
          endif
c  .......loop over components
          do j = 1,3
c  .........if geometrical tollerance is exceeded
            if (abs(rnorm2(j,i)-rnorm1(j,i)).gt.GEOM_TOL) then
c  ...........raise flag
              iflag = 1
              write(*,3334)nc,i,rnorm1(1:3,i)
 3334         format(' nc = ',i8,'; sub = ',i2,'; norm1 = ',3(e12.5,2x))
              write(*,3335)rnorm2(1:3,i)
 3335         format('                          norm2 = ',3(e12.5,2x))
cc              write(*,3333)nc,i,j,abs(rnorm2(j,i)-rnorm1(j,i))
cc 3333         format(' nc = ',i3,'; sub = ',i2,'; comp = ',i1,'; error
cc     .  = ',e12.5)
cc              write(*,*) 'ERROR IN COMP j = ',j,
cc     .                    abs(rnorm2(j,i)-rnorm1(j,i))
cc              write(*,8056) np,POINTS(np)%Rdata(4:6)
cc 8056         format('np,rn = ',i5,2x,3f8.3)
cc              write(*,8019) nr_save(1:2),nedg_save(1:2),nblock_save(1:2)
cc 8019         format('rectangles,edges,blocks = ',3(2i7,2x))
cc              write(*,8020) nc,i,rnorm1(1:3,i),rnorm2(1:3,i)
cc 8020         format('reconstruct: nc,i,rnorm1(1:3,i),
cc     . rnorm2(1:3,i) = ',i4,i2,2x,3f8.3,2x,3f8.3)
cc              write(*,8021) CURVES(nc)%EndPoNo(1:2)
cc 8021         format('             endpoints = ',2i5)
cc              call pause
cc              go to 50
            endif
c  .......end of loop over components
          enddo
c  .....end of loop over subdivision points
        enddo
c  .....increment counter and save curve number
        if (iflag.eq.1) then
          nfail = nfail + 1; curve_fail(nfail) = nc
        endif
c 50     continue
cc      call pause
c  ...end of loop through curves
      enddo
c
      if (nfail.eq.0) then
        write(*,*)'G1 OBJECT!!!'
      else
        write(*,*)'=================================================='
        write(*,555)NRCURVE,nfail
 555    format(' NRCURVE = ',i10,'; nfail = ',i10)
        do i = 1,nfail
          write(*,557)i,curve_fail(i),CURVES(curve_fail(i))%EndPoNo(1:2)
 557      format(' i = ',i2,'; nc = ',i7,', endpoints = ',i8,', ',i8)
        enddo
        write(*,*)'=================================================='
      endif

      write(*,*) 'reconstruct: HAVE CHECKED CONTINUITY'
c
 60   continue
c
c----------------------------------------------------------------------
c
c  ...Step 8: Check the 3D jacobians (solid angles)
c
      icheck=0
c  ...skip to next check if icheck if off
      if (icheck.eq.0) goto 90
c
      write(*,*)'reconstruct: checking Jacobians...'
c
      fmin = 1.0e30
      nfail = 0
      nsub = 5
c  ...loop over tetras
      do nh = 1,NRTETRA
        iflag = 0
c  .....loop over subdivision points
        do k = 0,nsub
          do j = 0,nsub-k
            do i = 0,nsub-k-j
              eta(1) = i*1.d0/nsub
              eta(2) = j*1.d0/nsub
              eta(3) = k*1.d0/nsub
              call tetra(nh,eta, void,dxdeta)
              do ivar = 1,3
                call normalize(dxdeta(1:3,ivar))
              enddo
              call mixed_product(dxdeta(1:3,1),dxdeta(1:3,2),
     .                           dxdeta(1:3,3), angle)
              if (angle.lt.GEOM_TOLL) then
                fmin = min(angle,fmin)
ccc              if ((angle.lt.angle_min).or.(angle.gt.angle_max)) then
                iflag = 1
ccc                write(*,6037)nh
 6037     format(' reconstruct: BAD SOLID ANGLE IN TET = ,',i4)
ccc                write(*,6038)eta
 6038     format('   Eta         = ',3(e12.5,2x))
ccc                write(*,6041)void
 6041     format('   X           = ',3(e12.5,2x))
                do ii = 1,3
ccc                  write(*,6039)ii,dxdeta(1:3,ii)
 6039     format('   dXdEta(:,',i1,') = ',3(e12.5,2x))
                enddo
ccc                write(*,6040)angle
 6040     format('   Jacobian    = ',e12.5)
cc                write(*,6036) nh,eta,angle
cc 6036  format(' reconstruct: BAD SOLID ANGLE IN BLOCK ntet = ',i7,
cc     . ' AT POINT eta,rjac = ',3(e12.5,2x)'; ',e12.5)
              endif
            enddo
          enddo
c  .....end of loop over subdivision points
        enddo
        if (iflag.eq.1) nfail = nfail + 1
c  ...end of loop over tetras
      enddo
c
      if (nfail.gt.0) then
        write(*,*)'=================================================='
        write(*,556)NRTETRA,nfail,fmin
 556    format(' NRTETRA = ',i10,'; nfail = ',i10,'; fmin = ',e12.5)
        call pause
        write(*,*)'=================================================='
      endif


      write(*,*)'reconstruct: Jacobians checked'
c
 90   continue
c
c
      end
