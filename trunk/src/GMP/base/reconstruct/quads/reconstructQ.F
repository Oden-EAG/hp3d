c----------------------------------------------------------------------
c
c   routine name       - reconstructQ
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine performs G^1 reconstruction of
c                        geometry for quads
c
c----------------------------------------------------------------------
c
      subroutine reconstruct
c
      use GMP
      use control
c
#include "syscom.blk"
      common /crhsub/ iprint_rhsub
      common /ctri/ iprint_tri
      common /creconstruct_curve/ iprint_reconstruct_curve
c
c  ...data for curve reconstruction, points and normals
      dimension x(3,2)
c
c  ...lists of curves adjacent to a point
      parameter(maxc=15)
      dimension ncurve_list1(maxc), ncurve_list2(maxc)
c
c  ...corresponding list of rectangles, vertex number, orientation
      dimension nrect_list(1:3,maxc), nrect_list2(1:3,maxc)
c
c  ...curvatures at a point
      dimension acoeff(3,3)
c
c  ...work space
      dimension temp(1:3),rn1(3),rn2(3),dnorm(1:3,2)
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
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
c  ...minimum admissible solid angle
      data pi /3.14159265358979312d0/
      angle_min = 1.d0/9.d0*pi/2.d0
      angle_max = pi - angle_min
c
c**********************************************************************
c
      iprint=0
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: DEBUGGING...'
        call pause
      endif
c
c*********************************************************************
c
c  ...Step 0: Verification of regularity conditions necessary for
c             geometry reconstruction
c
c  ...redefine rectangles on the reconstructed surface; a rectangle
c     is assumed to lie on the reconstructed surface if it satifies two
c     conditions: it is on A SURFACE, i.e. it has only one neigboring
c     block, and ALL of its four vertices have been marked as vertex-
c     normal points
c
      do nr=1,NRRECTA
        ndec=0
        do iv=1,4
          np = RECTANGLES(nr)%VertNo(iv)
          if (POINTS(np)%Type.eq.'CoorNrm') ndec=ndec+1
        enddo
        ilen=0
        do is=1,2
          if (RECTANGLES(nr)%BlockNo(is).ne.0) ilen=ilen+1
        enddo
        select case(ilen)
        case(0)
          write(*,7001)nr
 7001     format('reconstruct: HAVE ENCOUNTERED AN ISOLATED RECTANGLE',
     .           ' WITH NO NEIGHBORING BLOCKS nr = ',i6)
          call pause
        case(1)
          ndec = ndec+1
        case(2)
          if (ndec.eq.4) then
            write(*,7002) nr
 7002       format('reconstruct: RECTANGLE nr = ',i6,' SPANS ACROSS',
     .             ' THE DOMAIN ')
            call pause
          endif
        end select
c
        select case(ndec)
        case(2,3,4)
c
c  .......use the transfinite interpolation
          RECTANGLES(nr)%Type = 'TraQua'
        case(5)
c
c  .......use the G^1 rectangle
          if (iprint.eq.1) then
            write(*,7003) nr
 7003       format('reconstruct: G1 RECTANGLE nr = ',i6)
          endif
          RECTANGLES(nr)%Type = 'HermRec'
          allocate(RECTANGLES(nr)%Rdata(48))
          RECTANGLES(nr)%Rdata(1:48) = 0.d0
        end select
      enddo
c
c  ...redefine the curves on the reconstructed surface; a curve
c     is defined to be Hermitian if it has EXACTLY two neighboring
c     G1 rectangles, an additional check involves verifying that
c     both endpoints are point-normal points
      do nc=1,NRCURVE
c
c  .....collect the data for the curve
        nflag=0
        do i=1,2
          np = CURVES(nc)%EndPoNo(i)
          x(1:3,i) = POINTS(np)%Rdata(1:3)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag=nflag+1
          endif
        enddo
        do i=1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(i)), nr,lab)
          if (lab.eq.2) then
            if (RECTANGLES(nr)%Type.eq.'HermRec') nflag=nflag+1
          endif
        enddo
c
c  .....skip if not a Hermitian curve
        if (nflag.eq.4) then
          if (iprint.eq.1) then
            write(*,7004) nc
 7004       format('reconstruct: DEFINING HERMITIAN CURVE nc = ',i6)
          endif
          CURVES(nc)%Type='HermCur'
          allocate(CURVES(nc)%Rdata(12))
          CURVES(nc)%Rdata(1:3) = x(1:3,2) - x(1:3,1)
          CURVES(nc)%Rdata(4:6) = x(1:3,2) - x(1:3,1)
          CURVES(nc)%Rdata(7:12) = 0.d0
        elseif(nflag.gt.4) then
          write(*,7005) nc,nflag
 7005     format('reconstruct: WRONG TOPOLOGY ! nc,nflag = ',2i6)
          stop 1
        endif
      enddo
ccc      call print_GMP
ccc      stop
c
c*********************************************************************
c
c  ...Step 1: Construction of tangent bundles
c
c  ...loop through points and set up the tangents
      do np=1,NRPOINT
c
c  .....skip if not a point with the normal
        if (POINTS(np)%Type.ne.'CoorNrm') cycle
        if (POINTS(np)%NrCurv.eq.0) cycle
        if (iprint.eq.1) then
          write(*,8010)POINTS(np)%Rdata(1:6)
 8010     format('reconstruct: POINT, NORMAL = ',3f8.3,2x,3f8.3)
        endif
c
c  .....collect all Hermite curves connected to the point
        k=0
        do j=1,POINTS(np)%NrCurv
          nc = POINTS(np)%CurvNo(j)
          if (CURVES(nc)%Type.eq.'HermCur') then
            k=k+1
            if (k.gt.maxc) then
              write(*,*) 'reconstruct: INCREASE maxc'
              stop 1
            endif
            ncurve_list1(k) = nc
          endif
        enddo
        nr_connected_curves = k
        if (iprint.eq.1) then
          write(*,8005) np,ncurve_list1(1:nr_connected_curves)
 8005     format('reconstruct: np,ncurve_list1 = ',
     .            i5,2x,10i5)
          call pause
        endif
c
c  .....set the tangent bundles
        call set_tangents(np,ncurve_list1, nr_connected_curves,1,
     .                    ncurve_list2,nrect_list,
     .                    ifirst,isecond,rn1,rn2)
c
c  ...end of loop through points
      enddo
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: DONE WITH TANGENT BUNDLES'
      endif
c
c*********************************************************************
c
c  ...Step 2: Reconstruction of curves using the unit tangent
c             vectors and local curvatures
c
c  ...loop through all curves
      do nc=1,NRCURVE
        if (nc.eq.3214) then
          iprint_reconstruct_curve=1
        else
          iprint_reconstruct_curve=0
        endif
c
c  .....skip if not a Hermitian curve
        if (CURVES(nc)%Type.ne.'HermCur') cycle
c
c  .....collect the data for the curve
        do i=1,2
          np = CURVES(nc)%EndPoNo(i)
          x(1:3,i)  = POINTS(np)%Rdata(1:3)
        enddo
        e(1:3,1) = CURVES(nc)%Rdata(1:3)
        e(1:3,2) = CURVES(nc)%Rdata(4:6)
        dnorm(1:3,1) = CURVES(nc)%Rdata(7:9)
        dnorm(1:3,2) = CURVES(nc)%Rdata(10:12)
c
        call reconstruct_curve(x,e,dnorm,eprim,ii)
c
c  .....store the first derivatives
        CURVES(nc)%Rdata(1:3) = e(1:3,1)
        CURVES(nc)%Rdata(4:6) = e(1:3,2)
c
c  .....store the second derivatives
        CURVES(nc)%Rdata(7:9)   = eprim(1:3,1)
        CURVES(nc)%Rdata(10:12) = eprim(1:3,2)
        if (ii.ne.0) then
          write(*,7039) nc
 7039     format('reconstruct: REVERSED TANGENT FOR nc = ',i6)
          nrfig = CURVES(nc)%NrFig
          write(*,7040) abs(CURVES(nc)%FigNo(1:nrfig))/10
 7040     format('             ADJACENT RECTANGLES = ',10i6)
          call pause
        endif
c
        if (iprint.eq.2) then
          write(*,8001) nc
 8001     format('reconstruct: HAVE RECONSTRUCTED CURVE nc = ',i5)
          write(*,8002) x(1:3,1),x(1:3,2)
 8002     format('reconstruct: ENDPOINTS       = ',2(3f8.3,2x))
          write(*,8004) CURVES(nc)%Rdata(1:6)
 8004     format('reconstruct: TANGENT VECTORS = ',2(3f8.3,2x))
          write(*,8014) CURVES(nc)%Rdata(7:12)
 8014     format('reconstruct: SECOND DERIVATIVES VECTORS = ',
     .            2(3f8.3,2x))
          call pause
        endif
c
c  ...end of loop through curves
      enddo
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: HAVE RECONSTRUCTED CURVES'
      endif
      call graphg
c
c**********************************************************************
c
c  ...Step 3: Computation of normal components of mixed derivatives
c             for vertices of G^1 rectangles
c
c  ...loop through points
      do np=1,NRPOINT
c
c  .....skip if not a point with the normal
        if (POINTS(np)%Type.ne.'CoorNrm') cycle
        if (POINTS(np)%NrCurv.eq.0) cycle
        if (iprint.eq.3) then
          write(*,8010)POINTS(np)%Rdata(1:6)
        endif
c
c  .....collect all Hermite curves connected to the point
        k=0
        do j=1,POINTS(np)%NrCurv
          nc = POINTS(np)%CurvNo(j)
          if (CURVES(nc)%Type.eq.'HermCur') then
            k=k+1
            if (k.gt.maxc) then
              write(*,*) 'reconstruct: INCREASE maxc'
              stop 1
            endif
            ncurve_list1(k) = nc
          endif
        enddo
        nr_connected_curves = k
        if (iprint.eq.3) then
          write(*,8005) np,ncurve_list1(1:nr_connected_curves)
          call pause
        endif
c
c  .....get the patch data
        call set_tangents(np,ncurve_list1, nr_connected_curves,0,
     .                    ncurve_list2,nrect_list,
     .                    ifirst,isecond,rn1,rn2)
        if (ifirst.eq.0) then
c
          call reconstruct_patch(np,ncurve_list2,nr_connected_curves,0,
     .                           POINTS(np)%Rdata(4:6),nrect_list)
        else
          nr_curv = isecond-ifirst+1
          call reconstruct_patch(np,ncurve_list2(ifirst),nr_curv,1,
     .                           rn1,nrect_list(1,ifirst))
          ii=0
          do j=isecond,nr_connected_curves+ifirst
            i = imod(j,nr_connected_curves)
            ii=ii+1
            ncurve_list1(ii) = ncurve_list2(i)
            nrect_list2(1:3,ii) = nrect_list(1:3,i)
          enddo
          nr_curv = ii
          call reconstruct_patch(np,ncurve_list1,nr_curv,1,
     .                           rn2,nrect_list2)
        endif
c
c  ...end of loop through points
      enddo
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: HAVE COMPUTED SECOND DERIVATIVES'
      endif
ccc      return
c
c***********************************************************************
c
c  ...Step 4: Check the continuity of normals
c
c  ...loop through all curves
      do nc=1,NRCURVE
c
c  .....skip if not a Hermitian curve
        if (CURVES(nc)%Type.ne.'HermCur') cycle
c
c  .....skip if an edge curve
        call check_edge(nc, idec); if (idec.eq.2) cycle
c
        write(*,8062) nc
 8062   format('reconstruct: CHECKING NORMAL CONTINUITY FOR CURVE ',
     .         'nc = ',i5)
c
c  .....get the beginning point
        np = CURVES(nc)%EndPoNo(1)
c
c  .....loop through connected rectangles
        ifound=0
        do ir=1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(ir)), nr,lab)
          if (RECTANGLES(nr)%Type.ne.'HermRec') cycle
c
c  .......loop through edges of the rectangle
          do ie=1,4
            nc1 = abs(RECTANGLES(nr)%EdgeNo(ie))
            if (nc1.ne.nc) cycle
            ifound = ifound+1
            nr_save(ifound) = nr; nedg_save(ifound) = ie
            do is=1,2
              if (RECTANGLES(nr)%BlockNo(is).ne.0) then
                nblock_save(ifound) = RECTANGLES(nr)%BlockNo(is)
              endif
            enddo
            do i=0,10
              xi = i/10.d0
              select case(ie)
              case(1)
                if (RECTANGLES(nr)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = xi; eta(2) = 0.d0
              case(2)
                if (RECTANGLES(nr)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = 1.d0; eta(2) = xi
              case(3)
                if (RECTANGLES(nr)%EdgeNo(ie).gt.0) xi = 1.d0 - xi
                eta(1) = xi; eta(2) = 1.d0
              case(4)
                if (RECTANGLES(nr)%EdgeNo(ie).gt.0) xi = 1.d0 - xi
                eta(1) = 0.d0; eta(2) = xi
              end select
              call recta(nr,eta, void,dxdeta)
              call cross_product(dxdeta(1:3,1),dxdeta(1:3,2), temp)
              call normalize(temp)
              if (i.eq.0) then
c
c  .............determine the orientation
                call scalar_product(temp,POINTS(np)%Rdata(4:6), s)
                factor=1.d0
                rsign = sign(factor,s)
              endif
              if (iprint.eq.4) then
                write(*,*) 'ifound,rsign,i,rn = ',
     .                      ifound,rsign,i,temp(1:3)
              endif
              select case(ifound)
              case(1)
                rnorm1(1:3,i) = rsign*temp(1:3)
              case(2)
                rnorm2(1:3,i) = rsign*temp(1:3)
              case default
                write(*,*) 'reconstruct: INCONSISTENCY '
                stop 1
              end select
            enddo
          enddo
        enddo
        do i=0,10
          do j=1,3
            if (abs(rnorm2(j,i)-rnorm1(j,i)).gt.GEOM_TOL) then
              write(*,*) 'ERROR IN COMP j = ',j,
     .                    abs(rnorm2(j,i)-rnorm1(j,i))
              write(*,8056) np,POINTS(np)%Rdata(4:6)
 8056         format('np,rn = ',i5,2x,3f8.3)
              write(*,8019) nr_save(1:2),nedg_save(1:2),nblock_save(1:2)
 8019         format('rectangles,edges,blocks = ',3(2i7,2x))
              write(*,8020) nc,i,rnorm1(1:3,i),rnorm2(1:3,i)
 8020         format('reconstruct: nc,i,rnorm1(1:3,i),rnorm2(1:3,i) = ',
     .                             i4,i2,2x,3f8.3,2x,3f8.3)
              write(*,8021) CURVES(nc)%EndPoNo(1:2)
 8021         format('             endpoints = ',2i5)
              call pause
              go to 50
            endif
          enddo
        enddo
 50     continue
c
c  ...end of loop through curves
      enddo
      if (iprint.ge.1) then
        write(*,*) 'reconstruct: HAVE CHECKED CONTINUITY'
      endif
c
c***********************************************************************
c
c  ...Step 5: Check the 3D jacobians (solid angles)
c
      do nh=1,NRHEXAS
        do k=0,10
          eta(3) = k/10.d0
          do j=0,10
            eta(2) = j/10.d0
            do i=0,10
              eta(1) = i/10.d0
              call hexa(nh,eta, void,dxdeta)
              do ivar=1,3
                call normalize(dxdeta(1:3,ivar))
              enddo
              call mixed_product(dxdeta(1:3,1),dxdeta(1:3,2),
     .                           dxdeta(1:3,3), angle)
              if ((angle.lt.angle_min).or.(angle.gt.angle_max))
     .          go to 100
            enddo
          enddo
        enddo
        cycle
 100    write(*,6036) nh,eta,angle
 6036   format('reconstruct: BAD SOLID ANGLE IN BLOCK nh = ',i7,
     .         ' AT POINT eta,rjac = ',3f8.3,2x,e12.5)
      enddo
c
c
      end





cccc----------------------------------------------------------------------
cccc
cccc   routine name       - collect_tangent
cccc
cccc----------------------------------------------------------------------
cccc
cccc   computer           - machine independent
cccc
cccc   latest revision    - Feb 07
cccc
cccc   purpose            - routine collects tanget vector for a curve
cccc                        at one of its endpoints and its derivative
cccc                        wrt to the curve parameter
cccc
cccc   arguments
cccc     in:
cccc              Nc      - a Hermite curve
cccc              Np      - a point
cccc     out:
cccc              Tangent - the tangent vector
cccc
cccc----------------------------------------------------------------------
cccc
ccc      subroutine collect_tangent(Nc,Np, Tangent,Dtangent)
cccc
ccc      use GMP
cccc
ccc#include "syscom.blk"
ccc#include "cpmax.blk"
cccc
ccc      dimension Tangent(1:3), Dtangent(1:3)
cccc
cccc  ...locate the point
ccc      call locate(Np,CURVES(Nc)%EndPoNo(1:2),2, i)
ccc      if (i.eq.0) then
ccc        write(*,7001) Nc,Np
ccc 7001   format('collect_tangent: Nc,Np = ',2i5)
ccc        stop 1
ccc      endif
cccc
ccc      select case(i)
ccc      case(1)
ccc        Tangent(1:3)  = CURVES(Nc)%Rdata(1:3)
ccc        Dtangent(1:3) = CURVES(Nc)%Rdata(7:9)
ccc      case(2)
ccc        Tangent(1:3) = -CURVES(Nc)%Rdata(4:6)
cccc
cccc  .....second order derivatives DO NOT change the sign
cccc       with orientation !!!
ccc        Dtangent(1:3) = CURVES(Nc)%Rdata(10:12)
ccc      end select
cccc
ccc      end
ccc
cccc----------------------------------------------------------------------
cccc
cccc   routine name       - check_edge
cccc
cccc----------------------------------------------------------------------
cccc
cccc   computer           - machine independent
cccc
cccc   latest revision    - May 07
cccc
cccc   purpose            - routine identifies a curve on a sharp edge
cccc
cccc   arguments
cccc     in:
cccc              Nc      - a curve
cccc     out:
cccc              Idec    = 2 if the curve is on a sharp edge
cccc                      = 0,1 if not
cccc
cccc----------------------------------------------------------------------
cccc
ccc      subroutine check_edge(Nc, Idec)
cccc
ccc      use GMP
cccc
ccc#include "syscom.blk"
ccc#include "cpmax.blk"
cccc
ccc      dimension x(3),nedg_points(100),nface_points(1000)
ccc      dimension eta(2),dxeta(3,2)
ccc      save ivisit,nedg_points,nr_nedg_points
ccc      data ivisit /0/
cccc
ccc      Idec=0
ccc      return
cccc
cccc  ...SHORT CUT FOR THE HEAD PROBLEM !!!!
ccc      if (ivisit.eq.0) then
ccc        ivisit=1
ccc        eta(1:2) = 0.5d0
cccc
cccc  .....identify all points on the lower section
ccc        nr_face_points=0
ccc        do nr=1,NRRECTA
ccc          if (RECTANGLES(nr)%BlockNo(2).eq.0) then
ccc            call recta_linear(nr,eta, x,dxdeta)
ccc            if (x(3).le.-29.d0) then
ccc              do iv=1,4
ccc                np = RECTANGLES(nr)%VertNo(iv)
ccc                call locate(np,nface_points,nr_face_points, ii)
ccc                if (ii.eq.0) then
ccc                  nr_face_points = nr_face_points + 1
ccc                  if (nr_face_points.gt.1000) then
ccc                    write(*,*) 'check_edge: nr_face_points = ',
ccc     .                          nr_face_points
ccc                    stop 1
ccc                  endif
ccc                  nface_points(nr_face_points) = np
ccc                endif
ccc              enddo
ccc            endif
ccc          endif
ccc        enddo
cccc
cccc  .....select only those points that are connected to
cccc       other surface points then on the lower surface
ccc        ipe=0
ccc        do ip=1,nr_face_points
ccc          np = nface_points(ip)
ccc          nrcurv = POINTS(np)%NrCurv
ccc          do j=1,nrcurv
ccc            nc1 = POINTS(np)%CurvNo(j)
ccc            do i=1,2
ccc              np1 = CURVES(nc1)%EndPoNo(i)
ccc              if (np1.eq.np) cycle
ccc              call locate(np1,nface_points,nr_face_points, ii)
ccc              if ((ii.eq.0).and.(POINTS(np1)%Type.eq.'CoorNrm')) then
cccc
cccc  .............add to the list of edge points
ccc                call locate(np,nedg_points,ipe, jj)
ccc                if (jj.eq.0) then
ccc                  ipe=ipe+1
ccc                  if (ipe.gt.100) then
ccc                    write(*,*) 'check_edge: ipe = ',ipe
ccc                    stop 1
ccc                  endif
ccc                  nedg_points(ipe) = np
ccc                endif
ccc              endif
ccc            enddo
ccc          enddo
ccc        enddo
ccc        nr_nedg_points = ipe
cccc
ccc        write(*,7001) nedg_points(1:nr_nedg_points)
ccc 7001   format('check_edge: EDGE POINTS = ',20i5,/,
ccc     .         '                          ',20i5,/,
ccc     .         '                          ',20i5,/)
ccc        call pause
ccc      endif
cccc
cccc  ...check if both endpoints are on the list of edge points
ccc      Idec=0
ccc      do j=1,2
ccc        np = CURVES(Nc)%EndPoNo(j)
ccc        call locate(np,nedg_points,nr_nedg_points, ii)
ccc        if (ii.ne.0) Idec=Idec+1
ccc      enddo
cccc
ccc      end
