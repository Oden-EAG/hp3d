c----------------------------------------------------------------------
c   routine name       - split_surface_OLD
c----------------------------------------------------------------------
c   latest revision    - Dec 08
c
c   purpose            - Routine enables generation of a "membrane"
c                        by splitting points lying on a predefined
c                        surface. Triangles lying on the surface are
c                        duplicated and thin prisms, tetrahedra
c                        and pyramids are generated.
c
c   arguments
c     in:
c            Nsplit    - plane to be split and the twin plane resulting
c                        from splitting
c            Nr_bound  - number of surfaces bounding the area to be split
c            Ns_bound  - list of bounding surfaces
c            Nr_confm  - number of surfaces to conform to
c            Ns_confm  - list of surfaces to conform to
c            Dh1,Dh2   - shift to the left and right for the plane to
c                        be split
c     out:
c                        changes recorded in GMP data structure
c
c   REMARK: Routine operates only with curve to vertex, figure to
c           vertex, and block to vertex connectivities; complete
c           connectivities are determined in routine 'connect'
c           AFTERWARDS
c---------------------------------------------------------------------
c
      subroutine split_surface_OLD(Nsplit,Nr_bound,Ns_bound,
     .                                Nr_confm,Ns_confm,Dh1,Dh2)
c
      use kinds
      use control
      use GMP
      use parallelism
#include "syscom.blk"
#include "cinout.blk"
#include "mpif.h"
c
c  ...bounding surfaces and surfaces to conform to
      dimension Ns_bound(10),Ns_confm(10)
c
c  ...list of surfaces to conform to
      dimension nosrf(10)
c
c  ...list of figures to be split
      integer, dimension(:,:), pointer :: figs_split
c  ...list of duplicated points
      integer, dimension(:), pointer   :: new_point
c  ...list of new curves
      integer, dimension(:,:), pointer :: new_curves
c
c  ...arguments for 'surf'
      dimension xp(3),dfdx(3)
c
c  ...unit normal to the split plane
      dimension rnorm(3),rnorm1(3),rnorm2(3)
c
c  ...work space
      dimension void2(4),void3(4)
c
c  ...module function cycles btw 1,...,mod
      imod(j,mod) = j-(j-1)/mod*mod
c
c printing flag
#define I_PRINT
c======================================================================
c    PG, Jul 09                                                       |
c  REMARK: use sloppy geometric tolerance since geometry has not      |
c  been update to double precision yet!                               |
c----------------------------------------------------------------------
      real(DP) :: GEOM_TOL_ORIG                                       !
c  ...save original geometric tolerance                               !
      GEOM_TOL_ORIG = GEOM_TOL                                        !
c  ...set sloppy geometric tolerance                                  !
      GEOM_TOL = 0.1D-5                                               !
#ifdef I_PRINT
      write(*,*)'split_surface: setting GEOM_TOL = ',GEOM_TOL         !
#endif
c======================================================================
c
      iprint = 0
c
c  ...save current geometry parameters
      nrcurve_orig = NRCURVE
      nrtrian_orig = NRTRIAN
      nrrecta_orig = NRRECTA
      nrprism_orig = NRPRISM
      nrtetra_orig = NRTETRA
      nrpyram_orig = NRPYRAM
c
c----------------------------------------------------------------------
c  Step 0: Classify points into the following groups                  !
c         -4   points in the interior of the split surface            !

c              that have to conform to at least a listed surface      !
c         -3   points in the interior of the split surface            !
c          0   points on the intersection of the split and            !
c              bounding surfaces                                      !
c         -1   points on the negative side of the surface             !
c         -2   points on the positive side of the surface             !
c----------------------------------------------------------------------
c  ...allocate new_point and initialize it to 0 (needed in Step 1)
      allocate(new_point(MAXNP))
      new_point(1:MAXNP) = 0
c  ...loop through points
      do np = 1, NRPOINT
        xp(1:3) = POINTS(np)%Rdata(1:3)
        call surf(Nsplit,xp, fval,dfdx)
c  .....classify point
        if (abs(fval).le.GEOM_TOL) then
          new_point(np)=-3
          if (if_bound(xp,Nr_bound,Ns_bound).le.0) new_point(np)=0
        elseif (fval.lt.-GEOM_TOL) then
          new_point(np)=-1
        elseif (fval.gt.GEOM_TOL) then
          new_point(np)=-2
        endif
c  ...end of loop through points
      enddo
c
c----------------------------------------------------------------------
c  Step 1: Identify figures to be split; 2 iterations are needed:     |
c    1st iteration determines nr_figs_to_split                        |
c    2nd iteration fill array figs_to_split                           |
c----------------------------------------------------------------------
c  ...loop through iterations
      do iter = 1, 2
c  .....set figure counter to 0
        if = 0
c
c  *********************  T R I A N G L E S  **********************   |
c
c  .....loop through triangles
        do nt = 1, NRTRIAN
c======================================================================
c  ion:  number of vertices on INTERIOR of split surface              |
c  iout: number of vertices on POS/NEG SIDE, or INTERIOR_&_CONFORMING |
c  REMARK: vertices on BOUNDING_&_SPLIT are not considered!           |
c======================================================================
          ion = 0;  iout = 0
c  .......loop through vertices
          do iv = 1, 3
            np = TRIANGLES(nt)%VertNo(iv)
c  .........if vertex on INTERIOR
            if (new_point(np) .eq. -3) then
              ion = ion + 1
c  .........else vertex on POS/NEG SIDE, or INTERIOR_&_CONFORMING
            elseif (new_point(np) .lt. 0) then
              iout = iout + 1
            endif
c  .......end of loop through vertices
          enddo
c======================================================================
c  IF (at least one vertex is in INTERIOR) .AND.                      |
c        (no vertex is on POS/NEG SIDE or INTERIOR_&_CONFORMING) THEN |
c----------------------------------------------------------------------
c    1. figure needs to be split;                                     |
c    2. newly created edges will NOT have to conform to any surface;  |
c    3. newly created faces will NOT have to conform to any surface;  |
c    4. newly created block will be a prism or pyramid depending on   |
c       whether one vertex is on BOUNDING_&_SPLIT.                    |
c======================================================================
          if ((ion .ge. 1) .and. (iout .eq. 0)) then
c  .........increment figure counter
            if = if + 1
c  .........if 2nd iteration
            if (iter .eq. 2) then
c  ...........store figure type (1 = trian,  2 = recta)
              figs_split(1,if) = 1
c  ...........store triangle number
              figs_split(2,if) = nt
c  ...........check if triangle is on the split surface
c======================================================================
c  REMARK: when splitting two intersecting surfaces it may happen     |
c  that the triangle does not result on the split surface; disregard  |
c  the issue for now.                                                 |
c======================================================================
              if (.not.associated(TRIANGLES(nt)%Idata)) then
                write(*,6001) nt
ccc                call print_GMP
ccc                call pause
              elseif (TRIANGLES(nt)%Idata(1) .ne. Nsplit) then
                write(*,6001) nt
 6001           format(' split_surface: TRIANGLE nt = ',i5,
     .                 ' NOT ON THE SPLIT SURFACE')
ccc                call print_GMP
ccc                call pause
              endif
            endif
          endif
c  .....end of loop through triangles
        enddo
c
c  *******************  R E C T A N G L E S  ***********************  |
c
c  .....loop through rectangles
        do nr = 1, NRRECTA
          ion = 0;  iout = 0
c  .......loop through vertices
          do iv = 1, 4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np) .eq. -3) then
              ion=ion+1
            elseif (new_point(np).lt.0) then
              iout=iout+1
            endif
c  .......end of loop through vertices
          enddo
c  .......if rectangle needs to be split
          if ((ion.eq.2).and.(iout.eq.0)) then
            if = if + 1
c  .........if 2nd interation store figure type and number
            if (iter .eq. 2) then
              figs_split(1,if) = 2
              figs_split(2,if) = nr
            endif
          endif
c  .....end of loop through rectangles
        enddo
c
c  .....save number of figures to split
        nr_figs_to_split = if
        if (nr_figs_to_split.eq.0) then
          write(*,*) 'split_surface: nr_figs_to_split = ',
     .                                         nr_figs_to_split
          stop
        endif
c  .....allocate figs_to_split if 1st interations
        if (iter .eq. 1) allocate(figs_split(2,nr_figs_to_split))
c  ...end of loop through interations
      enddo

c  ...printing statement
      if (iprint .eq. 1) then
        write(*,8001) nr_figs_to_split
 8001   format('split_surface: nr_figs_to_split = ',i4)
        write(*,8002) ((figs_split(1,i),figs_split(2,i)),
     .                  i=1,nr_figs_to_split)
 8002   format(10(i2,':',i8))
        call pause
      endif
c
c----------------------------------------------------------------------
c  STEP 1.2: do a few things...
c
c  ...add one subdomain
      NRDOMAIN = NRDOMAIN + 1
      write(*,8003) NRDOMAIN
 8003 format(' split_surface: adding subdomain, NRDOMAINS = ',i2)
ccc      call pause
c  ...check if split surface is a plane
      if (SURFACES(Nsplit)%Type .ne. 'VecPt') then
        write(*,*) 'split_surface: CANNOT SPLIT ns = ',ns
        stop
      endif
c  ...determine normal of split plane
      call norm(SURFACES(Nsplit)%Rdata(4:6), s)
      rnorm(1:3) = SURFACES(Nsplit)%Rdata(4:6)/s
c  ...add FOUR planes to GMP data structure
      do i = 1, 4
c  .....add one surface
        NRSURFS = NRSURFS + 1
c  .....check if MAXSU is exceeded
        if (NRSURFS .gt. MAXSU) then
          write(*,*) 'split_surface: increase MAXSU.'
          stop
        endif
        SURFACES(NRSURFS)%Type = 'VecPt'
        allocate(SURFACES(NRSURFS)%Rdata(6))
c  .....create Rdata for first TWO planes only
        select case(i)
c  .......plane on NEG side of split plane
          case(1)
            SURFACES(NRSURFS)%Rdata(1:3) = SURFACES(Nsplit)%Rdata(1:3)
     .                                   - rnorm(1:3)*Dh1
            SURFACES(NRSURFS)%Rdata(4:6) = rnorm(1:3)
c  .......plane on POS side of split plane
          case(2)
            SURFACES(NRSURFS)%Rdata(1:3) = SURFACES(Nsplit)%Rdata(1:3)
     .                                   + rnorm(1:3)*Dh2
            SURFACES(NRSURFS)%Rdata(4:6) = rnorm(1:3)
        end select
c======================================================================
c  REMARK: other TWO planes are used shortly after, don't panic!      |
c======================================================================
      enddo
c
c----------------------------------------------------------------------
c  Step 2: Loop through figures to be split and duplicate their
c          vertex points except for those located on the terminating
c          surface; also generate new curves connecting the twin points
c
c  ...loop through figures
      do if = 1, nr_figs_to_split
        select case(figs_split(1,if))
          case(1);  nt = figs_split(2,if);  nvrt = 3
          case(2);  nr = figs_split(2,if);  nvrt = 4
        end select
c  .....loop through figure vertex points
        do iv = 1, nvrt
          select case(figs_split(1,if))
            case(1);  np = TRIANGLES(nt)%VertNo(iv)
            case(2);  np = RECTANGLES(nr)%VertNo(iv)
          end select
c  .......select point type (BOUND_&_SPLIT or INTERIOR)
          select case(new_point(np))
c  .........the point is not to be duplicated
            case(0)
c======================================================================
c  REMARK: changing classification for new_point(np):                 |
c    -4  old classification                                           |
c    -3  old classification                                           |
c    -2  old classification                                           |
c    -1  old classification                                           |
c    np  point needs to be duplicated                                 |
c======================================================================
              new_point(np) = np
c  .........duplicate point and create new curve
            case(-3)
              xp(1:3) = POINTS(np)%Rdata(1:3)
c  ...........determine surfaces VERTEX must conform to
              call give_surf(xp,Nr_confm,Ns_confm, nrsrf,nosrf)
c  ...........select number of surfaces to conform to
              select case(nrsrf)
c  .............if no surfaces to conform to set up 2 "ad hoc" planes
                case(0)
c  ...............rotate normal to split plane by +90 deg
                  do i = 1, 3
                    i1 = imod(i + 1,3);  i2 = imod(i + 2,3)
c  .................choose direction of rotation based on rnorm(i)
                    if (rnorm(i) .ne. 0d0) then
                      rnorm1(i1) =  rnorm(i)
                      rnorm1(i)  = -rnorm(i1)
                      rnorm1(i2) = 0.d0
                      exit
                    endif
                  enddo
c  ...............normalize rnorm1
                  call normalize(rnorm1)
c  ...............set up 1st plane
                  SURFACES(NRSURFS - 1)%Rdata(1:3) = xp(1:3)
                  SURFACES(NRSURFS - 1)%Rdata(4:6) = rnorm1(1:3)
c  ...............compute rnorm2 orthogonal to both rnorm and rnorm1
                  call cross_product(rnorm,rnorm1, rnorm2)
c  ...............set up 2nd plane
                  SURFACES(NRSURFS)%Rdata(1:3) = xp(1:3)
                  SURFACES(NRSURFS)%Rdata(4:6) = rnorm2(1:3)
c  ...............store planes number for Newton method
                  nosrf(1) = NRSURFS - 1;  nosrf(2) = NRSURFS
c  .............one surface to conform to
                case(1)
c  ...............set up 1 "ad hoc" plane
                  call surf(nosrf(1),xp, fval,dfdx)
                  call normalize(dfdx)
                  call cross_product(rnorm,dfdx, rnorm1)
                  SURFACES(NRSURFS - 1)%Rdata(1:3) = xp(1:3)
                  SURFACES(NRSURFS - 1)%Rdata(4:6) = rnorm1(1:3)
c  ...............store plane number for Newton method
                  nosrf(2) = NRSURFS - 1
c  .............two surfaces to conform to, nothing to do
                case(2)
c  .............three or more surfaces
                case default
                  write(*,*) 'split_surface: ????'
                  stop
c  ...........end select number of surfaces to conform to
              end select
c  ...........set third surface for Newton method to plane on NEG side
              nosrf(3) = NRSURFS - 3
c  ...........call Newton method using original point as initial guess
              call mnewt(1,nosrf,void1,void2,xp,void3,
     .                                        POINTS(np)%Rdata(1:3))
c  ...........generate new point
              NRPOINT = NRPOINT + 1
              if (NRPOINT .gt. MAXNP) then
                write(*,*) 'split_surface: INCREASE MAXNP'
                stop
              endif
              POINTS(NRPOINT)%Type  = 'Regular'
              allocate(POINTS(NRPOINT)%Rdata(1:3))
c  ...........set third surface for Newton method to plane on POS side
              nosrf(3) = NRSURFS - 2
c  ...........call Newton method using original point as initial guess
              call mnewt(1,nosrf,void1,void2,xp,void3,
     .                                   POINTS(NRPOINT)%Rdata(1:3))
              if (iprint .eq. 2) then
                write(*,*)'split_surface: OLD xp and NEW POINTS',
     .                                       ' np,nrsrf = ',np,nrsrf
                write(*,*) xp(1:3)
                write(*,*) POINTS(np)%Rdata(1:3)
                write(*,*) POINTS(NRPOINT)%Rdata(1:3)
                call pause
              endif
c  ...........store new point to its old twin point connectivity
c======================================================================
c  REMARK: changing classification for new_point(np):                 |
c    -4           old classification                                  |
c    -3           old classification                                  |
c    -2           old classification                                  |
c    -1           old classification                                  |
c    ntwin_point  point has been duplicated                           |
c======================================================================
              new_point(np) = NRPOINT
c  ...........generate new curve
              NRCURVE = NRCURVE + 1
              if (NRCURVE .gt. MAXNC) then
                write(*,*)'split_surface: INCREASE MAXNC'
                stop
              endif
              CURVES(NRCURVE)%Type       = 'Seglin'
              CURVES(NRCURVE)%EndPoNo(1) = np
              CURVES(NRCURVE)%EndPoNo(2) = NRPOINT
c  .......end select point type
          end select
c  .....end loop through figure vertices
        enddo
c  ...end of loop through figures
      enddo
c
c  ...eliminate the last two temporary planes in GMP data structure
      do i = 1, 2
        deallocate(SURFACES(NRSURFS)%Rdata)
        NRSURFS = NRSURFS - 1
      enddo
c
c-----------------------------------------------------------------------
c  Step 3: Loop through figures to be split and generate remaining
c          new CURVES, TRIANGLES, and RECTANGLES.
c
      max_new_curves = 1000
      allocate(new_curves(2,max_new_curves))
      new_curves(1:2,1:max_new_curves) = 0
c
c  ...initiate new curve counter
      ic_new = 0
c  ...loop through figures to be split
      do if = 1, nr_figs_to_split
c  .....select figure type
        select case(figs_split(1,if))
          case(1);  nt = figs_split(2,if);  nedg = 3
          case(2);  nr = figs_split(2,if);  nedg = 4
        end select
c
c------------------------------------------------------------------------
c  STEP 3.1: loop through the figure edges, duplicate curves and generate
c            new rectangles and triangles RESULTING FROM THE DUPLICATION
c            OF EDGE CURVES.
c------------------------------------------------------------------------
c  .....loop through figure edges
        do ie = 1, nedg
c  .......local edge endpoints
          iv = ie;  iv1 = imod(iv + 1,nedg)
c  .......global edge endpoints
          select case(figs_split(1,if))
            case(1)
              np  = TRIANGLES(nt)%VertNo(iv)
              np1 = TRIANGLES(nt)%VertNo(iv1)
            case(2)
              np  = RECTANGLES(nr)%VertNo(iv)
              np1 = RECTANGLES(nr)%VertNo(iv1)
          end select
c  .......if both vertices are on the bounding surface don't duplicate curve
          if ((new_point(np).eq.np).and.(new_point(np1).eq.np1)) cycle  ! <---------------- (*)
c  .......determine if curve has already been visited (curves are visited twice)
          call locate_curve(np,np1,new_curves,ic_new, ifound)
c  .......if the curve has not been found, i.e. first visit, duplicate curve
          if (ifound .eq. 0) then
c  .........add a curve
            NRCURVE = NRCURVE + 1
c  .........check if MAXNC is exceeded
            if (NRCURVE .gt. MAXNC) then
              write(*,*)'split_surface: INCREASE MAXNC'
              stop
            endif
c  .........store curve info
            CURVES(NRCURVE)%Type = 'Seglin'                             ! <-- duplicating CURVE
c======================================================================
c  REMARK: classification for new_point(np):                          |
c    -4           initial classification                              |
c    -3           initial classification                              |
c    -2           initial classification                              |
c    -1           initial classification                              |
c    ntwin_point  point has been duplicated                           |
c======================================================================
c  .........endpoints of new curve are twin points of old curve
            CURVES(NRCURVE)%EndPoNo(1) = new_point(np)                  ! <-- duplicating CURVE
            CURVES(NRCURVE)%EndPoNo(2) = new_point(np1)                 ! <-- duplicating CURVE
c  .........increase new curves counter
            ic_new = ic_new + 1
            if (ic_new .gt. max_new_curves) then
              write(*,*)'split_surface: INCREASE max_new_curves'
              stop
            endif
c  .........update list of new curves endpoints
            new_curves(1,ic_new) = np
            new_curves(2,ic_new) = np1
c======================================================================
c  REMARK: recall that we are looping over curves that need to be     |
c    duplicated, see (*) above, hence at least one curve endpoint     |
c    is not on the bounding surface!                                  |
c======================================================================
c           if 1st endpoint lays on bounding surface (hence 2nd doesn't)
            if (new_point(np) .eq. np) then
c  ...........generate a new triangle
              NRTRIAN = NRTRIAN + 1
              if (NRTRIAN .gt. MAXTR) then
                write(*,*)'split_surface: increase MAXTR.'
                stop
              endif
              TRIANGLES(NRTRIAN)%Type      = 'PlaneTri'                 ! <-- duplicating TRIANGLE
              TRIANGLES(NRTRIAN)%VertNo(1) = np                         ! <-- duplicating TRIANGLE
              TRIANGLES(NRTRIAN)%VertNo(2) = np1                        ! <-- duplicating TRIANGLE
              TRIANGLES(NRTRIAN)%VertNo(3) = new_point(np1)             ! <-- duplicating TRIANGLE
c  .........ONLY 2nd endpoint lays on bounding surface
            elseif (new_point(np1) .eq. np1) then
c  ...........generate a new triangle
              NRTRIAN = NRTRIAN + 1
              if (NRTRIAN .gt. MAXTR) then
                write(*,*)'split_surface: increase MAXTR.'
                stop
              endif
              TRIANGLES(NRTRIAN)%Type      = 'PlaneTri'                 ! <-- duplicating TRIANGLE
              TRIANGLES(NRTRIAN)%VertNo(1) = np                         ! <-- duplicating TRIANGLE
              TRIANGLES(NRTRIAN)%VertNo(2) = np1                        ! <-- duplicating TRIANGLE
              TRIANGLES(NRTRIAN)%VertNo(3) = new_point(np)              ! <-- duplicating TRIANGLE
c  .........neither endpoint lays on bounding surface
            else
c  ...........generate new rectangle
              NRRECTA = NRRECTA + 1
              if (NRRECTA .gt. MAXRE) then
                write(*,*) 'split_surface: INCREASE MAXRE'
                stop
              endif
              RECTANGLES(NRRECTA)%Type      = 'BilQua'                  ! <-- duplicating RECTANGLE
              RECTANGLES(NRRECTA)%VertNo(1) = np                        ! <-- duplicating RECTANGLE
              RECTANGLES(NRRECTA)%VertNo(2) = np1                       ! <-- duplicating RECTANGLE
              RECTANGLES(NRRECTA)%VertNo(3) = new_point(np1)            ! <-- duplicating RECTANGLE
              RECTANGLES(NRRECTA)%VertNo(4) = new_point(np)             ! <-- duplicating RECTANGLE
c
c  ...........determine if the rectangle lies on a surface to conform to
              do is=1,Nr_confm
                ns = Ns_confm(is)
                idec=0
                do iv=1,4
                  np = RECTANGLES(NRRECTA)%VertNo(iv)
                  xp(1:3) = POINTS(np)%Rdata(1:3)
                  call surf(ns,xp, fval,dfdx)
                  if (abs(fval).lt.GEOM_TOL) idec=idec+1
                enddo
                if (idec.eq.4) then
                  RECTANGLES(NRRECTA)%Type = 'PTIRec'                   ! <-- duplicating RECTANGLE
                  allocate(RECTANGLES(NRRECTA)%Idata(1))
                  RECTANGLES(NRRECTA)%Idata(1)=ns                       ! <-- duplicating RECTANGLE
ccc                  write(*,8038) NRRECTA,ns
 8038             format('split_surface: ATTACHING RECTANGLE ',i5,
     .                   ' TO SURFACE', i3)
                  exit
                endif
              enddo
            endif
c
c  .......end if 1st visit to curve
          endif
c  .....end of loop through figure edges
        enddo
c
c------------------------------------------------------------------------
c  STEP 3.2: split the figure
c------------------------------------------------------------------------
c  .....select figure type
        select case(figs_split(1,if))
c  .....figure is a triangle
        case(1)
c  .......determine whether the old and new triangle need to conform to new planes
          idec = 0
c  .......loop through triangle vertices
          do iv = 1, 3
            np = TRIANGLES(nt)%VertNo(iv)
c  .........if np is different from its twin point idec++
            if (new_point(np) .ne. np) idec = idec + 1
          enddo
c  .......add triangle
          NRTRIAN = NRTRIAN + 1
c  .......check if MAXTR is exceeded
          if (NRTRIAN .gt. MAXTR) then
            write(*,*)'split_surface: increase MAXTR.'
            stop
          endif
c  .......loop through vertices of triangle to be duplicated
          do iv = 1, 3
            np = TRIANGLES(nt)%VertNo(iv)
c  .........store vertices of new triangle
            TRIANGLES(NRTRIAN)%VertNo(iv) = new_point(np)               ! <-- duplicating TRIANGLE
          enddo
          select case(idec)
c  .......up to 2 vertices have been duplicated
          case(0,1,2)
c  .........detach old triangle from the original surface if dh1>0
c           and not a plane triangle (presence of a plane triangle
c           may result from successive splittings...)
            if ((Dh1 .gt. 0.d0) .and.
     .                     (TRIANGLES(nt)%Type .ne. 'PlaneTri')) then
              deallocate(TRIANGLES(nt)%Idata)
              TRIANGLES(nt)%Type = 'PlaneTri'
            endif
            TRIANGLES(NRTRIAN)%Type = 'PlaneTri'                        ! <-- duplicating TRIANGLE
c  .......all 3 vertices have been duplicated
          case(3)
c  .........attach old triangle to the plane on NEG side
            if (Dh1 .gt. 0.d0) TRIANGLES(nt)%Idata(1) = NRSURFS - 1
            TRIANGLES(NRTRIAN)%Type = 'PTITri'
            allocate(TRIANGLES(NRTRIAN)%Idata(1))
c  .........attach new triangle to the plane on POS side
            TRIANGLES(NRTRIAN)%Idata(1) = NRSURFS                       ! <-- duplicating TRIANGLE
          end select
c  .....figure is a rectangle
        case(2)
c  .......determine whether the old and new rectangles need to conform to new planes
          idec = 0
          do iv = 1, 4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np) .ne. np) idec = idec + 1
          enddo
c  .......add rectangle
          NRRECTA = NRRECTA + 1
c  .......check if MAXRE is exceeded
          if (NRRECTA .gt. MAXRE) then
            write(*,*)'split_surface: increase MAXRE.'
            stop
          endif
          do iv = 1, 4
            np = RECTANGLES(nr)%VertNo(iv)
            RECTANGLES(NRRECTA)%VertNo(iv) = new_point(np)              ! <-- duplicating RECTANGLE
          enddo
          select case(idec)
          case(0,1,2,3)
c
c  .........detach old rectangle from the original surface if dh1>0
            if (Dh1.gt.0.d0) then
              deallocate(RECTANGLES(nr)%Idata)
              RECTANGLES(nr)%Type = 'BilQua'
            endif
            RECTANGLES(NRRECTA)%Type = 'BilQua'                         ! <-- duplicating RECTANGLE
          case(4)
            if (Dh1.gt.0.d0) RECTANGLES(nr)%Idata(1) = NRSURFS-1
            RECTANGLES(NRRECTA)%Type = 'PTIRec'                         ! <-- duplicating RECTANGLE
            allocate(RECTANGLES(NRRECTA)%Idata(1))
            RECTANGLES(NRRECTA)%Idata(1) = NRSURFS                      ! <-- duplicating RECTANGLE
          end select
          if (iprint.eq.1) then
            write(*,7024) nr,NRRECTA
 7024       format('split_surface: HAVE DUPLICATED RECTANGLE',i5,
     .             ' NEW TWIN RECTANGLE = ',i5)
ccc            call print_GMP
          endif
c  .....end select figure type
        end select
c  ...end of loop through figures to split
      enddo
c
c---------------------------------------------------------------------
c  Step 4: Update connectivities to points for the
c          original mesh entities
c
c  ...loop through curves
      do nc=1,nrcurve_orig
        ileft=0; iright=0; ion=0
        do iv=1,2
          np = CURVES(nc)%EndPoNo(iv)
          if (new_point(np).gt.0) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.eq.1).and.(iright.eq.1)) then
          do iv=1,2
            np = CURVES(nc)%EndPoNo(iv)
            if (new_point(np).gt.0)
     .         CURVES(nc)%EndPoNo(iv) = new_point(np)
          enddo
        endif
      enddo
c
c  ...loop through triangles
      do nt=1,nrtrian_orig
        ileft=0; iright=0; ion=0
        do iv=1,3
          np = TRIANGLES(nt)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,3
            np = TRIANGLES(nt)%VertNo(iv)
            if (new_point(np).gt.0)
     .          TRIANGLES(nt)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: TRIANGLE CROSSING SPLIT SURFACE'
          write(*,*) 'nt = ',nt
          do iv=1,3
            np = TRIANGLES(nt)%VertNo(iv)
            write(*,*) 'iv,np,new_point(np) = ',iv,np,new_point(np)
          enddo
          call print_GMP
          stop 1
        endif
      enddo
c
c  ...loop through rectangles
      do nr=1,nrrecta_orig
        ileft=0; iright=0; ion=0
        do iv=1,4
          np = RECTANGLES(nr)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np).gt.0)
     .          RECTANGLES(nr)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: RECTANGLE CROSSING SPLIT SURFACE'
          write(*,*) 'nr = ',nr
          do iv=1,3
            np = RECTANGLES(nr)%VertNo(iv)
            write(*,*) 'iv,np,new_point(np) = ',iv,np,new_point(np)
          enddo
          call print_GMP
          stop 1
        endif
      enddo
c
c  ...loop through prisms
      do npri=1,nrprism_orig
        if (iprint.eq.1) then
          write(*,*) 'split_surface: RECONNECTING PRISM ',npri
        endif
        ileft=0; iright=0; ion=0
        do iv=1,6
          np = PRISMS(npri)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,6
            np = PRISMS(npri)%VertNo(iv)
            if (new_point(np).gt.0)
     .          PRISMS(npri)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: PRISM CROSSING SPLIT SURFACE'
          stop 1
        endif
      enddo
c
c  ...loop through tets
      do ntet=1,nrtetra_orig
        ileft=0; iright=0; ion=0
        do iv=1,4
          np = TETRAS(ntet)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,4
            np = TETRAS(ntet)%VertNo(iv)
            if (new_point(np).gt.0)
     .          TETRAS(ntet)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: TET CROSSING SPLIT SURFACE'
          stop 1
        endif
      enddo
c
c  ...loop through pyramids
      do npyr=1,nrpyram_orig
        ileft=0; iright=0; ion=0
        do iv=1,5
          np = PYRAMIDS(npyr)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,5
            np = PYRAMIDS(npyr)%VertNo(iv)
            if (new_point(np).gt.0)
     .          PYRAMIDS(npyr)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: PYRAMID CROSSING SPLIT SURFACE'
          stop 1
        endif
      enddo
c
c---------------------------------------------------------------------
c  Step 5: Generate new prisms, tets and pyramids
c
      do if=1,nr_figs_to_split
        select case(figs_split(1,if))
c
c  .....split triangle
        case(1)
          nt=figs_split(2,if)
c
c  .......determine how many points have been duplicated
          ile=0
          do iv=1,3
            np = TRIANGLES(nt)%VertNo(iv)
            if (new_point(np).ne.np) ile=ile+1
          enddo
c
          select case(ile)
c
c  .......generate a new tet
          case(1)
            NRTETRA = NRTETRA+1
            if (NRTETRA.gt.MAXTE) then
              write(*,*) 'split_surface: INCREASE MAXTE'
              stop 1
            endif
            TETRAS(NRTETRA)%Type = 'Linear'
            TETRAS(NRTETRA)%Domain = NRDOMAIN
            TETRAS(NRTETRA)%VertNo(1:3) = TRIANGLES(nt)%VertNo(1:3)
            do iv=1,3
              np = TRIANGLES(nt)%VertNo(iv)
              if (new_point(np).ne.np)
     .          TETRAS(NRTETRA)%VertNo(4) = new_point(np)
            enddo
            call check_orientation(3,NRTETRA)
            if (iprint.eq.1) then
              write(*,*) 'split_surface: GENERATED NEW TET = ',NRTETRA,
     .                    TETRAS(NRTETRA)%VertNo(1:4)
            endif

c
c  .......generate a new pyramid
          case(2)
            NRPYRAM = NRPYRAM+1
            if (NRPYRAM.gt.MAXPY) then
              write(*,*) 'split_surface: INCREASE MAXPY'
              stop 1
            endif
            PYRAMIDS(NRPYRAM)%Type = 'Linear'
            PYRAMIDS(NRPYRAM)%Domain = NRDOMAIN
            do ie=1,3
              iv=ie; iv1=imod(iv+1,3)
              np = TRIANGLES(nt)%VertNo(iv)
              np1 = TRIANGLES(nt)%VertNo(iv1)
              if ((new_point(np).ne.np).and.(new_point(np1).ne.np1))exit
            enddo
            PYRAMIDS(NRPYRAM)%VertNo(1) = np
            PYRAMIDS(NRPYRAM)%VertNo(2) = np1
            PYRAMIDS(NRPYRAM)%VertNo(3) = new_point(np1)
            PYRAMIDS(NRPYRAM)%VertNo(4) = new_point(np)
            iv2=imod(iv+2,3)
            PYRAMIDS(NRPYRAM)%VertNo(5) = TRIANGLES(nt)%VertNo(iv2)
            call check_orientation(4,NRPYRAM)
c
c  .......generate a new prism
          case(3)
            NRPRISM = NRPRISM+1
            if (NRPRISM.gt.MAXBT) then
              write(*,*) 'split_surface: INCREASE MAXBT'
              stop 1
            endif
            PRISMS(NRPRISM)%Type = 'Linear'
            PRISMS(NRPRISM)%Domain = NRDOMAIN
            do iv=1,3
              np = TRIANGLES(nt)%VertNo(iv)
              PRISMS(NRPRISM)%VertNo(iv) = np
              PRISMS(NRPRISM)%VertNo(3+iv) = new_point(np)
            enddo
            call check_orientation(1,NRPRISM)
c
          case default
            write(*,*) 'split_surface: INCONSISTENCY 1, ile = ',ile
            stop 1
          end select
c
c  .....split rectangle
        case(2)
          nr=figs_split(2,if)
c
c  .......determine how many points have been duplicated
          ile=0
          do iv=1,4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np).ne.np) ile=ile+1
          enddo
c
          select case(ile)
c
c  .......generate a new prism
          case(2)
            NRPRISM = NRPRISM+1
            if (NRPRISM.gt.MAXBT) then
              write(*,*) 'split_surface: INCREASE MAXBT'
              stop 1
            endif
            PRISMS(NRPRISM)%Type = 'Linear'
            PRISMS(NRPRISM)%Domain = NRDOMAIN
            do ie=1,4
              iv=ie; iv1=imod(iv+1,4)
              np = RECTANGLES(nr)%VertNo(iv)
              np1 = RECTANGLES(nr)%VertNo(iv1)
              if ((new_point(np).ne.np).and.(new_point(np1).eq.np1))exit
            enddo
            PRISMS(NRPRISM)%VertNo(1) = np
            PRISMS(NRPRISM)%VertNo(2) = new_point(np)
            PRISMS(NRPRISM)%VertNo(3) = np1
            iv2 = imod(iv+2,4);  iv3 = imod(iv+3,4)
            np2 = RECTANGLES(nr)%VertNo(iv2)
            np3 = RECTANGLES(nr)%VertNo(iv3)
            PRISMS(NRPRISM)%VertNo(4) = np3
            PRISMS(NRPRISM)%VertNo(5) = new_point(np3)
            PRISMS(NRPRISM)%VertNo(6) = np2
            call check_orientation(1,NRPRISM)
            if (iprint.eq.1) then
              write(*,7020) NRPRISM
 7020         format('split_surface: GENERATED NEW PRISM ',i6)
              write(*,7021) PRISMS(NRPRISM)%VertNo(1:6)
 7021         format('WITH VERTEX POINTS =',6i6)
ccc              call pause
            endif
c
          case default
            write(*,*) 'split_surface: INCONSISTENCY 2, ile = ',ile
            stop 1
          end select
        end select
c
c  ...end of loop through figures to split
      enddo
c
      deallocate(figs_split,new_point,new_curves)
c
c======================================================================
c     PG, Jul 09
c  ...reset original geometric tolerance
      GEOM_TOL = GEOM_TOL_ORIG
c======================================================================
c
      end subroutine split_surface_OLD
