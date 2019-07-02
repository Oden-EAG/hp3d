!-----------------------------------------------------------------------
!> Purpose : routine reads in the geometry data for the Geometry
!!           Modeling Package in the DEFAULT format
!!
!! @revision Nov 12
!-----------------------------------------------------------------------
!  DEFAULT file format:
!
!     NDIM (= 3)   MANDIM (= 3)
!
!     NRSURFS
!     SURFACES(1)%Type
!     SURFACES(1)%Rdata(:)   (as implied by Type)
!     ...
!
!     NRDOMAIN
!
!     NRPOINT
!     POINTS(1)%Type
!     POINTS(1)%Idata(:)   (as implied by Type)
!     POINTS(1)%Rdata(:)   (as implied by Type)
!     ...
!
!     NRCURVE
!     CURVES(1)%Type
!     CURVES(1)%EndPoNo(1:2)
!     CURVES(1)%Idata(:)   (as implied by Type)
!     CURVES(1)%Rdata(:)   (as implied by Type)
!     ...
!
!     NRTRIAN
!     TRIANGLES(1)%Type
!     TRIANGLES(1)%VertNo(1:3)
!     TRIANGLES(1)%Idata   (as implied by Type)
!     ...
!
!     NRRECTA
!     RECTANGLES(1)%Type
!     RECTANGLES(1)%VertNo(1:4)
!     RECTANGLES(1)%Idata   (as implied by Type)
!     ...
!
!     NRPRISM
!     PRISMS(1)%Type
!     PRISMS(1)%Domain     PRISMS(1)%VertNo(1:6)
!     ...
!
!     NRHEXAS
!     HEXAS(1)%Type
!     HEXAS(1)%Domain      HEXAS(1)%VertNo(1:8)
!     ...
!
!     NRTETRA
!     TETRAS(1)%Type
!     TETRAS(1)%Domain     TETRAS(1)%VertNo(1:4)
!     ...
!
!     NRPYRAM
!     PYRAMIDS(1)%Type
!     PYRAMIDS(1)%Domain   PYRAMIDS(1)%VertNo(1:5)
!     ...
!
!
!  Remarks:
!
!  1. Refer to "modules/GMP.F" for a detailed description of each
!     0,1,2,3D object.
!
!  2. Connectivities explicitly listed in the file are:
!
!       CURVES     --> ENDPOINTS ,
!       TRIANGLES  --> VERTICES  ,
!       RECTANGLES --> VERTICES  ,
!       3D BLOCKS  --> VERTICES  .
!
!     Remaining connectivities are automatically generated.
!
!  3. For each 3D block, the enumeration of the vertices dictates a
!     (local) system of coordinates for the block:
!
!        PRISM : Xi_1 = (v1,v2) ; Xi_2 = (v1,v3) ; Xi_3 = (v1,v4) ,
!        HEXA  : Xi_1 = (v1,v2) ; Xi_2 = (v1,v4) ; Xi_3 = (v1,v5) ,
!        TETRA : Xi_1 = (v1,v2) ; Xi_2 = (v1,v3) ; Xi_3 = (v1,v4) ,
!        PYRAM : Xi_1 = (v1,v2) ; Xi_2 = (v1,v4) ; Xi_3 = (v1,v5) .
!
!     If the implied system happens to be LEFT-oriented, it is
!     automatically modified into a RIGHT-oriented system by swapping
!     appropriate vertices.
!
!  4. Always punctuate formulas! (P. Gatto)
!-----------------------------------------------------------------------
subroutine input_DEFAULT(Fp)
!
      use GMP
      use error
      use environment , only : QUIET_MODE
!
      implicit none
      character(len=*), intent(in) :: Fp

      integer, parameter :: nin=16
      character(len=10) :: type
      integer :: ns,np,nc,j,k,nt,nr,npri,nh,ntet,npyr
      integer :: istat
      integer :: iprint
!-----------------------------------------------------------------------
!
      iprint=0
!
      open(unit=nin,file=Fp,form='formatted',access='sequential',status='unknown')
!
!  ...read in the dimension of the problem and the manifold
      read(nin,*) NDIM,MANDIM
!
!  ...allocate memory for GMP data structure
      call alloc_GMP
!
IF (.NOT. QUIET_MODE) write(*,*)'-- input_DEFAULT --'
!
!-----------------------------------------------------------------------
!  SURFACES                                                            |
!-----------------------------------------------------------------------
!
!  ...number of surfaces
      read(nin,*) NRSURFS
 IF (.NOT. QUIET_MODE)     write(*,1000) NRSURFS
 1000 format(' NRSURFS = ',i5,' ; reading surfaces...')
!
      if (MAXSU.lt.NRSURFS) then
        write(*,*)'input_DEFAULT: increase MAXSU!' ; stop
      endif
!
!  ...loop over surfaces
      do ns=1,NRSURFS
        read(nin,*) SURFACES(ns)%Type
        if (iprint.eq.1) then
          write(*,7003) ns,SURFACES(ns)%Type
 7003     format('   ns = ',i4,'; type = ',a10)
        endif
!
        select case(SURFACES(ns)%Type)
!
!  .....plane through Rdata(1:3), normal to Rdata(4:6)
        case('VecPt')
          allocate(SURFACES(ns)%Rdata(6))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4:6)
!
!  .....plane through Rdata(1:3), Rdata(4:6), Rdata(7:9)
        case('ThrPt')
          allocate(SURFACES(ns)%Rdata(9))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4:6)
          read(nin,*) SURFACES(ns)%Rdata(7:9)
!
!  .....plane through point Idata(1), normal to Rdata(1:3)
        case('PtNoVec')
          allocate(SURFACES(ns)%Idata(1))
          allocate(SURFACES(ns)%Rdata(3))
          read(nin,*) SURFACES(ns)%Idata(1)
          read(nin,*) SURFACES(ns)%Rdata(1:3)
!
!  .....plane through point Idata(1), and Rdata(1:3), Rdata(4:6)
        case('PtNo2Pt')
          allocate(SURFACES(ns)%Idata(1))
          allocate(SURFACES(ns)%Rdata(6))
          read(nin,*) SURFACES(ns)%Idata(1)
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4:6)
!
!  .....yz plane parametrized with cylindrical coordinates
        case('PPwCC')
          allocate(SURFACES(ns)%Rdata(3))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
!
!  .....sphere centered at Rdata(1:3), with radius Rdata(4)
        case('Sphere')
          allocate(SURFACES(ns)%Rdata(4))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4)
!
!    ...infinite cylinder with axis Rdata(4:6) through point Rdata(1:3),
!       and radius Rdata(7)
        case('Cylinder')
          allocate(SURFACES(ns)%Rdata(7))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4:6)
          read(nin,*) SURFACES(ns)%Rdata(7)
!
!  .....ellipsoid centered at Rdata(1:3), with x,y,z semiaxes of
!       lengths Rdata(4:6)
        case('Ellipsoid')
          allocate(SURFACES(ns)%Rdata(6))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4:6)
!
!  .....infinite cone with vertex at Rdata(1:3), axis Rdata(4:6),
!       aperture Rdata(7) [Rad]
        case('Cone')
          allocate(SURFACES(ns)%Rdata(7))
          read(nin,*) SURFACES(ns)%Rdata(1:3)
          read(nin,*) SURFACES(ns)%Rdata(4:6)
          read(nin,*) SURFACES(ns)%Rdata(7)
!
!  .....reconstructed surface
        case('RecSurf')
!
!  .....other cases
        case default
          write(*,7004) SURFACES(ns)%Type
 7004     format(' input_DEFAULT: unknown surface! Type = ',a10)
          stop
        endselect
      enddo
!
!---------------------------------------------------------------------
!  DOMAINS                                                           |
!---------------------------------------------------------------------
!
      read(nin,*) NRDOMAIN
!
!---------------------------------------------------------------------
!  POINTS                                                            |
!---------------------------------------------------------------------
!
!  ...read in number of points
      read(nin,*) NRPOINT
 IF (.NOT. QUIET_MODE)     write(*,1009) NRPOINT
 1009 format(' NRPOINT = ',i5,' ; reading points...')
!
      if (MAXNP.lt.NRPOINT) then
        write(*,*)'input_DEFAULT: increase MAXNP!' ; stop
      endif
!
!  ...loop over points
      do np=1,NRPOINT
        read(nin,*) POINTS(np)%Type
        if (iprint.eq.1) then
          write(*,1010) np,POINTS(np)%Type
 1010     format('   np = ',i4,'; type = ',a10)
        endif
!
        select case(POINTS(np)%Type)
!
!  .....regular point with coordinates Rdata(1:3)
        case('Regular')
          allocate(POINTS(np)%Rdata(3))
          read(nin,*) POINTS(np)%Rdata(1:3)
!
!  .....point defined implicitly as the intersecton of surfaces
!       Idata(1:3); Rdata(1:3) is initial guess for NR iterations
        case('Implicit')
          allocate(POINTS(np)%Idata(3))
          read(nin,*) POINTS(np)%Idata(1:3)
          allocate(POINTS(np)%Rdata(3))
          read(nin,*) POINTS(np)%Rdata(1:3)
!
!  .....point with coordinates Rdata(1:3) and normal Rdata(4:6)
        case('CoorNrm')
          allocate(POINTS(np)%Rdata(6))
          read(nin,*) POINTS(np)%Rdata(1:6)
!
!  .....point with coordiantes Rdata(1:3) on a sharp edge, with normals
!       Rdata(1:3) and Rdata(4:6)
        case('SharpPt')
          allocate(POINTS(np)%Rdata(9))
          read(nin,*) POINTS(np)%Rdata(1:3)
          read(nin,*) POINTS(np)%Rdata(4:6)
          read(nin,*) POINTS(np)%Rdata(7:9)
!
        case default
          write(*,1002) POINTS(np)%Type
 1002     format(' input_DEFAULT: unknown point type! Type = ',a10)
        endselect
!
      enddo
!
!---------------------------------------------------------------------
!  CURVES                                                            |
!---------------------------------------------------------------------

!  ...read in number of curves
      read(nin,*) NRCURVE
 IF (.NOT. QUIET_MODE)     write(*,1011) NRCURVE
 1011 format(' NRCURVE = ',i5,' ; reading curves...')
!
      if (MAXNC.lt.NRCURVE) then
        write(*,*) 'input_DEFAULT: increase MAXNC!' ; stop
        stop
      endif
!
!  ...loop over curves
      do nc=1,NRCURVE
        read(nin,*)  CURVES(nc)%Type
        read(nin,*) (CURVES(nc)%EndPoNo(j) , j=1,2)
!
        select case(CURVES(nc)%Type)
!
!  .....straight segment
        case('Seglin')
!
!  .....curve on algebraic surface Idata(1)
        case('1SurfsCur')
          allocate(CURVES(nc)%Idata(1), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) CURVES(nc)%Idata(1)
!
!  .....curve from intersection of 2 algebraic surfaces Idata(1:2)
        case('2SurfsCur')
          allocate(CURVES(nc)%Idata(2), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) (CURVES(nc)%Idata(k),k=1,2)
!
!  .....curve from intersection of 3 algebraic surfaces Idata(1:3)
        case('3SurfsCur')
          allocate(CURVES(nc)%Idata(3), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) (CURVES(nc)%Idata(k),k=1,3)
!
!  .....quarter of a circle centered at Rdata(1:3)
        case('QuaCir','QuaEl1','QuaEl2')
          allocate(CURVES(nc)%Rdata(NDIM), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) (CURVES(nc)%Rdata(k),k=1,NDIM)
!
!  .....quarter of a superellipse centered at Rdata(1:3) and powers 
!       R(4) and R(5) (regular ellipse corresponds to R(4)=R(5)=2)
        case('QuaSEl')
          allocate(CURVES(nc)%Rdata(NDIM+2), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) (CURVES(nc)%Rdata(k),k=1,NDIM+2)
!
!  .....implicit curve from the intersection of surfaces Idata(1:2),
!       and bounded by surfaces Idata(3:4)
        case('ImpCir')
          allocate(CURVES(nc)%Idata(4), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) (CURVES(nc)%Idata(k),k=1,4)
!
!  .....curve on the intersection of 1,2,3 reconstructed surfaces
!       (UNDEVELOPED! Straight segments used)
        case('1SurfrCur','2SurfrCur','3SurfrCur')
!
!  .....image of a straight line segment through a global system of
!       cylindrical coordinates
!        case('CylCur')
        case('CylCoord')
!
        endselect
!
      enddo
!
!---------------------------------------------------------------------
!  TRIANGLES                                                         |
!---------------------------------------------------------------------
!
      read(nin,*) NRTRIAN
 IF (.NOT. QUIET_MODE)     write(*,1012) NRTRIAN
 1012 format(' NRTRIAN = ',i5,' ; reading triangles...')
!
      if (MAXTR.lt.NRTRIAN) then
        write(*,*) 'MAXTR = ',MAXTR
        write(*,*) 'input_DEFAULT: increase MAXTR!' ; stop
      endif
!
!  ...loop over triangles
      do nt=1,NRTRIAN
        read(nin,*)  TRIANGLES(nt)%Type
        read(nin,*) (TRIANGLES(nt)%VertNo(j) , j=1,3)
!
        type=TRIANGLES(nt)%Type
        select case(type)
!
!  .....plane triangle, transfinite interpolation triangle
        case('PlaneTri','TransTri')
!
!  .....parametric transfinite interpolation triangle on surface Idata(1)
        case('PTITri')
          allocate(TRIANGLES(nt)%Idata(1), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) TRIANGLES(nt)%Idata(1)
!
!  .....triangle on G^1 reconstructed surface Idata(1)
        case('G1RecTri')
          TRIANGLES(nt)%Type='PlaneTri'
          allocate(TRIANGLES(nt)%Idata(1), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) TRIANGLES(nt)%Idata(1)
!
!  .....implicit triangle lying on surface Idata(1) and bounded by
!       surfaces Idata(2:4), listed counter clockwise
        case('ImpliTri')
          allocate(TRIANGLES(nt)%Idata(4), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) TRIANGLES(nt)%Idata(1:4)
!
!  .....image of a linear triangle through a global system of cylindrical
!       coordinates
        case('CylTri')
!
        case default
          write(*,1003)type
 1003     format(' input_DEFAULT: unknown triangle type! Type = ',a10)
          stop
        endselect
!
      enddo
!
!---------------------------------------------------------------------
!  RECTANGLES                                                        |
!---------------------------------------------------------------------
!
      read(nin,*) NRRECTA
 IF (.NOT. QUIET_MODE)     write(*,1013) NRRECTA
 1013 format(' NRRECTA = ',i5,' ; reading rectangles...')
!
      if (MAXRE.lt.NRRECTA) then
        write(*,*) 'input_DEFAULT: increase MAXRE!' ; stop
      endif
!
!  ...loop over rectangles
      do nr=1,NRRECTA
        read(nin,*)  RECTANGLES(nr)%Type
        read(nin,*) (RECTANGLES(nr)%VertNo(j) , j=1,4)
!
        select case(RECTANGLES(nr)%Type)
!
!  .....bilinear rectangle, transfinite interpolation rectangle
        case ('BilQua','TraQua')
!
!  .....parametric transfinite interpolation rectangle on surface Idata(1)
        case ('PTIRec')
          allocate(RECTANGLES(nr)%Idata(1), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) RECTANGLES(nr)%Idata(1)
!
!  .....implicit rectangle on surface Idata(1), bounded by surfaces Idata(2:4)
        case('ImpRec')
          allocate(RECTANGLES(nr)%Idata(5), stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
          endif
          read(nin,*) RECTANGLES(nr)%Idata(1:5)
!
!  .....image of a linear rectangle through a global system of cylindrical
!       coordinates
        case('CylRec')
!
        case default
          write(*,1004) RECTANGLES(nr)%Type
 1004     format(' input_DEFAULT: unknow rectangle type! Type = ',a10)
          stop
        endselect
      enddo
!
!---------------------------------------------------------------------
!  PRISMS
!---------------------------------------------------------------------
!
      read(nin,*) NRPRISM
 IF (.NOT. QUIET_MODE)     write(*,1014) NRPRISM
 1014 format(' NRPRISM = ',i5,' ; reading prisms...')
!
      if (MAXBT.lt.NRPRISM) then
        write(*,*) 'input_DEFAULT: increase MAXBT!' ; stop
      endif
!
!  ...loop over prisms
      do npri=1,NRPRISM
        read(nin,*) PRISMS(npri)%Type
        read(nin,*) PRISMS(npri)%domain,PRISMS(npri)%VertNo(1:6)
        call check_orientation(1,npri)
!
        select case(PRISMS(npri)%Type)
        case('Linear','TIprism')
        case default
          write(*,1005) PRISMS(npri)%Type
 1005     format(' input_DEFAULT: unknown prism type! Type = ',a10)
          stop
        endselect
      enddo
!
!---------------------------------------------------------------------
!  HEXAHEDRA
!---------------------------------------------------------------------
!
      read(nin,*) NRHEXAS
IF (.NOT. QUIET_MODE) write(*,1015) NRHEXAS
 1015 format(' NRHEXAS = ',i5,' ; reading hexas...')
!
      if (MAXHE.lt.NRHEXAS) then
        write(*,*) 'input_DEFAULT: increase MAXHE!' ; stop
      endif
!
      do nh=1,NRHEXAS
        read(nin,*) HEXAS(nh)%Type
        read(nin,*) HEXAS(nh)%domain,HEXAS(nh)%VertNo(1:8)
        call check_orientation(2,nh)
!
        select case(HEXAS(nh)%Type)
        case('Linear','TraHex','TrInHex','CylHex')
        case default
          write(*,1006) HEXAS(nh)%Type
 1006     format(' input_DEFAULT: unknown hexa type! Type = ',a10)
          stop
        endselect
      enddo
!
!---------------------------------------------------------------------
!  TETRAHEDRA                                                        |
!---------------------------------------------------------------------
!
      read(nin,*) NRTETRA
IF (.NOT. QUIET_MODE) write(*,1016) NRTETRA
 1016 format(' NRTETRA = ',i5,' ; reading tets...')
!
      if (MAXTE.lt.NRTETRA) then
        write(*,*) 'input_DEFAULT: increase MAXTE!'
        stop
      endif
      if (iprint.eq.1) then
        write(*,*) 'input_DEFAULT: READING TETRAS...'
        write(*,*) '           NRTETRA = ',NRTETRA
      endif
!
!
      do ntet=1,NRTETRA
        read(nin,*) TETRAS(ntet)%Type
        read(nin,*) TETRAS(ntet)%domain,TETRAS(ntet)%VertNo(1:4)
        call check_orientation(3,ntet)
!
        select case(TETRAS(ntet)%Type)
        case ('Linear','TraTet','CylTet')
        case default
          write(*,1007) type
 1007     format(' input_DEFAULT: unknown tet type! Type = ',a10)
          stop
        endselect
      enddo
!
!---------------------------------------------------------------------
!  PYRAMIDS                                                          |
!---------------------------------------------------------------------
!
      read(nin,*) NRPYRAM
IF (.NOT. QUIET_MODE) write(*,1017) NRPYRAM
 1017 format(' NRPYRAM = ',i5,' ; reading pyramids...')
!
      if (MAXPY.lt.NRPYRAM) then
        write(*,*) 'input_DEFAULT: increase MAXPY!'
        stop
      endif
      if (iprint.eq.1) then
        write(*,*) 'input_DEFAULT: READING PYRAMIDS...'
        write(*,*) '           NRPYRAM = ',NRPYRAM
      endif
!
!
      do npyr=1,NRPYRAM
        read(nin,*) PYRAMIDS(npyr)%Type
        read(nin,*) PYRAMIDS(npyr)%domain,PYRAMIDS(npyr)%VertNo(1:5)
        call check_orientation(4,npyr)
!
        select case(PYRAMIDS(npyr)%Type)
        case ('Linear','TIpyram')
        case default
          write(*,1008) PYRAMIDS(npyr)%Type
 1008     format(' input_DEFAULT: unknown pyramid type! Type = ',a10)
          stop
        endselect
      enddo
!
IF (.NOT. QUIET_MODE) write(*,*)''
!
      close(nin)
!
!///////////////////////////////////////////////////////////////////////
!  ...if needed, a geometry customization routine could be called HERE
!     call customize_geometry
!///////////////////////////////////////////////////////////////////////
!
!  ...complete full connectivities
      call connect
!
!
    endsubroutine input_DEFAULT
!
!
!
!-----------------------------------------------------------------------
!
!   routine name       - new_old
!
!-----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    - Aug 04
!
!   purpose            - routine reads in geometry data for GMP
!                        in the NEW FORMAT, allocates dynamically
!                        the necessary data structure arrays, and
!                        produces automatically a copy of the input
!                        file in the OLD FORMAT
!
!
!   arguments in       -
!
!   required  routines -
!
!-----------------------------------------------------------------------
!
      subroutine new_old
!
      use GMP
#include "syscom.blk"
#include "cinout.blk"
!
      character*10 type
!
!-----------------------------------------------------------------------
!
      iprint=1
!
!  ...read in the dimension of the problem and the manifold
      read(NWE,*)  NDIM,MANDIM
      write(NIN,*) NDIM,MANDIM,'             ...NDIM,MANDIM'
      write(NIN,*) ''
      write(NIN,*) ''
!
!  ...read in surfaces data
      read(NWE,*)  NRSURFS
      write(NIN,*) NRSURFS,   &
           '                         ...NUMBER OF SURFACES '
      write(NIN,*)''
!
!
!  ...allocate dynamically the memory
      allocate( SURFACES(NRSURFS), STAT=ic )
      if (ic.ne.0) then
        print*,'SURFACES ARE NOT ALLOCATED '
        stop
      endif
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING SURFACES...'
      endif
!
!  ...read in and store the data on surfaces...
      do isur=1,NRSURFS
!
!  .....transfer label of the points  to type of the points
        read(NWE,*) Type
        select case(Type)
!
!    ...plane with normal vector and one point
        case('VecPt')
          label= 1
!
!    ...plane with three points
        case('ThrPt')
          label = 2
!
!    ...sphere with radius and central point
        case('Sphere')
          label = 3
!
!    ...cylinder with vector of direction ,  point of begining and
!       radius
        case('Cylinder')
          label = 4
!
!  .....ellipsoid with a given center and axes aligned with
!       the global system of coordinates
        case('Ellipsoid')
          label = 7
!
!    ...other cases
        case default
          write(*,*)'input_geometry.f: WRONG TYPE FOR THE SURFACE!!'
        end select
        SURFACES(isur)%Type = Type
        write(NIN,*) label, &
        '                         ...LABEL OF SURFACE', isur
!
!
        select case(label)
!
!  .......plane which is normal to a given vector and passing through
!         a point, orientation specified by the normal
          case(1)
            allocate( SURFACES(isur)%Rdata(6), STAT = ic)
            if (ic.ne.0) then
              write(*,*) 'Rdate IS NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)SURFACES(isur)%Rdata(1:3), &
                 '          ...COORDINATES OF THE POINTS'
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3), &
                 '          ...COORDINATES OF THE NORMAL VECTOR'
!
!  .......plane passing through three given points, A,B,C
!         orientation specified by AB x AC
          case(2)
             allocate( SURFACES(isur)%Rdata(9), STAT = ic)
             if(ic.ne.0)then
                write(*,*) 'Rdate IS  NOT ALLOCATED', &
                     'SURFACE NUMBER = ',isur
                stop
            endif
            read(NWE,*)( SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3), &
                 '          ...COORDINATES OF THE FIRST POINT'
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3), &
                 '          ...COORDINATES OF THE SECOND POINT'
            read(NWE,*)(SURFACES(isur)%Rdata(6+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(6+j),j=1,3), &
                 '          ...COORDINATES OF THE THIRD POINT'
!
!  .......sphere
          case(3)
            allocate( SURFACES(isur)%Rdata(4), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'SURFACES COORDINATES ARE NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3),  &
                 '          ...COORDINATES OF THE CENTER'
            read(NWE,*)SURFACES(isur)%Rdata(4)
            write(NIN,100)SURFACES(isur)%Rdata(4),  &
                 '                            ...RADIUS OF THE SPHERE'
!
!  .......cylinder
          case(4)
            allocate( SURFACES(isur)%Rdata(7), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdate IS  NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3)
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3)
            read(NWE,*) SURFACES(isur)%Rdata(7)
            write(NIN,*) SURFACES(isur)%Rdata(7)
!
!  .......ellipsoid with semi-axes a,b,c along x,y,z and given central
!         point
          case(7)
            allocate( SURFACES(isur)%Rdata(6), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdata IS  NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
!
!  .........read coordinates of the central point
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
!
!  .........read semiaxes lengths
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
!
!  .........write out the data to the new format input file
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3)
!
!
!  .......other cases
          case default
            write(*,*)'new_old: WRONG SURFACE TYPE'
            stop
        end select
!
        write(NIN,*)''
      enddo
!
!----------------------------------------------------------------------
!
!                            P O I N T S
!
!  ...read in number of points
      read(NWE,*) NRPOINT
      write(NIN,*)''
      write(NIN,*)NRPOINT, &
      '                         ...NUMBER OF POINTS'
      write(NIN,*)''
!
!  ...allocate the memory dynamically
      allocate( POINTS(NRPOINT), STAT=ic )
      if (ic.ne.0) then
        write(*,*) 'new_old: POINTS ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING POINTS...'
      endif
!
!  ...read in and store the points data
      do ip=1,NRPOINT
        read(NWE,*)Type
        POINTS(ip)%Type  = Type
        POINTS(ip)%NrCurv = 0
!
        select case(Type)
!
!  .....read in and store the coordinates for regular points
        case('Regular')
          allocate( POINTS(ip)%Rdata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Rdata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Rdata(k),k=1,3)
!
        case('CoorNrm')
          allocate( POINTS(ip)%Rdata(6), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Rdata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Rdata(k),k=1,6)
!
!  .....read in and store the intersecting surfaces data and
!       starting point for implicit  points
        case('Implicit')
          allocate (POINTS(ip)%Idata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Idata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Idata(j),j=1,3)
          allocate( POINTS(ip)%Rdata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Rdata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Rdata(j),j=1,3)
!
!
        case default
           write(*,*)'new_old: WRONG POINT TYPE'
        end select
!
      enddo
!
!----------------------------------------------------------------------
!
!                        C U R V E S
!
!
!  ...read in number of curves
      read(NWE,*) NRCURVE
!
      allocate( CURVES(NRCURVE), STAT=ic )
      if( ic.ne.0 )then
        print*,'CURVES ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING CURVES...'
      endif
!
      do ic = 1, NRCURVE
        read(NWE,*) type
        CURVES(ic)%NrFig = 0
        CURVES(ic)%Type  = type
!
!  .....read in and store endpoints for the curve
        read(NWE,*) (CURVES(ic)%EndPoNo(j),j=1,2)
!
        select case (type)
        case('QuaCir','SegCir','QuaEl1','QuaEl2')
!
!  .......read in and store the coordinates of the center of a circle
          allocate(CURVES(ic)%Rdata(NDIM))
          read(NWE,*) (CURVES(ic)%Rdata(k),k=1,NDIM)
!
        case('QuaSEl')
!
!  .......read in and store the coordinates of the center of a 
!         superellipse and powers
          allocate(CURVES(ic)%Rdata(NDIM+2))
          read(NWE,*) (CURVES(ic)%Rdata(k),k=1,NDIM+2)
!
        case('ImpCir')
!
!  .......read in and store the intersecting surfaces
          allocate( CURVES(ic)%Idata(4), STAT = icc)
          if (icc.ne.0) then
            write(*,*) 'Idata IS NOT ALLOCATED', &
                       'THE CURVE NUMBER IS ', icur
            stop
          endif
          read(NWE,*) (CURVES(ic)%Idata(k),k=1,4)
        end select
!
      enddo
!
!---------------------------------------------------------------------
!
!                       T R I A N G L E S
!
!  ...read in number of triangles
      read(NWE,*) NRTRIAN
      if (NRTRIAN.ne.0) then
        write(*,*) 'new_old: UNFINISHED '
        stop
      endif
!
!
!  ...allocate dynamically the memory
      allocate( TRIANGLES(NRTRIAN), STAT=ic )
      if( ic.ne.0 )then
        print*,'TRIANGLES ARE NOT ALLOCATED'
        stop
      endif
!
!---------------------------------------------------------------------
!
!                        R E C T A N G L E S
!
!  ...read in and store the rectangles data
      read(NWE,*) NRRECTA
!
!  ...allocate dynamically the memory
      allocate( RECTANGLES(NRRECTA), STAT=ic )
      if (ic.ne.0) then
        print*,'RECTANGLES ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING RECTANGLES...'
      endif
!
      do irec=1,NRRECTA
        read(NWE,*) Type
        RECTANGLES(irec)%EdgeNo(1:4) = 0
        RECTANGLES(irec)%BlockNo(1:2) = 0
        RECTANGLES(irec)%Type  = Type
!
        read(NWE,*) RECTANGLES(irec)%VertNo(1:4)
!
!  .....read in and store surfaces that define the implicit rectangle
        if (type.eq.'ImpRec')then
          allocate(RECTANGLES(irec)%Idata(5), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Idata IS NOT ALLOCATED', &
                       'RECTANGLE NUMBER IS', irec
            stop
          endif
          read(NWE,*) (RECTANGLES(irec)%Idata(j),j=1,5)
        endif
!
      enddo
!
      read(NWE,*) NRPRISM
!
!
!  ...read in number of hexas
      read(NWE,*) NRHEXAS
!
!  ...allocate dynamically the memory
      allocate(HEXAS(NRHEXAS), STAT=ic )
      if (ic.ne.0) then
        print*,'new_old: HEXAS ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING HEXAS...'
      endif
!
      do ih=1,NRHEXAS
        read(NWE,*) HEXAS(ih)%Type
        read(NWE,*) HEXAS(ih)%VertNo(1:8)
!
!  ...end of loop through hexahedra
      enddo
!
!  ...determine missing connectivities and orientations
      call connect
!
!
!
!*************************************************************
!...write down the result into input file
!*************************************************************
!
      do i=1,NRPOINT
        select case(POINTS(i)%Type)
        case('Regular')
          write(NIN,*)1,  &
               '                         ...LABEL OF POINT',i
          write(NIN,*)POINTS(i)%NrCurv, &
               '                         ...NUMBER OF CURVES', &
               ' MEETING AT THE POINT '
          NrCurv = POINTS(i)%NrCurv
          select case(nrcurv)
            case(4)
               write(NIN,40)POINTS(i)%CurvNo(1:NrCurv), &
                    '             ...CURVES NUMBERS'
            case(5)
               write(NIN,50)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
            case(6)
               write(NIN,60)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
          end select

          write(NIN,300) (POINTS(i)%Rdata(k),k=1,3),'   ' &
               ,'       ...COORDINATES OF THE POINT'
          write(NIN,*)''
!
        case('CoorNrm')
          write(NOUT,*)30, &
               '                         ...LABEL OF POINT',i
          write(NIN,*)POINTS(i)%NrCurv, &
               '                   ...NUMBER OF CURVES', &
               ' MEETING AT THE POINT '
          NrCurv = POINTS(i)%NrCurv
          select case(nrcurv)
            case(4)
               write(NIN,40)POINTS(i)%CurvNo(1:NrCurv), &
                    '             ...CURVES NUMBERS'
            case(5)
               write(NIN,50)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
            case(6)
               write(NIN,60)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
          end select
          write(NIN,600) (POINTS(i)%Rdata(k),k=1,6),'   ' &
               ,'      ...COORDINATES AND NORMAL VECTORS'
          write(NIN,*)''
!
        case('Implicit')
          write(NIN,*) 50, &
               '                         ...LABEL OF POINT',i
          write(NIN,*) POINTS(i)%NrCurv, &
               '                         ...NUMBER OF CURVES', &
               ' MEETING AT THE POINT '
          NrCurv = POINTS(i)%NrCurv
          select case(nrcurv)
            case(4)
               write(NIN,40) POINTS(i)%CurvNo(1:NrCurv), &
                    '             ...CURVES NUMBERS'
            case(5)
               write(NIN,50) POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
            case(6)
               write(NIN,60) POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
          end select
            write(NIN,30)(POINTS(i)%Idata(j),j=1,3), &
                 '                   ...SURFACES THAT CONSTITUTE THE &
                 POINT'
          write(NIN,300)(POINTS(i)%Rdata(k),k=1,3), &
                 '          ...INITIAL APPROXIMATION POINT'
!
        case default
          write(*,*) 'new_old: WRONG POINT LABEL', POINTS(i)%Type
        end select
      enddo
!
!
      write(NIN,*)''
      write(NIN,*)NRCURVE  ,'                         ...NUMBER', &
           ' OF CURVES'
      write(NIN,*)''
      do i=1,NRCURVE
        if(CURVES(i)%Type.eq.'NormCoord')then
            write(NIN,*)30,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '              ...FIGURES NUMBERS'
            end select
            write(NIN,*)''
!
        elseif(CURVES(i)%Type.eq.'ImpCir')then
            write(NIN,*)50,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVES'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '              ...FIGURES NUMBERS'
            end select
            write(NIN,40)(CURVES(i)%Idata(k),k=1,4), &
                 '             ...SURFACES CONSTITUTING THE CURVE'
            write(NIN,*)''
!
        elseif(CURVES(i)%Type.eq.'Seglin')then
            write(NIN,*)1,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '             ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaCir')then
            write(NIN,*)-1,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

         elseif(CURVES(i)%Type.eq.'SegCir')then
            write(NIN,*)-2,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaEl1')then
            write(NIN,*)-3,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaEl2')then
            write(NIN,*)-4,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''
            
        elseif(CURVES(i)%Type.eq.'QuaSEl')then
            write(NIN,*)-5,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''
        else
            write(*,*)'INPUT_HEAD: WRONG CURVE LABEL', CURVES(i)%Type
        endif
      enddo
!
!
      write(NIN,*)''
      write(NIN,*) NRTRIAN,  &
           '                         ...NUMBER OF TRIANGLES'
!
!
      write(NIN,*)''
      write(NIN,*)''
      write(NIN,*)NRRECTA ,  &
           '                         ...NUMBER OF RECTANGLES'
      do i = 1,NRRECTA
         if(RECTANGLES(i)%Type.eq.'BilQua')then
            write(NIN,*)''
            write(NIN,*)1, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
                 '                         ...ADJACENT', &
                 'HEXAHEDRONS NUMBERS'
!
         else if(RECTANGLES(i)%Type.eq.'TraQua')then
            write(NIN,*)''
            write(NIN,*)3, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRONS NUMBERS'
!
          else if(RECTANGLES(i)%Type.eq.'RatioRec')then
            write(NIN,*)''
            write(NIN,*)30, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRON NUMBERS'
!
         else if(RECTANGLES(i)%Type.eq.'ImpRec')then
            write(NIN,*)''
            write(NIN,*)50, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRONS NUMBERS'
            write(NIN,50)RECTANGLES(i)%Idata(1:5), &
                 '       ...SURFACES THAT CONSISTUTE THE RECTANGLE'
!
         elseif(RECTANGLES(i)%Type.eq.'CylRec')then
            write(NIN,*)''
            write(NIN,*)-1, &
                 '                        ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRONS NUMBERS'
         else
            write(*,*)'INPUT_HEAD: WRONG RECTANGLE LABEL',  &
                 RECTANGLES(i)%Type
         endif
      enddo
!
!
      write(NIN,*)''
      write(NIN,*) NRPRISM  , &
           '                         ...NUMBER OF PRISMS'
!
!
      write(NIN,*)''
      write(NIN,*)''
      write(NIN,*) NRHEXAS  ,  &
           '                         ...NUMBER OF HEXAHEDRONS'
      do iii = 1,  NRHEXAS
         write(NIN,*)''
         write(NIN,*)1, &
              '                         ...TYPE OF HEXAHEDRON',III
         write(NIN,60) (HEXAS(iii)%FigNo(j), j = 1,6), &
              ' ...RECTANGULAR FACES NUMBERS'
       enddo
!
 100   format(1x,1f9.3,5a)
 300   format(1x,3f9.3,5a)
 600   format(1x,6f9.3,5a)
!
 20    format(1x,2i6,9a)
 30    format(1x,3i6,9a)
 40    format(1x,4i6,9a)
 50    format(1x,5i6,9a)
 60    format(1x,6i6,9a)
!
!
     end subroutine new_old
