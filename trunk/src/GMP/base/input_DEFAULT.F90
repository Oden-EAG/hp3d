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
!     isurf_flag     (surface domains flags
!                     =  0    no surface domains are read when inputing
!                             data for TRIANGLES and RECTANGLES
!                     =  1    surface domains are expected when inputing
!                             data for TRIANGLES and RECTANGLES  )
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
      integer :: ns,np,nc,j,k,nt,nr,npri,nh,ntet,npyr,isurf_flag
      integer :: istat
      integer :: iprint
!-----------------------------------------------------------------------
!
      iprint=0
!
      open(unit=nin,file=Fp, &
         form='formatted',access='sequential',status='old',action='read')
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
 1000 format(' NRSURFS = ',i7,' ; reading surfaces...')
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
!  POINTS  and surface domains flag                                  |
!---------------------------------------------------------------------
!
!  ...read in number of points
      read(nin,*) NRPOINT
!
!  ...most likely, we have read the surface domains flag
      if ((NRPOINT.eq.0).or.(NRPOINT.eq.1)) then
        isurf_flag = NRPOINT
!
!  .....proceed with reading the actual number of points
        read(nin,*) NRPOINT
!
!  ...most likely, the surface domains flag is missing, set it to 0
      else
        write(*,*) 'input_DEFAULT: surface domains flag missing, ',&
                   ' setting it to 0'
        call pause
        isurf_flag = 0
      endif 
!
      IF (.NOT. QUIET_MODE)     write(*,1009) NRPOINT
 1009 format(' NRPOINT = ',i7,' ; reading points...')
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
 1011 format(' NRCURVE = ',i7,' ; reading curves...')
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
 1012 format(' NRTRIAN = ',i7,' ; reading triangles...')
!
      if (MAXTR.lt.NRTRIAN) then
        write(*,*) 'MAXTR = ',MAXTR
        write(*,*) 'input_DEFAULT: increase MAXTR!' ; stop
      endif
!
!  ...loop over triangles
      do nt=1,NRTRIAN
        read(nin,*)  TRIANGLES(nt)%Type
        select case(isurf_flag)
        case(0); read(nin,*) (TRIANGLES(nt)%VertNo(j) , j=1,3)
        case(1); read(nin,*) TRIANGLES(nt)%Domain, (TRIANGLES(nt)%VertNo(j) , j=1,3)
        end select
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
!  .....triangle on G^1 reconstructed surface Idata(1) (LEGACY)
!        case('G1RecTri')
!          TRIANGLES(nt)%Type='PlaneTri'
!          allocate(TRIANGLES(nt)%Idata(1), stat=istat)
!          if (istat.ne.SUCCESS) then
!            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
!          endif
!          read(nin,*) TRIANGLES(nt)%Idata(1)
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
 1013 format(' NRRECTA = ',i7,' ; reading rectangles...')
!
      if (MAXRE.lt.NRRECTA) then
        write(*,*) 'input_DEFAULT: increase MAXRE!' ; stop
      endif
!
!  ...loop over rectangles
      do nr=1,NRRECTA
        read(nin,*)  RECTANGLES(nr)%Type
        select case(isurf_flag)
        case(0); read(nin,*) (RECTANGLES(nr)%VertNo(j) , j=1,4)
        case(1); read(nin,*) RECTANGLES(nr)%Domain, (RECTANGLES(nr)%VertNo(j) , j=1,4)
        end select
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
 1014 format(' NRPRISM = ',i7,' ; reading prisms...')
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
 1015 format(' NRHEXAS = ',i7,' ; reading hexas...')
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
 1016 format(' NRTETRA = ',i7,' ; reading tets...')
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
 1017 format(' NRPYRAM = ',i7,' ; reading pyramids...')
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
    end subroutine input_DEFAULT
