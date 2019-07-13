c----------------------------------------------------------------------
c   latest revision    - May 12
c
c   purpose            - module defines date structure arrays
c----------------------------------------------------------------------

      module data_structure3D
c
      use physics
      use parameters
      use element_data
c
c  ...parameters
      integer, parameter :: NHIST = 20
c
c  ...dirichlet data set up
      integer, save :: NR_DIRICHLET_LIST = 0
      integer, save, dimension(20) :: DIRICHLET_LIST
c  ...dirichlet homogenous data set up
      integer, save :: NR_DIRICHLET_HOMOGENEOUS_LIST = 0
      integer, save, dimension(20) :: DIRICHLET_HOMOGENEOUS_LIST
c
c  ...number of initial mesh elements, active elements, nodes
      integer, save :: NRELIS,NRELES,NRNODS
c
c  ...total number of active H1,H(curl),H(div),L2 dofs
      integer, save :: NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
c
c  ...maximum number of nodes, pointer to the first free entry in NODES
      integer, save :: MAXNODS,NPNODS
c
c----------------------------------------------------------------------
c  INITIAL MESH ELEMENT                                               |
c----------------------------------------------------------------------
      type element
c
c  .....type - 'bric','tetr','pris','pyra'
        character(len=4) :: type
c
c  .....number of physics attributes supported BY the element
        integer          :: nrphysics
c
c  .....list   of physics attributes supported BY the element
        character(len=5), dimension(:), pointer :: physics
c
c  .....boundary conditions nickname (1 BC flag per face) for each
c       attribute supported BY the element.
c       Reserved values (see subroutine meshgen/set_index):
c         H1 :
c           0 - no BC
c           1 - Dirichlet BC on all components
c           3 - Dirichlet BC on 2nd and 3rd component
c           4 - Dirichlet BC on 1st and 3rd component
c           5 - Dirichlet BC on 1st and 2nd component
c           6 - Dirichlet BC on 1st component
c           7 - Dirichlet BC on 2nd component
c           8 - Dirichlet BC on 3rd component
c         H(div) :
c           0 - no BC
c           1 - Dirichlet BC on all components
c           3 - Dirichlet BC on 1st component
c           4 - Dirichlet BC on 2nd component
c           5 - Dirichlet BC on 3rd component
c           6 - Dirichlet BC on 2nd and 3rd component
c           7 - Dirichlet BC on 1st and 3rd component
c           8 - Dirichlet BC on 1st and 2nd component
c         H(curl), L2 :
c           0 - no BC
c           1 - Dirichlet BC on all components
        integer, dimension(:), pointer :: bcond
c
c  .....element nodal connectivities: vertices,edges,faces, middle node
        integer, dimension(:), pointer :: nodes
c
c  .....corresponding orientations for edges
        integer          :: edge_orient
c
c  .....corresponding orientations for faces
        integer          :: face_orient
c
c  .....neighbors across faces
        integer, dimension(:), pointer :: neig
c
c  .....nickname of GMP block (1,2,3,4 - prism, hexa, tet, pyramid)
c       coinciding with the element
        integer          :: GMPblock
      endtype element
c
c----------------------------------------------------------------------
c  NODE                                                               |
c----------------------------------------------------------------------
      type node
c
c  .....type - 'vert','medg','mdlt','mdlq','mdlb','mdln','mdlp','mdld'
        character(len=4) :: type
c
c  .....case number indicating what physical attributes are supported
        integer          :: case
c
c  .....nickname storing info about supported variables:
c         0 - component does not exist
c         1 - H1 component with Dirichlet BC flag
c         2 - free H1 component
c         3 - H(curl) component with Dirichlet BC flag
c         4 - free H(curl) component
c         5 - H(div) component with Dirichlet BC flag
c         6 - free H(div) component
c         7 - L2 component with Dirichlet BC flag
c         8 - free L2 component
        integer*8        :: index
        integer          :: iback
c
c  .....order of approximation
        integer          :: order
c
c  .....boundary condition flag
        integer          :: bcond
c
c  .....refinement flag
        integer          :: ref_kind, ref_filter
c
c  .....father node
        integer          :: father
c
c  .....node sons
        integer, dimension(:), pointer :: sons
c
c  .....interface flag with GMP
        integer          :: geom_interf

c  .....visitation flag
        integer          :: visit
c
c  .....activation flag
        integer          :: act
c
c  .....locker number
        integer          :: lock
c
c  .....subdomain number (distributed mesh)
        integer          :: subd
c
c  .....geometry dof
        real*8, dimension(:,:), pointer :: coord
c
c  .....error
c       0   - scalar error
c       1-3 - gradient of any vector component-wise error
c       4   - additional info
        real*8, dimension(0:4,0:4) :: error
c
c  .....H1 solution dof
#if C_MODE
        complex*16, dimension(:,:), pointer :: zdofH
#else
        real*8    , dimension(:,:), pointer :: zdofH
#endif
c
c  .....H(curl) solution dof
#if C_MODE
        complex*16, dimension(:,:), pointer :: zdofE
#else
        real*8    , dimension(:,:), pointer :: zdofE
#endif
c
c  .....H(div) solution dof
#if C_MODE
        complex*16, dimension(:,:), pointer :: zdofV
#else
        real*8    , dimension(:,:), pointer :: zdofV
#endif
c
c  .....L^2 solution dof
#if C_MODE
        complex*16, dimension(:,:), pointer :: zdofQ
#else
        real*8    , dimension(:,:), pointer :: zdofQ
#endif
      endtype node
c
c-----------------------------------------------------------------------
c
c  ...data structure arrays
      type(element), allocatable, save  :: ELEMS(:)
      type(node)   , allocatable, save  :: NODES(:)
c
c-----------------------------------------------------------------------
c
      contains
c
c-----------------------------------------------------------------------
c
      subroutine open_history_file(fp)
      character(len=*) :: fp
      open(unit=NHIST,file=fp,
     .     form='formatted',access='sequential',status='unknown')
      end subroutine
c
c-----------------------------------------------------------------------
c
      subroutine close_history_file
      write(NHIST,*) '0 0'
      close(NHIST)
      end subroutine
c
c-----------------------------------------------------------------------
c
c  ...determine number of dof for a higher order node
      subroutine find_ndof(Nod, NdofH,NdofE,NdofV,NdofQ)
c
      integer Nod,NdofH,NdofE,NdofV,NdofQ
c
      call ndof_nod(NODES(Nod)%type,NODES(Nod)%order,
     .              NdofH,NdofE,NdofV,NdofQ)
c
      end subroutine find_ndof
c
c-----------------------------------------------------------------------
c
c  ...find number of H1,H(curl),H(div),L2 variables supported by a node
      subroutine find_nvar(Nod, NvarH,NvarE,NvarV,NvarQ)
c
      integer :: Nod,NvarH,NvarE,NvarV,NvarQ
      integer, dimension(NRINDEX) :: idx
      integer :: k
c
      NvarH=0 ; NvarE=0 ; NvarV=0 ; NvarQ=0
c
      call get_index(Nod, idx)
      do k=1,NRINDEX
        select case (idx(k))
        case (2); nvarH = nvarH + 1
        case (4); nvarE = nvarE + 1
        case (6); nvarV = nvarV + 1
        case (8); nvarQ = nvarQ + 1
        end select
      end do
c
      end subroutine find_nvar
c
c-----------------------------------------------------------------------
c
c  ...determine number of sons for a higher order node
      subroutine find_nsons(Nod, Nrsons)
c
      integer Nod,Nrsons
c
      if (associated(NODES(Nod)%sons)) then
        Nrsons = ubound(NODES(Nod)%sons,1)
      else
        Nrsons = 0
      endif
c
      end subroutine find_nsons
c
c-----------------------------------------------------------------------
c
c  ...get index for a node
      subroutine get_index(Nod, Indexd)
c
      integer Indexd(NRINDEX)
c
      call decodLong(NODES(Nod)%index,10,NRINDEX, Indexd)
c
      end subroutine get_index
c
c-----------------------------------------------------------------------
c
c  ...allocate memory for data structure
      subroutine allocds
c
      if (allocated(ELEMS)) then
        write(*,*) 'allocds: WARNING !! DATA STRUCTURE ARRAYS',
     .             ' HAVE NOT BEEN DEALLOCATED'
        call deallocds
      endif
c
      allocate(ELEMS(NRELIS))
      do nel=1,NRELIS
        ELEMS(nel)%type = 'none'
        ELEMS(nel)%nrphysics = 0
        nullify (ELEMS(nel)%physics)
        nullify (ELEMS(nel)%bcond)
        nullify (ELEMS(nel)%nodes)
        ELEMS(nel)%edge_orient = 0
        ELEMS(nel)%face_orient = 0
        nullify (ELEMS(nel)%neig)
        ELEMS(nel)%GMPblock = 0
      enddo
c
      allocate(NODES(MAXNODS))
      do nod=1,MAXNODS
        NODES(nod)%type = 'none'
        NODES(nod)%case = 0
        NODES(nod)%index = 0
        NODES(nod)%order = 0
        NODES(nod)%bcond = nod+1
        NODES(nod)%ref_kind = 0
        NODES(nod)%father = 0
        nullify(NODES(nod)%sons)
        NODES(nod)%geom_interf = 0
        nullify (NODES(nod)%coord)
        NODES(nod)%error = 0.d0
        nullify (NODES(nod)%zdofH)
        nullify (NODES(nod)%zdofE)
        nullify (NODES(nod)%zdofV)
        nullify (NODES(nod)%zdofQ)
      enddo
      NODES(MAXNODS)%bcond = 0
      NPNODS=1
c
c
      end subroutine allocds
c
c-----------------------------------------------------------------------
c
c  ...deallocate data structure arrays
      subroutine deallocds
c
      do nel=1,NRELIS
        if (associated(ELEMS(nel)%physics))
     .    deallocate(ELEMS(nel)%physics)
        if (associated(ELEMS(nel)%bcond))
     .    deallocate(ELEMS(nel)%bcond)
        if (associated(ELEMS(nel)%nodes))
     .    deallocate(ELEMS(nel)%nodes)
        if (associated(ELEMS(nel)%neig))
     .    deallocate(ELEMS(nel)%neig)
      enddo
      deallocate(ELEMS)
c
      do nod=1,MAXNODS
        if (Associated(NODES(nod)%sons))
     .    deallocate(NODES(nod)%sons)
        if (Associated(NODES(nod)%coord))
     .    deallocate(NODES(nod)%coord)
        if (Associated(NODES(nod)%zdofH))
     .    deallocate(NODES(nod)%zdofH)
        if (Associated(NODES(nod)%zdofE))
     .    deallocate(NODES(nod)%zdofE)
        if (Associated(NODES(nod)%zdofV))
     .    deallocate(NODES(nod)%zdofV)
        if (Associated(NODES(nod)%zdofQ))
     .    deallocate(NODES(nod)%zdofQ)
      enddo
      deallocate(NODES)
c
c
      end subroutine deallocds
c
c-----------------------------------------------------------------------
c
c  ...dump out hp3d data structure
      subroutine dumpout_hp3d(Dump_file)
c
      character(len=15) :: Dump_file
      ndump=31
cc    kyungjoo
cc      open(unit=ndump,file=Dump_file,
cc     .     buffered='yes',blocksize = 65536,
cc     .     form='formatted',access='sequential',status='unknown')
      open(unit=ndump,file=Dump_file,
     .     form='formatted',access='sequential',status='unknown')

c
      write(ndump,*) NRELIS,NRELES,NRNODS
      write(ndump,*) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
      write(ndump,*) MAXNODS,NPNODS
c
      do nel=1,NRELIS
        write(ndump,*) ELEMS(nel)%type
        write(ndump,*) ELEMS(nel)%nrphysics
        if (associated(ELEMS(nel)%physics)) then
          nn = ubound(ELEMS(nel)%physics,1)
          write(ndump,*) nn
          write(ndump,1010) ELEMS(nel)%physics
 1010     format(1x,20(a5,2x))
        else
          write(ndump,*) 0
        endif
        if (associated(ELEMS(nel)%bcond)) then
          nn = ubound(ELEMS(nel)%bcond,1)
          write(ndump,*) nn
          write(ndump,*) ELEMS(nel)%bcond
        else
          write(ndump,*) 0
        endif
        if (associated(ELEMS(nel)%nodes)) then
          nn = ubound(ELEMS(nel)%nodes,1)
          write(ndump,*) nn
          write(ndump,*) ELEMS(nel)%nodes
        else
          write(ndump,*) 0
        endif
        write(ndump,*) ELEMS(nel)%edge_orient
        write(ndump,*) ELEMS(nel)%face_orient
        if (associated(ELEMS(nel)%neig)) then
          nn = ubound(ELEMS(nel)%neig,1)
          write(ndump,*) nn
          write(ndump,*) ELEMS(nel)%neig
        else
          write(ndump,*) 0
        endif
        write(ndump,*) ELEMS(nel)%GMPblock
      enddo
c
      do nod=1,NRNODS
        write(ndump,*) NODES(nod)%type
        write(ndump,*) NODES(nod)%case
        write(ndump,*) NODES(nod)%index
        write(ndump,*) NODES(nod)%order
        write(ndump,*) NODES(nod)%bcond
        write(ndump,*) NODES(nod)%ref_kind
        write(ndump,*) NODES(nod)%ref_filter
        write(ndump,*) NODES(nod)%father
        if (associated(NODES(nod)%sons)) then
          nn = ubound(NODES(nod)%sons,1)
          write(ndump,*) nn
          write(ndump,*) NODES(nod)%sons
        else
          write(ndump,*) 0
        endif
        write(ndump,*) NODES(nod)%geom_interf
        write(ndump,*) NODES(nod)%visit
        write(ndump,*) NODES(nod)%act
        if (associated(NODES(nod)%coord)) then
          nn1 = ubound(NODES(nod)%coord,1)
          nn2 = ubound(NODES(nod)%coord,2)
          write(ndump,*) nn1, nn2
          write(ndump,*) NODES(nod)%coord
        else
          write(ndump,*) 0 , 0
        endif
        write(ndump,*) NODES(nod)%error
        if (associated(NODES(nod)%zdofH)) then
          nn1 = ubound(NODES(nod)%zdofH,1)
          nn2 = ubound(NODES(nod)%zdofH,2)
          write(ndump,*) nn1, nn2
          write(ndump,*) NODES(nod)%zdofH
        else
          write(ndump,*) 0 , 0
        endif
        if (associated(NODES(nod)%zdofE)) then
          nn1 = ubound(NODES(nod)%zdofE,1)
          nn2 = ubound(NODES(nod)%zdofE,2)
          write(ndump,*) nn1, nn2
          write(ndump,*) NODES(nod)%zdofE
        else
          write(ndump,*) 0 , 0
        endif
        if (associated(NODES(nod)%zdofV)) then
          nn1 = ubound(NODES(nod)%zdofV,1)
          nn2 = ubound(NODES(nod)%zdofV,2)
          write(ndump,*) nn1, nn2
          write(ndump,*) NODES(nod)%zdofV
        else
          write(ndump,*) 0 , 0
        endif
        if (associated(NODES(nod)%zdofQ)) then
          nn1 = ubound(NODES(nod)%zdofQ,1)
          nn2 = ubound(NODES(nod)%zdofQ,2)
          write(ndump,*) nn1, nn2
          write(ndump,*) NODES(nod)%zdofQ
        else
          write(ndump,*) 0 , 0
        endif
      enddo
c
      close(ndump)
c
      end subroutine dumpout_hp3d
c
c
c-----------------------------------------------------------------------
c
c  ...dump in hp3d data structure in LEGACY format, i.e
c     NODES%visit , NODES%act are not read from 'dumpc3Dhp'
c
      subroutine dumpin_hp3d_LEGACY(Dump_file)
c
      character(len=15) :: Dump_file
      character(len=20) :: type
      integer           :: npnods_loc
      integer,parameter :: iprint = 0
c
c  ...deallocate DS if necessary
      if (allocated(ELEMS).or.allocated(NODES)) call deallocds
c
      ndump=31
      open(unit=ndump,file=Dump_file,
     .     form='formatted',access='sequential',status='unknown')
c
      read(ndump,*) NRELIS,NRELES,NRNODS
      read(ndump,*) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
      read(ndump,*) MAXNODS,NPNODS
c
      if (iprint.eq.1) then
        write(*,*)'------------------------------------------'
        write(*,*)'dumpin_hp3d_LEGACY:'
        write(*,8000)NRELIS
 8000   format('     NRELIS  = ',i12)
        write(*,8001)NRELES
 8001   format('     NRELES  = ',i12)
        write(*,8002)NRNODS
 8002   format('     NRNODS  = ',i12)
        write(*,8003)NRDOFSH
 8003   format('     NRDOFSH = ',i12)
        write(*,8004)NRDOFSE
 8004   format('     NRDOFSE = ',i12)
        write(*,8005)NRDOFSV
 8005   format('     NRDOFSV = ',i12)
        write(*,8006)NRDOFSQ
 8006   format('     NRDOFSQ = ',i12)
        write(*,8007)MAXNODS
 8007   format('     MAXNODS = ',i12)
        write(*,8008)NPNODS
 8008   format('     NPNODS  = ',i12)
        write(*,*)'------------------------------------------'
        call pause
      endif
c
c  ...save NPNODS
      npnods_loc = NPNODS
c  ...allocate data structure
      call allocds
c  ...reset NPNODS
      NPNODS = npnods_loc
c
      do nel=1,NRELIS
        read(ndump,*) ELEMS(nel)%type
        read(ndump,*) ELEMS(nel)%nrphysics
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%physics(nn))
          read(ndump,*) (ELEMS(nel)%physics(i), i=1,nn)
        else
          nullify(ELEMS(nel)%physics)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%bcond(nn))
          read(ndump,*) (ELEMS(nel)%bcond(i), i=1,nn)
        else
          nullify(ELEMS(nel)%bcond)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%nodes(nn))
          read(ndump,*) (ELEMS(nel)%nodes(i), i=1,nn)
        else
          nullify(ELEMS(nel)%nodes)
        endif
        read(ndump,*) ELEMS(nel)%edge_orient
        read(ndump,*) ELEMS(nel)%face_orient
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%neig(nn))
          read(ndump,*) (ELEMS(nel)%neig(i), i=1,nn)
        else
          nullify(ELEMS(nel)%neig)
        endif
        read(ndump,*) ELEMS(nel)%GMPblock
      enddo
c
c  ...N O D E S........................................................
      do nod=1,NRNODS
        read(ndump,*) NODES(nod)%type
        read(ndump,*) NODES(nod)%case
        read(ndump,*) NODES(nod)%index
        read(ndump,*) NODES(nod)%order
        read(ndump,*) NODES(nod)%bcond
        read(ndump,*) NODES(nod)%ref_kind
        read(ndump,*) NODES(nod)%father
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(NODES(nod)%sons(nn))
          read(ndump,*) NODES(nod)%sons
        else
          nullify(NODES(nod)%sons)
        endif
        read(ndump,*) NODES(nod)%geom_interf
c  .....set undefined flags to 0
        NODES(nod)%visit=0; NODES(nod)%act=0
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%coord(nn1,nn2))
          read(ndump,*) NODES(nod)%coord
        else
          nullify(NODES(nod)%coord)
        endif
c
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofH(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofH
        else
          nullify(NODES(nod)%zdofH)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofE(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofE
        else
          nullify(NODES(nod)%zdofE)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofV(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofV
        else
          nullify(NODES(nod)%zdofV)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofQ(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofQ
        else
          nullify(NODES(nod)%zdofQ)
        endif
      enddo
c
      close(ndump)
c
      end subroutine dumpin_hp3d_LEGACY
c
c
c-----------------------------------------------------------------------
c
c  ...dump in hp3d data structure
      subroutine dumpin_hp3d(Dump_file)
c
      character(len=15) :: Dump_file
      character(len=20) :: type
      integer           :: npnods_loc
c
      if (allocated(ELEMS).or.allocated(NODES)) call deallocds
c
      ndump=31
      open(unit=ndump,file=Dump_file,
     .     form='formatted',access='sequential',status='unknown')
c
      read(ndump,*) NRELIS,NRELES,NRNODS
      read(ndump,*) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
      read(ndump,*) MAXNODS,NPNODS
      npnods_loc = NPNODS
      call allocds
      NPNODS = npnods_loc
c
      do nel=1,NRELIS
        read(ndump,*) ELEMS(nel)%type
        read(ndump,*) ELEMS(nel)%nrphysics
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%physics(nn))
          read(ndump,*) (ELEMS(nel)%physics(i), i=1,nn)
        else
          nullify(ELEMS(nel)%physics)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%bcond(nn))
          read(ndump,*) (ELEMS(nel)%bcond(i), i=1,nn)
        else
          nullify(ELEMS(nel)%bcond)
        endif
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%nodes(nn))
          read(ndump,*) (ELEMS(nel)%nodes(i), i=1,nn)
        else
          nullify(ELEMS(nel)%nodes)
        endif
        read(ndump,*) ELEMS(nel)%edge_orient
        read(ndump,*) ELEMS(nel)%face_orient
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(ELEMS(nel)%neig(nn))
          read(ndump,*) (ELEMS(nel)%neig(i), i=1,nn)
        else
          nullify(ELEMS(nel)%neig)
        endif
        read(ndump,*) ELEMS(nel)%GMPblock
      enddo
c
      do nod=1,NRNODS
        read(ndump,*) NODES(nod)%type
        read(ndump,*) NODES(nod)%case
        read(ndump,*) NODES(nod)%index
        read(ndump,*) NODES(nod)%order
        read(ndump,*) NODES(nod)%bcond
        read(ndump,*) NODES(nod)%ref_kind
        read(ndump,*) NODES(nod)%ref_filter
        read(ndump,*) NODES(nod)%father
        read(ndump,*) nn
        if (nn.gt.0) then
          allocate(NODES(nod)%sons(nn))
          read(ndump,*) NODES(nod)%sons
        else
          nullify(NODES(nod)%sons)
        endif
        read(ndump,*) NODES(nod)%geom_interf
        read(ndump,*) NODES(nod)%visit
        read(ndump,*) NODES(nod)%act
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%coord(nn1,nn2))
          read(ndump,*) NODES(nod)%coord
        else
          nullify(NODES(nod)%coord)
        endif
        read(ndump,*) NODES(nod)%error
c
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofH(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofH
        else
          nullify(NODES(nod)%zdofH)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofE(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofE
        else
          nullify(NODES(nod)%zdofE)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofV(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofV
        else
          nullify(NODES(nod)%zdofV)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%zdofQ(nn1,nn2))
          read(ndump,*) NODES(nod)%zdofQ
        else
          nullify(NODES(nod)%zdofQ)
        endif
      enddo
c
      close(ndump)
c
      end subroutine dumpin_hp3d
c
c-----------------------------------------------------------------------
      subroutine add_dirichlet_to_list(Iboundary)
      integer loc
      loc = 0
      call locate(Iboundary, DIRICHLET_LIST, NR_DIRICHLET_LIST, loc)
      if (loc.eq.0) then
        NR_DIRICHLET_LIST = NR_DIRICHLET_LIST + 1
        DIRICHLET_LIST(NR_DIRICHLET_LIST) = Iboundary
      end if
      end subroutine

c-----------------------------------------------------------------------
      subroutine add_dirichlet_homogeneous_to_list(Iboundary)
      integer loc1,loc2
      loc1 = 0
      loc2 = 0
      call locate(Iboundary, DIRICHLET_LIST, NR_DIRICHLET_LIST, loc1)
      call locate(Iboundary, DIRICHLET_HOMOGENEOUS_LIST,
     .            NR_DIRICHLET_HOMOGENEOUS_LIST, loc2)
      if (loc2.eq.0) then
        NR_DIRICHLET_HOMOGENEOUS_LIST = NR_DIRICHLET_LIST + 1
        DIRICHLET_HOMOGENEOUS_LIST(NR_DIRICHLET_HOMOGENEOUS_LIST)
     .    = Iboundary
c       any homogenous B.C flag is also a Dirichlet flag
        if (loc1.eq.0) then
          NR_DIRICHLET_LIST = NR_DIRICHLET_LIST + 1
          DIRICHLET_LIST(NR_DIRICHLET_LIST) = Iboundary
        end if ! loc1
      end if !loc2
      end subroutine
c

c-----------------------------------------------------------------------
c  ...reset visitation flags for all nodes
      subroutine reset_visit
c
      do i=1,NRNODS
        NODES(i)%visit = 0
      enddo
c
      end subroutine reset_visit
c
c
c  ...get visitation flag for a node
      subroutine get_visit(Nod, Vis)
c
      integer, intent(in)  :: Nod
      integer, intent(out) :: Vis
c
      Vis = NODES(Nod)%visit
c
      end subroutine

c  ...set visitation flag of a node
      subroutine set_visit(Nod)
c
      integer, intent(in) :: Nod
c
      NODES(Nod)%visit = 1
c
      end subroutine
c
c-----------------------------------------------------------------------
c
c  ...get current subdomain of a node
      subroutine get_subd(Nod, Subd)
c
      integer, intent(in)  :: Nod
      integer, intent(out) :: Subd
c
      Subd = NODES(Nod)%subd
c
      end subroutine
c
c  ...set new subdomain of a node
      subroutine set_subd(Nod,Subd)
c
      integer, intent(in) :: Nod,Subd
c
      NODES(Nod)%subd = Subd
c
      end subroutine
c
c-----------------------------------------------------------------------
      function Is_dirichlet(Nod)
      integer Nod
      integer ibc(NR_PHYSA), loc
      logical Is_dirichlet

      call decod(NODES(Nod)%bcond,10,NR_PHYSA, ibc)
      Is_dirichlet = .false.
      do iphys=1,NR_PHYSA
        if (ibc(iphys).eq.1) then
          Is_dirichlet = .true.
        else
          call locate(ibc(iphys),
     .      DIRICHLET_LIST, NR_DIRICHLET_LIST, loc)
          if (loc.ne.0) then
            Is_dirichlet = .true.
          endif
        endif
      enddo
c
      end function
c----------------------------------------------------------------------
      function Is_dirichlet_homogeneous(Nod)
      integer Nod
      integer ibc(NR_PHYSA), loc
      logical Is_dirichlet_homogeneous

      call decod(NODES(Nod)%bcond,10,NR_PHYSA, ibc)
      Is_dirichlet_homogeneous = .false.
      do iphys=1,NR_PHYSA
        call locate(ibc(iphys),DIRICHLET_HOMOGENEOUS_LIST,
     .              NR_DIRICHLET_HOMOGENEOUS_LIST,loc)
        if (loc.ne.0) then
          Is_dirichlet_homogeneous = .true.
        endif
      enddo
c
      end function
c----------------------------------------------------------------------
      subroutine check_dirichlet_homogeneous(istat)
      integer :: istat

c     local variables
      integer :: ibc_number,ii

c     we check if any flag number in list of homogeneous Dirichlet b.c.
c     is also in the list
      istat=0
      do ii=1,SIZE(DIRICHLET_HOMOGENEOUS_LIST)
        ibc_number=DIRICHLET_HOMOGENEOUS_LIST(ii)
          call locate(ibc_number,
     .      DIRICHLET_LIST, NR_DIRICHLET_LIST, loc)
          if (loc.eq.0) then
            istat = 1
            write(*,*) 'ERROR check_dirichlet_homogeneous: ',
     .       'ibc=',ibc_number,'is not in DIRICHLET_LIST array'
            return
          endif
      end do

      end subroutine
c ------------------------------------------------------------------
      function Is_right_handed(Mdle)
      integer :: Mdle
      logical :: Is_right_handed
      integer :: i, nod
      real*8  :: v(3,4), a(3,3), val
      select case(ELEMS(Mdle)%type)
      case('bric','pris','pyra')
        Is_right_handed = .true.
      case('tetr')
c
        do i=1,4
          nod = ELEMS(Mdle)%nodes(i)
          v(1:3,i) = NODES(nod)%coord(1:3,1)
        enddo
c
        do i=1,3
          a(1:3,i) = v(1:3,i+1) - v(1:3,1)
        enddo
c
        call mixed_product(a(1:3,1), a(1:3,2), a(1:3,3), val)
        Is_right_handed = (val > 0.d0)
      end select
      end function
c-----------------------------------------------------------------------
      function Is_active(Nod)
      integer Nod
      logical Is_active
      select case (NODES(Nod)%act)
      case (0)
        Is_active = .FALSE.
      case default
        Is_active = .TRUE.
      end select
      end function
c-----------------------------------------------------------------------
      function Is_leaf(Nod)
      integer Nod
      logical Is_leaf
      select case (NODES(Nod)%ref_kind)
      case (0)
         Is_leaf = .TRUE.
      case default
         Is_leaf = .FALSE.
      end select
      end function
c-----------------------------------------------------------------------
      function Is_root(Nod)
      integer Nod
      logical Is_root
      if (NODES(Nod)%father.lt.0) then
        Is_root = .TRUE.
      else
        Is_root = .FALSE.
      endif
      end function
c-----------------------------------------------------------------------
      function Is_middle(Nod)
      integer Nod
      logical Is_middle
      select case(NODES(Nod)%type)
      case('mdlb','mdlp','mdln','mdld')
        Is_middle = .TRUE.
      case default
        Is_middle = .FALSE.
      end select
      end function
c
      end module data_structure3D
