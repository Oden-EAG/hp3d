!----------------------------------------------------------------------
!   latest revision    - Aug 2019
!
!   purpose            - module defines date structure arrays
!----------------------------------------------------------------------
!
module data_structure3D
!
      use physics
      use parameters
      use element_data
      use mpi_param, only: RANK
!
!  ...parameters
      integer, parameter :: NHIST = 20
!
!  ...dirichlet data set up
      integer, save :: NR_DIRICHLET_LIST = 0
      integer, save, dimension(20) :: DIRICHLET_LIST
!  ...dirichlet homogenous data set up
      integer, save :: NR_DIRICHLET_HOMOGENEOUS_LIST = 0
      integer, save, dimension(20) :: DIRICHLET_HOMOGENEOUS_LIST
!
!  ...number of initial mesh elements, active elements (global,local), nodes
      integer, save :: NRELIS,NRELES,NRELES_SUBD,NRNODS
!
!  ...total number of active H1,H(curl),H(div),L2 dofs
      integer, save :: NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
!
!  ...maximum number of nodes, pointer to the first free entry in NODES
      integer, save :: MAXNODS=0,NPNODS=0
!
!----------------------------------------------------------------------
!  INITIAL MESH ELEMENT                                               |
!----------------------------------------------------------------------
      type element
!
!  .....type - 'bric','tetr','pris','pyra'
        character(len=4) :: type
!
!  .....number of physics attributes supported BY the element
        integer          :: nrphysics
!
!  .....list   of physics attributes supported BY the element
        character(len=5), dimension(:), pointer :: physics
!
!  .....boundary conditions nickname (1 BC flag per face) for each
!       attribute supported BY the element.
!       Reserved values (see subroutine meshgen/set_index):
!         H1 :
!           0 - no BC
!           1 - Dirichlet BC on all components
!           3 - Dirichlet BC on 2nd and 3rd component
!           4 - Dirichlet BC on 1st and 3rd component
!           5 - Dirichlet BC on 1st and 2nd component
!           6 - Dirichlet BC on 1st component
!           7 - Dirichlet BC on 2nd component
!           8 - Dirichlet BC on 3rd component
!         H(div) :
!           0 - no BC
!           1 - Dirichlet BC on all components
!           3 - Dirichlet BC on 1st component
!           4 - Dirichlet BC on 2nd component
!           5 - Dirichlet BC on 3rd component
!           6 - Dirichlet BC on 2nd and 3rd component
!           7 - Dirichlet BC on 1st and 3rd component
!           8 - Dirichlet BC on 1st and 2nd component
!         H(curl), L2 :
!           0 - no BC
!           1 - Dirichlet BC on all components
        integer, dimension(:), pointer :: bcond
!
!  .....element nodal connectivities: vertices,edges,faces, middle node
        integer, dimension(:), pointer :: nodes
!
!  .....corresponding orientations for edges
        integer          :: edge_orient
!
!  .....corresponding orientations for faces
        integer          :: face_orient
!
!  .....neighbors across faces
        integer, dimension(:), pointer :: neig
!
!  .....nickname of GMP block (1,2,3,4 - prism, hexa, tet, pyramid)
!       coinciding with the element
        integer          :: GMPblock
      endtype element
!
!----------------------------------------------------------------------
!  NODE                                                               |
!----------------------------------------------------------------------
      type node
!
!  .....type - 'vert','medg','mdlt','mdlq','mdlb','mdln','mdlp','mdld'
        character(4)     :: type
!
!  .....case number indicating what physical attributes are supported
        integer          :: case
!
!  .....nickname storing info about supported variables:
!         0 - component does not exist
!         1 - H1 component with Dirichlet BC flag
!         2 - free H1 component
!         3 - H(curl) component with Dirichlet BC flag
!         4 - free H(curl) component
!         5 - H(div) component with Dirichlet BC flag
!         6 - free H(div) component
!         7 - L2 component with Dirichlet BC flag
!         8 - free L2 component
        integer(8)       :: index
!
!  .....order of approximation
        integer          :: order
!
!  .....boundary condition flag
        integer          :: bcond
!
!  .....father node
        integer          :: father
!
!  .....node sons
        integer          :: first_son
        integer          :: nr_sons
        !integer(2)      :: nr_sons
!
!  .....refinement flag
        integer          :: ref_kind
        !integer(2)      :: ref_kind
!
!  .....interface flag with GMP
        integer          :: geom_interf
!
!  .....visitation flag
        integer          :: visit
!
!  .....activation flag
        integer          :: act
        !logical          :: act
!
!  .....subdomain number (distributed mesh)
        integer          :: subd
!
!  .....geometry dof
        real(8), dimension(:,:), pointer :: coord
!
#if DEBUG_MODE
!
        integer          :: iback
!
!  .....locker number
        integer          :: lock
!
!  .....error
!       0   - scalar error
!       1-3 - gradient of any vector component-wise error
!       4   - additional info
        real(8), dimension(0:4,0:4) :: error
#endif
!
!  .....H1 solution dof
#if C_MODE
        complex(8), dimension(:,:), pointer :: zdofH
#else
        real(8)   , dimension(:,:), pointer :: zdofH
#endif
!
!  .....H(curl) solution dof
#if C_MODE
        complex(8), dimension(:,:), pointer :: zdofE
#else
        real(8)   , dimension(:,:), pointer :: zdofE
#endif
!
!  .....H(div) solution dof
#if C_MODE
        complex(8), dimension(:,:), pointer :: zdofV
#else
        real(8)   , dimension(:,:), pointer :: zdofV
#endif
!
!  .....L2 solution dof
#if C_MODE
        complex(8), dimension(:,:), pointer :: zdofQ
#else
        real(8)   , dimension(:,:), pointer :: zdofQ
#endif
      endtype node
!
!-----------------------------------------------------------------------
!
!  ...data structure arrays
      type(element), allocatable, save  :: ELEMS(:)
      type(node)   , allocatable, save  :: NODES(:)
!
!  ...global and local (subdomain) order of elements
      integer      , allocatable, save  :: ELEM_ORDER(:)
      integer      , allocatable, save  :: ELEM_SUBD(:)
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!
      subroutine update_ELEM_ORDER()
         integer :: iel,mdle
         if (allocated(ELEM_ORDER)) deallocate(ELEM_ORDER)
         if (allocated(ELEM_SUBD))  deallocate(ELEM_SUBD)
         allocate(ELEM_ORDER(NRELES))
         allocate(ELEM_SUBD(NRELES))
         mdle = 0; NRELES_SUBD = 0
         do iel=1,NRELES
            call nelcon(mdle, mdle)
            ELEM_ORDER(iel) = mdle
            if (NODES(mdle)%subd .eq. RANK) then
               NRELES_SUBD = NRELES_SUBD + 1
               ELEM_SUBD(NRELES_SUBD) = mdle
            endif
         enddo
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine open_history_file(fp)
      character(len=*) :: fp
      open(unit=NHIST,file=fp,   &
           form='formatted',access='sequential',status='unknown')
      end subroutine
!
!-----------------------------------------------------------------------
!
      subroutine close_history_file
      write(NHIST,*) '0 0'
      close(NHIST)
      end subroutine
!
!-----------------------------------------------------------------------
!
!  ...determine number of dof for a higher order node
      subroutine find_ndof(Nod, NdofH,NdofE,NdofV,NdofQ)
!
      integer Nod,NdofH,NdofE,NdofV,NdofQ
!
      call ndof_nod(NODES(Nod)%type,NODES(Nod)%order, NdofH,NdofE,NdofV,NdofQ)
!
      end subroutine find_ndof
!
!-----------------------------------------------------------------------
!
!  ...find number of H1,H(curl),H(div),L2 variables supported by a node
      subroutine find_nvar(Nod, NvarH,NvarE,NvarV,NvarQ)
!
      integer :: Nod,NvarH,NvarE,NvarV,NvarQ
      integer, dimension(NRINDEX) :: idx
      integer :: k
!
      NvarH=0 ; NvarE=0 ; NvarV=0 ; NvarQ=0
!
      call get_index(Nod, idx)
      do k=1,NRINDEX
        select case (idx(k))
        case (2); nvarH = nvarH + 1
        case (4); nvarE = nvarE + 1
        case (6); nvarV = nvarV + 1
        case (8); nvarQ = nvarQ + 1
        end select
      end do
!
      end subroutine find_nvar
!
!-----------------------------------------------------------------------
!
!  ...determine number of sons for a higher order node
      subroutine find_nsons(Nod, Nrsons)
!
      integer Nod,Nrsons
!
      Nrsons = NODES(Nod)%nr_sons
!
      end subroutine find_nsons
!
!-----------------------------------------------------------------------
!
!  ...return i-th son of Nod
      function Son(Nod,I)
      integer Nod
      integer I
      integer Son
!
      Son = NODES(Nod)%first_son+I-1
!
      end function
!
!-----------------------------------------------------------------------
!
!  ...get index for a node
      subroutine get_index(Nod, Indexd)
!
      integer Indexd(NRINDEX)
!
      call decodLong(NODES(Nod)%index,10,NRINDEX, Indexd)
!
      end subroutine get_index
!
!-----------------------------------------------------------------------
!
!  ...allocate memory for data structure
      subroutine allocds
!
      if (allocated(ELEMS)) then
        write(*,*) 'allocds: WARNING !! DATA STRUCTURE ARRAYS', &
                   ' HAVE NOT BEEN DEALLOCATED'
        call deallocds
      endif
!
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
!
      allocate(NODES(MAXNODS))
!$OMP PARALLEL DO
      do nod=1,MAXNODS
        NODES(nod)%type = 'none'
        NODES(nod)%case = 0
        NODES(nod)%index = 0
        NODES(nod)%order = 0
        NODES(nod)%bcond = nod+1
        NODES(nod)%ref_kind = 0
        NODES(nod)%father = 0
        NODES(nod)%first_son = 0
        NODES(nod)%nr_sons = 0
        NODES(nod)%geom_interf = 0
#if DEBUG_MODE
        NODES(nod)%error = 0.d0
#endif
        nullify (NODES(nod)%coord)
        nullify (NODES(nod)%zdofH)
        nullify (NODES(nod)%zdofE)
        nullify (NODES(nod)%zdofV)
        nullify (NODES(nod)%zdofQ)
      enddo
!$OMP END PARALLEL DO
      NODES(MAXNODS)%bcond = 0
      NPNODS=1
!
!
      end subroutine allocds
!
!-----------------------------------------------------------------------
!
!  ...deallocate data structure arrays
      subroutine deallocds
!
      do nel=1,NRELIS
        if (associated(ELEMS(nel)%physics)) deallocate(ELEMS(nel)%physics)
        if (associated(ELEMS(nel)%bcond))   deallocate(ELEMS(nel)%bcond)
        if (associated(ELEMS(nel)%nodes))   deallocate(ELEMS(nel)%nodes)
        if (associated(ELEMS(nel)%neig))    deallocate(ELEMS(nel)%neig)
      enddo
      deallocate(ELEMS)
!
      do nod=1,MAXNODS
        if (Associated(NODES(nod)%coord))   deallocate(NODES(nod)%coord)
        if (Associated(NODES(nod)%zdofH))   deallocate(NODES(nod)%zdofH)
        if (Associated(NODES(nod)%zdofE))   deallocate(NODES(nod)%zdofE)
        if (Associated(NODES(nod)%zdofV))   deallocate(NODES(nod)%zdofV)
        if (Associated(NODES(nod)%zdofQ))   deallocate(NODES(nod)%zdofQ)
      enddo
      deallocate(NODES)
!
      if (allocated(ELEM_ORDER)) deallocate(ELEM_ORDER)
!
      end subroutine deallocds
!
!-----------------------------------------------------------------------
!
!  ...dump out hp3d data structure
      subroutine dumpout_hp3d(Dump_file)
!
      character(len=15) :: Dump_file
      ndump=31
!!    kyungjoo
!!      open(unit=ndump,file=Dump_file,
!!     .     buffered='yes',blocksize = 65536,
!!     .     form='formatted',access='sequential',status='unknown')
      open(unit=ndump,file=Dump_file,  &
           form='formatted',access='sequential',status='unknown')

!
      write(ndump,*) NRELIS,NRELES,NRNODS
      write(ndump,*) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
      write(ndump,*) MAXNODS,NPNODS
!
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
!
      do nod=1,NRNODS
        write(ndump,*) NODES(nod)%type
        write(ndump,*) NODES(nod)%case
        write(ndump,*) NODES(nod)%index
        write(ndump,*) NODES(nod)%order
        write(ndump,*) NODES(nod)%bcond
        write(ndump,*) NODES(nod)%ref_kind
        write(ndump,*) NODES(nod)%father
        write(ndump,*) NODES(nod)%first_son
        write(ndump,*) NODES(nod)%nr_sons
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
#if DEBUG_MODE
        write(ndump,*) NODES(nod)%error
#endif
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
!
      close(ndump)
!
      end subroutine dumpout_hp3d
!
!-----------------------------------------------------------------------
!  ...dump in hp3d data structure
      subroutine dumpin_hp3d(Dump_file)
!
      character(len=15) :: Dump_file
      character(len=20) :: type
      integer           :: npnods_loc
!
      if (allocated(ELEMS).or.allocated(NODES)) call deallocds
!
      ndump=31
      open(unit=ndump,file=Dump_file,  &
           form='formatted',access='sequential',status='unknown')
!
      read(ndump,*) NRELIS,NRELES,NRNODS
      read(ndump,*) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
      read(ndump,*) MAXNODS,NPNODS
      npnods_loc = NPNODS
      call allocds
      NPNODS = npnods_loc
!
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
!
      do nod=1,NRNODS
        read(ndump,*) NODES(nod)%type
        read(ndump,*) NODES(nod)%case
        read(ndump,*) NODES(nod)%index
        read(ndump,*) NODES(nod)%order
        read(ndump,*) NODES(nod)%bcond
        read(ndump,*) NODES(nod)%ref_kind
        read(ndump,*) NODES(nod)%father
        read(ndump,*) NODES(nod)%first_son
        read(ndump,*) NODES(nod)%nr_sons
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
#if DEBUG_MODE
        read(ndump,*) NODES(nod)%error
#endif
!
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
!
      close(ndump)
!
      end subroutine dumpin_hp3d
!
!-----------------------------------------------------------------------
      subroutine add_dirichlet_to_list(Iboundary)
      integer loc
      loc = 0
      call locate(Iboundary, DIRICHLET_LIST, NR_DIRICHLET_LIST, loc)
      if (loc.eq.0) then
        NR_DIRICHLET_LIST = NR_DIRICHLET_LIST + 1
        DIRICHLET_LIST(NR_DIRICHLET_LIST) = Iboundary
      end if
      end subroutine

!-----------------------------------------------------------------------
      subroutine add_dirichlet_homogeneous_to_list(Iboundary)
      integer loc1,loc2
      loc1 = 0
      loc2 = 0
      call locate(Iboundary,DIRICHLET_LIST,NR_DIRICHLET_LIST, loc1)
      call locate(Iboundary,DIRICHLET_HOMOGENEOUS_LIST, &
                  NR_DIRICHLET_HOMOGENEOUS_LIST, loc2)
      if (loc2.eq.0) then
        NR_DIRICHLET_HOMOGENEOUS_LIST = NR_DIRICHLET_LIST + 1
        DIRICHLET_HOMOGENEOUS_LIST(NR_DIRICHLET_HOMOGENEOUS_LIST) = Iboundary
!       any homogenous B.C flag is also a Dirichlet flag
        if (loc1.eq.0) then
          NR_DIRICHLET_LIST = NR_DIRICHLET_LIST + 1
          DIRICHLET_LIST(NR_DIRICHLET_LIST) = Iboundary
        end if ! loc1
      end if !loc2
      end subroutine
!

!-----------------------------------------------------------------------
!  ...reset visitation flags for all nodes
      subroutine reset_visit
!
!$OMP PARALLEL DO
      do i=1,NRNODS
        NODES(i)%visit = 0
      enddo
!$OMP END PARALLEL DO
!
      end subroutine reset_visit
!
!
!  ...get visitation flag for a node
      subroutine get_visit(Nod, Vis)
!
      integer, intent(in)  :: Nod
      integer, intent(out) :: Vis
!
      Vis = NODES(Nod)%visit
!
      end subroutine

!  ...set visitation flag of a node
      subroutine set_visit(Nod)
!
      integer, intent(in) :: Nod
!
      NODES(Nod)%visit = 1
!
      end subroutine
!
!-----------------------------------------------------------------------
!
!  ...get current subdomain of a node
      subroutine get_subd(Nod, Subd)
!
      integer, intent(in)  :: Nod
      integer, intent(out) :: Subd
!
      Subd = NODES(Nod)%subd
!
      end subroutine
!
!  ...set new subdomain of a node
      subroutine set_subd(Nod,Subd)
!
      integer, intent(in) :: Nod,Subd
!
      NODES(Nod)%subd = Subd
!
      end subroutine
!
!-----------------------------------------------------------------------
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
          call locate(ibc(iphys),DIRICHLET_LIST,NR_DIRICHLET_LIST, loc)
          if (loc.ne.0) then
            Is_dirichlet = .true.
          endif
        endif
      enddo
!
      end function
!----------------------------------------------------------------------
      function Is_dirichlet_homogeneous(Nod)
      integer Nod
      integer ibc(NR_PHYSA), loc
      logical Is_dirichlet_homogeneous

      call decod(NODES(Nod)%bcond,10,NR_PHYSA, ibc)
      Is_dirichlet_homogeneous = .false.
      do iphys=1,NR_PHYSA
        call locate(ibc(iphys),DIRICHLET_HOMOGENEOUS_LIST,  &
                    NR_DIRICHLET_HOMOGENEOUS_LIST, loc)
        if (loc.ne.0) then
          Is_dirichlet_homogeneous = .true.
        endif
      enddo
!
      end function
!----------------------------------------------------------------------
      subroutine check_dirichlet_homogeneous(istat)
      integer :: istat

!     local variables
      integer :: ibc_number,ii

!     we check if any flag number in list of homogeneous Dirichlet b.c.
!     is also in the list
      istat=0
      do ii=1,SIZE(DIRICHLET_HOMOGENEOUS_LIST)
        ibc_number=DIRICHLET_HOMOGENEOUS_LIST(ii)
          call locate(ibc_number,DIRICHLET_LIST,NR_DIRICHLET_LIST, loc)
          if (loc.eq.0) then
            istat = 1
            write(*,*) 'ERROR check_dirichlet_homogeneous: ',  &
             'ibc=',ibc_number,'is not in DIRICHLET_LIST array'
            return
          endif
      end do

      end subroutine
!----------------------------------------------------------------------
      function Is_right_handed(Mdle)
      integer :: Mdle
      logical :: Is_right_handed
      integer :: i, nod
      real*8  :: v(3,4), a(3,3), val
      select case(ELEMS(Mdle)%type)
      case('bric','pris','pyra')
        Is_right_handed = .true.
      case('tetr')
!
        do i=1,4
          nod = ELEMS(Mdle)%nodes(i)
          v(1:3,i) = NODES(nod)%coord(1:3,1)
        enddo
!
        do i=1,3
          a(1:3,i) = v(1:3,i+1) - v(1:3,1)
        enddo
!
        call mixed_product(a(1:3,1), a(1:3,2), a(1:3,3), val)
        Is_right_handed = (val > 0.d0)
      end select
      end function
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      function Is_root(Nod)
      integer Nod
      logical Is_root
      if (NODES(Nod)%father.lt.0) then
        Is_root = .TRUE.
      else
        Is_root = .FALSE.
      endif
      end function
!-----------------------------------------------------------------------
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
!
end module data_structure3D
