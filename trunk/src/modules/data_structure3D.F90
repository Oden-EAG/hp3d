!----------------------------------------------------------------------
!> @brief   Defines data structure arrays
!> @date    Feb 2023
!----------------------------------------------------------------------
module data_structure3D
!
      use physics
      use parameters
      use element_data
      use mpi_param, only: RANK
      implicit none
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
!  .....etype - BRIC, TETR, PRIS, PYRA
        integer :: etype
!
!  .....number of physics attributes supported BY the element
        integer :: nrphysics
!
!  .....list   of physics attributes supported BY the element
        character(len=5), dimension(:), pointer :: physics
!
!  .....array of boundary condition nicknames
!       each entry specifies boundary conditions for one particular component;
!       each entry is decimal-encoded per element face (1 BC flag per face)
!       Values (per component and face):
!         0   - No BC
!         1   - Dirichlet BC
!         2-9 - User-customizable BCs
        integer, dimension(:), pointer :: bcond
!
!  .....element nodal connectivities: vertices,edges,faces, middle node
        integer, dimension(:), pointer :: nodes
!
!  .....corresponding orientations for edges
        integer :: edge_orient
!
!  .....corresponding orientations for faces
        integer :: face_orient
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
!  .....ntype - VERT, MEDG, MDLT, MDLQ, MDLB, MDLN, MDLP, MDLD
        integer          :: ntype
!
!  .....case number indicating what physical attributes are supported
!       (binary-encoded per physics variable)
        integer          :: case
!
!  .....order of approximation (decimal-encoded per direction (x,y,z))
        integer          :: order
!
!  .....boundary condition flag (binary-encoded per component)
!       0: component DOFs treated as unknowns
!       1: component DOFs treated as Dirichlet DOFs
        integer          :: bcond
!
!  .....father node
        integer          :: father
!
!  .....node sons
        integer          :: first_son
        integer          :: nr_sons
!
!  .....refinement flag (decimal-encoded per direction (x,y,z))
        integer          :: ref_kind
!
!  .....visitation flag
        integer          :: visit
!
!  .....activation flag
        logical          :: act
!
!  .....subdomain number (distributed mesh)
        integer          :: subd
!
!  .....geometry and solution degrees of freedom
        type(dof_data), pointer :: dof
!
#if DEBUG_MODE
!
!  .....error
!       0   - scalar error
!       1-3 - gradient of any vector component-wise error
!       4   - additional info
        real(8), dimension(0:4,0:4) :: error
#endif
!
      endtype node
!
!----------------------------------------------------------------------
!  DOF DATA                                                            |
!----------------------------------------------------------------------
      type dof_data
!
!  .....geometry dof
        real(8), dimension(:,:), pointer :: coord
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
      endtype dof_data
!
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
      subroutine update_ELEM_ORDER
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
      end subroutine update_ELEM_ORDER
!
!-----------------------------------------------------------------------
!
      subroutine open_history_file(fp)
      character(len=*) :: fp
      open(unit=NHIST,file=fp,   &
           form='formatted',access='sequential',status='unknown')
      end subroutine open_history_file
!
!-----------------------------------------------------------------------
!
      subroutine close_history_file
      write(NHIST,*) '0 0'
      close(NHIST)
      end subroutine close_history_file
!
!-----------------------------------------------------------------------
!
!  ...determine number of dof for a higher order node
      subroutine find_ndof(Nod, NdofH,NdofE,NdofV,NdofQ)
!
      integer, intent(in)  :: Nod
      integer, intent(out) :: NdofH,NdofE,NdofV,NdofQ
!
      call ndof_nod(NODES(Nod)%ntype,NODES(Nod)%order, NdofH,NdofE,NdofV,NdofQ)
!
      end subroutine find_ndof
!
!-----------------------------------------------------------------------
!
!  ...find number of H1,H(curl),H(div),L2 variables supported by a node
      subroutine find_nvar(Nod, NvarH,NvarE,NvarV,NvarQ)
!
      integer, intent(in)  :: Nod
      integer, intent(out) :: NvarH,NvarE,NvarV,NvarQ
      integer :: idx(NRINDEX)
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
      integer, intent(in)  :: Nod
      integer, intent(out) :: Nrsons
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
      end function Son
!
!-----------------------------------------------------------------------
!
!  ...allocate memory for data structure
      subroutine allocds
!
      integer :: nel,nod
!
      if (allocated(ELEMS)) then
        write(*,*) 'allocds: WARNING !! DATA STRUCTURE ARRAYS', &
                   ' HAVE NOT BEEN DEALLOCATED'
        call deallocds
      endif
!
      allocate(ELEMS(NRELIS))
      do nel=1,NRELIS
        ELEMS(nel)%etype = 0
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
        NODES(nod)%ntype = 0
        NODES(nod)%case = 0
        NODES(nod)%order = 0
        NODES(nod)%act = .false.
        NODES(nod)%subd = -1
        NODES(nod)%bcond = nod+1
        NODES(nod)%ref_kind = 0
        NODES(nod)%father = 0
        NODES(nod)%first_son = 0
        NODES(nod)%nr_sons = 0
        nullify (NODES(nod)%dof)
#if DEBUG_MODE
        NODES(nod)%error = 0.d0
#endif
      enddo
!$OMP END PARALLEL DO
      NODES(MAXNODS)%bcond = 0
      NPNODS=1
!
      end subroutine allocds
!
!-----------------------------------------------------------------------
!
!  ...deallocate data structure arrays
      subroutine deallocds
!
      integer :: nel,nod
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
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%coord)) deallocate(NODES(nod)%dof%coord)
          if (associated(NODES(nod)%dof%zdofH)) deallocate(NODES(nod)%dof%zdofH)
          if (associated(NODES(nod)%dof%zdofE)) deallocate(NODES(nod)%dof%zdofE)
          if (associated(NODES(nod)%dof%zdofV)) deallocate(NODES(nod)%dof%zdofV)
          if (associated(NODES(nod)%dof%zdofQ)) deallocate(NODES(nod)%dof%zdofQ)
          deallocate(NODES(nod)%dof)
        endif
      enddo
      deallocate(NODES)
!
      if (allocated(ELEM_ORDER)) deallocate(ELEM_ORDER)
!
      end subroutine deallocds
!
!-----------------------------------------------------------------------
!
!  ...increase MAXNODS
      subroutine increase_MAXNODS
!
      type(node), allocatable :: NODES_NEW(:)
      integer :: MAXNODS_NEW,nod
!
      if (.not. allocated(NODES)) then
         write(*,*) 'increase_MAXNODS: NODES not allocated. returning...'
         return
      endif
      if (NPNODS .ne. 0 .or. NRNODS .lt. MAXNODS) then
         write(*,*) 'increase_MAXNODS: INCONSISTENCY. MAXNODS,NRNODS,NPNODS=',&
                                                      MAXNODS,NRNODS,NPNODS
         return
      endif
!
!  ...determine size of new NODES array
      MAXNODS_NEW = 2*MAXNODS
!
!  ...allocate new NODES array twice the size of the old
      allocate(NODES_NEW(MAXNODS_NEW))
!
!  ...copy data from old NODES array into the new one
      NODES_NEW(1:MAXNODS) = NODES(1:MAXNODS)
!
      !$OMP PARALLEL DO
      do nod=MAXNODS+1,MAXNODS_NEW
        NODES_NEW(nod)%ntype = 0
        NODES_NEW(nod)%case = 0
        NODES_NEW(nod)%order = 0
        NODES_NEW(nod)%act = .false.
        NODES_NEW(nod)%subd = -1
        NODES_NEW(nod)%bcond = nod+1
        NODES_NEW(nod)%ref_kind = 0
        NODES_NEW(nod)%father = 0
        NODES_NEW(nod)%first_son = 0
        NODES_NEW(nod)%nr_sons = 0
        nullify (NODES_NEW(nod)%dof)
#if DEBUG_MODE
        NODES_NEW(nod)%error = 0.d0
#endif
      enddo
      !$OMP END PARALLEL DO
      NODES_NEW(MAXNODS_NEW)%bcond = 0
!
!  ...move NODES pointer to the new array,
!     and deallocate the old array
      call move_alloc(NODES_NEW, NODES)
      NPNODS  = MAXNODS+1
      MAXNODS = MAXNODS_NEW
!
      end subroutine increase_MAXNODS
!
!-----------------------------------------------------------------------
!
!  ...dump out hp3d data structure
      subroutine dumpout_hp3d(Dump_file)
!
      character(len=15) :: Dump_file
      integer :: nel,nod,nn,nn1,nn2
      integer :: ndump
      ndump=31
!
      open(unit=ndump,file=Dump_file,  &
           form='formatted',access='sequential',status='unknown')
!
      write(ndump,*) NRELIS,NRELES,NRNODS
      write(ndump,*) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
      write(ndump,*) MAXNODS,NPNODS
!
      do nel=1,NRELIS
        write(ndump,*) ELEMS(nel)%etype
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
        write(ndump,*) NODES(nod)%ntype
        write(ndump,*) NODES(nod)%case
        write(ndump,*) NODES(nod)%order
        write(ndump,*) NODES(nod)%bcond
        write(ndump,*) NODES(nod)%ref_kind
        write(ndump,*) NODES(nod)%father
        write(ndump,*) NODES(nod)%first_son
        write(ndump,*) NODES(nod)%nr_sons
        write(ndump,*) NODES(nod)%visit
        write(ndump,*) NODES(nod)%act
        if (associated(NODES(nod)%dof)) then
          write(ndump,*) 1
        else
          write(ndump,*) 0
        endif
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%coord)) then
            nn1 = ubound(NODES(nod)%dof%coord,1)
            nn2 = ubound(NODES(nod)%dof%coord,2)
            write(ndump,*) nn1, nn2
            write(ndump,*) NODES(nod)%dof%coord
          else
            write(ndump,*) 0 , 0
          endif
        else
          write(ndump,*) 0 , 0
        endif
#if DEBUG_MODE
        write(ndump,*) NODES(nod)%error
#endif
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%zdofH)) then
            nn1 = ubound(NODES(nod)%dof%zdofH,1)
            nn2 = ubound(NODES(nod)%dof%zdofH,2)
            write(ndump,*) nn1, nn2
            write(ndump,*) NODES(nod)%dof%zdofH
          else
            write(ndump,*) 0 , 0
          endif
        else
          write(ndump,*) 0 , 0
        endif
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%zdofE)) then
            nn1 = ubound(NODES(nod)%dof%zdofE,1)
            nn2 = ubound(NODES(nod)%dof%zdofE,2)
            write(ndump,*) nn1, nn2
            write(ndump,*) NODES(nod)%dof%zdofE
          else
            write(ndump,*) 0 , 0
          endif
        else
          write(ndump,*) 0 , 0
        endif
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%zdofV)) then
            nn1 = ubound(NODES(nod)%dof%zdofV,1)
            nn2 = ubound(NODES(nod)%dof%zdofV,2)
            write(ndump,*) nn1, nn2
            write(ndump,*) NODES(nod)%dof%zdofV
          else
            write(ndump,*) 0 , 0
          endif
        else
          write(ndump,*) 0 , 0
        endif
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%zdofQ)) then
            nn1 = ubound(NODES(nod)%dof%zdofQ,1)
            nn2 = ubound(NODES(nod)%dof%zdofQ,2)
            write(ndump,*) nn1, nn2
            write(ndump,*) NODES(nod)%dof%zdofQ
          else
            write(ndump,*) 0 , 0
          endif
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
      integer :: npnods_loc,nel,nod,nn,nn1,nn2,i
      integer :: ndump
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
        read(ndump,*) ELEMS(nel)%etype
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
        read(ndump,*) NODES(nod)%ntype
        read(ndump,*) NODES(nod)%case
        read(ndump,*) NODES(nod)%order
        read(ndump,*) NODES(nod)%bcond
        read(ndump,*) NODES(nod)%ref_kind
        read(ndump,*) NODES(nod)%father
        read(ndump,*) NODES(nod)%first_son
        read(ndump,*) NODES(nod)%nr_sons
        read(ndump,*) NODES(nod)%visit
        read(ndump,*) NODES(nod)%act
        read(ndump,*) nn1
        if (nn1.eq.1) then
          allocate(NODES(nod)%dof)
        else
          nullify(NODES(nod)%dof)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%dof%coord(nn1,nn2))
          read(ndump,*) NODES(nod)%dof%coord
        else
          if(associated(NODES(nod)%dof)) nullify(NODES(nod)%dof%coord)
        endif
#if DEBUG_MODE
        read(ndump,*) NODES(nod)%error
#endif
!
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%dof%zdofH(nn1,nn2))
          read(ndump,*) NODES(nod)%dof%zdofH
        else
          if(associated(NODES(nod)%dof)) nullify(NODES(nod)%dof%zdofH)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%dof%zdofE(nn1,nn2))
          read(ndump,*) NODES(nod)%dof%zdofE
        else
          if(associated(NODES(nod)%dof)) nullify(NODES(nod)%dof%zdofE)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%dof%zdofV(nn1,nn2))
          read(ndump,*) NODES(nod)%dof%zdofV
        else
          if(associated(NODES(nod)%dof)) nullify(NODES(nod)%dof%zdofV)
        endif
        read(ndump,*) nn1, nn2
        if ((nn1.gt.0).and.(nn2.gt.0)) then
          allocate(NODES(nod)%dof%zdofQ(nn1,nn2))
          read(ndump,*) NODES(nod)%dof%zdofQ
        else
          if(associated(NODES(nod)%dof)) nullify(NODES(nod)%dof%zdofQ)
        endif
      enddo
!
      close(ndump)
!
      end subroutine dumpin_hp3d
!
!-----------------------------------------------------------------------
!
      subroutine add_dirichlet_to_list(Iboundary)
      integer, intent(in) :: Iboundary
      integer :: loc
      loc = 0
      call locate(Iboundary, DIRICHLET_LIST, NR_DIRICHLET_LIST, loc)
      if (loc.eq.0) then
        NR_DIRICHLET_LIST = NR_DIRICHLET_LIST + 1
        DIRICHLET_LIST(NR_DIRICHLET_LIST) = Iboundary
      end if
      end subroutine add_dirichlet_to_list
!
!-----------------------------------------------------------------------
!
      subroutine add_dirichlet_homogeneous_to_list(Iboundary)
      integer, intent(in) :: Iboundary
      integer :: loc1,loc2
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
      end subroutine add_dirichlet_homogeneous_to_list
!
!-----------------------------------------------------------------------
!
!  ...reset visitation flags for all nodes
      subroutine reset_visit
      integer :: i
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
      end subroutine get_visit

!  ...set visitation flag of a node
      subroutine set_visit(Nod)
!
      integer, intent(in) :: Nod
!
      NODES(Nod)%visit = 1
!
      end subroutine set_visit
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
      end subroutine get_subd
!
!  ...set new subdomain of a node
      subroutine set_subd(Nod,Subd)
!
      integer, intent(in) :: Nod,Subd
!
      NODES(Nod)%subd = Subd
!
      end subroutine set_subd
!
!-----------------------------------------------------------------------
!
      function Is_Dirichlet(Nod)
      logical Is_Dirichlet
      integer Nod
      integer ibc(NRINDEX), ic, iphys, ivar
!
      call decod(NODES(Nod)%bcond,2,NRINDEX, ibc)
      Is_Dirichlet = .false.
!
      ic = 0
!
!  ...check Dirichlet flags for all variable types
      do iphys=1,NR_PHYSA
         select case(DTYPE(iphys))
!     ...H1 checks all
         case('contin')
            continue
!     ...skip checking vertices of H(curl) variables
         case('tangen')
            if (NODES(nod)%ntype .eq. VERT) then
               ic = ic + NR_COMP(iphys)
               cycle
            endif
!     ...skip checking vertices and edges of H(div) variables
         case('normal')
            if (     NODES(nod)%ntype .eq. VERT   &
                .or. NODES(nod)%ntype .eq. MEDG ) then
               ic = ic + NR_COMP(iphys)
               cycle
            endif
!     ...skip checking all on L2
         case('discon')
            ic = ic + NR_COMP(iphys)
            cycle
         end select
!
         do ivar=1,NR_COMP(iphys)
            ic = ic + 1
            if (ibc(ic).eq.1) Is_Dirichlet = .true.
         enddo
!
      enddo
!
      end function Is_Dirichlet
!
!-----------------------------------------------------------------------
!
      function Is_Dirichlet_attr(Nod,D_type)
      logical Is_Dirichlet_attr
      integer Nod
      character(6) D_type
      integer ibc(NRINDEX), ic, iphys, ivar
!
      call decod(NODES(Nod)%bcond,2,NRINDEX, ibc)
      Is_Dirichlet_attr = .false.
      ic = 0
      do iphys=1,NR_PHYSA
        if (DTYPE(iphys).eq.D_type) then
          do ivar=1,NR_COMP(iphys)
            ic = ic + 1
            if (ibc(ic).eq.1) Is_Dirichlet_attr = .true.
          enddo
        else
          ic = ic + NR_COMP(iphys)
        endif
      enddo
!
      end function Is_Dirichlet_attr
!
!----------------------------------------------------------------------
!
      function Is_right_handed(Mdle)
      logical :: Is_right_handed
      integer :: Mdle
      integer :: i, nod
      real(8) :: v(3,4), a(3,3), val
!
      select case(ELEMS(Mdle)%etype)
      case(BRIC,PRIS,PYRA)
         Is_right_handed = .true.
      case(TETR)
         do i=1,4
            nod = ELEMS(Mdle)%nodes(i)
            v(1:3,i) = NODES(nod)%dof%coord(1:3,1)
         enddo
         do i=1,3
            a(1:3,i) = v(1:3,i+1) - v(1:3,1)
         enddo
         call mixed_product(a(1:3,1), a(1:3,2), a(1:3,3), val)
         Is_right_handed = (val > 0.d0)
      case default
         write(*,*) 'Is_right_handed'; stop
      end select
      end function Is_right_handed
!
!-----------------------------------------------------------------------
!
      function Is_active(Nod)
      logical Is_active
      integer Nod
      Is_active = NODES(Nod)%act
      end function Is_active
!
!-----------------------------------------------------------------------
!
      function Is_inactive(Nod)
      logical Is_inactive
      integer Nod
      Is_inactive = .not. NODES(Nod)%act
      end function Is_inactive
!
!-----------------------------------------------------------------------
!
      function Is_leaf(Nod)
      logical Is_leaf
      integer Nod
      select case(NODES(Nod)%ref_kind)
         case(0);       Is_leaf = .TRUE.
         case default;  Is_leaf = .FALSE.
      end select
      end function Is_leaf
!
!-----------------------------------------------------------------------
!
      function Is_root(Nod)
      logical Is_root
      integer Nod
      if (NODES(Nod)%father.lt.0) then
         Is_root = .TRUE.
      else
         Is_root = .FALSE.
      endif
      end function Is_root
!
!-----------------------------------------------------------------------
!
      function Is_middle(Nod)
      logical Is_middle
      integer Nod
      select case(NODES(Nod)%ntype)
         case(MDLB,MDLP,MDLN,MDLD); Is_middle = .TRUE.
         case default;              Is_middle = .FALSE.
      end select
      end function Is_middle
!
end module data_structure3D
