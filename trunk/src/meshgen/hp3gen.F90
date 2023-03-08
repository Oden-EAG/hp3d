!> @brief routine generates an initial FE mesh interfacing with GMP
!!
!! @param[in]  Fp - Path to input file (physics attributes)
!!
!> @date Feb 2023
subroutine hp3gen(Fp)
  !
  use error
  use GMP
  use element_data
  use data_structure3D
  !
  implicit none
  !
  ! input argument
  character(len=*) :: Fp
  !
  ! element order (temporary)
  integer, allocatable :: nelem_order(:)
  !
  ! node BC flags
  integer, allocatable :: ibc_nod(:)
  !
  ! physical attributes of a node
  character(len=5), allocatable :: phys_vect(:)
  !
  ! edge, face orientations
  integer :: nedge_orient(12),nface_orient(6)
  !
  ! neighbors of a point or curve
  integer, parameter :: max_neig=150
  integer :: neigbl(max_neig)
  !
  !  ...BC for element faces an attribute
  integer :: ibc_elem(6)
  !
  !  ...faces adjacent to a vertex
  integer :: nofaces(4)
  real(8) :: x(NDIMEN)
  !
  !  ...physical attributes for an element or node
  character(len=5) :: phys
  !
  !  ...element type
  integer :: etype
  !
  !  ...node type
  integer :: ntype
  !
  !  ...auxiliary
  integer :: nel, npri, nh, ntet, npyr, np, iv, is, ifc, ie, mdle
  integer :: nt, nrbl, nbl, nr, lab, nord, nod_new, nbcond, nod, nrfaces
  integer :: nb, nc, i, ib, iii, istat, icase, iphys, num, number, subd
  integer :: ivar,nvar
  logical :: iact
  !
#if DEBUG_MODE
  integer :: iprint = -1
  integer :: iprint_vert = 0
  integer :: iprint_edge = 0
  integer :: iprint_face = 0
  integer :: iprint_mdle = 0
#endif
  !----------------------------------------------------------------------
  !
#if DEBUG_MODE
  if (iprint .eq. 1) then
     write(*,*) 'hp3gen: DEBUGGING'
  endif
#endif
  !
  x(1:NDIMEN) = 0.d0
  !
  call read_input(Fp)

  !  ...total number of elements and middle nodes
  NRELIS = NRPRISM+NRHEXAS+NRTETRA+NRPYRAM
  NRELES = NRELIS
  !
  !  ...allocate physics and data structure arrays
  call allocds
  !
  !  ...initialize number of H1,H(curl),H(div),L2 dofs
  NRDOFSH=0 ; NRDOFSE=0 ; NRDOFSV=0 ; NRDOFSQ=0
  !
  allocate(nelem_order(NRELIS),ibc_nod(NRINDEX),phys_vect(NR_PHYSA), &
           stat=istat)
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  !
  !  ...elements are generated in the order listed above: prisms, hexas,
  !     tets, pyramids
  !
  !  ...nodes are enumerated in the order: middle nodes (to make the
  !     element and middle number for an initial mesh coincide),
  !     vertex nodes (points), edge nodes (curves), face nodes (rectangles
  !     and triangles next)
  !
  !  ...generate initial mesh elements.............
  nel=0
  !
  !  ...loop through prisms
  do npri=1,NRPRISM
     nel=nel+1
     ELEMS(nel)%etype = PRIS
     allocate(ELEMS(nel)%nodes(21))
     allocate(ELEMS(nel)%neig(5)); ELEMS(nel)%neig(1:5)=0
     do iv=1,6
        ELEMS(nel)%nodes(iv) = NRELIS+PRISMS(npri)%VertNo(iv)
     enddo
     nedge_orient(1:9)=0
     do ie=1,9
        if (PRISMS(npri)%EdgeNo(ie).lt.0) nedge_orient(ie)=1
        ELEMS(nel)%nodes(6+ie) = NRELIS+NRPOINT &
                               + abs(PRISMS(npri)%EdgeNo(ie))
     enddo
     do ifc=1,2
        call decode(PRISMS(npri)%FigNo(ifc), nt,nface_orient(ifc))
        ELEMS(nel)%nodes(15+ifc) = NRELIS+NRPOINT+NRCURVE + NRRECTA+nt
        do is=1,2
           if (TRIANGLES(nt)%BlockNo(is).ne.npri*10+1) then
              call decode(TRIANGLES(nt)%BlockNo(is), nbl,lab)
              select case(lab)
              case(1)
                 ELEMS(nel)%neig(ifc) = nbl
              case(3)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+nbl
              case(4)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+NRTETRA+nbl
              end select
           endif
        enddo
     enddo
     do ifc=3,5
        call decode(PRISMS(npri)%FigNo(ifc), nr,nface_orient(ifc))
        ELEMS(nel)%nodes(15+ifc) = NRELIS+NRPOINT+NRCURVE+nr
        do is=1,2
           if (RECTANGLES(nr)%BlockNo(is).ne.npri*10+1) then
              call decode(RECTANGLES(nr)%BlockNo(is), nbl,lab)
              select case(lab)
              case(1)
                 ELEMS(nel)%neig(ifc) = nbl
              case(2)
                 ELEMS(nel)%neig(ifc) = NRPRISM+nbl
              case(4)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+NRTETRA+nbl
              end select
           endif
        enddo
     enddo
     ELEMS(nel)%nodes(21) = nel
     call encodg(nedge_orient,2,9, ELEMS(nel)%edge_orient)
     call encodg(nface_orient,8,5, ELEMS(nel)%face_orient)
     ELEMS(nel)%GMPblock = npri*10+1
  enddo
#if DEBUG_MODE
  if ((iprint.eq.1).and.(NRPRISM.ne.0)) then
     write(*,*) 'hp3gen: HAVE GENERATED PRISMS'
  endif
#endif
  !
  !  ...loop through hexas
  do nh=1,NRHEXAS
     nel=nel+1
     ELEMS(nel)%etype = BRIC
     allocate(ELEMS(nel)%nodes(27))
     allocate(ELEMS(nel)%neig(6)); ELEMS(nel)%neig(1:6)=0
     do iv=1,8
        ELEMS(nel)%nodes(iv) = NRELIS+HEXAS(nh)%VertNo(iv)
     enddo
     nedge_orient(1:12)=0
     do ie=1,12
        if (HEXAS(nh)%EdgeNo(ie).lt.0) nedge_orient(ie)=1
        ELEMS(nel)%nodes(8+ie) = NRELIS+NRPOINT &
                               + abs(HEXAS(nh)%EdgeNo(ie))
     enddo
     do ifc=1,6
        call decode(HEXAS(nh)%FigNo(ifc), nr,nface_orient(ifc))
        ELEMS(nel)%nodes(20+ifc) = NRELIS+NRPOINT+NRCURVE+nr
        do is=1,2
           if (RECTANGLES(nr)%BlockNo(is).ne.nh*10+2) then
              call decode(RECTANGLES(nr)%BlockNo(is), nbl,lab)
              select case(lab)
              case(1)
                 ELEMS(nel)%neig(ifc) = nbl
              case(2)
                 ELEMS(nel)%neig(ifc) = NRPRISM+nbl
              case(4)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+NRTETRA+nbl
              end select
           endif
        enddo
     enddo
     ELEMS(nel)%nodes(27) = nel
     call encodg(nedge_orient,2,12, ELEMS(nel)%edge_orient)
     call encodg(nface_orient,8,6, ELEMS(nel)%face_orient)
     ELEMS(nel)%GMPblock = nh*10+2
  enddo
#if DEBUG_MODE
  if ((iprint.eq.1).and.(NRHEXAS.ne.0)) then
     write(*,*) 'hp3gen: HAVE GENERATED HEXAS'
  endif
#endif
  !
  !  ...loop through tets
  do ntet=1,NRTETRA
     nel=nel+1
     ELEMS(nel)%etype = TETR
     allocate(ELEMS(nel)%nodes(15))
     allocate(ELEMS(nel)%neig(4)); ELEMS(nel)%neig(1:4)=0
     do iv=1,4
        ELEMS(nel)%nodes(iv) = NRELIS+TETRAS(ntet)%VertNo(iv)
     enddo
     nedge_orient(1:6)=0
     do ie=1,6
        if (TETRAS(ntet)%EdgeNo(ie).lt.0) nedge_orient(ie)=1
        ELEMS(nel)%nodes(4+ie) = NRELIS+NRPOINT &
                               + abs(TETRAS(ntet)%EdgeNo(ie))
     enddo
     do ifc=1,4
        call decode(TETRAS(ntet)%FigNo(ifc), nt,nface_orient(ifc))
        ELEMS(nel)%nodes(10+ifc) = NRELIS+NRPOINT+NRCURVE+NRRECTA+nt
        do is=1,2
           if (TRIANGLES(nt)%BlockNo(is).ne.ntet*10+3) then
              call decode(TRIANGLES(nt)%BlockNo(is), nbl,lab)
              select case(lab)
              case(1)
                 ELEMS(nel)%neig(ifc) = nbl
              case(3)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+nbl
              case(4)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+NRTETRA+nbl
              end select
           endif
        enddo
     enddo
     ELEMS(nel)%nodes(15) = nel
     call encodg(nedge_orient,2,6, ELEMS(nel)%edge_orient)
     call encodg(nface_orient,8,4, ELEMS(nel)%face_orient)
     ELEMS(nel)%GMPblock = ntet*10+3
  enddo
#if DEBUG_MODE
  if ((iprint.eq.1).and.(NRTETRA.ne.0)) then
     write(*,*) 'hp3gen: HAVE GENERATED TETRAS'
  endif
#endif
  !
  !  ...loop through pyramids
  do npyr=1,NRPYRAM
     nel=nel+1
     ELEMS(nel)%etype = PYRA
     allocate(ELEMS(nel)%nodes(19))
     allocate(ELEMS(nel)%neig(5)); ELEMS(nel)%neig(1:5)=0
     do iv=1,5
        ELEMS(nel)%nodes(iv) = NRELIS+PYRAMIDS(npyr)%VertNo(iv)
     enddo
     nedge_orient(1:8)=0
     do ie=1,8
        if (PYRAMIDS(npyr)%EdgeNo(ie).lt.0) nedge_orient(ie)=1
        ELEMS(nel)%nodes(5+ie) = NRELIS+NRPOINT &
                               + abs(PYRAMIDS(npyr)%EdgeNo(ie))
     enddo
     call decode(PYRAMIDS(npyr)%FigNo(1), nr,nface_orient(1))
     ELEMS(nel)%nodes(14) = NRELIS+NRPOINT+NRCURVE+nr
     do is=1,2
        if (RECTANGLES(nr)%BlockNo(is).ne.npyr*10+4) then
           call decode(RECTANGLES(nr)%BlockNo(is), nbl,lab)
           select case(lab)
           case(1)
              ELEMS(nel)%neig(is) = nbl
           case(2)
              ELEMS(nel)%neig(is) = NRPRISM+nbl
           case(4)
              ELEMS(nel)%neig(is) = NRPRISM+NRHEXAS+NRTETRA+nbl
           end select
        endif
     enddo
     do ifc=2,5
        call decode(PYRAMIDS(npyr)%FigNo(ifc), nt,nface_orient(ifc))
        ELEMS(nel)%nodes(13+ifc) = NRELIS+NRPOINT+NRCURVE+NRRECTA+nt
        do is=1,2
           if (TRIANGLES(nt)%BlockNo(is).ne.npyr*10+4) then
              call decode(TRIANGLES(nt)%BlockNo(is), nbl,lab)
              select case(lab)
              case(1)
                 ELEMS(nel)%neig(ifc) = nbl
              case(3)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+nbl
              case(4)
                 ELEMS(nel)%neig(ifc) = NRPRISM+NRHEXAS+NRTETRA+nbl
              end select
           endif
        enddo
     enddo
     ELEMS(nel)%nodes(19) = nel
     call encodg(nedge_orient,2,8, ELEMS(nel)%edge_orient)
     call encodg(nface_orient,8,5, ELEMS(nel)%face_orient)
     ELEMS(nel)%GMPblock = npyr*10+4
  enddo
#if DEBUG_MODE
  if ((iprint.eq.1).and.(NRPYRAM.ne.0)) then
     write(*,*) 'hp3gen: HAVE GENERATED PYRAMIDS'
  endif
#endif
  !
  !----------------------------------------------------------------------
  !
  !  ...set physics flags, boundary conditions flags and order
  !     of approximation for initial mesh elements (user provided routine)
#if DEBUG_MODE
  if (iprint.ne.-1) then
     write(*,*)'CALLING set_initial_mesh'
  endif
#endif
  call set_initial_mesh(nelem_order)
  !
  !----------------------------------------------------------------------
  !
  !  ...initialize number of geometry dofs
  !
  NRNODS=0
  !
  !  ...generate middle nodes
  do nel=1,NRELIS
     nod_new = nel
     !
     !  .....identify the case number for the middle node...................
     call find_case(ELEMS(nel)%nrphysics,ELEMS(nel)%physics, icase)
     !
     nbcond=0
     nord = nelem_order(nel)
     !
     select case(ELEMS(nel)%etype)
     case(BRIC); ntype=MDLB
     case(TETR); ntype=MDLN
     case(PRIS); ntype=MDLP
     case(PYRA); ntype=MDLD
     end select
     !
     subd = -1; iact = .true.
     call nodgen(ntype,icase,nbcond,-nel,nord,subd,iact, nod)
     if (nod_new.ne.nod) then
        write(*,*) 'hp3gen: nod_new,nod = ',nod_new,nod
        stop 1
     endif
#if DEBUG_MODE
     if (iprint_mdle.eq.1) then
        write(*,7004) nod
7004    format('hp3gen: HAVE GENERATED MIDDLE NODE ',i6)
     endif
#endif
  enddo
  !
  deallocate(nelem_order , stat=istat)
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  !
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,*) 'hp3gen: HAVE GENERATED MIDDLE NODES'
  endif
#endif
  !
  !----------------------------------------------------------------------
  !
  !  ...generate vertex nodes
  do np=1,NRPOINT
     nod_new = NRELIS+np; ibc_nod=0
     !
     call find_point_to_block(np,max_neig, nrbl,neigbl)
     !
     !  .....identify the case number for the point........................
     num=0
     !
     !  .....loop through neighboring blocks
     do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)
        !
        !  .......determine element number
        select case(lab)
           !
           !  .......prism
        case(1); nel= nb; etype=PRIS
           !
           !  .......hexahedron
        case(2); nel= NRPRISM+nb; etype=BRIC
           !
           !  .......tetrahedron
        case(3); nel= NRPRISM+NRHEXAS+nb; etype=TETR
           !
           !  .......pyramid
        case(4); nel= NRPRISM+NRHEXAS+NRTETRA+nb; etype=PYRA
        end select
        !
        do iphys=1,ELEMS(nel)%nrphysics
           phys = ELEMS(nel)%physics(iphys)
           !
           !  .........check if on the list
           call locate_char(phys,phys_vect,num, number)
           if (number.eq.0) then
              num=num+1
              phys_vect(num) = phys
           endif
        enddo
        !
        !  .......loop through the element physics attributes and their components
        nvar=0
        do iphys=1,ELEMS(nel)%nrphysics
           phys = ELEMS(nel)%physics(iphys)
           call locate_char(phys,PHYSA,NR_PHYSA, iii)
           !
           !  .........loop through the variable components
           do ivar=1,NR_COMP(iii)
              nvar = nvar+1
              !
              !  .........decode the BC flags for the faces
              call decodg(ELEMS(nel)%bcond(nvar),10,nface(etype), ibc_elem)
              !
              !  .........determine faces adjacent to the vertex
              call locate(nod_new,ELEMS(nel)%nodes(1),nvert(etype), iv)
              call vert_to_faces(etype,iv, nrfaces,nofaces)
              !
              !  .........loop through the faces adjacent to the vertex and
              !           use element face BC flags to establish BC flag for the
              !           point
              do i=1,nrfaces
                 ifc = nofaces(i)
                 call copyBCflag(1,ibc_elem(ifc), ibc_nod(nvar))
              enddo
           enddo
        enddo
     enddo
     !
     call find_case(num,phys_vect, icase)
     !
     !  .....encode BC flags for the node into a single nickname
     call encod(ibc_nod,2,NRINDEX, nbcond)
     !
     subd = -1; iact = .true.
     call nodgen(VERT,icase,nbcond,-nel,1,subd,iact, nod)
     if (nod_new.ne.nod) then
        write(*,*) 'hp3gen: nod_new,nod = ',nod_new,nod
        stop 1
     endif
#if DEBUG_MODE
     if (iprint_vert.eq.1) then
        write(*,7001) nod
7001    format('hp3gen: HAVE GENERATED VERTEX NODE ',i6)
     endif
#endif
  enddo
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,*) 'hp3gen: HAVE GENERATED VERTEX NODES'
  endif
#endif
  !
  !----------------------------------------------------------------------
  !
  !  ...generate mid-edge nodes
  do nc=1,NRCURVE
     nod_new = NRELIS+NRPOINT+nc; ibc_nod=0; nord=MAXP
     !
     call find_curve_to_block(nc,max_neig, nrbl,neigbl)
     !
     !  .....identify the case number and the BC flag for the node
     num=0
!
!    loop through neighboring blocks
     do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)
!
!       determine element number
        select case(lab)
!
!       prism
        case(1); nel=                         nb; etype=PRIS
!
!       hexahedron
        case(2); nel= NRPRISM                +nb; etype=BRIC
!
!       tetrahedron
        case(3); nel= NRPRISM+NRHEXAS        +nb; etype=TETR
!
!       pyramid
        case(4); nel= NRPRISM+NRHEXAS+NRTETRA+nb; etype=PYRA
        endselect
!
!       loop through physical attributes
        do iphys=1,ELEMS(nel)%nrphysics
!
!          build up list of physical attributes supported by the edge node
!
!          select attribute
           phys = ELEMS(nel)%physics(iphys)
!
!          check if on the list
           call locate_char(phys,phys_vect,num, number)
!
!          add to the list if needed
           if (number == 0) then
              num=num+1
              phys_vect(num) = phys
           endif
!
!       loop through physical attributes
        enddo
        !
!       determine faces adjacent to the edge
        call locate(nod_new,ELEMS(nel)%nodes(nvert(etype)+1),nedge(etype), ie)
        call edge_to_faces(etype,ie, nofaces)
        !
        !  .......loop through the element physics attributes
        nvar=0
        do iphys=1,ELEMS(nel)%nrphysics
           phys = ELEMS(nel)%physics(iphys)
           call locate_char(phys,PHYSA,NR_PHYSA, iii)
           !
           !  .........loop through the variable components
           do ivar=1,NR_COMP(iii)
              nvar = nvar+1
              !
              !  .........decode the BC flags for the faces
              call decodg(ELEMS(nel)%bcond(nvar),10,nface(etype), ibc_elem)
              !
              !  .........loop through the faces adjacent to the edge and
              !           use element face BC flags to establish BC flag for the
              !           mid-edge node
              do i=1,2
                 ifc = nofaces(i)
                 call copyBCflag(2,ibc_elem(ifc), ibc_nod(nvar))
              enddo
           enddo
        enddo
        !
        mdle = nel
        call min_order(etype,2,ie,0,NODES(mdle)%order, nord)
     enddo
     !
     !  .....determine the case number
     call find_case(num,phys_vect, icase)
     !
     !  .....encode BC flags for the node into a single nickname
     call encod(ibc_nod,2,NRINDEX, nbcond)
     !
     subd = -1; iact = .true.
     call nodgen(MEDG,icase,nbcond,-nel,nord,subd,iact, nod)
     if (nod_new.ne.nod) then
        write(*,*) 'hp3gen: nod_new,nod = ',nod_new,nod
        stop 1
     endif
#if DEBUG_MODE
     if (iprint_edge.eq.1) then
        write(*,7002) nod
7002    format('hp3gen: HAVE GENERATED EDGE NODE ',i6)
     endif
#endif
  enddo
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,*) 'hp3gen: HAVE GENERATED EDGE NODES'
  endif
#endif
  !
  !----------------------------------------------------------------------
  !
  !  ...generate mid-face nodes for rectangles
  do nr=1,NRRECTA
     nod_new = NRELIS+NRPOINT+NRCURVE+nr; ibc_nod=0
     nord=MAXP*10+MAXP
     !
     nrbl=0
     do is=1,2
        if (RECTANGLES(nr)%BlockNo(is).ne.0) then
           nrbl=nrbl+1
           neigbl(nrbl) = RECTANGLES(nr)%BlockNo(is)
        endif
     enddo
     !
     !  .....identify the case number for the mid-face node.................
     num=0

     !  ...loop through neighboring blocks
     do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)

        !  ...determine element number
        select case(lab)
        case(1); nel=nb                         ; etype=PRIS
        case(2); nel=NRPRISM+nb                 ; etype=BRIC
        case(4); nel=NRPRISM+NRHEXAS+NRTETRA+nb ; etype=PYRA
        end select

        !  ...loop through neighbor's physical attributes
        do iphys=1,ELEMS(nel)%nrphysics
           phys = ELEMS(nel)%physics(iphys)

           !  ...check if on the list of attributes associated to the face node
           call locate_char(phys,phys_vect,num, number)

           !  ...build up a list of physical attributes for the face node
           if (number.eq.0) then
              num=num+1 ; phys_vect(num)=phys
           endif
        enddo

        !  ...determine face local number
        call locate(nod_new, &
             ELEMS(nel)%nodes(nvert(etype)+nedge(etype)+1), &
             nface(etype), ifc)

        !  ...loop through neighbor's physical attributes
        nvar=0
        do iphys=1,ELEMS(nel)%nrphysics
           phys = ELEMS(nel)%physics(iphys)
           call locate_char(phys,PHYSA,NR_PHYSA, iii)
           !
           !  .........loop through the variable components
           do ivar=1,NR_COMP(iii)
              nvar = nvar+1

           !  ...decode the BC flags for the faces
              call decodg(ELEMS(nel)%bcond(nvar),10,nface(etype), ibc_elem)

              !  ...copy face BC flag to GLOBAL list associated to the face node
              call copyBCflag(3,ibc_elem(ifc), ibc_nod(nvar))
           enddo
        enddo
        !
        mdle = nel
        call decodg(ELEMS(nel)%face_orient,8,nface(etype),nface_orient)
        call min_order(etype,3,ifc,nface_orient(ifc),NODES(mdle)%order, nord)

     !  ...end of loop through neighboring blocks
     enddo
     !
     call find_case(num,phys_vect, icase)
     !
     !  ...encode BC flags for the node into a single nickname
     call encod(ibc_nod,2,NRINDEX, nbcond)
     !
     subd = -1; iact = .true.
     call nodgen(MDLQ,icase,nbcond,-nel,nord,subd,iact, nod)
     if (nod_new.ne.nod) then
        write(*,*) 'hp3gen: nod_new,nod = ',nod_new,nod
        stop 1
     endif
#if DEBUG_MODE
     if (iprint_face.eq.1) then
        write(*,7003) nod
7003    format('hp3gen: HAVE GENERATED FACE NODE ',i6)
     endif
#endif
  enddo
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,*) 'hp3gen: HAVE GENERATED RECTANGULAR FACE NODES'
  endif
#endif
  !
  !----------------------------------------------------------------------
  !
  !  ...generate mid-face nodes for triangles
  do nt=1,NRTRIAN
     nod_new = NRELIS+NRPOINT+NRCURVE+NRRECTA+nt; ibc_nod=0
     nord=MAXP
     !
     nrbl=0
     do is=1,2
        if (TRIANGLES(nt)%BlockNo(is).ne.0) then
           nrbl=nrbl+1
           neigbl(nrbl) = TRIANGLES(nt)%BlockNo(is)
        endif
     enddo
     !
     !  .....identify the case number and BC flag for the mid-face node
     num=0
     !
     !  .....loop through neighboring blocks (at most two...)
     do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)
        !
        !  .......determine element number
        select case(lab)
           !
           !  .......prism
        case(1); nel= nb; etype=PRIS
           !
           !  .......tetrahedron
        case(3); nel= NRPRISM+NRHEXAS+nb; etype=TETR
           !
           !  .......pyramid
        case(4); nel= NRPRISM+NRHEXAS+NRTETRA+nb; etype=PYRA
        end select
        !
        ! loop over element physics attributes
        do iphys=1,ELEMS(nel)%nrphysics
           ! physics attribute label
           phys = ELEMS(nel)%physics(iphys)
           !
           !  .........check if on the list of the physics attributes for the node
           call locate_char(phys,phys_vect,num, number)
           if (number.eq.0) then
              num=num+1
              phys_vect(num) = phys
           endif
        enddo
        !
        !  .......determine adjacent face
        call locate(nod_new, &
             ELEMS(nel)%nodes(nvert(etype)+nedge(etype)+1),nface(etype), ifc)
        !
        !  .......loop through the element physics attributes
        nvar=0
        do iphys=1,ELEMS(nel)%nrphysics
           phys = ELEMS(nel)%physics(iphys)
           call locate_char(phys,PHYSA,NR_PHYSA, iii)
           !
           !  .........loop through the variable components
           do ivar=1,NR_COMP(iii)
              nvar = nvar+1
              !
              !  .........decode the BC flags for the faces
              call decodg(ELEMS(nel)%bcond(nvar),10,nface(etype), ibc_elem)
              call copyBCflag(3,ibc_elem(ifc), ibc_nod(nvar))
           enddo
        enddo
        !
        mdle = nel
        call decodg(ELEMS(nel)%face_orient,8,nface(etype), nface_orient)
        call min_order(etype,3,ifc,nface_orient(ifc),NODES(mdle)%order, nord)
     enddo
     !
     call find_case(num,phys_vect, icase)
     !
     !  .....encode BC flags for the node into a single nickname
     call encod(ibc_nod,2,NRINDEX, nbcond)
     !
     subd = -1; iact = .true.
     call nodgen(MDLT,icase,nbcond,-nel,nord,subd,iact, nod)
     !
     if (nod_new.ne.nod) then
        write(*,*) 'hp3gen: nod_new,nod = ',nod_new,nod
        stop 1
     endif
#if DEBUG_MODE
     if (iprint_face.eq.1) then
        write(*,7003) nod
     endif
#endif
  enddo
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,*) 'hp3gen: HAVE GENERATED TRIANGULAR FACE NODES'
  endif
#endif
  !
  deallocate(ibc_nod,phys_vect , stat=istat)
  if (istat.ne.SUCCESS) then
    call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  endif
  !
!----------------------------------------------------------------------
  !
  !  ...filling ELEM_ORDER array
  !
#if DEBUG_MODE
  if (iprint.ne.-1) then
     write(*,*)'CALLING update_ELEM_ORDER'
  endif
#endif
  call update_ELEM_ORDER
  !
  !  ...generate geometry and Dirichlet dof
  !
#if DEBUG_MODE
  if (iprint.ne.-1) then
     write(*,*)'CALLING update_gdof'
  endif
#endif
  !
  call update_gdof
#if DEBUG_MODE
  if (iprint.ne.-1) then
     write(*,2000)
2000 format('GEOMETRY DOFs HAVE BEEN UPDATED')
  endif
#endif
  !
  call update_Ddof
#if DEBUG_MODE
  if (iprint.ne.-1) then
     write(*,2001)
2001 format('DIRICHLET DOFs HAVE BEEN UPDATED')
  endif
#endif
  !
  !
end subroutine hp3gen
!
!----------------------------------------------------------------------
!   latest revision    - Sep 07
!
!   purpose            - auxiliary routine of hp3gen;
!                        it copies a BC flag from element to a node
!                        obeying hierarchy of Dirichlet BC flags
!
!   arguments
!     in:
!              Nflag   = 1 vertex node
!                      = 2 mid-edge node
!                      = 3 mid-face node
!              IBCelem - BC flag for an element face (coming from
!                        the element)
!     in/out:
!              IBCnod  - BC flag for a node of the element
!----------------------------------------------------------------------
!
subroutine copyBCflag(Nflag,IBCelem, IBCnod)
!
use data_structure3D
!
  implicit none
!
  integer, intent(in)    :: NFlag, IBCelem
  integer, intent(inout) :: IBCnod
!
  select case(Nflag)
!
!  ...vertex or edge node
  case(1,2)
!
!  ..copy only if a Dirichlet BC flag
     if (IBCelem.eq.1) IBCnod = IBCelem
!
!  ...face node
  case(3)
     select case(IBCnod)
!
!  .....zero BC flag, update if Dirichlet
     case(0)
        if (IBCelem.eq.1) IBCnod = IBCelem
!
!  .....non-zero flag, check compatibility
     case default
        if ((IBCelem.eq.1).and.(IBCnod.ne.IBCelem)) then
           write(*,7001) Nflag,IBCnod,IBCelem
7001       format(' copyBCflag: INCOMPATIBLE FACE FLAGS', &
                  ' Nflag, IBCnod,IBCelem = ',3i3)
           stop 1
        endif
     endselect
!
  endselect
!
!
end subroutine copyBCflag
