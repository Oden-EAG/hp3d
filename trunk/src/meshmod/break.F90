!-------------------------------------------------------------------------
!> Purpose : routine breaks an element middle node according to a given 
!!           refinement flag. Compatibility is assumed across the faces
!!           and the edges.
!!
!> @param[in] Mdle - middle node
!> @param[in] Kref - refinement flag
!!
!> rev@July 2019
!-------------------------------------------------------------------------
!
subroutine break(Mdle,Kref)
!
  use data_structure3D
  use element_data
  use refinements
  use par_mesh  , only: DISTRIBUTED
  use MPI_param , only: RANK
!  
  implicit none
!
  integer, intent(in) :: Mdle, Kref
!
  character(len=4)        :: type
  integer, dimension(27)  :: nodesl,norientl
  integer, dimension(8)   :: novert
  integer, dimension(6)   :: kreff
  integer, dimension(12)  :: krefe
  integer, dimension(4)   :: iv
  integer, dimension(4,6) :: neig
  integer :: i,j, is, iprint, iact, iface, ipass
  integer :: kref_face, kref_edge, nod, nr_vert, nrsons, subd
!
!-------------------------------------------------------------------------
!
  iprint = 0
!
! record on history file (needed for debugging)
  write(NHIST,8001) Mdle,Kref
8001 format(i6, 2x,i3,' # element number and refinement kind')
!
! nodal connectivities
  call elem_nodes(Mdle, nodesl,norientl)
!
! determine refinements for the element faces and edges
  type = NODES(Mdle)%type
  call find_face_ref_flags(type,Kref, kreff)
  call find_edge_ref_flags(type,Kref, krefe)
!
! printing
  if (iprint.eq.2) then
     call elem_show(Mdle, type, nodesl, norientl)
7000 format('break: Mdle = ',i6)
     write(*,7011) Kref
7011 format('       Kref = ',i2)
     write(*,7012) kreff(1:nface(type))
7012 format('       kreff = ',6i3)
     write(*,7013) krefe(1:nedge(type))
7013 format('       krefe = ',12i2)
  endif
!  
!
!=========================================================================
! E D G E S                                                              |
!=========================================================================
!
! loop over edges and generate INACTIVE son nodes
  do i=1,nedge(type)
     nod=nodesl(nvert(type)+i)
!    if edge is unrefined and needs to be refined
     if ( (krefe(i).ne.0).and.(is_leaf(nod)) ) then
        call edge_to_vert(type,i, iv(1),iv(2))
        nr_vert=2
!       collect edge vertex nodes        
        do j=1,2 ; novert(j) = nodesl(iv(j)) ; enddo
!       break edge
        iact = 0; kref_edge = 1
        call nodbreak(nod,kref_edge,iact,novert,nr_vert)
     endif
  enddo
!
!
!=========================================================================
! F A C E S                                                              |
!=========================================================================
!
! loop over faces and generate ACTIVE/INACTIVE son nodes
  do i=1,nface(type)
     iface=nvert(type)+nedge(type)+i
     nod=nodesl(iface)
!
!    face should be refined
     if (kreff(i).ne.0) then
!
!       collect vertices
        call face_to_vert(type,i, iv(1),iv(2),iv(3),iv(4))
        nr_vert=4
        do j=1,4 ; novert(j)=nodesl(iv(j)) ; enddo
!
!       modify refinement according to orientation
        call change_ref_flag('l2g',NODES(nod)%type,kreff(i), &
                             norientl(iface), kref_face)
!
!       check
        call check_ref(NODES(nod)%type,NODES(nod)%ref_kind, &
                            kref_face, ipass)
        if (ipass.eq.0) then
           write(*,7002) mdle,i,NODES(nod)%ref_kind,kref_face
7002       format('break: mdle,i,NODES(nod)%ref_kind,kref_face = ', &
                i6,2x,3(2x,i2))
           call result
           return
        endif
!
!       compare desired refinement "kref_face" to existing one "%ref_kind"
        kref_face=kref_face-NODES(nod)%ref_kind
!
!       by default, generate INACTIVE son nodes
        iact=0
!
!       existing refinement coincides with desired refinement
        if (kref_face.eq.0) then
!          do nothing, do not activate
!
!       existing refinement does NOT coincide with desired refinement
        else
           select case(NODES(nod)%type)
!           
!          triangular face           
           case('mdlt')
              nr_vert=0
              call nodbreak(nod,kref_face,iact,novert,nr_vert)
!              
!          quadrilateral face              
           case('mdlq')
              select case(NODES(nod)%ref_kind)
!             UNREFINED face : just break
              case(0)
                 call nodbreak(nod,kref_face,iact,novert,nr_vert)
!             ANISOTROPICALLY REFINED face : break son quads and son edge
              case(10,01)
!                loop over son quads and generate INACTIVE son nodes                 
                 do is=1,2
                    nr_vert=0
                    call nodbreak(NODES(nod)%sons(is),kref_face, &
                                  iact,novert,nr_vert)
                 enddo
!                break edge node                 
                 nr_vert=4
                 call nodbreak(NODES(nod)%sons(3),kref_edge,iact,novert,nr_vert)
              endselect
           endselect
        endif
     endif
  enddo
!
!
!=========================================================================
! M I D D L E                                                            |
!=========================================================================
!
  nr_vert=0
7100 format('break: Mdle ',i6,' ',a4, ' HAS BEEN BROKEN WITH Kref  ', i2)
  if ( (NODES(Mdle)%type.eq.'mdlb').and.(Kref.eq.111) ) then
     nr_vert=8
     do j=1,8 ; novert(j)=nodesl(j) ; enddo
  endif
!  
! break middle node
  iact = 0
  call nodbreak(Mdle,Kref,iact,novert,nr_vert)
!
  call activate_sons(Mdle)
!
  call deactivate(Mdle)
!
! update the number of active elements
  call nr_mdle_sons(type,Kref, nrsons)
  NRELES=NRELES+nrsons-1
!
! printing
  if (iprint.eq.1) then
     write(*,7100) Mdle, NODES(Mdle)%type, Kref
  endif
!
end subroutine break
!
!
!------------------------------------------------------------------------
!> Purpose : activate son nodes
!!
!> @param[in] Nod - node number
!------------------------------------------------------------------------
!
subroutine activate_sons(Nod)
  use data_structure3D
!
  implicit none
!
  integer, intent(in) :: Nod
  integer :: nrsons, i
!
!------------------------------------------------------------------------
!
  call find_nsons(Nod, nrsons)
! loop over sons
  do i=1,nrsons
     call activate(NODES(Nod)%sons(i))
  enddo
!
end subroutine activate_sons
