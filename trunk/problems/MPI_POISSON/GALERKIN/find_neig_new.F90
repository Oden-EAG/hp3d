
!--------------------------------------------------------------------------
!> Purpose : find element neighbors (4 per face) 
!            across faces of an ACTIVE element
!! @param[in]  Mdle - middle node
!! @param[out] Neig - neighbors
!              The four neighbors are returned in the FACE REFERENCE 
!
!              -------------             |\
!              |  4  |  3  |             |3\
!              -------------             |--
!              |  1  |  2  |             |\4|\
!              -------------             |1\|2\
!!                                       ------
!! @revision May 20
!--------------------------------------------------------------------------
!
      subroutine find_neig(Mdle, Neig_list)
!
      use element_data
      use data_structure3D
      use refinements
      implicit none
!
!  ...arguments
      integer,                 intent(in)  :: Mdle
      integer, dimension(4,6), intent(out) :: Neig_list
!
!  ...locals
      character(len=4) :: type
      integer, dimension(27) :: nodesl, norientl
      integer, parameter :: max_neig=40
      integer :: list_neig(max_neig)
      integer, dimension(27,max_neig)  :: nodesl_neig, norientl_neig
      integer, dimension(2)  :: neig, nsid_neig
      integer, dimension(4)  :: mface_sons
      integer :: nve, nrf, i, j, nod, nrneig, n, mdlen, loc, kref, nrsons, is,  &
                 iprint
!-------------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7010) Mdle
 7010   format('find_neig: Mdle = ',i6)
      endif
!
!  ...consistency check
      if (NODES(Mdle)%ref_kind.ne.0) then
        write(*,7020) Mdle
 7020   format('find_neig: Mdle = ',i6,' HAS BEEN REFINED') 
        stop 1
      endif
!
      type = NODES(Mdle)%type
      nve = nvert(type) + nedge(type); nrf = nface(type)
!
!  ...short cut for initial mesh elements only
      if (is_root(Mdle)) then
        do i=1,nrf
          Neig_list(1:4,i) = ELEMS(Mdle)%neig(i)
        enddo
        return
      endif
!
!  ...determine the element nodes
      call elem_nodes(Mdle, nodesl,norientl)
!
!  ...loop through element faces
      do i=1,nrf
        nod = nodesl(nve+i)   ! face node
!
!  .....determine element neighbors for the face
        call neig_face_extended(nod, nrneig,neig,nsid_neig,nodesl_neig,norientl_neig)
!    
        if (iprint.eq.1) then
          write(*,7030) i,nod,nrneig
 7030     format('find_neig: i = ',i1,' nod = ',i6' nrneig = ',i2)
          do j=1,nrneig
            nve = nvert(NODES(neig(j))%type) + nedge(NODES(neig(j))%type) 
            write(*,7040) neig(j),nsid_neig(j),nodesl_neig(nve+1:nve+nface(NODES(neig(j))%type),j)
 7040       format(' neig = ',i5,' nsid = ',i1,' nodesl = ',6i6)
          enddo
        endif 
        select case (nrneig)  ! nrneig = number of element neighbors for the face
        case(1)
!
!  .......element adjacent to the boundary; done
          Neig_list(1:4,i) = 0
          cycle
        case(2)
          if (Mdle.eq.neig(1)) then
!
!  .........copy the nodal connectivity info for the element neighbor to the first position
            neig(1) = neig(2)
            nodesl_neig(:,1) = nodesl_neig(:,2)
            norientl_neig(:,1) = norientl_neig(:,2)
          endif
        end select
!
!  .....initiate
        Neig_list(1:4,i) = neig(1)
!
!  .....the face node has not been refined; done
        if (NODES(nod)%ref_kind.eq.0) then
          cycle
        endif
!
!  .....determine mid-face node sons (and possibly grandsons)
        call find_mface_sons(nod, mface_sons)
        if (iprint.eq.1) then
          write(*,*) mface_sons(1:4)
 7050     format('find_neig: mface_sons = ',4i6)
        endif
!
!  .....initiate the neighbor list with the current neighbor
        n=1
        list_neig(1)  =  neig(1)
        do while (n.gt.0)
          mdlen = list_neig(1)
          type = NODES(mdlen)%type
          kref = NODES(mdlen)%ref_kind
!
!  .......if the element neighbor has been refined
          if (NODES(mdlen)%ref_kind.ne.0) then
            call nr_mdle_sons(type,kref, nrsons)
            if (iprint.eq.1) then
              write(*,*) 'mdlen, nrsons = ',mdlen, nrsons
            endif
!
!  .........loop through the sons of the element neigbor
            do is=1,nrsons
              n=n+1
              if (n.gt.max_neig) then
                write(*,*) 'find_neig: max_neig EXCEEDED !'
                stop 1
              endif
!
!  ...........determine and store the nodal connectivities for the son
              call elem_nodes_one(mdlen, &
                                  nodesl_neig(:,1), norientl_neig(:,1), is, &   
                                  list_neig(n), nodesl_neig(:,n), norientl_neig(:,n))
            enddo
          endif
!
!  .......look for any of the face node sons on the list of nodes of the element neighbor
          nve =nvert(type)+nedge(type)             ! combined number of vertices and edges
          nrf =nface(type)                         ! number of faces
          do j=1,4     
            call locate(mface_sons(j), nodesl_neig(nve+1:nve+nrf,1),nrf, loc)
            if (loc.ne.0) Neig_list(j,i) = mdlen
          enddo
          if (iprint.eq.1) then
            write(*,7060) n, Neig_list(1:4,i)
 7060       format('find_neig: n = ',i3,' Neig_list(1:4,i) = ',4i8)
            call pause
          endif
!
!  .......restack the lists
          n=n-1
          do j=1,n
            list_neig(j) = list_neig(j+1)
            nodesl_neig(:,j) = nodesl_neig(:,j+1)
            norientl_neig(:,j) = norientl_neig(:,j+1)
          enddo
!
!  .....end of loop through element neighbor list
        enddo
!
!  ...end of loop through element faces
      enddo
!
      end subroutine find_neig

!--------------------------------------------------------------------------
!> Purpose : find sons or grandsons of a mid-face element occupying 
!            the four positions on the face as indicated in routine find_neig
!! @param[in]  Mface - a mid-face node
!! @param[out] Mface_sons - its sons or grandsons
!!
!! @revision May 20
!--------------------------------------------------------------------------
!
      subroutine find_mface_sons(Mface, Mface_sons)
!
      use element_data
      use data_structure3D
      implicit none
!
!  ...arguments
      integer,               intent(in)  :: Mface
      integer, dimension(4), intent(out) :: Mface_sons
!
!  ...locals
      integer :: nson,i,j
!
!  ...initiate
      Mface_sons(1:4) = Mface
!
      select case(NODES(Mface)%type)
      case('mdlt')
        select case(NODES(Mface)%ref_kind)
        case(4)
          do j=1,4
            Mface_sons(j) = Son(Mface,j)
          enddo
        case default
          write(*,*) 'find_mface_sons: UNFINISHED 1'
          stop 1
        end select
      case('mdlq')
        select case(NODES(Mface)%ref_kind)
        case(11)
          do j=1,4
            Mface_sons(j) = son(Mface,j)
          enddo
        case(10)
          nson = Son(Mface,1)
          Mface_sons(1) = nson
          Mface_sons(4) = nson
          select case(NODES(nson)%ref_kind)
          case(0) ! nothing to do
          case(01)
            Mface_sons(1) = son(nson,1)
            Mface_sons(4) = son(nson,2)
          case default
            write(*,*) 'find_mface_sons: INCONSISTENCY 1'; stop 1
          end select
          nson = son(Mface,2)
          Mface_sons(2) = nson
          Mface_sons(3) = nson
          select case(NODES(nson)%ref_kind)
          case(0) ! nothing to do
          case(01)
            Mface_sons(2) = son(nson,1)
            Mface_sons(3) = son(nson,2)
          case default
            write(*,*) 'find_mface_sons: INCONSISTENCY 2'; stop 1
          end select
        case(01)
          nson = son(Mface,1)
          Mface_sons(1) = nson
          Mface_sons(2) = nson
          select case(NODES(nson)%ref_kind)
          case(0) ! nothing to do
          case(10)
            Mface_sons(1) = son(nson,1)
            Mface_sons(2) = son(nson,2)
          case default
            write(*,*) 'find_mface_sons: INCONSISTENCY 3'; stop 1
          end select
          nson = son(Mface,2)
          Mface_sons(3) = nson
          Mface_sons(4) = nson
          select case(NODES(nson)%ref_kind)
          case(0) ! nothing to do
          case(10)
            Mface_sons(4) = son(nson,1)
            Mface_sons(3) = son(nson,2)
          case default
            write(*,*) 'find_mface_sons: INCONSISTENCY 4'; stop 1
          end select
        case default
          write(*,*) 'find_mface_sons: UNFINISHED 2'
          stop 1
        end select
      end select
!
!  ...check consistency: all sons are supposed to be unrefined
      do i=1,4
        nson = Mface_sons(i)
        if (NODES(nson)%ref_kind.ne.0) then
          write(*,*) 'find_mface_sons: INCONSISTENVY 5'; stop 1
        endif
      enddo
!
      end subroutine find_mface_sons





