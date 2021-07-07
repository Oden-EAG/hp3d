!> Purpose - test routines determining neighbors
!> Last revision - July 21
!
      subroutine verify_neig
!
      use data_structure3D
      implicit none
      common /c_neig_edge/ iprint_neig_edge
!
!  ...workspace for elem_nodes
      integer, dimension(27) :: nodesl,norientl
!
!  ...workspace for neig_face
      integer, dimension(2)  :: nsid_list
!
!  ...workspace for neig-edge
      integer, parameter :: maxn=10
      integer :: neig(maxn), nedg_list(maxn), norient_list(maxn)
!
!  ...misc
      integer :: i,j, iface, mdle, nod, nrneig, loc, iprint, ie, iprint_neig_edge
      character(len=4) :: type
!
      iprint=0
      iprint_neig_edge=0
!
!  ...loop over active elements
      do i=1,NRELES
        mdle = ELEM_ORDER(i)
!
!  .....determine the element nodes and orientations
        call elem_nodes(mdle, nodesl,norientl)
!
!  .....loop over element's faces
        type=NODES(mdle)%type
        do iface=1,nface(type)
          j   = nvert(type) + nedge(type) + iface
          nod = nodesl(j)
!
          if (iprint.eq.1) then
            write(*,7010) mdle, iface, nod
 7010       format(' verify_neig: mdle, iface, nod    = ',3i10)
          endif
!
!  .......determine element neighbors for the face node
          call neig_face(nod, nrneig,neig,nsid_list,norient_list)
          if (iprint.eq.1) then
            write(*,7020) neig(1:nrneig)
 7020       format('verify_neig: mdle NODES NEIGHBORS = ',2i10)
          endif
!
!  .......look up the element on the list of the face neighbors
          call locate(mdle,neig,nrneig, loc)
!
!  .......check
          select case(loc)
!
!  .......mdle not found on list of neighbors
          case(0)
            write(*,*) nrneig, 'neig ', neig(1:nrneig)
            write(*,7030) mdle, iface, nod
 7030       format(' verify_neig: INCONSISTENCY, mdle,iface,nod = ',i6,i2,i6)
            call pause
!
!  .......check consistency of orientations
          case default
            if (norient_list(loc).ne.norientl(j)) then
              write(*,7040) mdle,iface,nod,norientl(j),norient_list(loc)
 7040         format(' verify_neig: INCONSISTENCY, ', &
                       'mdle,iface,nod,norient,norient_list = ', &
                       i6,i2,i6,2x,2i2)
              call pause
            endif
          endselect
!
!  .....end of loop through faces
        enddo
        if (iprint.eq.1) call pause
!
!  .....loop over element's edges
        do ie=1,nedge(type)
          j   = nvert(type) + ie
          nod = nodesl(j)
!
          if (iprint.eq.1) then
            write(*,7050) mdle, ie, nod
 7050       format(' verify_neig: mdle, ie, nod    = ',3i10)
          endif
!
!  .......determine element neighbors for the edge node
   30     call neig_edge(nod,maxn, nrneig,neig,nedg_list,norient_list)
!
!  .......look up the element on the list of the edge neighbors
          call locate(mdle,neig,nrneig, loc)
!
!  .......check
          select case(loc)
!
!  .......mdle not found on list of neighbors
          case(0)
            write(*,*) nrneig, 'neig ', neig(1:nrneig)
            write(*,7060) mdle, ie, nod
 7060       format(' verify_neig: INCONSISTENCY, mdle,ie,nod = ',i6,i2,i6)
            call pause
            iprint_neig_edge=1
            go to 30
!
!  .......check consistency of orientations
          case default
            if (norient_list(loc).ne.norientl(j)) then
              write(*,7070) mdle,NODES(mdle)%type,ie,nod,norientl(j),norient_list(loc)
 7070         format(' verify_neig: INCONSISTENCY, ', &
                       'mdle,type,ie,nod,norient,norient_list = ', &
                       i10,a4,i3,i10,2x,2i2)
              call pause
            endif
          endselect
!
!  .....end of loop through edges
        enddo
!
!  ...end of loop through elements
      enddo
!
!
      end subroutine verify_neig


