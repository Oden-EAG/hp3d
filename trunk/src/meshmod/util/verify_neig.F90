!> @brief   test routines determining neighbors
!> @date    Feb 2023
subroutine verify_neig
  use data_structure3D
  implicit none
  integer, dimension(27) :: nodesl, norientl
  integer, dimension(2)  :: neig, nsid_list, norient_list
  integer :: i, j, iface, mdle, nod, nrneig, loc, ntype

#if DEBUG_MODE
  integer :: iprint
  iprint=0
#endif

  !  ...loop over active elements
  do i=1,NRELES
     mdle = ELEM_ORDER(i)
     call elem_nodes(mdle, nodesl,norientl)

     !  ...loop over element's faces
     ntype=NODES(mdle)%ntype
     do iface=1,nface(ntype)
        j   = nvert(ntype) + nedge(ntype) + iface
        nod = nodesl(j)

#if DEBUG_MODE
        if (iprint.eq.1) then
           write(*,7002) mdle, iface, nod
7002       format(' verify_neig: mdle, iface, nod    = ',3i10)
        endif
#endif

        call neig_face(nod, nrneig,neig,nsid_list,norient_list)

#if DEBUG_MODE
        if (iprint.eq.1) then
          write(*,7003) neig(1:nrneig)
7003      format('verify_neig: mdle NODES NEIGHBORS = ',2i10)
        endif
#endif
        call locate(mdle,neig,nrneig, loc)
        
        !  ...check
        select case(loc)
        !  ...mdle not found on list of neighbors
        case(0)
           write(*,*) nrneig, 'neig ', neig(1:nrneig)
           write(*,7000) mdle, iface, nod
7000       format(' verify_neig: INCONSISTENCY, mdle,iface,nod = ',i6,i2,i6)
           call pause
        !  ...face orientation not consistent
        case default
           if (norient_list(loc).ne.norientl(j)) then
              write(*,7001) mdle,iface,nod,norientl(j),norient_list(loc)
7001          format(' verify_neig: INCONSISTENCY, ', &
                   'mdle,iface,nod,norient,norient_list = ', &
                   i6,i2,i6,2x,2i2)
              call pause
           endif
        endselect

     enddo
#if DEBUG_MODE
     if (iprint.eq.1) call pause
#endif
  enddo

end subroutine verify_neig
