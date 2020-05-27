c---------------------------------------------------------------------
c
c   routine name       - refresh
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine activates inactive nodes that need
c                        not be constrained
c
c   arguments          - none
c
c-----------------------------------------------------------------------
c
      subroutine refresh
c
      use data_structure3D
      use element_data
#include "syscom.blk"
c
      character(len=4) type
c
c  ...element nodes and their orientations
      dimension nodesl(27), norientl(27)
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7000)
 7000   format('refresh: DEBUGGING')
      endif
c
c  ...lower the visitation flags for all nodes
      call reset_visit
c
c  ...Step 1: loop through all active elements and flag their edge
c     and face nodes
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_nodes(mdle, nodesl, norientl)
        type = NODES(mdle)%type
        do j=1,nedge(type)+nface(type)
          i=nvert(type)+j
          NODES(nodesl(i))%visit = 1
        enddo
      enddo
c
c  ...Step 2: visit all inactive nodes and activate them if their
c     parent nodes are not flagged
      do nod=1,NRNODS
        if (NODES(nod)%act.eq.1) cycle
        nfath = NODES(nod)%father
        if (nfath.ne.0) then
          if (NODES(nfath)%visit.eq.0) then
            call activate(nod)
          endif
        endif
        if (iprint.eq.1) then
          write(*,7001) nod
 7001     format('refresh: ACTIVATING nod = ',i6)
        endif
      enddo
c
c  ...lower the visitation flags for the nodes
      call reset_visit
c
c
      end subroutine
c
c
c
c
c---------------------------------------------------------------------
c
c   routine name       - close
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - enforce 1-irregularity of the mesh
c
c-----------------------------------------------------------------------
c
      subroutine close
c
      use element_data
      use data_structure3D
#include "syscom.blk"
c
      character(len=4) type
c
c  ...list of elements to refine along with their refinement flags
      integer, allocatable ::  list(:,:)
c
c  ...element nodes and their orientation
      dimension nodesl(27),norientl(27)
c
c  ...refinement flags for element faces and edges
      dimension kreff(6),krefe(12)
c
      iprint=0
c
      allocate(list(2,NRELES))
c
      do
c
c  .....activate nodes that can be activated
        call refresh
c
c  .....flag edges and faces that must be refined
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          call get_connect_info(mdle, nodesl,norientl)
          call flag_constr_parents(mdle)
        enddo
c
c  .....identify elements to refine
        ic=0
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          call elem_nodes(mdle, nodesl,norientl)
          type = NODES(mdle)%type
          nflag=0
          krefe=0
          do j=1,nedge(type)
            i=nvert(type)+j
            nod = nodesl(i)
            if (iprint.eq.1) then
              write(*,7005) mdle,type,j,nod,NODES(nod)%visit
 7005         format('close: mdle,type,j,nod,NODES(nod)%visit = ',
     .               i6,2x,a5,2x,i2,i7,2x,i1)
            endif
            if (NODES(nod)%visit.eq.1) then
              krefe(j)=1; nflag=1
            endif
          enddo
          kreff=0
          do j=1,nface(type)
            i=nvert(type)+nedge(type)+1
            nod = nodesl(i)
            if (NODES(nod)%visit.eq.1) then
              select case(NODES(nod)%type)
              case('mdlt'); kreff(j)=1
              case('mdlq'); kreff(j)=11
              end select
              nflag=1
            endif
          enddo
          if (nflag.eq.1) then
            if (iprint.eq.1) then
              nre = nedge(NODES(mdle)%type)
              nrf = nface(NODES(mdle)%type)
              write(*,7003) mdle, krefe(1:nre)
 7003         format('close: mdle, krefe = ',i6,2x,12i2)
              write(*,7004) kreff(1:nrf)
 7004         format('             kreff = ',6i2)
              call pause
            endif
            ic=ic+1
            list(1,ic) = mdle
            call find_element_ref_flag2(type,kreff,krefe, kref)
            list(2,ic) = kref
          endif
        enddo
c
c  .....quit if there are no elements to refine
        if (ic.eq.0) exit
        if (iprint.ge.1) then
          write(*,7002) list(1:2,1:ic)
 7002     format('close: elements to refine = ',100(i6,':',i2,2x))
          call pause
        endif
c
c  .....refine the elements from the list
        do iel=1,ic
          mdle = list(1,iel); kref = list(2,iel)
          if (NODES(mdle)%ref_kind.eq.0) call refine(mdle,kref)
        enddo
      enddo
c
      call reset_visit
      deallocate(list)
c
c
      end subroutine
