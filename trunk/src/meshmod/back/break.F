c---------------------------------------------------------------------
c
c   routine name       - break
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine breaks an element
c
c   arguments :
c
c     in:
c               Mdle   - element middle node number
c               Kref   - refinement flag
c
c-----------------------------------------------------------------------
c
      subroutine break(Mdle,Kref)
c
      use data_structure3D
      use element_data
      use refinements
      use GMP
#include "syscom.blk"
c
      character(len=4) :: type
c
c  ...element nodes and their orientations
      dimension nodesl(27),norientl(27)
c
c  ...vertices
      dimension novert(8)
c
c  ...refinement flags for element faces and edges
      dimension kreff(6), krefe(12)
c
c  ...vertex nodes numbers
      dimension iv(4)
c
c  ...neighbors of the element
      dimension neig(4,6)
c
c  ...volume of the element sons
      dimension vol_son(10)
c
c
      iprint = 1
c
c Kyungjoo
c  ...record the refinement in the history file
c      write(NHIST,8001) Mdle,Kref
c 8001 format(i6, 2x,i3,' element number and refinement kind')
c
c  ...determine the element nodes
      call elem_nodes(Mdle, nodesl,norientl)
c
c  ...determine the element neighbors
      call find_neig(Mdle, neig)
      if (iprint.eq.1) then
        write(*,6900) neig(1:4,1:nface(NODES(Mdle)%type))
 6900   format('break: neig = ',6(4i4,2x))
      endif
c
c  ...determine refinements for the element faces and edges
      type = NODES(Mdle)%type
      call find_face_ref_flags(type,Kref, kreff)
      call find_edge_ref_flags(type,Kref, krefe)
c
      if (iprint.eq.1) then
        write(*,7000) Mdle
        call elem_show(Mdle, NODES(Mdle)%type,nodesl,norientl)
 7000   format('break: Mdle = ',i6)
        write(*,7011) Kref
 7011   format('       Kref = ',i2)
        write(*,7012) kreff(1:nface(type))
 7012   format('       kreff = ',6i3)
        write(*,7013) krefe(1:nedge(type))
 7013   format('       krefe = ',12i2)
      endif
c
c  ...loop through element edges
      do ie=1,nedge(type)
        if (krefe(ie).ne.0) then
          call edge_to_vert(type,ie, iv(1),iv(2))
          nod = nodesl(nvert(type)+ie)
          if (NODES(nod)%ref_kind.eq.0) then
            nr_vert=2
            do j=1,2
              novert(j) = nodesl(iv(j))
            enddo
            call nodbreak(nod,1,0,novert,nr_vert)
          endif
        endif
      enddo
c
c  ...loop through element faces
      do if=1,nface(type)
        i=nvert(type)+nedge(type)+if
        if (kreff(if).ne.0) then
          call face_to_vert(type,if, iv(1),iv(2),iv(3),iv(4))
          nr_vert=4
          do j=1,4
            novert(j) = nodesl(iv(j))
          enddo
          nod = nodesl(nvert(type)+nedge(type)+if)
c
c  .......transform face refinement flag to the face coordinates
          call change_ref_flag('l2g',NODES(nod)%type,kreff(if),
     .                          norientl(i), krefm)
c
c  .......check compatibility of existing and planned refinements
          ipass=1
          call check_ref(NODES(nod)%type,NODES(nod)%ref_kind,
     .                        krefm, ipass)
          if (ipass.eq.0) then
            write(*,7002) mdle,if,NODES(nod)%ref_kind,krefm
 7002       format('break: mdle,if,NODES(nod)%ref_kind,krefm = ',
     .              i6,2x,3(2x,i2))
            call result; call print_GMP
            return
cccc            stop 1
          endif
c
c  .......refinement to do
          kref2 = krefm-NODES(nod)%ref_kind
          if (neig(1,if).eq.0) then
            iact = 1
          else
            iact = 0
          endif
          if (iprint.eq.1) then
            write(*,7008) if,kreff(if),krefm,
     .                    nod,NODES(nod)%ref_kind,kref2,iact
 7008       format('break: if,kreff(if),krefm = ',3i3,
     .             ' nod,NODES(nod)%ref_kind,kref2,iact = ',i5,3i3)
          endif
          if (kref2.eq.0) then
            call activate_sons(nod)
          else
            select case(NODES(nod)%type)
            case('mdlt')
              nr_vert = 0
              call nodbreak(nod,kref2,iact,novert,nr_vert)
            case('mdlq')
              select case(NODES(nod)%ref_kind)
              case(0)
                nr_vert = 0
                call nodbreak(nod,kref2,iact,novert,nr_vert)
              case(10,01)
                call activate_sons(nod)
                do is=1,2
                  nod_son = NODES(nod)%sons(is)
                  nr_vert = 0
                  call nodbreak(nod_son,kref2,0,novert,nr_vert)
                enddo
                nod_son = NODES(nod)%sons(3)
                nr_vert = 0
                call nodbreak(nod_son,1,0,novert,nr_vert)
              end select
            end select
          endif
        endif
      enddo
c
c  ...break the middle node
      if ((NODES(Mdle)%type.eq.'mdlb').and.(Kref.eq.111)) then
        nr_vert=8
        do j=1,8
          novert(j) = nodesl(j)
        enddo
      endif
      call nodbreak(Mdle,Kref,1,novert,nr_vert)
c
c  ...update the number of active elements
      call nr_mdle_sons(type,Kref, nrsons)
      NRELES = NRELES + nrsons-1
c
ccc      call element_volume(Mdle, vol_fath)
ccc      call nr_mdle_sons(NODES(Mdle)%type,NODES(Mdle)%ref_kind,
ccc     .     nrsons)
ccc      vol_sons = 0.d0
ccc      do ison=1, nrsons
ccc        call element_volume(NODES(Mdle)%sons(ison), vol_son(ison))
ccc        vol_sons = vol_sons + vol_son(ison)
ccc      enddo
c
ccc      if (iprint.eq.1) then
ccc        if (abs(vol_sons - vol_fath).gt.(GEOM_TOL)) then
ccc           write(*,*) 'break: VOLUME MISMATCH, mdle, vol_fath,vol_sones'
ccc     .          ,mdle,vol_fath, vol_sons
ccc           do ison=1, nrsons
ccc             write(*,*) ison, ' son= ', NODES(Mdle)%sons(ison),
ccc     .             'volume= ', vol_son(ison)
ccc           enddo
ccc        endif
ccc      endif
c
      if (iprint.eq.1) then
        write(*,7100) Mdle, Kref
 7100   format('break: Mdle ',i6,' HAS BEEN BROKEN WITH Kref  ', i2)
      endif
c
c
      end
c
c---------------------------------------------------------------------
c
c   routine name       - activate_sons
c
c---------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine activates sons of a node
c
c   arguments :
c     in:
c               Nod    - node number
c
c---------------------------------------------------------------------
c
      subroutine activate_sons(Nod)
c
      use data_structure3D
c
      call find_nsons(Nod, nrsons)
c
      do is=1,nrsons
        nod_son = NODES(Nod)%sons(is)
        call activate(nod_son)
      enddo
c
c
      end subroutine


