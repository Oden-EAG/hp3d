c----------------------------------------------------------------------
c
c   routine name       - find_neig
c
c----------------------------------------------------------------------
c
c   latest revision    - Jul 10
c
c   purpose            - routine returns neighbors through faces
c                        for an element
c
c   arguments
c     in:
c            Mdle      - (middle node of an) element
c     out:
c
c            Neig      - neighbors
c
c----------------------------------------------------------------------
c
      subroutine find_neig(Mdle, Neig)
c
      use element_data
      use data_structure3D
#include "syscom.blk"
c
      dimension Neig(4,6)
c
      character(len=4) :: type
c
c  ...work space for neig_face
      dimension neig_loc(2),nsid_list_loc(2),norient_list_loc(2)
c
c  ...work space for elem_nodes
      dimension nodesl(27), norientl(27)
c
      Neig(1:4,1:6) = 0
c
      select case(Mdle)
      case(61)
        iprint=0
      case default
        iprint=0
      end select
      if (iprint.eq.1) then
        write(*,7000) Mdle
      endif
c
      type = NODES(Mdle)%type
c
c  Step 0: short cut for initial mesh elements only
c
      if (NODES(Mdle)%father.le.0) then
        nel = abs(NODES(Mdle)%father)
        do if=1,nface(type)
          Neig(1:4,if) = ELEMS(nel)%neig(if)
        enddo
        go to 999
      endif
c
c  Step 1: find neighbors for element mid-face nodes
c
      call elem_nodes(Mdle, nodesl,norientl)
      do if=1,nface(type)
        nod = nodesl(nvert(type)+nedge(type)+if)
        if (iprint.eq.1) then
          write(*,7020) if,nod
 7020     format('find_neig: if = ',i2,' MID_FACE NODE = ',i6)
        endif
        call neig_face(nod,
     .                 nrneig_loc,neig_loc,nsid_list_loc,norient_list)
        if (nrneig_loc.eq.2) then
          if (Mdle.eq.neig_loc(1)) then
            Neig(1:4,if) = neig_loc(2)
            if (iprint.eq.1) then
              write(*,7030) if, neig_loc(2)
 7030         format('find_neig: FOUND NEIGHBOR FOR FACE ',i2,i6)
            endif
          else
            Neig(1:4,if) = neig_loc(1)
            if (iprint.eq.1) then
              write(*,7030) if, neig_loc(1)
            endif
          endif
        endif
      enddo
c
 999  continue
c
      if (iprint.eq.1) then
        write(*,7000) Mdle
        write(*,7010) Neig(1:4,1:nface(type))
        call pause
 7000   format('find_neig: Mdle = ', i5)
 7010   format('           Neig = ',6(4i5,2x))
        write(*,*) 'find_neig: EXITING'
      endif
c
c
      end subroutine
