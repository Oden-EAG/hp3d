c-----------------------------------------------------------------------
c
c   routine name       - update_Ddof
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 11
c
c   purpose            - routine updates values of solution degrees
c                        of freedom for Dirichlet nodes
c
c-----------------------------------------------------------------------
c
      subroutine update_homogeneous_Ddof
c
      use data_structure3D
#include "syscom.blk"
c
c  ...element nodes
      dimension nodesl(27), norientl(27)
      character(len=4) :: type
c
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...lower the GMP interface flag for all nodes
      do nod=1,NRNODS
        NODES(nod)%geom_interf=0
      enddo
c
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_nodes(mdle, nodesl, norientl)
        type = NODES(mdle)%type
        nnod = nvert(type) + nedge(type) + nface(type) + 1
        do inod=1,nnod
          nod = nodesl(inod)
          if (NODES(nod)%geom_interf.eq.0) then
            if (associated(NODES(nod)%dof)) then
              if (associated(NODES(nod)%dof%zdofH)) then
                NODES(nod)%dof%zdofH = ZERO
              endif
              if (associated(NODES(nod)%dof%zdofE)) then
                NODES(nod)%dof%zdofE = ZERO
              endif
              if (associated(NODES(nod)%dof%zdofV)) then
                NODES(nod)%dof%zdofV = ZERO
              endif
              NODES(nod)%geom_interf=1
            endif
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
      if (iprint.ge.1) write(*,*) 'update_Ddof: DONE...'
c
      end
