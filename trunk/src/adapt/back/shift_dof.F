c----------------------------------------------------------------------
c
c   routine name       - shift_dof
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 09
c
c   purpose            - routine copies dof within the data structure
c                        arrays shifting them to further location
c
c   arguments:
c     in:
c         Nrcomp       - number of loads
c         Nshift       - shift
c
c----------------------------------------------------------------------
c
      subroutine shift_dof(Nrcomp,Nshift)
c
      use data_structure3D
      use element_data
#include "syscom.blk"
c
c----------------------------------------------------------------------
c
c  ...check input
      if (Nrcomp.gt.MAXNRHS) then
        write(*,7001) Nrcomp,MAXNRHS
 7001   format('shift_dof: Nrcomp,MAXNRHS = ',2i3)
        stop 1
      endif
      if (Nrcomp+Nshift.gt.NRCOMS) then
        write(*,7002) Nrcomp,Nshift,NRCOMS
 7002   format('shift_dof: Nrcomp,Nshift,NRCOMS = ',3i3)
        stop 1
      endif
c
c  ...loop through nodes and shift dof
      do nod=1,MAXNODS
        if (NODES(nod)%type.eq.'none') cycle
        icase = NODES(nod)%case
        call ndof_nod(NODES(nod)%type,NODES(nod)%order,
     .                ndofH,ndofE,ndofV,ndofQ)
        if (Associated(NODES(nod)%zdofH)) then
          nvar = NREQNH(Icase)*Nrcomp
          newbeg = NREQNH(Icase)*Nshift+1
          newend = NREQNH(Icase)*(Nshift+Nrcomp)
          do k=1,ndofH
            NODES(nod)%zdofH(newbeg:newend,k)
     .      = NODES(nod)%zdofH(1:nvar,k)
          enddo
        endif
        if (Associated(NODES(nod)%zdofE)) then
          nvar = NREQNE(Icase)*Nrcomp
          newbeg = NREQNE(Icase)*Nshift+1
          newend = NREQNE(Icase)*(Nshift+Nrcomp)
          do k=1,ndofE
            NODES(nod)%zdofE(newbeg:newend,k)
     .      = NODES(nod)%zdofE(1:nvar,k)
          enddo
        endif
        if (Associated(NODES(nod)%zdofV)) then
          nvar = NREQNV(Icase)*Nrcomp
          newbeg = NREQNV(Icase)*Nshift+1
          newend = NREQNV(Icase)*(Nshift+Nrcomp)
          do k=1,ndofV
            NODES(nod)%zdofV(newbeg:newend,k)
     .      = NODES(nod)%zdofV(1:nvar,k)
          enddo
        endif
        if (Associated(NODES(nod)%zdofQ)) then
          nvar = NREQNQ(Icase)*Nrcomp
          newbeg = NREQNQ(Icase)*Nshift+1
          newend = NREQNQ(Icase)*(Nshift+Nrcomp)
          do k=1,ndofQ
            NODES(nod)%zdofQ(newbeg:newend,k)
     .      = NODES(nod)%zdofQ(1:nvar,k)
          enddo
        endif
      enddo
c
c
      end

