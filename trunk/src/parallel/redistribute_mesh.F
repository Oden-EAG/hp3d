c----------------------------------------------------------------------
c
c   routine name       - redistribute_mesh
c
c----------------------------------------------------------------------
c
c   latest revision    - Jun 07
c
c   purpose            - assigns elements to processors,
c                        deallocates dofs assigned to other processors
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine redistribute_mesh

      use parallelism
      use data_structure3D
      use element_data
      use control
      use parallelism
      use frsolmod , ONLY: OLD_ELEM_ORDER,NEW_ELEM_ORDER,ELEM_CENTER
#include "syscom.blk"
#include "mpif.h"
      common /cnodcor/ iprint_nodcor
c
c  ...nodes for a modified element and the corresponding number
c     of H1,H(curl),H(div) and L2 dof
      dimension nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),
     .                        ndofmV(MAXNODM),ndofmQ(MAXNODM)
c
c  ...number of variables for each physics attribute for an element
      dimension nrdofs(NR_PHYSA)
c
c  ...decoded index for a node
      dimension index(NRINDEX)
c
c  ...element geometry dof, direction vector
      dimension xnod(3,MAXbrickH),xc(3),direction(3)

      double precision, allocatable :: weights(:)
      double precision :: total_weight, total_weight_per_proc, dnrproc
c
      iprint=0
ccc      iprint_nodcor=1
c
c-----------------------------------------------------------------------
c
c  ...reorder the elements in order to minimize the bandwidth...
c     OLD/NEW ELEM_ORDER are defined in frsolmod
      direction(1:3) = 0.d0; direction(1) = 1.d0;
c
c  ...for safety deallocate first
      if (allocated(OLD_ELEM_ORDER)) deallocate (OLD_ELEM_ORDER)
      if (allocated(NEW_ELEM_ORDER)) deallocate (NEW_ELEM_ORDER)
      if (allocated(ELEM_CENTER))    deallocate (ELEM_CENTER)
c
c  ...then allocate
      allocate(OLD_ELEM_ORDER(NRELES))
      allocate(NEW_ELEM_ORDER(NRELES))
      allocate(ELEM_CENTER(NRELES))
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        OLD_ELEM_ORDER(iel)=mdle
        call nodcor(mdle, xnod)
        xc(1:3) = 0.d0
        do iv=1,nvert(NODES(mdle)%type)
          xc(1:3) = xc(1:3) + xnod(1:3,iv)
        enddo
        xc(1:3) = xc(1:3)/nvert(NODES(mdle)%type)
        call scalar_product(xc,direction, ELEM_CENTER(iel))
      enddo
      call sortm(NRELES,NEW_ELEM_ORDER,ELEM_CENTER)
      do iel=1,NRELES
        iel1=NEW_ELEM_ORDER(iel)
        NEW_ELEM_ORDER(iel) = OLD_ELEM_ORDER(iel1)
        write(*,*) 'iel,NEW_ELEM_ORDER(iel) = ',iel,NEW_ELEM_ORDER(iel)
ccc        if (NEW_ELEM_ORDER(iel).eq.1200) stop 1
ccc        if (iel/20*20.eq.iel) call pause
      enddo
      deallocate(OLD_ELEM_ORDER,ELEM_CENTER)

      write(*,*)'redistribute mesh: ELEMS sorted'
      call pause
c
c-----------------------------------------------------------------------

      allocate(ELEMS_PARALLEL_DATA(NRELIS),STAT=i1)
      if(i1/=0)then
        if(i1/=0)then
          write(*,*)'redistribute_mesh:',
     .      'ELEMS_PARALLEL_DATA not allocated!'
          stop
        endif
      endif

      if(NRPROC.eq.1)then
        do iel=1,NRELIS
          ELEMS_PARALLEL_DATA(iel)%processor_owner=MYRANK
        enddo
        return
      endif

c  ...simple load balancing
      nrelem1 = NRELIS
      nrproc1 = NRPROC
      nrelem_per_proc = nrelem1 / nrproc1
      ireminder=0
      if(nrelem_per_proc*nrproc1.ne.nrelem1)ireminder=1
      iel_counter=1
      iel = NEW_ELEM_ORDER(iel_counter)
      icount=0
      iproc_owner=0
      if(iprint.eq.1)then
        write(*,*)'redistribute_mesh:NRELIS,NRPROC,nrelem_per_proc',
     .    NRELIS,NRPROC,nrelem_per_proc
      endif
      do
        ELEMS_PARALLEL_DATA(iel)%processor_owner=iproc_owner
        icount=icount+1
        if(icount.ge.nrelem_per_proc+ireminder)then
          nrelem1 = nrelem1-icount
          nrproc1 = nrproc1-1
          if(nrproc1.gt.0)then
            nrelem_per_proc = nrelem1 / nrproc1
          else
            exit
          endif
          ireminder=0
          if(nrelem_per_proc*nrproc1.ne.nrelem1)ireminder=1
          icount = 0
          iproc_owner = iproc_owner+1
        endif
        iel_counter = iel_counter+1
        if(iel_counter.gt.NRELIS)exit
        iel = NEW_ELEM_ORDER(iel_counter)
      enddo

      if (iprint.eq.1) then
       if(MYRANK.eq.0)then
       do iel=1,NRELIS
         write(*,*)'ELEMS(',iel,')%processor_owner',
     .    ELEMS_PARALLEL_DATA(iel)%processor_owner
       enddo
       endif
      endif

      end subroutine redistribute_mesh

      subroutine clean_after_mesh_redistribution

      use parallelism

      deallocate(ELEMS_PARALLEL_DATA,STAT=i1)
      if(i1/=0)then
        write(*,*)'error deallocation ELEMS_PARALLEL_DATA',i1
        stop
      endif

      end subroutine clean_after_mesh_redistribution
