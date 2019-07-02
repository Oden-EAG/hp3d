c----------------------------------------------------------------------
c
c   routine name       - refine
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - routine refines an element
c                        enforcing TWO mesh regularity rules:
c
c                        Rule 1: no element can be broken unless
c                                ALL its mid-face nodes are active
c                        Rule 2: an element refinement flag is
c                                always upgraded to accomodate
c                                existing refinements of faces
c
c   arguments :
c     in:
c             Mdle0    - middle node of an element, identified
c                        with the element
c             Kref0    - refinement flag as defined in the data
c                        structure
c                        hexas:   111, 110, 101, 011, 100, 010, 001
c                        prisms:  11, 10, 01
c                        tets:    11,12,13, 21,22,23,24,25,26
c                        pyramid: 10,01
c
c----------------------------------------------------------------------
c
      subroutine refine(Mdle0,Kref0)
c
      use element_data
      use data_structure3D
      use constrained_nodes
c
#include "syscom.blk"
      common /cneig_face/ iprint_face_neig
c
      character(len=4) :: type
c
c  ...list of elements to be refined, corresponding refinement flags
      parameter (maxel=20)
      dimension mdle_list(maxel), kref_list(maxel)
c
c  ...mid-face node middle node neighbors
      dimension neigf(2),nsid_list(2),norient_list(2)
c
c  ...an element nodes and their orientation
      dimension nodesl(27),norientl(27)
c
c  ...local face refinement flags
      dimension kreff(6)
c
      imod(j,mod) = j-(j-1)/mod*mod
c
c---------------------------------------------------------------------
c
      if (Mdle0.ge.2609) then
        iprint=0
      else
        iprint=0
      endif
      iprint_face_neig=0
c
    5 continue
      if (iprint.eq.1) then
        write(*,7002) Mdle0,Kref0
 7002   format('refine: DEBUGGING FOR Mdle0,Kref0 = ',i5,2x,i3)
      endif
c
c  ...put the element to be refined on the list
      mdle_list(1) = Mdle0
      kref_list(1) = Kref0
      n=1
c
      do while(n.gt.0)
 10     continue
c
c  .....pick the last element from the list
        mdle = mdle_list(n)
        kref = kref_list(n)
        type = NODES(mdle)%type
        if (iprint.ge.1) then
          write(*,*)'************************************************'
          write(*,7003) n,mdle,kref,type
 7003     format('refine: n,mdle,kref,type = ',i2,2x,i5,2x,i3,2x,a4)
          write(*,*)'************************************************'
        endif
c
        if (NODES(mdle)%ref_kind.ne.0) then
          write(*,7005) mdle,NODES(mdle)%ref_kind
 7005     format('refine: mdle,ref_kind = ',i5,2x,i3)
          stop 1
        endif
c
c  .....get the element nodes
        call get_connect_info(mdle, nodesl,norientl)
c
c  .....loop through element faces
        do if=1,nface(type)
          i = nvert(type)+nedge(type)+if
          nod = nodesl(i)
          select case(NODES(nod)%act)
c
c  .......the face is inactive (constrained); the element cannot be
c         refined, identify neighbor that has to be refined first
          case(0)
c
c  .........find the constraining face
            call decode2(NODES_CONSTR(i), nc,icase)
            nodp = NFACEC(nc)
c
c  .........determine neighbors of the mid-face node
 200        call neig_face(nod, nrneig,neigf,nsid_list,norient_list)
c
c  .........identify the neighbor across the face
            call locate(mdle,neigf,2, in)
            if (in.eq.0) then
              write(*,7006) mdle,if,nod,neigf
 7006         format('refine: ERROR !! mdle,if,nod, = ',3i6,
     .               ' neigf = ',2i6)
              stop 1
            endif
            inn = imod(in+1,2)
c
c  .........this should be a big neighbor indicated with a negative number
            if (neigf(inn).ge.0) then
              write(*,7006) mdle,if,nod,neigf
              stop 1
            endif
c
            mdlen = -neigf(inn)
            ifn = nsid_list(inn)
            norientn = norient_list(inn)
c
c  .........initiate face refinement flags for the neighbor with
c           the existing refinement of the face
            kreff=0
            call change_ref_flag('g2l',NODES(nodp)%type,
     .             NODES(nodp)%ref_kind,norientn, kreff(ifn))
c
c  .........find a minimum refinement for the neighbor to accomodate
c           the existing refinement of the face
            call find_element_ref_flag(NODES(mdlen)%type,0,kreff, krefm)
            if (krefm.lt.0) then
              write(*,*) 'NODES(mdlen)%type,kreff,krefm = ',
     .                    NODES(mdlen)%type,kreff,krefm
              write(*,*) 'nodp,NODES(nodp)%type,NODES(nodp)%ref_kind =',
     .                    nodp,NODES(nodp)%type,NODES(nodp)%ref_kind
              iprint_face_neig=1
              call pause
              go to 200
            endif
c
c  .........and the neighbor to the list of elements to be refined
            n=n+1
            mdle_list(n) = mdlen
            kref_list(n) = krefm
            go to 10
c
c  .......the face is active; record the existing refinement flag
          case(1)
            if (NODES(nod)%ref_kind.ne.0) then
              call change_ref_flag('g2l',NODES(nod)%type,
     .             NODES(nod)%ref_kind,norientl(i), kreff(if))
            else
              kreff(if)=0
            endif
          end select
        enddo
c
c  .....upgrade the element refinement flag to accomodate the existing
c       refinements of faces
        call find_element_ref_flag(type,kref,kreff, krefm)
        if (krefm.gt.0) then
c
c  .......proceed with the refinement
          if (iprint.ge.1) then
            write(*,7001) mdle,kref,krefm
 7001       format('refine: BREAKING mdle,kref,krefm = ',i5,2(2x,i3))
            call pause
          endif
          write(NHIST,*) mdle,krefm
          call break(mdle,krefm)
          n=n-1
c
c  .....the element cannot be refined, inconsistency
        else
          write(*,7010) mdle,type,kref,kreff(1:nface(type))
 7010     format('refine: mdle,type,kref,kreff = ',
     .            i5,2x,a4,2x,i2,2x,6(i2,2x))
          stop 1
        endif
      enddo
c
c
      end subroutine
