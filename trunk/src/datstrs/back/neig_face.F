c----------------------------------------------------------------------
c
c   routine name       - neig_face
c
c----------------------------------------------------------------------
c
c   latest revision    - Jul 10
c
c   purpose            - find middle node neigbors for a mid-face node
c                        NOTES: Every mid-face node number has at least
c                        one mdle neighbor that occupies the whole face.
c                        The second neighbor (if it exists) must occupy
c                        either the whole face as well, or the face
c                        corresponding to the parent mid-face node
c
c   arguments :
c     in:
c         Mfce0        - a mid-face node
c     out:
c         Nrneig       - number of neighbors (2 or 1, if the face is
c                        on the boundary of the domain)
c         Neig         - YOUNGEST mdle nodes neighbors for the two
c                        sides of the face;
c                          0 indicates no neighbor;
c                          negative mdle node number indicates a 'big'
c                            neighbor
c         Nsid_list    - local face numbers in neighbors' coordinates
c         Norient_list - orientations of the mid-face node wrt to the
c                        neighbors
c
c----------------------------------------------------------------------
c
      subroutine neig_face(Mfce0,
     .                     Nrneig,Neig,Nsid_list,Norient_list)
      use data_structure3D
      use element_data
      use refinements
#include "syscom.blk"
      common /cneig_face/ iprint
c
      dimension Neig(2),Nsid_list(2),Norient_list(2)
c
      character(len=4) :: type
c
      character(len=4) :: father_type, son_type,stype,ftype
      dimension nfather(maxgen),father_type(maxgen),
     .          nfather_ref_kind(maxgen),noson(maxgen)
c
c  ...decoded orientations for initial element faces
      dimension nface_ort(6)
c
c  ...face refinement flags
      dimension kreff(6)
c
c-----------------------------------------------------------------------
c
      select case(Mfce0)
      case(43)
        iprint=0
      case default
        iprint=0
      end select
      if (iprint.eq.1) then
        write(*,7000) Mfce0
 7000   format('neig_face: Mfce0 = ',i6)
      endif
c
c  ...check input
      select case(NODES(Mfce0)%type)
      case('mdlt','mdlq')
      case default
        write(*,7001) Mfce0,NODES(Mfce0)%type
 7001   format('neig_face: WRONG INPUT: Mfce, type = ',i6,2x,a5)
        stop 1
      end select
c
      Nrneig=0; Neig=0; Nsid_list=0; Norient_list=0
c
c-----------------------------------------------------------------------
c
c  ...Step 1: Go up the tree of nodes until you hit an initial mesh
c     mid-face node, or a mid-face node son of a middle node
c
      mfce = Mfce0
      nfath = NODES(mfce)%father
      igen=0
      do
        if (iprint.eq.1) then
          write(*,7002) igen,nfath
 7002     format('face_neig: igen,nfath = ',i3,i7)
        endif
c
c  .....initial mesh element mid-face node
        if (nfath.lt.0) then
c
c  .......initial mesh element number
          nel = -nfath
c
c  .......locate the face node on the element list of nodes and find
c         the face number
          type = ELEMS(nel)%type
          nn = nvert(type)+nedge(type); nrf = nface(type)
          call locate(mfce,ELEMS(nel)%nodes(nn+1:nn+nrf),nrf, if)
c
c  .......record the neighbor and the face number
          Nrneig = 1
          Neig(1) = nel     ! mdle=nel for initial mesh elements
          Nsid_list(1)=if
c
c  .......decode the element face orientations
          call decodg(ELEMS(nel)%face_orient,8,nrf, nface_ort)
          Norient_list(1) = nface_ort(if)
c
c  .......possible second neighbor
          neln = ELEMS(nel)%neig(if)
          if (neln.ne.0) then
            Nrneig = 2
            Neig(2) = neln
            type = ELEMS(neln)%type
            nn = nvert(type)+nedge(type); nrf = nface(type)
            call locate(mfce,ELEMS(neln)%nodes(nn+1:nn+nrf),nrf, ifn)
            Nsid_list(2)=ifn
            call decodg(ELEMS(neln)%face_orient,8,nrf, nface_ort)
            Norient_list(2) = nface_ort(ifn)
          endif
          exit
c
c  .....son of a middle node or a mid-face node
        else
          call nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind,
     .                      nr_msons)
          call nr_sons(NODES(nfath)%type,NODES(nfath)%ref_kind,
     .                       nrsons)
          call locate(mfce,NODES(nfath)%sons,nrsons, jss)
c
          ftype = NODES(nfath)%type
          select case(ftype)
c
c  .......son of of a middle node
          case('mdlp','mdln','mdld','mdlb')
            no_mdle = nvert(ftype)+nedge(ftype)+nface(ftype)+1
            call decode_ref(ftype,NODES(nfath)%ref_kind,
     .                      nref1,nref2,nref3)
c
c  .........loop through the middle sons only
            do is=1,nr_msons
              nson = NODES(nfath)%sons(is)
              stype = NODES(nson)%type
              jbeg = nvert(stype)+nedge(stype)
c
c  ...........loop through the face nodes only
              do j=jbeg+1,jbeg+nface(stype)
c
c  .............determine parent number and son number
                jp = npar_ref(ftype,j,is,nref1,nref2,nref3)
                js = nson_ref(ftype,j,is,nref1,nref2,nref3)
                if (iprint.eq.1) then
                  write(*,7003) nr_msons,nrsons,no_mdle,jss,is,j,jp,js
 7003             format('find_neig: nr_msons,nrsons,no_mdle,jss = ',
     .             4i3,/,'           is,j,jp,js = ',4i3)
                endif
c
c
c  .............if mdlp is the parent and son numbers match....
                if ((jp.eq.no_mdle).and.(js.eq.jss)) then
                  Nrneig=Nrneig+1
c
c  ...............brother is a neighbor
                  Neig(Nrneig) = NODES(nfath)%sons(is)
                  Nsid_list(Nrneig) = j-jbeg
                  Norient_list(Nrneig) =
     .            nort_ref(ftype,j,is,nref1,nref2,nref3)
                endif
              enddo
              if (iprint.eq.1) then
                if (is.eq.nr_msons) call pause
              endif
            enddo
            exit
c
c  .......son of a mid-face node
          case('mdlt','mdlq')
c
c  .........record the family information and pass to the father
            igen=igen+1
            father_type(igen) = NODES(nfath)%type
            nfather_ref_kind(igen) = NODES(nfath)%ref_kind
            call nr_face_sons(NODES(nfath)%type,NODES(nfath)%ref_kind,
     .                        nrsons)
            call locate(mfce,NODES(nfath)%sons,nrsons, noson(igen))
            mfce = nfath
            nfath = NODES(mfce)%father
          case default
            write(*,*) 'INCONSISTENCY', ftype
            write(*,8001) ftype, 1
 8001       format('neig_face: INCONSISTENCY ',a4,i2)
            stop 1
          end select
        endif
      enddo
c
c  ...save generation number for the mid-face node
      nrgen=igen
      if (iprint.eq.1) then
        do igen=nrgen,1,-1
          write(*,7010) igen,father_type(igen),
     .                  nfather_ref_kind(igen),noson(igen)
 7010     format('neig_face: igen = ',i2,' father_type = ',a5,
     .           ' ref_kind = ',i3,' noson = ',i3)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c  ...Step 2: Go down the tree of middle nodes, 'sticking to the face'
c
c  ...loop through the sides of the face
      do n=1,Nrneig
c
c  .....start with the neighbor for the mid-face ancestor found
c       in Step 1
        mdle = Neig(n)
        igen=nrgen
c
c  .....go down the tree of refinements of 'mdle'
        do
c
c  .......recall face number and its orientation
          if = Nsid_list(n); nort = Norient_list(n)
c
c  .......the flag will be raised if 'mdle' is replaced with one of its
c         sons
          nflag=0
          if (iprint.eq.1) then
            write(*,7012) n,igen,mdle,NODES(mdle)%ref_kind
 7012       format('neig_face: n,igen,mdle,NODES(mdle)%ref_kind = ',
     .              i2,i3,i6,i3)
            write(*,7013) if,nort
 7013       format('           if,nort = ',2i3)
          endif
c
c  .......the tree of middle node refinements has ended
          if (NODES(mdle)%ref_kind.eq.0) then
c
c  .........the tree of refinements of the mid-face node has ended,
c           'Mfce0' occupies a whole face of 'mdle'
            if (igen.eq.0) then
              Neig(n) =  mdle; go to 20
c
c  .........the tree of refinements of the mid-face node has NOT ended,
c           'mdle' is adjacent to an ancestor of the mid-face node and
c           'Mfce0' occupies only a portion of the element's face
            else
              Neig(n) = -mdle; go to 20
            endif
          endif
          call nr_mdle_sons(NODES(mdle)%type,NODES(mdle)%ref_kind,
     .                      nr_msons)
c
c  .......determine refinement flags for the element faces
          call find_face_ref_flags(NODES(mdle)%type,
     .                             NODES(mdle)%ref_kind, kreff)
c
          ftype = NODES(mdle)%type
          call decode_ref(ftype,NODES(mdle)%ref_kind,
     .                    nref1,nref2,nref3)
          nrve = nvert(ftype)+nedge(ftype)
c
c  .......loop through middle node sons
          do is=1,nr_msons
            nson = NODES(mdle)%sons(is)
            if (iprint.eq.1) then
              write(*,7004) is,nson
 7004         format('neig_face: is,nson = ',2i5)
            endif
            stype = NODES(nson)%type
c
c  .........loop through mid-face nodes of the son
            do if=1,nface(stype)
c
c  ...........local node number
              j=nvert(stype)+nedge(stype)+if
c
c  ...........determine the parent and son numbers
              jp = npar_ref(ftype,j,is,nref1,nref2,nref3)
              js = nson_ref(ftype,j,is,nref1,nref2,nref3)
              if (iprint.eq.1) then
                write(*,7005) if,j,jp,is0
 7005           format('neig_face: if,j,jp,is0 = ',4i3)
              endif
c
c  ...........the mid-face node is shared with the father; this is
c             a result of an anisotropic refinement of the element,
c             trade father for the son
              if (js.eq.0) then
                mdle = nson; nflag=1; exit
              endif
c
c  ...........quit if the tree of refinements for 'Mfce0' has ended
              if (igen.eq.0) cycle
c
c  ...........mid-face node numbers for elements must match
              if (jp-nrve.ne.Nsid_list(n)) cycle
c
c  ...........adjust for orientation
              new_ort = nort_ref(ftype,j,is,nref1,nref2,nref3)
              select case(type_nod(ftype,jp))
c
c  ...........triangular face
              case('mdlt')
                if (iprint.eq.1) then
                  write(*,7006) kreff(jp-nrve),nfather_ref_kind(igen),
     .                          nort,js,new_ort
 7006             format('neig_face: kreff(jp-nrve),',
     .            'nfather_ref_kind(igen),nort,js,new_ort = ',5i4)
                endif
                call rotate_trian(kreff(jp-nrve),nfather_ref_kind(igen),
ccc     .                            nort,js,iflag)
     .                            nort,js,new_ort)
                if (iprint.eq.1) then
                  write(*,7006) kreff(jp-nrve),nfather_ref_kind(igen),
     .                          nort,js,new_ort
                endif
c
c  .............face refinements must match
ccc                if (iflag.eq.100) cycle
                if (new_ort.eq.100) cycle
c
c  .............node son numbers must match
                if (js.ne.noson(igen)) cycle
c
c  .............switch to the son adjusting face number and orientation
                igen=igen-1; mdle = nson; nflag=1
                Nsid_list(n) = if
                Norient_list(n) = new_ort; exit
c
c  ...........quadrilateral face
              case('mdlq')
                call rotate_quad(kreff(if),nfather_ref_kind(igen),nort,
     .                           js,js1,nvoid)
c
c   ............two anistropic refinements of the face matching
c               isotropic refinement of the element...
                if (kreff(if)-nfather_ref_kind(igen).gt.0) then
c
c  ...............node sons must match
                  if (js.ne.noson(igen)) cycle
c
c  ...............two anisotropic face refinements needed....
                  if (igen.eq.1) then
                    Neig(n) = -mdle; go to 20
                  else
                    if (js1.ne.noson(igen-1)) cycle
                    igen=igen-2; mdle=NODES(mdle)%sons(is)
                    nflag=1;exit
                  endif
                else
c
c  ...............node sons must match
                  if (js.ne.noson(igen)) cycle
                  igen=igen-1; mdle = NODES(mdle)%sons(is)
                  nflag=1;exit
                endif
              end select
c
c  .........end of loop through mid-face nodes of the son
            enddo
            if (nflag.eq.1) exit
c
c  .......end of loop through middle node sons
          enddo
c
c  .......exit if nothing has changed
          if (nflag.eq.0) then
            Neig(n) = mdle; exit
          endif
c
c  .....end of loop through middle node generations
        enddo
 20     continue
        if (iprint.eq.1) then
          write(*,7020) n,igen,mdle,if,nort
 7020     format('neig_face: n,igen,mdle,if,nort = ',2i3,i6,2i3)
        endif
c
c  ...end of loop through the sides
      enddo
      if (iprint.eq.1) then
        do i=1,Nrneig
          write(*,7030) Neig(i),Nsid_list(i),Norient_list(i)
 7030     format('neig_face: Neig = ',i6,' FACE NO = ',i2,
     .           ' ORIENTATION = ',i1)
        enddo
        call pause
        call result
      endif
c
c
      end subroutine

