c----------------------------------------------------------------------
c
c   routine name       - result
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - (interactive) routine prints content
c                        of data structure arrays
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine result
c
      use data_structure3D
      use assembly
#include "syscom.blk"
c
c  ...list of elements connected to a node and corresponding local node
c     numbers
      parameter (max_conel=100)
      dimension list_elem(max_conel),nod_loc(max_conel)
c
c  ...nodes for a modified element and the corresponding number
c     of H1,H(curl),H(div) and L2 dof
      dimension nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),
     .                        ndofmV(MAXNODM),ndofmQ(MAXNODM)
      dimension nrdofs(NR_PHYSA)
c
c  ...index for a node
      dimension index(NRINDEX)
c
c  ...egde and face orientations decode
      dimension nedge_orient(12), nface_orient(6)
c  ...tetra elem print
      dimension nodesl(27),norientl(27)
c
c  ...work space for neig-face
      dimension neig(2), nsid_list(2),norient_list(2)
      dimension neig_mdle(4,6)
      character(len=4) :: type
c
c----------------------------------------------------------------------
c
c
 10   continue
      write(*,*) 'result: SELECT:'
      write(*,*) 'EXIT....................................0'
      write(*,*) 'GENERAL DATA STRUCTURE PARAMETRS........1'
      write(*,*) 'ELEMENT DATA............................2'
      write(*,*) 'NODE DATA...............................3'
      write(*,*) 'LIST ELEMENTS CONNECTED TO A NODE.......4'
      read(*,*) idec
c
      select case(idec)
c
      case(0)
        return
      case(1)
        write(*,6001) NRELIS,NRELES,NRNODS
 6001   format(' NRELIS,NRELES,NRNODS = ',3i10)
        write(*,6002) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
 6002   format(' NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ = ',4i10)
        write(*,6003) MAXNODS,NPNODS
 6003   format(' MAXNODS,NPNODS = ',2i10)
c
      case(2)
        write(*,*) 'SET INITIAL MESH ELEMENT NUMBER'
        read(*,*) nel
        if (nel.gt.NRELIS) then
          write(*,*) 'THE ELEMENT DOES NOT EXIST'
          goto 10
        endif
        write(*,7001) nel
 7001   format(' ELEMENT = ',i10)
        write(*,7002) ELEMS(nel)%type
 7002   format(' TYPE = ',a4)
        nphys=ELEMS(nel)%nrphysics
        write(*,7003) (ELEMS(nel)%physics(i),ELEMS(nel)%bcond(i),
     .                i=1,nphys)
 7003   format(' PHYSICS ATTRIBUTES WITH BC FLAGS = ',5(a5,2x,i6,2x))
        nb=0 ; ne=nvert(ELEMS(nel)%type)
        write(*,7004) ELEMS(nel)%nodes(nb+1:ne)
 7004   format(' VERTEX NODES = ',8i10)
        nb=ne; ne=nb+nedge(ELEMS(nel)%type)
        write(*,7005) ELEMS(nel)%nodes(nb+1:ne)
 7005   format(' EDGE NODES = ',12i10)
        nb=ne; ne=nb+nface(ELEMS(nel)%type)
        write(*,7006) ELEMS(nel)%nodes(nb+1:ne)
 7006   format(' FACE NODES = ',12i10)
        write(*,7007) ELEMS(nel)%nodes(ne+1)
 7007   format(' MIDDLE NODE = ',i10)
        call decodg(ELEMS(nel)%edge_orient,2,12, nedge_orient)
        write(*,7008) nedge_orient(1:nedge(ELEMS(nel)%type))
 7008   format(' EDGE ORIENTATIONS DECODED = ',12i2)
        call decodg(ELEMS(nel)%face_orient,8,6, Nface_orient)
        write(*,7009) nface_orient(1:nface(ELEMS(nel)%type))
 7009   format(' FACE ORIENTATIONS DECODED = ',6i2)
        write(*,7010) ELEMS(nel)%neig(1:nface(ELEMS(nel)%type))
 7010   format(' FACE NEIGHBORS = ',6i10)
        write(*,7011) ELEMS(nel)%GMPblock
 7011   format(' GMP BLOCK = ',6i10)
c
      case(3)
        write(*,*) 'SET NODE NUMBER'
        read(*,*) nod
        write(*,7021) nod, NODES(nod)%act
 7021   format(' NODE = ',i10, ' active = ', l2)
        write(*,7022) NODES(nod)%type
 7022   format(' TYPE = ',a4)
        write(*,7023) NODES(nod)%case
 7023   format(' CASE = ',i5)
        write(*,7024) NODES(nod)%index
 7024   format(' INDEX = ',i20)
        write(*,7025) NODES(nod)%order
 7025   format(' ORDER = ',i3)
        write(*,7026) NODES(nod)%bcond
 7026   format(' BC FLAG = ',i5)
        write(*,7426) NODES(nod)%visit
 7426   format(' VISITATION FLAG = ',i5)
        write(*,7027) NODES(nod)%ref_kind
 7027   format(' REFINEMENT KIND = ',i3)
        write(*,7028) NODES(nod)%father
 7028   format(' FATHER = ',i10)
        if (NODES(nod)%ref_kind.ne.0) then
           write(*,7029) NODES(nod)%first_son,NODES(nod)%nr_sons
 7029      format(' FIRST_SON = ',i10,', NR_SONS = ',i10)
        endif
c
        select case (NODES(nod)%type)
        case ('mdln','mdlp','mdlb','mdld')
           call find_domain(nod,ndom)
           write(*,8000) ndom
 8000      format(' DOMAIN = ',i3)
           call elem_nodes(nod, nodesl, norientl)
           type = NODES(nod)%type
c
           ibegin = 1
           iend   = nvert(type)
           write(*,7100) nodesl  (ibegin: iend)
           write(*,7200) norientl(ibegin: iend)
c
           ibegin = iend + 1
           iend   = iend + nedge(type)
           write(*,7110) nodesl  (ibegin: iend)
           write(*,7210) norientl(ibegin: iend)
c
           ibegin = iend + 1
           iend   = iend + nface(type)
           write(*,7120) nodesl  (ibegin: iend)
           write(*,7220) norientl(ibegin: iend)
        end select
c
ccc        call neig_face(nod, nrneig,neig,nsid_list,norient_list)
ccc        write(*,7320) neig(1:2)
c
 7100   format(' ELEM_VERT = ',8i6)
 7110   format(' ELEM_EDGE = ',12i6)
 7120   format(' ELEM_FACE = ',6i6)
 7200   format(' VERT_ORIT = ',8i6)
 7210   format(' EDGE_ORIT = ',12i6)
 7220   format(' FACE_ORIT = ',6i6)
 7310   format(' FACE_NEIG = ',6i6)
 7320   format(' FACE_NEIG = ',6i6)
c
        write(*,7030) NODES(nod)%geom_interf
 7030   format(' GMP INTERFACE FLAG = ',i2)
        call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
c
        if (associated(NODES(nod)%coord)) then
          write(*,7031)
 7031     format(' COORDINATES = ')
          do k=1,ndofH
            write(*,7032) NODES(nod)%coord(1:3,k)
 7032       format(10x,3f15.8)
          enddo
        endif
c
c  .....H1 dof
        if (associated(NODES(nod)%zdofH)) then
          write(*,7033)
 7033     format(' H1 DOF = ')
          nvar=NREQNH(NODES(nod)%case)
          do icomp=1,NRCOMS
            write(*,7035) icomp
 7035       format(' icomp = ',i3)
            nbeg =(icomp-1)*nvar; nend=nbeg+nvar
            do k=1,ndofH
              write(*,7034) NODES(nod)%zdofH(nbeg+1:nend,k)
#if C_MODE
 7034         format(10x,5(2e12.5,2x))
#else
 7034         format(10x,5e23.15)
#endif
            enddo
          enddo
        endif
c
c  .....H(curl) dof
        if (associated(NODES(nod)%zdofE)) then
          write(*,7036)
 7036     format(' H(curl) DOF = ')
          nvar=NREQNE(NODES(nod)%case)
          do icomp=1,NRCOMS
            write(*,7035) icomp
            nbeg =(icomp-1)*nvar; nend=nbeg+nvar
            do k=1,ndofE
              write(*,7034) NODES(nod)%zdofE(nbeg+1:nend,k)
            enddo
          enddo
        endif
c
c  .....H(div) dof
        if (associated(NODES(nod)%zdofV)) then
          write(*,7037)
 7037     format(' H(div) DOF = ')
          nvar=NREQNV(NODES(nod)%case)
          do icomp=1,NRCOMS
            write(*,7035) icomp
            nbeg =(icomp-1)*nvar; nend=nbeg+nvar
            do k=1,ndofV
              write(*,7034) NODES(nod)%zdofV(nbeg+1:nend,k)
            enddo
          enddo
        endif
c
c  .....L2 dof
        if (associated(NODES(nod)%zdofQ)) then
          write(*,7038)
 7038     format(' L2 DOF = ')
          nvar=NREQNQ(NODES(nod)%case)
          do icomp=1,NRCOMS
            write(*,7035) icomp
            nbeg =(icomp-1)*nvar; nend=nbeg+nvar
            do k=1,ndofQ
              write(*,7034) NODES(nod)%zdofQ(nbeg+1:nend,k)
            enddo
          enddo
        endif



c
      case(4)
        MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR
     .          + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
        allocate(NEXTRACT(MAXDOFM))
        allocate(IDBC(MAXDOFM))
        allocate(ZDOFD(MAXDOFM,NRCOMS))
        write(*,*) 'SET NODE NUMBER'
        read(*,*) nod
        call get_index(nod, index)
        write(*,7050) nod,NODES(nod)%case,index
 7050   format(' nod,case,index = ',i10,2x,i2,2x,20i2)
        l=0
        do iel=1,NRELES
          mdle = ELEM_ORDER(iel)
          call celem(mdle,1,
     .               nrdofs,nrdofm,nrdofc,
     .               nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .               zvoid,zvoid)
          call locate(nod,nodm,nrnodm, loc)
          if (loc.gt.0) then
            l=l+1
            if (l.le.max_conel) then
              list_elem(l) = mdle; nod_loc(l) = loc
            endif
          endif
        enddo
        if (l.gt.max_conel) then
          write(*,7051) l
 7051     format(' result: INCREASE max_conel to ',i3)
        else
          nr_conel=l
          do l=1,nr_conel
            mdle = list_elem(l)

            write(*,7052) l,mdle,NODES(mdle)%case,nod_loc(l)
 7052       format(' result: l, mdle,case,loc = ',
     .                      i3,2x,i10,2x,i2,2x,i3)
          enddo
        endif
        deallocate(NEXTRACT,IDBC,ZDOFD)
c
      endselect
      goto 10
c
c
      endsubroutine result
