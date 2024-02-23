!----------------------------------------------------------------------
!
!   routine name       - result
!
!----------------------------------------------------------------------
!
!   latest revision    - June 2021
!
!   purpose            - (interactive) routine prints content
!                        of data structure arrays
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine result
!
   use data_structure3D
   use assembly
!
   implicit none
!
!..list of elements connected to a node and corresponding local node numbers
   integer, parameter :: max_conel = 100
   integer :: list_elem(max_conel),nod_loc(max_conel)
!
!..nodes for a modified element and the corresponding number
!  of H1,H(curl),H(div) and L2 dof
   integer :: nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM), &
                            ndofmV(MAXNODM),ndofmQ(MAXNODM)
   integer :: nrdofs(NR_PHYSA)
!
!..index for a node
   integer :: index(NRINDEX)
!
!..decoded boundary flag and case for a node
   integer :: ibcnd(NRINDEX_HEV), icase(NR_PHYSA)
!
!..egde and face orientations decode
   integer :: nedge_orient(12), nface_orient(6)
!..tetra elem print
   integer :: nodesl(27),norientl(27)
!
!..work space for neig_face
!   dimension neig(2),nsid_list(2),norient_list(2)
!   dimension neig_mdle(4,6)
!
   integer :: ntype
!
!..auxiliary
   integer :: ne,nb,i,ii,k,l,iload,nbeg,nend,nel,loc,nvar,ivar,ibegin,iend
   integer :: ndofH,ndofE,ndofV,ndofQ,nod,iel,mdle,nr_conel,idec,ndom
   integer :: nrdofm,nrdofc,nrnodm
!
   VTYPE :: zvoid1(1),zvoid2(1)
!
!----------------------------------------------------------------------
!
 10  continue
   write(*,*) 'result: SELECT:'
   write(*,*) 'EXIT....................................0'
   write(*,*) 'GENERAL DATA STRUCTURE PARAMETERS.......1'
   write(*,*) 'ELEMENT DATA............................2'
   write(*,*) 'NODE DATA...............................3'
   write(*,*) 'LIST ELEMENTS CONNECTED TO A NODE.......4'
   read(*,*) idec
!
   select case(idec)
!
      case(0)
         return
      case(1)
         write(*,6001) NRELIS,NRELES,NRNODS
 6001    format(' NRELIS,NRELES,NRNODS            = ',3i10)
         write(*,6002) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
 6002    format(' NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ = ',4i10)
         write(*,6003) MAXNODS,NPNODS
 6003    format(' MAXNODS,NPNODS                  = ',2i10)
         write(*,6004) NRCOMS,N_COMS
 6004    format(' NRCOMS,N_COMS                   = ',2i10)
!
      case(2)
         write(*,*) 'SET INITIAL MESH ELEMENT NUMBER'
         read(*,*) nel
         if (nel.gt.NRELIS) then
            write(*,*) 'THE ELEMENT DOES NOT EXIST'
            goto 10
         endif
         write(*,7001) nel
 7001    format(' ELEMENT = ',i10)
         write(*,7002) S_Type(ELEMS(nel)%etype)
 7002    format(' TYPE = ',a4)
         ii=0
         do i=1,ELEMS(nel)%nrphysics
            call locate_char(ELEMS(nel)%physics(i),PHYSA,NR_PHYSA, loc)
            nvar = NR_COMP(loc)
            do ivar=1,nvar
               ii=ii+1
               write(*,7003) ELEMS(nel)%physics(i),ivar,ELEMS(nel)%bcond(ii)
 7003          format(' PHYSICS ATTRIBUTE = ',a5,' ivar = ',i2, &
                      ' BC FLAGS = ',i6)
            enddo
         enddo
         nb=0 ; ne=nvert(ELEMS(nel)%etype)
         write(*,7004) ELEMS(nel)%nodes(nb+1:ne)
 7004    format(' VERTEX NODES = ',8i10)
         nb=ne; ne=nb+nedge(ELEMS(nel)%etype)
         write(*,7005) ELEMS(nel)%nodes(nb+1:ne)
 7005    format(' EDGE NODES = ',12i10)
         nb=ne; ne=nb+nface(ELEMS(nel)%etype)
         write(*,7006) ELEMS(nel)%nodes(nb+1:ne)
 7006    format(' FACE NODES = ',12i10)
         write(*,7007) ELEMS(nel)%nodes(ne+1)
 7007    format(' MIDDLE NODE = ',i10)
         call decodg(ELEMS(nel)%edge_orient,2,12, nedge_orient)
         write(*,7008) nedge_orient(1:nedge(ELEMS(nel)%etype))
 7008    format(' EDGE ORIENTATIONS DECODED = ',12i2)
         call decodg(ELEMS(nel)%face_orient,8,6, Nface_orient)
         write(*,7009) nface_orient(1:nface(ELEMS(nel)%etype))
 7009    format(' FACE ORIENTATIONS DECODED = ',6i2)
         write(*,7010) ELEMS(nel)%neig(1:nface(ELEMS(nel)%etype))
 7010    format(' FACE NEIGHBORS = ',6i10)
         write(*,7011) ELEMS(nel)%GMPblock
 7011    format(' GMP BLOCK = ',6i10)
!
      case(3)
         write(*,*) 'SET NODE NUMBER'
         read(*,*) nod
         write(*,7021) nod, NODES(nod)%act
 7021    format(' NODE = ',i10, ' active = ', l2)
         write(*,7022) S_Type(NODES(nod)%ntype)
 7022    format(' TYPE = ',a4)
         call decod(NODES(nod)%case,2,NR_PHYSA, icase)
         write(*,7023) icase(1:NR_PHYSA)
 7023    format(' CASE = ',10i1)
         write(*,7025) NODES(nod)%order
 7025    format(' ORDER = ',i3)
         call decod(NODES(nod)%bcond,2,NRINDEX_HEV, ibcnd)
         write(*,7026) ibcnd(1:NRINDEX_HEV)
 7026    format(' BC FLAG = ',30i1)
         write(*,7027) NODES(nod)%visit
 7027    format(' VISITATION FLAG = ',i5)
         write(*,7028) NODES(nod)%subd
 7028    format(' SUBDOMAIN  FLAG = ',i5)
         write(*,7029) NODES(nod)%ref_kind
 7029    format(' REFINEMENT KIND = ',i3)
         write(*,7040) NODES(nod)%father
 7040    format(' FATHER = ',i10)
         if (NODES(nod)%ref_kind.ne.0) then
            write(*,7041) NODES(nod)%first_son,NODES(nod)%nr_sons
 7041       format(' FIRST_SON = ',i10,', NR_SONS = ',i10)
         endif
!
         select case (NODES(nod)%ntype)
            case (MDLB,MDLP,MDLN,MDLD)
               call find_domain(nod,ndom)
               write(*,8000) ndom
 8000          format(' DOMAIN = ',i3)
               call elem_nodes(nod, nodesl,norientl)
               ntype = NODES(nod)%ntype
!
               ibegin = 1
               iend   = nvert(ntype)
               write(*,7100) nodesl  (ibegin: iend)
               write(*,7200) norientl(ibegin: iend)
!
               ibegin = iend + 1
               iend   = iend + nedge(ntype)
               write(*,7110) nodesl  (ibegin: iend)
               write(*,7210) norientl(ibegin: iend)
!
               ibegin = iend + 1
               iend   = iend + nface(ntype)
               write(*,7120) nodesl  (ibegin: iend)
               write(*,7220) norientl(ibegin: iend)
         end select
!
!        call neig_face(nod, nrneig,neig,nsid_list,norient_list)
!        write(*,7310) neig(1:2)
!
 7100    format(' ELEM_VERT = ',8i6)
 7110    format(' ELEM_EDGE = ',12i6)
 7120    format(' ELEM_FACE = ',6i6)
 7200    format(' VERT_ORIT = ',8i6)
 7210    format(' EDGE_ORIT = ',12i6)
 7220    format(' FACE_ORIT = ',6i6)
!7310    format(' FACE_NEIG = ',6i6)
!
         call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
!
!     ...Geometry dof
         if (associated(NODES(nod)%dof)) then
            if (associated(NODES(nod)%dof%coord)) then
               write(*,7031)
 7031          format(' COORDINATES = ')
               do k=1,ndofH
                  write(*,7032) NODES(nod)%dof%coord(1:3,k)
 7032             format(10x,3f15.8)
               enddo
            endif
         endif
!
!     ...H1 dof
         if (associated(NODES(nod)%dof)) then
            if (associated(NODES(nod)%dof%zdofH)) then
               write(*,7033)
 7033          format(' H1 DOF = ')
               nvar=NREQNH(NODES(nod)%case)
               do iload=1,NRRHS
                  write(*,7035) iload
 7035             format(' iload = ',i3)
                  nbeg =(iload-1)*nvar; nend=nbeg+nvar
                  do k=1,ndofH
                     write(*,7034) NODES(nod)%dof%zdofH(nbeg+1:nend,k,N_COMS)
#if C_MODE
 7034                format(10x,5(2e12.5,2x))
#else
 7034                format(10x,5e23.15)
#endif
                  enddo
               enddo
            endif
         endif
!
!     ...H(curl) dof
         if (associated(NODES(nod)%dof)) then
            if (associated(NODES(nod)%dof%zdofE)) then
               write(*,7036)
 7036          format(' H(curl) DOF = ')
               nvar=NREQNE(NODES(nod)%case)
               do iload=1,NRRHS
                  write(*,7035) iload
                  nbeg =(iload-1)*nvar; nend=nbeg+nvar
                  do k=1,ndofE
                     write(*,7034) NODES(nod)%dof%zdofE(nbeg+1:nend,k,N_COMS)
                  enddo
               enddo
            endif
         endif
!
!     ...H(div) dof
         if (associated(NODES(nod)%dof)) then
            if (associated(NODES(nod)%dof%zdofV)) then
               write(*,7037)
 7037          format(' H(div) DOF = ')
               nvar=NREQNV(NODES(nod)%case)
               do iload=1,NRRHS
                  write(*,7035) iload
                  nbeg =(iload-1)*nvar; nend=nbeg+nvar
                  do k=1,ndofV
                     write(*,7034) NODES(nod)%dof%zdofV(nbeg+1:nend,k,N_COMS)
                  enddo
               enddo
            endif
         endif
!
!     ...L2 dof
         if (associated(NODES(nod)%dof)) then
            if (associated(NODES(nod)%dof%zdofQ)) then
               write(*,7038)
 7038          format(' L2 DOF = ')
               nvar=NREQNQ(NODES(nod)%case)
               do iload=1,NRRHS
                  write(*,7035) iload
                  nbeg =(iload-1)*nvar; nend=nbeg+nvar
                  do k=1,ndofQ
                     write(*,7034) NODES(nod)%dof%zdofQ(nbeg+1:nend,k,N_COMS)
                  enddo
               enddo
            endif
         endif
!
!
      case(4)
         MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR &
                 + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
         allocate(NEXTRACT(MAXDOFM))
         allocate(IDBC(MAXDOFM))
         allocate(ZDOFD(MAXDOFM,NRRHS))
         write(*,*) 'SET NODE NUMBER'
         read(*,*) nod
         call get_index(nod, index)
         write(*,7050) nod,NODES(nod)%case,index
 7050    format(' nod,case,index = ',i10,2x,i2,2x,20i2)
         l=0
         do iel=1,NRELES
            mdle = ELEM_ORDER(iel)
            call celem(mdle,1,nrdofs,nrdofm,nrdofc,              &
                       nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,  &
                       zvoid1,zvoid2)
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
 7051       format(' result: INCREASE max_conel to ',i3)
         else
            nr_conel=l
            do l=1,nr_conel
               mdle = list_elem(l)
!
               write(*,7052) l,mdle,NODES(mdle)%case,nod_loc(l)
 7052          format(' result: l,mdle,case,loc = ',i3,2x,i10,2x,i2,2x,i3)
            enddo
         endif
         deallocate(NEXTRACT,IDBC,ZDOFD)
!
   end select
   goto 10
!
!
end subroutine result
