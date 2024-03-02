!----------------------------------------------------------------------
!
!   routine name       - neig_edge
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - find neigbors for an element edge
!
!   arguments :
!     in:
!         Medge        - a mid-edge node
!         Maxn         - dimension of the arrays below
!                        (max number of neighbors)
!     out:
!         Neig         - list of neighbors of the edge
!         Nedg_list    - local edge numbers
!         Norient_list - orientations of the mid-edge node
!         Nrneig       - number of adjacent elements
!
!----------------------------------------------------------------------
!
   subroutine neig_edge(Medge,Maxn, &
                           Neig,Nedg_list,Norient_list,Nrneig)
!
      use GMP
      use data_structure3D
      use element_data
!
      implicit none
!
      integer, intent(in)  :: Medge
      integer, intent(in)  :: Maxn
      integer, intent(out) :: Neig(2,Maxn)
      integer, intent(out) :: Nedg_list(Maxn)
      integer, intent(out) :: Norient_list(Maxn)
      integer, intent(out) :: Nrneig
!
!  ...neighbors of a curve
      integer :: neigbl(100)
!
!  ...decoded orientations for edges
      integer :: nedge_orient(12)
!
!  ...element type
      integer :: etype
!
!  ...miscellanea
      integer :: nod,nc,nrbl,nb,lab,ib,ie,nel,ibeg,iend
!
#if HP3D_DEBUG
      integer :: iprint
      iprint = 0
#endif
!
!----------------------------------------------------------------------
!
!  ...initialize
      Neig = 0; Nedg_list = 0; Norient_list = 0; Nrneig = 0
!
      if (NODES(Medge)%ntype.ne.MEDG) then
        write(*,7001) Medge,S_Type(NODES(Medge)%ntype)
 7001   format('neig_edge: Medge,NODES(Medge)%type = ',i6,2x,a5)
        stop 1
      endif
!
!  ...a temporary version for the initial mesh only...
      if (NODES(Medge)%father.gt.0) then
        write(*,7002) Medge,NODES(Medge)%father
 7002   format('neig_edge: Medge,NODES(Medge)%father = ',i6,2x,i6)
        stop 1
      endif
!
!  ...the number of the GMP curve is (see hp3gen)
      nod = Medge
      nc = nod-(NRELIS+NRPOINT)
!
!  ...use GMP utility to find the adjacent GMP blocks
      call find_curve_to_block(nc,100, nrbl,neigbl)
      if (nrbl.gt.Maxn) then
        write(*,7003) Maxn,nrbl
 7003   format('neig_edge: Maxn,nrbl = ',2i4)
        stop 1
      endif
      Nrneig = nrbl
!
!  ...loop through neighboring blocks
      do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)
!
!  .....determine element number
        select case(lab)
!
!  .....prism
        case(1); nel = nb
!
!  .....hexahedron
        case(2); nel = NRPRISM+nb
!
!  .....tetrahedron
        case(3); nel = NRPRISM+NRHEXAS+nb
!
!  .....pyramid
        case(4); nel = NRPRISM+NRHEXAS+NRTETRA+nb
        
        case default; write(*,*) 'neig_edge'; stop
        end select
!
!  .....store the neighbor, for the initial mesh
!       middle node numbers coincide with element numbers...
        Neig(1,ib) = nel; Neig(2,ib) = 0
        etype = ELEMS(nel)%etype
        ibeg = nvert(etype)+1
        iend = nvert(etype)+nedge(etype)
        call locate(nod,ELEMS(nel)%nodes(ibeg:iend),nedge(etype), ie)
        if (ie.eq.0) then
          write(*,7004)
 7004     format('neig_edge: INCONSISTENCY')
          stop 1
        endif
        Nedg_list(ib) = ie
        call decodg(ELEMS(nel)%edge_orient,2,12, nedge_orient)
        Norient_list(ib) = nedge_orient(ie)
!
!  ...end of loop through adjacent blocks
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7011) Medge
 7011   format('neig_edge: FIRST NEIGHBORS OF Medge = ',i5)
        write(*,7012) Neig(1,1:Nrneig)
 7012   format(10i8)
        write(*,7012) Nedg_list(1:Nrneig)
        write(*,7012) Norient_list(1:Nrneig)
        call pause
      endif
#endif
!
   end subroutine neig_edge
