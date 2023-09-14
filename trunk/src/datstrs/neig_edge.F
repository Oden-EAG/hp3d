c----------------------------------------------------------------------
c
c   routine name       - neig_edge
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - find neigbors for an element edge
c
c   arguments :
c     in:
c         Medge        - a mid-edge node
c         Maxn         - dimension of the arrays below
c                        (max number of neighbors)
c     out:
c         Neig         - list of neighbors of the edge
c         Nedg_list    - local edge numbers
c         Norient_list - orientations of the mid-edge node
c         Nrneig       - number of adjacent elements
c
c----------------------------------------------------------------------
c
      subroutine neig_edge(Medge,Maxn,
     .                     Neig,Nedg_list,Norient_list,Nrneig)
c
      use GMP
      use data_structure3D
      use element_data
c
      implicit none
c
      integer, intent(in)  :: Medge
      integer, intent(in)  :: Maxn
      integer, intent(out) :: Neig(2,Maxn)
      integer, intent(out) :: Nedg_list(Maxn)
      integer, intent(out) :: Norient_list(Maxn)
      integer, intent(out) :: Nrneig
c
c  ...neighbors of a curve
      integer :: neigbl(100)
c
c  ...decoded orientations for edges
      integer :: nedge_orient(12)
c
c  ...element type
      integer :: etype
c
c  ...miscellanea
      integer :: nod,nc,nrbl,nb,lab,ib,ie,nel,ibeg
c
#if DEBUG_MODE
      integer :: iprint
      iprint = 0
#endif
c
c----------------------------------------------------------------------
c
c  ...initialize
      Neig = 0; Nedg_list = 0; Norient_list = 0; Nrneig = 0
c
      if (NODES(Medge)%ntype.ne.MEDG) then
        write(*,7001) Medge,S_Type(NODES(Medge)%ntype)
 7001   format('neig_edge: Medge,NODES(Medge)%type = ',i6,2x,a5)
        stop 1
      endif
c
c  ...a temporary version for the initial mesh only...
      if (NODES(Medge)%father.gt.0) then
        write(*,7002) Medge,NODES(Medge)%father
 7002   format('neig_edge: Medge,NODES(Medge)%father = ',i6,2x,i6)
        stop 1
      endif
c
c  ...the number of the GMP curve is (see hp3gen)
      nod = Medge
      nc = nod-(NRELIS+NRPOINT)
c
c  ...use GMP utility to find the adjacent GMP blocks
      call find_curve_to_block(nc,100, nrbl,neigbl)
      if (nrbl.gt.Maxn) then
        write(*,7003) Maxn,nrbl
 7003   format('neig_edge: Maxn,nrbl = ',2i4)
        stop 1
      endif
      Nrneig = nrbl
c
c  ...loop through neighboring blocks
      do ib=1,nrbl
        call decode(neigbl(ib), nb,lab)
c
c  .....determine element number
        select case(lab)
c
c  .....prism
        case(1); nel = nb
c
c  .....hexahedron
        case(2); nel = NRPRISM+nb
c
c  .....tetrahedron
        case(3); nel = NRPRISM+NRHEXAS+nb
c
c  .....pyramid
        case(4); nel = NRPRISM+NRHEXAS+NRTETRA+nb
        
        case default; write(*,*) 'neig_edge'; stop
        end select
c
c  .....store the neighbor, for the initial mesh
c       middle node numbers coincide with element numbers...
        Neig(1,ib) = nel; Neig(2,ib) = 0
        etype = ELEMS(nel)%etype
        ibeg = NVERT(etype)+1
        call locate(nod, ELEMS(nel)%nodes(ibeg),NEDGE(etype), ie)
        if (ie.eq.0) then
          write(*,7004)
 7004     format('neig_edge: INCONSISTENCY')
          stop 1
        endif
        Nedg_list(ib) = ie
        call decodg(ELEMS(nel)%edge_orient,2,12, nedge_orient)
        Norient_list(ib) = nedge_orient(ie)
c
c  ...end of loop through adjacent blocks
      enddo
c
#if DEBUG_MODE
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
c
      end subroutine neig_edge
