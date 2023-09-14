!> @brief Return nodal connectivity
!! @param[in]  Mdle     - middle node number
!! @param[out] Nodesl   - element nodes
!! @param[out] Norientl - their orientations
!> @date Feb 2023
!-----------------------------------------------------------------------
subroutine elem_nodes(Mdle, Nodesl,Norientl)
!
   use data_structure3D
   use refinements
   use constrained_nodes
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Nodesl(27),Norientl(27)
!
!..history information nfathl, isonl
   integer :: nfathl(MAXGEN),isonl(MAXGEN)
   integer :: nodesl_fath(27),norientl_fath(27)
   integer :: igen,nfath,nod,nson,nrgen,n_nodes,nrsons
!
#if DEBUG_MODE
   integer :: iprint
   iprint=0
#endif
!
!-----------------------------------------------------------------------
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) '------------------------------------------------'
      write(*,*) 'elem_nodes: Collecting ancestors FOR Mdle = ', Mdle
   endif
#endif
!
!..initialize
   Nodesl(1:27) = 0; Norientl(1:27) = 0
!
!-------------------------------------------------------------
!  Step 1 : Short cut for initial mesh
!-------------------------------------------------------------
   if (Is_root(Mdle)) then
      call elem_dump(Mdle, Nodesl,Norientl)
      return
   endif
!-------------------------------------------------------------
!  Step 2: Go up the tree to the initial mesh ancestor
!-------------------------------------------------------------
   nfath = NODES(Mdle)%father
   nson  = Mdle
   igen  = 0
   do while(nfath.gt.0)
      igen=igen+1
      nfathl(igen) = nfath
      call nr_mdle_sons(NODES(nfath)%ntype,NODES(nfath)%ref_kind, nrsons)
!     call locate(nson,NODES(nfath)%sons,nrsons, isonl(igen))
      isonl(igen) = nson - NODES(nfath)%first_son + 1
!     if (isonl(igen)<0 .or. isonl(igen)>nrsons) call pause
      nson = nfath
      nfath = NODES(nson)%father
   enddo
   nrgen = igen
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      do igen=1,nrgen
         write(*,7011) igen, nfathl(igen), isonl(igen)
7011     format('igen = ',i2,', father = ',i9,', son number = ',i2)
      enddo
   endif
#endif
!-------------------------------------------------------------
!  Step 3: Dump elem information
!-------------------------------------------------------------
!  nson is the top level element
   call elem_dump(nson, nodesl_fath,norientl_fath)
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      call elem_show(nson)
   endif
#endif
!-------------------------------------------------------------
! Step 4: Reconstruct
!-------------------------------------------------------------
   do igen=nrgen,1,-1
      call elem_nodes_one( &
            nfathl(igen),nodesl_fath,norientl_fath,isonl(igen), &
            nod,Nodesl,Norientl)
!
      if ((igen.eq.1).and.(INFO_CONSTRAINTS.eq.1)) then
         FATH_NODES  = nodesl_fath
         FATH_ORIENT = norientl_fath
         FATH_TYPE   = NODES(nfathl(1))%ntype
         SON_NUM     = isonl(1)
      endif
      nodesl_fath = Nodesl; norientl_fath = Norientl
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,*) 'igen = ', igen
         call elem_show(nod,NODES(nod)%ntype,Nodesl,Norientl)
      endif
#endif
!
  enddo
!
#if DEBUG_MODE
   if (iprint.eq.1) call pause
   if ((iprint.eq.2).and.(INFO_CONSTRAINTS.eq.1)) then
      write(*,7200) Mdle
7200  format('elem_nodes: FATHER INFO FOR Mdle = ',i6)
      n_nodes = NVERT(FATH_TYPE)+NEDGE(FATH_TYPE)+NFACE(FATH_TYPE)+1
      write(*,7201) FATH_NODES(1:n_nodes)
7201  format('FATHER NODES = ',27i6)
      write(*,7202) FATH_ORIENT(1:n_nodes)
7202  format('NODES ORIENT = ',27i6)
      write(*,7203) S_Type(FATH_TYPE),SON_NUM
7203  format('FATH_TYPE = ',a5,' SON_NUM = ',i1)
      call pause
   endif
#endif
!
end subroutine elem_nodes
