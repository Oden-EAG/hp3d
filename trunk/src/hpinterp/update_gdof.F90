!-----------------------------------------------------------------------
!
!    routine name       - update_gdof
!
!-----------------------------------------------------------------------
!
!    latest revision    - Jul 18
!
!    purpose            - routine updates values of geometry degrees
!                         of freedom, using the GMP parametrizations
!
!    arguments          - none
!
!------------------------------------------------------------------------
!
subroutine update_gdof()
!
   use constrained_nodes
   use data_structure3D
   use element_data
   use environment, only : QUIET_MODE
   use GMP
!
   implicit none
!
   character(len=4) :: ntype
   character(len=6) :: mdltype
!
!..orientation of element nodes
   integer, dimension(12) :: nedge_orient
   integer, dimension(6)  :: nface_orient
!
!..element nodes
   integer, dimension(27)      :: nodesl, norientl
   integer, dimension(MAXNODM) :: nodm
!
!..order of approximation for an element
   integer, dimension(19) :: norder
!
!..reference coordinates for an element
   real*8, dimension(3,8) :: xsub
!
!..geometry dofs for an element
   real*8, dimension(3,MAXbrickH) :: xnod
!
!..auxiliary variables for timing
   real*8 :: start, OMP_get_wtime
!
!..auxiliary array with active mdle nodes
   integer, dimension(NRELES) :: mdlel
!
!..auxiliary variables
   integer :: iel, iv, ie, ifc, ind, iflag, i, k, loc
   integer :: mdle, nf, no, nod, nr_elem_nodes, nrnodm, nr_up_elem
!
!-----------------------------------------------------------------------
!
   start = OMP_get_wtime()
!
!..lower the GMP interface flag for all nodes
   do nod=1,NRNODS
      NODES(nod)%geom_interf=0
   enddo
!
!-----------------------------------------------------------------------  
!
!..fetch active elements
   mdle = 0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdlel(iel) = mdle
   enddo
!
!-----------------------------------------------------------------------
!
!..multiple passes
   do
!
!  ...initiate number of updated elements
      nr_up_elem=0
!
!  ...loop through active elements
      do 100 iel=1,NRELES
!
         mdle = mdlel(iel)
         call find_elem_type(mdle, mdltype)
!
!     ...skip if the element has already been processed
         if (NODES(mdle)%geom_interf.eq.1) cycle
         call refel(mdle, iflag,no,xsub)
!
!     ...determine nodes for the element (active and constrained)
!           and build the data base for the constrained nodes
!           (module constrained_nodes)
         call get_connect_info(mdle, nodesl,norientl)
!
!     ...use the info to determine all nodes of the corresponding
!           modified element
         call logic_nodes(mdle,nodesl, nodm,nrnodm)
!
!     ...check if all parent nodes have been updated
         ntype = NODES(mdle)%type
         nr_elem_nodes = nvert(ntype)+nedge(ntype)+nface(ntype)+1
         do i=1,nrnodm
            nod = nodm(i)
            call locate(nod,nodesl,nr_elem_nodes, loc)
!
!        ...if not a regular node of the element
            if (loc.eq.0) then
!
!           ...check if the node has been updated
               if (NODES(nod)%geom_interf.eq.0) go to 100
            endif
         enddo
!
!     ...update the number of processed elements
         nr_up_elem = nr_up_elem+1
!
!----------------------------------------------------------------------- 
!
!     ...loop through the element vertex nodes
         do iv=1,nvert(ntype)
            nod = nodesl(iv)
            if (NODES(nod)%geom_interf.eq.1) cycle
            call hpvert(iflag,no,xsub(1:3,iv), NODES(nod)%coord)
            NODES(nod)%geom_interf=1
!
!     ...end of loop through element vertices
         enddo
!
!-----------------------------------------------------------------------      
!
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
         call nodcor(mdle, xnod)
!
!     ...loop through element edge nodes
         do ie=1,nedge(ntype)
            ind = nvert(ntype)+ie
            nod = nodesl(ind)
!
!        ...if no gdof, mark as processed
            if (.not.associated(NODES(nod)%coord)) then
               NODES(nod)%geom_interf=1
            endif        
            if (NODES(nod)%geom_interf.eq.1) cycle
            if (NODES(nod)%act        .eq.0) cycle
!
            if (mdltype .ne. 'Linear') then
               call hpedge(mdle,iflag,no,xsub,ntype,             &
                           nedge_orient,nface_orient,norder,ie,  &
                           xnod,NODES(nod)%coord)
            endif
!
            NODES(nod)%geom_interf=1
!
!     ...end of loop through element edges
         enddo
!
!-----------------------------------------------------------------------      
!
         call nodcor(mdle, xnod)
!
!     ...loop through element face nodes
         do ifc=1,nface(ntype)
            ind = nvert(ntype)+nedge(ntype)+ifc
            nod = nodesl(ind)
!        ...if no gdof, mark as processed
            if (.not.associated(NODES(nod)%coord)) then
               NODES(nod)%geom_interf=1
            endif        
            if (NODES(nod)%geom_interf.eq.1) cycle
            if (NODES(nod)%act        .eq.0) cycle
!
            if (mdltype .ne. 'Linear') then
               call hpface(mdle,iflag,no,xsub,ntype,             &
                           nedge_orient,nface_orient,norder,ifc, &
                           xnod,NODES(nod)%coord)
            endif
            NODES(nod)%geom_interf=1
!
!     ...end of loop through element faces
         enddo
!
         NODES(mdle)%geom_interf=1
!
!  ...end of loop through elements
 100  continue
!
      if (nr_up_elem.eq.0) exit
!
!..end of infinite loop through multiple passes
   enddo
!
!-----------------------------------------------------------------------
!
!..update the middle node
!
!$OMP PARALLEL DO                                           &
!$OMP PRIVATE(iflag,mdle,mdltype,nedge_orient,nface_orient, &
!$OMP         no,norder,ntype,xnod,xsub)                    &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = mdlel(iel)
      if (.not.associated(NODES(mdle)%coord)) cycle
      call find_elem_type(mdle, mdltype)
      if (mdltype .ne. 'Linear') then
         ntype = NODES(mdle)%type
         call refel(mdle, iflag,no,xsub)
         call find_orient(mdle, nedge_orient,nface_orient)
         call find_order(mdle, norder)
         call nodcor(mdle, xnod)
         call hpmdle(mdle,iflag,no,xsub,ntype,          &
                     nedge_orient,nface_orient,norder,  &
                     xnod,NODES(mdle)%coord)
      endif
   enddo
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!
   if (.not. QUIET_MODE) then
      write(*,8004) OMP_get_wtime()-start
 8004 format(' update_gdof: ',f12.5,'  seconds')
   endif
!
end subroutine update_gdof
