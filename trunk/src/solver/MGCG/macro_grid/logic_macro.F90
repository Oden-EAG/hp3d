!------------------------------------------------------------------------
!
!    routine name      - logic_macro
!
!------------------------------------------------------------------------
!
!    latest revision   - Jan 2018
!
!    purpose           - routine establishes dof connectivities
!                        between a macro-element and the corresponding
!                        coarse grid element
!
!   arguments :
!     in:
!                Igrid - grid index
!                MdleC - middle node of a coarse element
!            Nod_macro - macro-element nodes
!          Nrnod_macro - # of macro-element nodes
!         NrdofH_macro - # of H1      macro element dofs
!         NrdofE_macro - # of H(curl) macro element dofs
!         NrdofV_macro - # of H(div)  macro element dofs
!     out:
!               Nrnodm - # of modified coarse grid element nodes
!                        (excluding the middle node)
!                 Nodm - nodes of the modified coarse grid element
!                        (excluding the middle node)
!               NdofmH - the corresponding # of H1 single component dof
!               NdofmE - the corresponding # of H(curl) single component dof
!               NdofmV - the corresponding # of H(div)  single component dof
!          NdofH_macro - # of H1      dof for the macro-element nodes
!          NdofE_macro - # of H(curl) dof for the macro-element nodes
!          NdofV_macro - # of H(div)  dof for the macro-element nodes
!
!                        For the k-th dof of the MACRO-ELEMENT element
!      NrconH_macro(k) - number of connected coarse grid element H1 dof
!      NacH_macro(j,k) - list of the connected dof.   j=1,...,NrconH(k)
!   ConstrH_macro(j,k) - constraints coefficients.    j=1,...,NrconH(k)
!      NrconE_macro(k) - number of connected coarse grid element H(curl) dof
!      NacE_macro(j,k) - list of the connected dof.   j=1,...,NrconE(k)
!   ConstrE_macro(j,k) - constraints coefficients.    j=1,...,NrconE(k)
!      NrconV_macro(k) - number of connected coarse grid element H(div) dof
!      NacV_macro(j,k) - list of the connected dof.   j=1,...,NrconV(k)
!   ConstrV_macro(j,k) - constraints coefficients.    j=1,...,NrconV(k)
!
!-----------------------------------------------------------------------
!
   subroutine logic_macro(Igrid, MdleC, Nod_macro,Nrnod_macro,        &
                          NrdofH_macro, NrdofE_macro, NrdofV_macro,   &
                          Nrnodm,Nodm,NdofmH, NdofmE, NdofmV,         &
                          NdofH_macro, NdofE_macro, NdofV_macro,      &
                          NrconH_macro, NacH_macro, ConstrH_macro,    &
                          NrconE_macro, NacE_macro, ConstrE_macro,    &
                          NrconV_macro, NacV_macro, ConstrV_macro)
!
   use data_structure3D
   use refinements
   use prolongation
   use mg_data_structure
!
   implicit none
!
!..globals in
   integer, intent(in)  :: Igrid, MdleC, Nod_macro(Nrnod_macro), Nrnod_macro
!..total number of dof for a macro-element
   integer, intent(in)  :: NrdofH_macro, NrdofE_macro, NrdofV_macro
!
!..globals out
   integer, intent(out) :: Nrnodm
   integer, intent(out) :: Nodm(MAXNODM),   NdofmH(MAXNODM),      &
                           NdofmE(MAXNODM), NdofmV(MAXNODM)
   integer, intent(out) :: NdofH_macro(Nrnod_macro),              &
                           NdofE_macro(Nrnod_macro),   NdofV_macro(Nrnod_macro)
   integer, intent(out) :: NrconH_macro(NrdofH_macro), NacH_macro(NACDIM,NrdofH_macro),   &
                           NrconE_macro(NrdofE_macro), NacE_macro(NACDIM,NrdofE_macro),   &
                           NrconV_macro(NrdofV_macro), NacV_macro(NACDIM,NrdofV_macro)
   real*8,  intent(out) :: ConstrH_macro(NACDIM,NrdofH_macro),  &
                           ConstrE_macro(NACDIM,NrdofE_macro),  &
                           ConstrV_macro(NACDIM,NrdofV_macro)
!
!..locals
!..element characteristics
   character(len=4)     :: type, ftype, facetype
   integer              :: nrv, nre, nrf
!
!..workspace for elem_nodes
   integer              :: nodesl(27),norientl(27),norder(19)
!
!..offsets for first dof for coarse (local) element nodes
   integer              :: naHl(26), naEl(26), naVl(26)
!
!..offsets for first dof for coarse (coarse) element nodes
   integer              :: naH(MAXNODM),naE(MAXNODM),naV(MAXNODM)
!
!..offsets for first dof for macro-element element nodes
   integer              :: naH_macro(Nrnod_macro),naE_macro(Nrnod_macro),  &
                           naV_macro(Nrnod_macro)
!
!..workspace for logicC
   integer              :: idec
   integer              :: nrconH(MAXbrickH), nacH(NACDIM,MAXbrickH),   &
                           nrconE(MAXbrickE), nacE(NACDIM,MAXbrickE),   &
                           nrconV(MAXbrickV), nacV(NACDIM,MAXbrickV)
   real*8               :: ConstrH(NACDIM,MAXbrickH),                   &
                           ConstrE(NACDIM,MAXbrickE),                   &
                           ConstrV(NACDIM,MAXbrickV)
!
!..local offsets and dof counters
   integer              :: icH, icE, icV, ndofH, ndofE, ndofV, nvoid
!..local counters
   integer              :: i, j, nod, loc, loc1, locC, iv, nodl
   integer              :: kH, kE, kV, kHC, kEC, kVC
!
!..maps establishing injections between edge or face dof and element dof
   integer              :: mapH(MAXquadH),mapE(MAXquadE),mapV(MAXquadQ)
!
!..local for edges
   integer              :: ie, je, jp, kp, iedg(2), nodes_edge(3)
   integer              :: kH_macro, kE_macro, kV_macro
!
!..# of connected modified coarse grid element parent dof
   integer              :: ndof_macroH, ndof_macroE, ndof_macroV
   real*8               :: coeff
!
!..local for faces
   integer              :: if, nface_or, nrvf, nrfn, nod_face(9), nodesl_face(9)
!
!..workspace for face_to_vert_nos and face_to_vert_orient
   integer              :: nedge_or(4), nfedg_or(4) , nfver(4), nfedg(4)
!
   integer              :: nord, nordh, nordv, nordp(5), nordh_new
!
!..injection of a coarse grid node dof into dof after p-refinement
   integer              :: InjH(MAXquadH), InjE(MAXquadE), InjV(MAXquadV)
!
!..printing flag
   integer              :: iprint, inod
!
!----------------------------------------------------------------------------------------
!
   select case(MdleC)
   case(1)
      iprint = 0
   case default
      iprint = 0
   end select

!..initialization
!..modified coarse element
   Nrnodm = 0; Nodm = 0
   NdofmH = 0; NdofmE = 0; NdofmV = 0
!..macro-element
   NdofH_macro   = 0;    NdofE_macro   = 0;     NdofV_macro   = 0
!
!----------------------------------------------------------------------------------------
!
!..element type
   type = NODES(MdleC)%type
   nrv = nvert(type); nre = nedge(type); nrf = nface(type)

!
!..get element nodes and build the local database on constraints for the element
   call get_connect_infoC(Igrid,MdleC, nodesl,norientl)
!
!..get nodes of the modified coarse element and the corresponding
!  constraints coefficients
!..initialize arrays
   idec = 2
   call logicC(Igrid,MdleC,idec,                  &
               Nodm,NdofmH,NdofmE,NdofmV,Nrnodm,  &
               nrconH,nacH,constrH,               &
               nrconE,nacE,constrE,               &
               nrconV,nacV,constrV)
!
!..remove the middle node from the list
   Nodm(Nrnodm) = 0; Nrnodm = Nrnodm-1
!
!..establish offsets for the coarse element nodal dof
   icH=0; icE=0; icV=0
   do i = 1,nrv+nre+nrf
      naHl(i)=icH; naEl(i)=icE; naVl(i)=icV
      nod = nodesl(i)
      call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                    ndofH,ndofE,ndofV,nvoid)
      icH=icH+ndofH; icE=icE+NdofE; icV=icV+NdofV
   enddo
!
!..establish offsets for the modified coarse element nodal dof
   icH=0; icE=0; icV=0
   do i = 1,Nrnodm
      naH(i)=icH; naE(i)=icE; naV(i)=icV
      nod = Nodm(i)
      call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                    NdofmH(i),NdofmE(i),NdofmV(i),nvoid)
      icH=icH+NdofmH(i); icE=icE+NdofmE(i); icV=icV+NdofmV(i)
   enddo
!
!..establish offsets for the macro-element nodal dof
   icH=0; icE=0; icV=0
   do i=1,Nrnod_macro
      naH_macro(i)=icH; naE_macro(i)=icE; naV_macro(i)=icV
      nod = Nod_macro(i)
      call ndof_nod(NODES(nod)%type,NODES(nod)%order,  &
                    NdofH_macro(i),NdofE_macro(i),NdofV_macro(i),nvoid)
      icH=icH+NdofH_macro(i); icE=icE+NdofE_macro(i); icV=icV+NdofV_macro(i)
   enddo
   ! nrdofH_macro = icH; nrdofE_macro = icE; nrdofV_macro = icV
!
!
   if (iprint.eq.1) then
 7000 format(' i   = ',i3,' nod = ',i6,' TYPE = ',a6,' ORDER = ',i3, &
             ' naH = ',i3,' NdofH = ',i3, &
             ' naE = ',i3,' NdofE = ',i3, &
             ' naV = ',i3,' NdofV = ',i3)
!
      write(*,7001) MdleC
 7001 format('logic_macro: COARSE LOCAL ELEMENT NODES FOR MdleC = ',i6)
      do i = 1,nrv+nre+nrf
         nod = nodesl(i)
         call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                       ndofH,ndofE,ndofV,nvoid)
         write(*,7000) i,nod,NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                       naHl(i),ndofH, &
                       naEl(i),ndofE, &
                       naVl(i),ndofV
      enddo

      write(*,7002) MdleC
 7002 format('logic_macro: COARSE MODIFIED ELEMENT NODES FOR MdleC = ',i6)
      do i = 1,Nrnodm
         nod = Nodm(i)
         write(*,7000) i,nod,NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                       naH(i),NdofmH(i), &
                       naE(i),NdofmE(i), &
                       naV(i),NdofmV(i)
      enddo

      write(*,7003) MdleC
 7003 format('logic_macro: MACRO-ELEMENT NODES FOR MdleC = ',i6)
      do i=1,Nrnod_macro
         nod = Nod_macro(i)
         write(*,7000) i,nod,NODES(nod)%type,NODES(nod)%order, &
                       naH_macro(i),NdofH_macro(i), &
                       naE_macro(i),NdofE_macro(i), &
                       naV_macro(i),NdofV_macro(i)
      enddo
      call pause
   endif
!
!---------------------------------------------------------------------------------------
!
!..initiate the connectivities arrays
   NrconH_macro = 0; NacH_macro = 0; ConstrH_macro = 0.d0
   NrconE_macro = 0; NacE_macro = 0; ConstrE_macro = 0.d0
   NrconV_macro = 0; NacV_macro = 0; ConstrV_macro = 0.d0
!
!---------------------------------------------------------------------------------------
!..Step 0: Modified coarse grid nodes
!
   do i=1,Nrnodm
      nod = Nodm(i)
      call locate (nod, Nod_macro,Nrnod_macro, loc)
!
!  ...if the node belongs also to the macro-element then connect it to itself
      if (loc.ne.0) then
         kH = naH_macro(loc); kE = naE_macro(loc); kV = naV_macro(loc)
         call inject_dof(NODES(nod)%type, NODES_MG(nod)%orderC(Igrid),NODES(nod)%order,   &
                         InjH, InjE, InjV)
         do j=1,NdofmH(i)
            NrconH_macro(kH+injH(j))=1
            NacH_macro(1,kH+injH(j)) = naH(i)+j
            ConstrH_macro(1,kH+injH(j)) = 1.d0
         enddo
         do j=1,NdofmE(i)
            NrconE_macro(kE+injE(j)) = 1
            NacE_macro(1,kE+injE(j)) = naE(i)+j
            ConstrE_macro(1,kE+injE(j)) = 1.d0
         enddo
         do j=1,NdofmV(i)
            NrconV_macro(kV+injV(j)) = 1
            NacV_macro(1,kV+injV(j)) = naV(i)+j
            ConstrV_macro(1,kV+injV(j)) = 1.d0
         enddo
      endif
   enddo
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
!..Step 1: Vertices
!---------------------------------------------------------------------------------------
!
!..loop through the coarse grid element vertices (they may also belong
!  to macro-element)
   do iv=1,nrv
      nod = nodesl(iv)
!
!  ...skip if not on the list of macro-element nodes
      call locate(nod,Nod_macro,Nrnod_macro, loc)
      if (loc.eq.0) cycle
!
!  ...skip if on the list of modified coarse element node
      call locate(nod,Nodm,Nrnodm, locC)
      if (locC.ne.0) cycle
!
!  ...the vertex is a constrained vertex of the coarse grid element,
!     copy prolongation coefficients from those delivered by logicC
      kH = naH_macro(loc)+1  !  dof # within the macro-element
      kHC = naH(iv)+1        !  dof # within the coarse grid element
      NrconH_macro(kH) = NrconH(kHC)
      NacH_macro   (1:NrconH(kHC),kH) = nacH   (1:NrconH(kHC),kHC)
      ConstrH_macro(1:NrconH(kHC),kH) = constrH(1:NrconH(kHC),kHC)
!
!..end of loop through coarse grid element vertices
   enddo
!
!---------------------------------------------------------------------------------------
!..Step 2: Edges
!---------------------------------------------------------------------------------------
!
!..loop through the macro-element (coarse grid element) edges
   do ie=1,nre
!
!  ...local mid-edge node number
      i = nrv + ie
!
!  ...mid-edge node
      nodes_edge(3) = nodesl(i)

!  ...vertex nodes (account for edge orientation)
      call edge_to_vert(type,ie, iedg(1),iedg(2))
!  ...establish the injections between edge and element dof
      if (norientl(i).eq. 1) call swap(iedg(1),iedg(2))
      nodes_edge(1:2) = nodesl(iedg(1:2))
!
!  ...establish the injections between edge and element dof
      mapH(1:2) = naHl(iedg(1:2))+1
      nod = nodes_edge(3)
      call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                    ndofH,ndofE,nvoid,nvoid)
      do kH=1,ndofH
        mapH(2+kH) = naHl(i)+kH
      enddo
      do kE=1,ndofE
        mapE(kE) = naEl(i)+kE
      enddo
!
      call build_edge_tree(nodes_edge,NODES_MG(nod)%orderC(Igrid))
!
!  ...loop through local nodes in the tree
      do nodl=1,NR_NODES_PR
!
!     ...global node number
         nod = NODES_PR(nodl)%nod
!
!     ...skip if not on the list of macro-element nodes
         call locate(nod,Nod_macro,Nrnod_macro, loc)
         if (loc.eq.0) cycle
!
!     ...skip if on the list of modified coarse element node
         call locate(nod,Nodm,Nrnodm, locC)
         if (locC.ne.0) cycle
!
!     ...skip if the node is not in the interior of the edge
!        (i.e. it has been considered in the loop through vertices)
         if (.not. Edge_interior_node(nodl)) cycle
!
         call ndof_nod(NODES(nod)%type,NODES_PR(nodl)%order, &
                        ndofH,ndofE,nvoid,nvoid)
!
!     ...loop through nodal H1 dof
         do j=1,ndofH
!
!        ...dof number in the macro-element
            kH_macro = naH_macro(loc)+j
!
!        ...initiate # of connected modified coarse grid element parent dof
            ndof_macroH = 0
!
!        ...loop through connected parent edge dof
            do je=1,MDOFH
!
!           ...find the corresponding local element dof number
               kH = mapH(je)
!
!           ...loop through the modified coarse grid element parent dof
               do jp=1,nrconH(kH)
                  kp = nacH(jp,kH)
                  coeff = constrH(jp,kH)*NODES_PR(nodl)%coeffH(je,j)
!
                  if (coeff .eq. 0.d0) cycle

!              ...adjust # of connected parent dof and accumulate for
!                 the prolongation coefficient
                  call add_dof(kp,coeff, &
                  NacH_macro(:,kH_macro),ndof_macroH,ConstrH_macro(:,kH_macro))
!
!           ...end of loop through the ultimate modified coarse grid element
!              parent dof
               enddo
!
!        ...end of loop through connected parent edge dof
            enddo
!
!        ...save the # of connected parent dof
            NrconH_macro(kH_macro) = ndof_macroH
!
!     ...loop through nodal H1 dof
         enddo
!
!     ...loop through nodal H(curl) dof
         do j=1,ndofE
!
!        ...dof number in the macro-element
            kE_macro = naE_macro(loc)+j
!
!        ...initiate # of connected modified coarse grid element parent dof
            ndof_macroE = 0
!
!        ...loop through connected parent edge dof
            do je=1,MDOFE
!
!           ...find the corresponding local element dof number
               kE = mapE(je)
!
!           ...loop through the modified coarse grid element parent dof
               do jp=1,nrconE(kE)
                  kp = nacE(jp,kE)
                  coeff = constrE(jp,kE)*NODES_PR(nodl)%coeffE(je,j)
!
                  if (coeff .eq. 0.d0) cycle

!              ...adjust # of connected parent dof and accumulate for
!                 the prolongation coefficient
                  call add_dof(kp,coeff, &
                  NacE_macro(:,kE_macro),ndof_macroE,ConstrE_macro(:,kE_macro))
!
!           ...end of loop through the ultimate modified coarse grid element
!              parent dof
               enddo
!
!        ...end of loop through connected parent edge dof
            enddo
!
!        ...save the # of connected parent dof
            NrconE_macro(kE_macro) = ndof_macroE
!
!     ...loop through nodal H(curl) dof
         enddo
!
!  ...end of loop through nodes in the tree
      enddo
!
!..end of loop through edges
   enddo
!
!---------------------------------------------------------------------------------------
!..Step 3: Faces
!---------------------------------------------------------------------------------------
!
!..loop through the macro-element faces
   do if=1,nrf
!
      ftype = facetype(type,if); nface_or = norientl(nrv+nre+if)
!
!  ...number of vertices on the face (3 or 4)
      nrvf = nvert(ftype)
!
!  ...pick up the local face nodes numbers
      call face_nodes(type,if, nod_face,nrfn)
!
!  ...pick up vertex (and edge) numbers in the face global coordinates
      call face_to_vert_nos(ftype,nface_or, nfver)
!
!  ...pick up edge numbers in the face global coordinates
      call face_to_edge_nos(ftype,nface_or, nfedg)
!
!  ...collect element orientations for the face edges
      do i=1,nrvf
        ie = nod_face(nrvf+i)
        nedge_or(i) = norientl(ie)
      enddo
!
!  ...pick up edge orientations in the face global coordinates
      call face_to_edge_orient(ftype,nface_or,nedge_or, nfedg_or)
!
!  ...pick up the global face nodes numbers in global face coordinates
      do i=1,nrvf
         nodesl_face(i) = nodesl(nod_face(nfver(i)))
         nodesl_face(nrvf+i) = nodesl(nod_face(nrvf+nfedg(i)))
      enddo
      nodesl_face(nrvf+nrvf+1) = nodesl(nrv+nre+if)
!
!  ...determine the face parent nodes order for the prolongation operator
      do i=1,nrvf+1
        nordp(i) = NODES_MG(nodesl_face(nrvf+i))%orderC(Igrid)
      enddo
!
!  ...establish the injections between face and element dof
!  ...vertices of the face
      do i=1,nrvf
         j = nod_face(nfver(i))
         mapH(i) = naHl(j)+1
      enddo
      icH=nrvf; icE=0
!
!  ...edges of the face
      do i=1,nrvf
         j = nod_face(nrvf+nfedg(i))
         nod = nodesl_face(nrvf+i)
         call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                       ndofH,ndofE,nvoid,nvoid)
         do kH=1,ndofH
            mapH(icH+kH) = naHl(j)+kH
         enddo
         do kE=1,ndofE
            mapE(icE+kE) = naEl(j)+kE
         enddo
         icH = icH + ndofH
         icE = icE + ndofE
      enddo
      nod = nodesl(nrv+nre+if)
!
      call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC(Igrid), &
                    ndofH,ndofE,ndofV,nvoid)
!
!  ...the mid-face node
      j = nod_face(nrvf+nrvf+1)
      do kH=1,ndofH
         mapH(icH+kH) = naHl(j)+kH
      enddo
      do kE=1,ndofE
         mapE(icE+kE) = naEl(j)+kE
      enddo
      do kV=1,ndofV
         mapV(kV) = naVl(j)+kV
      enddo
!
      icH = icH + ndofH; icE = icE + ndofE
!
      call build_face_tree(ftype,nordp,nfedg_or,nodesl_face)
!
!  ...loop through nodes in the tree
      do nodl=1,NR_NODES_PR
!
!     ...global node number
         nod = NODES_PR(nodl)%nod
!
!     ...skip if not on the list of macro-element nodes
         call locate(nod,Nod_macro,Nrnod_macro, loc)
         if (loc.eq.0) cycle
!
!     ...skip if on the list of modified coarse element node
         call locate(nod,Nodm,Nrnodm, locC)
         if (locC.ne.0) cycle
!
!     ...skip if the node is not in the interior of the face
!     ...(i.e. it has been considered in the loop through edges)
         if(.not. Face_interior_node(nodl)) cycle
!
!     ...get the # number of dof implied by the parent (coarse) order
         call ndof_nod(NODES(nod)%type,NODES_PR(nodl)%order, &
                       ndofH,ndofE,ndofV,nvoid)

         call inject_dof(NODES(nod)%type, NODES_PR(nodl)%order,NODES(nod)%order,   &
                         InjH, InjE, InjV)
!
!     ...loop through nodal H1 dof
         do j=1,ndofH
!
!        ...dof number in the macro-element
            kH_macro = naH_macro(loc)+injH(j)
!
!        ...initiate # of connected modified coarse grid element parent dof
            ndof_macroH = 0
!
!        ...loop through connected parent dof
            do je=1,MDOFH
!
!           ...find the corresponding local element dof number
               kH = mapH(je)
!
!           ...loop through the modified coarse grid element parent dof
               do jp=1,nrconH(kH)
                  kp = nacH(jp,kH)
                  coeff = constrH(jp,kH)*NODES_PR(nodl)%coeffH(je,j)
!
                  if (coeff .eq. 0.d0) cycle

!              ...adjust # of connected parent dof and accumulate for
!                 the prolongation coefficient
                  call add_dof(kp,coeff, &
                  NacH_macro(:,kH_macro),ndof_macroH,ConstrH_macro(:,kH_macro))
!
!           ...end of loop through the ultimate modified coarse grid element
!              parent dof
               enddo
!
!        ...end of loop through connected parent edge dof
            enddo
!
!        ...save the # of connected parent dof
            NrconH_macro(kH_macro) = ndof_macroH
!
!     ...end of loop through nodal H1 dof
         enddo
!
!     ...loop through nodal H(curl) dof
         do j=1,ndofE
!
!        ...dof number in the macro-element
            kE_macro = naE_macro(loc)+injE(j)
!
!        ...initiate # of connected modified coarse grid element parent dof
            ndof_macroE = 0
!
!        ...loop through connected parent edge dof
            do je=1,MDOFE
!
!           ...find the corresponding local element dof number
               kE = mapE(je)
!
!           ...loop through the modified coarse grid element parent dof
               do jp=1,nrconE(kE)
                  kp = nacE(jp,kE)
                  coeff = constrE(jp,kE)*NODES_PR(nodl)%coeffE(je,j)
!
                  if (coeff .eq. 0.d0) cycle
!
!              ...adjust # of connected parent dof and accumulate for
!                 the prolongation coefficient
                  call add_dof(kp,coeff, &
                  NacE_macro(:,kE_macro),ndof_macroE,ConstrE_macro(:,kE_macro))
!
!           ...end of loop through the ultimate modified coarse grid element
!              parent dof
               enddo
!
!        ...end of loop through connected parent edge dof
            enddo
!
!        ...save the # of connected parent dof
            NrconE_macro(kE_macro) = ndof_macroE
!
!     ...end of loop through nodal H(curl) dof
         enddo
!
!     ...loop through nodal H(div) dof
         do j=1,ndofV
!
!        ...dof number in the macro-element
            kV_macro = naV_macro(loc)+injV(j)
!
!        ...initiate # of connected modified coarse grid element parent dof
            ndof_macroV = 0
!
!        ...loop through connected parent edge dof
            do je=1,MDOFV
!
!           ...find the corresponding local element dof number
               kV = mapV(je)
!
!           ...loop through the modified coarse grid element parent dof
               do jp=1,nrconV(kV)
                  kp = nacV(jp,kV)
                  coeff = constrV(jp,kV)*NODES_PR(nodl)%coeffV(je,j)

                  if (coeff .eq. 0.d0) cycle
!
!              ...adjust # of connected parent dof and accumulate for
!                 the prolongation coefficient
                  call add_dof(kp,coeff, &
                  NacV_macro(:,kV_macro),ndof_macroV,ConstrV_macro(:,kV_macro))
!
!           ...end of loop through the ultimate modified coarse grid element
!              parent dof
               enddo
!
!        ...end of loop through connected parent edge dof
            enddo
!
!        ...save the # of connected parent dof
            NrconV_macro(kV_macro) = ndof_macroV
!
!     ...loop through nodal H(div) dof
         enddo
!
!  ...end of loop through nodes in the tree
      enddo
!..end of loop through faces
   enddo
!
!
   end subroutine logic_macro
!
!
!------------------------------------------------------------------------
!
!    routine name      - add_dof
!
!------------------------------------------------------------------------
!
!    latest revision   - Jan 2018
!
!    purpose           - Routine updates prolongation (constraint)
!                        coefficients for the macro-element dofs
!
!   arguments :
!     in:
!         Kp           - a modified coarse grid element dof
!         Coeff        - an update to the constraint coefficient
!   in/out:
!         Nac_macro    - current list of modified coarse grid element
!                        dof connected to a macro-element dof
!         Ndof_macro   - current # of entries in Nac_macro
!       Constr_macro   - the corresponding list of constraint
!                        coefficients
!
!-----------------------------------------------------------------------
!
   subroutine add_dof(Kp,Coeff, Nac_macro,Ndof_macro,Constr_macro)
!
   implicit none
   integer, intent(in)    :: Kp
   integer, intent(inout) :: Nac_macro(*),Ndof_macro
   real*8,  intent(in)    :: Coeff
   real*8,  intent(inout) :: Constr_macro(*)
!
!..locals
   integer :: loc
!
!-----------------------------------------------------------------------
!
!..locate the parent dof on the list of connected parent dof so far
   call locate(Kp,Nac_macro,Ndof_macro, loc)
   if (loc.eq.0) then
      Ndof_macro = Ndof_macro+1
      loc = Ndof_macro
      Nac_macro(loc) = Kp
      Constr_macro(loc) = 0.d0
   endif
   Constr_macro(loc) = Constr_macro(loc) + Coeff

!
   end subroutine add_dof
!

