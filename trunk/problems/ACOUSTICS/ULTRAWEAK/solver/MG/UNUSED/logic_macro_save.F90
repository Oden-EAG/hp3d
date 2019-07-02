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
!           MdleC      - middle node of a coarse element
!           Nod_macro  - macro-element nodes
!           Nr_macro   - # of macro-element nodes
!             
!     out:     
!           Nrnodm     - number of modified coarse grid element nodes
!           Nodm       - nodes of the modified coarse grid element
!                        excluding the middle node
!           NdofmH,NdofmE,NdofmV - the corresponding number of
!                        H1, H(curl),H(div) single component dof
!
!-----------------------------------------------------------------------
!
   subroutine logic_macro(MdleC,&
                          Nod_macro,Nr_macro, &
                          Nrnodm,Nodm,NdofmH, NdofmE, NdofmV)
!
   use data_structure3D
   use refinements
   use prolongation

   IMPLICIT NONE

   integer, intent(in)  :: MdleC,Nod_macro(Nr_macro),Nr_macro
   integer, intent(out) :: Nrnodm
   integer, intent(out) :: Nodm(MAXNODM),   NdofmH(MAXNODM),  &
                           NdofmE(MAXNODM), NdofmV(MAXNODM)

!..element type
   character(len=4) :: type, facetype, ftype
!..workspace for elem_nodes
   integer          :: nodesl(27),norientl(27),norder(19)
!..workspace for face_to_vert_nos and face_to_vert_orient
   integer          :: nedge_or(4), nfedg_or(4) , nfver(4)
!
   integer          :: icH,icE,icV,im,nod,nvoid
!
!..offsets for first dof for coarse element nodes
   integer          :: naHl(26),naEl(26),naVl(26)
!
!..offsets for first dof for modified coarse element nodes
   integer          :: naH(MAXNODM),naE(MAXNODM),naV(MAXNODM)
!
!..offsets for first dof for macro-element element nodes
   integer          :: naH_macro(Nr_macro),naE_macro(Nr_macro),naV_macro(Nr_macro)   
!   
!..number of H1, H(div), H(curl), single component dof for a macro-element
   integer          :: ndofH_macro(Nr_macro),     &
                       ndofE_macro(Nr_macro),ndofV_macro(Nr_macro)
!
!..total number of dof for a macro-element
   integer          :: nrdofH_macro, nrdofE_macro, nrdofV_macro
!
!..workspace for logicC
   integer   :: idec, NrconH,NrconE,NrconV,NacH,NacE,NacV
   real*8    :: ConstrH,ConstrE,ConstrV
   dimension :: nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV),   &
                nacH(NACDIM,MAXbrickH),constrH(NACDIM,MAXbrickH),        &
                nacE(NACDIM,MAXbrickE),constrE(NACDIM,MAXbrickE),        &
                nacV(NACDIM,MAXbrickV),constrV(NACDIM,MAXbrickV)
!   
!..connectivity arrays to be returned
   integer, parameter   :: Max_macroH = 500, Max_macroE = 500, Max_macroV = 500
   integer   :: NrconH_macro, NrconE_macro, NrconV_macro
   integer   :: NacH_macro, NacE_macro, NacV_macro
   real*8    :: ConstrH_macro, ConstrE_macro, ConstrV_macro, coeff
   dimension :: NrconH_macro(Max_macroH),                                       &
                NrconE_macro(Max_macroE), NrconV_macro(Max_macroV),             &
                NacH_macro(NACDIM,Max_macroH),ConstrH_macro(NACDIM,Max_macroH), &
                NacE_macro(NACDIM,Max_macroE),ConstrE_macro(NACDIM,Max_macroE), &
                NacV_macro(NACDIM,Max_macroV),ConstrV_macro(NACDIM,Max_macroV)
!
!..maps establishing injections between edge or face dof and element dof
   integer   :: mapH(MAXquadH),mapE(MAXquadE),mapV(MAXquadQ)
!
   integer   :: kH, kE, kV, j, loc, loc1, locC
   integer   :: kH_macro, kE_macro , kV_macro
   integer   :: ndof_macroH, ndof_macroE , ndof_macroV
   integer   :: nrv, nre, nrf, nrvf, nrfn, iv
   integer   :: ie, i, iedg(2), nodes_edge(3), nodl
   integer   :: if, ief, ivf(4), nod_face(9), nodesl_face(9)
   integer   :: nord, nordh, nordv
   integer   :: kHC, ndofH, ndofE, ndofV
   integer   :: je, jp, kp
   integer   :: nface_or

   integer   :: iprint              
!
! 
   iprint = 0
!
!-----------------------------------------------------------------------
!
   type = NODES(MdleC)%type
   nrv = nvert(type); nre = nedge(type); nrf = nface(type)
!   
!..get element nodes and build the local data base on constrained
!..nodes for the element
   call get_connect_infoC(MdleC, nodesl,norientl)
!
!..get nodes of the modified coarse element and the corresponding
!  constraints coefficients
   idec = 1
   call logicC(MdleC,idec,                        &
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
      call ndof_nod(NODES(nod)%type,NODES(nod)%orderC, &
                    ndofH,ndofE,ndofV,nvoid)
      icH=icH+ndofH; icE=icE+NdofE; icV=icV+NdofV
   enddo
!
!..establish offsets for the modified coarse element nodal dof 
!  USING THE REVERSE ORDER
   icH=0; icE=0; icV=0
   do im = Nrnodm,1,-1
      naH(im)=icH; naE(im)=icE; naV(im)=icV
      nod = Nodm(im)
      call ndof_nod(NODES(nod)%type,NODES(nod)%orderC, &
                    NdofmH(im),NdofmE(im),NdofmV(im),nvoid)
      icH=icH+NdofmH(im); icE=icE+NdofmE(im); icV=icV+NdofmV(im)
   enddo
!
   if (iprint.eq.1) then
      write(*,7001) MdleC
 7001 format('logic_macro: MODIFIED ELEMENT NODES FOR MdleC = ',i6)
      do im = Nrnodm,1,-1
         nod = Nodm(im)
         write(*,7002) im,nod,NODES(nod)%type,NODES(nod)%orderC, &
                       naH(im),NdofmH(im), &
                       naE(im),NdofmE(im), &
                       naV(im),NdofmV(im)
 7002 format(' im  = ',i3,' nod = ',i6,' TYPE = ',a6,' ORDER = ',i3, &
             ' naH = ',i3,' NdofmH = ',i3, &
             ' naE = ',i3,' NdofmE = ',i3, &
             ' naV = ',i3,' NdofmV = ',i3)
      enddo
   endif
!
!..establish offsets for the macro-element nodal dof
   icH=0; icE=0; icV=0
   do im=1,Nr_macro
      naH_macro(im)=icH; naE_macro(im)=icE; naV_macro(im)=icV
      nod = Nod_macro(im)
      call ndof_nod(NODES(nod)%type,NODES(nod)%order,  &
                    NdofH_macro(im),NdofE_macro(im),NdofV_macro(im),nvoid)
      icH=icH+NdofH_macro(im); icE=icE+NdofE_macro(im); icV=icV+NdofV_macro(im)
   enddo
   nrdofH_macro = icH; nrdofE_macro = icE; nrdofV_macro = icV
   if (iprint.eq.1) then
      write(*,8012) MdleC
 8012 format('logic_macro: MACRO-ELEMENT NODES FOR MdleC = ',i6)
      do im=1,Nr_macro
         nod = Nod_macro(im)
         write(*,7002) im,nod,NODES(nod)%type,NODES(nod)%order, &
                       naH_macro(im),NdofH_macro(im), &
                       naE_macro(im),NdofE_macro(im), &
                       naV_macro(im),NdofV_macro(im)
      enddo
      call pause
   endif
!
!-----------------------------------------------------------------------
!
!..initiate the connectivities arrays
   NrconH_macro = 0; NacH_macro = 0; ConstrH_macro = 0.d0
   NrconE_macro = 0; NacE_macro = 0; ConstrE_macro = 0.d0
   NrconV_macro = 0; NacV_macro = 0; ConstrV_macro = 0.d0
!
!-----------------------------------------------------------------------
!
!..loop through the coarse grid element  vertices (they may also belong
!  to macro-element)
   do iv=1,nrv
      nod = nodesl(iv)
!
!  ...skip if not on the list of macro-element nodes
      call locate(nod,Nod_macro,Nr_macro, loc)
      if (loc.eq.0) cycle
      kH = naH_macro(loc)+1  !  dof # within the macro-element
!
!  ...skip if not on the list of modified element nodes
      call locate(nod,Nodm,Nrnodm, locC)
      if (locC.eq.0) cycle
!
!  ...copy prolongation coefficients from those delivered by logicC
      kHC = naH(locC)+1  ! dof # within the coarse grid element
      NrconH_macro(kH) = NrconH(kHC)
      NacH_macro(1:NrconH(kHC),kH) = nacH(1:NrconH(kHC),kHC)
      ConstrH_macro(1:NrconH(kHC),kH) = constrH(1:NrconH(kHC),kHC)
!
!..end of loop through coarse grid element vertices
   enddo
!
!-----------------------------------------------------------------------
!
!..loop through the macro-element edges
   do ie=1,nre
!
!  ...local mid-edge node number
      i = nrv + ie
!
!  ...mid-edge node     
      nodes_edge(3) = nodesl(i)
!      
!  ...vertex nodes (account for edge orientation)
      call edge_to_vert(type,ie, iedg(1),iedg(2))
      if (norientl(i).eq. 1) call swap(iedg(1),iedg(2))
      nodes_edge(1:2) = nodesl(iedg(1:2))
!
!  ...establish the injections between edge and element dof
      mapH(1:2) = naHl(iedg(1:2))+1
      nod = nodes_edge(3)
      call ndof_nod(NODES(nod)%type,NODES(nod)%orderC, &
                    ndofH,ndofE,nvoid,nvoid)
      do kH=1,ndofH
        mapH(2+kH) = naHl(i)+kH
      enddo
      do kE=1,ndofE
        mapE(kE) = naEl(i)+kE
      enddo
!
      call build_edge_tree(nodes_edge)
!
!  ...loop through local nodes in the tree
      do nodl=1,NR_NODES_PR
!
!     ...global node number
         nod = NODES_PR(nodl)%nod
         call locate(nod,Nod_macro,Nr_macro, loc)
!
!     ...skip if not on the list of macro-element nodes
         if (loc.eq.0) cycle
!    
!     ...skip if the node is not in the interior of the edge
!     ...(i.e. it has been considered in the loop through vertices)
         if (.not. Edge_interior_node(nodl)) cycle         
!
         call ndof_nod(NODES(nod)%type,NODES(nod)%orderC, &
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
!     ...end of loop through nodes in the tree
      enddo
!
!..end of loop through edges
   enddo
!
!-----------------------------------------------------------------------
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
!  ...collect element orientations for the face edges
      do i=1,nrvf
        ie = nod_face(nrvf+i)
        nedge_or(i) = norientl(ie)
      enddo
!
!  ...pick up vertex (and edge) numbers in the face global coordinates
      call face_to_vert_nos(ftype,nface_or, nfver)
!
!  ...pick up edge orientations in the face global coordinates
      call face_to_edge_orient(ftype,nface_or,nedge_or, nfedg_or)
!
!  ...pick up the global face nodes numbers in global face coordinates
      do i=1,nrvf
         nodesl_face(i) = nodesl(nod_face(nfver(i)))
         nodesl_face(nrvf+i) = nodesl(nod_face(nrvf+nfver(i)))
      enddo
      nodesl_face(nrvf+nrvf+1) = nodesl(nrv+nre+if)
!
!  ...determine the face order for the prolongation operator
      nod = nodesl(nrv+nre+if)
      select case(facetype(type,if))
      case('mdlt')
         nord = NODES(nod)%orderC
         do ie=1,3
            nod = nodesl_face(3+ie)
            nord = max(nord,NODES(nod)%orderC)
         enddo
      case('mdlq')
         call decode(NODES(nod)%orderC, nordh,nordv)
         do ie=1,3,2
            nod = nodesl_face(4+ie)
            nordh = max(nordh,NODES(nod)%orderC)
         enddo
         do ie=2,4,2
            nod = nodesl_face(4+ie)
            nordv = max(nordv,NODES(nod)%orderC)
         enddo
         nord = nordh*10+nordv
      end select
!
!  ...establish the injections between face and element dof
      do i=1,nrvf
         mapH(i) = naHl(nod_face(nfver(i)))+1
      enddo
      icH=nrvf; icE=0
      do i=1,nrvf
         j = nod_face(nrvf+nfver(i))
         nod = nodesl(j)
         call ndof_nod(NODES(nod)%type,NODES(nod)%orderC, &
                       ndofH,ndofE,nvoid,nvoid)
         do kH=1,ndofH
            mapH(icH+kH) = naHl(j)+kH
         enddo
         do kE=1,ndofE
            mapE(kE) = naEl(j)+kE
         enddo
         icH = icH + ndofH
         icE = icE + ndofE
      enddo
      nod = nodesl(nrv+nre+if)
      call ndof_nod(NODES(nod)%type,NODES(nod)%orderC, &
                    ndofH,ndofE,ndofV,nvoid)
      do kH=1,ndofH
         mapH(icH+kH) = naHl(j)+kH
      enddo
      do kE=1,ndofE
         mapE(kE) = naEl(j)+kE
      enddo
      do kV=1,ndofV
         mapV(kV) = naVl(j)+kV
      enddo
!
      call build_face_tree(ftype,nord,nfedg_or,nodesl_face)
!
!  ...loop through nodes in the tree
      do nodl=1,NR_NODES_PR
!
!     ...global node number
         nod = NODES_PR(nodl)%nod
!
!     ...skip if not on the list of macro-element nodes
         call locate(nod,Nod_macro,Nr_macro, loc)
         if (loc.eq.0) cycle
!
!     ...skip if the node has already been considered in the first step
!     ...of the algorithm (i.e., it is common to both coarse-grid element
!     ...and the macro-element)
         call locate(nod,Nodm,Nrnodm, loc1)
         if (loc1.ne.0) cycle
!    
!     ...skip if the node is not in the interior of the face 
!     ...(i.e. it has been considered in the loop through edges)
         if(.not. Face_interior_node(nodl)) cycle         
!
!     ...end of loop through nodes in the tree
      enddo
!    ...end of loop through faces
   enddo      
!
!
!-----------------------------------------------------------------------
!  PRINT
!-----------------------------------------------------------------------
!
   if (iprint.ge.2) then
      write(*,7011) MdleC
 7011 format('logic: MdleC = ',i6)
      do i=1,Nr_macro
         nod = Nod_macro(i)
         write(*,7012) i,nod
 7012 format('logicC: i = ',i3,' nod = ',i5)
         do j=1,NdofH_macro(i)
            kH = naH_macro(i)+j
            write(*,7013) kH,NrconH(kH)
 7013 format('logicC: kH = ',i4,' NrconH = ',i3)
            write(*,7014) NacH(1:NrconH(kH),kH)
 7014 format(20i6)
            write(*,7015) ConstrH(1:NrconH(kH),kH)
 7015 format(20f6.3)
         enddo
         do j=1,NdofE_macro(i)
            kE = naE_macro(i)+j
            write(*,7022) kE,NrconE(kE)
 7022 format('logicC: kE = ',i4,' NrconE = ',i3)
            write(*,7023) NacE(1:NrconE(kE),kE)
 7023 format(20i6)
            write(*,7024) ConstrE(1:NrconE(kE),kE)
 7024 format(20f6.3)
         enddo
      enddo
      call pause
   endif



   end subroutine logic_macro
!
!-----------------------------------------------------------------------
!
!    routine name      - print_nodesPR
!
!-----------------------------------------------------------------------
!
!    latest revision   - Jan 2018
!
!    purpose           - 
!
!   arguments :
!
!----------------------------------------------------------------------
!
   subroutine print_nodesPR
!
   use prolongation

   implicit none 

   integer :: i, nrs,j, sonl(9), iH, jH, iE, jE, iV, jV, temp(2)


!..print out NODES_PR
   do i=1,NR_NODES_PR
      write (*,6001) i,NODES_PR(i)%type,NODES_PR(i)%nod,       &
                       NODES_PR(i)%ref_kind,NODES_PR(i)%father
 6001 format('  i     = ',i4,'    type = ',a4,'    nod = ',i5,     &
                      '    ref_kind = ',i2,'    father = ',i4)
!
 !      if (NODES_PR(i)%ref_kind .ne. 0) then
 !         call nr_sons(NODES_PR(i)%type,NODES_PR(i)%ref_kind, nrs)
 !         write(*,6002) NODES_PR(i)%sons(1:nrs)
 ! 6002    format('sons      = ',9i6)
 !         write(*,6003) NODES_PR(NODES_PR(i)%sons(1:nrs))%nod
 ! 6003    format('sons nod  = ',9i6)
 !         write(*,*)
 !      endif

      temp =  ubound(NODES_PR(i)%coeffH); iH = temp(1); jH = temp(2)
      temp =  ubound(NODES_PR(i)%coeffE); iE = temp(1); jE = temp(2) 
      temp =  ubound(NODES_PR(i)%coeffV); iV = temp(1); jV = temp(2)
! 
      do j = 1, jH
         write(*,6004) j, NODES_PR(i)%coeffH(1:iH,j)
 6004    format(' jH     = ', i4,/,' coeffH = ', 5e13.4,/,5(10x,5e13.4,/))
      enddo        
      do j = 1, jE
         write(*,6005) j, NODES_PR(i)%coeffE(1:iE,j)
 6005    format(' jE     = ', i4,/,' coeffE = ', 5e13.4,/,5(10x,5e13.4,/))
      enddo
      do j = 1, jV
         write(*,6006) j, NODES_PR(i)%coeffV(1:iV,j)
 6006    format(' jV     = ', i4,/,' coeffV = ', 5e13.4,/,5(10x,5e13.4,/))
      enddo
      write(*,*)
      call pause
   enddo
!
!
   end subroutine print_nodesPR



   function Facetype(Type,If)
!
   character(len=4) :: Facetype,Type
!
   select case(Type)
   case('mdlp')
      select case(If)
      case(1,2)
         Facetype = 'mdlt'
      case(3,4,5)
         Facetype = 'mdlq'
      end select
   case('mdlb')
      Facetype = 'mdlq'
   case('mdln')
      Facetype = 'mdlt'
   case('mdld')
      select case(If)
      case(1)
         Facetype = 'mdlq'
      case(2,3,4,5)
         Facetype = 'mdlt' 
      end select
   end select
   end function Facetype

!..Routine updates prolongation (constraint) coefficients for
!  the macro-element dofs
!  in :
!       Kp           - a modified coarse grid element dof
!       Coeff        - an update to the constraint coefficient
!  in/out:
!       Nac_macro    - current list of modified coarse grid element
!                      dof connected to a macro-element dof
!       Ndof_macro   - current # of entries in Nac_macro
!       Constr_macro - the corresponding list of constraint 
!                      coefficients

   subroutine add_dof(Kp,Coeff, Nac_macro,Ndof_macro,Constr_macro)
!
   implicit none
   integer :: Kp,Nac_macro(*),Constr_macro(*),Ndof_macro
   real*8  :: Coeff
!
!..locals
   integer :: loc
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