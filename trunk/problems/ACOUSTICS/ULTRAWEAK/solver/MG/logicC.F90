!----------------------------------------------------------------------
!
!   routine name       - logicC
!
!----------------------------------------------------------------------
!
!   latest revision    - Jan 18
!
!
!   purpose            - routine establishes list of nodes
!                        for a modified element, and the corresponding
!                        constraints coefficients 
!                        (taking to account the master flag)
!
!   arguments :
!     in:
!            Mdle      - an element number, same as the middle node
!                        number
!            Idec      = 1 if only the element nodes are to be
!                          determined
!     out:
!            Nodm      - actual (unconstrained) nodes returned in the
!                        standard order: vertex, mid-edge, mid-face
!                        and middle nodes
!            NdofmH,NdofmE,NdofmV - the corresponding number of
!                        H1, H(curl) and H(div) dof
!            Nrnodm    - number of modified element nodes
!
!            NrconH    - number of parent dof's for an i-th local H1 dof
!            NacH(j,i),j=1,...,NrconH(i) - list of parent dof's
!                        for i-th local dof
!            ConstrH(j,i),j=1,...,nrcon(i) - the corresponding
!                        constraints coefficients
!
!            NrconE,NacE,ConstrE - same for H(curl) dof
!            NrconV,NacV,ConstrV - same for H(div) dof
!
!----------------------------------------------------------------------
!
   subroutine logicC(Mdle,Idec,                         &
                    Nodm,NdofmH,NdofmE,NdofmV,Nrnodm,  &
                    NrconH,NacH,ConstrH,               &
                    NrconE,NacE,ConstrE,               &
                    NrconV,NacV,ConstrV)
!
   use element_data
   use constraints
   use data_structure3D
   use constrained_nodes
   use refinements
   use mg_data_structure
#include "syscom.blk"
!
   dimension Nodm(MAXNODM),NdofmH(MAXNODM),NdofmE(MAXNODM),    &
                           NdofmV(MAXNODM),NdofmQ(MAXNODM)
!
   dimension NrconH(MAXbrickH),NrconE(MAXbrickE),NrconV(MAXbrickV),   &
             NacH(NACDIM,MAXbrickH),ConstrH(NACDIM,MAXbrickH),        &
             NacE(NACDIM,MAXbrickE),ConstrE(NACDIM,MAXbrickE),        &
             NacV(NACDIM,MAXbrickV),ConstrV(NACDIM,MAXbrickV)
!
!..middle node type
   character(len=4) :: etype
 7002     format('im = ',i3,' nod = ',i6,' TYPE = ',a6,' ORDER = ',i3,   &
           ' naH = ',i3,' NdofmH = ',i3,  &
           ' naE = ',i3,' NdofmE = ',i3,  &
           ' naV = ',i3,' NdofmV = ',i3)
 7030       format('logic: ACTIVE NODE      i = ',i3,' nod = ',i6,  &
                   ' TYPE = ',a6,' ORDER = ',i2)
 7031       format('logic: CONSTRAINED NODE i = ',i3,' nod = ',i6,  &
                   ' TYPE = ',a6,' ORDER = ',i2,' icase = ',i2,     &
                   ' EDGE NODES = ',i5,2x,2i5)
 7032       format('logic: CONSTRAINED NODE i = ',i3,' nod = ',i6,  &
                   ' TYPE = ',a6,' ORDER = ',i2,' icase = ',i2,     &
                   ' FACE NODES = ',i5,2x,4i5,2x,4i5)
!
!..element local nodes, orientations
   dimension nodesl(27),norientl(27)
!
!..offsets for first dof for modified element nodes
   dimension naH(MAXNODM),naE(MAXNODM),naV(MAXNODM)
!
!..order for the element
   dimension norder(19)
!
!..addresses and number of dof for parent mid-edge and vertex nodes
   dimension locp(8),ndofHp(8),ndofEp(8)
!
!..unconstrained/constrained edge flags
!  = 0          if edge is unconstrained (active)
!  = son number if edge is constrained (inactive)
   dimension ncflag(4)
!
!----------------------------------------------------------------------
!
!..initialize constraints arrays
   if (.NOT. INITIALIZED_CONSTR) call init_cnstr
!
   select case(Mdle)
   case(1084)
     if (Idec.eq.2) then
       iprint=0
     else
       iprint=0
     endif
   case default
     iprint=0
   endselect
!
   if (iprint.ge.1) then
      write(*,7000) Mdle
 7000 format(' logicC: Mdle = ',i6,' --')
   endif
!
!..get element nodes and build the local data base on constrained
!..nodes for the element
   call get_connect_infoC(Mdle, nodesl,norientl)
   if (iprint.eq.1) then
      call elem_show(Mdle, NODES(Mdle)%type,nodesl,norientl)
   endif
!
!..get nodes of the modified element
   call logic_nodesC(Mdle,nodesl, Nodm,Nrnodm)
!
!
!..establish offsets for the modified element nodal dof
!
!..initialize
   icH=0 ; icE=0 ; icV=0
!
!..loop over nodes of modified element
   do im=1,Nrnodm
!
!  ...number of H1, H(curl), H(div) dofs stored so far
      naH(im)=icH ; naE(im)=icE ; naV(im)=icV
!
!  ...number of H1, H(curl), H(div) dofs associated to node
      nod=Nodm(im)
      call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC,         &
                    NdofmH(im),NdofmE(im),NdofmV(im),NdofmQ(im))
!
!  ...update H1,H(curl),H(div) dofs counters
      icH=icH+NdofmH(im) ; icE=icE+NdofmE(im) ; icV=icV+NdofmV(im)
   enddo
!
!  ...printing
   if (iprint.eq.1) then
      write(*,*) 'logicC: MODIFIED ELEMENT NODES:'
      do im=1,Nrnodm
         nod=Nodm(im)
         write(*,7002) im,nod,NODES(nod)%type,NODES_MG(nod)%orderC,    &
                       naH(im),NdofmH(im),naE(im),NdofmE(im),  &
                       naV(im),NdofmV(im)
      enddo
   endif
!
!----------------------------------------------------------------------
   if (Idec.eq.1) go to 100
!
!..initialize
   NrconH=0 ; NacH=0 ; ConstrH=0.d0
   NrconE=0 ; NacE=0 ; ConstrE=0.d0
   NrconV=0 ; NacV=0 ; ConstrV=0.d0
!
!..construct the constrained approximation info
!
!..initialize local H1, H(curl), H(div) dof counters
   kH=0 ; kE=0 ; kV=0
!
!..number of element's nodes
   etype=NODES(Mdle)%type
   nrnodes=nvert(etype)+nedge(etype)+nface(etype)+1
!
!..loop over element nodes
   do i=1,nrnodes
      nod=nodesl(i)
!
!======================================================================
!       NODE IS UNCONSTRAINED                                         |
!======================================================================
      if (NODES_MG(nod)%master .eq. 1) then
!
!     ...node location on the list of modified element's nodes
         call locate(nod,Nodm,Nrnodm, loc)
!
!     ...H1-loop through the nodal dof
         do j=1,NdofmH(loc)
!
!        ...advance dof counter
            kH=kH+1
!
!        ...number of parent dof's
            NrconH(kH)=1
!
!        ...parent dof location
            NacH(1,kH)=naH(loc)+j
!
!        ...constraint coefficient
            ConstrH(1,kH)=1.d0
         enddo
!
!     ...H(curl)-loop through the nodal dof
         do j=1,NdofmE(loc)
            kE=kE+1
            NrconE(kE)=1
            NacE(1,kE)=naE(loc)+j
            ConstrE(1,kE)=1.d0
         enddo
!
!     ...H(div)-loop through the nodal dof
         do j=1,NdofmV(loc)
            kV=kV+1
            NrconV(kV)=1
            NacV(1,kV)=naV(loc)+j
            ConstrV(1,kV)=1.d0
         enddo
!
!     ...printing
         if (iprint.eq.1) then
            write(*,7030) i,nod,NODES(nod)%type,NODES_MG(nod)%orderC
         endif
!
!======================================================================
!       NODE IS CONSTRAINED                                           |
!======================================================================
      else
!
!     ...compute number of dof for the node
         call ndof_nod(NODES(nod)%type,NODES_MG(nod)%orderC,    &
                             ndofH,ndofE,ndofV,ndofQ)
!
!     ...identify the node case
         call decode2(NODES_CONSTR(i), nc,icase)
!
!     ...printing
         if (iprint.eq.1) then
            select case(icase)
!
!        ...node constrained by an edge
            case(11,12,13,37,38,39,47,48,49)
               write(*,7031) i,nod,NODES(nod)%type,NODES_MG(nod)%orderC,icase,   &
                             NEDGC(nc),NEDG_CONS(1:2,nc)
!
!           ...node constrained by a face
            case default
               write(*,7032) i,nod,NODES(nod)%type,NODES_MG(nod)%orderC,icase,   &
                             NFACEC(nc),NFACE_CONS(1:8,nc)
            end select
         endif
!
         select case(icase)
!
!-----------------------------------------------------------------------
!        H2-REFINED EDGE
!-----------------------------------------------------------------------
!
!     -- FIRST AND SECOND MID-EDGE NODE CONSTRAINED BY AN EDGE --
!
         case(11,12, 37,38, 47,48)
!
!        ...local son number
            select case(icase)
            case(11,37,47) ; is=1
            case(12,38,48) ; is=2
            endselect
!
!        ...parent mid-edge node
            nodp = NEDGC(nc)
!
!        ...parent node location on the list of modified element's nodes
            call locate(nodp,Nodm,Nrnodm, loc)
!
!           H1-loop through nodal dof's
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
!
!           ...advance dof counter
               kH=kH+1
!
!           ...number of parent dof's
               NrconH(kH) = ndofH
!
!           ...loop through parent node dof's
               do jp=1,NdofmH(loc)
!
!              ...parent dof location
                  NacH(jp,kH)=naH(loc)+jp
!
!              ...constraint coefficient
                  ConstrH(jp,kH)=RRRH(1,jp,is,l)
               enddo
            enddo
!
!           H(curl)-loop through nodal dof's
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
               kE=kE+1
               NrconE(kE) = ndofE
!
!           ...loop through parent node dof's
               do jp=1,NdofmE(loc)
                  NacE(   jp,kE)=naE(loc)+jp
                  ConstrE(jp,kE)=RRRE(1,jp,is,l)
               enddo
            enddo
!
!           NO H(div) constraints across edges
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!     -- VERTEX NODE CONSTRAINED BY AN EDGE --
!
         case(13,39,49)
!
!           H1 constraint ONLY
!           ~~~~~~~~~~~~~~~~~~
!
!        ...advance dof counter
            kH=kH+1
!
!        ...no parent dof's
            NrconH(kH)=0
!
!        ...parent mid-edge node
            nodp = NEDGC(nc)
!
!        ...number of parent node dof's
            ndofHpp = NODES_MG(nodp)%orderC-1
!
!        ...number of parent dof's
            NrconH(kH) = NrconH(kH)+ndofHpp
!
!        ...parent node location on the list of modified element's nodes
            call locate(nodp,Nodm,Nrnodm, loc)
!
!        ...initialize parent dof local number
            jp=0
!
!        ...STEP 1 : loop through parent mid-edge node dofs
            do lp=1,ndofHpp
!
!           ...advance parent dof local number
               jp=jp+1
!
!           ...parent dof location
               NacH(jp,kH)=naH(loc)+lp
!
!           ...constraint coefficient
               ConstrH(jp,kH)=RRRH(1,lp,3,1)
            enddo
!
!        ...STEP 2 : loop through parent edge vertex dofs                        
            do iv=1,2
!
!           ...parent vertex node
               nodp = NEDG_CONS(iv,nc)
!
!           ...active (unconstrained) vertex node
               if (nodp.gt.0) then
!
!              ...update number of parent dof's
                  NrconH(kH) = NrconH(kH)+1
!
!              ...vertex node location on list of modified element's nodes
                  call locate(nodp,Nodm,Nrnodm, loc)
!
!              ...advance parent dof local number
                  jp=jp+1
!
!              ...parent dof location
                  NacH(jp,kH)=naH(loc)+1
!
!              ...constraint coefficient
                  ConstrH(jp,kH)=RRRH(1+iv,1,3,1)
!
!              ...inactive (constrained) vertex node
               else
!
!              ...grandparent edge node
                  call decode2(-nodp, nce,nvoid)
                  nodp = NEDGC(nce)
!
!              ...number of grandparent node dof's
                  ndofHpp = NODES_MG(nodp)%orderC-1
!
!              ...update number of parent dof's
                  NrconH(kH) = NrconH(kH)+ndofHpp
!
!              ...grandparent node location on list of modified element's nodes
                  call locate(nodp,Nodm,Nrnodm, loc)
!
!              ...STEP 2.1 : loop through grandparent mid-edge node dofs
                  do lp=1,ndofHpp
!
!                 ...advance parent dof local number
                     jp=jp+1
!
!                 ...parent dof location
                     NacH(   jp,kH)=naH(loc)+lp
!
!                 ...constraint coefficient
                     ConstrH(jp,kH)=RRRH(1,lp,3,1)*0.5d0
                  enddo
!
!              ...STEP 2.2 : loop through grandparent vertex nodes dof
                  do ivp=1,2
!
!                 ...grandparent vertex node
                     nodp = NEDG_CONS(ivp,nce)
!
!                 ...update number of parent dof's
                     NrconH(kH) = NrconH(kH)+1
!
!                 ...vertex node location on list of modified element's nodes
                     call locate(nodp,Nodm,Nrnodm, loc)
!
!                 ...advance parent dof local number
                     jp=jp+1
!
!                 ...parent dof location
                     NacH(jp,kH)=naH(loc)+1
!
!                 ...constraint coefficient
                     ConstrH(jp,kH)=0.25d0
                  enddo
               endif
            enddo
!
!
!-----------------------------------------------------------------------
!     H4-REFINED QUAD FACE
!-----------------------------------------------------------------------
!
!     -- MID-FACE NODE CONSTRAINED BY AN H4-REFINED FACE --
!
         case(21,22,23,24)
!
!        ...determine horizontal and vertical son number
            select case(icase)
            case(21) ; ish=1 ; isv=1
            case(22) ; ish=2 ; isv=1
            case(23) ; ish=2 ; isv=2
            case(24) ; ish=1 ; isv=2
            end select
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
!
!        ...parent node location on list of modified element's nodes
            call locate(nodp,Nodm,Nrnodm, loc)
!
!        ...number of parent node horizontal and vertical dof's
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh = nordh-1 ; ndofEh = nordh
            ndofHv = nordv-1 ; ndofEv = nordv
!
!           H1-loop through the constrained mid-face node dof's
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do lv=1,ndofHv
               do lh=1,ndofHh
!
!              ...advance dof counter
                  kH=kH+1
!
!              ...number of parent dof's
                  NrconH(kH) = ndofHv*ndofHh
!
!              ...initialize parent dof local number
                  jp=0
!
!              ...loop through the parent node dof's
                  do lpv=1,ndofHv
                     do lph=1,ndofHh
!
!                    ...advance parent dof local number
                        jp=jp+1
!
!                    ...parent dof location
                        NacH(   jp,kH)=naH(loc)+(lpv-1)*ndofHh+lph
!
!                    ...constraint coefficient
                        ConstrH(jp,kH)=RRRH(1,lph,ish,lh)*RRRH(1,lpv,isv,lv)
                     enddo
                  enddo
               enddo
            enddo
!
!           H(curl)-loop through the constrained mid-face node dof's
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ...STEP 1 : horizontal dof's first
            do lv=1,ndofHv
               do lh=1,ndofEh
                  kE=kE+1; jp=0
                  NrconE(kE) = ndofHv*ndofEh
                  do lpv=1,ndofHv
                     do lph=1,ndofEh
                        jp=jp+1
                        NacE(jp,kE)   =naE(loc)+(lpv-1)*ndofEh+lph
                        ConstrE(jp,kE)=RRRE(1,lph,ish,lh)*RRRH(1,lpv,isv,lv)
                     enddo
                  enddo
               enddo
            enddo
!
!        ...STEP 2 : vertical dof's next
            do lv=1,ndofEv
               do lh=1,ndofHh
                  kE=kE+1; jp=0
                  NrconE(kE) = ndofEv*ndofHh
                  do lpv=1,ndofEv
                     do lph=1,ndofHh
                        jp=jp+1
                        NacE(jp,kE)  = naE(loc) + ndofEh*ndofHv    &
                                     + (lpv-1)*ndofHh+lph
                        ConstrE(jp,kE)=RRRH(1,lph,ish,lh)*RRRE(1,lpv,isv,lv)
                     enddo
                  enddo
               enddo
            enddo
!
!           H(div)-loop through the constrained mid-face node dof's
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ...horizontal first then vertical
            do lv=1,ndofEv
               do lh=1,ndofEh
                  kV=kV+1; jp=0
                  NrconV(kV) = ndofEv*ndofEh
                  do lpv=1,ndofEv
                     do lph=1,ndofEh
                        jp=jp+1
                        NacV(jp,kV)   =naV(loc)+(lpv-1)*ndofEh+lph
                        ConstrV(jp,kV)=RRRE(1,lph,ish,lh)*RRRE(1,lpv,isv,lv)
                     enddo
                  enddo
               enddo
            enddo
!
!
!     -- HORIZONTAL MID-EDGE NODE CONSTRAINED BY AN H4-REFINED FACE --
!
         case(26,28)
!
!        ...determine horizontal son number
            select case(icase)
            case(28) ; is=1
            case(26) ; is=2
            endselect
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
!
!        ...location of parent node on list of modified element's nodes
            call locate(nodp,Nodm,Nrnodm, loc)
!
!        ...number of parent node horizontal and vertical dof's
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh=nordh-1 ; ndofEh=nordh
            ndofHv=nordv-1 ; ndofEv=nordv
!
!        ...parent mid-edge nodes (south,north)
            do ip=1,3,2
               nodp = iabs(NFACE_CONS(ip,nc))
               call locate(nodp,Nodm,Nrnodm, locp(ip))
               ndofHp(ip) = NODES_MG(nodp)%orderC-1
               ndofEp(ip) = NODES_MG(nodp)%orderC
            enddo
!
!           H1::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
               kH=kH+1
               NrconH(kH) = ndofHh*ndofHv + ndofHp(1) + ndofHp(3)
!
!           ...loop through the parent mid-face node dof
               jp=0
               do lpv=1,ndofHv
                  do lph=1,ndofHh
                     jp=jp+1
                     NacH(jp,kH)    = naH(loc)+(lpv-1)*ndofHh+lph
                     ConstrH(jp,kH) = RRRH(1,lph,is,l)*RRRH(1,lpv,3,1)
                  enddo
               enddo
!
!           ...loop through the parent mid-edge nodes dofs
               do ip=1,3,2
                  do lp=1,ndofHp(ip)
                     jp=jp+1
                     NacH(jp,kH) = naH(locp(ip))+lp
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrH(jp,kH)=.5d0*RRRH(1,lp,is,l)
                     else
                        ConstrH(jp,kH)=.5d0*(-1.d0)**(lp+1)*RRRH(1,lp,is,l)
                     endif
                  enddo
               enddo
            enddo
!
!           Hcurl::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
               kE=kE+1
               NrconE(kE) = ndofEh*ndofHv + ndofEp(1) + ndofEp(3)
!
!           ...loop through the parent mid-face node dof
               jp=0
               do lpv=1,ndofHv
                  do lph=1,ndofEh
                     jp=jp+1
                     NacE(jp,kE)    = naE(loc)+(lpv-1)*ndofEh+lph
                     ConstrE(jp,kE) = RRRE(1,lph,is,l)*RRRH(1,lpv,3,1)
                  enddo
               enddo
!
!           ...loop through the parent mid-edge nodes dof
               do ip=1,3,2
                  do lp=1,ndofEp(ip)
                     jp=jp+1
                     NacE(jp,kE) = naE(locp(ip))+lp
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrE(jp,kE) = .5d0*RRRE(1,lp,is,l)
                     else
                        ConstrE(jp,kE) = .5d0*RRRE(1,lp,is,l)*(-1.d0)**(lp)
                     endif
                  enddo
               enddo
            enddo
!
!
!     -- VERTICAL MID-EDGE NODE CONSTRAINED BY AN H4-REFINED FACE --
!
         case(25,27)
!
!        ...determine vertical son number
            select case(icase)
            case(25); is=1
            case(27); is=2
            endselect
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
!
!        ...parent mid-edge nodes (east,west)
            do ip=2,4,2
               nodp = iabs(NFACE_CONS(ip,nc))
               call locate(nodp,Nodm,Nrnodm, locp(ip))
               ndofHp(ip) = NODES_MG(nodp)%orderC-1
               ndofEp(ip) = NODES_MG(nodp)%orderC
            enddo
!
!           H1::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
               kH=kH+1
               NrconH(kH) = ndofHh*ndofHv + ndofHp(2) + ndofHp(4)
!
!           ...loop through the parent mid-face node dofs
               jp=0
               do lpv=1,ndofHv
                  do lph=1,ndofHh
                     jp=jp+1
                     NacH(jp,kH) = naH(loc)+(lpv-1)*ndofHh+lph
                     ConstrH(jp,kH) = RRRH(1,lph,3,1)*RRRH(1,lpv,is,l)
                  enddo
               enddo
!
!           ...loop through the parent mid-edge nodes dof
               do ip=2,4,2
                  do lp=1,ndofHp(ip)
                     jp=jp+1
                     NacH(jp,kH)    = naH(locp(ip))+lp
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrH(jp,kH)=.5d0*RRRH(1,lp,is,l)
                     else
                        ConstrH(jp,kH)=.5d0*(-1.d0)**(lp+1)*RRRH(1,lp,is,l)
                     endif
                  enddo
               enddo
            enddo
!
!           Hcurl::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
               kE=kE+1
               NrconE(kE) = ndofHh*ndofEv + ndofEp(2) + ndofEp(4)
!
!           ...loop through the parent mid-face node dof
               jp=0
               do lpv=1,ndofEv
                  do lph=1,ndofHh
                     jp=jp+1
                     NacE(jp,kE) = naE(loc) + ndofEh*ndofHv    &
                                 + (lpv-1)*ndofHh+lph
                     ConstrE(jp,kE) = RRRH(1,lph,3,1)*RRRE(1,lpv,is,l)
                  enddo
               enddo
!
!           ...loop through the parent mid-edge nodes dof
               do ip=2,4,2
                  do lp=1,ndofEp(ip)
                     jp=jp+1
                     NacE(jp,kE) = naE(locp(ip))+lp
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrE(jp,kE) = .5d0*RRRE(1,lp,is,l)
                     else
                        ConstrE(jp,kE) = .5d0*RRRE(1,lp,is,l)*(-1.d0)**(lp)
                     endif
                  enddo
               enddo
            enddo
!
!
!     -- VERTEX NODE CONSTRAINED BY AN H4-REFINED FACE --
!
         case(29)
!
!        ...parent mid-face node
         nodp = NFACEC(nc)
         call locate(nodp,Nodm,Nrnodm, loc)
         call decode(NODES_MG(nodp)%orderC, nordh,nordv)
         ndofHh = nordh-1
         ndofHv = nordv-1
!
!        ...parent mid-edge nodes (south,east,north,west)
            do ip=1,4
               nodp = iabs(NFACE_CONS(ip,nc))
               call locate(nodp,Nodm,Nrnodm, locp(ip))
               ndofHp(ip) = NODES_MG(nodp)%orderC-1
            enddo
!
!        ...only one dof
            kH=kH+1
            NrconH(kH) = ndofHh*ndofHv      &
                       + ndofHp(1)+ndofHp(2)+ndofHp(3)+ndofHp(4)+4
!
!        ...loop through the parent mid-face node dof
            jp=0
            do lpv=1,ndofHv
               do lph=1,ndofHh
                  jp=jp+1
                  NacH(jp,kH)    = naH(loc)+(lpv-1)*ndofHh+lph
                  ConstrH(jp,kH) = RRRH(1,lph,3,1)*RRRH(1,lpv,3,1)
               enddo
            enddo
!
!        ...loop through the parent mid-edge node dof
            do ip=1,4
               do lp=1,ndofHp(ip)
                  jp=jp+1
                  NacH(jp,kH)    = naH(locp(ip))+lp
                  if (NFACE_CONS(ip,nc).gt.0) then
                     ConstrH(jp,kH) = .5d0*RRRH(1,lp,3,1)
                  else
                     ConstrH(jp,kH) = .5d0*(-1.d0)**(lp+1)*RRRH(1,lp,3,1)
                  endif
               enddo
            enddo
!
!        ...loop through the parent vertex dof
            do ip=1,4
               nodp = NFACE_CONS(4+ip,nc)
               call locate(nodp,Nodm,Nrnodm, loc)
               jp=jp+1
               NacH(jp,kH)  = naH(loc)+1
               ConstrH(jp,kH) = .25d0
            enddo
!
!
!-----------------------------------------------------------------------
!        HORIZONTALLY H2-REFINED QUAD FACE
!-----------------------------------------------------------------------
!
!     -- MID-FACE NODE CONSTRAINED BY A HORIZONTALLY H2-REFINED FACE --
!
         case(31,32, 34,35, 61,62)
!
!        ...determine vertical son number
            select case(icase)
            case(31,34,61); is=1
            case(32,35,62); is=2
            end select
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
!
!           H1::loop through the constrained mid-face node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do lv=1,ndofHv
               do lh=1,ndofHh
                  kH=kH+1
                  NrconH(kH) = ndofHv
   !
   !              ...loop through the parent mid-face node dof
                  jp=0
                  do lpv=1,ndofHv
                     jp=jp+1
                     NacH(jp,kH) = naH(loc)+(lpv-1)*ndofHh+lh
                     ConstrH(jp,kH) = RRRH(1,lpv,is,lv)
                  enddo
               enddo
            enddo
!
!           Hcurl::loop through the constrained mid-face node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!        ...horizontal dof first
            do lv=1,ndofHv
               do lh=1,ndofEh
                  kE=kE+1; jp=0
                  NrconE(kE) = ndofHv
                  do lpv=1,ndofHv
                     jp=jp+1
                     NacE(jp,kE) = naE(loc)+(lpv-1)*ndofEh+lh
                     ConstrE(jp,kE) = RRRH(1,lpv,is,lv)
                  enddo
               enddo
            enddo
!
!        ...vertical dof next
            do lv=1,ndofEv
               do lh=1,ndofHh
                  kE=kE+1; jp=0
                  NrconE(kE) = ndofEv
                  do lpv=1,ndofEv
                     jp=jp+1
                     NacE(jp,kE) = naE(loc)+ndofEh*ndofHv+(lpv-1)*ndofHh+lh
                     ConstrE(jp,kE) = RRRE(1,lpv,is,lv)
                  enddo
               enddo
            enddo
!
!           Hdiv(L2)::loop through the constrained mid-face node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ...horizontal first then vertical
            do lv=1,ndofEv
               do lh=1,ndofEh
                  kV=kV+1; jp=0
                  NrconV(kV) = ndofEv
                  do lpv=1,ndofEv
                     jp=jp+1
                     NacV(jp,kV)   =naV(loc)+(lpv-1)*ndofEh+lh
                     ConstrV(jp,kV)=RRRE(1,lpv,is,lv)
                  enddo
               enddo
            enddo
!
!
!     -- HORIZONTAL MID-EDGE NODE CONSTRAINED BY A HORIZONTALLY H2-REFINED FACE --
!
         case(33,36,63)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
!
!        ...parent mid-edge nodes (south,north)
            do ip=1,3,2
               nodp = iabs(NFACE_CONS(ip,nc))
               select case(NODES_MG(nodp)%master)
               case(1)
                  ncflag(ip)=0
                  call locate(nodp,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = 1
                  ndofEp(ip) = 1
               case(0)
                  nfath = NODES(nodp)%father
                  call locate(nodp,NODES(nfath)%sons(1:2),2, ncflag(ip))
                  call locate(nfath,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = NODES_MG(nfath)%orderC-1
                  ndofEp(ip) = NODES_MG(nfath)%orderC
               end select
            enddo
!
!           H1::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
               kH=kH+1
               NrconH(kH) = ndofHv+ndofHp(1)+ndofHp(3)
!
!           ...loop through the parent mid-face node dof
               jp=0
               do lpv=1,ndofHv
                  jp=jp+1
                  NacH(jp,kH)    = naH(loc)+(lpv-1)*ndofHh+l
                  ConstrH(jp,kH) = RRRH(1,lpv,3,1)
               enddo
!
!           ...loop through the parent mid-edge nodes
               do ip=1,3,2
                  select case(ncflag(ip))
!
!              ...unconstrained mid-edge node, only one parent dof
                  case(0)
                     jp=jp+1
                     NacH(jp,kH) = naH(locp(ip))+l
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrH(jp,kH) = .5d0
                     else
                        ConstrH(jp,kH) = .5d0*(-1.d0)**(l-1)
                     endif
!
!              ...constrained mid-edge node
                  case(1,2)
                     do lp=1,ndofHp(ip)
                        jp=jp+1
                        NacH(jp,kH) = naH(locp(ip))+lp
                        if (NFACE_CONS(ip,nc).gt.0) then
                           ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)
                        else
                           ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)  &
                                          *(-1.d0)**(l-1)
                        endif
                     enddo
                  end select
               enddo
            enddo
!
!           Hcurl::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
               kE=kE+1; jp=0
               NrconE(kE) = ndofHv+ndofEp(1)+ndofEp(3)
!
!           ...loop through the parent mid-face node dof
               do lpv=1,ndofHv
                  jp=jp+1
                  NacE(jp,kE)    = naE(loc)+(lpv-1)*ndofEh+l
                  ConstrE(jp,kE) = RRRH(1,lpv,3,1)
               enddo
!
!           ...loop through the parent mid-edge nodes
               do ip=1,3,2
                  select case(ncflag(ip))
!
!              ...unconstrained mid-edge node, only one parent dof
                  case(0)
                     jp=jp+1
                     NacE(jp,kE) = naE(locp(ip))+l
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrE(jp,kE) = .5d0
                     else
                        ConstrE(jp,kE) = .5d0*(-1.d0)**(l)
                     endif
!
!              ...constrained mid-edge node
                  case(1,2)
                     do lp=1,ndofEp(ip)
                        jp=jp+1
                        NacE(jp,kE) = naE(locp(ip))+lp
                        if (NFACE_CONS(ip,nc).gt.0) then
                           ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)
                        else
                           ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)   &
                                          *(-1.d0)**(l)
                        endif
                     enddo
                  end select
               enddo
            enddo
!
!
!-----------------------------------------------------------------------
!     VERTICALLY H2-REFINED QUAD FACE
!-----------------------------------------------------------------------
!
!     -- MID-FACE NODE CONSTRAINED BY A VERTICALLY H2-REFINED FACE --
!
         case(41,42, 44,45, 51,52)
!
!        ...determine vertical son number
            select case(icase)
            case(41,44,51); is=1
            case(42,45,52); is=2
            end select
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
!
!        ...H1:loop through the constrained mid-face node dof
            do lv=1,ndofHv
               do lh=1,ndofHh
                  kH=kH+1
                  NrconH(kH) = ndofHh
!
!              ...loop through the parent mid-face node dof
                  jp=0
                  do lph=1,ndofHh
                     jp=jp+1
                     NacH(jp,kH) = naH(loc)+(lv-1)*ndofHh+lph
                     ConstrH(jp,kH) = RRRH(1,lph,is,lh)
                  enddo
               enddo
            enddo
!
!           Hcurl::loop through the constrained mid-face node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!        ...horizontal dof first
            do lv=1,ndofHv
               do lh=1,ndofEh
                  kE=kE+1; jp=0
                  NrconE(kE) = ndofEh
                  do lph=1,ndofEh
                     jp=jp+1
                     NacE(jp,kE) = naE(loc)+(lv-1)*ndofEh+lph
                     ConstrE(jp,kE) = RRRE(1,lph,is,lh)
                  enddo
               enddo
            enddo
!
!        ...vertical dof next
            do lv=1,ndofEv
               do lh=1,ndofHh
                  kE=kE+1; jp=0
                  NrconE(kE) = ndofHh
                  do lph=1,ndofHh
                     jp=jp+1
                     NacE(jp,kE) = naE(loc)+ndofEh*ndofHv+(lv-1)*ndofHh+lph
                     ConstrE(jp,kE) = RRRH(1,lph,is,lh)
                  enddo
               enddo
            enddo
!
!           Hdiv(L2)::loop through the constrained mid-face node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ...horizontal first then vertical
            do lv=1,ndofEv
               do lh=1,ndofEh
                  kV=kV+1; jp=0
                  NrconV(kV) = ndofEh
                  do lph=1,ndofEh
                     jp=jp+1
                     NacV(jp,kV)   =naV(loc)+(lv-1)*ndofEh+lph
                     ConstrV(jp,kV)=RRRE(1,lph,is,lh)
                  enddo
               enddo
            enddo
!
!
!     -- VERTICAL MID-EDGE NODE CONSTRAINED BY A VERTICALLY H2-REFINED FACE --
!
         case(43,46,53)
!
!        ...parent mid-face node
            nodp = NFACEC(nc)
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES_MG(nodp)%orderC, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
!
!        ...parent mid-edge nodes (east,west)
            do ip=2,4,2
               nodp = iabs(NFACE_CONS(ip,nc))
               select case(NODES_MG(nodp)%master)
               case(1)
                  ncflag(ip)=0
                  call locate(nodp,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = 1
                  ndofEp(ip) = 1
               case(0)
                  nfath = NODES(nodp)%father
                  call locate(nodp,NODES(nfath)%sons(1:2),2, ncflag(ip))
                  call locate(nfath,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = NODES_MG(nfath)%orderC-1
                  ndofEp(ip) = NODES_MG(nfath)%orderC
               end select
            enddo
!
!           ...H1::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
               kH=kH+1
               NrconH(kH) = ndofHh+ndofHp(2)+ndofHp(4)
!
!           ...loop through the parent mid-face node dof
               jp=0
               do lph=1,ndofHh
                  jp=jp+1
                  NacH(jp,kH)    = naH(loc)+(l-1)*ndofHh+lph
                  ConstrH(jp,kH) = RRRH(1,lph,3,1)
               enddo
!
!           ...loop through the parent mid-edge nodes
               do ip=2,4,2
                  select case(ncflag(ip))
!
!              ...unconstrained mid-edge node, only one parent dof
                  case(0)
                     jp=jp+1
                     NacH(jp,kH)    = naH(locp(ip))+l
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrH(jp,kH) = .5d0
                     else
                        ConstrH(jp,kH) = .5d0*(-1.d0)**(l-1)
                     endif
!
!              ...constrained mid-edge node
                  case(1,2)
                     do lp=1,ndofHp(ip)
                        jp=jp+1
                        NacH(jp,kH) = naH(locp(ip))+lp
                        if (NFACE_CONS(ip,nc).gt.0) then
                           ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)
                        else
                           ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)   &
                                          *(-1.d0)**(l-1)
                        endif
                     enddo
                  end select
               enddo
            enddo
!
!           Hcurl::loop through the constrained mid-edge node dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
               kE=kE+1; jp=0
               NrconE(kE) = ndofHh+ndofEp(2)+ndofEp(4)
!
!           ...loop through the parent mid-face node dof
               do lph=1,ndofHh
                  jp=jp+1
                  NacE(jp,kE)  = naE(loc)+ndofEh*ndofHv  &
                               + (l-1)*ndofHh+lph
                  ConstrE(jp,kE) = RRRH(1,lph,3,1)
               enddo
!
!           ...loop through the parent mid-edge nodes
               do ip=2,4,2
                  select case(ncflag(ip))
!
!              ...unconstrained mid-edge node, only one parent dof
                  case(0)
                     jp=jp+1
                     NacE(jp,kE) = naE(locp(ip))+l
                     if (NFACE_CONS(ip,nc).gt.0) then
                        ConstrE(jp,kE) = .5d0
                     else
                        ConstrE(jp,kE) = .5d0*(-1.d0)**(l)
                     endif
!
!               ...constrained mid-edge node
                   case(1,2)
                     do lp=1,ndofEp(ip)
                        jp=jp+1
                        NacE(jp,kE) = naE(locp(ip))+lp
                        if (NFACE_CONS(ip,nc).gt.0) then
                           ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)
                        else
                           ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)  &
                                          *(-1.d0)**(l)
                        endif
                     enddo
                  end select
               enddo
            enddo
!
!
!-----------------------------------------------------------------------
!     H4-REFINED TRIANGULAR FACE
!-----------------------------------------------------------------------
!
!     -- MID-FACE NODE CONSTRAINED BY AN H4-REFINED FACE --
!
         case(71,72,73,74)
            call decode(icase, nvoid,is)
!
!           H1::loop through the nodal dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
               kH=kH+1
!
!           ...initiate parent dof counter
               jp=0
!
!           ...parent mid-face node
               nodp = NFACEC(nc)
               call locate(nodp,Nodm,Nrnodm, loc)
!
!           ...loop through parent node dof
               do lp=1,NdofmH(loc)
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+lp
                  ConstrH(jp,kH) = RRTH(1,lp,is,l)
               enddo
!
!           ...loop through parent medg nodes
               do ie=1,3
                  nodp = iabs(NFACE_CONS(ie,nc))
                  call locate(nodp,Nodm,Nrnodm, loc)
!
!              ...loop through parent node dof
                  do lp=1,NdofmH(loc)
                     jp=jp+1
                     NacH(jp,kH) = naH(loc)+lp
                     if (NFACE_CONS(ie,nc).gt.0) then
                        ConstrH(jp,kH) = RRTH(1+ie,lp,is,l)
                     else
                        ConstrH(jp,kH) = RRTH(1+ie,lp,is,l)*(-1.d0)**(lp+1)
                     endif
                  enddo
               enddo
!
!           ...save the number of constraining dof
               NrconH(kH) = jp
            enddo
!
!           Hcurl::loop through the nodal dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
               kE=kE+1
!
!           ...initiate parent dof counter
               jp=0
!
!           ...parent mid-face node
               nodp = NFACEC(nc)
               call locate(nodp,Nodm,Nrnodm, loc)
!
!           ...loop through parent node dof
               do lp=1,NdofmE(loc)
                  jp=jp+1
                  NacE(jp,kE) = naE(loc)+lp
                  ConstrE(jp,kE) = RRTE(1,lp,is,l)
               enddo
!
!           ...loop through parent medg nodes
               do ie=1,3
                  nodp = iabs(NFACE_CONS(ie,nc))
                  call locate(nodp,Nodm,Nrnodm, loc)
!
!              ...loop through parent node dof
                  do lp=1,NdofmE(loc)
                     jp=jp+1
                     NacE(jp,kE) = naE(loc)+lp
                     if (NFACE_CONS(ie,nc).gt.0) then
                        ConstrE(jp,kE) = RRTE(1+ie,lp,is,l)
                     else
                        ConstrE(jp,kE) = RRTE(1+ie,lp,is,l)*(-1.d0)**(lp)
                     endif
                  enddo
               enddo
!
!           ...save the number of constraining dof
               NrconE(kE) = jp
            enddo
!
!           H(div)::loop through the nodal dof
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofV
               kV=kV+1
!
!           ...initiate parent dof counter
               jp=0
!
!           ...parent mid-face node
               nodp = NFACEC(nc)
               call locate(nodp,Nodm,Nrnodm, loc)
!
!           ...loop through parent node dof
               do lp=1,NdofmV(loc)
                  jp=jp+1
                  NacV(jp,kV) = naV(loc)+lp
                  ConstrV(jp,kV) = RRTQ(lp,is,l)
               enddo
!
!           ...save the number of constraining dof
               NrconV(kV) = jp
            enddo
!
!
!     -- MID-EDGE NODE CONSTRAINED BY AN H4-REFINED FACE --
!
         case(75,76,77)
            call decode(icase, nvoid,is)
!
!        ...H1:loop through the nodal dof
            do l=1,ndofH
               kH=kH+1
!
!           ...initiate parent dof counter
               jp=0
!
!           ...parent mid-face node
               nodp = NFACEC(nc)
               call locate(nodp,Nodm,Nrnodm, loc)
!
!           ...loop through parent mdlt node dof
               do lp=1,NdofmH(loc)
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+lp
                  ConstrH(jp,kH) = RRTH(1,lp,is,l)
               enddo
!
!           ...loop through parent medg nodes
               do ie=1,3
                  nodp = iabs(NFACE_CONS(ie,nc))
                  call locate(nodp,Nodm,Nrnodm, loc)
                  if (iprint.eq.1) then
                     write(*,*) 'logic: ie,nodp,loc = ',ie,nodp,loc
                   endif
!
!              ...loop through parent medg node dof
                  do lp=1,NdofmH(loc)
                     jp=jp+1
                     NacH(jp,kH) = naH(loc)+lp
                     if (NFACE_CONS(ie,nc).gt.0) then
                        ConstrH(jp,kH) = RRTH(1+ie,lp,is,l)
                     else
                        ConstrH(jp,kH) = RRTH(1+ie,lp,is,l)*(-1.d0)**(lp+1)
                     endif
                  enddo
               enddo
               if (iprint.eq.1) then
                   write(*,*) 'l = ',l
                   write(*,7013) NacH(1:jp,kH)
                   write(*,7014) ConstrH(1:jp,kH)
               endif
!
!           ...save the number of constraining dof
               NrconH(kH) = jp
            enddo
!
!        ...Hcurl::loop through the nodal dof
            do l=1,ndofE
               kE=kE+1
!
!           ...initiate parent dof counter
               jp=0
!
!           ...parent mid-face node
               nodp = NFACEC(nc)
               call locate(nodp,Nodm,Nrnodm, loc)
!
!           ...loop through parent mdlt node dof
               do lp=1,NdofmE(loc)
                  jp=jp+1
                  NacE(jp,kE) = naE(loc)+lp
                  ConstrE(jp,kE) = RRTE(1,lp,is,l)
               enddo
!
!           ...loop through parent medg nodes
               do ie=1,3
                  nodp = iabs(NFACE_CONS(ie,nc))
                  call locate(nodp,Nodm,Nrnodm, loc)
                  if (iprint.eq.1) then
                     write(*,*) 'logic: ie,nodp,loc = ',ie,nodp,loc
                  endif
!
!              ...loop through parent medg node dof
                  do lp=1,NdofmE(loc)
                     jp=jp+1
                     NacE(jp,kE) = naE(loc)+lp
                     if (NFACE_CONS(ie,nc).gt.0) then
                        ConstrE(jp,kE) = RRTE(1+ie,lp,is,l)
                     else
                        ConstrE(jp,kE) = RRTE(1+ie,lp,is,l)*(-1.d0)**(lp)
                     endif
                  enddo
               enddo
               if ((iprint.eq.1).and.(Idec.gt.1)) then
                  write(*,*) '75,76,77::Mdle,l,jp = ',Mdle,l,jp
                  write(*,7013) NacE(1:jp,kE)
                  write(*,7014) ConstrE(1:jp,kE)
               endif
!
!           ...save the number of constraining dof
               NrconE(kE) = jp
            enddo
!
!
!-----------------------------------------------------------------------
!  H2-REFINED TRIANGULAR FACE (QUAD + TRIANGLE)
!-----------------------------------------------------------------------
!
!     -- MDLQ NODE CONSTRAINED BY AN H2-REFINED TRIANGULAR FACE --
!
         case(82,83,84)
            call decode(icase, nvoid,iref)
!
!        ...loop through the nodal dof
            do l=1,ndofH
               kH=kH+1
!
!           ...initiate parent dof counter
               jp=0
!
!           ...parent mdlt node
               nodp = NFACEC(nc)
               call locate(nodp,Nodm,Nrnodm, loc)
!
!           ...loop through parent mdlt node dof
               do lp=1,NdofmH(loc)
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+lp
                  np = NODES_MG(nodp)%orderC
                  ConstrH(jp,kH) = get_rrqh(iref,np,1,lp,l)
               enddo
!
!           ...loop through parent medg nodes
               do ie=1,3
                  nodp = iabs(NFACE_CONS(ie,nc))
                  call locate(nodp,Nodm,Nrnodm, loc)
                  if (iprint.eq.1) then
                     write(*,*) 'logic: ie,nodp,loc = ',ie,nodp,loc
                  endif
!
!              ...loop through parent medg node dof
                  do lp=1,NdofmH(loc)
                     jp=jp+1
                     NacH(jp,kH) = naH(loc)+lp
                     if (NFACE_CONS(ie,nc).gt.0) then
                        ConstrH(jp,kH) = get_rrqh(iref,np,1+ie,lp,l)
                     else
                        ConstrH(jp,kH) = get_rrqh(iref,np,1+ie,lp,l)   &
                                       *(-1.d0)**(lp+1)
                     endif
                  enddo
               enddo
               if (iprint.eq.1) then
                  write(*,*) 'l = ',l
                  write(*,7013) NacH(1:jp,kH)
                  write(*,7014) ConstrH(1:jp,kH)
               endif
!
!           ...save the number of constraining dof
               NrconH(kH) = jp
!
            enddo
!         
!-----------------------------------------------------------------------
!
         end select
!
!     ...if a constrained node
      endif
!
!  ...printing
      if (iprint.eq.1) then
         write(*,9000) i,nod,kH,kE,kV
 9000    format(' i = ',i2,' ; nod = ',i5,' ; kH = ',i3,    &
                ' ; kE = ',i3,' ; kV = ',i3 )
      endif
!
! ...end of loop through local nodes
   enddo
!
!..sanity check
   kEtot = kE
   do kE=1,kEtot
      do kp=1,NrconE(kE)
         if ((NacE(kp,kE).le.0).or.(NacE(kp,kE).gt.icE)) then
            write(*,7200) Mdle,kE,kp,NacE(kp,kE)
 7200       format('logic: Mdle,kE,kp,NacE(kp,kE) = ',4i5)
            stop 1
         endif
      enddo
   enddo
!
 100  continue
!
   if ((iprint.ge.1).and.(Idec.ne.1)) then
      write(*,7011) Mdle
!
      write(*,*) ' -- H1 --'
      nrdofH_loc=kH
      do kH=1,nrdofH_loc
         write(*,7012) kH,NrconH(kH)
         write(*,7013) NacH(1:NrconH(kH),kH)
         write(*,7014) ConstrH(1:NrconH(kH),kH)
      enddo
!
      write(*,*) ' '
      write(*,*) ' -- H(curl) --'
      nrdofE_loc=kE
      do kE=1,nrdofE_loc
         write(*,7012) kE,NrconE(kE)
         write(*,7013) NacE(1:NrconE(kE),kE)
         write(*,7014) ConstrE(1:NrconE(kE),kE)
      enddo
!
      write(*,*) ' -- H(div) --'
      nrdofV_loc=kV
      do kV=1,nrdofV_loc
         write(*,7012) kV,NrconV(kV)
         write(*,7013) NacV(1:NrconV(kV),kV)
         write(*,7014) ConstrV(1:NrconV(kV),kV)
      enddo
      call pause
!
 7011 format(' logic: Mdle = ',i6)
 7012 format(' logic: k = ',i4,' Nrcon = ',i3)
 7013 format(20i6)
 7014 format(20f6.3)

   endif
!
!
   end subroutine logicC
