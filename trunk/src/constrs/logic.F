c----------------------------------------------------------------------
c
c   routine name       - logic
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine establishes list of nodes
c                        for a modified element, and the corresponding
c                        constraints coefficients
c   remark: this routine must be OMP thread-safe
c
c   arguments :
c     in:
c            Mdle      - an element number, same as the middle node
c                        number
c            Idec      = 1 if only the element nodes are to be
c                          determined
c     out:
c            Nodm      - actual (unconstrained) nodes returned in the
c                        standard order: vertex, mid-edge, mid-face
c                        and middle nodes
c            NdofmH    - the corresponding number of H1 dof
c            NdofmE    - the corresponding number of H(curl) dof
c            NdofmV    - the corresponding number of H(div) dof
c            Nrnodm    - number of modified element nodes
c
c            NrconH    - number of parent dof's for an i-th local H1 dof
c            NacH(j,i),j=1,...,NrconH(i) - list of parent dof's
c                        for i-th local dof
c            ConstrH(j,i),j=1,...,nrcon(i) - the corresponding
c                        constraints coefficients
c
c            NrconE,NacE,ConstrE - same for H(curl) dof
c            NrconV,NacV,ConstrV - same for H(div) dof
c
c----------------------------------------------------------------------
c
      subroutine logic(Mdle,Idec,
     .                 Nodm,NdofmH,NdofmE,NdofmV,Nrnodm,
     .                 NrconH,NacH,ConstrH,
     .                 NrconE,NacE,ConstrE,
     .                 NrconV,NacV,ConstrV)
c
      use element_data
      use constraints
      use data_structure3D
      use constrained_nodes
      use refinements
c
      implicit none
c
      integer, intent(in)  :: Mdle,Idec
c
      integer, intent(out) :: Nodm(MAXNODM)
c
      integer, intent(out) :: NdofmH(MAXNODM),
     .                        NdofmE(MAXNODM),
     .                        NdofmV(MAXNODM)
c
      integer, intent(out) :: Nrnodm
c
      integer, intent(out) :: NrconH(MAXbrickH),
     .                        NrconE(MAXbrickE),
     .                        NrconV(MAXbrickV)
c
      integer, intent(out) :: NacH(NACDIM,MAXbrickH),
     .                        NacE(NACDIM,MAXbrickE),
     .                        NacV(NACDIM,MAXbrickV)
c
      real(8), intent(out) :: ConstrH(NACDIM,MAXbrickH),
     .                        ConstrE(NACDIM,MAXbrickE),
     .                        ConstrV(NACDIM,MAXbrickV)
c
c  ...dummy argument for ndof_nod
      integer :: NdofmQ(MAXNODM)
c
c  ...middle node type
      integer :: ntype
c
 7002     format('im = ',i3,' nod = ',i6,' TYPE = ',a6,' ORDER = ',i3,
     .      ' naH = ',i3,' NdofmH = ',i3,
     .      ' naE = ',i3,' NdofmE = ',i3,
     .      ' naV = ',i3,' NdofmV = ',i3)
 7030       format('logic: ACTIVE NODE      i = ',i3,' nod = ',i6,
     .             ' TYPE = ',a6,' ORDER = ',i2)
 7031       format('logic: CONSTRAINED NODE i = ',i3,' nod = ',i6,
     .             ' TYPE = ',a6,' ORDER = ',i2,' icase = ',i2,
     .             ' EDGE NODES = ',i5,2x,2i5)
 7032       format('logic: CONSTRAINED NODE i = ',i3,' nod = ',i6,
     .             ' TYPE = ',a6,' ORDER = ',i2,' icase = ',i2,
     .             ' FACE NODES = ',i5,2x,4i5,2x,4i5)
c
c  ...element local nodes, orientations
      integer :: nodesl(27),norientl(27)
c
c  ...offsets for first dof for modified element nodes
      integer :: naH(MAXNODM),naE(MAXNODM),naV(MAXNODM)
c
c  ...order for the element
      integer :: norder(19)
c
c  ...addresses and number of dof for parent mid-edge and vertex nodes
      integer :: locp(8),ndofHp(8),ndofEp(8)
c
c  ...unconstrained/constrained edge flags
c     = 0          if edge is unconstrained (active)
c     = son number if edge is constrained (inactive)
      integer :: ncflag(4)
c
c  ...H(div) constrained approximation flag
      logical :: CONSTRAINED_HDIV = .TRUE.
c
c  ...miscellanea
      integer :: i,ie,ip,iref,is,ish,isv,iv,ivp,jp,k,kp
      integer :: icH,icE,icV,im,kH,kE,kV,j,icase
      integer :: kEtot,l,lh,loc,lp,lph,lpv,lv
      integer :: nc,nce,nfath,nod,nodp,nordh,nordv,np,nrnodes,nvoid
      integer :: ndofE,ndofEh,ndofEv
      integer :: ndofH,ndofHh,ndofHv
      integer :: ndofV,ndofQ,ndofHpp
      integer :: nrdofH,nrdofE,nrdofV,nrdofQ
c
#if HP3D_DEBUG
      integer :: nrdofH_loc,nrdofE_loc,nrdofV_loc
      integer :: iprint
      iprint=0
#endif
c
c----------------------------------------------------------------------
c
      Nodm   = 0
      NdofmH = 0; NrconH = 0; NacH = 0; ConstrH = 0.d0
      NdofmE = 0; NrconE = 0; NacE = 0; ConstrE = 0.d0
      NdofmV = 0; NrconV = 0; NacV = 0; ConstrV = 0.d0
c
c     initialize constraints arrays
      if (.not. INITIALIZED_CONSTR) call init_cnstr
c
#if HP3D_DEBUG
      if (iprint.ge.1) then
        write(*,7000) Mdle
 7000   format(' -- logic: Mdle = ',i6,' --')
      endif
#endif
c
c     get element nodes and build the local data base on constrained
c     nodes for the element
      call get_connect_info(Mdle, nodesl,norientl)
c
#if HP3D_DEBUG
      if (iprint.eq.1) then
        call elem_show(Mdle, NODES(Mdle)%ntype,nodesl,norientl)
      endif
#endif
c
c     get nodes of the modified element
      call logic_nodes(Mdle,nodesl, Nodm,Nrnodm)
cc
cccc     store location of the node on the list of modified nodes
ccc      do im=1,Nrnodm
ccc        NODES(Nodm(im))%visit=im
ccc      enddo
c
c
c     establish offsets for the modified element nodal dof
c
c     initialize
      icH=0 ; icE=0 ; icV=0
c
c     loop over nodes of modified element
      do im=1,Nrnodm
c
c       number of H1, H(curl), H(div) dofs stored so far
        naH(im)=icH ; naE(im)=icE ; naV(im)=icV
c
c       number of H1, H(curl), H(div) dofs associated to node
        nod=Nodm(im)
        call ndof_nod(NODES(nod)%ntype,NODES(nod)%order, NdofmH(im),
     .                                                   NdofmE(im),
     .                                                   NdofmV(im),
     .                                                   NdofmQ(im))
c
c       update H1,H(curl),H(div) dofs counters
        icH=icH+NdofmH(im) ; icE=icE+NdofmE(im) ; icV=icV+NdofmV(im)
      enddo
c
#if HP3D_DEBUG
c     printing
      if (iprint.eq.1) then
        write(*,*) 'MODIFIED ELEMENT NODES:'
        do im=1,Nrnodm
          nod=Nodm(im)
          write(*,7002) im,nod,S_Type(NODES(nod)%ntype),
     .                  NODES(nod)%order,
     .                  naH(im),NdofmH(im),
     .                  naE(im),NdofmE(im),
     .                  naV(im),NdofmV(im)
        enddo
      endif
#endif
c
c----------------------------------------------------------------------
      if (Idec.eq.1) go to 100
c
c     initialize
      NrconH=0 ; NacH=0 ; ConstrH=0.d0
      NrconE=0 ; NacE=0 ; ConstrE=0.d0
      NrconV=0 ; NacV=0 ; ConstrV=0.d0
c
c     construct the constrained approximation info
c
c     initialize local H1, H(curl), H(div) dof counters
      kH=0 ; kE=0 ; kV=0
c
c     number of element's nodes
      ntype=NODES(Mdle)%ntype
      nrnodes=nvert(ntype)+nedge(ntype)+nface(ntype)+1
c
c     loop over element nodes
      do i=1,nrnodes
        nod=nodesl(i)
c
c======================================================================
c       NODE IS UNCONSTRAINED                                         |
c======================================================================
C         if (Is_active(nod)) then
        if ((Is_active(nod)) .or.
     .      (Is_inactive(nod) .and. NODES(nod)%ref_kind /= 0)) then
c
c         node location on the list of modified element's nodes
ccc       loc=NODES(nod)%visit
          call locate(nod,Nodm,Nrnodm, loc)
c
c         H1-loop through the nodal dof
          do j=1,NdofmH(loc)
c
c           advance dof counter
            kH=kH+1
c
c           number of parent dof's
            NrconH(   kH)=1
c
c           parent dof location
            NacH(   1,kH)=naH(loc)+j
c
c           constraint coefficient
            ConstrH(1,kH)=1.d0
          enddo
c
c         H(curl)-loop through the nodal dof
          do j=1,NdofmE(loc)
            kE=kE+1
            NrconE(   kE)=1
            NacE(   1,kE)=naE(loc)+j
            ConstrE(1,kE)=1.d0
          enddo
c
c         H(div)-loop through the nodal dof
          do j=1,NdofmV(loc)
            kV=kV+1
            NrconV(   kV)=1
            NacV(   1,kV)=naV(loc)+j
            ConstrV(1,kV)=1.d0
          enddo
c
#if HP3D_DEBUG
c         printing
          if (iprint.eq.1) then
            write(*,7030) i,nod,S_Type(NODES(nod)%ntype),
     .                    NODES(nod)%order
          endif
#endif
c
c======================================================================
c       NODE IS CONSTRAINED                                           |
c======================================================================
        else
c
c         compute number of dof for the node
          call ndof_nod(NODES(nod)%ntype,NODES(nod)%order,
     .                  ndofH,ndofE,ndofV,ndofQ)
c
c         identify the node case
          call decode2(NODES_CONSTR(i), nc,icase)
c
#if HP3D_DEBUG
c         printing
          if (iprint.eq.1) then
            select case(icase)
c
c           node constrained by an edge
            case(11,12,13,37,38,39,47,48,49)
            write(*,7031) i,nod,NODES(nod)%ntype,NODES(nod)%order,icase,
     .                    NEDGC(nc),NEDG_CONS(1:2,nc)
c
c           node constrained by a face
            case default
            write(*,7032) i,nod,NODES(nod)%ntype,NODES(nod)%order,icase,
     .                    NFACEC(nc),NFACE_CONS(1:8,nc)
            end select
          endif
#endif
c
c
          select case(icase)
c
c
c-----------------------------------------------------------------------
c  H2-REFINED EDGE
c-----------------------------------------------------------------------
c
c  -- FIRST AND SECOND MID-EDGE NODE CONSTRAINED BY AN EDGE --
c
          case(11,12, 37,38, 47,48)
c
c           local son number
            select case(icase)
            case(11,37,47) ; is=1
            case(12,38,48) ; is=2
            endselect
c
c           parent mid-edge node
            nodp = NEDGC(nc)
c
c           parent node location on the list of modified element's nodes
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
c
c           H1-loop through nodal dof's
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
c
c             advance dof counter
              kH=kH+1
c
c             number of parent dof's
              NrconH(kH) = ndofH
c
c             loop through parent node dof's
              do jp=1,NdofmH(loc)
c
c               parent dof location
                NacH(   jp,kH)=naH(loc)+jp
c
c               constraint coefficient
                ConstrH(jp,kH)=RRRH(1,jp,is,l)
              enddo
            enddo
c
c           H(curl)-loop through nodal dof's
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
              kE=kE+1
              NrconE(kE) = ndofE
c
c             loop through parent node dof's
              do jp=1,NdofmE(loc)
                NacE(   jp,kE)=naE(loc)+jp
                ConstrE(jp,kE)=RRRE(1,jp,is,l)
              enddo
            enddo
c
c           NO H(div) constraints across edges
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
c  -- VERTEX NODE CONSTRAINED BY AN EDGE --
c
          case(13,39,49)
c
c           H1 constraint ONLY
c           ~~~~~~~~~~~~~~~~~~
c
c           advance dof counter
            kH=kH+1
c
c           no parent dof's
            NrconH(kH)=0
c
c           parent mid-edge node
            nodp = NEDGC(nc)
c
c           number of parent node dof's
            ndofHpp = NODES(nodp)%order-1
c
c           number of parent dof's
            NrconH(kH) = NrconH(kH)+ndofHpp
c
c           parent node location on the list of modified element's nodes
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
c
c           initialize parent dof local number
            jp=0
c
c           STEP 1 : loop through parent mid-edge node dofs
            do lp=1,ndofHpp
c
c             advance parent dof local number
              jp=jp+1
c
c             parent dof location
              NacH(   jp,kH)=naH(loc)+lp
c
c             constraint coefficient
              ConstrH(jp,kH)=RRRH(1,lp,3,1)
            enddo
c
c           STEP 2 : loop through parent edge vertex dofs
            do iv=1,2
c
c             parent vertex node
              nodp = NEDG_CONS(iv,nc)
c
c             active (unconstrained) vertex node
              if (nodp.gt.0) then
c
c               update number of parent dof's
                NrconH(kH) = NrconH(kH)+1
c
c               vertex node location on list of modified element's nodes
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
c
c               advance parent dof local number
                jp=jp+1
c
c               parent dof location
                NacH(   jp,kH)=naH(loc)+1
c
c               constraint coefficient
                ConstrH(jp,kH)=RRRH(1+iv,1,3,1)
c
c             inactive (constrained) vertex node
              else
c
c               grandparent edge node
                call decode2(-nodp, nce,nvoid)
                nodp = NEDGC(nce)
c
c               number of grandparent node dof's
                ndofHpp = NODES(nodp)%order-1
c
c               update number of parent dof's
                NrconH(kH) = NrconH(kH)+ndofHpp
c
c               grandparent node location on list of modified element's nodes
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
c
c               STEP 2.1 : loop through grandparent mid-edge node dofs
                do lp=1,ndofHpp
c
c                 advance parent dof local number
                  jp=jp+1
c
c                 parent dof location
                  NacH(   jp,kH)=naH(loc)+lp
c
c                 constraint coefficient
                  ConstrH(jp,kH)=RRRH(1,lp,3,1)*0.5d0
                enddo
c
c               STEP 2.2 : loop through grandparent vertex nodes dof
                do ivp=1,2
c
c                 grandparent vertex node
                  nodp = NEDG_CONS(ivp,nce)
c
c                 update number of parent dof's
                  NrconH(kH) = NrconH(kH)+1
c
c                 vertex node location on list of modified element's nodes
ccc               loc = NODES(nodp)%visit
                  call locate(nodp,Nodm,Nrnodm, loc)
c
c                 advance parent dof local number
                  jp=jp+1
c
c                 parent dof location
                  NacH(   jp,kH)=naH(loc)+1
c
c                 constraint coefficient
                  ConstrH(jp,kH)=0.25d0
                enddo
              endif
            enddo
c
c
c-----------------------------------------------------------------------
c  H4-REFINED QUAD FACE
c-----------------------------------------------------------------------
c
c  -- MID-FACE NODE CONSTRAINED BY AN H4-REFINED FACE --
c
          case(21,22,23,24)
c
c           determine horizontal and vertical son number
            select case(icase)
            case(21) ; ish=1 ; isv=1
            case(22) ; ish=2 ; isv=1
            case(23) ; ish=2 ; isv=2
            case(24) ; ish=1 ; isv=2
            end select
c
c           parent mid-face node
            nodp = NFACEC(nc)
c
c           parent node location on list of modified element's nodes
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
c
c           number of parent node horizontal and vertical dof's
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1 ; ndofEh = nordh
            ndofHv = nordv-1 ; ndofEv = nordv
c
c           H1-loop through the constrained mid-face node dof's
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do lv=1,ndofHv
              do lh=1,ndofHh
c
c               advance dof counter
                kH=kH+1
c
c               number of parent dof's
                NrconH(kH) = ndofHv*ndofHh
c
c               initialize parent dof local number
                jp=0
c
c               loop through the parent node dof's
                do lpv=1,ndofHv
                  do lph=1,ndofHh
c
c                   advance parent dof local number
                    jp=jp+1
c
c                   parent dof location
                    NacH(   jp,kH)=naH(loc)+(lpv-1)*ndofHh+lph
c
c                   constraint coefficient
                    ConstrH(jp,kH)=RRRH(1,lph,ish,lh)*RRRH(1,lpv,isv,lv)
                  enddo
                enddo
              enddo
            enddo
c
c           H(curl)-loop through the constrained mid-face node dof's
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           STEP 1 : horizontal dof's first
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
c
c           STEP 2 : vertical dof's next
            do lv=1,ndofEv
              do lh=1,ndofHh
                kE=kE+1; jp=0
                NrconE(kE) = ndofEv*ndofHh
                do lpv=1,ndofEv
                  do lph=1,ndofHh
                    jp=jp+1
                    NacE(jp,kE)  = naE(loc) + ndofEh*ndofHv
     .                           + (lpv-1)*ndofHh+lph
                    ConstrE(jp,kE)=RRRH(1,lph,ish,lh)*RRRE(1,lpv,isv,lv)
                  enddo
                enddo
              enddo
            enddo
c
c           H(div)-loop through the constrained mid-face node dof's
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           horizontal first then vertical
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
c
c
c  -- HORIZONTAL MID-EDGE NODE CONSTRAINED BY AN H4-REFINED FACE --
c
          case(26,28)
c
c           determine horizontal son number
            select case(icase)
            case(28) ; is=1
            case(26) ; is=2
            endselect
c
c           parent mid-face node
            nodp = NFACEC(nc)
c
c           location of parent node on list of modified element's nodes
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
c
c           number of parent node horizontal and vertical dof's
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh=nordh-1 ; ndofEh=nordh
            ndofHv=nordv-1 ; ndofEv=nordv
c
c           parent mid-edge nodes (south,north)
            do ip=1,3,2
              nodp = abs(NFACE_CONS(ip,nc))
ccc           locp(ip) = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, locp(ip))
              ndofHp(ip) = NODES(nodp)%order-1
              ndofEp(ip) = NODES(nodp)%order
            enddo
c
c           H1::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
              kH=kH+1
              NrconH(kH) = ndofHh*ndofHv + ndofHp(1) + ndofHp(3)
c
c             loop through the parent mid-face node dof
              jp=0
              do lpv=1,ndofHv
                do lph=1,ndofHh
                  jp=jp+1
                  NacH(jp,kH)    = naH(loc)+(lpv-1)*ndofHh+lph
                  ConstrH(jp,kH) = RRRH(1,lph,is,l)*RRRH(1,lpv,3,1)
                enddo
              enddo
c
c             loop through the parent mid-edge nodes dofs
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
c
c           Hcurl::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
              kE=kE+1
              NrconE(kE) = ndofEh*ndofHv + ndofEp(1) + ndofEp(3)
c
c             loop through the parent mid-face node dof
              jp=0
              do lpv=1,ndofHv
                do lph=1,ndofEh
                  jp=jp+1
                  NacE(jp,kE)    = naE(loc)+(lpv-1)*ndofEh+lph
                  ConstrE(jp,kE) = RRRE(1,lph,is,l)*RRRH(1,lpv,3,1)
                enddo
              enddo
c
c             loop through the parent mid-edge nodes dof
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
c
c
c  -- VERTICAL MID-EDGE NODE CONSTRAINED BY AN H4-REFINED FACE --
c
          case(25,27)
c
c           determine vertical son number
            select case(icase)
            case(25); is=1
            case(27); is=2
            endselect
c
c           parent mid-face node
            nodp = NFACEC(nc)
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
c
c           parent mid-edge nodes (east,west)
            do ip=2,4,2
              nodp = abs(NFACE_CONS(ip,nc))
ccc           locp(ip) = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, locp(ip))
              ndofHp(ip) = NODES(nodp)%order-1
              ndofEp(ip) = NODES(nodp)%order
            enddo
c
c           H1::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
              kH=kH+1
              NrconH(kH) = ndofHh*ndofHv + ndofHp(2) + ndofHp(4)
c
c             loop through the parent mid-face node dofs
              jp=0
              do lpv=1,ndofHv
                do lph=1,ndofHh
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+(lpv-1)*ndofHh+lph
                  ConstrH(jp,kH) = RRRH(1,lph,3,1)*RRRH(1,lpv,is,l)
                enddo
              enddo
c
c             loop through the parent mid-edge nodes dof
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
c
c           Hcurl::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
              kE=kE+1
              NrconE(kE) = ndofHh*ndofEv + ndofEp(2) + ndofEp(4)
c
c             loop through the parent mid-face node dof
              jp=0
              do lpv=1,ndofEv
                do lph=1,ndofHh
                  jp=jp+1
                  NacE(jp,kE) = naE(loc) + ndofEh*ndofHv
     .                        + (lpv-1)*ndofHh+lph
                  ConstrE(jp,kE) = RRRH(1,lph,3,1)*RRRE(1,lpv,is,l)
                enddo
              enddo
c
c             loop through the parent mid-edge nodes dof
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
c
c
c  -- VERTEX NODE CONSTRAINED BY AN H4-REFINED FACE --
c
          case(29)
c
c           parent mid-face node
            nodp = NFACEC(nc)
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1
            ndofHv = nordv-1
c
c           parent mid-edge nodes (south,east,north,west)
            do ip=1,4
              nodp = abs(NFACE_CONS(ip,nc))
ccc           locp(ip) = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, locp(ip))
              ndofHp(ip) = NODES(nodp)%order-1
            enddo
c
c           only one dof
            kH=kH+1
            NrconH(kH) = ndofHh*ndofHv
     .                 + ndofHp(1)+ndofHp(2)+ndofHp(3)+ndofHp(4)+4
c
c           loop through the parent mid-face node dof
            jp=0
            do lpv=1,ndofHv
              do lph=1,ndofHh
                jp=jp+1
                NacH(jp,kH)    = naH(loc)+(lpv-1)*ndofHh+lph
                ConstrH(jp,kH) = RRRH(1,lph,3,1)*RRRH(1,lpv,3,1)
              enddo
            enddo
c
c           loop through the parent mid-edge node dof
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
c
c           loop through the parent vertex dof
            do ip=1,4
              nodp = NFACE_CONS(4+ip,nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
              jp=jp+1
              NacH(jp,kH)  = naH(loc)+1
              ConstrH(jp,kH) = .25d0
            enddo
c
c
c-----------------------------------------------------------------------
c  HORIZONTALLY H2-REFINED QUAD FACE
c-----------------------------------------------------------------------
c
c  -- MID-FACE NODE CONSTRAINED BY A HORIZONTALLY H2-REFINED FACE --
c
          case(31,32, 34,35, 61,62)
c
c           determine vertical son number
            select case(icase)
            case(31,34,61); is=1
            case(32,35,62); is=2
            end select
c
c           parent mid-face node
            nodp = NFACEC(nc)
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
c
c           H1::loop through the constrained mid-face node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do lv=1,ndofHv
              do lh=1,ndofHh
                kH=kH+1
                NrconH(kH) = ndofHv
c
c               loop through the parent mid-face node dof
                jp=0
                do lpv=1,ndofHv
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+(lpv-1)*ndofHh+lh
                  ConstrH(jp,kH) = RRRH(1,lpv,is,lv)
                enddo
              enddo
            enddo
c
c           Hcurl::loop through the constrained mid-face node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c           horizontal dof first
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
c
c           vertical dof next
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
c
c           Hdiv(L2)::loop through the constrained mid-face node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           horizontal first then vertical
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
c
c
c  -- HORIZONTAL MID-EDGE NODE CONSTRAINED BY A HORIZONTALLY H2-REFINED FACE --
c
          case(33,36,63)
c
c           parent mid-face node
            nodp = NFACEC(nc)
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
c
c           parent mid-edge nodes (south,north)
            do ip=1,3,2
              nodp = abs(NFACE_CONS(ip,nc))
              if (Is_active(nodp)) then
                ncflag(ip)=0
                call locate(nodp,Nodm,Nrnodm, locp(ip))
                ndofHp(ip) = 1
                ndofEp(ip) = 1
              else
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            ...added by S.P Mar 07 2018
                select case(NODES(nodp)%ref_kind)
                case default
                  ncflag(ip)=0
                  call locate(nodp,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = 1
                  ndofEp(ip) = 1
                case(0)
                  nfath = NODES(nodp)%father
c                  call locate(nodp,NODES(nfath)%sons(1:2),2, ncflag(ip))
                  ncflag(ip) = nodp - NODES(nfath)%first_son + 1
                  call locate(nfath,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = NODES(nfath)%order-1
                  ndofEp(ip) = NODES(nfath)%order
                end select
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif
            enddo
c
c           H1::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
              kH=kH+1
              NrconH(kH) = ndofHv+ndofHp(1)+ndofHp(3)
c
c             loop through the parent mid-face node dof
              jp=0
              do lpv=1,ndofHv
                jp=jp+1
                NacH(jp,kH)    = naH(loc)+(lpv-1)*ndofHh+l
                ConstrH(jp,kH) = RRRH(1,lpv,3,1)
              enddo
c
c             loop through the parent mid-edge nodes
              do ip=1,3,2
                select case(ncflag(ip))
c
c               unconstrained mid-edge node, only one parent dof
                case(0)
                  jp=jp+1
                  NacH(jp,kH) = naH(locp(ip))+l
                  if (NFACE_CONS(ip,nc).gt.0) then
                    ConstrH(jp,kH) = .5d0
                  else
                    ConstrH(jp,kH) = .5d0*(-1.d0)**(l-1)
                  endif
c
c               constrained mid-edge node
                case(1,2)
                  do lp=1,ndofHp(ip)
                    jp=jp+1
                    NacH(jp,kH) = naH(locp(ip))+lp
                    if (NFACE_CONS(ip,nc).gt.0) then
                      ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)
                    else
                      ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)
     .                               *(-1.d0)**(l-1)
                    endif
                  enddo
                end select
              enddo
            enddo
c
c           Hcurl::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
              kE=kE+1; jp=0
              NrconE(kE) = ndofHv+ndofEp(1)+ndofEp(3)
c
c             loop through the parent mid-face node dof
              do lpv=1,ndofHv
                jp=jp+1
                NacE(jp,kE)    = naE(loc)+(lpv-1)*ndofEh+l
                ConstrE(jp,kE) = RRRH(1,lpv,3,1)
              enddo
c
c             loop through the parent mid-edge nodes
              do ip=1,3,2
                select case(ncflag(ip))
c
c               unconstrained mid-edge node, only one parent dof
                case(0)
                  jp=jp+1
                  NacE(jp,kE) = naE(locp(ip))+l
                  if (NFACE_CONS(ip,nc).gt.0) then
                    ConstrE(jp,kE) = .5d0
                  else
                    ConstrE(jp,kE) = .5d0*(-1.d0)**(l)
                  endif
c
c               constrained mid-edge node
                case(1,2)
                  do lp=1,ndofEp(ip)
                    jp=jp+1
                    NacE(jp,kE) = naE(locp(ip))+lp
                    if (NFACE_CONS(ip,nc).gt.0) then
                      ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)
                    else
                      ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)
     .                  *(-1.d0)**(l)
                    endif
                  enddo
                end select
              enddo
            enddo
c
c
c-----------------------------------------------------------------------
c  VERTICALLY H2-REFINED QUAD FACE
c-----------------------------------------------------------------------
c
c  -- MID-FACE NODE CONSTRAINED BY A VERTICALLY H2-REFINED FACE --
c
          case(41,42, 44,45, 51,52)
c
c           determine vertical son number
            select case(icase)
            case(41,44,51); is=1
            case(42,45,52); is=2
            end select
c
c           parent mid-face node
            nodp = NFACEC(nc)
ccc         loc = NODES(nodp)%visit
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
c
c           H1::loop through the constrained mid-face node dof
            do lv=1,ndofHv
              do lh=1,ndofHh
                kH=kH+1
                NrconH(kH) = ndofHh
c
c           ....loop through the parent mid-face node dof
                jp=0
                do lph=1,ndofHh
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+(lv-1)*ndofHh+lph
                  ConstrH(jp,kH) = RRRH(1,lph,is,lh)
                enddo
              enddo
            enddo
c
c           Hcurl::loop through the constrained mid-face node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c           horizontal dof first
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
c
c           vertical dof next
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
c
c           Hdiv(L2)::loop through the constrained mid-face node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           horizontal first then vertical
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
c
c
c  -- VERTICAL MID-EDGE NODE CONSTRAINED BY A VERTICALLY H2-REFINED FACE --
c
          case(43,46,53)
c
c           parent mid-face node
            nodp = NFACEC(nc)
            call locate(nodp,Nodm,Nrnodm, loc)
            call decode(NODES(nodp)%order, nordh,nordv)
            ndofHh = nordh-1; ndofEh = nordh
            ndofHv = nordv-1; ndofEv = nordv
c
c           parent mid-edge nodes (east,west)
            do ip=2,4,2
              nodp = abs(NFACE_CONS(ip,nc))
              if (Is_active(nodp)) then
                ncflag(ip)=0
                call locate(nodp,Nodm,Nrnodm, locp(ip))
                ndofHp(ip) = 1
                ndofEp(ip) = 1
              else
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            ...added by S.P Mar 07 2018
                select case(NODES(nodp)%ref_kind)
                case default
                  ncflag(ip)=0
                  call locate(nodp,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = 1
                  ndofEp(ip) = 1
                case(0)
                  nfath = NODES(nodp)%father
c                  call locate(nodp,NODES(nfath)%sons(1:2),2, ncflag(ip))
                  ncflag(ip) = nodp - NODES(nfath)%first_son + 1
                  call locate(nfath,Nodm,Nrnodm, locp(ip))
                  ndofHp(ip) = NODES(nfath)%order-1
                  ndofEp(ip) = NODES(nfath)%order
                end select
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif
            enddo
c
c           H1::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
              kH=kH+1
              NrconH(kH) = ndofHh+ndofHp(2)+ndofHp(4)
c
c             loop through the parent mid-face node dof
              jp=0
              do lph=1,ndofHh
                jp=jp+1
                NacH(jp,kH)    = naH(loc)+(l-1)*ndofHh+lph
                ConstrH(jp,kH) = RRRH(1,lph,3,1)
              enddo
c
c             loop through the parent mid-edge nodes
              do ip=2,4,2
                select case(ncflag(ip))
c
c               unconstrained mid-edge node, only one parent dof
                case(0)
                  jp=jp+1
                  NacH(jp,kH)    = naH(locp(ip))+l
                  if (NFACE_CONS(ip,nc).gt.0) then
                    ConstrH(jp,kH) = .5d0
                  else
                    ConstrH(jp,kH) = .5d0*(-1.d0)**(l-1)
                  endif
c
c               constrained mid-edge node
                case(1,2)
                  do lp=1,ndofHp(ip)
                    jp=jp+1
                    NacH(jp,kH) = naH(locp(ip))+lp
                    if (NFACE_CONS(ip,nc).gt.0) then
                      ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)
                    else
                      ConstrH(jp,kH) = .5d0*RRRH(1,lp,ncflag(ip),l)
     .                               *(-1.d0)**(l-1)
                    endif
                  enddo
                end select
              enddo
            enddo
c
c           Hcurl::loop through the constrained mid-edge node dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
              kE=kE+1; jp=0
              NrconE(kE) = ndofHh+ndofEp(2)+ndofEp(4)
c
c             loop through the parent mid-face node dof
              do lph=1,ndofHh
                jp=jp+1
                NacE(jp,kE)    = naE(loc)+ndofEh*ndofHv
     .                         + (l-1)*ndofHh+lph
                ConstrE(jp,kE) = RRRH(1,lph,3,1)
              enddo
c
c             loop through the parent mid-edge nodes
              do ip=2,4,2
                select case(ncflag(ip))
c
c               unconstrained mid-edge node, only one parent dof
                case(0)
                  jp=jp+1
                  NacE(jp,kE) = naE(locp(ip))+l
                  if (NFACE_CONS(ip,nc).gt.0) then
                    ConstrE(jp,kE) = .5d0
                  else
                    ConstrE(jp,kE) = .5d0*(-1.d0)**(l)
                  endif
c
c               constrained mid-edge node
                case(1,2)
                  do lp=1,ndofEp(ip)
                    jp=jp+1
                    NacE(jp,kE) = naE(locp(ip))+lp
                    if (NFACE_CONS(ip,nc).gt.0) then
                      ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)
                    else
                      ConstrE(jp,kE) = .5d0*RRRE(1,lp,ncflag(ip),l)
     .                               *(-1.d0)**(l)
                    endif
                  enddo
                end select
              enddo
            enddo
c
c
c-----------------------------------------------------------------------
c  H4-REFINED TRIANGULAR FACE
c-----------------------------------------------------------------------
c
c  -- MID-FACE NODE CONSTRAINED BY AN H4-REFINED FACE --
c
          case(71,72,73,74)
            call decode(icase, nvoid,is)
c
c           H1::loop through the nodal dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofH
              kH=kH+1
c
c             initiate parent dof counter
              jp=0
c
c             parent mid-face node
              nodp = NFACEC(nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
c
c             loop through parent node dof
              do lp=1,NdofmH(loc)
                jp=jp+1
                NacH(jp,kH) = naH(loc)+lp
                ConstrH(jp,kH) = RRTH(1,lp,is,l)
              enddo
c
c             loop through parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
c
c               loop through parent node dof
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
c
c             save the number of constraining dof
              NrconH(kH) = jp
            enddo
c
c           Hcurl::loop through the nodal dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofE
              kE=kE+1
c
c           ..initiate parent dof counter
              jp=0
c
c           ..parent mid-face node
              nodp = NFACEC(nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
c
c           ..loop through parent node dof
              do lp=1,NdofmE(loc)
                jp=jp+1
                NacE(jp,kE) = naE(loc)+lp
                ConstrE(jp,kE) = RRTE(1,lp,is,l)
              enddo
c
c           ..loop through parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
c
c           ....loop through parent node dof
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
c
c           ..save the number of constraining dof
              NrconE(kE) = jp
            enddo
c
c           H(div)::loop through the nodal dof
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do l=1,ndofV
              kV=kV+1
c
c           ..initiate parent dof counter
              jp=0
c
c           ..parent mid-face node
              nodp = NFACEC(nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
c
c           ..loop through parent node dof
              do lp=1,NdofmV(loc)
                jp=jp+1
                NacV(jp,kV) = naV(loc)+lp
                ConstrV(jp,kV) = RRTQ(lp,is,l)
              enddo
c
c           ..save the number of constraining dof
              NrconV(kV) = jp
            enddo
c
c
c  -- MID-EDGE NODE CONSTRAINED BY AN H4-REFINED FACE --
c
          case(75,76,77)
            call decode(icase, nvoid,is)
c
c           H1::loop through the nodal dof
            do l=1,ndofH
              kH=kH+1
c
c             initiate parent dof counter
              jp=0
c
c             parent mid-face node
              nodp = NFACEC(nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
c
c             loop through parent mdlt node dof
              do lp=1,NdofmH(loc)
                jp=jp+1
                NacH(jp,kH) = naH(loc)+lp
                ConstrH(jp,kH) = RRTH(1,lp,is,l)
              enddo
c
c             loop through parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
c
#if HP3D_DEBUG
                if (iprint.eq.1) then
                  write(*,*) 'logic: ie,nodp,loc = ',ie,nodp,loc
                endif
#endif
c
c               loop through parent medg node dof
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
#if HP3D_DEBUG
              if (iprint.eq.1) then
                write(*,*) 'l = ',l
                write(*,7013) NacH(1:jp,kH)
                write(*,7014) ConstrH(1:jp,kH)
              endif
#endif
c
c             save the number of constraining dof
              NrconH(kH) = jp
            enddo
c
c           Hcurl::loop through the nodal dof
            do l=1,ndofE
              kE=kE+1
c
c             initiate parent dof counter
              jp=0
c
c             parent mid-face node
              nodp = NFACEC(nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
c
c             loop through parent mdlt node dof
              do lp=1,NdofmE(loc)
                jp=jp+1
                NacE(jp,kE) = naE(loc)+lp
                ConstrE(jp,kE) = RRTE(1,lp,is,l)
              enddo
c
c             loop through parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
#if HP3D_DEBUG
                if (iprint.eq.1) then
                  write(*,*) 'logic: ie,nodp,loc = ',ie,nodp,loc
                endif
#endif
c
c               loop through parent medg node dof
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
#if HP3D_DEBUG
              if ((iprint.eq.1).and.(Idec.gt.1)) then
                write(*,*) '75,76,77::Mdle,l,jp = ',Mdle,l,jp
                write(*,7013) NacE(1:jp,kE)
                write(*,7014) ConstrE(1:jp,kE)
              endif
#endif
c
c             save the number of constraining dof
              NrconE(kE) = jp
            enddo
c
c
c-----------------------------------------------------------------------
c  H2-REFINED TRIANGULAR FACE (QUAD + TRIANGLE)
c-----------------------------------------------------------------------
c
c  -- MDLQ NODE CONSTRAINED BY AN H2-REFINED TRIANGULAR FACE --
c
          case(82,83,84)
            call decode(icase, nvoid,iref)
c
c           loop through the nodal dof
            do l=1,ndofH
              kH=kH+1
c
c             initiate parent dof counter
              jp=0
c
c             parent mdlt node
              nodp = NFACEC(nc)
ccc           loc = NODES(nodp)%visit
              call locate(nodp,Nodm,Nrnodm, loc)
c
c             loop through parent mdlt node dof
              do lp=1,NdofmH(loc)
                jp=jp+1
                NacH(jp,kH) = naH(loc)+lp
                np = NODES(nodp)%order
                ConstrH(jp,kH) = get_rrqh(iref,np,1,lp,l)
              enddo
c
c             loop through parent medg nodes
              do ie=1,3
                nodp = abs(NFACE_CONS(ie,nc))
ccc             loc = NODES(nodp)%visit
                call locate(nodp,Nodm,Nrnodm, loc)
#if HP3D_DEBUG
                if (iprint.eq.1) then
                  write(*,*) 'logic: ie,nodp,loc = ',ie,nodp,loc
                endif
#endif
c
c               loop through parent medg node dof
                do lp=1,NdofmH(loc)
                  jp=jp+1
                  NacH(jp,kH) = naH(loc)+lp
                  if (NFACE_CONS(ie,nc).gt.0) then
                    ConstrH(jp,kH) = get_rrqh(iref,np,1+ie,lp,l)
                  else
                    ConstrH(jp,kH) = get_rrqh(iref,np,1+ie,lp,l)
     .                               *(-1.d0)**(lp+1)
                  endif
                enddo
              enddo
#if HP3D_DEBUG
              if (iprint.eq.1) then
                write(*,*) 'l = ',l
                write(*,7013) NacH(1:jp,kH)
                write(*,7014) ConstrH(1:jp,kH)
              endif
#endif
c
c             save the number of constraining dof
              NrconH(kH) = jp
c
            enddo
c-----------------------------------------------------------------------
c
          end select
c
c  .....if a constrained node
        endif
c
#if HP3D_DEBUG
c  .....printing
        if (iprint.eq.1) then
          write(*,9000) i,nod,kH,kE,kV
 9000     format(' i = ',i2,' ; nod = ',i5,' ; kH = ',i3,
     .                                     ' ; kE = ',i3,
     .                                     ' ; kV = ',i3 )
        endif
#endif
c
c  ...end of loop through local nodes
      enddo
c
c  ...sanity check
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
c
c  ...clean visitation flag
 100  continue
ccc    do im=1,Nrnodm
ccc      NODES(Nodm(im))%visit=0
ccc    enddo
c
c-----------------------------------------------------------------------
c
c  UNCONSTRAINED VERSION FOR H(div)
c
      if (.not.CONSTRAINED_HDIV) then
        call find_order(Mdle, norder)
        call celndof(NODES(Mdle)%ntype,norder,
     .               nrdofH,nrdofE,nrdofV,nrdofQ)
c
        NrconV=0 ; NacV=0 ; ConstrV=0.d0
        do k=1,nrdofV
          NrconV(k)=1 ; NacV(1,k)=k ; ConstrV(1,k)=1.d0
        enddo
      endif
c
c-----------------------------------------------------------------------
c
#if HP3D_DEBUG
      if ((iprint.ge.1).and.(Idec.ne.1)) then
        write(*,7011) Mdle
c
        write(*,*) ' -- H1 --'
        nrdofH_loc=kH
        do kH=1,nrdofH_loc
          write(*,7012) kH,NrconH(kH)
          write(*,7013) NacH(1:NrconH(kH),kH)
          write(*,7014) ConstrH(1:NrconH(kH),kH)
        enddo
c
        write(*,*) ' '
        write(*,*) ' -- H(curl) --'
        nrdofE_loc=kE
        do kE=1,nrdofE_loc
          write(*,7012) kE,NrconE(kE)
          write(*,7013) NacE(1:NrconE(kE),kE)
          write(*,7014) ConstrE(1:NrconE(kE),kE)
        enddo
c
        write(*,*) ' -- H(div) --'
        nrdofV_loc=kV
        do kV=1,nrdofV_loc
          write(*,7012) kV,NrconV(kV)
          write(*,7013) NacV(1:NrconV(kV),kV)
          write(*,7014) ConstrV(1:NrconV(kV),kV)
        enddo
        call pause
c
 7011   format(' logic: Mdle = ',i6)
 7012   format(' logic: k = ',i4,' Nrcon = ',i3)
 7013   format(20i6)
 7014   format(20f6.3)
c
      endif
#endif
c
c
      end subroutine logic
