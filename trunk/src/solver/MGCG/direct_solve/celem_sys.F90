!-----------------------------------------------------------------------
!
!   routine name       - celem_system
!
!-----------------------------------------------------------------------
!
!   latest revision    - May 10
!
!   purpose            - routine calculates load vector and
!                        stiffness matrix for a modified element
!                        for a coupled problem
!
!            ATTENTION: the modified element arrays correspond
!            to the order of dofs induced by different type
!            of discretization: H1, H(curl), H(div) and L2
!            and by the nodes:
!              middle node first, then
!              mid-face nodes,
!              mid-edge nodes,
!              vertex nodes;
!            within each node, the corresponding dof corresponding
!            to the same kind of discretization
!            are enumerated in the order dictated by index
!
!   arguments :
!     in:
!             Mdle     - middle node of an element
!             Idec     = 1 return only the info on nodes
!                      = 2 compute also the element matrices
!     out:
!             Nrdofs   - number of element local dof for each physical
!                        attribute (needed to dimension work arrays
!                        in calling routines)
!             Nrdofm   - number of modified element dof in the expanded
!                        mode (needed to dimension work arrays in calling
!                        routines)
!             Nrdofc   - number of modified element dof for the coupled
!                        problem after compression
!
!             Nodm     - actual (unconstrained) nodes in the order:
!                        middle,mid-face,mid-edge, vertex nodes
!             NdofmH   - the corresponding number of H1 dof
!             NdofmE   - the corresponding number of H(curl) dof
!             NdofmV   - the corresponding number of H(div) dof
!             NdofmQ   - the corresponding number of L2 dof
!             Nrnodm   - number of the modified element nodes
!
!             Zbload   - 1D array containing the modified load vector
!             Zastif   - 1D array containing the modified stiffness
!                        matrix
!
!----------------------------------------------------------------------
!
subroutine celem_sys(Mdle,Idec,                                &
                        Nrdofs,Nrdofm,Nrdofc,                     &
                        Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,  &
                        Zbload,Zastif)
!
      use data_structure3D
      use assembly
      use control
#include "syscom.blk"
!
      common /ccelem/ iprint
!
      dimension Nrdofs(NR_PHYSA)
      dimension Nodm(MAXNODM),NdofmH(MAXNODM),NdofmE(MAXNODM), &
                              NdofmV(MAXNODM),NdofmQ(MAXNODM)
      dimension Zbload(*), Zastif(*)
!
!  ...element order
      dimension norder(19)
!
!  ...index for a node
      dimension index(NRINDEX)
!
!  ...number of single H1 and H(curl) or H(div) dof for a higher order
!     node
      dimension ndofmHl(MAXNODM),ndofmEl(MAXNODM),ndofmVl(MAXNODM)
!
!  ...location of first dof for a modified element node and
!     discretization type
      dimension namH(MAXNODM),namE(MAXNODM),namV(MAXNODM)
!
!  ...constrained approximation coefficients for a scalar-valued
!     element
      dimension nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV), &
                nacH(NACDIM,MAXbrickH),constrH(NACDIM,MAXbrickH),      &
                nacE(NACDIM,MAXbrickE),constrE(NACDIM,MAXbrickE),      &
                nacV(NACDIM,MAXbrickV),constrV(NACDIM,MAXbrickV)
!
!  ...test and trial flags
      dimension itest(NR_PHYSA), jtrial(NR_PHYSA)
!
!  ...for supporting non-conforming H(div)
!!      logical :: CONSTRAINED_HDIV = .FALSE.
      logical :: CONSTRAINED_HDIV = .TRUE.
!
!----------------------------------------------------------------------
!
      iprint=0
!
!  ...initialize
      Nrdofs(1:NR_PHYSA)=0
      Nodm(  1:MAXNODM )=0
      NdofmH(1:MAXNODM )=0
      NdofmE(1:MAXNODM )=0
      NdofmV(1:MAXNODM )=0
      NdofmQ(1:MAXNODM )=0
!
!  ...number of dofs for a SINGLE H1,H(curl),H(div),L2 component
      call find_order(Mdle, norder)
      call celndof(NODES(Mdle)%type,norder, nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
      if (iprint.eq.1) then
        write(*,5001) Mdle
 5001   format(' celem_system: DEBUGGING FOR Mdle = ',i8)
        call print_order(NODES(Mdle)%type,norder)
      endif
!
!  ...number of dofs for ALL H1,H(curl),H(div),L2 components (EXPANDED mode)
      do i=1,NR_PHYSA
        select case(DTYPE(i))
        case('contin') ; Nrdofs(i)=nrdoflH*NR_COMP(i)
        case('tangen') ; Nrdofs(i)=nrdoflE*NR_COMP(i)
        case('normal') ; Nrdofs(i)=nrdoflV*NR_COMP(i)
        case('discon') ; Nrdofs(i)=nrdoflQ*NR_COMP(i)
        endselect
      enddo
!
!  ...determine the constraints coefficients
      call logic(Mdle,Idec,                             &
                 Nodm,ndofmHl,ndofmEl,ndofmVl,Nrnodm,   &
                 nrconH,nacH,constrH,                   &
                 nrconE,nacE,constrE,                   &
                 nrconV,nacV,constrV)
!
!  ...printing
      if ((iprint.ge.1).and.(Idec.ne.1)) then
        write(*,*) ' - H1 - '
        do k=1,nrdoflH
          write(*,2000) nrconH(k),k
          write(*,2100) nacH(1:nrconH(k),k)
          write(*,2200) constrH(1:nrconH(k),k)
 2000     format(' nrconH(k),k  = ',i2,3x,i10)
 2100     format('    nacH(:,k) = ',20(i10,2x))
 2200     format(' constrH(:,k) = ',20(e12.5,2x))
        enddo
!
        write(*,*) ' - Hcurl - '
        do k=1,nrdoflE
          write(*,1000) nrconE(k),k
          write(*,1100) nacE(1:nrconE(k),k)
          write(*,1200) constrE(1:nrconE(k),k)
 1000     format(' nrconE(k),k  = ',i2,3x,i10)
 1100     format('    nacE(:,k) = ',20(i10,2x))
 1200     format(' constrE(:,k) = ',20(e12.5,2x))
        enddo
        call pause
!
        write(*,*) ' - Hdiv - '
        do k=1,nrdoflV
          write(*,3000) nrconV(k),k
          write(*,3100) nacV(1:nrconV(k),k)
          write(*,3200) constrV(1:nrconV(k),k)
 3000     format(' nrconV(k),k  = ',i2,3x,i10)
 3100     format('    nacV(:,k) = ',20(i10,2x))
 3200     format(' constrV(:,k) = ',20(e12.5,2x))
        enddo
        call pause
      endif
!
!  ...number of H1 dofs for the modified element (EXPANDED mode)
      k=0
      do i=1,Nrnodm
        namH(i) = k
        k = k + ndofmHl(i)
      enddo
      nrdofmH = k*NRHVAR
!
!  ...number of H(curl) dofs for the modified element (EXPANDED mode)
      k=0
      do i=1,Nrnodm
        namE(i) = k
        k = k + ndofmEl(i)
      enddo
      nrdofmE = k*NREVAR
!
!  ...number of H(div) dofs for the modified element (EXPANDED mode)
      k=0
      do i=1,Nrnodm
        namV(i) = k
        k = k + ndofmVl(i)
      enddo
      nrdofmV = k*NRVVAR
!
IF (CONSTRAINED_HDIV) THEN
!
      nrdofmHE  = nrdofmH   + nrdofmE
      nrdofmHEV = nrdofmHE  + nrdofmV
      Nrdofm    = nrdofmHEV + nrdoflQ*NRQVAR
!
ELSE
!
      nrdofmHE  = nrdofmH   + nrdofmE
      nrdofmHEV = nrdofmHE  + nrdoflV*NRVVAR
      Nrdofm    = nrdofmHEV + nrdoflQ*NRQVAR
!
ENDIF
!
      if (iprint .eq. 1) then
        write(*,6001) nrdoflH,nrdoflE,nrdoflV,nrdoflQ
 6001   format('celem: nrdoflH,nrdoflE,nrdoflV,nrdoflQ = ',4i4)
        write(*,6002) nrdofmH,nrdofmHE,nrdofmHEV,Nrdofm
 6002   format('celem: nrdofmH,nrdofmHE,nrdofmHEV,Nrdofm = ',4i4)
      endif
!
!***********************************************************************
!
!  ...establish extraction vector and Dirichlet data
      ZDOFD = ZERO
      IDBC  = 0
!
!  ...initiate global dof counter for the compressed mode
      k=0
!
!---------------------------- H1 dof -----------------------------------
!
!  ...loop through nodes in the reversed order
      do i=1,Nrnodm
        nod = Nodm(i)
        nvarHt = NREQNH(NODES(nod)%case) ! total number of H1 components for a single load
        l = namH(i)*NRHVAR
        call get_index(nod, index)
        kold = k
!
!  .....loop through H1 dof for the node
        do j=1,ndofmHl(i)
          ivar=0
          do ii=1,NRHVAR
            select case(index(ii))
!
!  .........dof present but known from Dirichlet BC, save the BC data
            case(1)
              l=l+1; ivar=ivar+1
              IDBC(l)=1
              do iload=1,NR_RHS
                ZDOFD(l,iload) = NODES(nod)%zdofH((iload-1)*nvarHt+ivar,j)
              enddo
!
!  .........dof present and active
            case(2)
              k=k+1; l=l+1; ivar=ivar+1
              NEXTRACT(k) = l
!
!  .........dof absent, skip it but update the expanded dof counter
            case(0)
              l=l+1
            case default
              write(*,*) 'celem: nod,ii,index(ii) = ',nod,ii,index(ii)
              stop
            end select
          enddo
        enddo
!
!  .....store the actual number of dof
        NdofmH(i) = k-kold
!
!  ...end of loop through nodes
      enddo
!
!---------------------------- H(curl) dof ------------------------------
!
!  ...loop through nodes in reversed order
      do i=1,Nrnodm
        nod = Nodm(i)
        nvarEt = NREQNE(NODES(nod)%case) ! total number of H(curl) components for a single load
        l = nrdofmH + namE(i)*NREVAR
        call get_index(nod, index)
        kold = k
!
!  .....loop through H(curl) dof for the node
        do j=1,ndofmEl(i)
          ivar=0
          do ii=NRHVAR+1,NRHVAR+NREVAR
            select case(index(ii))
!
!  .........dof present but known from Dirichlet BC, save the BC data
            case(3)
              l=l+1; ivar=ivar+1
              IDBC(l)=1
              do iload=1,NR_RHS
                ZDOFD(l,iload) = NODES(nod)%zdofE((iload-1)*nvarEt+ivar,j)
              enddo
!
!  .........dof present and active
            case(4)
              k=k+1; l=l+1; ivar=ivar+1
              NEXTRACT(k) = l
!
!  .........dof absent, skip it but update the local dof counter
            case(0)
              l=l+1
            case default
              write(*,*) 'celem: nod,ii,index(ii) = ',nod,ii,index(ii)
              stop
            end select
          enddo
        enddo
!
!  .....compute the actual number of dof
        NdofmE(i) = k-kold
!
!  ...end of loop through nodes
      enddo
!
!---------------------------- H(div) dof -------------------------------
!
IF (CONSTRAINED_HDIV) THEN
!
!  ...loop through nodes in the reversed order
      do i=1,Nrnodm
        nod = Nodm(i)
        nvarVt = NREQNV(NODES(nod)%case) ! total number of H(div) components for a single load
        l = nrdofmHE + namV(i)*NRVVAR
        call get_index(nod, index)
        kold = k
!
!  .....loop through H(div) dof for the node
        do j=1,ndofmVl(i)
          ivar=0
          do ii=NRHVAR+NREVAR+1,NRHVAR+NREVAR+NRVVAR
            select case(index(ii))
!
!  .........dof present but known from Dirichlet BC, save the BC data
            case(5)
              l=l+1; ivar=ivar+1
              IDBC(l)=1
              do iload=1,NR_RHS
                ZDOFD(l,iload) = NODES(nod)%zdofV((iload-1)*nvarVt+ivar,j)
              enddo
!
!  .........dof present and active
            case(6)
              k=k+1; l=l+1; ivar=ivar+1
              NEXTRACT(k) = l
!
!  .........dof absent, skip it but update the local dof counter
            case(0)
              l=l+1
            case default
              write(*,*) 'celem: nod,ii,index(ii) = ',nod,ii,index(ii)
              stop
            end select
          enddo
        enddo
!
!  .....compute the actual number of dof
        NdofmV(i) = k-kold
!
!  ...end of loop through nodes
      enddo
!
ELSE
!
      nod = Mdle
      nvarVt = NREQNV(NODES(nod)%case) ! total number of H(div) components for a single load
      call get_index(nod, index)
      kold = k
      call find_ndof(nod, nvoid,nvoid,ndofV,nvoid)
      l = nrdofmHE
!
!  ...loop through H(div) dof for the node
      do j=1,ndofV
        ivar=0
        do ii=NRHVAR+NREVAR+1,NRHVAR+NREVAR+NRVVAR
          select case(index(ii))
!
!  .......dof present but known from Dirichlet BC, save the BC data
          case(5)
            l=l+1; ivar=ivar+1
            IDBC(l)=1
            do iload=1,NR_RHS
              ZDOFD(l,iload) = NODES(nod)%zdofV((iload-1)*nvarVt+ivar,j)
            enddo
!
!  .......dof present and active
          case(6)
            k=k+1; l=l+1; ivar=ivar+1
            NEXTRACT(k) = l
!
!  .......dof absent, skip it but update the local dof counter
          case(0)
            l=l+1
          case default
            write(*,*) 'celem: nod,ii,index(ii) = ',nod,ii,index(ii)
            stop
          end select
        enddo
      enddo
!
!  ...compute the actual number of dof
      NdofmV(1:Nrnodm-1) = 0
      NdofmV(Nrnodm) = k-kold
!
ENDIF
!
!---------------------------- L2 dof -----------------------------------
!
!  ...loop through higher order nodes, middle node first
      nod = Mdle
      nvarQt = NREQNQ(NODES(nod)%case) ! total number of L2 components for a single load
      call get_index(nod, index)
      kold = k
      call find_ndof(nod, nvoid,nvoid,nvoid,ndofQ)
      l = nrdofmHEV
!
!  ...loop through L2 dof for the node
      do j=1,ndofQ
        ivar=0
        do ii=NRHVAR+NREVAR+NRVVAR+1,NRINDEX
          select case(index(ii))
!
!  .......dof present but known from Dirichlet BC, save the BC data
          case(7)
            l=l+1; ivar=ivar+1
            IDBC(l)=1
            do iload=1,NR_RHS
              ZDOFD(l,iload) = NODES(nod)%zdofQ((iload-1)*nvarQt+ivar,j)
            enddo
!
!  .......dof present and active
          case(8)
            k=k+1; l=l+1; ivar=ivar+1
            NEXTRACT(k) = l
!
!  .......dof absent, skip it but update the local dof counter
          case(0)
            l=l+1
          case default
            write(*,*) 'celem: nod,ii,index(ii) = ',nod,ii,index(ii)
            stop
          end select
        enddo
      enddo
!
!  ...compute the actual number of dof
      NdofmQ(1:Nrnodm-1) = 0
      NdofmQ(Nrnodm) = k-kold
!
!-----------------------------------------------------------------------
!
!  ...save number of dof for the modified element in the compressed
!     mode
      Nrdofc = k
!
!  ...printing
      if (iprint.eq.1) then
        write(*,7002) (Nodm(i),NdofmH(i),NdofmE(i),NdofmV(i), &
                        NdofmQ(i),i=1,Nrnodm)
 7002   format('       Nodm,NdofmH,NdofmE,NdofmV,NdofmQ = ',  &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x),           &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x),           &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x),           &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x))
        write(*,7003) Nrdofc
 7003   format('celem: Nrdofc = ',i4)
        write(*,7004) NEXTRACT(1:Nrdofc)
 7004   format('celem: NEXTRACT = ',10(/,20i5))
        write(*,7005) IDBC(1:Nrdofm)
 7005   format('celem: IDBC = ',10(4i2,1x),                   &
             10(/,'              ',10(4i2,1x)))
        call pause
      endif
!
!  ...return if only INFO is needed
      if (Idec.eq.1) return
!
!***********************************************************************
!
!  ...determine the element local load vector and stiffness matrix
      call elem(Mdle, itest,jtrial)
!
!  ...initiate modified load vector and half-way modified
!     auxiliary stiffness matrices
      ZBMOD(1:Nrdofm,1:NR_RHS) = ZERO
      do i=1,NR_PHYSA
        AAUX(i)%array(1:Nrdofm,1:Nrdofs(i)) = ZERO
      enddo
!
!  ...compute the modified element matrices in the expanded form
!
!  ...first loop through physical attributes
      do iphys1=1,NR_PHYSA
!
!  .....skip if the attribute is absent
        if (itest(iphys1).eq.0) cycle
!
        select case(DTYPE(iphys1))
        case('contin')
!
!  .......loop through element test functions
          do k=1,nrdoflH
!
!  .........loop through the connected real dofs
            do kp=1,nrconH(k)
!
!  ...........modified matrix location
              l = nacH(kp,k)
!
!  ...........loop through components
              do ivar=1,NR_COMP(iphys1)
                ll = (l-1)*NRHVAR+ADRES(iphys1)+ivar
                kk = (k-1)*NR_COMP(iphys1)+ivar
                ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)   &
                        + BLOC(iphys1)%array(kk,1:NR_RHS)*constrH(kp,k)
!
!  .............second loop through physical attributes
                do iphys2=1,NR_PHYSA
!
!  ...............skip if the attribute is absent
                  if (jtrial(iphys2).eq.0) cycle
!
!  ...............transform the stiffness matrix
                  AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                  = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                  + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))  &
                    *constrH(kp,k)
!
!  .............end of second loop through physical attributes
                enddo
!
!  ...........end of loop through components
              enddo
!
!  .........end of loop through connected parent dof
            enddo
!
!  .......end of loop through element test functions
          enddo
!
       case('tangen')
!
!  .......loop through element test functions
          do k=1,nrdoflE
!
!  .........loop through the connected real dofs
            do kp=1,nrconE(k)
!
!  ...........modified matrix location
              l = nacE(kp,k)
!
!  ...........loop through components
              do ivar = 1, NR_COMP(iphys1)
                ll = nrdofmH + (l - 1)*NREVAR + ADRES(iphys1) + ivar
                kk = (k - 1)*NR_COMP(iphys1) + ivar
                ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)  &
                  + BLOC(iphys1)%array(kk,1:NR_RHS)*constrE(kp,k)
!
!  .............second loop through physical attributes
                do iphys2=1,NR_PHYSA
!
!  ...............skip if the attribute is absent
                  if (jtrial(iphys2).eq.0) cycle
!
!  ...............transform the stiffness matrix
                  AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                  = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                  + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))  &
                    *constrE(kp,k)
!
!  .............end of second loop through physical attributes
                enddo
!
!  ...........end of loop through components
              enddo
!
!  .........end of loop through connected parent dof
            enddo
!
!  .......end of loop through element test functions
          enddo
!
       case('normal')
!
IF (CONSTRAINED_HDIV) THEN
!
!  .......loop through element test functions
          do k=1,nrdoflV
!
!  .........loop through the connected real dofs
            do kp=1,nrconV(k)
!
!  ...........modified matrix location
              l = nacV(kp,k)
!
!  ...........loop through components
              do ivar=1,NR_COMP(iphys1)
                ll = nrdofmHE + (l-1)*NRVVAR+ADRES(iphys1)+ivar
                kk = (k-1)*NR_COMP(iphys1)+ivar
                ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)              &
                      + BLOC(iphys1)%array(kk,1:NR_RHS)*constrV(kp,k)
!
!  .............second loop through physical attributes
                do iphys2=1,NR_PHYSA
!
!  ...............skip if the attribute is absent
                  if (jtrial(iphys2).eq.0) cycle
!
!  ...............transform the stiffness matrix
                  AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                  = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                  + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))  &
                    *constrV(kp,k)
!
!  .............end of second loop through physical attributes
                enddo
!
!  ...........end of loop through components
              enddo
!
!  .........end of loop through connected parent dof
            enddo
!
!  .......end of loop through element test functions
          enddo
!
ELSE
!
!  .......loop through element test functions
          do k=1,nrdoflV
!
!  .........loop through components
            do ivar=1,NR_COMP(iphys1)
              ll = nrdofmHE + (k-1)*NRVVAR+ADRES(iphys1)+ivar
              kk = (k-1)*NR_COMP(iphys1)+ivar
              ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)              &
                    + BLOC(iphys1)%array(kk,1:NR_RHS)
!
!  ...........second loop through physical attributes
              do iphys2=1,NR_PHYSA
!
!  .............skip if the attribute is absent
                if (jtrial(iphys2).eq.0) cycle
!
!  .............transform the stiffness matrix
                AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))
!
!  ...........end of second loop through physical attributes
              enddo
!
!  .........end of loop through components
            enddo
!
!  .......end of loop through element test functions
          enddo
!
ENDIF
!
       case('discon')
!
!  .......loop through element test functions
          do k=1,nrdoflQ
!
!  .........loop through components
            do ivar=1,NR_COMP(iphys1)
              ll = nrdofmHEV + (k-1)*NRQVAR+ADRES(iphys1)+ivar
              kk = (k-1)*NR_COMP(iphys1)+ivar
              ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)      &
                        + BLOC(iphys1)%array(kk,1:NR_RHS)
!
!  ...........second loop through physical attributes
              do iphys2=1,NR_PHYSA
!
!  .............skip if the attribute is absent
                if (jtrial(iphys2).eq.0) cycle
!
!  .............transform the stiffness matrix
                AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))     &
                = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))   &
                + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))
!
!  ...........end of second loop through physical attributes
              enddo
!
!  .........end of loop through components
            enddo
!
!  .......end of loop through element test functions
          enddo
!
        end select
!
!  ...end of first loop through physical attributes
      enddo
!
!
!***********************************************************************
!
!  ...initiate modified stiffness matrix
      ZAMOD(1:Nrdofm,1:Nrdofm) = ZERO
!
!  ...loop through physical attributes
      do iphys=1,NR_PHYSA
!
!  .....skip if the attribute is absent
        if (jtrial(iphys).eq.0) cycle
!
        select case(DTYPE(iphys))
        case('contin')
!
!  .......loop through element trial functions
          do k=1,nrdoflH
!
!  .........loop through the connected real dofs
            do kp=1,nrconH(k)
!
!  ...........modified matrix location
              l = nacH(kp,k)
!
!  ...........loop through components
              do ivar=1,NR_COMP(iphys)
                ll = (l-1)*NRHVAR+ADRES(iphys)+ivar
                kk = (k-1)*NR_COMP(iphys)+ivar
!
!  .............accumulate
                ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
                + AAUX(iphys)%array(1:Nrdofm,kk)*constrH(kp,k)
!
!  ...........end of loop through components
              enddo
!
!  .........end of loop through connected parent dof
            enddo
!
!  .......end of loop through element trial functions
          enddo
!
       case('tangen')
!
!  .......loop through element trial functions
          do k=1,nrdoflE
!
!  .........loop through the connected real dofs
            do kp=1,nrconE(k)
!
!  ...........modified matrix location
              l = nacE(kp,k)
!
!  ...........loop through components
              do ivar=1,NR_COMP(iphys)
                ll = nrdofmH + (l-1)*NREVAR+ADRES(iphys)+ivar
                kk = (k-1)*NR_COMP(iphys)+ivar
!
!  .............accumulate
                ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
                + AAUX(iphys)%array(1:Nrdofm,kk)*constrE(kp,k)
!
!  ...........end of loop through components
              enddo
!
!  .........end of loop through connected parent dof
            enddo
!
!  .......end of loop through element trial functions
          enddo
!
       case('normal')
!
IF (CONSTRAINED_HDIV) THEN
!
!  .......loop through element trial functions
          do k=1,nrdoflV
!
!  .........loop through the connected real dofs
            do kp=1,nrconV(k)
!
!  ...........modified matrix location
              l = nacV(kp,k)
!
!  ...........loop through components
              do ivar=1,NR_COMP(iphys)
                ll = nrdofmHE + (l-1)*NRVVAR+ADRES(iphys)+ivar
                kk = (k-1)*NR_COMP(iphys)+ivar
!
!  .............accumulate
                ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
                + AAUX(iphys)%array(1:Nrdofm,kk)*constrV(kp,k)
!
!  ...........end of loop through components
              enddo
!
!  .........end of loop through connected parent dof
            enddo
!
!  .......end of loop through element trial functions
          enddo
!
ELSE
!
!  .......loop through element trial functions
          do k=1,nrdoflV
!
!  .........loop through components
            do ivar=1,NR_COMP(iphys)
              ll = nrdofmHE + (k-1)*NRVVAR+ADRES(iphys)+ivar
              kk = (k-1)*NR_COMP(iphys)+ivar
!
!  ...........accumulate
              ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
              + AAUX(iphys)%array(1:Nrdofm,kk)
!
!  .........end of loop through components
            enddo
!
!  .......end of loop through element trial functions
          enddo
!
ENDIF
!
       case('discon')
!
!  .......loop through element trial functions
          do k=1,nrdoflQ
!
!  .........loop through components
            do ivar=1,NR_COMP(iphys)
              ll = nrdofmHEV + (k-1)*NRQVAR+ADRES(iphys)+ivar
              kk = (k-1)*NR_COMP(iphys)+ivar
!
!  ...........accumulate
              ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
              + AAUX(iphys)%array(1:Nrdofm,kk)
!
!  .........end of loop through components
            enddo
!
!  .......end of loop through element trial functions
          enddo
!
        end select
!
!  ...end of loop through physical attributes
      enddo
!
!-----------------------------------------------------------------------
!
!  ...account for Dirichlet BCs, modify the load vector
!
!  ...loop through Dirichlet trial functions
      do k2=1,Nrdofm
        if (IDBC(k2).eq.1) then
!
!  .......loop through all test functions
          do k1=1,Nrdofm
!
!  .........modify the load vector
            ZBMOD(k1,1:NR_RHS) = ZBMOD(k1,1:NR_RHS)  &
                               - ZAMOD(k1,k2)*ZDOFD(k2,1:NR_RHS)
          enddo
        endif
      enddo
!
!-----------------------------------------------------------------------
!
!  ...compress the element matrices
      do l1=1,Nrdofc
        k1=NEXTRACT(l1)
        do load=1,NR_RHS
          ll1 = (load-1)*Nrdofc+l1
          Zbload(ll1) = ZBMOD(k1,load)
        enddo
        select case(ISYM_FLAG)
!
!  .....symmetric case
        case(1)
          do l2=1,l1
            k2=NEXTRACT(l2)
            k = (l1-1)*l1/2 + l2
            Zastif(k) = (ZAMOD(k1,k2)+ZAMOD(k2,k1))/2.d0
          enddo
!
!  .....unsymmetric case: row-major
        case(2)
          do l2=1,Nrdofc
            k2=NEXTRACT(l2)
            k = (l1-1)*Nrdofc + l2
            Zastif(k) = ZAMOD(k1,k2)
          enddo
!
!  .....unsymmetric case: col-major
        case(3)
          do l2=1,Nrdofc
            k2=NEXTRACT(l2)
            k = (l2-1)*Nrdofc + l1
            Zastif(k) = ZAMOD(k1,k2)
          enddo
        endselect
      enddo
!
!-----------------------------------------------------------------------
!     PRINTING STATEMENTS                                              |
!-----------------------------------------------------------------------
      if (Mdle.eq.1) then
        iprint=0
      else
        iprint=0
      endif
!
      if (iprint.ge.2) then
 20     continue
        write(*,6666) NEXTRACT(1:Nrdofc)
 6666   format(' NEXTRACT = ',20(i3))
        write(*,*)'--------------------------------------------------'
        write(*,*) 'celem: Mdle = ',Mdle
        write(*,*) ' SET LOAD CASE'
        read(*,*) load
        write(*,*) 'SET PHYSICAL ATTRIBUTE'
        read(*,*) iphys1
        write(*,8001) iphys1,load
 8001   format(' LOCAL LOAD VECTOR FOR ATTRIBUTE = ',i2,' AND LOAD = ',i2)
        do ivar1=1,NR_COMP(iphys1)
          write(*,*) 'COMPONENT = ',ivar1
          write(*,8000) ( BLOC(iphys1)%array((k1-1)*NR_COMP(iphys1)+ivar1,    &
                                            load                         ),   &
                          k1=1,Nrdofs(iphys1)/NR_COMP(iphys1)               )
#if HP3D_COMPLEX
 8000     format(24(2e11.4,1x))
#else
 8000     format(16(e12.5,1x))
#endif
        enddo
!
!  .....loop through physical attributes
        do iphys2 = 1, NR_PHYSA
          if (jtrial(iphys2).eq.0) cycle
          write(*,8002) iphys1,iphys2
 8002     format(' LOCAL STIFFNESS MATRIX  FOR ATTRIBUTES = ',2i2)
          do ivar1=1,NR_COMP(iphys1)
          do ivar2=1,NR_COMP(iphys2)
            write(*,*) 'COMPONENTS = ',ivar1,ivar2
            do k1=1,Nrdofs(iphys1)/NR_COMP(iphys1)
              write(*,8000) ( ALOC(iphys1,iphys2)%array((k1-1)*NR_COMP(iphys1)+ivar1,     &
                                                        (k2-1)*NR_COMP(iphys2)+ivar2),    &
                              k2=1,Nrdofs(iphys2)/NR_COMP(iphys2)                      )
            enddo
          enddo
          enddo
        enddo
        write(*,*) 'celem: CONTINUE ? 0-Yes / 1-No'
        read(*,*) ians
        if (ians.eq.1) goto 20
!
 30     write(*,*)'MODIFIED STIFFNESS MATRIX, Nrdofm = ',Nrdofm
        write(*,*) 'ibeg,iend,jbeg,jend (ibeg = 0 - Exit)'
        read(*,*) ibeg,iend,jbeg,jend
        if (ibeg.eq.0) goto 40
        do i=ibeg,iend
!cc          write(*,*) 'i,ZBMOD(i) = ',i,ZBMOD(i,load)
          write(*,8000) ZAMOD(i,jbeg:jend)
        enddo
        goto 30
 40     continue
!
 50     write(*,*)'COMPRESSED STIFFNESS MATRIX, Nrdofc = ',Nrdofc
        select case(ISYM_FLAG)
!  .....SYMMETRIC CASE
        case(1)
          write(*,*)'SYMMETRIC CASE'
          write(*,*)'ibeg,iend,jbeg (ibeg = 0 - Exit)'
          read(*,*) ibeg,iend,jbeg
!  .....SYMMETRIC CASE
        case(2)
          write(*,*)'UNSYMMETRIC CASE'
          write(*,*)'ibeg,iend,jbeg,jend (ibeg = 0 - Exit)'
          read(*,*) ibeg,iend,jbeg,jend
        end select
        if (ibeg.eq.0) goto 60
        do i=ibeg,iend
!cc          write(*,*) 'i,Zbload = ',i,Zbload((load-1)*Nrdofc+i)
          select case(ISYM_FLAG)
!  .......SYMMETRIC CASE
          case(1)
            kbeg=(i-1)*i/2+jbeg; kend=i*(i+1)/2
!  .......UNSYMMETRIC CASE
          case(2)
            kbeg=(i-1)*Nrdofc+jbeg; kend=(i-1)*Nrdofc+jend
          endselect
          write(*,8000) Zastif(kbeg:kend)
        enddo
        goto 50
 60     continue

      endif
!
!
   end subroutine celem_sys
