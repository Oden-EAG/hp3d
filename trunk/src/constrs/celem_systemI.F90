!-----------------------------------------------------------------------
!
!   routine name       - celem_systemI
!
!-----------------------------------------------------------------------
!
!   latest revision    - Sept 18
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
!             Iel      - global element id (used in stc)
!             Mdle     - middle node of an element
!             Idec     = 1 return only the info on nodes
!                      = 2 compute also the element matrices
!     out:
!             Nrdofs   - number of element local dof for each physical
!                        attribute (needed to dimension work arrays)
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
#include "typedefs.h"
subroutine celem_systemI(Iel,Mdle,Idec,                            &
                         Nrdofs,Nrdofm,Nrdofc,                     &
                         Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,  &
                         Zbload,Zastif)
!
   use data_structure3D
   use assembly
   use control
   use stc, only: stc_fwd_wrapper
!
   implicit none
!
   integer,                     intent(in)  :: Iel,Mdle,Idec
   integer,                     intent(out) :: Nrdofs(NR_PHYSA)
   integer,                     intent(out) :: Nrdofm,Nrdofc,Nrnodm
   integer, dimension(MAXNODM), intent(out) :: Nodm,NdofmH,NdofmE,NdofmV,NdofmQ
   VTYPE,                       intent(out) :: Zbload(MAXDOFM)
   VTYPE,                       intent(out) :: Zastif(MAXDOFM**2)
!
!..element type
   integer :: ntype
!
!..element order
   integer :: norder(19),norderi(19)
!
!..index for a node
   integer :: index(NRINDEX)
!
!..number of single H1 and H(curl) or H(div) dof for a higher order node
   integer, dimension(MAXNODM) :: ndofmHl,ndofmEl,ndofmVl
!
!..location of first dof for a modified element node and discretization type
   integer, dimension(MAXNODM) :: namH,namE,namV
!
!..constrained approximation coefficients for a scalar-valued element
   integer :: nrconH(MAXbrickH),nacH(NACDIM,MAXbrickH),  &
              nrconE(MAXbrickE),nacE(NACDIM,MAXbrickE),  &
              nrconV(MAXbrickV),nacV(NACDIM,MAXbrickV)
   real(8) :: constrH(NACDIM,MAXbrickH),  &
              constrE(NACDIM,MAXbrickE),  &
              constrV(NACDIM,MAXbrickV)
!
!..test and trial flags
   integer, dimension(NR_PHYSA) :: itest,jtrial
!
!..auxiliary variables
   integer :: i,ii,il,j,k,k1,k2,kk,kold,kp,l,l1,l2,ll,ll1,load,nvoid
   integer :: iphys,iphys1,iphys2,icomp,iload,ivar
   integer :: nod,nvarHt,nvarEt,nvarVt
   integer :: nrPhysH,nrPhysE,nrPhysV,nrPhysHE,nrPhysHEV
   integer :: nrdofmH,nrdofmE,nrdofmV,nrdofmHE,nrdofmHEV
   integer :: nrdoflH,nrdoflE,nrdoflV,nrdoflQ
   integer :: nrdoflHi,nrdoflEi,nrdoflVi
!
#if DEBUG_MODE
   integer :: ians,ibeg,iend,jbeg,jend,kbeg,kend,ivar1,ivar2
   integer :: iprint
   iprint=0
#endif
!
!----------------------------------------------------------------------
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
      norderi = norder
      ntype = NODES(Mdle)%ntype
      select case(ntype)
        case(MDLB); norderi(19) = 111
        case(MDLP); norderi(15) = 11
        case(MDLN); norderi(11) = 1
        case(MDLD); norderi(14) = 1
      end select
      call celndof(ntype,norder, nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
      call celndof(ntype,norderi, nrdoflHi,nrdoflEi,nrdoflVi,nvoid)
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,5001) Mdle
 5001   format(' celem_systemI: DEBUGGING for Mdle = ',i8)
        call print_order(ntype,norder)
      endif
#endif
!
!  ...number of dofs for ALL H1,H(curl),H(div),L2 components (EXPANDED mode)
      do i=1,NR_PHYSA
        if (PHYSAi(i)) then
          select case(DTYPE(i))
            case('contin') ; Nrdofs(i)=nrdoflHi*NR_COMP(i)
            case('tangen') ; Nrdofs(i)=nrdoflEi*NR_COMP(i)
            case('normal') ; Nrdofs(i)=nrdoflVi*NR_COMP(i)
          endselect
        else
          select case(DTYPE(i))
            case('contin') ; Nrdofs(i)=nrdoflH*NR_COMP(i)
            case('tangen') ; Nrdofs(i)=nrdoflE*NR_COMP(i)
            case('normal') ; Nrdofs(i)=nrdoflV*NR_COMP(i)
            case('discon') ; Nrdofs(i)=nrdoflQ*NR_COMP(i)
          endselect
        endif
      enddo
!
!  ...determine the constraints coefficients
      call logic(Mdle,Idec,                             &
                 Nodm,ndofmHl,ndofmEl,ndofmVl,Nrnodm,   &
                 nrconH,nacH,constrH,                   &
                 nrconE,nacE,constrE,                   &
                 nrconV,nacV,constrV)
!
!  ...eliminate the middle node from the list of modified element nodes
      Nrnodm = Nrnodm-1
!
!  ...printing
#if DEBUG_MODE
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
#endif
!
!  ...number of H1 dofs for the modified element w/o middle node (EXPANDED mode)
      k=0
      do i=1,Nrnodm
        namH(i) = k
        k = k + ndofmHl(i)
      enddo
      nrdofmH = k*NRHVAR
!
!  ...number of H(curl) dofs for the modified element w/o middle node  (EXPANDED mode)
      k=0
      do i=1,Nrnodm
        namE(i) = k
        k = k + ndofmEl(i)
      enddo
      nrdofmE = k*NREVAR
!
!  ...number of H(div) dofs for the modified element w/o middle node (EXPANDED mode)
      k=0
      do i=1,Nrnodm
        namV(i) = k
        k = k + ndofmVl(i)
      enddo
      nrdofmV = k*NRVVAR
!
      nrdofmHE  = nrdofmH   + nrdofmE
      nrdofmHEV = nrdofmHE  + nrdofmV
      Nrdofm    = nrdofmHEV
!
#if DEBUG_MODE
      if (iprint .eq. 1) then
        write(*,6001) nrdoflH,nrdoflE,nrdoflV,nrdoflQ
 6001   format('celem_systemI: nrdoflH,nrdoflE,nrdoflV,nrdoflQ = ',4i4)
        write(*,6002) nrdoflHi,nrdoflEi,nrdoflVi
 6002   format('celem_systemI: nrdoflHi,nrdoflEi,nrdoflVi      = ',3i4)
       write(*,6003) nrdofmH,nrdofmHE,nrdofmHEV,Nrdofm
 6003   format('celem_systemI: nrdofmH,nrdofmHE,nrdofmHEV,Nrdofm = ',4i4)
      endif
#endif
!
!***********************************************************************
!
!..establish extraction vector and Dirichlet data
   if (Idec.eq.2) then
      ZDOFD = ZERO
      IDBC  = 0
   endif
!
!..count number of variables (not components) for each physics type
   nrPhysH=0; nrPhysE=0; nrPhysV=0
   do iphys=1,NR_PHYSA
      select case(DTYPE(iphys))
         case('contin')
            nrPhysH=nrPhysH+1
         case('tangen')
            nrPhysE=nrPhysE+1
         case('normal')
            nrPhysV=nrPhysV+1
         case default
      end select
   enddo
#if DEBUG_MODE
   if (iprint .eq. 1) then
      write(*,*) 'Physics variables:'
      write(*,*) 'nrPhysH = ', nrPhysH
      write(*,*) 'nrPhysE = ', nrPhysE
      write(*,*) 'nrPhysV = ', nrPhysV
   endif
#endif
   nrPhysHE  = nrPhysH+nrPhysE
   nrPhysHEV = nrPhysH+nrPhysE+nrPhysV
!
!..initiate global dof counter for the compressed mode
   k=0
!
!---------------------------- H1 dof -----------------------------------
!
!..skip H1 dof if no H1 physics present
   if (nrPhysH .eq. 0) goto 100
!
!..loop through nodes in the reversed order
   do i=Nrnodm,1,-1
      nod = Nodm(i)
!  ...total number of H1 components for a single load
      nvarHt = NREQNH(NODES(nod)%case)
      l = namH(i)*NRHVAR
      call get_index(nod, index)
      kold = k
!
!  ...loop through H1 dof for the node
      do j=1,ndofmHl(i)
         ivar=0; ii=0
!
!     ...loop through physics variables
         do iphys=1,nrPhysH
            il = NR_COMP(iphys)
            if (index(ii+1) .eq. 0) then
               l=l+il; ii=ii+il
               cycle
            endif
            if (.not. PHYSAm(iphys)) then
               l=l+il; ii=ii+il; ivar=ivar+il
               cycle
            endif
!
!        ...loop through components of this physics variable
            do icomp=1,il
               l=l+1; ii=ii+1; ivar=ivar+1
               select case(index(ii))
!
!           ...dof present but known from Dirichlet BC, save the BC data
               case(1)
                  if (Idec.eq.2) then
                     IDBC(l)=1
                     do iload=1,NR_RHS
                        ZDOFD(l,iload) = NODES(nod)%dof%zdofH((iload-1)*nvarHt+ivar,j)
                     enddo
                  endif
!
!           ...dof present and active
               case(2)
                  k=k+1
                  if (Idec.eq.2) NEXTRACT(k) = l
#if DEBUG_MODE
               case default
                  write(*,*) 'celem_systemI: nod,ii,index(ii) = ',nod,ii,index(ii)
                  write(*,*) '                    in H1: Mdle = ', Mdle, '. stop.'
                  stop
#endif
               end select
!        ...end loop through components
            enddo
!     ...end loop through H1 physics variables
         enddo
!  ...end loop through H1 dof for node
      enddo
!
!  ..store the actual number of dof
     NdofmH(i) = k-kold
!
!..end of loop through nodes
   enddo
!
!---------------------------- H(curl) dof ------------------------------
!
!  ...skip H(curl) dof if no H(curl) physics present
  100 if (nrPhysE .eq. 0) goto 200
!
!..loop through nodes in the reversed order
   do i=Nrnodm,1,-1
      nod = Nodm(i)
!  ...total number of H(curl) components for a single load
      nvarEt = NREQNE(NODES(nod)%case)
      l = nrdofmH + namE(i)*NREVAR
      call get_index(nod, index)
      kold = k
!
!  ...loop through H(curl) dof for the node
      do j=1,ndofmEl(i)
         ivar=0; ii=NRHVAR
!
!     ...loop through physics variables
         do iphys=nrPhysH+1,nrPhysHE
            il = NR_COMP(iphys)
            if (index(ii+1) .eq. 0) then
               l=l+il; ii=ii+il
               cycle
            endif
            if (.not. PHYSAm(iphys)) then
               l=l+il; ii=ii+il; ivar=ivar+il
               cycle
            endif
!
!        ...loop through components of this physics variable
            do icomp=1,il
               l=l+1; ii=ii+1; ivar=ivar+1
               select case(index(ii))
!
!           ...dof present but known from Dirichlet BC, save the BC data
               case(3)
                  if (Idec.eq.2) then
                     IDBC(l)=1
                     do iload=1,NR_RHS
                        ZDOFD(l,iload) = NODES(nod)%dof%zdofE((iload-1)*nvarEt+ivar,j)
                     enddo
                  endif
!
!           ...dof present and active
               case(4)
                  k=k+1
                  if (Idec.eq.2) NEXTRACT(k) = l
#if DEBUG_MODE
               case default
                  write(*,*) 'celem_systemI: nod,ii,index(ii) = ',nod,ii,index(ii)
                  write(*,*) '               in H(curl): Mdle = ', Mdle, '. stop.'
                  stop
#endif
               end select
!        ...end loop through components
            enddo
!     ...end loop through H(curl) physics variables
         enddo
!  ...end loop through H(curl) dof for node
      enddo
!
!  ...compute the actual number of dof
      NdofmE(i) = k-kold
!
!..end of loop through nodes
   enddo
!
!---------------------------- H(div) dof -------------------------------
!
!  ...skip Hdiv if no H(div) physics present
  200 if (nrPhysV .eq. 0) goto 300
!
!..loop through nodes in the reversed order
   do i=Nrnodm,1,-1
      nod = Nodm(i)
!  ...total number of H(div) components for a single load
      nvarVt = NREQNV(NODES(nod)%case)
      l = nrdofmHE + namV(i)*NRVVAR
      call get_index(nod, index)
      kold = k
!
!  ...loop through H(curl) dof for the node
      do j=1,ndofmVl(i)
         ivar=0; ii=NRHVAR+NREVAR
!
!     ...loop through physics variables
         do iphys=nrPhysHE+1,nrPhysHEV
            il = NR_COMP(iphys)
            if (index(ii+1) .eq. 0) then
               l=l+il; ii=ii+il
               cycle
            endif
            if (.not. PHYSAm(iphys)) then
               l=l+il; ii=ii+il; ivar=ivar+il
               cycle
            endif
!
!        ...loop through components of this physics variable
            do icomp=1,il
               l=l+1; ii=ii+1; ivar=ivar+1
               select case(index(ii))
!
!           ...dof present but known from Dirichlet BC, save the BC data
               case(5)
                  if (Idec.eq.2) then
                     IDBC(l)=1
                     do iload=1,NR_RHS
                        ZDOFD(l,iload) = NODES(nod)%dof%zdofV((iload-1)*nvarVt+ivar,j)
                     enddo
                  endif
!
!           ...dof present and active
               case(6)
                  k=k+1
                  if (Idec.eq.2) NEXTRACT(k) = l
#if DEBUG_MODE
               case default
                  write(*,*) 'celem_systemI: nod,ii,index(ii) = ',nod,ii,index(ii)
                  write(*,*) '                in H(div): Mdle = ', Mdle, '. stop.'
                  stop
#endif
               end select
!        ...end loop through components
            enddo
!     ...end loop through H(div) physics variables
         enddo
!  ...end loop through H(div) dof for node
      enddo
!
!  ...compute the actual number of dof
      NdofmV(i) = k-kold
!
!..end of loop through nodes
   enddo
!
!-----------------------------------------------------------------------
!
  300 continue
!
!..save number of dof for the modified element in the compressed mode
   Nrdofc = k
!
!  ...printing
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7002) (Nodm(i),NdofmH(i),NdofmE(i),NdofmV(i), i=1,Nrnodm)
 7002   format('       Nodm,NdofmH,NdofmE,NdofmV = ',         &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x),           &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x),           &
               /,7x,4(i8,4i4,2x),/,7x,4(i8,4i4,2x))
        write(*,7003) Nrdofc
 7003   format('celem_systemI: Nrdofc = ',i4)
        if (Idec .eq. 2) then
          write(*,7004) NEXTRACT(1:Nrdofc)
 7004     format('celem_systemI: NEXTRACT = ',10(/,20i5))
          write(*,7005) IDBC(1:Nrdofm)
 7005     format('celem_systemI: IDBC = ',10(4i2,1x),                   &
             10(/,'              ',10(4i2,1x)))
          call pause
        endif
      endif
#endif
!
!..return if only INFO is needed
   if (Idec.eq.1) return
!
!***********************************************************************
!
!***********************************************************************
!
!..determine the element local load vector and stiffness matrix
   call elem(Mdle, itest,jtrial)
#if DEBUG_MODE
   do i=1,NR_PHYSA
     if((.not. PHYSAm(i)) .and. ((itest(i).eq.1) .or. (jtrial(i).eq.1))) then
       write(*,*) 'celem_systemI: INCONSISTENCY of PHYSAm and itest/jtrial. stop.'
       stop
     endif
   enddo
#endif
!
!..perform static condensation of middle node dof
   call stc_fwd_wrapper(Iel,Mdle)
!
!***********************************************************************
!
!***********************************************************************
!
!..initiate modified load vector and half-way modified
!  auxiliary stiffness matrices
   ZBMOD(1:Nrdofm,1:NR_RHS) = ZERO
   do i=1,NR_PHYSA
      AAUX(i)%array(1:Nrdofm,1:Nrdofs(i)) = ZERO
   enddo
!
!  ...compute the modified element matrices in the expanded form
!
!..first loop through physical attributes
   do iphys1=1,NR_PHYSA
!
!  ...skip if the attribute is absent
      if (itest(iphys1).eq.0) cycle
!
      select case(DTYPE(iphys1))
      case('contin')
!
!     ...loop through element test functions excluding middle node
         do k=1,nrdoflHi
!
!        ...loop through the connected real dofs
            do kp=1,nrconH(k)
!
!           ...modified matrix location
               l = nacH(kp,k)
!
!           ...loop through components
               do ivar=1,NR_COMP(iphys1)
                  ll = (l-1)*NRHVAR+ADRES(iphys1)+ivar
                  kk = (k-1)*NR_COMP(iphys1)+ivar
                  ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)              &
                     + BLOC(iphys1)%array(kk,1:NR_RHS)*constrH(kp,k)
!
!              ...second loop through physical attributes
                  do iphys2=1,NR_PHYSA
!
!                 ...skip if the attribute is absent
                     if (jtrial(iphys2).eq.0) cycle
!
!                 ...transform the stiffness matrix
                     AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                     = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                     + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))  &
                     *constrH(kp,k)
!
!              ...end of second loop through physical attributes
                  enddo
!
!           ...end of loop through components
               enddo
!
!        ...end of loop through connected parent dof
            enddo
!
!     ...end of loop through element test functions
         enddo
!
      case('tangen')
!
!     ...loop through element test functions excluding middle node
         do k=1,nrdoflEi
!
!        ...loop through the connected real dofs
            do kp=1,nrconE(k)
!
!           ...modified matrix location
               l = nacE(kp,k)
!
!           ...loop through components
               do ivar = 1, NR_COMP(iphys1)
                  ll = nrdofmH + (l - 1)*NREVAR + ADRES(iphys1) + ivar
                  kk = (k - 1)*NR_COMP(iphys1) + ivar
                  ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)              &
                     + BLOC(iphys1)%array(kk,1:NR_RHS)*constrE(kp,k)
!
!              ...second loop through physical attributes
                  do iphys2=1,NR_PHYSA
!
!                 ...skip if the attribute is absent
                     if (jtrial(iphys2).eq.0) cycle
!
!                 ...transform the stiffness matrix
                     AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                     = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                     + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))  &
                     *constrE(kp,k)
!
!              ...end of second loop through physical attributes
                  enddo
!
!           ...end of loop through components
               enddo
!
!        ...end of loop through connected parent dof
            enddo
!
!     ...end of loop through element test functions
         enddo
!
      case('normal')
!
!     ...loop through element test functions excluding middle node
         do k=1,nrdoflVi
!
!        ...loop through the connected real dofs
            do kp=1,nrconV(k)
!
!           ...modified matrix location
               l = nacV(kp,k)
!
!           ...loop through components
               do ivar=1,NR_COMP(iphys1)
                  ll = nrdofmHE + (l-1)*NRVVAR+ADRES(iphys1)+ivar
                  kk = (k-1)*NR_COMP(iphys1)+ivar
                  ZBMOD(ll,1:NR_RHS) = ZBMOD(ll,1:NR_RHS)              &
                     + BLOC(iphys1)%array(kk,1:NR_RHS)*constrV(kp,k)
!
!              ...second loop through physical attributes
                  do iphys2=1,NR_PHYSA
!
!                 ...skip if the attribute is absent
                     if (jtrial(iphys2).eq.0) cycle
!
!                 ...transform the stiffness matrix
                     AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))           &
                     = AAUX(iphys2)%array(ll,1:Nrdofs(iphys2))         &
                     + ALOC(iphys1,iphys2)%array(kk,1:Nrdofs(iphys2))  &
                     *constrV(kp,k)
!
!              ...end of second loop through physical attributes
                  enddo
!
!           ...end of loop through components
               enddo
!
!        ...end of loop through connected parent dof
            enddo
!
!     ...end of loop through element test functions
         enddo
!
      end select
!
!..end of first loop through physical attributes
   enddo
!
!
!***********************************************************************
!
!..initiate modified stiffness matrix
   ZAMOD(1:Nrdofm,1:Nrdofm) = ZERO
!
!..loop through physical attributes
   do iphys=1,NR_PHYSA
!
!  ...skip if the attribute is absent
      if (jtrial(iphys).eq.0) cycle
!
      select case(DTYPE(iphys))
      case('contin')
!
!     ...loop through element trial functions
         do k=1,nrdoflHi
!
!        ...loop through the connected real dofs
            do kp=1,nrconH(k)
!
!           ...modified matrix location
               l = nacH(kp,k)
!
!           ...loop through components
               do ivar=1,NR_COMP(iphys)
                  ll = (l-1)*NRHVAR+ADRES(iphys)+ivar
                  kk = (k-1)*NR_COMP(iphys)+ivar
!
!              ...accumulate
                  ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
                  + AAUX(iphys)%array(1:Nrdofm,kk)*constrH(kp,k)
!
!           ...end of loop through components
               enddo
!
!        ...end of loop through connected parent dof
            enddo
!
!     ...end of loop through element trial functions
         enddo
!
      case('tangen')
!
!     ...loop through element trial functions
         do k=1,nrdoflEi
!
!        ...loop through the connected real dofs
            do kp=1,nrconE(k)
!
!           ...modified matrix location
               l = nacE(kp,k)
!
!           ...loop through components
               do ivar=1,NR_COMP(iphys)
                  ll = nrdofmH + (l-1)*NREVAR+ADRES(iphys)+ivar
                  kk = (k-1)*NR_COMP(iphys)+ivar
!
!              ...accumulate
                  ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
                  + AAUX(iphys)%array(1:Nrdofm,kk)*constrE(kp,k)
!
!           ...end of loop through components
               enddo
!
!        ...end of loop through connected parent dof
            enddo
!
!     ...end of loop through element trial functions
         enddo
!
      case('normal')
!
!     ...loop through element trial functions
         do k=1,nrdoflVi
!
!        ...loop through the connected real dofs
            do kp=1,nrconV(k)
!
!           ...modified matrix location
               l = nacV(kp,k)
!
!           ...loop through components
               do ivar=1,NR_COMP(iphys)
                  ll = nrdofmHE + (l-1)*NRVVAR+ADRES(iphys)+ivar
                  kk = (k-1)*NR_COMP(iphys)+ivar
!
!              ...accumulate
                  ZAMOD(1:Nrdofm,ll) = ZAMOD(1:Nrdofm,ll) &
                  + AAUX(iphys)%array(1:Nrdofm,kk)*constrV(kp,k)
!
!           ...end of loop through components
               enddo
!
!        ...end of loop through connected parent dof
            enddo
!
!     ...end of loop through element trial functions
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
!..loop through Dirichlet trial functions
   do k2=1,Nrdofm
      if (IDBC(k2).eq.1) then
!
!     ...loop through all test functions
         do k1=1,Nrdofm
!
!        ...modify the load vector
            ZBMOD(k1,1:NR_RHS) = ZBMOD(k1,1:NR_RHS)  &
                               - ZAMOD(k1,k2)*ZDOFD(k2,1:NR_RHS)
         enddo
      endif
   enddo
!
!-----------------------------------------------------------------------
!
!..compress the element matrices
   do l1=1,Nrdofc
      k1=NEXTRACT(l1)
      do load=1,NR_RHS
         ll1 = (load-1)*Nrdofc+l1
         Zbload(ll1) = ZBMOD(k1,load)
      enddo
      select case(ISYM_FLAG)
!
!  ...symmetric case
      case(1)
         do l2=1,l1
            k2=NEXTRACT(l2)
            k = (l1-1)*l1/2 + l2
            Zastif(k) = (ZAMOD(k1,k2)+ZAMOD(k2,k1))/2.d0
         enddo
!
!  ...unsymmetric case: row-major
      case(2)
         do l2=1,Nrdofc
            k2=NEXTRACT(l2)
            k = (l1-1)*Nrdofc + l2
            Zastif(k) = ZAMOD(k1,k2)
         enddo
!
!  ...unsymmetric case: col-major
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
#if DEBUG_MODE
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
        write(*,*) 'celem_systemI: Mdle = ',Mdle
        write(*,*) ' SET LOAD CASE'
        read(*,*) load
        write(*,*) 'SET PHYSICAL ATTRIBUTE'
        read(*,*) iphys1
        if (DTYPE(iphys1).eq.'discon') goto 20
        write(*,8001) iphys1,load
 8001   format(' LOCAL LOAD VECTOR FOR ATTRIBUTE = ',i2,' AND LOAD = ',i2)
        do ivar1=1,NR_COMP(iphys1)
          write(*,*) 'COMPONENT = ',ivar1
          write(*,8000) ( BLOC(iphys1)%array((k1-1)*NR_COMP(iphys1)+ivar1,    &
                                            load                         ),   &
                          k1=1,Nrdofs(iphys1)/NR_COMP(iphys1)               )
#if C_MODE
 8000     format(24(2e11.4,1x))
#else
 8000     format(16(e12.5,1x))
#endif
        enddo
!
!  .....loop through physical attributes
        do iphys2 = 1, NR_PHYSA
          if (DTYPE(iphys2).eq.'discon') cycle
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
        write(*,*) 'celem_systemI: CONTINUE ? 0-Yes / 1-No'
        read(*,*) ians
        if (ians.eq.1) goto 20
!
 30     write(*,*)'MODIFIED STIFFNESS MATRIX, Nrdofm = ',Nrdofm
        write(*,*) 'ibeg,iend,jbeg,jend (ibeg = 0 - Exit)'
        read(*,*) ibeg,iend,jbeg,jend
        if (ibeg.eq.0) goto 40
        do i=ibeg,iend
!         write(*,*) 'i,ZBMOD(i) = ',i,ZBMOD(i,load)
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
        case default
          write(*,*) ' celem_systemI: ISYM_FLAG = ',ISYM_FLAG
        end select
        if (ibeg.eq.0) goto 60
        do i=ibeg,iend
!         write(*,*) 'i,Zbload = ',i,Zbload((load-1)*Nrdofc+i)
          select case(ISYM_FLAG)
!  .......SYMMETRIC CASE
          case(1)
            kbeg=(i-1)*i/2+jbeg; kend=i*(i+1)/2
!  .......UNSYMMETRIC CASE
          case(2)
            kbeg=(i-1)*Nrdofc+jbeg; kend=(i-1)*Nrdofc+jend
          case default
            kbeg=0; kend = 0;
            write(*,*) ' celem_systemI: ISYM_FLAG = ',ISYM_FLAG
          endselect
          write(*,8000) Zastif(kbeg:kend)
        enddo
        goto 50
 60     continue
      endif
#endif
!
!
end subroutine celem_systemI
