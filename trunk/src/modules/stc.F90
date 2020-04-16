!
#include "implicit_none.h"
!
!-----------------------------------------------------------------------
!
!    module name        - stc
!
!-----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - module sets up required workspace for the
!                         static condensation, and it contains routines
!                         for computing the condensed system and the
!                         Schur complements, as well as to perform
!                         backsubstitution to solve for bubble dof.
!
!    contains
!        subroutines:   - stc_alloc
!                       - stc_dealloc
!                       - stc_get_nrdof
!                       - stc_fwd_wrapper
!                       - stc_fwd_herm
!                       - stc_fwd_gen
!                       - stc_bwd_wrapper
!                       - stc_bwd
!                       - stc_solout
!
!----------------------------------------------------------------------!
module stc
!
   use assembly, only: NR_RHS
   use data_structure3D
   use par_mesh, only: DISTRIBUTED
!
!..storing Schur complements
   logical, save :: STORE_STC   = .true.
!
!..optimization flag: real symmetric/complex Hermitian element matrices
   logical, save :: HERM_STC    = .true.
!
!..workspace for storing element matrices
   type super_mat
      integer :: mdle
      integer :: ni
!
!  ...Schur complement of matrix A
      VTYPE, allocatable :: ASchur(:,:)
!  ...and the corresponding right-hand side B
      VTYPE, allocatable :: BSchur(:,:)
!  ...local connectivity map
      integer, allocatable :: con(:)
   endtype super_mat
!
!..element workspace for static condensation
   type(super_mat), allocatable :: CLOC(:)
!
contains
!
!
subroutine stc_alloc
   if (.not. DISTRIBUTED) NRELES_SUBD = NRELES
   allocate(CLOC(NRELES_SUBD))
end subroutine stc_alloc
!
!
subroutine stc_dealloc
   if (allocated(CLOC)) deallocate(CLOC)
end subroutine stc_dealloc
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_get_nrdof
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine computes the number of interface and
!                         bubble dof for each physics variable of this 
!                         element. It accounts only for variables that
!                         are supported on the middle node of the element.
!
!    arguments :
!      in:
!              Mdle         - Global middle node id (used in NODES)
!     out:
!              Nrdofi       - Number of interface dofs
!              Nrdofb       - Number of bubble dofs
!
!----------------------------------------------------------------------
subroutine stc_get_nrdof(Mdle, Nrdofi,Nrdofb)
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Nrdofi(NR_PHYSA),Nrdofb(NR_PHYSA)
!
!..element type
   character(len=4) :: etype
!
!..element order
   integer :: nord,norderi(19)
!
   integer :: index(NRINDEX),ncase(NR_PHYSA)
!
   integer :: nrdofHmdl,nrdofEmdl,nrdofVmdl,nrdofQmdl
   integer :: nrdoflHi,nrdoflEi,nrdoflVi
!
   integer :: i,nvoid
!
!----------------------------------------------------------------------
!
   Nrdofi(1:NR_PHYSA)=0
   Nrdofb(1:NR_PHYSA)=0
!
!..set node to be middle node number (global id in NODES)
   call decod(NODES(Mdle)%case,2,NR_PHYSA, ncase)
!
!..number of dofs for a SINGLE H1,H(curl),H(div),L2 component
   call find_order(Mdle, norderi)
   etype = NODES(Mdle)%type
   select case(etype)
      case('mdlb'); nord = norderi(19); norderi(19) = 111
      case('mdlp'); nord = norderi(15); norderi(15) = 11
      case('mdln'); nord = norderi(11); norderi(11) = 1
      case('mdld'); nord = norderi(14); norderi(14) = 1
   end select
   call celndof(etype,norderi, nrdoflHi,nrdoflEi,nrdoflVi,nvoid)
   call ndof_nod(etype,nord, nrdofHmdl,nrdofEmdl,nrdofVmdl,nrdofQmdl)
!
!..number of dofs for ALL H1,H(curl),H(div),L2 components
   do i=1,NR_PHYSA
!  ...stop if variable is not supported on middle node
      if (ncase(i) .eq. 0) then
         write(*,*) 'stc_get_nrdof: static condensation requires ', &
                    'elements to support all physics variables. stop.'
         stop
      endif
!  ...skip if variable is deactivated
      if (.not. PHYSAm(i)) cycle
!
      select case(DTYPE(i))
         case('contin') ; Nrdofi(i)=nrdoflHi*NR_COMP(i)
         case('tangen') ; Nrdofi(i)=nrdoflEi*NR_COMP(i)
         case('normal') ; Nrdofi(i)=nrdoflVi*NR_COMP(i)
      end select
      if (.not. PHYSAi(i)) then
         select case(DTYPE(i))
            case('contin') ; Nrdofb(i)=nrdofHmdl*NR_COMP(i)
            case('tangen') ; Nrdofb(i)=nrdofEmdl*NR_COMP(i)
            case('normal') ; Nrdofb(i)=nrdofVmdl*NR_COMP(i)
            case('discon') ; Nrdofb(i)=nrdofQmdl*NR_COMP(i)
         end select
      endif
   enddo
!
end subroutine stc_get_nrdof
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_fwd_wrapper
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine performs static condensation for an
!                         element based on the ALOC, BLOC arrays.
!                         It writes the condensed system back into
!                         corresponding ALOC, BLOC locations.
!
!    arguments :
!      in:
!              Iel      - Global element id (used in CLOC)
!              Mdle     - Global middle node id (used in NODES)
!
!----------------------------------------------------------------------
subroutine stc_fwd_wrapper(Iel,Mdle)
!
   use assembly, only: ALOC, BLOC
!
   implicit none
!
   integer, intent(in) :: Iel,Mdle
!
!..workspace for interface and bubble system
   VTYPE, allocatable :: Aii(:,:),Abi(:,:),Aib(:,:),Abb(:,:),Bi(:,:),Bb(:,:)
!
!..number of interface and bubble dof for each variable
   integer :: nrdofi(NR_PHYSA),nrdofb(NR_PHYSA)
!
!..aux variables
   integer :: i,j,ii,ib,ji,jb,ki,kb,ni,nb
!
!--------------------------------------------------------------
!
   call stc_get_nrdof(Mdle, nrdofi,nrdofb)
!
   CLOC(Iel)%mdle = Mdle
!
!..initialize interface and bubble dof counter
   ni = 0; nb = 0
   do j=1,NR_PHYSA
      ni = ni + nrdofi(j)
      nb = nb + nrdofb(j)
   enddo
!
!..nothing to do if there are no bubble dofs
   if (nb .eq. 0) goto 500
   if (ni .eq. 0) then
      write(*,*) 'stc_fwd_wrapper: no interface dof. stop.'
      stop
   endif
!
   allocate(Aii(ni,ni))
   allocate(Abi(nb,ni))
   allocate(Aib(ni,nb))
   allocate(Abb(nb,nb))
   allocate(Bi(ni,NR_RHS))
   allocate(Bb(nb,NR_RHS))
!
!..construct interface and bubble system
   ni = 0; nb = 0
   do j=1,NR_PHYSA
      ki = 0; kb = 0
      ji = nrdofi(j)
      jb = nrdofb(j)
      do i=1,NR_PHYSA
!     ...construct interface and bubble stiffness
         ii = nrdofi(i)
         ib = nrdofb(i)
         if (ii > 0 .and. ji > 0) then
            Aii(ki+1:ki+ii,ni+1:ni+ji) = ALOC(i,j)%array(1:ii,1:ji)
         endif
         if (ib > 0 .and. ji > 0) then
            Abi(kb+1:kb+ib,ni+1:ni+ji) = ALOC(i,j)%array(ii+1:ii+ib,1:ji)
         endif
         if (ii > 0 .and. jb > 0) then
            Aib(ki+1:ki+ii,nb+1:nb+jb) = ALOC(i,j)%array(1:ii,ji+1:ji+jb)
         endif
         if (ib > 0 .and. jb > 0) then
            Abb(kb+1:kb+ib,nb+1:nb+jb) = ALOC(i,j)%array(ii+1:ii+ib,ji+1:ji+jb)
         endif
         ki = ki + ii
         kb = kb + ib
      enddo
!
!  ...construct interface and bubble load
      if (ji > 0) then
         Bi(ni+1:ni+ji,1:NR_RHS) = BLOC(j)%array(1:ji,1:NR_RHS)
      endif
      if (jb > 0) then
         Bb(nb+1:nb+jb,1:NR_RHS) = BLOC(j)%array(ji+1:ji+jb,1:NR_RHS)
      endif
      ni = ni + ji
      nb = nb + jb
   enddo
!
!..perform static condensation of bubble matrix
   if (HERM_STC) then
      call stc_fwd_herm(ni,nb, Aii,Abi,Aib,Abb,Bi,Bb)
   else
      call stc_fwd_gen(ni,nb, Aii,Abi,Aib,Abb,Bi,Bb)
   endif
!
!..Store Schur complements if STORE_STC
!..   inv(Abb)*Bb
!..   inv(Abb)*Abi
   if (STORE_STC) then
      allocate(CLOC(Iel)%ASchur(nb,ni))
      allocate(CLOC(Iel)%BSchur(nb,NR_RHS))
      CLOC(Iel)%ASchur = Abi
      CLOC(Iel)%BSchur = Bb
   else
      write(*,*) 'stc_fwd_wrapper: ', &
                 'recomputing schur complement not yet implemented. stop.'
      stop
   endif
!
   deallocate(Abi,Aib,Abb,Bb)
!
!..extract condensed system
   ni = 0
   do j=1,NR_PHYSA
      ki = 0
      ji = nrdofi(j)
      do i=1,NR_PHYSA
!     ...extract condensed stiffness into ALOC
         ii = nrdofi(i)
         if (ii > 0 .and. ji > 0) then
            ALOC(i,j)%array(1:ii,1:ji) = Aii(ki+1:ki+ii,ni+1:ni+ji) 
         endif
         ki = ki + ii
      enddo
!
!  ...extract condensed load into BLOC
      if (ji > 0) then
         BLOC(j)%array(1:ji,1:NR_RHS) = Bi(ni+1:ni+ji,1:NR_RHS)
      endif
      ni = ni + ji
   enddo
!
   deallocate(Aii,Bi)
!
  500 continue
!
end subroutine stc_fwd_wrapper
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_fwd_herm
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine performs static condensation for the
!                         real symmetric or complex Hermitian case
!                         (eliminates: H1,H(curl),H(div) and L2 bubbles)
!
!    arguments :
!      in:
!              Ni           - total number of interface dof
!              Nb           - total number of bubble dof
!              Aii          - stiffness: interface-interface
!              Abi          - stiffness: bubble-interface
!              Aib          - stiffness: interface-bubble
!              Abb          - stiffness: bubble-bubble
!              Bi           - load: interface
!              Bb           - load: bubble
!
!----------------------------------------------------------------------
subroutine stc_fwd_herm(Ni,Nb, Aii,Abi,Aib,Abb,Bi,Bb)
!
   implicit none
!
   integer, intent(in)    :: Ni,Nb
   VTYPE  , intent(inout) :: Aii(Ni,Ni),Abi(Nb,Ni),Aib(Ni,Nb),Abb(Nb,Nb), &
                             Bi(Ni,NR_RHS),Bb(Nb,NR_RHS)
!
!..workspace for Hermitian bubble matrix in packed format
   VTYPE, dimension(Nb*(Nb+1)/2) :: AP
!
   integer  :: i,j,k,nk,k1,k2
   integer  :: info
!
!..function for vector storage of a symmetric matrix
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!----------------------------------------------------------------------
!
!..Compute Cholesky factorization: Abb = LL^*
!..first rewrite in packed form
   do i=1,Nb
      do j=i,Nb
         k = nk(i,j)
         AP(k) = Abb(i,j)
      enddo
   enddo
!
!..Cholesky factorization
#if C_MODE
   call ZPPTRF('U', Nb, AP, info)
#else
   call DPPTRF('U', nb, AP, info)
#endif
   if (info.ne.0) then
      write(*,*) 'stc_fwd: PPTRF: info = ', info
      stop
   endif
!
!..Compute inv(Abb)*Bb
#if C_MODE
   call ZPPTRS('U', Nb, NR_RHS, AP, Bb, Nb, info)
#else
   call DPPTRS('U', Nb, NR_RHS, AP, Bb, Nb, info)
#endif
   if (info.ne.0) then
      write(*,*) 'stc_fwd: PPTRS: info = ', info
      stop
   endif
!
!..Compute inv(Abb)*Abi
#if C_MODE
   call ZPPTRS('U', Nb, Ni, AP, Abi, Nb, info)
#else
   call DPPTRS('U', Nb, Ni, AP, Abi, Nb, info)
#endif
   if (info.ne.0) then
      write(*,*) 'stc_fwd: PPTRS: info = ', info
      stop
   endif
!
!..Compute condensed RHS Bi (store it in Bi)
!  Bi = Bi - Aib*inv(Abb)*Bb
!  GEMM: Bi = -1.0 * Aib * (inv(Abb)*Bb) + 1.0 * Bi
#if C_MODE
   call ZGEMM('N','N',Ni,NR_RHS,Nb,-ZONE,Aib,Ni,Bb,Nb,ZONE,Bi,Ni)
#else
   call DGEMM('N','N',Ni,NR_RHS,Nb,-ZONE,Aib,Ni,Bb,Nb,ZONE,Bi,Ni)
#endif
!
!..Compute condensed system Aii (store it in Aii)
!  Aii = Aii - Aib*inv(Abb)*Abi
!  GEMM: Aii = -1.0 * Aib * (inv(Abb)*Abi) + 1.0 * Aii
#if C_MODE
   call ZGEMM('N','N',Ni,Ni,Nb,-ZONE,Aib,Ni,Abi,Nb,ZONE,Aii,Ni)
#else
   call DGEMM('N','N',Ni,Ni,Nb,-ZONE,Aib,Ni,Abi,Nb,ZONE,Aii,Ni)
#endif
!
end subroutine stc_fwd_herm
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_fwd_gen
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine performs static condensation for the
!                         general case (non symmetric/hermitian)
!                         (eliminates: H1,H(curl),H(div) and L2 bubbles)
!
!    arguments :
!      in:
!              Ni           - total number of interface dof
!              Nb           - total number of bubble dof
!              Aii          - stiffness: interface-interface
!              Abi          - stiffness: bubble-interface
!              Aib          - stiffness: interface-bubble
!              Abb          - stiffness: bubble-bubble
!              Bi           - load: interface
!              Bb           - load: bubble
!
!----------------------------------------------------------------------
subroutine stc_fwd_gen(Ni,Nb, Aii,Abi,Aib,Abb,Bi,Bb)
!
   implicit none
!
   integer, intent(in)    :: Ni,Nb
   VTYPE  , intent(inout) :: Aii(Ni,Ni),Abi(Nb,Ni),Aib(Ni,Nb),Abb(Nb,Nb), &
                             Bi(Ni,NR_RHS),Bb(Nb,NR_RHS)
!
!..workspace for pivot matrix
   integer :: piv(Nb)
!
   integer  :: info
!
!----------------------------------------------------------------------
!
!..Compute pivoted LU factorization: Abb = PLU
#if C_MODE
      call ZGETRF(Nb,Nb,Abb,Nb,piv,info)
#else
      call DGETRF(Nb,Nb,Abb,Nb,piv,info)
#endif
   if (info .ne. 0) then
      write(*,*) 'stc_gen: GETRF: info = ', info
      stop
   endif
!
!..Compute inv(Abb)*Bb
#if C_MODE      
      call ZGETRS('N',Nb,NR_RHS,Abb,Nb,piv,Bb,Nb,info)
#else
      call DGETRS('N',Nb,NR_RHS,Abb,Nb,piv,Bb,Nb,info)
#endif      
   if (info.ne.0) then
      write(*,*) 'stc_gen: GETRS: info = ', info
      stop
   endif
!
!..Compute inv(Abb)*Abi
#if C_MODE      
      call ZGETRS('N',Nb,Ni,Abb,Nb,piv,Abi,Nb,info)
#else
      call DGETRS('N',Nb,Ni,Abb,Nb,piv,Abi,Nb,info)
#endif 
   if (info.ne.0) then
      write(*,*) 'stc_gen: GETRS: info = ', info
      stop
   endif
!
!..Compute condensed RHS Bi (store it in Bi)
!  Bi = Bi - Aib*inv(Abb)*Bb
!  GEMM: Bi = -1.0 * Aib * (inv(Abb)*Bb) + 1.0 * Bi
#if C_MODE
   call ZGEMM('N','N',Ni,NR_RHS,Nb,-ZONE,Aib,Ni,Bb,Nb,ZONE,Bi,Ni)
#else
   call DGEMM('N','N',Ni,NR_RHS,Nb,-ZONE,Aib,Ni,Bb,Nb,ZONE,Bi,Ni)
#endif
!
!..Compute condensed system Aii (store it in Aii)
!  Aii = Aii - Aib*inv(Abb)*Abi
!  GEMM: Aii = -1.0 * Aib * (inv(Abb)*Abi) + 1.0 * Aii
#if C_MODE
   call ZGEMM('N','N',Ni,Ni,Nb,-ZONE,Aib,Ni,Abi,Nb,ZONE,Aii,Ni)
#else
   call DGEMM('N','N',Ni,Ni,Nb,-ZONE,Aib,Ni,Abi,Nb,ZONE,Aii,Ni)
#endif
!
end subroutine stc_fwd_gen
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_bwd_wrapper
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - loop through elements and perform back-
!                         substitution using stored Schur complements.
!
!    arguments:
!      in:
!              Iel      - Global element id (used in CLOC)
!                         Subdomain id (if mesh is distributed)
!
!----------------------------------------------------------------------
subroutine stc_bwd_wrapper(Iel)
!
   use assembly, ONLY:  MAXbrickH,MAXmdlbH,MAXEQNH, &
                        MAXbrickE,MAXmdlbE,MAXEQNE, &
                        MAXbrickV,MAXmdlbV,MAXEQNV
!
   implicit none
!
   integer, intent(in) :: Iel
!
   integer :: mdle
!
!..number of interface and bubble dof for each variable
   integer :: nrdofi(NR_PHYSA),nrdofb(NR_PHYSA)
!
!..workspace for local solution dofs
   VTYPE :: zdofH(MAXEQNH,MAXbrickH-MAXmdlbH)
   VTYPE :: zdofE(MAXEQNE,MAXbrickE-MAXmdlbE)
   VTYPE :: zdofV(MAXEQNV,MAXbrickV-MAXmdlbV)
!
!..workspace for interface and bubble solution
   integer            :: ni,nb
   VTYPE, allocatable :: xi(:,:),xb(:,:)
!
!..aux variables
   integer :: i,j,k,l,ii,load,lH,lE,lV
!
!--------------------------------------------------------------
!
   mdle = CLOC(iel)%mdle
!
!..get interface and bubble dof for each variable
   call stc_get_nrdof(mdle, nrdofi,nrdofb)
!
!..initialize interface and bubble dof counter
   ni = 0; nb = 0
   do j=1,NR_PHYSA
      ni = ni + nrdofi(j)
      nb = nb + nrdofb(j)
   enddo
!
!..nothing to do if there are no bubble dofs
   if (nb .eq. 0) goto 500
!
!..get local solution dofs
!..these are returned in the EXPANDED mode
!..this includes even components of variables that
!  are not supported in this element
   call solelmI(mdle, zdofH,zdofE,zdofV)
!
   allocate(xi(ni,NR_RHS))
!
!..collapse interface solution dofs into one vector
   lH = 0; lE = 0; lV = 0;
   do load=1,NR_RHS
      ni = 0
      do j=1,NR_PHYSA
         l = NR_COMP(j)
         if (.not. PHYSAm(j)) goto 100
         ii = nrdofi(j)
         if (ii .eq. 0) cycle
         k = ii/l
!     ...copy data
         select case(DTYPE(j))
            case('contin')
               xi(ni+1:ni+ii,load) = RESHAPE(zdofH(lH+1:lH+l,1:k),(/ii/))
            case('tangen')
               xi(ni+1:ni+ii,load) = RESHAPE(zdofE(lE+1:lE+l,1:k),(/ii/))
            case('normal')
               xi(ni+1:ni+ii,load) = RESHAPE(zdofV(lV+1:lV+l,1:k),(/ii/))
#if DEBUG_MODE
            case default
               write(*,*) 'stc_bwd_wrapper: INCONSISTENCY. stop.'
               stop
#endif
         end select
!
         ni = ni+ii
 100     continue
         select case(DTYPE(j))
            case('contin'); lH = lH + l
            case('tangen'); lE = lE + l
            case('normal'); lV = lV + l
         end select
      enddo
   enddo
!
!..perform backsubstitution
   allocate(xb(nb,NR_RHS))

   if (STORE_STC) then
      xb = CLOC(Iel)%BSchur
      deallocate(CLOC(Iel)%BSchur)
      call stc_bwd(ni,nb,CLOC(Iel)%ASchur,xi, xb)
      deallocate(CLOC(Iel)%ASchur)
   else
      write(*,*) 'stc_bwd_wrapper: recomputing schur complement not yet implemented. stop.'
      stop
   endif
!
!..add bubble solution to data structure
   call stc_solout(CLOC(Iel)%mdle,nb,xb,nrdofb)
!
   deallocate(xi,xb)
!
  500 continue
!
!   write(*,*) 'stc_bwd_wrapper finished. stop.'
!   stop
!
end subroutine stc_bwd_wrapper
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_bwd
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine perform backsubstitution to obtain
!                         bubble solution dof. It makes use of the 
!                         Schur complements stored in CLOC.
!
!    arguments :
!      in:
!              Ni           - total number of interface dof
!              Nb           - total number of bubble dof
!              Abi          - left-hand side Schur complement
!              Xi           - solution dof : interface
!      inout:
!              Xb           - in: right-hand side Schur complement
!                           - out: solution dof : bubble
!
!----------------------------------------------------------------------
subroutine stc_bwd(Ni,Nb,Abi,Xi, Xb)
!
   implicit none
!
   integer, intent(in)    :: Ni,Nb
   VTYPE  , intent(in)    :: Abi(Nb,Ni)
   VTYPE  , intent(in)    :: Xi(Ni,NR_RHS)
   VTYPE  , intent(inout) :: Xb(Nb,NR_RHS)
!
!----------------------------------------------------------------------
!
!..GEMM: Xb = -1.0 * Abi * Xi + 1.0 * Xb
#if C_MODE
   call ZGEMM('N','N',Nb,NR_RHS,Ni,-ZONE,Abi,Nb,Xi,Ni,ZONE,Xb,Nb)
#else
   call DGEMM('N','N',Nb,NR_RHS,Ni,-ZONE,Abi,Nb,Xi,Ni,ZONE,Xb,Nb)
#endif
!
end subroutine stc_bwd
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_solout
!
!----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - routine writes bubble solution dof to data
!                         structure
!
!    arguments :
!      in:
!              Mdle         - global middle node number
!              Nb           - total number of bubble dof
!              Xb           - bubble solution dof
!              Nrdofb       - number of bubble solution dof per
!                             physics variable
!
!----------------------------------------------------------------------
subroutine stc_solout(Mdle,Nb,Xb,Nrdofb)
!
   implicit none
!
   integer, intent(in) :: Mdle,Nb
   VTYPE  , intent(in) :: Xb(Nb,NR_RHS)
   integer, intent(in) :: Nrdofb(NR_PHYSA)
!
!..decoded index for a node
   integer :: index(NRINDEX), ncase(NR_PHYSA)
   integer :: nn,nod,nvarH,nvarE,nvarV,nvarQ
   integer :: ivar,i,j,k,l,ik,il,load
!
!..component counters for the nodes (use in case of multiple loads)
   integer :: mvarH,mvarE,mvarV,mvarQ
!
#if DEBUG_MODE
   integer :: iprint=0
#endif
!
!----------------------------------------------------------------------
!
!..initiate component counters
   mvarH = 0; mvarE = 0; mvarV = 0; mvarQ = 0
!
!..set node to be middle node number (global id in NODES)
   nod = Mdle
   call decod(NODES(nod)%case,2,NR_PHYSA, ncase)
!
!-----------------------------------------------------------------------
!
!..get index for each component of middle node
   call get_index(nod, index)
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7001) nod,index
 7001 format('stc_solout: nod = ',i5,' index = ',10i2)
   endif
#endif
!
!..compute the number of active non-flux H,E,V,Q variables for the node
   nvarH=0; nvarE=0; nvarV=0; nvarQ=0; k=0
   do i=1,NR_PHYSA
      il = NR_COMP(i)
      if ((.not. PHYSAm(i)) .or. PHYSAi(i)) then
         k=k+il
         cycle
      endif
      do l=1,il
         k=k+1
         select case(index(k))
            case(2); nvarH=nvarH+1
            case(4); nvarE=nvarE+1
            case(6); nvarV=nvarV+1
            case(8); nvarQ=nvarQ+1
         end select
      enddo
   enddo
!
!..initialize the number of H1,H(curl),H(div),L2 components stored so far
   ivar=0
!
!..loop through right-hand sides (loads).............................
   do load=1,NR_RHS
!
!  ...initiate the bubble dof counter
      nn=0
!
!  ...H1 dof .................................
      if (nvarH .eq. 0) goto 200
!  ...initiate index offset ik = 0
      ik = 0
!  ...loop over physics variables
      do i=1,NR_PHYSA
         il = NR_COMP(i)
         if (DTYPE(i) .ne. 'contin')  goto 110
         if (ncase(i) .eq. 0)         goto 110
         if ((.not. PHYSAm(i)) .or. PHYSAi(i)) then
            mvarH=mvarH+il; goto 110
         endif
!     ...loop over nodal dofs (per component)
         do j=1,nrdofb(i)/il
            ivar = mvarH
!        ...loop over components of this physics variable
            do l=1,il
!           ...get index of this component for this node
               k = ik+l
               select case(index(k))
                  case(2)
                     ivar = ivar+1
                     nn = nn+1
!                 ...copy the dof
                     NODES(nod)%dof%zdofH(ivar,j) = Xb(nn,load)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Xb(nn,load)
 7006                   format('stc_solout: nn,load,Xb(nn,load) = ',i4,i2,x,2e13.5)
                        write(*,7007) nod,j,ivar,NODES(nod)%dof%zdofH(ivar,j)
 7007                   format('stc_solout: nod,j,ivar,NODES(nod)%dof%zdofH(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
                  case default
!                 ...No other case should ever occur on middle node
                     write(*,*) 'stc_solout: INCONSISTENCY. stop.'
                     stop
#endif
               end select
            enddo
!     ...end loop over dofs
         enddo
!     ...update the number of components stored so far
         mvarH = ivar
  110    continue
         ik = ik+il
!  ...end loop over physics variables
      enddo
!
!-----------------------------------------------------------------------
!
!  ...H(curl) dof .................................
  200 continue
      if (nvarE .eq. 0) goto 300
!  ...initiate index offset ik = 0
      ik = 0
!  ...loop over physics variables
      do i=1,NR_PHYSA
         il = NR_COMP(i)
         if (DTYPE(i) .ne. 'tangen')  goto 210
         if (ncase(i) .eq. 0)         goto 210
         if ((.not. PHYSAm(i)) .or. PHYSAi(i)) then
            mvarE=mvarE+il; goto 210
         endif
!     ...loop over nodal dofs (per component)
         do j=1,nrdofb(i)/il
            ivar = mvarE
!        ...loop over components of this physics variable
            do l=1,il
!           ...get index of this component for this node
               k = ik+l
               select case(index(k))
                  case(4)
                     ivar = ivar+1
                     nn = nn+1
!                 ...copy the dof
                     NODES(nod)%dof%zdofE(ivar,j) = Xb(nn,load)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Xb(nn,load)
                        write(*,7009) nod,j,ivar,NODES(nod)%dof%zdofE(ivar,j)
 7009                   format('stc_solout: nod,j,ivar,NODES(nod)%dof%zdofE(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
                  case default
!                 ...No other case should ever occur on middle node
                     write(*,*) 'stc_solout: INCONSISTENCY. stop.'
                     stop
#endif
               end select
            enddo
!     ...end loop over dofs
         enddo
!     ...update the number of components stored so far
         mvarE = ivar
  210    continue
         ik = ik+il
!  ...end loop over physics variables
      enddo
!
!-----------------------------------------------------------------------
!
!  ...H(div) dof .................................
  300 continue
      if (nvarV .eq. 0) goto 400
!  ...initiate index offset ik = 0
      ik = 0
!  ...loop over physics variables
      do i=1,NR_PHYSA
         il = NR_COMP(i)
         if (DTYPE(i) .ne. 'normal')  goto 310
         if (ncase(i) .eq. 0)         goto 310
         if ((.not. PHYSAm(i)) .or. PHYSAi(i)) then
            mvarV=mvarV+il; goto 310
         endif
!     ...loop over nodal dofs (per component)
         do j=1,nrdofb(i)/il
            ivar = mvarV
!        ...loop over components of this physics variable
            do l=1,il
!           ...get index of this component for this node
               k = ik+l
               select case(index(k))
                  case(6)
                     ivar = ivar+1
                     nn = nn+1
!                 ...copy the dof
                     NODES(nod)%dof%zdofV(ivar,j) = Xb(nn,load)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Xb(nn,load)
                        write(*,7010) nod,j,ivar,NODES(nod)%dof%zdofV(ivar,j)
 7010                   format('stc_solout: nod,j,ivar,NODES(nod)%dof%zdofV(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
                  case default
!                 ...No other case should ever occur on middle node
                     write(*,*) 'stc_solout: INCONSISTENCY. stop.'
                     stop
#endif
               end select
            enddo
!     ...end loop over dofs
         enddo
!     ...update the number of components stored so far
         mvarV = ivar
  310    continue
         ik = ik+il
!  ...end loop over physics variables
      enddo
!
!-----------------------------------------------------------------------
!
!  ...L2 dof .................................
  400 continue
      if (nvarQ .eq. 0) goto 500
!
!  ...initiate index offset ik = 0
      ik = 0
!  ...loop over physics variables
      do i=1,NR_PHYSA
         il = NR_COMP(i)
         if (DTYPE(i) .ne. 'discon')  goto 410
         if (ncase(i) .eq. 0)         goto 410
         if (.not. PHYSAm(i)) then
            mvarQ=mvarQ+il; goto 410
         endif
!     ...loop over nodal dofs (per component)
         do j=1,nrdofb(i)/il
!        ...ivar is component counter for this physics type ('discon')
!           it is used for addressing zdofQ in the data structure
!           ivar counts differently depending on which physics variables are supported
!        ...mvarQ is needed if using multiple loads
            ivar = mvarQ
!        ...loop over components of this physics variable
            do l=1,il
!           ...get index of this component for this node
               k = ik+l
               select case(index(k))
                  case(8)
                     ivar = ivar+1
                     nn = nn+1
!                 ...copy the dof
                     NODES(nod)%dof%zdofQ(ivar,j) = Xb(nn,load)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Xb(nn,load)
                        write(*,7011) nod,j,ivar,NODES(nod)%dof%zdofQ(ivar,j)
 7011                   format('stc_solout: nod,j,ivar,NODES(nod)%dof%zdofQ(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
                  case default
!                 ...No other case should ever occur on middle node
                     write(*,*) 'stc_solout: INCONSISTENCY. stop.'
                     stop
#endif
               end select
            enddo
!     ...end loop over dofs
         enddo
!     ...update the number of components stored so far
         mvarQ = ivar
  410    continue
         ik = ik+il
!  ...end loop over physics variables
      enddo
!
 500  continue
!
#if DEBUG_MODE
      if (iprint.eq.1) call pause
#endif
!
!..end of loop through loads
   enddo       
!
#if DEBUG_MODE
   if (iprint.eq.1) write(*,*) 'stc_solout finished.'
#endif
!
end subroutine
!
!
end module stc
