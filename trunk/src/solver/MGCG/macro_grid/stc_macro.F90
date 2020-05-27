!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!    routine name       - stc_fwd_wrapper_mg
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
!              Ielc     - Global coarse grid element id
!              Iel      - Local fine grid element id
!              Mdle     - Global middle node id (used in NODES)
!
!----------------------------------------------------------------------
subroutine stc_fwd_wrapper_mg(Ielc,Iel,Mdle)
!
   use assembly,           only: ALOC, BLOC, NR_RHS, NR_PHYSA
   use macro_grid_info,    only: A_MACRO
   use mg_data_structure,  only: ISTORE, ISTORE_YES
   use stc,                only: stc_get_nrdof,stc_fwd_herm
!
   implicit none
!
   integer, intent(in) :: Ielc,Iel,Mdle
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
   A_MACRO(Ielc)%GLOC(Iel)%mdle = Mdle
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
   call stc_fwd_herm(ni,nb, Aii,Abi,Aib,Abb,Bi,Bb)
!
!..Store Schur complements if STORE_STC
!..   inv(Abb)*Bb
!..   inv(Abb)*Abi
   if (ISTORE .eq. ISTORE_YES) then
      allocate(A_MACRO(Ielc)%GLOC(Iel)%array(nb,ni))
      allocate(A_MACRO(Ielc)%GLOC(Iel)%vect(nb))
!
      A_MACRO(Ielc)%GLOC(Iel)%array = Abi
      A_MACRO(Ielc)%GLOC(Iel)%vect  = Bb(:,1)
   else
      write(*,*) 'stc_fwd_wrapper_mg: ', &
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
end subroutine stc_fwd_wrapper_mg
!
!
!----------------------------------------------------------------------
!
!    routine name       - stc_bwd_wrapper_mg
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
!              Ielc     - Global coarse grid element id
!              Iel      - Local fine grid element id
!
!----------------------------------------------------------------------
subroutine stc_bwd_wrapper_mg(Ielc,Iel)
!
   use assembly, ONLY:  MAXbrickH,MAXmdlbH,MAXEQNH, &
                        MAXbrickE,MAXmdlbE,MAXEQNE, &
                        MAXbrickV,MAXmdlbV,MAXEQNV, &
                        NR_RHS, NR_PHYSA, DTYPE, PHYSAm, NR_COMP
   use macro_grid_info,    only: A_MACRO
   use mg_data_structure,  only: ISTORE, ISTORE_YES
   use stc,                only: stc_get_nrdof,stc_bwd,stc_solout
!
   implicit none
!
   integer, intent(in) :: Ielc,Iel
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
   mdle = A_MACRO(Ielc)%GLOC(Iel)%mdle
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
               write(*,*) 'stc_bwd_wrapper_mg: INCONSISTENCY. stop.'
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

   if (ISTORE .eq. ISTORE_YES) then
      xb(:,1) = A_MACRO(Ielc)%GLOC(Iel)%vect
      deallocate(A_MACRO(Ielc)%GLOC(Iel)%vect)
      call stc_bwd(ni,nb,A_MACRO(Ielc)%GLOC(Iel)%array,xi, xb)
      deallocate(A_MACRO(Ielc)%GLOC(Iel)%array)
   else
      write(*,*) 'stc_bwd_wrapper_mg: recomputing schur complement not yet implemented. stop.'
      stop
   endif
!
!..add bubble solution to data structure
   call stc_solout(A_MACRO(Ielc)%GLOC(Iel)%mdle,nb,xb,nrdofb)
!
   deallocate(xi,xb)
!
  500 continue
!
end subroutine stc_bwd_wrapper_mg
!
