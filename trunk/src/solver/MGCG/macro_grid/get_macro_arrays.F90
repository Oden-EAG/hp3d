!-----------------------------------------------------------------------
!
!    routine name       - get_macro_arrays
!
!-----------------------------------------------------------------------
!
!    latest revision    - July 2018
!
!    purpose            - returns Schur complement matrices for
!                         computation of statically condensed unknowns
!
!
!   arguments :
!     in:
!              Ielc     - element number
!     out:
!
!             Abi       - Schur Complement matrix
!             bb        - Schur Complement vector
!
!-----------------------------------------------------------------------
!
#include "implicit_none.h"
!
   subroutine get_macro_arrays(Igrid, Ielc, Abi, bb,nb,ni)
!
   use data_structure3D, ONLY: NODES,MAXNODM, NRNODS,NR_PHYSA
   use physics,          ONLY: NRHVAR, NREVAR, NRVVAR, NRQVAR
   use assembly,         ONLY: NR_RHS, NEXTRACT, IDBC, ZDOFD,   &
                               BLOC, AAUX, ALOC, ZBMOD, ZAMOD
   use parameters,       ONLY: ZERO, ZONE,      &
                               MAXbrickH, MAXbrickE, MAXbrickV, MAXbrickQ
   use macro_grid_info,  ONLY: A_MACRO
   use mg_data_structure,ONLY: GRID, NODES_MG, ISTORE, ISTORE_YES
!
   implicit none
!
!-----------------------------------------------------------------------
!
!..globals
   integer, intent( in) :: Igrid, Ielc, nb, ni
   VTYPE,   intent(out) :: Abi(nb,ni), bb(nb)
!
!..dof offsets
   integer, allocatable :: mdle_macro(:)
!
!..workspace for celem
   integer, allocatable :: maxdofs(:)
   integer :: maxdofm, nrdofs(NR_PHYSA), nrdofm, nrdofc, nrnodm
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM),  &
              ndofmE(MAXNODM),ndofmV(MAXNODM), ndofmQ(MAXNODM)
!
!..dof counters (i for interface, b for bubble)
   integer :: nrdof_Hi, nrdof_Ei, nrdof_Vi, nrdof_Qi
   integer :: nrdof_Hb, nrdof_Eb, nrdof_Vb, nrdof_Qb
   integer :: nrdofi, nrdofb
   integer :: nHb, nEb, nVb, nQ, nHc, nEc, nVc

!
   integer, allocatable :: naH(:),naE(:), naV(:)
!..various
   integer :: mdle, nrmdle, iel, nod, master, ndof, ndofi, ndofb, idec
   integer :: l, lb, j, i, k1, k2, k, ii, nnz_tot, nz_b, nz_ib,nsum_b, nsum_i
   integer, allocatable :: lcon(:)
   character(len=4) :: type

   VTYPE :: zvoid(1)
!
!..required structure for efficient assembly of global stiffness
!..matrix to CSR format
   type sparse_matrix
!  ...number of non-zero with duplicates
      integer :: nz_old
!  ...number of non-zero after removing the duplicates
      integer :: nz
!  ...row indices for each column
      integer, allocatable :: col(:)
!  ...values of the matrix corresponding to ia ordering
      VTYPE,  allocatable  :: val(:)
!
   end type sparse_matrix
!
   type(sparse_matrix), allocatable :: SKbb(:)
!
!..global matrices

   VTYPE, allocatable :: Aii(:,:), Atemp(:,:), btemp(:)
   VTYPE, allocatable :: bi(:)
   VTYPE, allocatable :: SAbb(:)
   VTYPE, allocatable :: Abi_ext(:,:)
   VTYPE, allocatable :: zload(:), ztemp(:), zloadc(:), ztempc(:)
!
   integer, allocatable :: IAbb(:), JAbb(:), IAib(:), JAib(:)
!
   character(len=1)     :: transa
   character(len=1)     :: matdescra(4)
   integer              :: n1, n2
!
!-----------------------------------------------------------------------
!
   maxdofm = MAXbrickH*NRHVAR + MAXbrickE*NREVAR  &
           + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
!
!..allocate required variables for celem
   allocate(NEXTRACT(maxdofm), IDBC(maxdofm))
   allocate(ZDOFD(maxdofm,NR_RHS))
!
   allocate(maxdofs(NR_PHYSA)); maxdofs = 0
!
!..initialize offsets
   allocate(naH(NRNODS)); naH = -1
   allocate(naE(NRNODS)); naE = -1
   allocate(naV(NRNODS)); naV = -1
!
!..Initialize dof counters
   nrdof_Hi = 0 ; nrdof_Ei = 0; nrdof_Vi = 0; nrdof_Qi = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0; nrdof_Qb = 0
!
   nz_b = 0; nz_ib = 0
!
   nrmdle = A_MACRO(Ielc)%nrmdle
!
   allocate(mdle_macro(nrmdle))
!
   idec = 1
!..loop through the sub-element within the macro elements
   do iel = 1,nrmdle
!
      mdle =  A_MACRO(ielc)%GLOC(iel)%mdle
!
      call celem1(mdle,idec,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                  ndofmE,ndofmV,ndofmQ,nrnodm, zvoid,zvoid)
!
!  ...update the maximum number of local dof
      do i=1,NR_PHYSA
         maxdofs(i) = max0(maxdofs(i),nrdofs(i))
      enddo
!
!  ...update the maximum number of modified element dof in the expanded mode
      maxdofm = max0(maxdofm,nrdofm)
!
!-------------------------------------------------------------------------------
!
!  ...count non zero for sparse assembly of Abb
      nsum_b = 0; nsum_i = 0
      do i =1, nrnodm-1
         nod = nodm(i)
         call find_master(Igrid,nod, master)
         type = NODES(master)%type
!
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
            nsum_b = nsum_b + ndofmH(i) + ndofmE(i) + ndofmV(i)
         case default
            nsum_i = nsum_i + ndofmH(i) + ndofmE(i) + ndofmV(i)
         end select
!
      enddo
!  ...non zero for sparse assembly
      nz_b  = nz_b  + nsum_b*(nsum_b+1)/2
      nz_ib = nz_ib + nsum_i*nsum_b
!
!-------------------------------------------------------------------------------
!
!  ...compute offsets for H1 (excluding interior dof)
      do i = 1, nrnodm-1
         nod = nodm(i)
         if (naH(nod).ge.0) cycle
!
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!
!     ...check if the the type of the master node
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!
!        ...store the first dof offset
            naH(nod) = nrdof_Hb
!        ...update the H1 dof counter
            nrdof_Hb = nrdof_Hb + ndofmH(i)
!
         case default
!        ...store the first dof offset
            naH(nod) = nrdof_Hi
!        ...update the H1 dof counter
            nrdof_Hi = nrdof_Hi + ndofmH(i)
         end select
      enddo
!
!  ...compute offsets for H(curl) (excluding interior dof)
      do i = 1, nrnodm-1
         nod = nodm(i)
         if (naE(nod).ge.0) cycle
!
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!
!     ...check if the the type of the master node
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!
!        ...store the first dof offset
            naE(nod) = nrdof_Eb
!        ...update the H(curl) dof counter
            nrdof_Eb = nrdof_Eb + ndofmE(i)
!
         case default
!
!        ...store the first dof offset
            naE(nod) = nrdof_Ei
!        ...update the H(curl) dof counter
            nrdof_Ei = nrdof_Ei + ndofmE(i)
         end select
      enddo
!
!  ...compute offsets for H(div) (excluding interior dof)
      do i = 1, nrnodm-1
         nod = nodm(i)
         if (naV(nod).ge.0) cycle
!
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!
!     ...check if the the type of the master node
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!
!        ...store the first dof offset
            naV(nod) = nrdof_Vb
!
!        ...update the H(div) dof counter
            nrdof_Vb = nrdof_Vb + ndofmV(i)
!
         case default
!
!        ...store the first dof offset
            naV(nod) = nrdof_Vi
!        ...update the H(div) dof counter
            nrdof_Vi = nrdof_Vi + ndofmV(i)
         end select
      enddo
!
!..end of loop through elements
   enddo
!
   deallocate(NEXTRACT,IDBC,ZDOFD)
!
!..macro dof counters
   nrdofi = nrdof_Hi + nrdof_Ei + nrdof_Vi +  nrdof_Qi
   nrdofb = nrdof_Hb + nrdof_Eb + nrdof_Vb +  nrdof_Qb
!
!..build and store connectivity maps
   if (nrdofb .gt. 0) then
      allocate(SKbb(nrdofb))
!  ...initialize
      do j=1,nrdofb
         SKbb(j)%nz_old = 0
      enddo
   endif
!
!----------------------------------------------------------------------
!                              STEP 2
!                  BUILD AND STORE CONNECTIVITY MAPS
!----------------------------------------------------------------------
!
!..allocate required variables for celem
   allocate(NEXTRACT(maxdofm))
   allocate(IDBC(maxdofm))
   allocate(ZDOFD(maxdofm,NR_RHS))
   allocate(lcon(maxdofm))

!
!..loop through elements
   do iel = 1, nrmdle
!
      nHc = A_MACRO(Ielc)%GLOC(iel)%nHc
      nEc = A_MACRO(Ielc)%GLOC(iel)%nEc
      nVc = A_MACRO(Ielc)%GLOC(iel)%nVc
      nHb = A_MACRO(Ielc)%GLOC(iel)%nHb
      nEb = A_MACRO(Ielc)%GLOC(iel)%nEb
      nVb = A_MACRO(Ielc)%GLOC(iel)%nVb
      nQ  = A_MACRO(Ielc)%GLOC(iel)%nQ
!
      ndof = nHc + nEc + nVc
      ndofb = nHb + nEb + nVb + nQ
      ndofi = l-lb
!
      lcon(1:ndof) = A_MACRO(ielc)%GLOC(iel)%lcon(1:ndof)
!
!  ...loop through number of dof for the element
      if (ndofb .eq. 0) cycle
      do k1=1,ndof
!
!     ...global dof is:
         i = lcon(k1)
         if (i .lt. 0 ) then
            i = abs(i)
!
!        ...loop through dof
            do k2=1,ndof
!
!           ...global dof is:
               j = lcon(k2)
               if (j .lt. 0 ) then
                  j = abs(j)
                  if (i.le.j) then
!
!                 ...count non-zero values with duplicates
                     SKbb(i)%nz_old = SKbb(i)%nz_old + 1
                  endif
               endif
            enddo
         endif

!  ...end of loop through number of dof for the element
      enddo
!
!..end of loop through elements
   enddo
!
   if (nrdofb .gt. 0 ) then
      do i = 1,nrdofb
         allocate(SKbb(i)%col(SKbb(i)%nz_old))
         allocate(SKbb(i)%val(SKbb(i)%nz_old))
!     ...reinitialize counters
         SKbb(i)%nz_old = 0
      enddo
!
      allocate(Abi_ext(nrdofb,1+nrdofi)); Abi_ext = ZERO
   endif
!
!..allocate required arrays for celem
   allocate(BLOC(NR_PHYSA))
   allocate(AAUX(NR_PHYSA))
   allocate(ALOC(NR_PHYSA,NR_PHYSA))
   do i=1,NR_PHYSA
      BLOC(i)%nrow = maxdofs(i)
      BLOC(i)%ncol = NR_RHS
      allocate(BLOC(i)%array(maxdofs(i),NR_RHS))
      do j=1,NR_PHYSA
         ALOC(i,j)%nrow = maxdofs(i)
         ALOC(i,j)%ncol = maxdofs(j)
         allocate(ALOC(i,j)%array(maxdofs(i),maxdofs(j)))
      enddo
      AAUX(i)%nrow = maxdofm
      AAUX(i)%ncol = maxdofs(i)
      allocate(AAUX(i)%array(maxdofm,maxdofs(i)))
   enddo
   allocate(ZBMOD(maxdofm,NR_RHS))
   allocate(ZAMOD(maxdofm,maxdofm))
   allocate(zload(maxdofm),ztemp(maxdofm**2))
!
!..loop through elements
   idec = 2
!..reset counter
   nz_ib = 0

   do iel=1,nrmdle
!  ...pick up the mdle number
      mdle =  A_MACRO(ielc)%GLOC(iel)%mdle
!
      call celem1(mdle,idec,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                  ndofmE,ndofmV,ndofmQ,nrnodm,zload,ztemp)
!

      nHc = A_MACRO(Ielc)%GLOC(iel)%nHc
      nEc = A_MACRO(Ielc)%GLOC(iel)%nEc
      nVc = A_MACRO(Ielc)%GLOC(iel)%nVc
      nHb = A_MACRO(Ielc)%GLOC(iel)%nHb
      nEb = A_MACRO(Ielc)%GLOC(iel)%nEb
      nVb = A_MACRO(Ielc)%GLOC(iel)%nVb
      nQ  = A_MACRO(Ielc)%GLOC(iel)%nQ
!
      n1 = nHc + nEc + nVc
      n2 = nHb + nEb + nVb + nQ
      ndof = n1 + n2
!
      allocate(ztempc(n1**2))
      allocate(zloadc(n1))
!
      if (n2 .gt. 0) then
         call stc_macro_out_nostore(Ielc,iel,ztemp(1:ndof**2),zload(1:ndof),ndof,n1,n2,  &
                            nHc,nHb,nEc,nEb,nVc,nVb,nQ,ztempc,zloadc)
      else
         ztempc(1:n1**2) = ztemp(1:n1**2)
         zloadc(1:n1) = zload(1:n1)
      endif
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      ndof = n1
      do k1=1,ndof
!     ...global dof is:
         i = A_MACRO(ielc)%GLOC(iel)%lcon(k1)
         if (i .lt. 0 ) then
            i = abs(i)
!        ...Assemble global load vector
            Abi_ext(i,nrdofi+1) = Abi_ext(i,nrdofi+1) + zloadc(k1)
!        ...loop through dof
            do k2=1,ndof
!           ...global dof is:
               j = A_MACRO(ielc)%GLOC(iel)%lcon(k2)
               if (j .lt. 0 ) then
                  j = abs(j)
                  if (i.le.j) then
                     k = (k1-1)*ndof + k2
!                    ...sparse storage
                        SKbb(i)%nz_old = SKbb(i)%nz_old + 1
                        SKbb(i)%col(SKbb(i)%nz_old) = j
                        SKbb(i)%val(SKbb(i)%nz_old) = ztempc(k)
                  endif
               else
                  k = (k1-1)*ndof + k2
                  Abi_ext(i,j) = Abi_ext(i,j) + ztempc(k)
               endif
            enddo
         endif
      enddo
      deallocate(ztempc)
      deallocate(zloadc)
! ...end of loop through the elements
   enddo
!
!..deallocate everything that is not required further
   do i=1,NR_PHYSA
      deallocate(BLOC(i)%array)
      do j=1,NR_PHYSA
         deallocate(ALOC(i,j)%array)
      enddo
      deallocate(AAUX(i)%array)
   enddo

   deallocate(zload,ztemp)
   deallocate(lcon,NEXTRACT,IDBC,ZDOFD,BLOC,AAUX,ALOC,ZBMOD,ZAMOD)
   deallocate(mdle_macro)
   deallocate(naH, naE, naV)
   deallocate(maxdofs)
!
   if (nrdofb .eq. 0) then
!  ...nothing to condense out. The element is not refined
      return
   endif

!..convert to compressed row format
   nnz_tot = 0
   do i = 1, nrdofb
      j = SKbb(i)%nz_old
      call assemble_double(SKbb(i)%col(1:j),SKbb(i)%val(1:j),j,SKbb(i)%nz)
      nnz_tot = nnz_tot + SKbb(i)%nz
   enddo
!
   allocate(IAbb(nnz_tot))
   allocate(JAbb(nnz_tot))
   allocate(SAbb(nnz_tot))
!
   k = 0
   do i = 1, nrdofb
      do j =  1, SKbb(i)%nz
         k = k + 1
         IAbb(k) = i
         JAbb(k) = SKbb(i)%col(j)
         SAbb(k) = SKbb(i)%val(j)
      enddo
   enddo
!
!..free memory
   do i = 1,nrdofb
      deallocate(SKbb(i)%col)
      deallocate(SKbb(i)%val)
   enddo
   deallocate(SKbb)
!
   call get_pointers(IAbb(1:nnz_tot), nnz_tot)

!
!..call pardiso to solve
   call pardiso_solve(IAbb(1:nrdofb+1),JAbb(1:nnz_tot),SAbb(1:nnz_tot),'H', &
                      nnz_tot,nrdofb,nrdofi+1,Abi_ext)

   Abi(1:nrdofb,1:nrdofi) = Abi_ext(1:nrdofb,1:nrdofi)
   bb(1:nrdofb) = Abi_ext(1:nrdofb,nrdofi+1)


   deallocate(Abi_ext)
   deallocate(SAbb,IAbb,JAbb)
!
!
   end subroutine get_macro_arrays
