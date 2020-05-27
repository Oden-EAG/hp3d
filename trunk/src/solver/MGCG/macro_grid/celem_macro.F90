!-----------------------------------------------------------------------
!
!    routine name       - celem_macro
!
!-----------------------------------------------------------------------
!
!    latest revision    - Feb 2018
!
!    purpose            - returns the macro-element stiffness matrix
!                         and load vector
!
!   arguments :
!     in:
!              Ielc     - element number
!             MdleC     - middle node of a coarse element
!     out:

!        Nrdof_macro    - number of condensed element dof
!
!        Nrnod_macro    - number of nodes
!        Nod_macro      - list of macro-element nodes nodes
!        NdofH_macro    - the corresponding number of H1 dof
!        NdofE_macro    - the corresponding number of H(curl) dof
!        NdofV_macro    - the corresponding number of H(div) dof
!
!             Zbload   - 1D array containing the condensed load vector
!             Zastif   - 1D array containing the condensed stiffness
!                        matrix
!
!-----------------------------------------------------------------------
!
#include "typedefs.h"
!
   subroutine celem_macro(Igrid,Idec_sch,Ielc,MdleC, Nrdof_macro,Nrnod_macro,Nod_macro,  &
                          NdofH_macro,NdofE_macro,NdofV_macro,Zbload,Zastif)
!
   use data_structure3D, ONLY: NODES,MAXNODM, NRNODS,NR_PHYSA
   use physics,          ONLY: NRHVAR, NREVAR, NRVVAR, NRQVAR
   use assembly,         ONLY: NR_RHS, NEXTRACT, IDBC, ZDOFD,   &
                               BLOC, AAUX, ALOC, ZBMOD, ZAMOD,  &
                               ZERO, ZONE,                      &
                               MAXbrickH, MAXmdlbH, NRHVAR,     &
                               MAXbrickE, MAXmdlbE, NREVAR,     &
                               MAXbrickV, MAXmdlbV, NRVVAR
   use macro_grid_info,  ONLY: A_MACRO
   use mg_data_structure,ONLY: GRID, NODES_MG, ISTORE, ISTORE_YES
!
   implicit none
!
!-----------------------------------------------------------------------
!
!..globals
   integer, intent( in) :: Igrid,Idec_sch,Ielc,MdleC
   integer, intent(out) :: Nrdof_macro, Nrnod_macro
   integer, intent(out) :: Nod_macro(*),   NdofH_macro(*), &
                           NdofE_macro(*), NdofV_macro(*)
   VTYPE,   intent(out) :: Zbload(*), Zastif(*)
!
!..locals
!..dof offsets
   integer, allocatable :: naH(:),naE(:), naV(:), nvisit(:)
   integer, allocatable :: mdle_macro(:)
!
!..workspace for celem
   integer, allocatable :: maxdofs(:)
   integer :: maxdofm, nrdofs(NR_PHYSA), nrdofm, nrdofc, nrnodm
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM),  &
              ndofmE(MAXNODM),ndofmV(MAXNODM), ndofmQ(MAXNODM)
!
!
!..dof counters (i for interface, b for bubble)
   integer :: nrdof_Hi, nrdof_Ei, nrdof_Vi
   integer :: nrdof_Hb, nrdof_Eb, nrdof_Vb
   integer :: nrdofi, nrdofb
!
!..various
   integer :: mdle, nrmdle, iel, nod, master, ndof, ndofi, ndofb, idec
   integer :: i, ii, j, k1, k2, k, l, lb
   integer, allocatable :: lcon(:)
   character(len=4) :: type

   VTYPE :: zvoid(1)
!
!..for sparse assembly
   integer  :: nsum_b, nsum_i, nz_b, nz_ib, nnz_tot
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

   VTYPE, allocatable :: Abi(:,:), Aii(:,:), Atemp(:,:), btemp(:)
   VTYPE, allocatable :: bb(:), bi(:)
   VTYPE, allocatable :: SAbb(:),SAib(:)
   VTYPE, allocatable :: DA(:),Abi_ext(:,:)
   VTYPE, allocatable :: zload(:), ztemp(:), zloadc(:), ztempc(:)
!
   integer, allocatable :: IAbb(:), JAbb(:), IAib(:), JAib(:)
!
   character(len=1)     :: transa
   character(len=1)     :: matdescra(4)
   integer              :: nk, n1, n2
   nk(n1,n2) = (n2-1)*n2/2+n1
!
!-----------------------------------------------------------------------
!
   maxdofm  = (MAXbrickH-MAXmdlbH)*NRHVAR   &
            + (MAXbrickE-MAXmdlbE)*NREVAR   &
            + (MAXbrickV-MAXmdlbV)*NRVVAR
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
   allocate(nvisit(NRNODS)); nvisit = 0
!
!..Initialize dof counters
   nrdof_Hi = 0 ; nrdof_Ei = 0; nrdof_Vi = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0
!
   mdle=0 ;   Nrnod_macro = 0 ; nrmdle = 0
!
   nz_b = 0; nz_ib = 0
   ! m_nz = 0
!
!..count number of sub-elements within the macro-element
   call nelcon_macro(MdleC,mdle, mdle)
   do while (mdle.ne.0)
      nrmdle = nrmdle +1
      call nelcon_macro(MdleC,mdle, mdle)
   enddo
!
   allocate(mdle_macro(nrmdle))
!
   mdle = 0; idec = 1
!..loop through the sub-element within the macro elements
   do iel = 1,nrmdle
!
      call nelcon_macro(MdleC,mdle, mdle)
!
      mdle_macro(iel) = mdle
!
      call celem_mg(Ielc,iel,mdle,idec,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
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
      do i =1, nrnodm
         nod = nodm(i)
         call find_master(Igrid,nod, master)
         type = NODES(master)%type
!
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
            nsum_b = nsum_b + ndofmH(i) + ndofmE(i) + ndofmV(i)
         case default
            nsum_i = nsum_i + ndofmH(i) + ndofmE(i) + ndofmV(i)
            if (nvisit(nod).eq.0) then
               Nrnod_macro = Nrnod_macro + 1
               Nod_macro(Nrnod_macro) = nod
               ndofH_macro(Nrnod_macro) = ndofmH(i)
               ndofE_macro(Nrnod_macro) = ndofmE(i)
               ndofV_macro(Nrnod_macro) = ndofmV(i)
               nvisit(nod) = 1
            endif
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
      do i = 1, nrnodm
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
!
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
      do i = 1, nrnodm
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
!
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
      do i = 1, nrnodm
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
   nrdofi = nrdof_Hi + nrdof_Ei + nrdof_Vi
   nrdofb = nrdof_Hb + nrdof_Eb + nrdof_Vb
   Nrdof_macro = nrdofi
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
   allocate(A_MACRO(ielc)%GLOC(nrmdle))
   A_MACRO(Ielc)%nrmdle = nrmdle
!
!
   if (Idec_sch .eq. 1) then
      allocate(GRID(igrid)%sch(ielc)%gloc(nrmdle))
      GRID(igrid)%sch(ielc)%nrmdle = nrmdle
   endif
!
!..allocate required variables for celem
   allocate(NEXTRACT(maxdofm))
   allocate(IDBC(maxdofm))
   allocate(ZDOFD(maxdofm,NR_RHS))
   allocate(lcon(maxdofm)); lcon = 0
!
!..loop through elements
   do iel = 1, nrmdle
!  ...save the mdle node number
      A_MACRO(ielc)%GLOC(iel)%mdle = mdle_macro(iel)
      if (Idec_sch .eq. 1) then
         GRID(igrid)%sch(ielc)%gloc(iel)%mdle = mdle_macro(iel)
         GRID(igrid)%sch(ielc)%gloc(iel)%iel  = NODES_MG(mdle_macro(iel))%visit
      endif
      call celem_mg(Ielc,iel,mdle_macro(iel),idec,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                    ndofmE,ndofmV,ndofmQ,nrnodm, zvoid,zvoid)
!  ...construct local to global connectivities
      l = 0; lb = 0
! ....H1 dof (excluding interior dof)
      do i = 1,nrnodm
         nod = nodm(i)
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check the type of the master node
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
            do j=1,ndofmH(i)
               l=l+1; lb = lb+1
               lcon(l) = -(naH(nod)+j)
            enddo
         case default
            do j=1,ndofmH(i)
               l=l+1
               lcon(l) = naH(nod)+j
            enddo
         end select
      enddo
!
!  ...H(curl) dof (excluding interior dof)
      do i = 1,nrnodm
         nod = nodm(i)
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check if the the type of the master node
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
            do j=1,ndofmE(i)
               l=l+1; lb=lb+1
               lcon(l) = -(nrdof_Hb + naE(nod)+j)
            enddo
         case default
            do j=1,ndofmE(i)
               l=l+1
               lcon(l) = nrdof_Hi + naE(nod)+j
            enddo
         end select
      enddo
!
!  ...H(div) dof (excluding interior dof)
      do i = 1,nrnodm
         nod = nodm(i)
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check if the the type of the master node
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
            do j=1,ndofmV(i)
               l=l+1; lb=lb+1
               lcon(l) = -(nrdof_Hb + nrdof_Eb + naV(nod)+j)
            enddo
         case default
            do j=1,ndofmV(i)
               l=l+1 ; lb=lb+1
               lcon(l) = nrdof_Hi + nrdof_Ei + naV(nod)+j
            enddo
         end select
      enddo
!
!  ...local number of dof is ndof_tot
      ndof = l
      ndofb = lb
      ndofi = l-lb

      A_MACRO(Ielc)%GLOC(iel)%ni = ndof

      allocate(A_MACRO(ielc)%GLOC(iel)%lcon(ndof))
      A_MACRO(ielc)%GLOC(iel)%lcon(1:ndof) = lcon(1:ndof)

      if (Idec_sch .eq. 1) then
         allocate(GRID(igrid)%sch(ielc)%gloc(iel)%lcon(ndof))
         GRID(igrid)%sch(ielc)%gloc(iel)%lcon = lcon
         GRID(igrid)%sch(ielc)%gloc(iel)%ndof_c = ndof
      endif

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
      allocate(Abi(nrdofb,nrdofi)); Abi = ZERO;
      allocate(bb(nrdofb)); bb = ZERO
      allocate(SAib(nz_ib),IAib(nz_ib),JAib(nz_ib));
      allocate(Abi_ext(nrdofb,1+nrdofi)); Abi_ext = ZERO
   endif
!
   allocate(Aii(nrdofi,nrdofi)); Aii = ZERO
   allocate(bi(nrdofi)); bi = ZERO;
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
!  ...save the mdle number
      call celem_mg(Ielc,iel,mdle_macro(iel),idec, nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                    ndofmE,ndofmV,ndofmQ,nrnodm,zload,ztemp)
!
      ndof = nrdofc
!
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      do k1=1,ndof
!     ...global dof is:
         i = A_MACRO(ielc)%GLOC(iel)%lcon(k1)
         if (i .lt. 0 ) then
            i = abs(i)
!        ...Assemble global load vector
            Abi_ext(i,nrdofi+1) = Abi_ext(i,nrdofi+1) + zload(k1)
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
                        SKbb(i)%val(SKbb(i)%nz_old) = ztemp(k)
                  endif
               else
                  k = (k1-1)*ndof + k2
                  Abi_ext(i,j) = Abi_ext(i,j) + ztemp(k)
               endif
            enddo
         else
!        ...Assemble global load vector
            bi(i) = bi(i) + zload(k1)
!        ...loop through dof
            do k2=1,ndof
               j = A_MACRO(ielc)%GLOC(iel)%lcon(k2)
               if (j .lt. 0 ) then
                  j = abs(j)
                  k = (k1-1)*ndof + k2
                  nz_ib = nz_ib + 1
                  SAib(nz_ib) = ztemp(k)
                  IAib(nz_ib) = i
                  JAib(nz_ib) = j
               else
                  k = (k1-1)*ndof + k2
                  Aii(i,j) = Aii(i,j) + ztemp(k)
               endif
            enddo
         endif
      enddo
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
   deallocate(naH, naE, naV, nvisit)
   deallocate(maxdofs)
!
   A_MACRO(ielc)%ni = nrdofi
   A_MACRO(ielc)%nb = nrdofb
   if (Idec_sch .eq. 1) then
      GRID(igrid)%sch(ielc)%n1 =  nrdofi
      GRID(igrid)%sch(ielc)%n2 =  nrdofb
   endif

   if (nrdofb .eq. 0) then
!
!  ...nothing to condense out. The element is not refined
      zbload(1:nrdofi) = bi
!
      do j = 1,nrdofi
         do i = 1,j
            k = nk(i,j)
            zastif(k) = Aii(i,j)
         enddo
      enddo
      deallocate(Aii, bi)
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
   call get_pointers(IAbb(1:nnz_tot),nnz_tot)


   if (Idec_sch .eq. 1) then
!  ...store the matrix A22 in sparse form
      GRID(Igrid)%sch(ielc)%inz = nnz_tot
      GRID(Igrid)%sch(ielc)%n1 =  nrdofi
      GRID(Igrid)%sch(ielc)%n2 =  nrdofb
!
      allocate(GRID(Igrid)%sch(ielc)%IA(nrdofb+1))
      allocate(GRID(Igrid)%sch(ielc)%JA(nnz_tot))
      allocate(GRID(Igrid)%sch(ielc)%SA(nnz_tot))
      GRID(Igrid)%sch(ielc)%IA(1:nrdofb+1) = IAbb(1:nrdofb+1)
      GRID(Igrid)%sch(ielc)%JA(1:nnz_tot) = JAbb(1:nnz_tot)
      GRID(Igrid)%sch(ielc)%SA(1:nnz_tot) = SAbb(1:nnz_tot)
   endif
!
!..call pardiso to solve
   call pardiso_solve(IAbb(1:nrdofb+1),JAbb(1:nnz_tot),SAbb(1:nnz_tot),'H', &
                      nnz_tot,nrdofb,nrdofi+1,Abi_ext)

   Abi = Abi_ext(:,1:nrdofi)
   bb = Abi_ext(:,nrdofi+1)


   deallocate(Abi_ext)
   deallocate(SAbb,IAbb,JAbb)



   if (Idec_sch .eq. 1) then
      allocate(GRID(igrid)%sch(ielc)%A21(nrdofb, nrdofi))
      GRID(igrid)%sch(ielc)%A21 = Abi
   endif

   allocate(Atemp(nrdofi,nrdofi))

   transa = 'N'
   matdescra(1) = 'G' ; matdescra(2) = 'L'
   matdescra(3) = 'N' ; matdescra(4) = 'F'
!
   call mkl_zcoomm('N',nrdofi,nrdofi,nrdofb,ZONE,matdescra,SAib,IAib,  &
                   JAib,nz_ib,Abi,nrdofb,ZERO,Atemp,nrdofi)
!
   Atemp = Aii - Atemp
!


   if (ISTORE .eq. ISTORE_YES) then
      allocate(A_MACRO(ielc)%array(nrdofb,nrdofi))
      A_MACRO(ielc)%array = Abi
!
      allocate(A_MACRO(ielc)%vect(nrdofb))
      A_MACRO(ielc)%vect = bb
   endif
!
   allocate(btemp(nrdofi))
!
   call mkl_zcoogemv(transa,nrdofi,SAib,IAib,JAib,nz_ib,bb,btemp)
!
   deallocate(SAib,IAib,JAib)
!
   btemp = bi - btemp
!
   zbload(1:nrdofi) = btemp(1:nrdofi)
!
   do j = 1,nrdofi
      do i = 1,j
         k = nk(i,j)
         zastif(k) = Atemp(i,j)
      enddo
   enddo

   deallocate(Abi, Aii, bb,bi, btemp, Atemp)
!
!
   end subroutine celem_macro
