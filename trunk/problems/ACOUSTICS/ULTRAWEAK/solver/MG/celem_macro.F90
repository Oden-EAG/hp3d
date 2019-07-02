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
   subroutine celem_macro(Ielc, MdleC, Nrdof_macro, Nrnod_macro, Nod_macro,  &
                          NdofH_macro, NdofE_macro, NdofV_macro, Zbload, Zastif)
!
   use data_structure3D, ONLY: NODES,MAXNODM, NRNODS,NR_PHYSA
   use physics,          ONLY: NRHVAR, NREVAR, NRVVAR, NRQVAR
   use assembly,         ONLY: NR_RHS, NEXTRACT, IDBC, ZDOFD,   &
                               BLOC, AAUX, ALOC, ZBMOD, ZAMOD
   use parameters,       ONLY: ZERO, ZONE,      &
                               MAXbrickH, MAXbrickE, MAXbrickV, MAXbrickQ
   use macro_grid_info,  ONLY: A_MACRO
!
   IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!..globals  
   integer, intent( in)    :: Ielc, MdleC
   integer, intent(out)    :: Nrdof_macro, Nrnod_macro 
   integer, intent(out)    :: Nod_macro(*),   NdofH_macro(*), &
                              NdofE_macro(*), NdofV_macro(*)
#if C_MODE      
   complex*16, intent(out) :: Zbload(*), Zastif(*)                     
#else
   real*8,     intent(out) :: Zbload(*), Zastif(*)
#endif
!   
!..locals   
!..dof offsets
   integer, allocatable    :: naH(:),naE(:), naV(:), naQ(:)
   integer, allocatable    :: mdle_macro(:)
!
!..workspace for celem
   integer, allocatable :: maxdofs(:)
   integer :: maxdofm, nrdofs(NR_PHYSA), nrdofm, nrdofc, nrnodm
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM),  &
              ndofmE(MAXNODM),ndofmV(MAXNODM), ndofmQ(MAXNODM)
!
!   
!..dof counters (i for interface, b for bubble)
   integer :: nrdof_Hi, nrdof_Ei, nrdof_Vi, nrdof_Qi
   integer :: nrdof_Hb, nrdof_Eb, nrdof_Vb, nrdof_Qb
   integer :: nrdofi, nrdofb
!
!..various
   integer :: mdle, nrmdle, iel, nod, master, ndof, ndofi, ndofb, idec
   integer :: l, lb, j, i, k1, k2, k, ii
   integer, allocatable :: lcon(:)
   character(len=4) :: type

#if C_MODE      
   complex*16 :: zvoid(1)
#else
   real*8     :: zvoid(1)
#endif
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
#if C_MODE
      complex*16, allocatable  :: val(:)
#else      
      real*8, allocatable      :: val(:)
#endif      
!
   end type sparse_matrix
!
   type(sparse_matrix), allocatable :: SKbb(:)
!
!..global matrices

#if C_MODE
   complex*16, allocatable  :: Abi(:,:), Aii(:,:), Atemp(:,:), btemp(:)
   complex*16, allocatable  :: bb(:), bi(:)
   complex*16, allocatable  :: SAbb(:),SAib(:)
   complex*16, allocatable  :: DA(:),Abi_ext(:,:)
   complex*16, allocatable  :: zload(:), ztemp(:)
#else      
   real*8,     allocatable  :: Abi(:,:), Aii(:,:), Atemp(:,:), btemp(:)
   real*8,     allocatable  :: bb(:), bi(:)
   real*8,     allocatable  :: SAbb(:), SAib(:)
   real*8,     allocatable  :: DA(:),Abi_ext(:)
   real*8,     allocatable  :: zload(:), ztemp(:)
#endif      
   integer,    allocatable  :: IAbb(:), JAbb(:), IAib(:), JAib(:)
!
   character(len=1)         :: transa
   character(len=1)         :: matdescra(4)
   integer                  :: nk, n1, n2
   nk(n1,n2) = (n2-1)*n2/2+n1
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
   allocate(naQ(NRNODS)); naQ = -1
!
!..Initialize dof counters 
   nrdof_Hi = 0 ; nrdof_Ei = 0; nrdof_Vi = 0; nrdof_Qi = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0; nrdof_Qb = 0
!   
   mdle=0 ;   Nrnod_macro = 0 ; nrmdle = 0
!
   nz_b = 0; nz_ib = 0
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
      do i =1, nrnodm
         nod = nodm(i)
         call find_master(nod, master)
         type = NODES(master)%type
!         
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
            nsum_b = nsum_b + ndofmH(i) + ndofmE(i) + ndofmV(i)+ ndofmQ(i)
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
!  ...compute offsets for H1     
      do i = 1, nrnodm
         nod = nodm(i)
         if (naH(nod).ge.0) cycle
!         
!     ...find its coarse master
         call find_master(nod, master)
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
            Nrnod_macro = Nrnod_macro + 1 
            Nod_macro(Nrnod_macro) = nod
            ndofH_macro(Nrnod_macro) = ndofmH(i)
            ndofE_macro(Nrnod_macro) = ndofmE(i)
            ndofV_macro(Nrnod_macro) = ndofmV(i)
!        ...store the first dof offset
            naH(nod) = nrdof_Hi
!        ...update the H1 dof counter
            nrdof_Hi = nrdof_Hi + ndofmH(i)
         end select
      enddo

!  ...compute offsets for H(curl)     
      do i = 1, nrnodm
         nod = nodm(i)
         if (naE(nod).ge.0) cycle
!         
!     ...find its coarse master
         call find_master(nod, master)
!         
!     ...check if the the type of the master node 
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!            
!        ...store the first dof offset
            naE(nod) = nrdof_Hb
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

!  ...compute offsets for H(div)     
      do i = 1, nrnodm
         nod = nodm(i)
         if (naV(nod).ge.0) cycle
!         
!     ...find its coarse master
         call find_master(nod, master)
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
      nod = nodm(nrnodm)
!  ...store the first dof offset
      naQ(nod) = nrdof_Qb
!  ...L2 dof are all bubbles
      nrdof_Qb = nrdof_Qb + ndofmQ(nrnodm)
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
   Nrdof_macro = nrdofi   
!
!..build and store connectivity maps
   allocate(SKbb(nrdofb))
!..initialize
   do j=1,nrdofb   
      SKbb(j)%nz_old = 0
   enddo 
!
   allocate(A_MACRO(ielc)%GLOC(nrmdle))
   A_MACRO(Ielc)%nrmdle = nrmdle
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
      call celem1(mdle_macro(iel),idec,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                  ndofmE,ndofmV,ndofmQ,nrnodm, zvoid,zvoid)
!  ...construct local to global connectivities
      l=0; lb = 0
! ....H1 dof
      do i =1, nrnodm
         nod = nodm(i)
!     ...find its coarse master
         call find_master(nod, master)
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
!  ...H(curl) dof 
      do i = 1,nrnodm
         nod = nodm(i)
!     ...find its coarse master
         call find_master(nod, master)
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
!  ...H(div) dof 
      do i =1, nrnodm
         nod = nodm(i)
!     ...find its coarse master
         call find_master(nod, master)
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
!  ...L2 dof (all bubbles)
      do j=1,ndofmQ(nrnodm)
         l=l+1; lb=lb+1
         lcon(l) = -(nrdof_Hb + nrdof_Eb + nrdof_Vb + naQ(nod)+j)
      enddo 
!
!  ...local number of dof is ndof_tot
      ndof = l
      ndofb = lb
      ndofi = l-lb
      allocate(A_MACRO(ielc)%GLOC(iel)%lcon(ndof))
      A_MACRO(ielc)%GLOC(iel)%lcon(1:ndof) = lcon(1:ndof)
      A_MACRO(ielc)%GLOC(iel)%ndof = ndof
!
!  ...loop through number of dof for the element
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
   do i = 1,nrdofb
      allocate(SKbb(i)%col(SKbb(i)%nz_old))
      allocate(SKbb(i)%val(SKbb(i)%nz_old))
!  ...reinitialize counters      
      SKbb(i)%nz_old = 0
   enddo
!
!..assemble matrices
   allocate(Abi(nrdofb,nrdofi), Aii(nrdofi,nrdofi))
   allocate(bb(nrdofb),bi(nrdofi)) 
   allocate(SAib(nz_ib),IAib(nz_ib),JAib(nz_ib));
   allocate(Abi_ext(nrdofb,1+nrdofi)) 
   ! allocate(DA(nrdofb)) 
!
!..initialize
   Abi = ZERO; Aii = ZERO ; bb = ZERO; bi = ZERO; 
   Abi_ext = ZERO 
   ! DA = ZERO
!
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
      call celem1(mdle_macro(iel),idec,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                  ndofmE,ndofmV,ndofmQ,nrnodm,zload,ztemp)
!      
!  ...assemble the global stiffness matrix and load vector
!  ...loop through element dof
      ndof = A_MACRO(ielc)%GLOC(iel)%ndof
      do k1=1,ndof
!     ...global dof is:
         i = A_MACRO(ielc)%GLOC(iel)%lcon(k1)
         if (i .lt. 0 ) then
            i = abs(i)
!        ...Assemble global load vector
            ! ii = nrdofb*nrdofi+i
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
                     ! if (i .eq. j) then   
                     !    DA(i) = DA(i) + ztemp(k)
                     ! endif   
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
   deallocate(naH, naE, naV, naQ)
   deallocate(maxdofs)

!
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

   ! call diag_scaling_s(nrdofb,nrdofi+1,nnz_tot,IAbb,JAbb,  &
   !                           DA, SAbb,Abi_ext)
! 
   call coo2csr_new(IAbb, JAbb, SAbb, nrdofb, nnz_tot)
!
! 
!..call pardiso to solve
   call pardiso_solve(IAbb(1:nrdofb+1),JAbb(1:nnz_tot),SAbb(1:nnz_tot),'H', &
                      nnz_tot,nrdofb,nrdofi+1,Abi_ext)

   ! do i = 1, nrdofb
   !    Abi_ext(i,1:nrdofi+1) = Abi_ext(i,1:nrdofi+1)/sqrt(DA(i))
   ! enddo

   ! deallocate(DA)

   Abi = Abi_ext(:,1:nrdofi)
   bb = Abi_ext(:,nrdofi+1)


   deallocate(Abi_ext)
   deallocate(SAbb,IAbb,JAbb)

   allocate(Atemp(nrdofi,nrdofi))

   transa = 'N'
   matdescra(1) = 'G' ; matdescra(2) = 'L'
   matdescra(3) = 'N' ; matdescra(4) = 'F'
!
   call mkl_zcoomm('N',nrdofi,nrdofi,nrdofb,ZONE,matdescra,SAib,IAib,  &
                   JAib,nz_ib,Abi,nrdofb,ZERO,Atemp,nrdofi)
!
   Atemp = Aii - Atemp

   allocate(A_MACRO(ielc)%array(nrdofb,nrdofi))
   A_MACRO(ielc)%array = Abi
!   
   allocate(A_MACRO(ielc)%vect(nrdofb))
   A_MACRO(ielc)%vect = bb
   A_MACRO(ielc)%ni = nrdofi
   A_MACRO(ielc)%nb = nrdofb
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
