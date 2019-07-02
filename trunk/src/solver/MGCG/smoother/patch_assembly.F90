!
! -----------------------------------------------------------------------
!
!    routine name       - patch_assembly
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - routine assembles global matrix blocks 
!                         for the additive Schwarz preconditioner
!                         LU (or Cholesky) decompositions of blocks
!                         are stored in module patch_info
!
! ----------------------------------------------------------------------
!
   subroutine patch_assembly(Igrid)
!
   use data_structure3D,  only: NRNODS
   ! use patch_info     
   use mg_data_structure, only: GRID, NODES_MG
   use parameters,        only: ZERO


   IMPLICIT NONE
!
   integer,intent(in) :: Igrid 
   integer :: ip, inod,jel, iel, i,j, nloc, mdle, nod, l, lp
   integer :: k1, k2, k, ij, kk, nk
   integer :: nrnodp, nrmdle, ndof_macro, nrnod_macro
   integer :: nrdof_H, nrdof_E, nrdof_V, nrdof
   integer, allocatable :: npH(:), npE(:), npV(:), nvisit(:)
! $omp threadprivate(npH,npE,npV,nvisit)   
   integer, allocatable :: nodl_aux(:), nod_macro(:)
! $omp threadprivate(nodl_aux, nod_macro)   
   integer, allocatable :: ndofH_macro(:), ndofE_macro(:), ndofV_macro(:)
! $omp threadprivate(ndofH_macro,ndofE_macro,ndofV_macro)
   integer, allocatable :: lcon(:)
! $omp threadprivate(lcon)
!
   integer :: iprint, info
!..function for vector storage of a symmetric (Hermitian) matrix
   nk(k1,k2) = (k2-1)*k2/2+k1


!----------------------------------------------------------------------------
!
   iprint = 0
!
!..loop through the patches
!$omp parallel default(shared)   &
!$omp private(ip, nrnodp, inod,nrdof_H, nrdof_E, nrdof_V, nrmdle, jel, i, j, &
!$omp         mdle, iel, ndof_macro, nrnod_macro, nod, nloc, nrdof, l, lp, &
!$omp         npH,npE,npV,nvisit,nodl_aux,nod_macro,lcon, &
!$omp         ndofH_macro, ndofE_macro, ndofV_macro) 
!..allocate and initialize offsets
   allocate(npH(NRNODS))
   allocate(npE(NRNODS))
   allocate(npV(NRNODS))
   allocate(nvisit(NRNODS))
!$omp do schedule(guided)
   do ip = 1, GRID(Igrid)%nrpatch
!  ...get number of nodes in the patch
      nrnodp = GRID(Igrid)%ptch(ip)%nsz      
!  ...allocate auxiliary nod list for each patch
      allocate(nodl_aux(nrnodp))   
!  ...nod counter
      inod = 0
!  ...initialize offsets for the patch
      npH = -1; npE = -1; npV = -1; nvisit = 0
!  ...initialize dof counter 
      nrdof_H  = 0 ; nrdof_E  = 0; nrdof_V  = 0 
!  ...number of mdle nodes contributing to the patch
      nrmdle = GRID(Igrid)%ptch(ip)%nrmdle
!      
!  ...loop through the mdle nodes within the patch   
      do jel = 1,nrmdle
!     ...pick up the mdle node number
         mdle = GRID(Igrid)%ptch(ip)%mdlel(jel)%mdle
!         
!     ...get info for this macro-element
         iel = NODES_MG(mdle)%iel(Igrid)
         ndof_macro  = GRID(Igrid)%macro_elem(iel)%ndof_macro
         nrnod_macro = GRID(Igrid)%macro_elem(iel)%nrnod_macro
!
         allocate(nod_macro(nrnod_macro))   
         nod_macro   = GRID(Igrid)%macro_elem(iel)%nod_macro
!         
         allocate(ndofH_macro(nrnod_macro)) 
         ndofH_macro = GRID(Igrid)%macro_elem(iel)%ndofH_macro
!         
         allocate(ndofE_macro(nrnod_macro)) 
         ndofE_macro = GRID(Igrid)%macro_elem(iel)%ndofE_macro
!         
         allocate(ndofV_macro(nrnod_macro)) 
         ndofV_macro = GRID(Igrid)%macro_elem(iel)%ndofV_macro
!         
!         
!     ...construct auxiliary nod list that respects the node ordering
         do i =1, nrnod_macro
            nod = nod_macro(i)
!        ...avoid repetition
            if (nvisit(nod) .gt. 0) cycle
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
               inod = inod + 1
               nodl_aux(inod) = nod
               nvisit(nod) = 1
            else
               nvisit(nod) = 1
            endif         
         enddo

!     ...construct offsets for the patch
!     ...H1 dof
         do i =1, nrnod_macro
            nod = nod_macro(i)
!        ...avoid repetition
            if (npH(nod).ge.0) cycle
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
!           ...store the first dof offset
               npH(nod) = nrdof_H
!           ...update the H1 dof counter
               nrdof_H = nrdof_H + ndofH_macro(i)
            endif         
         enddo
!
!     ...H(curl) dof
         do i = 1,nrnod_macro
            nod = nod_macro(i)
!        ...avoid repetition
            if (npE(nod).ge.0) cycle
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
!           ...store the first dof offset
               npE(nod) = nrdof_E
!           ...update the H(curl) dof counter
               nrdof_E = nrdof_E + ndofE_macro(i)
            endif   
         enddo
!
!     ...H(div) dof
         do i=1,nrnod_macro
            nod = nod_macro(i)
!        ...avoid repetition
            if (npV(nod).ge.0) cycle
!        ...store the first dof offset
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
               npV(nod) = nrdof_V
!           ...update the H(div) dof counter
               nrdof_V = nrdof_V + ndofV_macro(i)
            endif   
         enddo
!
         deallocate(nod_macro)
         deallocate(ndofH_macro, ndofE_macro, ndofV_macro)
!
!     ...end of loop through mdle nodes within the patch         
      enddo
      nrdof = nrdof_H + nrdof_E + nrdof_V
!  ...total dof for the patch       
      GRID(Igrid)%ptch(ip)%nrdof = nrdof
!      
!  ...overwrite nod list for the patch to change order
!  ...the node ordering has to agree with local dof 
      GRID(Igrid)%ptch(ip)%nodl(1:nrnodp) = nodl_aux(1:nrnodp)
!
      deallocate(nodl_aux)
!
!  ...construct and store connectivity arrays for each element in the patch
!
!  ...loop through the mdle nodes within the patch
      do jel = 1,nrmdle
!     ...pick up the mdle node number
         mdle = GRID(Igrid)%ptch(ip)%mdlel(jel)%mdle
!         
!     ...get info for this macro-element
         iel = NODES_MG(mdle)%iel(Igrid)
         ndof_macro  = GRID(Igrid)%macro_elem(iel)%ndof_macro
         nrnod_macro = GRID(Igrid)%macro_elem(iel)%nrnod_macro
!
         allocate(nod_macro(nrnod_macro))   
         nod_macro   = GRID(Igrid)%macro_elem(iel)%nod_macro
!         
         allocate(ndofH_macro(nrnod_macro)) 
         ndofH_macro = GRID(Igrid)%macro_elem(iel)%ndofH_macro
!         
         allocate(ndofE_macro(nrnod_macro)) 
         ndofE_macro = GRID(Igrid)%macro_elem(iel)%ndofE_macro
!         
         allocate(ndofV_macro(nrnod_macro)) 
         ndofV_macro = GRID(Igrid)%macro_elem(iel)%ndofV_macro
!         
!
!     ...allocate connectivity array
         allocate(lcon(ndof_macro))
         l = 0; lp = 0
!     ...H1 dof
         do i =1, nrnod_macro
            nod = nod_macro(i)
!        ...check if the nod is in the patch            
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
               do j = 1, ndofH_macro(i)
                  l = l + 1; lp = lp+1
                  lcon(l) = npH(nod) + j
               enddo
            else
               do j = 1, ndofH_macro(i)
                  l = l + 1
                  lcon(l) = -1
               enddo
            endif
         enddo   
!         
!     ...H(curl) dof
         do i =1, nrnod_macro
            nod = nod_macro(i)
!        ...check if the nod is in the patch            
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
               do j = 1, ndofE_macro(i)
                  l = l + 1; lp = lp+1
                  lcon(l) = nrdof_H + npE(nod) + j
               enddo
            else
               do j = 1, ndofE_macro(i)
                  l = l + 1
                  lcon(l) = -1
               enddo
            endif
         enddo   
!     ...H(div) dof
         do i =1, nrnod_macro
            nod = nod_macro(i)
!        ...check if the nod is in the patch            
            call locate(nod,GRID(Igrid)%ptch(ip)%nodl,nrnodp,nloc)
            if (nloc .ne. 0) then
               do j = 1, ndofV_macro(i)
                  l = l + 1; lp = lp+1
                  lcon(l) = nrdof_H + nrdof_E + npV(nod) + j
               enddo
            else
               do j = 1, ndofV_macro(i)
                  l = l + 1
                  lcon(l) = -1
               enddo
            endif
         enddo   

         GRID(Igrid)%ptch(ip)%mdlel(jel)%ndof = lp
!     ...store connectivity map for this mdle in the patch 
         allocate(GRID(Igrid)%ptch(ip)%mdlel(jel)%lcon(l))
         GRID(Igrid)%ptch(ip)%mdlel(jel)%lcon = lcon
         deallocate(lcon)
         deallocate(nod_macro)
         deallocate(ndofH_macro, ndofE_macro, ndofV_macro)
!  ...end of loop through mdle nodes with in the patch
      enddo
!
!..end of loop through the patches
   enddo
!$omp end do
  deallocate(npH,npE,npV,nvisit)
!$omp end parallel

!
!..finally the assembly of the patches
!$omp parallel default(shared)   &
!$omp private(ip, nrdof, nrmdle, jel, mdle, iel, ndof_macro, k1, k2, i, j,  &
!$omp         kk, ij, info,lcon ) 
!$omp do schedule(guided)
!..loop through the patches
   do ip = 1, GRID(Igrid)%nrpatch
!      
!  ...total number of dof in the patch      
      nrdof = GRID(Igrid)%ptch(ip)%nrdof
! 
      allocate(GRID(Igrid)%ptch(ip)%zAp(nrdof*(nrdof+1)/2)) 
      GRID(Igrid)%ptch(ip)%zAp = ZERO
!      
!  ...number of mdle nodes contributing to the patch
      nrmdle = GRID(Igrid)%ptch(ip)%nrmdle
!  ...loop through the mdle nodes within the patch   
      do jel = 1,nrmdle
!
!     ...pick up the mdle node number   
         mdle = GRID(Igrid)%ptch(ip)%mdlel(jel)%mdle
!         
!     ...get info for this macro-element
         iel = NODES_MG(mdle)%iel(Igrid)
         ndof_macro  = GRID(Igrid)%macro_elem(iel)%ndof_macro
!
!     ...allocate connectivity array
         allocate(lcon(ndof_macro))
!
         lcon = GRID(Igrid)%ptch(ip)%mdlel(jel)%lcon
         deallocate(GRID(Igrid)%ptch(ip)%mdlel(jel)%lcon)
!         
!     ...loop through dof in the matrix and assemble the ones corresponding to 
!     ...the specific patch
!--------------------------------------------------------------------------------
!--------->>>>>>FOR OPRIMIZATION: ORDERING OF LOOPS HAS TO CHANGE<<<<<<----------
!--------------------------------------------------------------------------------
!     ...loop through macro element dof
         do k1 = 1,ndof_macro
            i = lcon(k1)
!        ...check the global index (if it is negative then it is not in the patch)            
            if (i .gt. 0 ) then
!           ...loop through macro element dof
               do k2=1,ndof_macro
                  j = lcon(k2)
!              ...check the global index (if it is negative then it is not in the patch)            
                  if (j .gt. 0 ) then
                     if (k2 .ge. k1 ) then
                        kk = nk(k1,k2)       
                     else
                        kk = nk(k2,k1)
                     endif
                     if (i .le. j) then
                        ij = nk(i,j)
#if C_MODE
                        if (k2 .ge. k1 ) then
                           GRID(Igrid)%ptch(ip)%zAp(ij) = GRID(Igrid)%ptch(ip)%zAp(ij) &
                                                        + GRID(Igrid)%dloc(iel)%zstiff(kk)
                        else
                           GRID(Igrid)%ptch(ip)%zAp(ij) = GRID(Igrid)%ptch(ip)%zAp(ij) &
                                                        + conjg(GRID(Igrid)%dloc(iel)%zstiff(kk))
                        endif
#else                                                  
                        GRID(Igrid)%ptch(ip)%zAp(ij) = GRID(Igrid)%ptch(ip)%zAp(ij) &
                                                     + GRID(Igrid)%dloc(iel)%zstiff(kk)
#endif                                                                
                     endif
                  endif
               enddo
            endif
         enddo   
!
         deallocate(lcon)
!
!  ...end of loop through mdle nodes with in the patch
      enddo
!
#if C_MODE      
      call ZPPTRF('U', nrdof, GRID(Igrid)%ptch(ip)%zAp, info)
#else
      call DPPTRF('U', nrdof, GRID(Igrid)%ptch(ip)%zAp, info)
#endif      
      if (info .ne. 0) then 
         write(*,*) 'patch_assembly: info = ', info
         write(*,*) 'patch_assembly: ip   = ', ip
         ! write(*,*) 'mdle nodes  = ', PTCH(ip)%mdlel(1:nrmdle)%mdle
         ! write(*,*) 'nodl = ', PTCH(ip)%nodl(1:PTCH(ip)%nsz)
         ! write(*,*) 'patch list    = ', PTCH(ip)%nodl(1:PTCH(ip)%nsz)
         stop 
         ! call graphb
      endif
!..end of loop through the patches
   enddo
!$omp end do   
!$omp end parallel

!
!..create connectivity arrays between patch and global dof
   call patch_connect_maps(Igrid)
!
!
   end subroutine patch_assembly
!
!
!-----------------------------------------------------------------------
!
!    routine name       - patch_connect_maps
!
!-----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - routine computes and stores connectivity 
!                         maps between local and global level for 
!                         each patch
!
!----------------------------------------------------------------------
!
   subroutine patch_connect_maps(Igrid)

   use data_structure3D,   only: NRNODS
   use mg_data_structure,  only: NODES_MG, GRID
   use patch_info
!
   IMPLICIT NONE
!
   integer, intent(in)  :: Igrid
   integer, allocatable :: nod_macro(:)
   integer, allocatable :: ndofH_macro(:), ndofE_macro(:), ndofV_macro(:)
   integer, allocatable :: npH(:), npE(:), npV(:)
   integer              :: nrdof_H, nrdof_E, nrdof_V
   integer              :: iel,i, k, l, j, nod
   integer              :: nrdof, nrnod_macro, nrnodp
!
!-----------------------------------------------------------------------
!
!..Step 1 : compute offsets
!..allocate and initialize offsets
   allocate (npH(NRNODS)); npH = -1
   allocate (npE(NRNODS)); npE = -1
   allocate (npV(NRNODS)); npV = -1
!..determine the first dof offsets for active nodes
   nrdof_H = 0 ; nrdof_E = 0; nrdof_V = 0;
!
   do i = 1, NRNODS
      allocate(NODES_MG(i)%nod_ndof(3))
   enddo
!  
!..loop through the elements
   do iel=1,GRID(Igrid)%nreles
!
      nrnod_macro = GRID(Igrid)%macro_elem(iel)%nrnod_macro
!
      allocate(nod_macro(nrnod_macro))  
      nod_macro   = GRID(Igrid)%macro_elem(iel)%nod_macro
!      
      allocate(ndofH_macro(nrnod_macro))
      ndofH_macro = GRID(Igrid)%macro_elem(iel)%ndofH_macro
!      
      allocate(ndofE_macro(nrnod_macro))
      ndofE_macro = GRID(Igrid)%macro_elem(iel)%ndofE_macro
!      
      allocate(ndofV_macro(nrnod_macro))
      ndofV_macro = GRID(Igrid)%macro_elem(iel)%ndofV_macro
!      
!
!  ...H1 dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!         
         NODES_MG(nod)%nod_ndof(1) = ndofH_macro(i)
         NODES_MG(nod)%nod_ndof(2) = ndofE_macro(i)
         NODES_MG(nod)%nod_ndof(3) = ndofV_macro(i)
!
!     ...avoid repetition
         if (npH(nod).ge.0) cycle
!     ...store the first dof offset
         npH(nod) = nrdof_H
!     ...update the H1 dof counter
         nrdof_H = nrdof_H + ndofH_macro(i)
      enddo
!
!  ...H(curl) dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (npE(nod).ge.0) cycle
!     ...store the first dof offset
         npE(nod) = nrdof_E
!     ...update the H(curl) dof counter
         nrdof_E = nrdof_E + ndofE_macro(i)
      enddo
!
!  ...H(div) dof
      do i=1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (npV(nod).ge.0) cycle
!     ...store the first dof offset
         npV(nod) = nrdof_V
!     ...update the H(div) dof counter
         nrdof_V = nrdof_V + ndofV_macro(i)
      enddo
!
      deallocate(nod_macro, ndofH_macro, ndofE_macro, ndofV_macro)
!..end of loop through the elements
   enddo
!
!..Step 2: Compute Destination vectors for each patch
   do k = 1,GRID(Igrid)%nrpatch
      nrnodp = GRID(Igrid)%ptch(k)%nsz
      nrdof = GRID(Igrid)%ptch(k)%nrdof
      allocate(GRID(Igrid)%ptch(k)%lcon(nrdof))
      l=0
!  ...H1 dof 
      do i = 1,nrnodp
         nod = GRID(Igrid)%ptch(k)%nodl(i)
         do j=1,NODES_MG(nod)%nod_ndof(1)
            l=l+1
            GRID(Igrid)%ptch(k)%lcon(l) = npH(nod)+j
         enddo
      enddo
!  ...H(curl) dof 
      do i =1,nrnodp
         nod = GRID(Igrid)%ptch(k)%nodl(i)
         do j=1,NODES_MG(nod)%nod_ndof(2)
            l=l+1
            GRID(Igrid)%ptch(k)%lcon(l) = nrdof_H + npE(nod)+j
         enddo
      enddo
!        
!  ...H(div) dof 
      do i = 1,nrnodp
         nod = GRID(Igrid)%ptch(k)%nodl(i)
         do j=1,NODES_MG(nod)%nod_ndof(3)
            l=l+1
            GRID(Igrid)%ptch(k)%lcon(l) = nrdof_H + nrdof_E + npV(nod)+j
         enddo
      enddo
!  ...number of dof in a patch is ndof
      GRID(Igrid)%ptch(k)%nrdof = l
   enddo
!
   do i = 1, NRNODS
      deallocate(NODES_MG(i)%nod_ndof)
   enddo
!
   deallocate(npH,npE,npV)
!
!
   end subroutine patch_connect_maps




