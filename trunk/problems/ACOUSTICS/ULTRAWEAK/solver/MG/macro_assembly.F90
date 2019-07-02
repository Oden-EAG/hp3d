! -----------------------------------------------------------------------
!
!    routine name       - macro_assembly
!
! -----------------------------------------------------------------------
!
!    latest revision    - Feb 2018
!
!    purpose            - routine computes global stiffness matrix
!                         and global load vector corresponding to a
!                         macro-grid 
!    out :              - the global stiffness matrix and load vector 
!                         are stored in sparse coordinate form in the
!                         module sparse_solve
!
!-----------------------------------------------------------------------
!
   subroutine macro_assembly
!
   use data_structure3D,   ONLY: NRNODS, NRELES
   use m_assembly,         ONLY: NFIRSTH, NFIRSTE, NFIRSTV, NRDOF_TOT, NRDOF_CON
   use macro_grid_info,    ONLY: MDLE_MACRO, NRELES_COARSE, MACRO_ELEM,  &
                                 DLOC, A_MACRO, ZSOL_C, NRDOF_MACRO
   use parameters,         ONLY: ZERO
   use assembly,           ONLY: NR_RHS
   use pcg_info,           ONLY: RHS, XSOL, TOL, MAX_ITER
   use patch_info,         ONLY: deallocate_patches

   IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!..dof counters
   integer :: nrdof_H, nrdof_E, nrdof_V, ndof_macro
   integer :: iel,i,l,j
!
!..locals for with macro_elem_info
   integer :: mdle, nrnod_macro, maxndof, nod
    
   integer, allocatable, save :: nod_macro(:),   ndofH_macro(:),  &
                                 ndofE_macro(:), ndofV_macro(:) 
!$omp threadprivate(nod_macro,ndofH_macro,ndofE_macro,ndofV_macro)

   integer, allocatable, save :: lcon(:)
!$omp threadprivate(lcon)                                 
#if C_MODE
   complex*16, allocatable, save :: zsol_macro(:)
#else
   real*8,     allocatable, save :: zsol_macro(:)
#endif
!$omp threadprivate(zsol_macro)                                 


   integer :: k, k1, k2, ip, ndof, nz
   integer :: index_nz(NRELES_COARSE), index_inz(NRELES_COARSE)
   real*8  :: finish, start, omp_get_wtime
!
!-----------------------------------------------------------------------
!
   start = omp_get_wtime()

!..allocate and initialize offsets
   allocate(NFIRSTH(NRNODS)); NFIRSTH = -1
   allocate(NFIRSTE(NRNODS)); NFIRSTE = -1
   allocate(NFIRSTV(NRNODS)); NFIRSTV = -1
!
!..initialize dof counters
   nrdof_H = 0 ; nrdof_E = 0; nrdof_V = 0;
   maxndof = 0; nz = 0
!
!..loop through the coarse elements
   do iel = 1, NRELES_COARSE
!
      ndof_macro = MACRO_ELEM(iel)%ndof_macro
      nrnod_macro = MACRO_ELEM(iel)%nrnod_macro
      allocate(nod_macro(nrnod_macro))  ; nod_macro = MACRO_ELEM(iel)%nod_macro
      allocate(ndofH_macro(nrnod_macro)); ndofH_macro = MACRO_ELEM(iel)%ndofH_macro
      allocate(ndofE_macro(nrnod_macro)); ndofE_macro = MACRO_ELEM(iel)%ndofE_macro
      allocate(ndofV_macro(nrnod_macro)); ndofV_macro = MACRO_ELEM(iel)%ndofV_macro
!
      maxndof = max0(maxndof,ndof_macro)
!
!  ...counting nz for OMP
      index_nz(iel)  = ndof_macro**2
      nz = nz + index_nz(iel)
      index_inz(iel) = sum(index_nz(1:iel-1))
!
!  ...compute offsets for H1 dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (NFIRSTH(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTH(nod) = nrdof_H
!     ...update the H1 dof counter
         nrdof_H = nrdof_H + ndofH_macro(i)
      enddo
!
!  ...compute offsets for H(curl) dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (NFIRSTE(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTE(nod) = nrdof_E
!     ...update the H(curl) dof counter
         nrdof_E = nrdof_E + ndofE_macro(i)
      enddo
!
!  ...compute offsets for H(div) dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (NFIRSTV(nod).ge.0) cycle
!     ...store the first dof offset
         NFIRSTV(nod) = nrdof_V
!     ...update the H(div) dof counter
         nrdof_V = nrdof_V + ndofV_macro(i)
      enddo
!
      deallocate(nod_macro,ndofH_macro, ndofE_macro, ndofV_macro)
!..end of loop through the coarse elements
   enddo
!
   NRDOF_MACRO = nrdof_H + nrdof_E + nrdof_V
!
!..allocate memory for sparse matrix storage
   allocate(RHS(NRDOF_MACRO)) ; RHS = ZERO
!
   allocate(A_MACRO(NRELES_COARSE))
   allocate(DLOC(NRELES_COARSE))
!
!
!$omp parallel default(shared)    &
!$omp private(iel,mdle,nrnod_macro, ndof_macro,l,i,nod,j,ndof)
!$omp do schedule(guided)
   do iel=1,NRELES_COARSE
      mdle = MDLE_MACRO(iel)
!
      nrnod_macro = MACRO_ELEM(iel)%nrnod_macro
      ndof_macro  = MACRO_ELEM(iel)%ndof_macro
      allocate(nod_macro(nrnod_macro))  
      allocate(ndofH_macro(nrnod_macro))
      allocate(ndofE_macro(nrnod_macro))
      allocate(ndofV_macro(nrnod_macro))
      allocate(DLOC(iel)%zstiff(ndof_macro*(ndof_macro+1)/2))
      allocate(DLOC(iel)%zbload(ndof_macro))

!    
      call celem_macro(iel,mdle,ndof_macro,nrnod_macro,nod_macro, ndofH_macro,  &
                       ndofE_macro, ndofV_macro, DLOC(iel)%zbload, DLOC(iel)%zstiff)    
!
      DLOC(iel)%ndof = ndof_macro
      allocate(DLOC(iel)%lcon(ndof_macro))

!  ...create connectivity maps
      l=0
!  ...H1 dof 
      do i = 1,nrnod_macro
         nod = nod_macro(i)
         do j=1,ndofH_macro(i)
            l=l+1
            DLOC(iel)%lcon(l) = NFIRSTH(nod)+j
         enddo
      enddo
!  ...H(curl) dof 
      do i =1,nrnod_macro
         nod = nod_macro(i)
         do j=1,ndofE_macro(i)
            l=l+1
            DLOC(iel)%lcon(l) = nrdof_H + NFIRSTE(nod)+j
         enddo
      enddo
!  ...H(div) dof 
      do i = 1,nrnod_macro
         nod = nod_macro(i)
         do j=1,ndofV_macro(i)
            l=l+1
            DLOC(iel)%lcon(l) = nrdof_H + nrdof_E + NFIRSTV(nod)+j
         enddo
      enddo
!      
      deallocate(nod_macro,ndofH_macro, ndofE_macro, ndofV_macro)
!      
!..assembly finished
   enddo
!$omp end do   
!$omp end parallel
!
   deallocate(NFIRSTH, NFIRSTE, NFIRSTV)

   finish = omp_get_wtime()
   write(*,1002) finish-start
1002 format(' Assemble macro-grid system:t=     ', es13.4)    

!
   write(*,*) 'Assemble patches...'
   start = omp_get_wtime()
   call patch_assembly
   finish = omp_get_wtime()
   write(*,1003) finish-start
1003 format(' Assemble patches;t=               ', es13.4)    
!
   TOL = 1e-2 ;  MAX_ITER = 20
   allocate(XSOL(NRDOF_MACRO)) ; XSOL = 0
!
   call prolong_solution
!
   do iel = 1, NRELES_COARSE
      ndof = DLOC(iel)%ndof
!
      do k1=1,ndof
         i = DLOC(iel)%lcon(k1)
         XSOL(i) = ZSOL_C(iel)%macro(k1)
      enddo
!..end of loop through macro elements
   enddo
!
!
   write(*,*) 'Solve the system using PCG...'
   start = omp_get_wtime()
   call mg_pcg_solve

   finish = omp_get_wtime()
   write(*,1004) finish-start
1004 format(' Solve the system using PCG:t=     ', es13.4)    
! 
   call deallocate_patches
!
!..end of solution step
!
!..store the solution into the data structure
!
!..loop through macro elements
!$omp parallel default(shared)    &
!$omp private(iel,ndof,k1,i)
!$omp do schedule(guided)
   do iel = 1, NRELES_COARSE
      ndof = DLOC(iel)%ndof
!
      allocate(zsol_macro(ndof)) 
!
      do k1=1,ndof
         i = DLOC(iel)%lcon(k1)
         zsol_macro(k1) = XSOL(i)
      enddo
!
      ! zsol_macro = ZSOL_C(iel)%macro
!
      call solout_macro(iel,zsol_macro,ndof)
      deallocate(A_MACRO(iel)%GLOC)
      deallocate(ZSOL_C(iel)%macro)
      deallocate(zsol_macro) 
      deallocate(DLOC(iel)%lcon)
      deallocate(DLOC(iel)%zstiff)
      deallocate(DLOC(iel)%zbload)
!..end of loop through macro elements
   enddo
!$omp end do   
!$omp end parallel
!
   deallocate(RHS,XSOL,A_MACRO,DLOC)
!  
!
   end subroutine macro_assembly
