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
!    in :             
!               Igrid   - Mesh index 
!
!           Idec_sch    - flag to store or not store data for the 
!                         extension operator (Schur complements)
!                         0: Compute A_MACRO structure for the finest Mesh
!                         1: Compute Schur complement too
!                         
!
!-----------------------------------------------------------------------
!
#include 'implicit_none.h'
!
   subroutine macro_assembly(Igrid, Idec_sch)
!
   use data_structure3D,   ONLY: NRNODS, NRELES
   use mg_data_structure,  ONLY: GRID, NODES_MG, mg_reset_visit
   use macro_grid_info,    ONLY: NRDOF_MACRO, A_MACRO
   use assembly_sc,        ONLY: NRDOF_TOT, NRDOF_CON,IPRINT_TIME


   IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
   integer, intent(in) :: Igrid, Idec_sch
!..dof counters
   integer :: nrdof_H, nrdof_E, nrdof_V, ndof_macro
   integer :: iel,i,l,j
!
!..locals for with macro_elem_info
   integer :: mdle, nrnod_macro, maxndof, nod
   integer :: iprint_temp
!..offsets   
   integer, allocatable :: naH(:), naE(:), naV(:)
   integer, allocatable, save :: nod_macro(:),   ndofH_macro(:),  &
                                 ndofE_macro(:), ndofV_macro(:) 
!$omp threadprivate(nod_macro,ndofH_macro,ndofE_macro,ndofV_macro)
!
!-----------------------------------------------------------------------
!
   if (Idec_sch .eq. 1) then
      call mg_reset_visit

!  ...save the natural order of elements to the visitation flag  &
      mdle = 0
      do iel=1,nreles
         call nelcon(mdle,mdle)
         NODES_MG(mdle)%visit = iel
      enddo   
   endif   

!..allocate and initialize offsets
   allocate(naH(NRNODS)); naH = -1
   allocate(naE(NRNODS)); naE = -1
   allocate(naV(NRNODS)); naV = -1
!
!..initialize dof counters
   nrdof_H = 0 ; nrdof_E = 0; nrdof_V = 0;
   maxndof = 0
!
!..loop through the coarse elements
   do iel = 1, GRID(Igrid)%nreles
!
      ndof_macro = GRID(Igrid)%macro_elem(iel)%ndof_macro
      nrnod_macro = GRID(Igrid)%macro_elem(iel)%nrnod_macro
      allocate(nod_macro(nrnod_macro),ndofH_macro(nrnod_macro))  ; 
      allocate(ndofE_macro(nrnod_macro),ndofV_macro(nrnod_macro)); 
!
      nod_macro   = GRID(Igrid)%macro_elem(iel)%nod_macro
      ndofH_macro = GRID(Igrid)%macro_elem(iel)%ndofH_macro
      ndofE_macro = GRID(Igrid)%macro_elem(iel)%ndofE_macro
      ndofV_macro = GRID(Igrid)%macro_elem(iel)%ndofV_macro

      maxndof = max0(maxndof,ndof_macro)
!
!  ...compute offsets for H1 dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (naH(nod).ge.0) cycle
!     ...store the first dof offset
         naH(nod) = nrdof_H
!     ...update the H1 dof counter
         nrdof_H = nrdof_H + ndofH_macro(i)
      enddo
!
!  ...compute offsets for H(curl) dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (naE(nod).ge.0) cycle
!     ...store the first dof offset
         naE(nod) = nrdof_E
!     ...update the H(curl) dof counter
         nrdof_E = nrdof_E + ndofE_macro(i)
      enddo
!
!  ...compute offsets for H(div) dof
      do i = 1,nrnod_macro
         nod = nod_macro(i)
!     ...avoid repetition
         if (naV(nod).ge.0) cycle
!     ...store the first dof offset
         naV(nod) = nrdof_V
!     ...update the H(div) dof counter
         nrdof_V = nrdof_V + ndofV_macro(i)
      enddo
!
      deallocate(nod_macro,ndofH_macro, ndofE_macro, ndofV_macro)
!..end of loop through the coarse elements
   enddo
!
   NRDOF_MACRO(Igrid) = nrdof_H + nrdof_E + nrdof_V
   NRDOF_CON   = NRDOF_MACRO(Igrid)
   NRDOF_TOT   = NRDOF_TOT + NRDOF_CON
!
!..allocate memory for sparse matrix storage
   allocate(GRID(Igrid)%dloc(GRID(Igrid)%nreles))
!
   allocate(A_MACRO(GRID(Igrid)%nreles))
   if(Idec_sch .eq. 1) then   
      allocate(GRID(Igrid)%sch(GRID(igrid)%nreles))
   endif   
!
!
!..supress printing statements 
   iprint_temp = IPRINT_TIME
   IPRINT_TIME = 0
!   
!$omp parallel default(shared)    &
!$omp private(iel,mdle,nrnod_macro, ndof_macro,l,i,nod,j)
!$omp do schedule(guided)
   do iel=1,GRID(Igrid)%nreles
      mdle = GRID(Igrid)%mdlel(iel)
!
      nrnod_macro = GRID(Igrid)%macro_elem(iel)%nrnod_macro
      ndof_macro  = GRID(Igrid)%macro_elem(iel)%ndof_macro
      allocate(nod_macro(nrnod_macro))  
      allocate(ndofH_macro(nrnod_macro))
      allocate(ndofE_macro(nrnod_macro))
      allocate(ndofV_macro(nrnod_macro))
      allocate(GRID(Igrid)%dloc(iel)%zstiff(ndof_macro*(ndof_macro+1)/2))
      allocate(GRID(Igrid)%dloc(iel)%zbload(ndof_macro))
!    
      call celem_macro(Igrid,Idec_sch,iel,mdle,ndof_macro,nrnod_macro,nod_macro,  &
                       ndofH_macro, ndofE_macro, ndofV_macro,      &
                       GRID(Igrid)%dloc(iel)%zbload, GRID(Igrid)%dloc(iel)%zstiff)    
!
      GRID(Igrid)%dloc(iel)%ndof = ndof_macro
      allocate(GRID(Igrid)%dloc(iel)%lcon(ndof_macro))
!
!  ...create connectivity maps
      l=0
!  ...H1 dof 
      do i = 1,nrnod_macro
         nod = nod_macro(i)
         do j=1,ndofH_macro(i)
            l=l+1
            GRID(Igrid)%dloc(iel)%lcon(l) = naH(nod)+j
         enddo
      enddo
!  ...H(curl) dof 
      do i =1,nrnod_macro
         nod = nod_macro(i)
         do j=1,ndofE_macro(i)
            l=l+1
            GRID(Igrid)%dloc(iel)%lcon(l) = nrdof_H + naE(nod)+j
         enddo
      enddo
!  ...H(div) dof 
      do i = 1,nrnod_macro
         nod = nod_macro(i)
         do j=1,ndofV_macro(i)
            l=l+1
            GRID(Igrid)%dloc(iel)%lcon(l) = nrdof_H + nrdof_E + naV(nod)+j
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
   IPRINT_TIME = iprint_temp
!
   deallocate(naH, naE, naV)
!
!
   end subroutine macro_assembly
