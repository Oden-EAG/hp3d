!-----------------------------------------------------------------------
!
!    routine name       - solout_macro
!
!-----------------------------------------------------------------------
!
!    latest revision    - SEP 16
!
!    purpose            - performs backward substitution for the
!                         interior dofs of a macro element and
!                         store the solution in the data structure
!
!   arguments :
!     in:
!             ielc    - coarse element number
!             Zsol    - coarse skeleton solution dof
!             ndof    - size of macro_element solution vector
! ----------------------------------------------------------------------
!
#include "typedefs.h"
!
   subroutine solout_macro(Igrid, Ielc,Zsol_globi,Ndof)
!
   use assembly,          ONLY: assembly_begin_par, assembly_end_par
   use parameters,        ONLY: ZERO, ZONE
   use macro_grid_info,   ONLY: A_MACRO
   use mg_data_structure, ONLY: ISTORE, ISTORE_YES
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer :: Igrid
   VTYPE              :: Zsol_globi(Ndof)
   VTYPE, allocatable :: macro_array(:,:), macro_vector(:)
   VTYPE, allocatable :: Zsol_globb(:), zsol_loci(:), Zsol_locb(:), zsol_loc(:)
   integer            :: jb, ji, mdle, i, iel, k1, ielc
   integer            :: Ndof, ndof_c,ndof_tot, maxdofm
   integer            :: ni
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   jb = A_MACRO(ielc)%nb
   ji = A_MACRO(ielc)%ni
!
   allocate(Zsol_globb(jb))
!
   if (ISTORE .eq. ISTORE_YES) then
      if (jb .ne. 0) then
!
#if C_MODE
         call ZGEMV('N',jb,ji,ZONE,A_MACRO(ielc)%array,jb,Zsol_globi,1,ZERO,Zsol_globb,1)
#else
         call DGEMV('N',jb,ji,ZONE,A_MACRO(ielc)%array,jb,Zsol_globi,1,ZERO,Zsol_globb,1)
#endif
         Zsol_globb = A_MACRO(ielc)%vect - Zsol_globb
!
         deallocate(A_MACRO(ielc)%array)
         deallocate(A_MACRO(ielc)%vect)
      endif

   else
!  ...recompute matrices
      allocate(macro_array(jb, ji))
      allocate(macro_vector(jb))

      !call get_macro_arrays(Igrid, Ielc, macro_array, macro_vector, jb, ji)
      write(*,*) 'solout_macro: ISTORE_YES required. stop.'
      stop
      if (jb .ne. 0) then

#if C_MODE
         call ZGEMV('N',jb,ji,ZONE,macro_array,jb,Zsol_globi,1,ZERO,Zsol_globb,1)
#else
         call DGEMV('N',jb,ji,ZONE,macro_array,jb,Zsol_globi,1,ZERO,Zsol_globb,1)
#endif
         Zsol_globb = macro_vector - Zsol_globb
!
         deallocate(macro_array, macro_vector)

      endif
!
   endif
!
!..allocation of arrays required from solout_macro
   call assembly_begin_par
!
   do iel = 1, A_MACRO(ielc)%nrmdle
!
      ni = A_MACRO(ielc)%GLOC(iel)%ni
!
      allocate(zsol_loci(ni))
!
      do k1=1,ni
!     ...global dof is:
         i = A_MACRO(ielc)%GLOC(iel)%lcon(k1)
         if (i .lt. 0) then
            i = abs(i)
            zsol_loci(k1) = Zsol_globb(i)
         else
            zsol_loci(k1) = Zsol_globi(i)
         endif
      enddo
      deallocate(A_MACRO(ielc)%GLOC(iel)%lcon)
!
      mdle = A_MACRO(ielc)%GLOC(iel)%mdle
      call solout_mg(ielc,iel,mdle,ni,1,zsol_loci)
      deallocate(zsol_loci)
   enddo
!
   call assembly_end_par
!
   deallocate(Zsol_globb)
!
!
   end subroutine solout_macro

