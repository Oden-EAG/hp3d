! -----------------------------------------------------------------------
!
!    routine name       - schur_compl_restriction
!
! -----------------------------------------------------------------------
!
!    latest revision    - MAY 18
!
!    purpose            - routine applies the Schur complement extension
!                         operator P to a macro-grid solution vector
!                              --             --
!                             |        I        |
!                         P = |                 |
!                             |  -A22^(-1)*A21  |
!                              --             --
!    arguments
!           in
!               Igrid   - grid index
!                 Iel   - element identifier corresponding to
!                         the natural ordering of elements
!                  n1   - size of the solution vector
!                  z1   - coarse grid solution vector
!          out
!               the extension is stored in the module mg_data_structure
!
!
! ----------------------------------------------------------------------
!
!
#include "typedefs.h"

   subroutine schur_compl_restriction(Igrid,Rsd,n)
!
   use mg_data_structure,only: GRID
   use parameters, only: ZERO
!
   implicit none
!
   integer, intent(in)  :: Igrid,n
   VTYPE,  intent(out)  :: Rsd(n)
!
!..locals
   integer :: iel, ndof_m, ndof_c, nrmdlem, n1, n2
   integer :: ielm, iel_index,k,ic
   VTYPE, allocatable :: r_loc_1(:), r_loc_2(:)

!
!---------------------------------------------------
!
!..restrict local residuals to the coarser grid
   do iel=1,GRID(Igrid+1)%nreles
      ndof_m = GRID(Igrid+1)%dloc(iel)%ndof
      ndof_c = GRID(Igrid+1)%dloc(iel)%ndof_c
!
      call macro2coarse(Igrid+1,iel,GRID(Igrid+1)%loc(iel)%r, &
                        ndof_m,GRID(Igrid+1)%loc(iel)%r_c,ndof_c)
   enddo
!
!..step 2: Schur complement restriction
   do iel=1,GRID(Igrid)%nreles
      nrmdlem = GRID(Igrid)%sch(iel)%nrmdle
      n1 = GRID(Igrid)%sch(iel)%n1
      n2 = GRID(Igrid)%sch(iel)%n2
      allocate(r_loc_1(n1)); r_loc_1 = ZERO
      allocate(r_loc_2(n2)); r_loc_2 = ZERO
      do ielm = 1, nrmdlem
         ndof_c = GRID(Igrid)%sch(iel)%gloc(ielm)%ndof_c
         iel_index = GRID(Igrid)%sch(iel)%gloc(ielm)%iel
!     ...calculate residuals corresponding to both bubble and skeleton dof
         do k = 1,ndof_c
            ic = GRID(Igrid)%sch(iel)%gloc(ielm)%lcon(k)
            if (ic .lt. 0 ) then
               ic = abs(ic)
               r_loc_2(ic) = r_loc_2(ic) + GRID(Igrid+1)%loc(iel_index)%r_c(k)
            else
               r_loc_1(ic) = r_loc_1(ic) + GRID(Igrid+1)%loc(iel_index)%r_c(k)
            endif
         enddo
      enddo
!
      call local_schur_compl_restriction(Igrid,iel,n1,n2,r_loc_1,r_loc_2)
!
!  ...store r_loc_1 for each macro element
      GRID(Igrid)%loc(iel)%r = r_loc_1

!  ...assemble global residual for the sub-master macro mesh
      ndof_m = GRID(Igrid)%dloc(iel)%ndof
      do k=1,ndof_m
!     ...global coarse dof is:
         ic = GRID(Igrid)%dloc(iel)%lcon(k)
!     ...assemble global load vector
         Rsd(ic) =  Rsd(ic) + r_loc_1(k)
         enddo
      deallocate(r_loc_1, r_loc_2)
   enddo




   end subroutine schur_compl_restriction


! -----------------------------------------------------------------------
!
!    routine name       - local_schur_compl_restriction
!
! -----------------------------------------------------------------------
!
!    latest revision    - MAY 18
!
!    purpose            - routine applies the Schur complement restriction
!                         operator P* to a fine-grid residual vector
!
!                               --                    --
!                         P* = |   I    -A22^(-1)*A21   |
!                               --                    --
!    arguments
!           in
!               Igrid   - grid index
!                 Iel   - element identifier corresponding to
!                         the natural ordering of elements
!                  n1   - size of the residual vector for the interface dof
!                  n2   - size of the residual vector for the interior dof
!                  r2   - residual vector for for the interior dof
!
!          inout
!                  r1   - residual vector for for the interface dof
!
!
! ----------------------------------------------------------------------
!
!

   subroutine local_schur_compl_restriction(Igrid,Iel,n1,n2,r1,r2)
!
   use mg_data_structure,  only: GRID
   use parameters,         only: ZERO, ZONE
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer,    intent(in) :: Igrid,Iel, n1, n2
   VTYPE,      intent(in) :: r2(n2)
   VTYPE,   intent(inout) :: r1(n1)
   VTYPE,     allocatable :: A12(:,:)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!..check if there is something to condense out. If not return
   if (n2 .eq. 0) return

!..save r2
   allocate(GRID(Igrid)%sch(Iel)%r2(n2)) ; GRID(Igrid)%sch(iel)%r2 = r2

   allocate(A12(n1,n2))
!..compute r1 = r1 -A12*A22^(-1)*r2
#if HP3D_COMPLEX
   A12 = transpose(conjg(GRID(Igrid)%sch(iel)%A21))
   call ZGEMV('N',n1,n2,-ZONE,A12,n1,r2,1,ZONE,r1,1)
#else
   A12 = transpose(GRID(Igrid)%sch(iel)%A21)
   call DGEMV('N',n1,n2,-ZONE,A12,n1,r2,1,ZONE,r1,1)
#endif
!
   deallocate(A12)
!
!
   end subroutine local_schur_compl_restriction


