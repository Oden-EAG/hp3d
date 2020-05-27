! -----------------------------------------------------------------------
!
!    routine name       - local_schur_compl_extension
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
#include "implicit_none.h"
! 
   subroutine schur_compl_extension(Igrid,z,n)

   use mg_data_structure,only: GRID
   use assembly_sc,      only: IPRINT_TIME
!
   implicit none

   integer, intent(in) :: Igrid,n
   VTYPE,   intent(in) :: z(n)
!
!..locals
   integer :: iel, ndof_m,ndof_c,i,k, iprint_temp
!-----------------------------------------------------------------------
!
   iprint_temp = IPRINT_TIME
   IPRINT_TIME = 0
!
   do iel=1,GRID(Igrid)%nreles
      ndof_m = GRID(Igrid)%dloc(iel)%ndof ! macro dof
      do k=1,ndof_m
         i = GRID(Igrid)%dloc(iel)%lcon(k)
         GRID(Igrid)%loc(iel)%z(k) = z(i)
      enddo
!  ...compute Schur complement extension
      call local_schur_compl_extension(Igrid,iel,ndof_m,GRID(Igrid)%loc(iel)%z)
   enddo   
! 
!..prolong the correction to the finest level
   do iel=1,GRID(Igrid+1)%nreles
      ndof_m = GRID(Igrid+1)%dloc(iel)%ndof ! macro dof
      ndof_c = GRID(Igrid+1)%dloc(iel)%ndof_c  ! coarse dof
!  ...apply the prolongation operator
      call coarse2macro(Igrid+1,iel,GRID(Igrid+1)%loc(iel)%z_c, &
                        ndof_c,GRID(Igrid+1)%loc(iel)%z,ndof_m)
!      
   enddo

   IPRINT_TIME = iprint_temp
!
   end subroutine schur_compl_extension
!
! -----------------------------------------------------------------------
!
!    routine name       - local_schur_compl_extension
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
! 
   subroutine local_schur_compl_extension(Igrid,Iel,n1,z1)
! 
   use mg_data_structure,  only: GRID
   use parameters,         only: ZERO, ZONE
!  
   implicit none
!
!-----------------------------------------------------------------------

   integer, intent(in)  :: Igrid, Iel, n1
   VTYPE,   intent(in)  :: z1(n1)
!
!..locals
   integer              :: i, inz, k, ielm, iel_index, n2
   integer              :: ndof_c, ndof_m
   integer, allocatable :: IA(:), JA(:)
   VTYPE,   allocatable :: Zsol(:), z2(:), SA(:)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 

   n2 = GRID(igrid)%sch(Iel)%n2
!
!..check if there is anything to put in. Otherwise return
   if (n2 .eq. 0) then 

      do ielm = 1, GRID(Igrid)%sch(Iel)%nrmdle
!         
!     ...find the correct mdle number according to the natural order
         iel_index = GRID(Igrid)%sch(Iel)%GLOC(ielm)%iel


         ndof_c = GRID(Igrid)%sch(Iel)%gloc(ielm)%ndof_c
         allocate(Zsol(ndof_c))
!
         do k=1,ndof_c
!        ...global dof is:
            i = GRID(Igrid)%sch(iel)%gloc(ielm)%lcon(k)
            if (i .lt. 0) then
               write(*,*) 'schur_complement_in: INCONSISTENCY'
            else    
               Zsol(k) = z1(i)
            endif   
         enddo
!
!     ...store the solution in the solution arrays
         GRID(Igrid+1)%loc(iel_index)%z_c = Zsol
         deallocate(Zsol)
      enddo
!
      return
!      
   endif

   inz = GRID(Igrid)%sch(Iel)%inz


   allocate(z2(n2)) ; z2 = GRID(Igrid)%sch(Iel)%r2
   deallocate(GRID(Igrid)%sch(Iel)%r2)

   allocate(IA(n2+1))
   IA = GRID(Igrid)%sch(iel)%IA

   allocate(JA(inz)) 
   JA = GRID(Igrid)%sch(iel)%JA

   allocate(SA(inz)) 
   SA = GRID(Igrid)%sch(iel)%SA


   call pardiso_solve(IA,JA,SA,'H',inz,n2,1,z2)

!..compute -A22^-1*A21*z1
#if C_MODE
   call ZGEMV('N',n2,n1,-ZONE,GRID(Igrid)%sch(iel)%A21,n2,z1,1,ZERO,z2,1)
#else
   call DGEMV('N',n2,n1,-ZONE,GRID(Igrid)%sch(iel)%A21,n2,z1,1,ZERO,z2,1)
#endif   
!  
!..now restrict solutions z1 and z2 each element
!
   do ielm = 1, GRID(Igrid)%sch(Iel)%nrmdle
!      
!  ...find the correct mdle number according to the natural order
      iel_index = GRID(Igrid)%sch(iel)%gloc(ielm)%iel
      ndof_c = GRID(Igrid)%sch(iel)%gloc(ielm)%ndof_c
      allocate(Zsol(ndof_c))
!
      do k=1,ndof_c
!    ...global dof is:
         i = GRID(Igrid)%sch(iel)%gloc(ielm)%lcon(k)
         if (i .lt. 0) then
            i = abs(i)
            Zsol(k) = z2(i)
         else    
            Zsol(k) = z1(i)
         endif   
      enddo
!
!  ...store the solution in the solution arrays
      GRID(Igrid+1)%loc(iel_index)%z_c = Zsol
      deallocate(Zsol)
   enddo
!
   deallocate(z2, IA, JA, SA)
!
!
   end subroutine local_schur_compl_extension

