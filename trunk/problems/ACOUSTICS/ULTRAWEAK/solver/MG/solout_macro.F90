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
   subroutine solout_macro(Ielc,Zsoli,Ndof)
! 
   use assembly,         ONLY: assembly_begin_par, assembly_end_par
   use parameters,       ONLY: ZERO, ZONE
   use macro_grid_info,  ONLY: A_MACRO
!  
   IMPLICIT NONE
!
!-----------------------------------------------------------------------

#if C_MODE
   complex*16              :: Zsoli(Ndof)
   complex*16, allocatable :: Zsolb(:), Zsol_glob(:)
#else
   real*8                  :: Zsoli(Ndof)
   real*8,     allocatable :: Zsolb(:), Zsol_glob(:)
#endif
   integer                 :: jb, ji, mdle, i, iel, k1, ielc
   integer                 :: Ndof, ndof_tot, maxdofm
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
   jb = A_MACRO(ielc)%nb
   ji = A_MACRO(ielc)%ni
!
   allocate(Zsolb(jb))
!
#if C_MODE
   call ZGEMV('N',jb,ji,ZONE,A_MACRO(ielc)%array,jb,Zsoli,1,ZERO,Zsolb,1)
#else   
   call DGEMV('N',jb,ji,ZONE,A_MACRO(ielc)%array,jb,Zsoli,1,ZERO,Zsolb,1)
#endif

   Zsolb = A_MACRO(ielc)%vect - Zsolb
!
   deallocate(A_MACRO(ielc)%array)
   deallocate(A_MACRO(ielc)%vect)
! 
!..allocation of arrays required from solout_macro  
   call assembly_begin_par
!
   do iel = 1, A_MACRO(ielc)%nrmdle
      ndof_tot = A_MACRO(ielc)%GLOC(iel)%ndof
      allocate(Zsol_glob(ndof_tot))
!  
      do k1=1,ndof_tot
!     ...global dof is:
         i = A_MACRO(ielc)%GLOC(iel)%lcon(k1)
         if (i .lt. 0) then
            i = abs(i)
            Zsol_glob(k1) = Zsolb(i)
         else    
            Zsol_glob(k1) = Zsoli(i)
         endif    
      enddo
      deallocate(A_MACRO(ielc)%GLOC(iel)%lcon)
!
!  ...call solout(mdle)
      mdle = A_MACRO(ielc)%GLOC(iel)%mdle
      call solout1(mdle,Ndof_tot,1,1,Zsol_glob)
      deallocate(Zsol_glob)
   enddo
!
   call assembly_end_par
!
   deallocate(Zsolb)
!

   end subroutine solout_macro

