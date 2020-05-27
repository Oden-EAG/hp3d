!------------------------------------------------------------------------
!
!    routine name      - inject_dof
!
!------------------------------------------------------------------------
!
!    latest revision   - Apr 2018
!
!    purpose           - routine establishes dof injections
!                        after a p-refinement
!
!
!   arguments :
!     in:
!                 Type - type of the nod
!             Nord_old - old order of the nod
!                 Nord - new order of the nod
!     out:
!                 InjH - Injection of the H1 dof
!                 InjE - Injection of the H(curl) dof
!                 InjV - Injection of the H(div) dof
!
!-----------------------------------------------------------------------
!
!
   subroutine inject_dof(Type,Nord_old,Nord, InjH, InjE, InjV)
!
   use data_structure3D, only: ndof_nod
   use parameters
!
   implicit none
!
   character(len=4), intent(in)  :: Type
   integer,          intent(in)  :: Nord_old, Nord
   integer,          intent(out) :: InjH(MAXquadH), InjE(MAXquadE), InjV(MAXquadV)
!
   integer :: i, j, kH, kE, kV, mE
   integer :: nordh, nordv, nordh_old, nordv_old
   integer :: ndofH_old, ndofE_old, ndofV_old, nvoid
!
!------------------------------------------------------------------------
!
   InjH = 0; InjE = 0; InjV = 0
!
   call ndof_nod(Type,Nord_old, ndofH_old,ndofE_old,ndofV_old,nvoid)
!
   select case(Type)
   case('vert','medg','mdlt')
      do j=1,ndofH_old
         InjH(j) = j
      enddo
      do j=1,ndofE_old
         InjE(j) = j
      enddo
      do j=1,ndofV_old
         InjV(j) = j
      enddo
   case('mdlq')
      call decode(Nord, nordh,nordv)
      call decode(Nord_old, nordh_old,nordv_old)
!
!  ...H1
      kH=0
      do j=1,nordv_old-1
         do i=1,nordh_old-1
            kH=kH+1
            InjH(kH) = (j-1)*(nordh-1)+i
         enddo
      enddo
!
!  ...H(curl)
!  ...NOTE: The ordering of loops for the two families are different
      kE = 0
!  ...first family
      do j=1,nordh_old-1
         do i=1,nordv_old
            kE=kE+1
            InjE(kE) = (j-1)*(nordv)+i
         enddo
      enddo
!
      mE = (nordh-1)*nordv
!
!  ...second family
      do j=1,nordh_old
         do i=1,nordv_old-1
            kE=kE+1
            InjE(kE) = mE + (j-1)*(nordv-1)+i
         enddo
      enddo
!
!  ...H(div)
      kV = 0
      do j=1,nordv_old
         do i=1,nordh_old
            kV=kV+1
            InjV(kV) = (j-1)*(nordh)+i
         enddo
      enddo
   end select


   end subroutine inject_dof
