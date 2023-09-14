!-------------------------------------------------------------------------------
!>@brief       Set isotropic polynomial order for all initial mesh elements
!>@param[in]   Order     - polynomial order: 1,...,MAXP
!>@param[out]  ElemOrder - list of "anisotropic" element orders:
!!                            TETR:  1*Order; PYRA:   1*Order;
!!                            PRIS: 11*Order; HEXA: 111*Order.
!>@date        Sep 2023
!-------------------------------------------------------------------------------
subroutine set_order(Order, ElemOrder)
!
   use data_structure3D, only: NRELIS
!
   implicit none
!
   integer, intent(in)  :: Order
   integer, intent(out) :: ElemOrder(NRELIS)
!
   integer :: iel
!
!$OMP PARALLEL DO
!..loop through initial mesh elements
   do iel=1,NRELIS
      call set_order_elem(iel,Order, ElemOrder(iel))
   enddo
!$OMP END PARALLEL DO
!
end subroutine set_order


!-------------------------------------------------------------------------------
!>@brief       Set isotropic polynomial order for an initial mesh element
!>@param[in]   Iel       - initial mesh element: 1,...,NRELIS
!>@param[in]   Order     - polynomial order: 1,...,MAXP
!>@param[out]  ElemOrder - "anisotropic" element order:
!!                            TETR:  1*Order; PYRA:   1*Order;
!!                            PRIS: 11*Order; BRIC: 111*Order.
!>@date        Sep 2023
!-------------------------------------------------------------------------------
subroutine set_order_elem(Iel,Order, ElemOrder)
!
   use data_structure3D, only: ELEMS,NRELIS
   use parameters      , only: MAXP
   use node_types      , only: TETR,PYRA,PRIS,BRIC
!
   implicit none
!
   integer, intent(in)  :: Iel,Order
   integer, intent(out) :: ElemOrder
!
   integer :: etype
!
!..verifying input arguments
   if (Iel.lt.1 .or. Iel.gt.NRELIS) then
      write(*,1000) 'Iel',Iel
      stop
   endif
   if (Order.lt.1 .or. Order.gt.MAXP) then
      write(*,1000) 'Order',Order
      write(*,1000) 'MAXP ',MAXP
      stop
   endif
   1000 format('set_order_elem: invalid input: ',A,' = ',I9)
!
!..set isotropic order of approximation
   etype = ELEMS(iel)%etype
   select case(etype)
      case(TETR); ElemOrder =   1*Order
      case(PYRA); ElemOrder =   1*Order
      case(PRIS); ElemOrder =  11*Order
      case(BRIC); ElemOrder = 111*Order
   end select
!
end subroutine set_order_elem
