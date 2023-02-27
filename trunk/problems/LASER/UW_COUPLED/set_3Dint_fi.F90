!----------------------------------------------------------------------
!
!     routine name      - set_3Dint_fi
!
!----------------------------------------------------------------------
!
!     latest revision:  - Apr 2019
!
!     purpose:          - routine sets up quadrature data for a 3D
!                         element, accounting for different element
!                         types and orders of approximation
!
!     arguments:
!
!     in:
!             EType     - element type - IT ONLY SUPPORTS BRICK ELEMENTS
!             Norder    - order of approximation
!
!     out:
!             intx      - number of integration points for parameter 1
!             inty      - number of integration points for parameter 2
!             intz      - number of integration points for parameter 3
!             Xiloc     - integration points
!             Waloc     - weights
!
!----------------------------------------------------------------------
!
subroutine set_3Dint_fi(EType,Norder, nordx,nordy,nordz,nintx,ninty, &
                                      nintz,Xiloc,Waloc)
!
   use parametersDPG    , only : MAXNINT3ADD, MAXPP
   use control          , only : INTEGRATION
   use gauss_quadrature , only : INITIALIZED, XIGAUS1, WAGAUS1
   use node_types
!
   implicit none
!
   integer, intent(in)  :: EType
   integer, intent(in)  :: Norder(19)
   real(8), intent(out) :: Xiloc(3,MAXNINT3ADD)
   real(8), intent(out) :: Waloc(3,MAXNINT3ADD)
   integer, intent(out) :: nordx,nordy,nordz,nintx,ninty,nintz
!
   integer :: i,l,l1,l2,l3,nordh,nordv,nord1,nord2,nord3
!
#if DEBUG_MODE
   integer :: iprint
#endif
!
!----------------------------------------------------------------------
!
!..initialize if needed
   if (.not. INITIALIZED) call init_gauss_quadrature
!
#if DEBUG_MODE
   iprint=0
   if (iprint.eq.1) then
      write(*,7001) S_Type(etype),Norder
 7001 format('set_3Dint_fi: Type, Norder = ',a4,2x,19i4)
   endif
#endif
!
   select case(EType)
   case(MDLB,BRIC)
!
!  ...determine order of approximation
      call decode(Norder(19), nordh,nord3)
      call decode(nordh     , nord1,nord2)
      nordx = nord1
      nordy = nord2
      nordz = nord3
!
!  ...account for overintegration
      nordx=min(nordx+INTEGRATION,MAXPP)
      nordy=min(nordy+INTEGRATION,MAXPP)
      nordz=min(nordz+INTEGRATION,MAXPP)
!
!  ...compute number of integration points
      nintx=nordx+1
      ninty=nordy+1
      nintz=nordz+1
!
!  ...compute integration points and weights
      l=0
      do l3=1,nintz
         do l2=1,ninty
            do l1=1,nintx
               l=l+1
               Xiloc(1,l)=XIGAUS1(l1,nintx)
               Xiloc(2,l)=XIGAUS1(l2,ninty)
               Xiloc(3,l)=XIGAUS1(l3,nintz)
               Waloc(1,l)=WAGAUS1(l1,nintx)
               Waloc(2,l)=WAGAUS1(l2,ninty)
               Waloc(3,l)=WAGAUS1(l3,nintz)
            enddo
         enddo
      enddo
!
   case default
      write(*,*) 'set_3Dint_fi: Type = ', S_Type(Etype)
      stop
   endselect
!
!
end subroutine set_3Dint_fi
