!----------------------------------------------------------------------
!                                                                     
!     routine name      - set_3Dint_fi
!                                                                     
!----------------------------------------------------------------------
!                                                                     
!     latest revision:  - 10/27/2016
!                                                                     
!     purpose:          - routine sets up quadrature data for a 3D
!                         element, accouting for different element
!                         types and orders of approximation
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             Type      - element type - IT ONLY SUPPORTS BRICK ELEMENTS
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
subroutine set_3Dint_fi(Type,Norder,nordx,nordy,nordz,nintx,ninty, &
                                       nintz,Xiloc,Waloc)
!
   use parametersDPG    , only : MAXNINT3ADD, MAXPP
   use control          , only : INTEGRATION
   use gauss_quadrature , only : INITIALIZED, XIGAUS1, WAGAUS1
!      
!
   implicit none
!
   integer :: iprint1,i,l,l1,l2,l3,nintx,ninty,nintz, &
              nordx,nordy,nordz,nordh,nordv,nord1,nord2,nord3
   character(len=4),                         intent(in)  :: Type
   integer,         dimension(19),           intent(in)  :: Norder
   real*8,          dimension(3,MAXNINT3ADD),intent(out) :: Xiloc,Waloc
!
!----------------------------------------------------------------------
!
!  ...initialize if needed
   if (.NOT. INITIALIZED) call init_gauss_quadrature
!      
   iprint1=0
   if (iprint1.eq.1) then
      write(*,7001) Type,Norder
 7001 format('set_3Dint_fi: Type, Norder = ',a4,2x,19i4)
   endif
!
   select case(Type)
!
!======================================================================      
!  BRICK                                                              |
!======================================================================      
   case('mdlb','bric')
!
!  ...determine order of approximation
      ! nordx=0 ; nordy=0 ; nordz=0
      ! do i=1,19
      !    select case(i)
      !    case(1,3,5,7)
      !       nordx = max(nordx,Norder(i))
      !    case(2,4,6,8)
      !       nordy = max(nordy,Norder(i))
      !    case(9,10,11,12)
      !       nordz = max(nordz,Norder(i))
      !    case(13,14)
      !       call decode(Norder(i), nordh,nordv)
      !       nordx = max(nordx,nordh)
      !       nordy = max(nordy,nordv)
      !    case(15,17)
      !       call decode(Norder(i), nordh,nordv)
      !       nordx = max(nordx,nordh)
      !       nordz = max(nordz,nordv)
      !    case(16,18)
      !       call decode(Norder(i), nordh,nordv)
      !       nordy = max(nordy,nordh)
      !       nordz = max(nordz,nordv)
      !    case(19)
            call decode(Norder(19), nordh,nord3)
            call decode(nordh, nord1,nord2)
            ! nordx = max(nordx,nord1)
            ! nordy = max(nordy,nord2)
            ! nordz = max(nordz,nord3)
            nordx = nord1
            nordy = nord2
            nordz = nord3
         ! end select
      ! enddo
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
      write(*,*) 'set_3Dint_fi:Type = ',Type
      stop
   endselect
!
!
end subroutine set_3Dint_fi
