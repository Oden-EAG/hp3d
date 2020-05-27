!----------------------------------------------------------------------
!> Purpose  - routine uses the "minimum rule" to update order of
!!            approximation of an initial mesh edge-node or face-node.
!!
!> @param[in]    Etype   - element type
!> @param[in]    Nflag   = 2 mid-edge node
!> @param[in]            = 3 mid-face node
!> @param[in]    Ient    - edge or face number (entities)
!> @param[in]    Norient - face orientation (for faces only)
!> @param[in]    Norder  - element (middle node) order
!> @param[inout] Nord    - order of the node
!
!> @date Dec 14
!----------------------------------------------------------------------
!
subroutine min_order(Etype,Nflag,Ient,Norient,Norder, Nord)
!
      implicit none
      character(len=4), intent(in)    :: Etype
      integer         , intent(in)    :: Nflag,Ient,Norient,Norder
      integer         , intent(inout) :: Nord
!
      integer :: nordx,nordy,nordz, nordh,nordv
!
!----------------------------------------------------------------------
!
      select case(Etype)
!
!
!     -- PRISM --
      case('pris')
         call decode(Norder, nordx,nordz)
!
!        -- EDGES --
         select case(Nflag)
         case(2)
            select case(Ient)
!
!           horizontal
            case(1,2,3,4,5,6) ; Nord = min(Nord,nordx)
!
!           vertical
            case(7,8,9)       ; Nord = min(Nord,nordz)
            endselect
!
!        -- FACES --
         case(3)
            select case(Ient)
!
!           horizontal (triangle)
            case(1,2) ; Nord = min(Nord,nordx)
!
!           vertical (quad)
            case(3,4,5)
!
!              decode face ...
               call decode(Nord, nordh,nordv)
               select case(Norient)
!
!              non-flipping orientation
               case(0,2,5,7) ; nordh=min(nordh,nordx) ; nordv=min(nordv,nordz)
!
!              flipping orientation
               case(1,3,4,6) ; nordh=min(nordh,nordz) ; nordv=min(nordv,nordx)
               endselect
!
!              ... encode face
               Nord = nordh*10+nordv
            endselect
         endselect
!
!
!     -- BRICK --
      case('bric')
!
!        decode brick
         call decode(Norder, nordh,nordz)
         call decode(nordh, nordx,nordy)
!
!        -- EDGES --
         select case(Nflag)
         case(2)
            select case(Ient)
!
!           parallel to x-axis
            case(1,3,5,7)    ; Nord = min(Nord,nordx)
!
!           parallel to y-axis
            case(2,4,6,8)    ; Nord = min(Nord,nordy)
!
!           parallel to z-axis
            case(9,10,11,12) ; Nord = min(Nord,nordz)
            endselect
!
!        -- FACES --
         case(3)
!
!           decode face ...
            call decode(Nord, nordh,nordv)
!
            select case(Ient)
!
!           parallel to xy-plane
            case(1,2)
               select case(Norient)
!
!              non-flipping orientation
               case(0,2,5,7) ; nordh=min(nordh,nordx) ; nordv=min(nordv,nordy)
!
!              flipping orientation
               case(1,3,4,6) ; nordh=min(nordh,nordy) ; nordv=min(nordv,nordx)
               endselect
!
!
!           parallel to xz-plane
            case(3,5)
               select case(Norient)
!
!              non-flipping orientation
               case(0,2,5,7) ; nordh=min(nordh,nordx) ; nordv=min(nordv,nordz)
!
!              flipping orientation
               case(1,3,4,6) ; nordh=min(nordh,nordz) ; nordv=min(nordv,nordx)
               endselect
!
!           parallel to yz-plane
            case(4,6)
               select case(Norient)
!
!              non-flipping orientation
               case(0,2,5,7) ; nordh=min(nordh,nordy) ; nordv=min(nordv,nordz)
!
!              flipping orientation
               case(1,3,4,6) ; nordh=min(nordh,nordz) ; nordv=min(nordv,nordy)
               endselect
            endselect
!
!           ... encode face
            Nord = nordh*10+nordv
         endselect
!
!
!     -- TET --
      case('tetr')
         Nord = min(Nord,Norder)
!
!
!     -- PYRAMID --
      case('pyra')
!
!        -- EDGES --
         select case(Nflag)
         case(2)
            Nord = min(Nord,Norder)
!
!        -- FACES --
         case(3)
            select case(Ient)
!
!           bottom (quad)
            case(1)
!
!              decode face ...
               call decode(Nord, nordh,nordv)
!
!              face has ISOTROPIC order
               nordh=min(nordh,Norder) ; nordv=min(nordv,Norder)
!
!              ... encode face
               Nord = nordh*10+nordv
!
!           lateral(triangle)
            case(2,3,4,5)
               Nord = min(Nord,Norder)
            endselect
         endselect
!
!
      case default
         write(*,*) 'min_order: Etype = ',Etype
         stop
      endselect
!
!
endsubroutine min_order
