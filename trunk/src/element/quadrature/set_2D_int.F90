!----------------------------------------------------------------------
!
!     routine name      - set_2Dint
!
!     DEPRECATED, kept for backward compatibility (Oct 2019)
!                 use set_2D_int instead
!----------------------------------------------------------------------
subroutine set_2Dint(Ntype,Norder, Nint,Xiloc,Waloc)
!
      use parameters, only : MAX_NINT2, MAXP
      implicit none
!
      integer                        , intent(in)  :: Ntype
      integer, dimension(5)          , intent(in)  :: Norder
      integer                        , intent(out) :: Nint
      real(8), dimension(2,MAX_NINT2), intent(out) :: Xiloc
      real(8), dimension(  MAX_NINT2), intent(out) :: Waloc
!
      call set_2Dint_aux(Ntype,Norder,MAXP,MAX_NINT2, Nint,Xiloc,Waloc)
!
end subroutine set_2Dint
!
!----------------------------------------------------------------------
!     routine name: set_2D_int
!
!     latest rev  : Feb 2023
!
!     purpose     : routine sets up quadrature data for a standard
!                   2D element, accounting for different element
!                   types and orders of approximation in the element
!                   system of coordinates
!
!     arguments:
!
!     in:
!          Ntype        - element type
!          Norder       - order of approximation
!          Norient_face - face orientation
!
!     out:
!          Nint         - number of integration points
!          Xiloc        - integration points
!          Waloc        - weights
!
!----------------------------------------------------------------------
subroutine set_2D_int(Ntype,Norder,Norient_face, Nint,Xiloc,Waloc)
!
      use element_data, only : NFAXES
      use parameters  , only : MAX_NINT2
      use node_types
      implicit none
!
      integer                        , intent(in)  :: Ntype
      integer, dimension(5)          , intent(in)  :: Norder
      integer                        , intent(in)  :: Norient_face
      integer                        , intent(out) :: Nint
      real(8), dimension(2,MAX_NINT2), intent(out) :: Xiloc
      real(8), dimension(  MAX_NINT2), intent(out) :: Waloc
!
      integer, dimension(5) :: norder_loc
      integer               :: nordh,nordv
!
      norder_loc = Norder
      select case(Ntype)
      case(RECT,MDLQ,QUAD)
        if (NFAXES(3,Norient_face).eq.1) then
          call decode(Norder(5), nordh,nordv)
          norder_loc(5) = nordv*10+nordh
        endif
      end select
!
      call set_2Dint(Ntype,norder_loc, Nint,Xiloc,Waloc)
!
end subroutine set_2D_int
!
!----------------------------------------------------------------------
!
!     routine name      - set_2Dint_DPG
!
!     DEPRECATED, kept for backward compatibility (Oct 2019)
!                 use set_2D_int_DPG instead
!----------------------------------------------------------------------
subroutine set_2Dint_DPG(Ntype,Norder, Nint,Xiloc,Waloc)
!
      use parametersDPG, only : MAXNINT2ADD, MAXPP
      implicit none
!
      integer                          , intent(in)  :: Ntype
      integer, dimension(5)            , intent(in)  :: Norder
      integer                          , intent(out) :: Nint
      real(8), dimension(2,MAXNINT2ADD), intent(out) :: Xiloc
      real(8), dimension(  MAXNINT2ADD), intent(out) :: Waloc
!
      call set_2Dint_aux(Ntype,Norder,MAXPP,MAXNINT2ADD, Nint,Xiloc,Waloc)
!
end subroutine set_2Dint_DPG
!
!----------------------------------------------------------------------
!     routine name: set_2D_int_DPG
!
!     latest rev  : Feb 2023
!
!     purpose     : routine sets up quadrature data for a DPG
!                   2D element, accounting for different element
!                   types and orders of approximation in the element
!                   system of coordinates
!
!     arguments:
!
!     in:
!          Ntype        - element type
!          Norder       - order of approximation
!          Norient_face - face orientation
!
!     out:
!          Nint         - number of integration points
!          Xiloc        - integration points
!          Waloc        - weights
!
!----------------------------------------------------------------------
subroutine set_2D_int_DPG(Ntype,Norder,Norient_face, Nint,Xiloc,Waloc)
!
      use element_data , only : NFAXES
      use parametersDPG, only : MAXNINT2ADD
      use node_types
      implicit none
!
      integer                          , intent(in)  :: Ntype
      integer, dimension(5)            , intent(in)  :: Norder
      integer                          , intent(in)  :: Norient_face
      integer                          , intent(out) :: Nint
      real(8), dimension(2,MAXNINT2ADD), intent(out) :: Xiloc
      real(8), dimension(  MAXNINT2ADD), intent(out) :: Waloc
!
      integer, dimension(5) :: norder_loc
      integer               :: nordh,nordv
!
      norder_loc = Norder
      select case(Ntype)
      case(RECT,MDLQ,QUAD)
        if (NFAXES(3,Norient_face).eq.1) then
          call decode(Norder(5), nordh,nordv)
          norder_loc(5) = nordv*10+nordh
        endif
      end select
!
      call set_2Dint_DPG(Ntype,norder_loc, Nint,Xiloc,Waloc)
!
end subroutine set_2D_int_DPG
!
!----------------------------------------------------------------------
!
!     routine name      - set_2Dint_aux
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine sets up quadrature data for a 2D
!                         element, accounting for different element
!                         types and orders of approximation
!
!     arguments:
!
!     in:
!             Ntype     - element type
!             Norder    - order of approximation
!             Maxp      - maximum p
!             Max_int2  - maximum number of integration points
!
!     out:
!             Nint      - number of integration points
!             Xiloc     - integration points
!             Waloc     - weights
!
!----------------------------------------------------------------------
subroutine set_2Dint_aux(Ntype,Norder,Maxp,Max_nint2, Nint,Xiloc,Waloc)
!
      use parameters       , only : MODORDER
      use control          , only : INTEGRATION
      use gauss_quadrature , only : INITIALIZED,NSELECT,NRGAUPO,    &
                                    XIGAUSS,WAGAUSS,XIGAUS1,WAGAUS1
      use node_types
      implicit none
!
      integer                        , intent(in)  :: Ntype
      integer, dimension(5)          , intent(in)  :: Norder
      integer                        , intent(in)  :: Maxp,Max_nint2
      integer                        , intent(out) :: Nint
      real(8), dimension(2,Max_nint2), intent(out) :: Xiloc
      real(8), dimension(  Max_nint2), intent(out) :: Waloc
!
      integer :: nordxy(2)
      integer :: kint,l,l1,l2,nord,nordx,nordy,nintx,ninty
!
!  ...initialize if needed
      if (.not. INITIALIZED) call init_gauss_quadrature
!
      select case(Ntype)
!
!  ...triangle
      case(MDLT,TRIA)
        nord = max(Norder(1),Norder(2),Norder(3),Norder(4))
!
!  .....set the quadrature
        nord = nord + INTEGRATION
!
        nord = min(nord,Maxp)
        kint = NSELECT(nord)
        Nint = NRGAUPO(kint)
        do l=1,Nint
          Xiloc(1:2,l) = XIGAUSS(2:3,l,kint)
          Waloc(l) = WAGAUSS(l,kint)/2.d0
        enddo
!
      nordx = 0; nordy = 0
!  ...quad
      case(MDLQ,QUAD,RECT)
        call decod(Norder(5),MODORDER,2, nordxy)
        nordx = max(Norder(1),Norder(3),nordxy(1))
        nordy = max(Norder(2),Norder(4),nordxy(2))
!
!  .....set the quadrature
        nordx=nordx+INTEGRATION
        nordy=nordy+INTEGRATION
!
!  .....limit nordx,y to Maxp
        nordx=min(nordx,Maxp) ; nordy=min(nordy,Maxp)
        nintx=nordx+1         ; ninty=nordy+1
!
        Nint=nintx*ninty
        l=0
        do l2=1,ninty
          do l1=1,nintx
            l=l+1
            Xiloc(1,l) = XIGAUS1(l1,nintx)
            Xiloc(2,l) = XIGAUS1(l2,ninty)
            Waloc(  l) = WAGAUS1(l1,nintx)*WAGAUS1(l2,ninty)
          enddo
        enddo
!
      case default
        write(*,*) 'set_2Dint: Type = ',S_Type(Ntype)
        stop
      end select
!
end subroutine set_2Dint_aux
