!----------------------------------------------------------------------
!
!     routine name      - set_1D_int
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine sets up quadrature data for a 2D
!                         standard element, accounting for different
!                         element types and orders of approximation
!
!     arguments:
!
!     in:
!             Nord      - order of approximation
!
!     out:
!             Nint      - number of integration points
!             Xiloc     - integration points
!             Waloc     - weights
!
!----------------------------------------------------------------------
subroutine set_1Dint(Nord, Nint,Xiloc,Waloc)
!
      use parameters       , only : MAXP
      use control          , only : INTEGRATION
      use gauss_quadrature , only : INITIALIZED,XIGAUS1,WAGAUS1
      implicit none
!
      integer,                    intent(in)  :: Nord
      integer,                    intent(out) :: Nint
      real(8), dimension(MAXP+1), intent(out) :: Xiloc,Waloc
!
      integer :: l,nord_aux
!
!  ...initialized integration constants, if needed
      if (.not. INITIALIZED) call init_gauss_quadrature
!
!  ...collect integration points and weights
      nord_aux = Nord + INTEGRATION
      nord_aux = min(nord_aux,MAXP)
      Nint = nord_aux + 1
      do l=1,Nint
        Xiloc(l) = XIGAUS1(l,Nint)
        Waloc(l) = WAGAUS1(l,Nint)
      enddo
!
end subroutine set_1Dint
!
!
!----------------------------------------------------------------------
!
!     routine name      - set_1Dint_DPG
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine sets up quadrature data for a 1D
!                         DPG element, accounting for different
!                         element types and orders of approximation
!
!     arguments:
!
!     in:
!             Nord      - order of approximation
!
!     out:
!             Nint      - number of integration points
!             Xiloc     - integration points
!             Waloc     - weights
!
!----------------------------------------------------------------------
subroutine set_1Dint_DPG(Nord, Nint,Xiloc,Waloc)
!
      use parametersDPG    , only : MAXPP
      use control          , only : INTEGRATION
      use gauss_quadrature , only : INITIALIZED,XIGAUS1,WAGAUS1
      implicit none
!
      integer,                     intent(in)  :: Nord
      integer,                     intent(out) :: Nint
      real(8), dimension(MAXPP+1), intent(out) :: Xiloc,Waloc
!
      integer :: l,nord_aux
!
!  ...initialized integration constants, if needed
      if (.not. INITIALIZED) call init_gauss_quadrature
!
!  ...collect integration points and weights
      nord_aux = Nord + INTEGRATION
      nord_aux = min(nord_aux,MAXPP)
      Nint = nord_aux + 1
      do l=1,Nint
        Xiloc(l) = XIGAUS1(l,Nint)
        Waloc(l) = WAGAUS1(l,Nint)
      enddo
!
end subroutine set_1Dint_DPG
