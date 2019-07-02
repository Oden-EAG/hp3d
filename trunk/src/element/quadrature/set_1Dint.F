c--------------------------------------------------------------------
c                                                                     
c     routine name      - set_1Dint
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Feb 11
c                                                                     
c     purpose:          - routine sets up quadrature data for a 2D
c                         element, accouting for different element
c                         types and orders of approximation
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Norder    - order of approximation
c        
c     out:    
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c             
c-----------------------------------------------------------------------
c
      subroutine set_1Dint(Nord, Nint,Xiloc,Waloc)
c
      use parameters       , only : MAXP
      use control          , only : INTEGRATION
      use gauss_quadrature , only : INITIALIZED,XIGAUS1,WAGAUS1
      implicit none 
c
      integer,                   intent(in)  :: Nord
      integer,                   intent(out) :: Nint
      real*8, dimension(MAXP+1), intent(out) :: Xiloc,Waloc
c
      integer :: l,nord_aux      
c
c  ...initialized integration constants, if needed      
      if (.NOT. INITIALIZED) call init_gauss_quadrature
c
c  ...collect integration points and weights      
      nord_aux = Nord + INTEGRATION
      nord_aux = min(nord_aux,MAXP)
      Nint = nord_aux + 1
      do l=1,Nint
        Xiloc(l) = XIGAUS1(l,Nint)
        Waloc(l) = WAGAUS1(l,Nint)
      enddo
c
c      
      end subroutine
