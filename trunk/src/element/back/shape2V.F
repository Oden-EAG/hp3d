c--------------------------------------------------------------------
c
c     routine name      - shape2V_1
c
c--------------------------------------------------------------------
c
c     latest revision:  - Feb 14
c
c     purpose:          - routine returns values of a 2D element
c                         H(div) shape functions and their divergence
c
c     arguments:
c
c     in:
c             Type      - element type
c             Xi        - master coordinates
c             Norder    - polynomial order
c             Norient   - edge orientation
c
c     out:
c             NrdofV    - number of dof
c             ShapV     - values of the shape functions at the point
c             DivV      - divergence of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2V_1(Type,Xi,Norder,Norient, NrdofV,ShapV,DivV)
c
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Xi(2),Norder(5),Norient(4),
     .          ShapV(2,MAXtriaV),DivV(MAXtriaV)
c
c  ...H(div) shape functions are obtained by "rotating" H(curl) shape
c     functions
      call shape2E(Type,Xi,Norder,Norient, NrdofV,ShapV,VdivV)
c
c  ...'rotate'
      do k=1,NrdofV
        val = ShapV(1,k)
        ShapV(1,k) =  ShapV(2,k)
        ShapV(2,k) = -val
      enddo
c
      end subroutine shape2V_1
