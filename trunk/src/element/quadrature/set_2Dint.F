c--------------------------------------------------------------------
c
c     routine name      - set_2Dint
c
c--------------------------------------------------------------------
c
c     latest revision:  - Sep 17
c
c     purpose:          - routine sets up quadrature data for a standard 
c                         2D element, accounting for different element
c                         types and orders of approximation
c
c     arguments:
c
c     in:
c             Type      - element type
c             Norder    - order of approximation
c
c     out:
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c
c-----------------------------------------------------------------------
c
      subroutine set_2Dint(Type,Norder,Nint,Xiloc,Waloc)
c
      use parameters, only : MAXquadH, MAXP
c
#include "syscom.blk"
      character(len=4) :: Type
      dimension Norder(5), Xiloc(2,MAXquadH),Waloc(MAXquadH)
C
      call set_2Dint_aux(Type,Norder,MAXP,MAXquadH,Nint,Xiloc,Waloc)       
C
      end subroutine set_2Dint
c      
c--------------------------------------------------------------------
c
c     routine name      - set_2Dint_DPG
c
c--------------------------------------------------------------------
c
c     latest revision:  - Sep 17
c
c     purpose:          - routine sets up quadrature data for a DPG 
c                         2D element, accounting for different element
c                         types and orders of approximation
c
c     arguments:
c
c     in:
c             Type      - element type
c             Norder    - order of approximation
c
c     out:
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c
c-----------------------------------------------------------------------
c
      subroutine set_2Dint_DPG(Type,Norder,Nint,Xiloc,Waloc)
c
      use parametersDPG, only : MAXquadHH, MAXPP
c
#include "syscom.blk"
      character(len=4) :: Type
      dimension Norder(5), Xiloc(2,MAXquadHH),Waloc(MAXquadHH)
C
      call set_2Dint_aux(Type,Norder,MAXPP,MAXquadHH,Nint,Xiloc,Waloc)       
C
      end subroutine set_2Dint_DPG      

c--------------------------------------------------------------------
c
c     routine name      - set_2Dint_aux
c
c--------------------------------------------------------------------
c
c     latest revision:  - Sep 17
c
c     purpose:          - routine sets up quadrature data for a 2D
c                         element, accounting for different element
c                         types and orders of approximation
c
c     arguments:
c
c     in:
c             Type      - element type
c             Norder    - order of approximation
c             Maxp      - maximum p
c             Max_int2  - maximum number of integration points
c
c     out:
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c
c-----------------------------------------------------------------------
c
      subroutine set_2Dint_aux(Type,Norder,Maxp,Max_nint2,
     .                         Nint,Xiloc,Waloc)
c
      use parameters       , only : MODORDER
      use control          , only : INTEGRATION
      use gauss_quadrature , only : INITIALIZED,NSELECT,NRGAUPO,
     .                              XIGAUSS,WAGAUSS,
     .                              XIGAUS1,WAGAUS1
c
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Norder(5),Xiloc(2,Max_nint2),Waloc(Max_nint2)
      integer :: nordxy(2)
c
c  ...initialize if needed
      if (.NOT. INITIALIZED) call init_gauss_quadrature
c
      select case(Type)
c
c  ...triangle
      case('mdlt','tria')
        nord = max(Norder(1),Norder(2),Norder(3),Norder(4))
c
c  .....set the quadrature
        nord = nord + INTEGRATION
c
c  .....limit nord to Maxp, so that dimension Max_nint2 = (Maxp + 1)^2 
c       is not exceeded
        nord = min(nord,Maxp)
        kint = NSELECT(nord)
        Nint = NRGAUPO(kint)
        do l=1,Nint
          Xiloc(1:2,l) = XIGAUSS(2:3,l,kint)
          Waloc(l) = WAGAUSS(l,kint)/2.d0
        enddo
c
      nordx = 0; nordy=0
c  ...quad
      case('mdlq','quad','rect')
        call decod(Norder(5),MODORDER,2, nordxy)
        nordx = max(Norder(1),Norder(3),nordxy(1))
        nordy = max(Norder(2),Norder(4),nordxy(2))
c
c  .....set the quadrature
        nordx=nordx+INTEGRATION
        nordy=nordy+INTEGRATION
c
c  .....limit nordx,y to Maxp
        nordx=min(nordx,Maxp) ; nordy=min(nordy,Maxp)
        nintx=nordx+1         ; ninty=nordy+1
c
        Nint=nintx*ninty
        l=0
        do l2=1,ninty; do l1=1,nintx
          l=l+1
          Xiloc(1,l) = XIGAUS1(l1,nintx)
          Xiloc(2,l) = XIGAUS1(l2,ninty)
          Waloc(  l) = WAGAUS1(l1,nintx)*WAGAUS1(l2,ninty)
        enddo ; enddo
c
      case default
        write(*,*) 'set_2Dint: Type = ',Type
        stop 1
      end select
c
c
      end subroutine set_2Dint_aux
