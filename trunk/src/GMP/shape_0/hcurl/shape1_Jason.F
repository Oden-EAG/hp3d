
c----------------------------------------------------------------------
c
c   routine name       - shape1
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 05
c
c   purpose            - routine evaluates 1D hierarchical shape
c                        functions with their derivatives
c
c   arguments :
c     in:
c              Xi      - coordinate between 0 and 1
c              Norder  - order of approximation
c     out:
c              Nrdof   - number of dof (shape functions)
c              Vshap   - values of the shape functions at the point
c              Dvshap  - their derivatives
c
c----------------------------------------------------------------------
c
      subroutine shape1_Jason(Xi,Norder,  Nrdof,Vshap,Dvshap)
c
      use control
      use parameters, only:  NPMAX => MAXP
c
#include "syscom.blk"
c
      dimension Vshap(*),Dvshap(*)
c
      dimension d(1:2,3:(NPMAX+2)),v(3:(NPMAX+2))
      save d,v
      save nflag
      data nflag /0/
c
      Nrdof = Norder+1
      Vshap(1) = 1.d0-Xi; Dvshap(1) = -1.d0
      Vshap(2) = Xi;      Dvshap(2) =  1.d0
      if (Norder==1) return
c
      select case(SHAPE_FLAG)
      case(1)
        write(*,*)'PEANO SHAPE FUNCTIONS, REALLY?'
        call pause
c
c       PEANO SHAPE FUNCTIONS
c
        Vshap(3) = Vshap(1)*Vshap(2); Dvshap(3) = Vshap(1)-Vshap(2)
        x = Xi-0.5d0
        do i=4,Nrdof
          Vshap(i) = x*Vshap(i-1); Dvshap(i) = x*Dvshap(i-1)+Vshap(i-1)
        enddo
c
      case(2)
c
c       INTEGRATED LEGENDRE POLYNOMIALS
c
c       on first visit, set coefficients for 3-term recurrence relation
        if (nflag==0) then
          nflag = 1
          do i=3,NPMAX+2
            d(1,i) = dble(2*i-5)/dble(i-2)
            d(2,i) = dble(i-3)/dble(i-2)
            v(i)   = dble(1)/dble(2*(i-1))
          enddo
        endif
c
        x = 2.d0*Xi-1.d0
        do i=3,Nrdof
          Dvshap(i) = d(1,i)*x*Dvshap(i-1) - d(2,i)*Dvshap(i-2)
          Vshap(i)  = v(i)*(x*Dvshap(i)-Dvshap(i-1))
        enddo
c
      case default
        write(*,*)'shape1: INVALID SHAPE_FLAG = ',SHAPE_FLAG
        stop
      end select
c
      end subroutine shape1_Jason


      



