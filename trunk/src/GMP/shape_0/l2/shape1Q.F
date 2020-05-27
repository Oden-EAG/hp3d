c----------------------------------------------------------------------
c
c   routine name       - shape1Q_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 2012
c
c   purpose            - evaluate L2 1D shape functions
c
c   arguments :
c     in:
c        T             - local edge coordinate [0,1]
c        Norder        - order of approximation for the nodes (in the H1 sense)
c     out:
c        Nrdof         - number of shape functions
c        Shap1Q         - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1Q_0(T,Nord,Nrdof,Shap1Q)
c
      use control, only: SHAPE_FLAG
      use parameters

      implicit none
c
c  ... subroutine arguments
      double precision, intent(in) :: T
      integer,intent(in) :: Nord ! Order p of the element understood in the H1 sense
      integer, intent(out) :: Nrdof ! Number of L2 dof
      double precision,dimension(MAXP),intent(inout) :: Shap1Q


c     ...local variables
      double precision :: x
      double precision,dimension(MAXP-1) :: fake_dshapH ! needed to call H1 routines
      integer :: nordL2         ! Polynomial order of L2 space
      integer :: iprint, k
c
      iprint=0
c$$$      iprint_shapke=iprint
c
      nordL2 = Nord-1           ! WE ASSUME DUMMY ARGUMENT Nord REFERS TO H1 ORDER

      if (nordL2.lt.0) then
        Nrdof=0; Shap1Q=0.d0
        return
      endif
      if (nordL2.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format('shape1Q: nordL2,MAXP = ',2i2)
        stop
      endif
c

      select case(SHAPE_FLAG)
c
c  ...PEANO SHAPE FUNCTION.............................................
      case(1)
        write(*,*)'PEANO SHAPE FUNCTIONS, REALLY?'
        call pause
c     ..we simly call an existing routine for 1D
c       that have PEANO basis functions coded
        call shapHke(T,nordL2+2,0, Nrdof,Shap1Q,
     .               fake_dshapH)


c  ...LEGENDRE SHAPE FUNCTION.............................................
      case(2)

        x = 2.d0*T-1.d0         ! change to x in [-1,1]

c     ..LEGENDRE POLYNOMIALS....(NOT INTEGRATED LEGENDRE POLINOMIALS!!)....
c     .. Note that "shapHke" routine computes INTEGRATED legendre polynomials
        call shap1Q_legendre(T,nordL2,Nrdof,Shap1Q)

c
      case default
        write(*,*)'shape1Q: INVALID SHAPE_FLAG = ',SHAPE_FLAG
        stop
      end select
c


      if (iprint.eq.1) then
        write(*,7002) T,nordL2
 7002   format('shape1Q: T,nordL2 = ',f8.3,2x,i2,' Shap1Q = ')
        do k=1,Nrdof
          write(*,7003) k,Shap1Q(k)
 7003     format('         k = ',i2,f8.3)
        enddo
        call pause
      endif
c
      end subroutine shape1Q_0
