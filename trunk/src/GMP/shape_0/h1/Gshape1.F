c----------------------------------------------------------------------
c
c   routine name       - Gshape1
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine evaluates 1D hierarchical C2-conforming
c                        shape functions, the first six are determined
c                        by boundary conditions
c
c  arguments
c              Nord    - order of approximation
c              Xi      - master element coordinate
c     out:
c              Vshap   - values of the shape functions
c              Dvshap  - first derivatives
c              Ddvshap - second derivatives
c
c----------------------------------------------------------------------
c
      subroutine Gshape1(Nord,Xi, Vshap,Dvshap,Ddvshap)
c
      use parameters
c
      implicit none
c
      integer,intent(in)                             :: Nord
      double precision,intent(in)                    :: Xi
      double precision,dimension(MAXP+1),intent(out) :: Vshap, Dvshap,
     .                                                  Ddvshap
c
c  ...coefficients defining the shape functions in terms of monomials
      double precision,dimension(6,6),save :: ashap
c
c  ...visitation flag
      integer,save :: nflag_Gshape1 = 0
c
c  ...auxiliary matrix and vector
      double precision,dimension(6,6) :: a
      double precision,dimension(6)   :: b
c
      integer :: j,k
c
c----------------------------------------------------------------------
c
      integer :: iprint
      iprint=0
c
c  ...check data
      if ((Nord.lt.5).or.(Nord.gt.5)) then
        write(*,7001) Nord
 7001   format('Gshape1: Nord = ',i3)
        stop 1
      endif
      if ((Xi.lt.0.d0).or.(Xi.gt.1.d0)) then
        write(*,7002)
 7002   format('Gshape1: WARNING! Xi = ',f8.3)
        call pause
      endif
c
c  ...if 1st visit determine coefficients
      if (nflag_Gshape1.eq.0) then
c  .....raise visitation flag
        nflag_Gshape1=1
        a(1:6,1:6) = 0.d0
        a(1,1) = 1.d0
        a(2,1:6) = 1.d0
        a(3,2) = 1.d0
        a(4,2)=1.d0; a(4,3)=2.d0; a(4,4)=3.d0; a(4,5)=4.d0; a(4,6)=5.d0
        a(5,3) = 2.d0
        a(6,3)=2.d0; a(6,4)=6.d0; a(6,5)=12.d0; a(6,6)=20.d0
        call tri(a,6,6)
        do j=1,6
          b(1:6) = 0.d0
          b(j) = 1.d0
          call rhsub(a,ashap(1:6,j),b,6,6)
        enddo
        if (iprint.eq.1) then
          do j=1,6
            write(*,7003) j
 7003       format('Gshape1: COEFFICIENTS FOR j= ',i1)
            write(*,7004) ashap(1:6,j)
 7004       format(6e12.5)
          enddo
          call pause
        endif
      endif
c
c  ...compute shape functions and their derivatives
      do k=1,6
        Vshap(k)=0.d0; Dvshap(k)=0.d0; Ddvshap(k) = 0.d0
        do j=1,6
          Vshap(k) = Vshap(k) + ashap(j,k)*Xi**(j-1)
        enddo
        do j=2,6
          Dvshap(k) = Dvshap(k) + ashap(j,k)*(j-1)*Xi**(j-2)
        enddo
        do j=3,6
          Ddvshap(k) = Ddvshap(k)
     .               + ashap(j,k)*(j-1)*(j-2)*Xi**(j-3)
        enddo
      enddo
c
c
      end subroutine Gshape1
c
c
c----------------------------------------------------------------------
c
c   routine name       - Gshap1
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine evaluates 1D hierarchical C2-conforming
c                        shape functions, the first six are determined
c                        by boundary conditions
c
c  arguments
c              Nord    - order of approximation
c              Xi      - master element coordinate
c     out:
c              Vshap   - values of the shape functions
c              Dvshap  - first derivatives
c              Ddvshap - second derivatives
c              D3vshap - third derivatives
c
c----------------------------------------------------------------------
c
      subroutine Gshap1(Nord,Xi, Vshap,Dvshap,Ddvshap,D3vshap)
c
      use parameters
c
      implicit none
c
      integer,intent(in)                             :: Nord
      double precision,intent(in)                    :: Xi
      double precision,dimension(MAXP+1),intent(out) :: Vshap,
     .                                                  Dvshap,
     .                                                  Ddvshap,
     .                                                  D3vshap
c
c  ...coefficients defining the shape functions in terms of monomials
      double precision,dimension(6,6),save :: ashap
c
c  ...auxiliary matrices
      double precision,dimension(6,6) :: a
      double precision,dimension(6)   :: b
c
c  ...visitation flag
      integer,save :: nflag_Gshap1 = 0
c
      integer :: j,k
c
c----------------------------------------------------------------------
c
      integer :: iprint
      iprint=0
c
c  ...check data
      if ((Nord.lt.5).or.(Nord.gt.5)) then
        write(*,7001) Nord
 7001   format('Gshap1: Nord = ',i3)
        stop 1
      endif
      if ((Xi.lt.0.d0).or.(Xi.gt.1.d0)) then
        write(*,7002)
 7002   format('Gshap1: WARNING! Xi = ',f8.3)
        call pause
      endif
c
c  ...determine coefficients
      if (nflag_Gshap1.eq.0) then
        nflag_Gshap1=1
        a(1:6,1:6) = 0.d0
        a(1,1) = 1.d0
        a(2,1:6) = 1.d0
        a(3,2) = 1.d0
        a(4,2)=1.d0; a(4,3)=2.d0; a(4,4)=3.d0; a(4,5)=4.d0; a(4,6)=5.d0
        a(5,3) = 2.d0
        a(6,3)=2.d0; a(6,4)=6.d0; a(6,5)=12.d0; a(6,6)=20.d0
        call tri(a,6,6)
        do j=1,6
          b(1:6) = 0.d0
          b(j) = 1.d0
          call rhsub(a,ashap(1:6,j),b,6,6)
        enddo
        if (iprint.eq.1) then
          do j=1,6
            write(*,7003) j
 7003       format('Gshap1: COEFFICIENTS FOR j= ',i1)
            write(*,7004) ashap(1:6,j)
 7004       format(6e12.5)
          enddo
          call pause
        endif
      endif
c
c  ...compute shape functions and their derivatives
      do k=1,6
        Vshap(k)=0.d0; Dvshap(k)=0.d0;
        Ddvshap(k) = 0.d0; D3vshap(k) = 0.d0
        do j=1,6
          Vshap(k) = Vshap(k) + ashap(j,k)*Xi**(j-1)
        enddo
        do j=2,6
          Dvshap(k) = Dvshap(k) + ashap(j,k)*(j-1)*Xi**(j-2)
        enddo
        do j=3,6
          Ddvshap(k) = Ddvshap(k)
     .               + ashap(j,k)*(j-1)*(j-2)*Xi**(j-3)
        enddo
        do j=4,6
          D3vshap(k) = D3vshap(k)
     .               + ashap(j,k)*(j-1)*(j-2)*(j-3)*Xi**(j-4)
        enddo
      enddo
c
c
      end subroutine Gshap1
