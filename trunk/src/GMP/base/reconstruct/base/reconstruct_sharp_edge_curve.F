c----------------------------------------------------------------------
c
c   routine name       - reconstruct_sharp_edge_curve
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine reconstructs a curve
c
c   arguments
c     in:
c            Xpoint    - endpoint coordinates
c            Evect     - edpoint unit vectors
c     out:
c            Evect     - endpoint tangent vectors
c            Eprim     - endpoint second derivatives
c            Nerr_flag - error flag
c                        = 0 if no error
c                        = i(=1,2) if i-th tangent vector has been
c                                  reversed
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_sharp_edge_curve(Xpoint,Evect,Eprim,
     .                                        Nerr_flag)
c
      use control
c
      implicit none
c
      integer :: Nerr_flag
      real(8) :: Xpoint(3,2),Evect(3,2),Eprim(3,2)
c
c  ...matrices for a system of linear equations
      real(8) :: aa(2,2),bb(2),uu(2)
c
c  ...secant vector
      real(8) :: dx(3)
c
      integer :: i
      real(8) :: s,s1,s2,s3
c
      integer :: iprint
c**********************************************************************
c
      iprint=0
      Nerr_flag=0
c
c  ...check input
      do i=1,2
        call norm(Evect(1:3,i), s)
        if (abs(s-1.d0).gt.GEOM_TOL) then
          write(*,7001)  Evect(1:3,i)
 7001     format('reconstruct_curve: Evect = ',3f8.3)
          stop 1
        endif
      enddo
c
c  ...set up the system of 2x2 equations for the coefficients
c     (length of tangent vectors)
      aa(1,1) = 4.d0; aa(2,2) = 4.d0
      call scalar_product(Evect(1:3,1),Evect(1:3,2), s1)
      aa(1,2) = 2.d0*s1; aa(2,1) = aa(1,2)
c
      dx(1:3) = Xpoint(1:3,2)-Xpoint(1:3,1)
      call scalar_product(dx,Evect(1:3,1), s2)
      bb(1) = 6.d0*s2
      call scalar_product(dx,Evect(1:3,2), s3)
      bb(2) = 6.d0*s3
c
      call gausse(aa,2,bb, uu,2)
c
c  ...check consistency and redefine unit vectors
      do i=1,2
        if (uu(i).le.0.d0) then
          Nerr_flag = i
          write(*,7002) s1,s2,s3
 7002     format('             e_1 o e_2, e_1 o dr, e_2 o dr = ',
     .           3f8.3)
          stop 1
        else
          Evect(1:3,i) =  uu(i)*Evect(1:3,i)
        endif
      enddo
c
c  ...compute the second derivatives
      Eprim(1:3,1) = -6.d0*Xpoint(1:3,1) + 6.d0*Xpoint(1:3,2)
     .                - 4.d0*Evect(1:3,1) - 2.d0*Evect(1:3,2)
      Eprim(1:3,2) =  6.d0*Xpoint(1:3,1) - 6.d0*Xpoint(1:3,2)
     .                + 2.d0*Evect(1:3,1) + 4.d0*Evect(1:3,2)
c
c
      end subroutine reconstruct_sharp_edge_curve
