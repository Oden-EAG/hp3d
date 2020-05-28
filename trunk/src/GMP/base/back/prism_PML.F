c----------------------------------------------------------------------
c
c   routine name       - prism_PML
c
c---------------------------------------------------------------------
c
c
c   computer           - machine independent
c
c   latest revision    - Nov 08
c
c   purpose            - generic defines parametrization for a PML
c                        prism interpolating between two surfaces
c
c   arguments :
c     in:
c               No     - prism number
c               Eta    - reference coordinates of a point
c                        in the reference prism
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine prism_PML(No,Eta, X,Dxdeta)
c
c---------------------------------------------------------------------
      use GMP
c---------------------------------------------------------------------
      implicit none
      common /cprism_PML/ iprint
c---------------------------------------------------------------------
c     DUMMY ARGUMENTS
      integer, intent(in) :: No
      double precision, dimension(3), intent(in) :: Eta
      double precision, dimension(3), intent(out) :: X
      double precision, dimension(3,3), intent(out) :: Dxdeta
c---------------------------------------------------------------------
c     LOCAL VARIABLES
      integer :: nt1,nt2,lab1,lab2
      double precision, dimension(2) :: zetaB                            ! local coordinates of bottom face of master prism
      double precision, dimension(2) :: zetaT                            ! local coordinates of top face of master prism
      double precision, dimension(2,2) :: dzetaB                         ! derivatives
      double precision, dimension(2,2) :: dzetaT                         ! derivatives
      double precision, dimension(3) :: x1                               ! physical coordinates of bottom face of prism
      double precision, dimension(3) :: x2                               ! physical coordinates of top face of prism
      double precision, dimension(3,2) :: dx1, dx2                       ! derivatives
      integer :: j,iprint
c---------------------------------------------------------------------
c     PARAMETERS
      double precision, parameter, dimension(2,2,6) :: dzeta =           ! derivatives of local coordinates of triangular faces...
     .  reshape((/ 1.d0,  0.d0,  0.d0,  1.d0,                            ! ...with respect to master triangle coordinates (eta1,eta2)...
     .             0.d0, -1.d0,  1.d0, -1.d0,                            ! ...third index is orientation
     .            -1.d0,  1.d0, -1.d0,  0.d0,
     .             0.d0,  1.d0,  1.d0,  0.d0,
     .            -1.d0,  0.d0, -1.d0,  1.d0,
     .             1.d0, -1.d0,  0.d0, -1.d0/), (/2,2,6/))
c---------------------------------------------------------------------
c
ccc      iprint = 0
c
      if (PRISMS(No)%Type .ne. 'PMLprism') then
        write(*,7001) No
 7001   format('prism_PML: WRONG INPUT, No = ',i4)
        call print_GMP
        stop
      endif
c
      call decode(PRISMS(No)%FigNo(1), nt1,lab1)                         ! decode orientations of triangular faces
      call decode(PRISMS(No)%FigNo(2), nt2,lab2)

      if (iprint .eq. 1) then
        write(*,7002) nt1,nt2,lab1,lab2
 7002   format('prism_PML: nt1,nt2 = ',2i8,' lab1,lab2 = ',2i2)
      endif
c
      select case(lab1)                                                  ! account for orientation of bottom face
      case (0); zetaB(1:2) = Eta(1:2)
      case (1); zetaB(1) = Eta(2); zetaB(2) = 1.d0 - Eta(1) - Eta(2)
      case (2); zetaB(1) = 1.d0 - Eta(1) - Eta(2); zetaB(2) = Eta(1)
      case (3); zetaB(1) = Eta(2); zetaB(2) = Eta(1)
      case (4); zetaB(1) = 1.d0 - Eta(1) - Eta(2); zetaB(2) = Eta(2)
      case (5); zetaB(1) = Eta(1); zetaB(2) = 1.d0 - Eta(1) - Eta(2)
      case default
        write(*,7003)
 7003   format('prism_PML: unknown orientation of bottom face')
        stop
      end select
      dzetaB(1:2,1:2) = dzeta(1:2,1:2,lab1 + 1)
c
      select case(lab2)                                                  ! account for orientation of top face
      case (0); zetaT(1:2) = Eta(1:2)
      case (1); zetaT(1) = Eta(2); zetaT(2) = 1.d0 - Eta(1) - Eta(2)
      case (2); zetaT(1) = 1.d0 - Eta(1) - Eta(2); zetaT(2) = Eta(1)
      case (3); zetaT(1) = Eta(2); zetaT(2) = Eta(1)
      case (4); zetaT(1) = 1.d0 - Eta(1) - Eta(2); zetaT(2) = Eta(2)
      case (5); zetaT(1) = Eta(1); zetaT(2) = 1.d0 - Eta(1) - Eta(2)
      case default
        write(*,7004)
 7004   format('prism_PML: unknown orientation of top face')
        stop
      end select
      dzetaT(1:2,1:2) = dzeta(1:2,1:2,lab2 + 1)
c
      call trian(nt1,zetaB, x1,dx1)                                      ! physical coordinates of bottom and top face
      call trian(nt2,zetaT, x2,dx2)
      if (iprint.eq.1) then
        write(*,7005) zetaB
 7005   format('prism_PML: zetaB = ',2e12.5)
        do j=1,3
          write(*,7006) j,x1(j),dx1(j,1:2)
 7006     format('           ivar,x,dx = ',i1,2x,e12.5,2x,2e12.5)
        enddo
        write(*,7007) zetaT
 7007   format('prism_PML: zetaT = ',2e12.5)
        do j=1,3
          write(*,7006) j,x2(j),dx2(j,1:2)
        enddo
        call pause
      endif
c
      X(1:3) = x1(1:3)*(1.d0 - Eta(3)) + x2(1:3)*Eta(3)                  ! interpolate linearly between bottom and top face
      Dxdeta(1:3,1:2) = matmul(dx1,dzetaB)*(1.d0 - Eta(3)) +             ! compute derivative wrt eta1,eta2
     .                  matmul(dx2,dzetaT)*Eta(3)
      Dxdeta(1:3,3) = x2(1:3) - x1(1:3)                                  ! compute derivative wrt eta3
c
      end subroutine prism_PML


