!----------------------------------------------------------------------
!> @Purpose : Routine checks orientation of a GMP block by computing
!>            Jacobian of the linear parametrization at its 1st vertex,
!>            and changes enumeration of vertices in the case of a
!>            negative Jacobian
!
!> @param[in] Lab  = 1 (prism), 2(hexa), 3(tet), 4(pyramid)
!> @param[in] Nb   - block number
!
!> @revision Feb 13
!---------------------------------------------------------------------
subroutine check_orientation(Lab,Nb)
!
      use GMP
      implicit none
      integer, intent(in)     :: Lab,Nb
!
      real(8), dimension(3,3) :: vect
      integer                 :: np1,np2,np3,np4,iprint
      real(8)                 :: det
!---------------------------------------------------------------------
!
      iprint=0
!
      select case(Lab)
      case(1)
        np1=PRISMS(Nb)%VertNo(1); np2=PRISMS(Nb)%VertNo(2)
        np3=PRISMS(Nb)%VertNo(3); np4=PRISMS(Nb)%VertNo(4)
      case(2)
        np1=HEXAS(Nb)%VertNo(1); np2=HEXAS(Nb)%VertNo(2)
        np3=HEXAS(Nb)%VertNo(4); np4=HEXAS(Nb)%VertNo(5)
      case(3)
        np1=TETRAS(Nb)%VertNo(1); np2=TETRAS(Nb)%VertNo(2)
        np3=TETRAS(Nb)%VertNo(3); np4=TETRAS(Nb)%VertNo(4)
      case(4)
        np1=PYRAMIDS(Nb)%VertNo(1); np2=PYRAMIDS(Nb)%VertNo(2)
        np3=PYRAMIDS(Nb)%VertNo(4); np4=PYRAMIDS(Nb)%VertNo(5)
      endselect
!
      vect(1:3,1) = POINTS(np2)%Rdata(1:3) - POINTS(np1)%Rdata(1:3)
      vect(1:3,2) = POINTS(np3)%Rdata(1:3) - POINTS(np1)%Rdata(1:3)
      vect(1:3,3) = POINTS(np4)%Rdata(1:3) - POINTS(np1)%Rdata(1:3)
      call mixed_product(vect(1:3,1),vect(1:3,2),vect(1:3,3), det)
!
      if (iprint == 1) then
        write(*,7001) Lab,Nb,det
 7001   format(' check_orientation: Lab,Nb,det = ',i2,i6,e12.5)
      endif
!
      if (det <= 0.d0) then
        select case(Lab)
        case(1)
          call swap(PRISMS(Nb)%VertNo(2),PRISMS(Nb)%VertNo(3))
          call swap(PRISMS(Nb)%VertNo(5),PRISMS(Nb)%VertNo(6))
        case(2)
          call swap(HEXAS(Nb)%VertNo(2),HEXAS(Nb)%VertNo(4))
          call swap(HEXAS(Nb)%VertNo(6),HEXAS(Nb)%VertNo(8))
        case(3)
          call swap(TETRAS(Nb)%VertNo(2),TETRAS(Nb)%VertNo(3))
        case(4)
          call swap(PYRAMIDS(Nb)%VertNo(1),PYRAMIDS(Nb)%VertNo(2))
          call swap(PYRAMIDS(Nb)%VertNo(3),PYRAMIDS(Nb)%VertNo(4))
        endselect
      endif
!
!
end subroutine check_orientation
