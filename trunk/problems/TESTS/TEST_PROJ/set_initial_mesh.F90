!----------------------------------------------------------------------
!> Purpose : set problem dependent attributes for the initial mesh
!            (multiphysics data, boundary conditions, and initial 
!            order of approximation.)
!
!  @param[out]  Nelem_order - order of approximation for initial mesh
!                             elements
!----------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
!
      use PROJ             , only : IORDER_PRIS, &
                                    IORDER_BRIC, &
                                    IORDER_TETR, &
                                    IORDER_PYRA
      use parameters       , only : MAXP
      use data_structure3D , only : ELEMS,NRELIS
!
      implicit none
      integer,dimension(NRELIS),intent(out) :: Nelem_order
!
      integer :: nel
      integer,dimension(6) :: ibc
!----------------------------------------------------------------------
!
!     loop through initial mesh elements
      do nel=1,NRELIS
!
!       all elements support 4 attributes
        ELEMS(nel)%nrphysics=4
!
!       assign attributes
        allocate(ELEMS(nel)%physics(4)) 
        ELEMS(nel)%physics(1)='projH'
        ELEMS(nel)%physics(2)='projE'
        ELEMS(nel)%physics(3)='projV'
        ELEMS(nel)%physics(4)='projQ'
!
!       set up order of approximationn
        select case(ELEMS(nel)%Type)
        case('pris') ; Nelem_order(nel)=IORDER_PRIS*11
        case('bric') ; Nelem_order(nel)=IORDER_BRIC*111
        case('tetr') ; Nelem_order(nel)=IORDER_TETR*1
        case('pyra') ; Nelem_order(nel)=IORDER_PYRA*1
        endselect
!
!       no BC's on all faces
        ibc(1:6)=0
!
!       assign BC's for all attributes
        allocate(ELEMS(nel)%bcond(4)) 
        call encodg(ibc,10,6, ELEMS(nel)%bcond(1))
        call encodg(ibc,10,6, ELEMS(nel)%bcond(2))
        call encodg(ibc,10,6, ELEMS(nel)%bcond(3))
        call encodg(ibc,10,6, ELEMS(nel)%bcond(4))
!
      enddo
!
!
endsubroutine set_initial_mesh
