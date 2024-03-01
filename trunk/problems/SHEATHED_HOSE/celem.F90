!> @brief calculate modified stiffness matrix and load vector
!> @param[in]  Mdle   - middle node of an element
!> @param[in]  Idec   - 1 nodes, 2 matrix
!!
!> @param[out] Nrdofs - # of element local dof for each physics
!> @param[out] Nrdofm - # of modified element dof in the expanded
!> @param[out] Nrdofc - # of modified element dof after compression
!!
!> @param[out] Nodm   - actual nodes in the order
!!
!> @param[out] NdofmH - the number of H1 dof
!> @param[out] NdofmE - the number of H(curl) dof
!> @param[out] NdofmV - the number of H(div) dof
!> @param[out] NdofmQ - the number of L2 dof
!!
!> @param[out] Nrnodm - number of the modified element nodes
!!
!> @param[out] Bload  - 1D array containing the modified load vector
!> @param[out] Astif  - 1D array containing the modified stiffness matrix
!!
subroutine celem( &
     Mdle,Idec, &
     Nrdofs,Nrdofm,Nrdofc, &
     Nodm, &
     NdofmH,NdofmE,NdofmV,NdofmQ, &
     Nrnodm, &
     Bload,Astif)
  use physics
  use data_structure3D
  ! use parameters
!--------------------------------------------------------------------------
  implicit none
  integer,                      intent(in)  :: Mdle,Idec
  integer, dimension(NR_PHYSA), intent(out) :: Nrdofs
  integer,                      intent(out) :: Nrdofm,Nrdofc
  integer, dimension(MAXNODM),  intent(out) :: Nodm
  integer, dimension(MAXNODM),  intent(out) :: NdofmH,NdofmE,NdofmV,NdofmQ
  integer,                      intent(out) :: Nrnodm
  real(8),                       intent(out) :: Bload(*),Astif(*)
!--------------------------------------------------------------------------
  integer, dimension(NRINDEX_HEV) :: nbcond
!--------------------------------------------------------------------------

! Physics attributes:               Components:
! 1 TrDis  contin  (3 components)    1-3
! 2 TrStr  normal  (3 components)    4-6
! 3 Displ  discon  (3 components)    7-9
! 4 Stres  discon  (6 components)   10-15
! 5 Omega  discon  (3 components)   16-18

  select case(NODES(Mdle)%case)
! PRIMAL
  case(24)  ! 2^4+2^3

    ! eliminate middle node H(div) dof in celem_system by using BC flag
    nbcond(1:3) = 0
    nbcond(4:6) = 1
    call encod(nbcond,2,NRINDEX_HEV, NODES(Mdle)%bcond)

! ULTRA-WEAK
  case(31)  ! 2^4+2^3+2^2+2^1+2^0

    ! eliminate middle node H^1 and H(div) dof in celem_system by using BC flag
    nbcond(1:6) = 1
    call encod(nbcond,2,NRINDEX_HEV, NODES(Mdle)%bcond)

  case default

    write(*,*) 'ERROR: unknown formulation'
    stop 1

  end select

! redirect to the system routine
  call celem_system(Mdle,Idec, &
                    Nrdofs,Nrdofm,Nrdofc, &
                    Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm, &
                    Bload,Astif)

!
! reset the BC flag to zero
  select case(NODES(Mdle)%case)
  case(24,31)
    NODES(Mdle)%bcond = 0
  end select

end subroutine celem
