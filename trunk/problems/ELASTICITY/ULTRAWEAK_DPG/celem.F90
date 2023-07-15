!--------------------------------------------------------------------------
!> @breif      Calculate modified stiffness matrix and load vector
!!
!! @param[in]  Mdle   - middle node of an element
!! @param[in]  Idec   - 1 nodes, 2 matrix
!!
!! @param[out] Nrdofs - # of element local dof for each physics
!! @param[out] Nrdofm - # of modified element dof in the expanded
!! @param[out] Nrdofc - # of modified element dof after compression
!!
!! @param[out] Nodm   - actual nodes in the order
!!
!! @param[out] NdofmH - the number of H1 dof
!! @param[out] NdofmE - the number of H(curl) dof
!! @param[out] NdofmV - the number of H(div) dof
!! @param[out] NdofmQ - the number of L2 dof
!!
!! @param[out] Nrnodm - number of the modified element nodes
!!
!! @param[out] Bload  - 1D array containing the modified load vector
!! @param[out] Astif  - 1D array containing the modified stiffness matrix
!!
!> @date       July 2023
!--------------------------------------------------------------------------
   subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,   &
                 Nodm,NdofmH,NdofmE,NdofmV,NdofmQ, &
                 Nrnodm, Bload,Astif)
!
      use physics
      use data_structure3D
!      use parameters
!
      implicit none
      integer,                      intent(in)  :: Mdle,Idec
      integer, dimension(NR_PHYSA), intent(out) :: Nrdofs
      integer,                      intent(out) :: Nrdofm,Nrdofc
      integer, dimension(MAXNODM),  intent(out) :: Nodm
      integer, dimension(MAXNODM),  intent(out) :: NdofmH,NdofmE,NdofmV,NdofmQ
      integer,                      intent(out) :: Nrnodm
      real*8,                       intent(out) :: Bload(*),Astif(*)
!
      integer :: nbcond(NRINDEX)
!
!--------------------------------------------------------------------------
!
!  ...this is a hack to eliminate the bubbles in the trace physics variables
!
!  ...physics attributes components (3+3+3+6 = 15 components)
!     phys|  comp   | description
!     (1) |  1 - 3  | H1 diplacement trace
!     (2) |  4 - 6  | Hdiv traction trace
!     (3) |  7 - 9  | L2 displacement field
!     (4) | 10 - 15 | L2 traction field
!
      nbcond = 0; nbcond(1:6) = 1
      call encod(nbcond,2,NRINDEX, NODES(Mdle)%bcond)
!
!  ...redirect to the system routine
      call celem_system(Mdle,Idec,                           &
                        Nrdofs,Nrdofm,Nrdofc,Nodm,           &
                        NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,  &
                        Bload,Astif)
!
!  ...reset the BC flags back to zero (unknown)
      NODES(Mdle)%bcond = 0
!
end subroutine celem
