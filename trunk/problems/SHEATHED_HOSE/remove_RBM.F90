!> Purpose : eliminate infinitesimal rigid body motions from problem
!!
subroutine remove_RBM
  use physics
  use data_structure3D
!
!--------------------------------------------------------------------------
  implicit none
!--------------------------------------------------------------------------
  integer :: pt
  integer :: nbcond(NRINDEX)
!--------------------------------------------------------------------------
!
! This is a hack to eliminate vertex DOF in the trace variables.
! Six rigid body motion DOFs are assigned to be removed during assembly
!       ___ ___ _______________
!      /8__ __6\               \
!     / /7_|_5\ \               \
!   C/ /_/   \_\ \B_ _ _ _ _ _ _ \      x--->
!    \ \4\_ _/2/ /               /
!     \3\__|__/1/               /
!      \___ ___/_______________/
!          A
!  Elements 1,3,6,8 : Steel
!  Elements 2,4,5,7 : Rubber
!  Vertex DOFs are removed at A,B,C
!
!  A
  pt=1
  nbcond = (/1,1,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0/) ! 3 H1 + 3 H(div) + 12 L2
  call encod(nbcond,2,NRINDEX, NODES(NRELIS+pt)%bcond)
!
!  B
  pt=2
  nbcond = (/1,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0/)
  call encod(nbcond,2,NRINDEX, NODES(NRELIS+pt)%bcond)
!
!  C
  pt=11
  nbcond = (/1,0,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0/)
  call encod(nbcond,2,NRINDEX, NODES(NRELIS+pt)%bcond)

end subroutine remove_RBM
