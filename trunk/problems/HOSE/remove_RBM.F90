!> Purpose : eliminate infinitesimal rigid body motions from problem
!!
subroutine remove_RBM
  use physics
  use data_structure3D
  ! use parameters
!--------------------------------------------------------------------------
  implicit none
!--------------------------------------------------------------------------
  integer :: pt
  integer :: nbcondtmp
  integer, dimension(NRINDEX) :: nbcond
!--------------------------------------------------------------------------
!
! This is a hack to eliminate vertex DOF in the trace variables.
! The 'known' dof are eliminated by static condensation locally
!       ___ ___ _____________________________
!      /8__C__6\                             \
!     / /7_|_5\ \                             \
!    / /_/   \_\B\  _  _  _  _  _  _  _  _  _  \       x--->
!    \ \4\_ _/2/ /                             /
!     \3\__|__/1/                             /
!      \___A___/_____________________________/
!
!  Vertex DOF are introduced at A,B,C
!
!  A
  pt=1
  nbcond = (/1,1,1, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0/) ! 3 H1 + 3 H(div) + 12 L2
  call encod(nbcond,2,NRINDEX, NODES(NRELIS+pt)%bcond)
!
!  B
  pt=2
  nbcond = (/1,0,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0/)
  call encod(nbcond,2,NRINDEX, NODES(NRELIS+pt)%bcond)
!
!  C
  pt=12
  nbcond = (/1,1,0, 0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0/)
  call encod(nbcond,2,NRINDEX, NODES(NRELIS+pt)%bcond)

end subroutine remove_RBM
