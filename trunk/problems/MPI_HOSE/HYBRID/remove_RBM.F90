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
  integer, dimension(NR_PHYSA) :: nbcond
!--------------------------------------------------------------------------
!
  ! call add_dirichlet_to_list(5)
  ! call add_dirichlet_to_list(6)
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
!  Vertex DOF are introduced at A,B,C,D
!
!  A
  pt=1
  nbcond = (/1,0,0,0,0/)
  ! nbcondtmp = NODES(NRELIS+pt)%bcond
  call encod(nbcond,10,NR_PHYSA, NODES(NRELIS+pt)%bcond)
  call set_index(NODES(NRELIS+pt)%case,NODES(NRELIS+pt)%bcond, NODES(NRELIS+pt)%index)
  ! NODES(NRELIS+pt)%bcond = nbcondtmp
!
!  B
  pt=2
  nbcond = (/6,6,0,0,0/)
  ! nbcondtmp = NODES(NRELIS+pt)%bcond
  call encod(nbcond,10,NR_PHYSA, NODES(NRELIS+pt)%bcond)
  call set_index(NODES(NRELIS+pt)%case,NODES(NRELIS+pt)%bcond, NODES(NRELIS+pt)%index)
  ! NODES(NRELIS+pt)%bcond = nbcondtmp
!
!  C
  pt=12
  nbcond = (/5,5,0,0,0/)
  ! nbcondtmp = NODES(NRELIS+pt)%bcond
  call encod(nbcond,10,NR_PHYSA, NODES(NRELIS+pt)%bcond)
  call set_index(NODES(NRELIS+pt)%case,NODES(NRELIS+pt)%bcond, NODES(NRELIS+pt)%index)
  ! NODES(NRELIS+pt)%bcond = nbcondtmp


! !  A
!   pt=1
!   nbcond = (/1,0,0,0,0/)
!   call encod(nbcond,10,NR_PHYSA, NODES(NRELIS+pt)%bcond)
!   call set_index(NODES(NRELIS+pt)%case,NODES(NRELIS+pt)%bcond, NODES(NRELIS+pt)%index)

! !  B
!   pt=24
!   nbcond = (/1,0,0,0,0/)
!   call encod(nbcond,10,NR_PHYSA, NODES(NRELIS+pt)%bcond)
!   call set_index(NODES(NRELIS+pt)%case,NODES(NRELIS+pt)%bcond, NODES(NRELIS+pt)%index)


end subroutine remove_RBM