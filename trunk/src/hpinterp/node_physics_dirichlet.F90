! routine that returns a logical vector corresponding to whether
! the physical attribute is Dirichlet at this node
! (similar to is_dirichlet() but expanding the physics list)
!
! Aug. 2020
! Jaime
!
! inputs:
! Mdle   Middle element node id
! Lnod   local node whose bc we want to check
! 
! Output
! Dflag  Vector of logical variables corresponding 
!        to a Dirichlet physical attribute at the node
! 
subroutine node_physics_dirichlet(Mdle,Lnod,Dflag)
use data_structure3D
implicit none
integer :: Mdle,Lnod
logical :: Dflag(NR_PHYSA)
!---------
integer :: nodesl(27),norientl(27),iphys,ibc(NR_PHYSA),loc,nod


Dflag = .false.
call get_connect_info(Mdle, nodesl,norientl)

nod = nodesl(Lnod)

call decod(NODES(Nod)%bcond,10,NR_PHYSA, ibc)

do iphys=1,NR_PHYSA
   if (ibc(iphys).eq.1) then
      Dflag(iphys) = .true.
   else
      call locate(ibc(iphys),DIRICHLET_LIST,NR_DIRICHLET_LIST, loc)
      if (loc.ne.0) then
         Dflag(iphys) = .true.
      endif
   endif
enddo


end subroutine