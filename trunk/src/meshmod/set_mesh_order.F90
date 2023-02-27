!--------------------------------------------------------------------
!> @brief modify the order of approximation of the mesh
!!
!> @param[in] P - new order of approximation in the x1-direction
!!
!> @date  Feb 2023
!--------------------------------------------------------------------
subroutine set_mesh_order(P)
!
   use data_structure3D , only : NRNODS,Is_inactive
!
   implicit none
   integer, intent(in) :: P
!
   integer :: i
!
!--------------------------------------------------------------------
!TODO: parallelize with OMP, if this routine is called
!
!..loop over all nodes and modify order
   do i=1,NRNODS
!
!  ...skip inactive nodes
        if (Is_inactive(i)) cycle
!
!  ...modify order
      call set_p(i, P,P,P)
!
   enddo
!
!..update geometry and Dirichlet dofs
   call update_gdof
   call update_ddof
!
end subroutine set_mesh_order
!
!--------------------------------------------------------------------
!> @brief modify the order of approximation of a node.
!!
!! @param[in] Nod - node number
!! @param[in] P   - new order of approximation in the x1-direction
!! @param[in] Q   - new order of approximation in the x2-direction
!! @param[in] R   - new order of approximation in the x3-direction
!!
!> @date Feb 2023
!--------------------------------------------------------------------
subroutine set_p(Nod,P,Q,R)
!
   use data_structure3D , only : NODES,MAXP
   use node_types
!
   implicit none
   integer, intent(in) :: Nod,P,Q,R
!
   integer :: ntype,nord
!
!!!      write(*,*)'Nod,P,Q,R [0] = ',Nod,P,Q,R
!!!!     check upper bound
!!!      P=min(P,MAXP)
!!!      Q=min(Q,MAXP)
!!!      R=min(R,MAXP)
!!!!
!!!      write(*,*)'Nod,P,Q,R [1] = ',Nod,P,Q,R
!!!
!!!!     check lower bound
!!!      P=max(P,1)
!!!      Q=max(Q,1)
!!!      R=max(R,1)
!
!..node type and order of approximation
   ntype=NODES(Nod)%ntype
!
!..encode new order of approximation
   select case (ntype)
      case(MEDG,MDLT,MDLN,MDLD); nord =     P
      case(MDLQ,MDLP);           nord =  10*P +    Q
      case(MDLB);                nord = 100*P + 10*Q + R
   end select
!
!..modify node order of approximation (and everything related)
   call nodmod(Nod, nord)
!
end subroutine set_p
