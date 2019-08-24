!--------------------------------------------------------------------
!> Purpose - modify the order of approximation of the mesh
!!
!! @param[in] P   - new order of approximation in the x1-direction
!!
!> @date Dec 14
!--------------------------------------------------------------------
!
subroutine set_mesh_order(P)
!      
      use data_structure3D , only : NODES,NRNODS,Is_inactive
!
      implicit none      
      integer,intent(in) :: P
!      
      integer :: i
!
!--------------------------------------------------------------------
!
!     loop over all nodes and modify order
      do i=1,NRNODS
!
!       skip inactive nodes
        if (Is_inactive(i)) cycle
!
!       modify order
        call set_p(i, P,P,P)
!        
      enddo
!
!     update geometry and Dirichlet dofs
      call update_gdof
      call update_ddof
!
!
endsubroutine set_mesh_order
!
!
!
!--------------------------------------------------------------------
!> Purpose - modify the order of approximation of a node.
!!
!! @param[in] Nod - node number
!! @param[in] P   - new order of approximation in the x1-direction
!! @param[in] Q   - new order of approximation in the x2-direction
!! @param[in] R   - new order of approximation in the x3-direction
!!
!> @date Dec 14
!--------------------------------------------------------------------
!
subroutine set_p(Nod,P,Q,R)
!
      use data_structure3D , only : NODES,MAXP
!      
      implicit none
      integer, intent(in) :: Nod
      integer, intent(in) :: P,Q,R
!      
      character(len=4) :: ntype
      integer :: nord
!    
!--------------------------------------------------------------------
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
!     node type and order of approximation
      ntype=NODES(Nod)%type
!
!     encode new order of approximation
      select case (ntype)
      case('medg','mdlt','mdln','mdld'); nord =     P
      case('mdlq','mdlp');               nord =  10*P +    Q
      case('mdlb');                      nord = 100*P + 10*Q + R
      end select
!
!     modify node order of approximation (and everything related)
      call nodmod(Nod, nord)
!
!
endsubroutine set_p
