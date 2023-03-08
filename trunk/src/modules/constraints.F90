!----------------------------------------------------------------------
!
!   module name        - constraints
!
!----------------------------------------------------------------------
!> @brief   contains info on constrained approximation;
!!          see routines constrs_util/setcontr_edge,
!!                       constrs_util/setcontr_trian,
!!          for the definition of arrays RRRH,RRTH
!> @date    Feb 2023
!----------------------------------------------------------------------
module constraints
!
   use node_types
   use parameters
   implicit none
!
   save
!
   logical :: INITIALIZED_CONSTR = .false.
!
!..constrained approximation coefficients for the master edge
   real(8) :: RRRH(3,MAXP-1,3,MAXP-1), &
              RRRE(1,MAXP  ,2,MAXP  )
!
!..constrained approximation coefficients for the master triangle
   real(8) :: RRTH(4,MAXmdltH,7,MAXmdltH), &
              RRTE(4,MAXmdltE,7,MAXmdltE), &
              RRTQ(  MAXmdltQ,4,MAXmdltQ)
!
!..constrained approximation coefficients for the master triangle
!  constraining a triangle and a quad
!  iref/p_nodes/p_dof/child_dof
!  for now, H(curl) does not support these anisotropic constraints
   real(8) :: RRQH(2:4,1:4,MAXmdltH,MAXmdlqH)
!
!
contains
!
!----------------------------------------------------------------------
!> @date  Feb 2023
   subroutine map_quad(Iref,Xi, X)
!
      integer, intent(in)  :: Iref
      real(8), intent(in)  :: Xi(1:2)
      real(8), intent(out) :: X (1:2)
      real(8)              :: pts(1:2,1:4,2:4)
!
      pts = reshape(                                           &
        (/0.5d0,0.d0, 1.d0 ,0.d0, 0.d0 ,1.d0, 0.d0,0.5d0,      &
          0.d0 ,0.d0, 0.5d0,0.d0, 0.5d0,0.5d0, 0.d0,1.d0,      &
          0.d0 ,0.d0, 1.d0 ,0.d0, 0.5d0,0.5d0, 0.d0,0.5d0/),   &
        (/2,4,3/) )
!
      X(1:2) = pts(1:2,1,Iref)*(1.d0-Xi(1))*(1.d0-Xi(2))    &
             + pts(1:2,2,Iref)*Xi(1)*(1.d0-Xi(2))           &
             + pts(1:2,3,Iref)*Xi(1)*Xi(2)                  &
             + pts(1:2,4,Iref)*(1.d0-Xi(1))*Xi(2)
!
   end subroutine map_quad
!
!----------------------------------------------------------------------
!  ...in
!       Iref :: refinement flag
!       Np   :: the order of approximation
!       Ip   ::
!       J    :: index of parent dof
!       I    :: index of constrained dof
!> @date  Feb 2023
   function Get_rrqh(Iref, Np, Ip, J, I)
!
      integer Iref, Np, Ip, J, I
      real(8) Get_rrqh
      integer i1, i2, ii
!
      i2 = (I-1)/(Np-1)+1
      i1 = I - (i2-1)*(Np-1)
!
      ii = (i2-1)*(MAXP-1)+i1
      Get_rrqh = RRQH(Iref, Ip, J, ii);
!
   end function Get_rrqh
!
end module constraints
