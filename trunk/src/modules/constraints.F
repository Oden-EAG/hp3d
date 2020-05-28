c----------------------------------------------------------------------
c
c   module name        - module constraints
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 12
c
c   purpose            - contains info on constrained approximation;
c                        see routines constrs_util/setcontr_edge,
c                                     constrs_util/setcontr_trian,
c                        for the definition of arrays RRRH,RRTH
c
c
c----------------------------------------------------------------------

      module constraints
c
      use parameters
c
      save
      logical :: INITIALIZED_CONSTR = .FALSE.
c
c  ...constrained approximation coefficients for the master edge
      real(8) :: RRRH(3,MAXP-1,3,MAXP-1),
     .           RRRE(1,MAXP  ,2,MAXP  )
c
c  ...constrained approximation coefficients for the master triangle
      real(8) :: RRTH(4,MAXmdltH,7,MAXmdltH),
     .           RRTE(4,MAXmdltE,7,MAXmdltE),
     .           RRTQ(MAXmdltQ,4,MAXmdltQ)
c
c  ...constrained approximation coefficients for the master triangle
c     constraining a triangle and a quad
c     iref/p_nodes/p_dof/child_dof
c     for now, hcurl does not support this anisotropic constraints
      real(8) :: RRQH(2:4,1:4,MAXmdltH,MAXmdlqH)
c
c
      contains
c
      subroutine map_quad(Iref, Xi, X)
c
      integer :: Iref
      real(8), dimension(1:2,1:4,2:4) :: pts
c
      real(8), dimension(1:2) :: X, Xi
c
      pts = reshape(
     .  (/0.5d0,0.d0, 1.d0 ,0.d0, 0.d0 ,1.d0, 0.d0,0.5d0,
     .    0.d0 ,0.d0, 0.5d0,0.d0, 0.5d0,0.5d0, 0.d0,1.d0,
     .    0.d0 ,0.d0, 1.d0 ,0.d0, 0.5d0,0.5d0, 0.d0,0.5d0/),
     .  (/2,4,3/) )
c
      X(1:2) = pts(1:2,1,Iref)*(1.d0-Xi(1))*(1.d0-Xi(2))
     .       + pts(1:2,2,Iref)*Xi(1)*(1.d0-Xi(2))
     .       + pts(1:2,3,Iref)*Xi(1)*Xi(2)
     .       + pts(1:2,4,Iref)*(1.d0-Xi(1))*Xi(2)
c
      end subroutine map_quad
c
      function Get_rrqh(Iref, Np, Ip, J, I)
c
c  ...in
c       Iref :: refinement flag
c       I    :: index of constrained dof
c       J    :: index of parent dof
c       Np   :: the order of approximation
c
      integer Iref, Np, Ip, J, I, i1, i2, ii
      real(8) Get_rrqh
c
      i2 = (I-1)/(Np-1)+1
      i1 = I - (i2-1)*(Np-1)
c
      ii = (i2-1)*(MAXP-1)+i1
      Get_rrqh = RRQH(Iref, Ip, J, ii);
c
      end function Get_rrqh
c
      end module constraints
