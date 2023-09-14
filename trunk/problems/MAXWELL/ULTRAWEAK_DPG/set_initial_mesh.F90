!-------------------------------------------------------------------------------
!> @brief       Define problem dependent data (multiphysics, BC, approximation)
!!
!> @param[out]  Nelem_order  - order of initial mesh elements
!!
!> @date        July 2023
!-------------------------------------------------------------------------------
subroutine set_initial_mesh(Nelem_order)
!
   use data_structure3D, only: NRELIS
   use physics         , only: NR_PHYSA
   use commonParam
!
   implicit none
!
!..polynomial order for initial mesh elements
   integer, intent(out) :: Nelem_order(NRELIS)
!
!..misc
   integer :: attr,bdom,comp,flag,i
   integer :: nr_attr,attr_list(NR_PHYSA)
!
!-------------------------------------------------------------------------------
!
!  STEP 1 : set initial element order
!
!..setting uniform isotropic polynomial order "IP"
   call set_order(IP, Nelem_order)
!
!  STEP 2 : set up physics
!
!..set up number of physical attributes supported by the element
   nr_attr = NR_PHYSA
   attr_list = [(i, i=1,NR_PHYSA)]
   call set_attr(nr_attr,attr_list)
!
!..BC flag
!  0 - no BC
!  1 - Dirichlet BC
!  2 - impedance BC via penalty term
!  3 - impedance BC via elimination
!
!..physics attributes components
!  attr | comp | index | description
!     1 |  1-2 |   1-2 | Hcurl for Maxwell trace (\hat E,\hat H) (2 components)
!     2 |  1-6 |   3-8 | L2 field for Maxwell (E,H) (6 components)
   attr = 1 ! set BC for Hcurl variables
!
!..boundary domain "0" (Dirichlet BC on E-trace)
   bdom = 0 ! set on all exterior faces with boundary domain "0" (default domain)
   comp = 1 ! E-trace
   flag = 1 ! Dirichlet BC flag
   call set_bcond(bdom,attr,comp,flag)
!
!..boundary domain "1" (BC depends on IBCFLAG)
   bdom = 1 ! set on all exterior faces with boundary domain "1"
   if (IBCFLAG.eq.2..or.IBCFLAG.eq.3) then
!  ...impedance BC on H-trace
      comp = 2       ! H-trace
      flag = IBCFLAG ! impedance BC flag
   else
!  ...Dirichlet BC on E-trace
      comp = 1 ! E-trace
      flag = 1 ! Dirichlet BC flag
   endif
   call set_bcond(bdom,attr,comp,flag)
!
end subroutine set_initial_mesh
