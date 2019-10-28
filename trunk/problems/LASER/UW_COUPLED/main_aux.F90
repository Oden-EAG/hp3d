!
!-----------------------------------------------------------------------
!                   subroutine set_ANISO_FLAG
!-----------------------------------------------------------------------
!
subroutine set_ANISO_FLAG(ref_xyz)
!
   use commonParam
!
   implicit none
!
   integer, intent(in) :: ref_xyz
!
   select case(ref_xyz)
      case(1); ANISO_FLAG = IREFINE_X
      case(2); ANISO_FLAG = IREFINE_Y
      case(3); ANISO_FLAG = IREFINE_Z
      case(4); ANISO_FLAG = IREFINE_XY
      case(5); ANISO_FLAG = IREFINE_XZ
      case(6); ANISO_FLAG = IREFINE_YZ
      case default
         write(*,*) 'set_ANISO_FLAG: invalid ref_xyz param. stop.'
         stop
   end select
!
end subroutine set_ANISO_FLAG
!
!
!-----------------------------------------------------------------------
!                   subroutine select_phys_problem
!-----------------------------------------------------------------------
!
subroutine select_phys_problem(NO_PROBLEM)
!
   implicit none
!
   integer, intent(out) :: NO_PROBLEM
   integer              :: numProb
!
   write(*,*) 'Select problem:  1 = Heat Step'
   write(*,*) '                 2 = Multi-step Heat equation'
   write(*,*) '                 3 = Time harmonic Maxwell (signal)'
   write(*,*) '                 4 = Time harmonic Maxwell (pump)'
   read(*,*) numProb
   NO_PROBLEM = numProb
!
end subroutine select_phys_problem
!
!
!-----------------------------------------------------------------------
!                   subroutine select_phys_problem_maxwell
!-----------------------------------------------------------------------
!
subroutine select_phys_problem_maxwell(NO_PROBLEM)
!
   implicit none
!
   integer, intent(out) :: NO_PROBLEM
   integer              :: numProb
!
   write(*,*) 'Select problem:  3 = Time harmonic Maxwell (signal)'
   write(*,*) '                 4 = Time harmonic Maxwell (pump)'
   read(*,*) numProb
   NO_PROBLEM = numProb
!
end subroutine select_phys_problem_maxwell
!
!
!-----------------------------------------------------------------------
!                   subroutine select_phys_problem_heat
!-----------------------------------------------------------------------
!
subroutine select_phys_problem_heat(NO_PROBLEM)
!
   implicit none
!
   integer, intent(out) :: NO_PROBLEM
   integer              :: numProb
!
   write(*,*) 'Select problem:  1 = Heat Step'
   write(*,*) '                 2 = Multi-step Heat equation'
   read(*,*) numProb
   NO_PROBLEM = numProb
!
end subroutine select_phys_problem_heat
!
!
!-----------------------------------------------------------------------
!                       subroutine set_physAm
!-----------------------------------------------------------------------
!
subroutine set_physAm(NO_PROBLEM, PhysNick,Flag)
!
   use physics, only: PHYSAm
!
   implicit none
!
   integer,               intent(in)  :: NO_PROBLEM
   integer,               intent(out) :: PhysNick
   integer, dimension(6), intent(out) :: Flag
!
   select case(NO_PROBLEM)
      case(1,2)
         PhysNick = 1000
         Flag=0; Flag(1)=1
         PHYSAm(1:6) = (/.true.,.false.,.false.,.true.,.false.,.false./)
      case(3)
         PhysNick = 1
         Flag=0; Flag(5)=1
         PHYSAm(1:6) = (/.false.,.true.,.false.,.false.,.true.,.false./)
      case(4)
         PhysNick = 1
         Flag=0; Flag(6)=1
         PHYSAm(1:6) = (/.false.,.false.,.true.,.false.,.false.,.true./)
      case default
      write(*,*) 'set_physAm: invalid NO_PROBLEM param. stop.'
         stop
   end select
!
end subroutine set_physAm
!
!
!-------------------------
! subroutine ref_core
! refine fiber core in xy
!-------------------------
subroutine ref_core()
   use data_structure3D
   use assembly
!
   implicit none
!
   integer, parameter :: kref_mdlb = 110
   integer, parameter :: kref_mdlp = 10
!
   integer :: mdle, iel, ndom, nr_elem_ref
!
!..geometry dof (work space for nodcor)
   integer :: n_elem(NRELES)
!
   nr_elem_ref = 0
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
!  ...check if mdle is in the fiber core 
      call find_domain(mdle, ndom)
      select case(ndom)
         case(1,2)
            nr_elem_ref = nr_elem_ref + 1
            n_elem(nr_elem_ref) = mdle
         case default
            ! do nothing
      end select
   enddo
!
   write(*,2410) nr_elem_ref
 2410 format(' found ', i6, ' elements to refine.')
!
   do iel=1,nr_elem_ref
      mdle = n_elem(iel)
      write(*,2400) mdle
!  ...refine element mdle
      select case(NODES(mdle)%type)
         case('mdlb')
            call refine(mdle,kref_mdlb)
         case('mdlp')
            call refine(mdle,kref_mdlp)
         case default
         write(*,*) 'ref_core: unsupported node type. stop.'
         stop
      end select
   enddo
!
 2400 format(' refining mdle = ', i6)
!
end subroutine ref_core
!
!
!-----------------------------------------------------------------------
!                       subroutine ref_layer
!-----------------------------------------------------------------------
!
subroutine ref_layer(Kref)
!
   use data_structure3D
   use assembly
!
   implicit none
!
   integer, intent(in) :: Kref
!
   integer :: mdle, iel, nr_elem_ref
!
!..geometry dof (work space for nodcor)
   real*8, dimension(3,MAXbrickH) :: xnod
   integer :: n_elem(NRELES)
!
   nr_elem_ref = 0
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
!  ...check if mdle is at the fiber input
      call nodcor(mdle, xnod)
      ! For HEXA only
      if (minval(xnod(3,1:8)) .lt. 0.5d0) then
         nr_elem_ref = nr_elem_ref + 1
         n_elem(nr_elem_ref) = mdle
      endif
   enddo
!
   write(*,2510) nr_elem_ref
 2510 format(' found ', i6, ' elements to refine.')
!
   do iel=1,nr_elem_ref
      mdle = n_elem(iel)
! ...refine element mdle
      write(*,2500) mdle
      call refine(mdle,Kref)
   enddo
!
 2500 format(' refining mdle = ', i6)
!
end subroutine ref_layer

