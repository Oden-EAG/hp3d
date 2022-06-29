!------------------------------------------------------------------------------------
!> Purpose : return stiffness matrix and load vector for element
!!
!! @param[in]  Mdle   - an element (middle node) number
!! @param[out] Nrdof  - number of dof for a single component
!! @param[out] Itest  - index for assembly
!! @param[out] Itrial - index for assembly
!------------------------------------------------------------------------------------
!
subroutine elem(Mdle, Itest,Itrial)
  use parameters, only : ZERO
  use physics   , only : NR_PHYSA
  use assembly  , only : ALOC,BLOC,NR_RHS
  use data_structure3D
!------------------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!
!------------------------------------------------------------------------------------
!     I N I T I A L I Z E
!------------------------------------------------------------------------------------
!
  Itest (1:NR_PHYSA) = 0
  Itrial(1:NR_PHYSA) = 0
!
!------------------------------------------------------------------------------------
!     C O N S I D E R    C A S E S
!------------------------------------------------------------------------------------
!
  select case(NODES(Mdle)%case)
! PRIMAL
  case(24)  ! 2^4+2^3
    Itest(1:2)=1; Itrial(1:2)=1
    call elem_DPG_PRIMAL(Mdle)
! ULTRA-WEAK
  case(31)  ! 2^4+2^3+2^2+2^1+2^0
    Itest(1:5)=1; Itrial(1:5)=1
    call elem_DPG_UWEAK(Mdle)
  case default
    write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
               Mdle,NODES(Mdle)%case
    call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
  end select
!
end subroutine elem
