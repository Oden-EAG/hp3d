!----------------------------------------------------------------------------
!> Purpose : return stiffness matrix and load vector for element
!!
!! @param[in]  Mdle   - element (middle node) number
!! @param[out] Nrdof  - number of dof for a single component
!! @param[out] Itest  - index for assembly 
!! @param[out] Itrial - index for assembly 
!----------------------------------------------------------------------------
!
subroutine elem(Mdle, Itest,Itrial)
!
  use error
  use data_structure3D , only : NODES
  use parameters       , only : ZERO
  use physics          , only : NR_PHYSA
  use assembly         , only : ALOC, BLOC
!----------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!----------------------------------------------------------------------------
  !
  Itest (1:NR_PHYSA) = 0
  Itrial(1:NR_PHYSA) = 0
  Itest(1)=1; Itrial(1)=1
  !
  ALOC(1,1)%array = ZERO ; BLOC(1)%array = ZERO
  !
  select case(NODES(Mdle)%case)
  ! electromagnetic
  case(1)
     call elem_em(Mdle,BLOC(1  )%array,BLOC(1  )%nrow,  &
                       ALOC(1,1)%array,ALOC(1,1)%ncol)
  case default
     call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
  end select
!
end subroutine elem
