!> Purpose : with given element refinement flag, return the face refinements
!! @param[in]  Type     - middle node type
!! @param[in]  Kref     - refinement flag for the middle node
!! @param[out] Kreff    - refinement flags for the element faces
subroutine find_face_ref_flags(Type,Kref, Kreff)
  character(len=4),      intent(in)  :: Type
  integer,               intent(in)  :: Kref
  integer, dimension(6), intent(out) :: Kreff

  Kreff(1:6)=0
  select case(Type)
  case('mdln')
     select case(Kref)
     case(11,12,13);        Kreff(1:4)       = 1
     case(24);              Kreff(1:4)       = (/1,3,1,3/)
     case(32);              Kreff((/1,2,3/)) = (/3,3,2/)
     case default;          write(*,7001) Type, Kref; call pause
     end select
  case('mdlp')
     select case(Kref)
     case(11);              Kreff(1:2) = 1; Kreff(3:5) = 11
     case(10);              Kreff(1:2) = 1; Kreff(3:5) = 10
     case(01);              Kreff(1:2) = 0; Kreff(3:5)=01
     case default;          write(*,7001) Type, Kref; call pause
     end select
  case('mdld')
     select case(Kref)
     case(10);              Kreff(1:5) = (/10,1,4,1,4/)
     case default;          write(*,7001) Type, Kref; call pause
     end select
  case('mdlb')
     select case(Kref)
     case(111);             Kreff(1:6)=11
     case(110);             Kreff(1:2)=11; Kreff(3:6)=10
     case(101);             Kreff(1:6) = (/10,10,11,01,11,01/)
     case(011);             Kreff(1:6) = (/01,01,01,11,01,11/)
     case(100);             Kreff( (/1,2,3,5/) ) = 10
     case(010);             Kreff( (/1,2,4,6/) ) = (/01,01,10,10/)
     case(001);             Kreff(3:6)=01
     case default;          write(*,7001) Type, Kref; call pause
     end select
  end select

7001 format('find_face_ref_flags: NOT SUPPORT Kref ',a4, i6)

end subroutine find_face_ref_flags
