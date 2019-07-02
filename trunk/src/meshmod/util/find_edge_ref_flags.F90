!> Purpose : with given element refinement flag, return the edge refinements
!! @param[in]  Type     - middle node type
!! @param[in]  Kref     - refinement flag for the middle node 
!! @param[out] Krefe    - refinement flags for the edges
subroutine find_edge_ref_flags(Type,Kref, Krefe)
  character(len=4),       intent(in)  :: Type
  integer,                intent(in)  :: Kref
  integer, dimension(12), intent(out) :: Krefe

  Krefe=0
  select case(Type)
  case('mdln')
     select case(Kref)
     case(11,12,13);          Krefe(1:6) = 1
     case(24); Krefe(1:6) = 1; Krefe(kref-20)=0
     case(32);                Krefe( (/1,2,5/) ) = 1
     case default;            write(*,7001) Type, Kref; call pause
     end select
  case('mdlp')
     select case(Kref)
     case(11);                Krefe(1:9)=1
     case(10);                Krefe(1:6)=1
     case(01);                Krefe(7:9)=1
     case default;            write(*,7001) Type, Kref; call pause
     end select
  case('mdld')
     select case(Kref)
     case(10);                Krefe(1:8)=(/1,0,1,0,1,1,1,1/)
     case default;            write(*,7001) Type, Kref; call pause
     end select
  case('mdlb')
     select case(Kref)
     case(111);               Krefe(1:12)=1
     case(110);               Krefe(1:8)=1
     case(101);               Krefe( (/1,3,5,7/) ) = 1; Krefe(9:12)=1
     case(011);               Krefe( (/2,4,6,8/) ) = 1; Krefe(9:12)=1
     case(100);               Krefe( (/1,3,5,7/) ) = 1
     case(010);               Krefe( (/2,4,6,8/) ) = 1
     case(001);               Krefe(9:12)=1
     case default;            write(*,7001) Type, Kref; call pause
     end select
  end select

7001 format('find_edge_ref_flags: NOT SUPPORT Kref ',a4, i3)

end subroutine find_edge_ref_flags
