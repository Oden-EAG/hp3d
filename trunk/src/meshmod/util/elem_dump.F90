!> Purpose : dump the initial element's nodal connectivity
!! @param[in]  Mdle     - middle node
!! @param[out] Nodesl   - nodes
!! @param[out] Norientl - orients
subroutine elem_dump(Mdle, Nodesl, Norientl)
  use data_structure3D
  implicit none
  integer,                intent(in)  :: Mdle
  integer, dimension(27), intent(out) :: Nodesl, Norientl

  if (Mdle.gt.NRELIS) then
     write(*,*) 'elem_dump: MDLE, NRELIS = ', Mdle, NRELIS
     write(*,*) 'elem_dump: OUT OF RANGE'
  endif

  select case(ELEMS(Mdle)%type)
#ifdef _PYHP3D
  case(NODE_TYPE_TETRA)
#else
  case('tetr')
#endif
     Nodesl(1:15) = ELEMS(Mdle)%nodes(1:15)
     call decodg(ELEMS(Mdle)%edge_orient,2,6,  Norientl(5:10))
     call decodg(ELEMS(Mdle)%face_orient,8,4,  Norientl(11:14))
#ifdef _PYHP3D
  case(NODE_TYPE_PRISM)
#else
  case('pris')
#endif
     Nodesl(1:21) = ELEMS(Mdle)%nodes(1:21)
     call decodg(ELEMS(Mdle)%edge_orient,2,9,  Norientl(7:15))
     call decodg(ELEMS(Mdle)%face_orient,8,5,  Norientl(16:20))
#ifdef _PYHP3D
  case(NODE_TYPE_PYRAMID)
#else
  case('pyra')
#endif
     Nodesl(1:19) = ELEMS(Mdle)%nodes(1:19)
     call decodg(ELEMS(Mdle)%edge_orient,2,8,  Norientl(6:13))
     call decodg(ELEMS(Mdle)%face_orient,8,5,  Norientl(14:18))
#ifdef _PYHP3D
  case(NODE_TYPE_HEXA)
#else
  case('bric')
#endif
     Nodesl(1:27) = ELEMS(Mdle)%nodes(1:27)
     call decodg(ELEMS(Mdle)%edge_orient,2,12, Norientl(9:20))
     call decodg(ELEMS(Mdle)%face_orient,8,6,  Norientl(21:26))
  end select

end subroutine elem_dump
