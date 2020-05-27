!> Purpose : Write middle node connectivity
subroutine mdle_with_visit_numbering(Nout, Nbegin, Iflag)
  use data_structure3D
  use element_data
  implicit none
  integer, intent(in) :: Nout, Iflag, Nbegin
  character(len=4) :: type
  integer :: iv, iel, mdle, ndom, nr_vert, nverl(8), ntmp(3), nodesl(27), norientl(27)

  mdle = 0
  do iel=1,NRELES

     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl,norientl)

     type    = NODES(mdle)%type
     nr_vert = nvert(type)

     if (type.ne.'mdln') then
        cycle
     end if

     nverl = 0
     do iv=1,nr_vert
        nverl(iv) = abs(NODES(nodesl(iv))%visit) - Nbegin
     end do

     call find_domain(mdle, ndom)

     if (Iflag.eq.1) then
        write(Nout,7000) nr_vert, nverl(1:nr_vert)
     else
        write(Nout,7001) ndom, nverl(1:nr_vert)
     end if

  enddo

7000 format(i10, 8i10, '    # element id, vertices list')
7001 format(i10, 8i10, '    # element id, vertices list')

end subroutine mdle_with_visit_numbering
