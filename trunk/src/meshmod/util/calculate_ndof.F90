!> Purpose - calculate total degree of freedoms
subroutine calculate_ndof(NrdofH,NrdofE,NrdofV,NrdofQ)
  use data_structure3D
  implicit none
  integer, intent(out) :: NrdofH,NrdofE,NrdofV,NrdofQ

  integer, dimension(27) :: nodesl,norientl
  integer :: i,iel,mdle,nod,ndofH,ndofE,ndofV,ndofQ,nrnodes
  character(len=4) :: type

  call reset_visit

  NrdofH = 0;   NrdofE = 0;   NrdofV = 0;   NrdofQ = 0
  mdle = 0
  do iel=1, NRELES
     call nelcon(mdle, mdle)
     call elem_nodes(mdle, nodesl, norientl)

     type = NODES(mdle)%type
     nrnodes = nvert(type) + nedge(type) + nface(type) + 1
     do i=1,nrnodes
        nod = nodesl(i)

        if (NODES(nod)%visit.eq.0) then
           call ndof_nod(NODES(nod)%type, NODES(nod)%order, &
                ndofH,ndofE,ndofV,ndofQ)

           NrdofH = NrdofH + ndofH
           NrdofE = NrdofE + ndofE
           NrdofV = NrdofV + ndofV
           NrdofQ = NrdofQ + ndofQ

           NODES(nod)%visit = 1
        endif
     enddo
  enddo
end subroutine calculate_ndof


