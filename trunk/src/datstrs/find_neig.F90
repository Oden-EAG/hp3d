!--------------------------------------------------------------------------
!> Purpose : find neighbors (up to 4, for an h4 refined face) across faces
!! @param[in]  Mdle - middle node
!! @param[out] Neig - neighbors
!!
!! @revision May 12
!--------------------------------------------------------------------------
!
subroutine find_neig(Mdle, Neig_list)
  use element_data
  use data_structure3D
  implicit none
  ! ** Arguements
  !---------------------------------------------------
  integer,                 intent(in)  :: Mdle
  integer, dimension(4,6), intent(out) :: Neig_list

  ! ** Locals
  !---------------------------------------------------
  character(len=4) :: type
  integer, dimension(27) :: nodesl,norientl
  integer, dimension(2)  :: neig, nsid_list,norient_list
  integer :: iprint, i, nod, nrneig
  !---------------------------------------------------
  iprint = 0

  !  ...initialize
  Neig_list(1:4,1:6)=0
  type=NODES(Mdle)%type
  if (iprint.eq.1) then
     write(*,*) 'find_neig: Mdle, type = ', Mdle, type
  endif
  
  !---------------------------------------------------
  ! Step 0: short cut for initial mesh elements only
  !---------------------------------------------------
  if (is_root(Mdle)) then
     do i=1,nface(type)
        Neig_list(1:4,i) = ELEMS(Mdle)%neig(i)
     enddo
     return
  endif

  !---------------------------------------------------
  ! Step 1: use neig_face
  !---------------------------------------------------
  call elem_nodes(Mdle, nodesl,norientl) 
  do i=1,nface(type)
     nod = nodesl(nvert(type)+nedge(type)+i)
     call neig_face(nod, nrneig,neig,nsid_list,norient_list)
     select case (nrneig)
     case(1)
        Neig_list(1:4,i) = neig(1)
     case(2)
        ! pick the other one
        if (Mdle.eq.neig(1)) then
           Neig_list(1:4,i) = neig(2)
        else
           Neig_list(1:4,i) = neig(1)
        endif
     end select
  enddo
  
end subroutine find_neig
