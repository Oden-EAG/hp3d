!--------------------------------------------------------------------------
!> @brief      find neighbors (up to 4 for an h4 refined face) across faces
!> @param[in]  Mdle      - middle node
!> @param[out] Neig_list - neighbors
!!
!> @date       Feb 2023
!--------------------------------------------------------------------------
!
subroutine find_neig(Mdle, Neig_list)
  use element_data
  use data_structure3D
  implicit none
  ! ** Arguments
  !---------------------------------------------------
  integer, intent(in)  :: Mdle
  integer, intent(out) :: Neig_list(4,6)

  ! ** Locals
  !---------------------------------------------------
  integer, dimension(27) :: nodesl,norientl
  integer, dimension(2)  :: neig,nsid_list,norient_list
  integer :: i,nod,nrneig,ntype
  !
#if HP3D_DEBUG
  integer :: iprint
  iprint=0
#endif
  !---------------------------------------------------

  !  ...initialize
  Neig_list(1:4,1:6)=0
  ntype=NODES(Mdle)%ntype
  !
#if HP3D_DEBUG
  if (iprint.eq.1) then
     write(*,*) 'find_neig: Mdle, type = ', Mdle, S_Type(ntype)
  endif
#endif
  !---------------------------------------------------
  ! Step 0: short cut for initial mesh elements only
  !---------------------------------------------------
  if (is_root(Mdle)) then
     do i=1,nface(ntype)
        Neig_list(1:4,i) = ELEMS(Mdle)%neig(i)
     enddo
     return
  endif

  !---------------------------------------------------
  ! Step 1: use neig_face
  !---------------------------------------------------
  call elem_nodes(Mdle, nodesl,norientl)
  do i=1,nface(ntype)
     nod = nodesl(nvert(ntype)+nedge(ntype)+i)
     call neig_face(nod, nrneig,neig,nsid_list,norient_list)
     select case (nrneig)
     case(1)
        ! Neig_list(1:4,i) = neig(1)
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
