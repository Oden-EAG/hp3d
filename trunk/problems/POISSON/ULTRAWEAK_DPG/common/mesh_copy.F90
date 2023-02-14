


subroutine Nodes_copy(NODES_cp)

    use control
    use data_structure3D
    use element_data
    use parametersDPG
    use mpi_param, only: RANK

    implicit none



    type(node), intent(out) :: NODES_cp(MAXNODS)
    
    ! integer,    intent(out) :: ELEM_ORDER_cp(NRELES)
    ! integer,    intent(out) :: ELEM_SUBD_cp(NRELES)

    integer :: iel,mdle,NRELES_SUBD_cp,nod
    integer :: NdofH,NdofE,NdofV,NdofQ
    integer :: nn1c,nn2c,ic
    integer :: nnc(2,5), alloc_ind(5)

    do nod = 1,MAXNODS

        NODES_cp(nod)%type = NODES(NOD)%type
        NODES_cp(nod)%case = NODES(nod)%case
        NODES_cp(nod)%order = NODES(nod)%order
        NODES_cp(nod)%bcond = NODES(nod)%bcond
        NODES_cp(nod)%father = NODES(nod)%father
        NODES_cp(nod)%first_son = NODES(nod)%first_son
        NODES_cp(nod)%nr_sons = NODES(nod)%nr_sons
        NODES_cp(nod)%ref_kind = NODES(nod)%ref_kind
        NODES_cp(nod)%geom_interf = NODES(nod)%geom_interf
        NODES_cp(nod)%visit = NODES(nod)%visit
        NODES_cp(nod)%act = NODES(nod)%act
        NODES_cp(nod)%subd = NODES(nod)%subd
        nullify(NODES_cp(nod)%dof)
      !   #if DEBUG_MODE
      !    NODES_cp(nod)%error = NODES(nod)%error
      !   #endif
    enddo


    do nod = 1,MAXNODS
        if(associated(NODES(nod)%dof)) then
            
            allocate(NODES_cp(nod)%dof)
            nullify(NODES_cp(Nod)%dof%coord)
            nullify(NODES_cp(Nod)%dof%zdofH)
            nullify(NODES_cp(Nod)%dof%zdofE)
            nullify(NODES_cp(Nod)%dof%zdofV)
            nullify(NODES_cp(Nod)%dof%zdofQ)

            if(associated(NODES(nod)%dof%coord)) then
                nn1c = ubound(NODES(nod)%dof%coord,1)
                nn2c = ubound(NODES(nod)%dof%coord,2)
                allocate(NODES_cp(nod)%dof%coord(nn1c,nn2c))
                NODES_cp(nod)%dof%coord = ZERO
                
            endif

            if(associated(NODES(nod)%dof%zdofH)) then
                nn1c = ubound(NODES(nod)%dof%zdofH,1)
                nn2c = ubound(NODES(nod)%dof%zdofH,2)
                allocate(NODES_cp(nod)%dof%zdofH(nn1c,nn2c))
                NODES_cp(nod)%dof%zdofH = ZERO
            endif


            if(associated(NODES(nod)%dof%zdofE)) then
                nn1c = ubound(NODES(nod)%dof%zdofE,1)
                nn2c = ubound(NODES(nod)%dof%zdofE,2)
                allocate(NODES_cp(nod)%dof%zdofE(nn1c,nn2c))
                NODES_cp(nod)%dof%zdofE= ZERO
            endif

            if(associated(NODES(nod)%dof%zdofV)) then
                nn1c = ubound(NODES(nod)%dof%zdofV,1)
                nn2c = ubound(NODES(nod)%dof%zdofV,2)
                allocate(NODES_cp(nod)%dof%zdofV(nn1c,nn2c))
                NODES_cp(nod)%dof%zdofV = ZERO  
            endif

            if(associated(NODES(nod)%dof%zdofQ)) then
                nn1c = ubound(NODES(nod)%dof%zdofQ,1)
                nn2c = ubound(NODES(nod)%dof%zdofQ,2)
                allocate(NODES_cp(nod)%dof%zdofQ(nn1c,nn2c))
                NODES_cp(nod)%dof%zdofQ = ZERO
            endif


        endif
    enddo

end subroutine Nodes_copy



subroutine Nodes_replace(NODES_cp)

    use control
    use data_structure3D
    use element_data
    use parametersDPG
    use mpi_param, only: RANK

    implicit none

    integer :: nod
    integer :: nn1c,nn2c
    type(node), intent(out) :: NODES_cp(MAXNODS)

    do nod=1,MAXNODS
        if (associated(NODES(nod)%dof)) then
          if (associated(NODES(nod)%dof%coord)) deallocate(NODES(nod)%dof%coord)
          if (associated(NODES(nod)%dof%zdofH)) deallocate(NODES(nod)%dof%zdofH)
          if (associated(NODES(nod)%dof%zdofE)) deallocate(NODES(nod)%dof%zdofE)
          if (associated(NODES(nod)%dof%zdofV)) deallocate(NODES(nod)%dof%zdofV)
          if (associated(NODES(nod)%dof%zdofQ)) deallocate(NODES(nod)%dof%zdofQ)
          deallocate(NODES(nod)%dof)
        endif
    enddo
    
!    deallocate(NODES)
!
!    if (allocated(ELEM_ORDER)) deallocate(ELEM_ORDER)
!    if (allocated(ELEM_SUBD)) deallocate(ELEM_SUBD)




end subroutine Nodes_replace


subroutine update_ELEM_ORDER_new()
    use control
    use data_structure3D
    use element_data
    use parametersDPG
    use mpi_param, only: RANK

    integer :: iel,mdle
    if (allocated(ELEM_ORDER)) deallocate(ELEM_ORDER)
    if (allocated(ELEM_SUBD))  deallocate(ELEM_SUBD)
    allocate(ELEM_ORDER(NRELES))
    allocate(ELEM_SUBD(NRELES))
    ELEM_ORDER = ZERO
    ELEM_SUBD = ZERO
    mdle = 0; NRELES_SUBD = 0
    do iel=1,NRELES
       call nelcon(mdle, mdle)
       write(*,*) mdle,NODES(mdle)%subd
       ELEM_ORDER(iel) = mdle
       if (NODES(mdle)%subd .eq. RANK) then
          NRELES_SUBD = NRELES_SUBD + 1
          ELEM_SUBD(NRELES_SUBD) = mdle
       endif
    enddo
 end subroutine update_ELEM_ORDER_new
