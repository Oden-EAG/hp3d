


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


    ! mdle = 0;NRELES_SUBD_cp = 0
    ! do iel = 1,NRELES
    !     call nelcon(mdle,mdle)
    !     ELEM_ORDER_cp(iel) = mdle
    !     if (NODES(mdle)%subd .eq. RANK) then
    !         NRELES_SUBD_cp = NRELES_SUBD_cp + 1
    !         ELEM_SUBD_cp(NRELES_SUBD_cp) = mdle
    !     endif
    ! enddo

    !copying the old nodes array in the new copy array
    ! NODES_cp(1:MAXNODS) = NODES(1:MAXNODS)

    ! !checking whether dofs are allocated for NODES array
    ! !and if yes then allocating the smae amount in NODES_cp
    ! nnc = ZERO
    ! do nod = 1,NRNODS

    !     if(associated(NODES(nod)%dof)) then
    !         call find_ndof(nod, NdofH,NdofE,NdofV,NdofQ)

    !         if (associated(NODES(nod)%dof%coord)) then
    !             ! allocate(NODES_cp(nod)%dof%coord(3,NdofH))
    !             ! write(*,*) 1
    !             nn1c = ubound(NODES(nod)%dof%coord,1)
    !             nn2c = ubound(NODES(nod)%dof%coord,2)
    !             alloc_ind(1) = 1;nnc(1,5) = nn1c;nnc(2,5) = nn2c
    !         endif

    !         if(associated(NODES(nod)%dof%zdofH)) then
    !             nn1c = ubound(NODES(nod)%dof%zdofH,1)
    !             nn2c = ubound(NODES(nod)%dof%zdofH,2)
                
    !             alloc_ind(2) = 1;nnc(1,1) = nn1c;nnc(2,1) = nn2c
                
    !             ! allocate(NODES_cp(nod)%dof%zdofH(nn1c,nn2c))
    !         endif


    !         if(associated(NODES(nod)%dof%zdofE)) then
    !             nn1c = ubound(NODES(nod)%dof%zdofE,1)
    !             nn2c = ubound(NODES(nod)%dof%zdofE,2)

    !             alloc_ind(3) = 1;nnc(1,2) = nn1c;nnc(2,2) = nn2c
    !             ! allocate(NODES_cp(nod)%dof%zdofE(nn1c,nn2c))
    !         endif


    !         if(associated(NODES(nod)%dof%zdofV)) then
    !             nn1c = ubound(NODES(nod)%dof%zdofV,1)
    !             nn2c = ubound(NODES(nod)%dof%zdofV,2)

    !             alloc_ind(4) = 1;nnc(1,3) = nn1c; nnc(2,3) = nn2c
    !             ! allocate(NODES_cp(nod)%dof%zdofV(nn1c,nn2c))
    !         endif


    !         if(associated(NODES(nod)%dof%zdofQ)) then
    !             nn1c = ubound(NODES(nod)%dof%zdofQ,1)
    !             nn2c = ubound(NODES(nod)%dof%zdofQ,2)

    !             alloc_ind(5) = 1; nnc(1,4) = nn1c; nnc(2,4) = nn2c
    !             ! allocate(NODES_cp(nod)%dof%zdofQ(nn1c,nn2c))
    !         endif

    !         ! allocating nodes to NODES_cp
    !         ! tomorrow fix this
    !         if(alloc_ind(1) .eq. 1) then
    !             ! nullify(NODES_cp(nod)%dof%coord)
    !             ! allocate(NODES_cp(nod)%dof%coord(nnc(1,5),nnc(2,5)))
    !             write(*,*) (NODES_cp(nod)%dof%coord .eq. NODES(nod)%dof%coord)
    !             ! write(*,*) (NODES_cp(nod)%dof%coord(1,1)), (NODES(nod)%dof%coord(1,1))
            
    !         endif

    !         ! if(alloc_ind(2) .eq. 1) then

    !         !     allocate(NODES_cp(nod)%dof%zdofH(nnc(1,1),nnc(2,1)))
    !         !     ! write(*,*) associated(NODES_cp(nod)%dof%zdofH)
    !         !     ! write(*,*) (NODES_cp(nod)%dof%zdofH(1,1)), (NODES(nod)%dof%zdofH(1,1))
    !         ! endif

    !         ! if(alloc_ind(3) .eq. 1) then

    !         !     allocate(NODES_cp(nod)%dof%zdofE(nnc(1,2),nnc(2,2)))

    !         ! endif

    !         ! if(alloc_ind(4) .eq. 1) then

    !         !     allocate(NODES_cp(nod)%dof%zdofV(nnc(1,3),nnc(2,3)))

    !         ! endif

    !         ! if(alloc_ind(5) .eq. 1) then

    !         !     allocate(NODES_cp(nod)%dof%zdofQ(nnc(1,4),nnc(2,4)))

    !         ! endif

    !     endif

    ! enddo

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



subroutine Nodes_replace()

    use control
    use data_structure3D
    use element_data
    use parametersDPG
    use mpi_param, only: RANK

    implicit none

    integer :: nod


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
    
    deallocate(NODES)
!
    if (allocated(ELEM_ORDER)) deallocate(ELEM_ORDER)
    if (allocated(ELEM_SUBD)) deallocate(ELEM_SUBD)
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
       ELEM_ORDER(iel) = mdle
       if (NODES(mdle)%subd .eq. RANK) then
          NRELES_SUBD = NRELES_SUBD + 1
          ELEM_SUBD(NRELES_SUBD) = mdle
       endif
    enddo
 end subroutine