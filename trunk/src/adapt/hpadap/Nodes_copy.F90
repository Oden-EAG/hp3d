!
!-----------------------------------------------------------------------
!> @brief Routine copies the current Nodes array into NODES_cp
!> @param[in]   NODES_cp - Secondary nodes data structure required to
!!                         store the current mesh, we don't copy the
!!                         associated with the current mesh
!> @date May 2024
!-----------------------------------------------------------------------
subroutine Nodes_copy(NODES_cp)
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use mpi_param, only: RANK
!
   implicit none
!
   type(node), intent(out) :: NODES_cp(MAXNODS)
!
   integer :: iel,mdle,NRELES_SUBD_cp,nod
   integer :: NdofH,NdofE,NdofV,NdofQ
   integer :: nn1c,nn2c,nn3c,ic
   integer :: nnc(2,5), alloc_ind(5)
!
   do nod = 1,MAXNODS
      NODES_cp(nod)%ntype = NODES(NOD)%ntype
      NODES_cp(nod)%case = NODES(nod)%case
      NODES_cp(nod)%order = NODES(nod)%order
      NODES_cp(nod)%bcond = NODES(nod)%bcond
      NODES_cp(nod)%father = NODES(nod)%father
      NODES_cp(nod)%first_son = NODES(nod)%first_son
      NODES_cp(nod)%nr_sons = NODES(nod)%nr_sons
      
      NODES_cp(nod)%ref_kind = NODES(nod)%ref_kind
      NODES_cp(nod)%visit = NODES(nod)%visit
      NODES_cp(nod)%act = NODES(nod)%act
      NODES_cp(nod)%subd = NODES(nod)%subd
      nullify(NODES_cp(nod)%dof)
   enddo
!
   do nod = 1,MAXNODS
      if(associated(NODES(nod)%dof)) then
!           
         allocate(NODES_cp(nod)%dof)
         nullify(NODES_cp(Nod)%dof%coord)
         nullify(NODES_cp(Nod)%dof%zdofH)
         nullify(NODES_cp(Nod)%dof%zdofE)
         nullify(NODES_cp(Nod)%dof%zdofV)
         nullify(NODES_cp(Nod)%dof%zdofQ)
!
         if(associated(NODES(nod)%dof%coord)) then
               nn1c = ubound(NODES(nod)%dof%coord,1)
               nn2c = ubound(NODES(nod)%dof%coord,2)
               allocate(NODES_cp(nod)%dof%coord(nn1c,nn2c))
               NODES_cp(nod)%dof%coord = ZERO
         endif
!
         if(associated(NODES(nod)%dof%zdofH)) then
               nn1c = ubound(NODES(nod)%dof%zdofH,1)
               nn2c = ubound(NODES(nod)%dof%zdofH,2)
               nn3c = ubound(NODES(nod)%dof%zdofH,3)
               allocate(NODES_cp(nod)%dof%zdofH(nn1c,nn2c,nn3c))
               NODES_cp(nod)%dof%zdofH = ZERO
         endif
!
         if(associated(NODES(nod)%dof%zdofE)) then
               nn1c = ubound(NODES(nod)%dof%zdofE,1)
               nn2c = ubound(NODES(nod)%dof%zdofE,2)
               nn3c = ubound(NODES(nod)%dof%zdofE,3)
               allocate(NODES_cp(nod)%dof%zdofE(nn1c,nn2c,nn3c))
               NODES_cp(nod)%dof%zdofE= ZERO
         endif
!
         if(associated(NODES(nod)%dof%zdofV)) then
               nn1c = ubound(NODES(nod)%dof%zdofV,1)
               nn2c = ubound(NODES(nod)%dof%zdofV,2)
               nn3c = ubound(NODES(nod)%dof%zdofV,3)
               allocate(NODES_cp(nod)%dof%zdofV(nn1c,nn2c,nn3c))
               NODES_cp(nod)%dof%zdofV = ZERO  
         endif
!
         if(associated(NODES(nod)%dof%zdofQ)) then
               nn1c = ubound(NODES(nod)%dof%zdofQ,1)
               nn2c = ubound(NODES(nod)%dof%zdofQ,2)
               nn3c = ubound(NODES(nod)%dof%zdofQ,3)
               allocate(NODES_cp(nod)%dof%zdofQ(nn1c,nn2c,nn3c))
               NODES_cp(nod)%dof%zdofQ = ZERO
         endif
!
      endif
   enddo
end subroutine Nodes_copy
!
!-----------------------------------------------------------------------
!> @brief Routine deallocates the copy of Nodes array
!> @param[in]   NODES_cp - Secondary nodes data structure required to
!!                         store the current mesh, we don't copy the
!!                         associated with the current mesh
!> @date May 2024
!-----------------------------------------------------------------------
subroutine Nodes_dealloc(NODES_cp)
!
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use mpi_param, only: RANK
!
   implicit none
!
   integer :: nod
   integer :: nn1c,nn2c
   type(node), intent(out) :: NODES_cp(MAXNODS)
!
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
!
end subroutine Nodes_dealloc

