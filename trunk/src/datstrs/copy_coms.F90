!-------------------------------------------------------------------------------
!> @brief      copies all solution dofs from one component set to another
!!
!> @param[in]  Src: solution component set to copy from = 1,...,NRCOMS
!> @param[in]  Dst: solution component set to copy to   = 1,...,NRCOMS
!!
!> @date       Sep 2023
!-------------------------------------------------------------------------------
subroutine copy_coms(Src,Dst)
!
   use parameters
   use data_structure3D
!
   implicit none
!
   integer, intent(in)  :: Src,Dst
!
   integer :: nod,ndof,nvar
!
!-------------------------------------------------------------------------------
!
!..check consistency
   if ((Src < 1).or.(Src > NRCOMS)) then
      write(*,1000) 'Src',Src
      stop
   endif
   if ((Dst < 1).or.(Dst > NRCOMS)) then
      write(*,1000) 'Dst',Dst
      stop
   endif
   1000 format('copy_coms: Invalid parameter:',A,' = ',I5)
!
   if (Src.eq.Dst) then
      write(*,*) 'copy_coms: Src = Dst = ',Src
      return
   endif
!
!$OMP PARALLEL DO          &
!$OMP PRIVATE(ndof,nvar)   &
!$OMP SCHEDULE(DYNAMIC)
!..loop through active nodes
   do nod=1,NRNODS
      if (Is_inactive(nod)) cycle
      if (.not. associated(NODES(nod)%dof)) cycle
!
!  ...H1 dof
      if (.not. associated(NODES(nod)%dof%zdofH)) goto 10
      nvar = ubound(NODES(nod)%dof%zdofH,1)
      ndof = ubound(NODES(nod)%dof%zdofH,2)
      if (nvar>0 .and. ndof>0) then
         NODES(nod)%dof%zdofH(1:nvar,1:ndof,Dst) = &
         NODES(nod)%dof%zdofH(1:nvar,1:ndof,Src)
      endif
  10  continue
!
!  ...H(curl) dof
      if (.not. associated(NODES(nod)%dof%zdofE)) goto 20
      nvar = ubound(NODES(nod)%dof%zdofE,1)
      ndof = ubound(NODES(nod)%dof%zdofE,2)
      if (nvar>0 .and. ndof>0) then
         NODES(nod)%dof%zdofE(1:nvar,1:ndof,Dst) = &
         NODES(nod)%dof%zdofE(1:nvar,1:ndof,Src)
      endif
  20  continue
!
!  ...H(div) dof
      if (.not. associated(NODES(nod)%dof%zdofV)) goto 30
      nvar = ubound(NODES(nod)%dof%zdofV,1)
      ndof = ubound(NODES(nod)%dof%zdofV,2)
      if (nvar>0 .and. ndof>0) then
         NODES(nod)%dof%zdofV(1:nvar,1:ndof,Dst) = &
         NODES(nod)%dof%zdofV(1:nvar,1:ndof,Src)
      endif
  30  continue
!
!  ...L2 dof
      if (.not. associated(NODES(nod)%dof%zdofQ)) goto 40
      nvar = ubound(NODES(nod)%dof%zdofQ,1)
      ndof = ubound(NODES(nod)%dof%zdofQ,2)
      if (nvar>0 .and. ndof>0) then
         NODES(nod)%dof%zdofQ(1:nvar,1:ndof,Dst) = &
         NODES(nod)%dof%zdofQ(1:nvar,1:ndof,Src)
      endif
  40  continue
!
!..end of loop through nodes
   enddo
!$OMP END PARALLEL DO
!
end subroutine copy_coms
