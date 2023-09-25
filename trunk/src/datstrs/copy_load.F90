!-------------------------------------------------------------------------------
!> @brief      copies solution dofs from one RHS solution set to another
!!
!> @param[in]  Src: RHS solution set to copy from = 1,...,NRRHS
!> @param[in]  Dst: RHS solution set to copy to   = 1,...,NRRHS
!!
!> @date       Sep 2023
!-------------------------------------------------------------------------------
subroutine copy_load(Src,Dst)
!
   use parameters
   use data_structure3D
!
   implicit none
!
   integer, intent(in)  :: Src,Dst
!
   integer :: nod,nsrc,ndst,ndof,idof
!
!-------------------------------------------------------------------------------
!
!..check consistency
   if ((Src < 1).or.(Src > NRRHS)) then
      write(*,1000) 'Src',Src
      stop
   endif
   if ((Dst < 1).or.(Dst > NRRHS)) then
      write(*,1000) 'Dst',Dst
      stop
   endif
   1000 format('copy_load: Invalid parameter:',A,' = ',I5)
!
   if (Src.eq.Dst) then
      write(*,*) 'copy_load: Src = Dst = ',Src
      return
   endif
!
!$OMP PARALLEL DO                  &
!$OMP PRIVATE(nsrc,ndst,ndof,idof) &
!$OMP SCHEDULE(DYNAMIC)
!..loop through active nodes
   do nod=1,NRNODS
      if (Is_inactive(nod)) cycle
      if (.not. associated(NODES(nod)%dof)) cycle
!
!  ...H1 dof
      if (.not. associated(NODES(nod)%dof%zdofH)) goto 10
      nsrc = (Src-1)*NRHVAR
      ndst = (Dst-1)*NRHVAR
      ndof = ubound(NODES(nod)%dof%zdofH,2)
      do idof=1,ndof
         NODES(nod)%dof%zdofH(ndst+1:ndst+NRHVAR,idof,N_COMS) = &
         NODES(nod)%dof%zdofH(nsrc+1:nsrc+NRHVAR,idof,N_COMS)
      enddo
  10  continue
!
!  ...H(curl) dof
      if (.not. associated(NODES(nod)%dof%zdofE)) goto 20
      nsrc = (Src-1)*NREVAR
      ndst = (Dst-1)*NREVAR
      ndof = ubound(NODES(nod)%dof%zdofE,2)
      do idof=1,ndof
         NODES(nod)%dof%zdofE(ndst+1:ndst+NREVAR,idof,N_COMS) = &
         NODES(nod)%dof%zdofE(nsrc+1:nsrc+NREVAR,idof,N_COMS)
      enddo
  20  continue
!
!  ...H(div) dof
      if (.not. associated(NODES(nod)%dof%zdofV)) goto 30
      nsrc = (Src-1)*NRVVAR
      ndst = (Dst-1)*NRVVAR
      ndof = ubound(NODES(nod)%dof%zdofV,2)
      do idof=1,ndof
         NODES(nod)%dof%zdofV(ndst+1:ndst+NRVVAR,idof,N_COMS) = &
         NODES(nod)%dof%zdofV(nsrc+1:nsrc+NRVVAR,idof,N_COMS)
      enddo
  30  continue
!
!  ...L2 dof
      if (.not. associated(NODES(nod)%dof%zdofQ)) goto 40
      nsrc = (Src-1)*NRQVAR
      ndst = (Dst-1)*NRQVAR
      ndof = ubound(NODES(nod)%dof%zdofQ,2)
      do idof=1,ndof
         NODES(nod)%dof%zdofQ(ndst+1:ndst+NRQVAR,idof,N_COMS) = &
         NODES(nod)%dof%zdofQ(nsrc+1:nsrc+NRQVAR,idof,N_COMS)
      enddo
  40  continue
!
!..end of loop through nodes
   enddo
!$OMP END PARALLEL DO
!
end subroutine copy_load
