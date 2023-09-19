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
   integer :: nod, nsrc, ndst, ndof, ivar
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
   1000 format('copy_coms: Invalid parameter:',A,' = ',I5)
!
   if (Src.eq.Dst) then
      write(*,*) 'copy_coms: Src = Dst = ',Src
      return
   endif
!
!$OMP PARALLEL DO                  &
!$OMP PRIVATE(nsrc,ndst,ndof,ivar) &
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
      if (ndof > 0) then
         do ivar=1,NRHVAR
            NODES(nod)%dof%zdofH(ndst+ivar,1:ndof,N_COMS) = &
            NODES(nod)%dof%zdofH(nsrc+ivar,1:ndof,N_COMS)
         enddo
      endif
  10  continue
!
!  ...H(curl) dof
      if (.not. associated(NODES(nod)%dof%zdofE)) goto 20
      nsrc = (Src-1)*NREVAR
      ndst = (Dst-1)*NREVAR
      ndof = ubound(NODES(nod)%dof%zdofE,2)
      if (ndof > 0) then
         do ivar=1,NREVAR
            NODES(nod)%dof%zdofE(ndst+ivar,1:ndof,N_COMS) = &
            NODES(nod)%dof%zdofE(nsrc+ivar,1:ndof,N_COMS)
         enddo
      endif
  20  continue
!
!  ...H(div) dof
      if (.not. associated(NODES(nod)%dof%zdofV)) goto 30
      nsrc = (Src-1)*NRVVAR
      ndst = (Dst-1)*NRVVAR
      ndof = ubound(NODES(nod)%dof%zdofV,2)
      if (ndof > 0) then
         do ivar=1,NRVVAR
            NODES(nod)%dof%zdofV(ndst+ivar,1:ndof,N_COMS) = &
            NODES(nod)%dof%zdofV(nsrc+ivar,1:ndof,N_COMS)
         enddo
      endif
  30  continue
!
!  ...L2 dof
      if (.not. associated(NODES(nod)%dof%zdofQ)) goto 40
      nsrc = (Src-1)*NRQVAR
      ndst = (Dst-1)*NRQVAR
      ndof = ubound(NODES(nod)%dof%zdofQ,2)
      if (ndof > 0) then
         do ivar=1,NRQVAR
            NODES(nod)%dof%zdofQ(ndst+ivar,1:ndof,N_COMS) = &
            NODES(nod)%dof%zdofQ(nsrc+ivar,1:ndof,N_COMS)
         enddo
      endif
  40  continue
!
!..end of loop through nodes
   enddo
!$OMP END PARALLEL DO
!
end subroutine copy_load
