!-------------------------------------------------------------------------------
!> @brief      copies one physics variable component set into another, all loads
!!
!> @param[in]  Src: attribute component set to copy from = 1,...,NR_PHYSA
!> @param[in]  Dst: attribute component set to copy to   = 1,...,NR_PHYSA
!!
!> @date       Sep 2023
!-------------------------------------------------------------------------------
subroutine copy_attr(Src,Dst)
!
   use parameters
   use data_structure3D
!
   implicit none
!
   integer, intent(in)  :: Src,Dst
!
   integer :: nod,src_off,dst_off,src_comp,dst_comp,ncomp,ndof,idof,iload
!
!-------------------------------------------------------------------------------
!
!..check consistency
   if ((Src < 1).or.(Src > NR_PHYSA)) then
      write(*,1000) 'Src',Src
      stop
   endif
   if ((Dst < 1).or.(Dst > NR_PHYSA)) then
      write(*,1000) 'Dst',Dst
      stop
   endif
   1000 format('copy_attr: Invalid parameter:',A,' = ',I5)
!
   if (D_TYPE(Src).ne.D_TYPE(Dst)) then
      write(*,2000) 'D_TYPE(Src)',D_TYPE(Src),'D_TYPE(Dst)',D_TYPE(Dst)
      stop
   endif
   if (NR_COMP(Src).ne.NR_COMP(Dst)) then
      write(*,2000) 'NR_COMP(Src)',NR_COMP(Src),'NR_COMP(Dst)',NR_COMP(Dst)
      stop
   endif
   2000 format('copy_attr: Invalid parameters:',A,' = ',I5,', ',A,' = ',I5)
!
   if (Src.eq.Dst) then
      write(*,*) 'copy_attr: Src = Dst = ',Src
      return
   endif
!
!..first component of Src and Dst variable for the attribute type
   src_off = ADRES(Src)
   dst_off = ADRES(Dst)
!
!..number of attribute components to copy (per load)
   ncomp = NR_COMP(Src)
!
!$OMP PARALLEL DO                                  &
!$OMP PRIVATE(src_comp,dst_comp,ndof,idof,iload)   &
!$OMP SCHEDULE(DYNAMIC)
!..loop through active nodes
   do nod=1,NRNODS
      if (Is_inactive(nod)) cycle
      if (.not. associated(NODES(nod)%dof)) cycle
!
      select case(D_TYPE(Src))
!
      case(CONTIN)
!
!     ...H1 dof
         if (.not. associated(NODES(nod)%dof%zdofH)) goto 10
         ndof = ubound(NODES(nod)%dof%zdofH,2)
         do idof=1,ndof
            do iload=1,NRRHS
               src_comp = (iload-1)*NRHVAR + src_off
               dst_comp = (iload-1)*NRHVAR + dst_off
               NODES(nod)%dof%zdofH(dst_comp+1:dst_comp+ncomp,idof,N_COMS) = &
               NODES(nod)%dof%zdofH(src_comp+1:src_comp+ncomp,idof,N_COMS)
            enddo
         enddo
     10  continue
!
      case(TANGEN)
!
!     ...H(curl) dof
         if (.not. associated(NODES(nod)%dof%zdofE)) goto 20
         ndof = ubound(NODES(nod)%dof%zdofE,2)
         do idof=1,ndof
            do iload=1,NRRHS
               src_comp = (iload-1)*NREVAR + src_off
               dst_comp = (iload-1)*NREVAR + dst_off
               NODES(nod)%dof%zdofE(dst_comp+1:dst_comp+ncomp,idof,N_COMS) = &
               NODES(nod)%dof%zdofE(src_comp+1:src_comp+ncomp,idof,N_COMS)
            enddo
         enddo
     20  continue
!
      case(NORMAL)
!
!     ...H(div) dof
         if (.not. associated(NODES(nod)%dof%zdofV)) goto 30
         ndof = ubound(NODES(nod)%dof%zdofV,2)
         do idof=1,ndof
            do iload=1,NRRHS
               src_comp = (iload-1)*NRVVAR + src_off
               dst_comp = (iload-1)*NRVVAR + dst_off
               NODES(nod)%dof%zdofV(dst_comp+1:dst_comp+ncomp,idof,N_COMS) = &
               NODES(nod)%dof%zdofV(src_comp+1:src_comp+ncomp,idof,N_COMS)
            enddo
         enddo
     30  continue
!
      case(DISCON)
!
!     ...L2 dof
         if (.not. associated(NODES(nod)%dof%zdofQ)) goto 40
         ndof = ubound(NODES(nod)%dof%zdofQ,2)
         do idof=1,ndof
            do iload=1,NRRHS
               src_comp = (iload-1)*NRQVAR + src_off
               dst_comp = (iload-1)*NRQVAR + dst_off
               NODES(nod)%dof%zdofQ(dst_comp+1:dst_comp+ncomp,idof,N_COMS) = &
               NODES(nod)%dof%zdofQ(src_comp+1:src_comp+ncomp,idof,N_COMS)
            enddo
         enddo
     40  continue
!
     end select
!
!..end of loop through nodes
   enddo
!$OMP END PARALLEL DO
!
end subroutine copy_attr
