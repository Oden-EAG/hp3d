!-------------------------------------------------------------------------------
!> @brief      copies one physics variable component into another, all loads
!!
!> @param[in]  SrcAttr: physics attribute to copy from   = 1,...,NR_PHYSA
!> @param[in]  SrcComp: attribute component to copy from = 1,...,NR_COMP(SrcAttr)
!> @param[in]  DstAttr: physics attribute to copy to     = 1,...,NR_PHYSA
!> @param[in]  DstComp: attribute component to copy to   = 1,...,NR_COMP(DstAttr)
!!
!> @date       Sep 2023
!-------------------------------------------------------------------------------
subroutine copy_comp(SrcAttr,SrcComp,DstAttr,DstComp)
!
   use parameters
   use data_structure3D
!
   implicit none
!
   integer, intent(in)  :: SrcAttr,SrcComp
   integer, intent(in)  :: DstAttr,DstComp
!
   integer :: nod,src_comp,dst_comp,ndof
   integer :: idof,iload,icomp_src,icomp_dst
!
!-------------------------------------------------------------------------------
!
!..check consistency
   if ((SrcAttr < 1).or.(SrcAttr > NR_PHYSA)) then
      write(*,1000) 'SrcAttr',SrcAttr
      stop
   endif
   if ((SrcComp < 1).or.(SrcComp > NR_COMP(SrcAttr))) then
      write(*,1000) 'SrcComp',SrcComp
      stop
   endif
   if ((DstAttr < 1).or.(DstAttr > NR_PHYSA)) then
      write(*,1000) 'DstAttr',DstAttr
      stop
   endif
   if ((DstComp < 1).or.(DstComp > NR_COMP(DstAttr))) then
      write(*,1000) 'DstComp',DstComp
      stop
   endif
   1000 format('copy_comp: Invalid parameter:',A,' = ',I5)
!
   if (D_TYPE(SrcAttr).ne.D_TYPE(DstAttr)) then
      write(*,2000) 'D_TYPE(SrcAttr)',D_TYPE(SrcAttr), &
                    'D_TYPE(DstAttr)',D_TYPE(DstAttr)
      stop
   endif
   2000 format('copy_comp: Invalid parameters:',A,' = ',I5,', ',A,' = ',I5)
!
   if ((SrcAttr.eq.DstAttr).and.(SrcComp.eq.DstComp)) then
      write(*,*) 'copy_comp: SrcAttr = DstAttr = ',SrcAttr, &
                          ', SrcComp = DstComp = ',SrcComp
      return
   endif
!
!..component index of Src and Dst attribute component in dof array
   icomp_src = ADRES(SrcAttr) + SrcComp
   icomp_dst = ADRES(DstAttr) + DstComp
!
!$OMP PARALLEL DO                                 &
!$OMP PRIVATE(src_comp,dst_comp,ndof,idof,iload)  &
!$OMP SCHEDULE(DYNAMIC)
!..loop through active nodes
   do nod=1,NRNODS
      if (Is_inactive(nod)) cycle
      if (.not. associated(NODES(nod)%dof)) cycle
!
      select case(D_TYPE(SrcAttr))
!
      case(CONTIN)
!
!     ...H1 dof
         if (.not. associated(NODES(nod)%dof%zdofH)) goto 10
         ndof = ubound(NODES(nod)%dof%zdofH,2)
         do idof=1,ndof
            do iload=1,NRRHS
               src_comp = (iload-1)*NRHVAR + icomp_src
               dst_comp = (iload-1)*NRHVAR + icomp_dst
               NODES(nod)%dof%zdofH(dst_comp,idof,N_COMS) = &
               NODES(nod)%dof%zdofH(src_comp,idof,N_COMS)
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
               src_comp = (iload-1)*NREVAR + icomp_src
               dst_comp = (iload-1)*NREVAR + icomp_dst
               NODES(nod)%dof%zdofE(dst_comp,idof,N_COMS) = &
               NODES(nod)%dof%zdofE(src_comp,idof,N_COMS)
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
               src_comp = (iload-1)*NRVVAR + icomp_src
               dst_comp = (iload-1)*NRVVAR + icomp_dst
               NODES(nod)%dof%zdofV(dst_comp,idof,N_COMS) = &
               NODES(nod)%dof%zdofV(src_comp,idof,N_COMS)
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
               src_comp = (iload-1)*NRQVAR + icomp_src
               dst_comp = (iload-1)*NRQVAR + icomp_dst
               NODES(nod)%dof%zdofQ(dst_comp,idof,N_COMS) = &
               NODES(nod)%dof%zdofQ(src_comp,idof,N_COMS)
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
end subroutine copy_comp
