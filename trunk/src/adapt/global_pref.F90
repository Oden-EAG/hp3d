!--------------------------------------------------------------------
!> @brief   routine performs a global p-enrichment
!> @date    Feb 2023
!--------------------------------------------------------------------
subroutine global_pref
!
   use parameters      , only: MAXP
   use data_structure3D, only: NRELES,NODES,ELEM_ORDER
   use node_types
!
   implicit none
!
   integer :: mdle,p,iel,nord,nordx,nordy,nordz,naux
!
!$OMP PARALLEL DO  &
!$OMP PRIVATE(mdle,nord,nordx,nordy,nordz,naux,p)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      nord = NODES(mdle)%order
      select case(NODES(mdle)%ntype)
         case(MDLB)
            nord = nord+111
            call decode(nord, naux ,nordz)
            call decode(naux, nordx,nordy)
            p = MAX(nordx,nordy,nordz)
         case(MDLP)
            nord = nord+11
            call decode(nord, nordx,nordz)
            p = MAX(nordx,nordz)
         case(MDLN,MDLD)
            nord = nord+1
            p = nord
      end select
      if (p .gt. MAXP) then
         write(*,100) 'global_pref: mdle,p,MAXP = ',mdle,p,MAXP,'. stop.'
         stop
     100 format(A,I7,I3,I3,A)
      endif
!
      call nodmod(mdle,nord)
   enddo
!$OMP END PARALLEL DO
!
!..raise order of approximation on non-middle nodes by enforcing minimum rule
   call enforce_min_rule
!
end subroutine global_pref
!
!
!--------------------------------------------------------------------
!> @brief   routine performs a global p-unrefinement
!> @date    Feb 2023
!--------------------------------------------------------------------
subroutine global_punref
!
   use data_structure3D, only: NRELES,NODES,ELEM_ORDER
   use par_mesh        , only: DISTRIBUTED
   use node_types
!
   implicit none
!
   integer :: mdle,iel,nord
!
   if (DISTRIBUTED) then
      write(*,*) 'global_punref: not implemented for distributed mesh. stop.'
      stop
   endif
!
!..loop over active elements
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      nord = NODES(mdle)%order
      select case(NODES(mdle)%ntype)
         case(MDLB); nord = nord-111
         case(MDLP); nord = nord-11
         case(MDLN,MDLD); nord = nord-1
      end select
      call nodmod(mdle, nord)
   enddo
!
end subroutine global_punref
