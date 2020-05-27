!--------------------------------------------------------------------
!> Purpose : routine performs a global p-enrichment
!!
!> @date Aug 2019
!--------------------------------------------------------------------
subroutine global_pref
!
   use parameters      , only: MAXP
   use data_structure3D, only: NRELES,NRNODS,NODES,MAXNODM, &
                               ELEM_ORDER,get_subd,set_subd
   use par_mesh        , only: DISTRIBUTED,get_elem_nodes
   use mpi_param       , only: RANK
!
   implicit none
!
   integer :: mdle,p,i,iel,nord,nordx,nordy,nordz,naux,nrnodm,subd
   integer :: nodm(MAXNODM)
!
!
!$OMP PARALLEL DO  &
!$OMP PRIVATE(mdle,nord,nordx,nordy,nordz,naux,p)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      nord = NODES(mdle)%order
      select case(NODES(mdle)%type)
         case('mdln','mdld')
            nord = nord+1
            p = nord
         case('mdlp')
            nord = nord+11
            call decode(nord, nordx,nordz)
            p = MAX(nordx,nordz)
         case('mdlb')
            nord = nord+111
            call decode(nord, naux ,nordz)
            call decode(naux, nordx,nordy)
            p = MAX(nordx,nordy,nordz)
      end select
      if (p .gt. MAXP) then
         write(*,100) 'global_pref: mdle,p,MAXP = ',mdle,p,MAXP,'. stop.'
         stop
     100 format(A,I7,I3,I3,A)
      endif
!     -------- begin setting subdomain for distributed mesh
      !..should not be necessary anymore after recent changes in distr_mesh,
      !  because subd values should already be set correctly within subdomain
      !  so that nodmod will work fine.
!      call get_subd(mdle, subd)
!      if ((DISTRIBUTED) .and. (subd .eq. RANK)) then
!         call get_elem_nodes(mdle, nodm,nrnodm)
!         do i=1,nrnodm
!            call set_subd(nodm(i),subd)
!         enddo
!      endif
!     -------- end setting subdomain for distributed mesh
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
subroutine global_punref
!
   use data_structure3D, only: NRELES,NRNODS,NODES,ELEM_ORDER
   use par_mesh        , only: DISTRIBUTED
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
      select case(NODES(mdle)%type)
         case('mdln','mdld'); nord = nord-1
         case('mdlp'); nord = nord-11
         case('mdlb'); nord = nord-111
      end select
      call nodmod(mdle, nord)
   enddo
!
end subroutine global_punref
