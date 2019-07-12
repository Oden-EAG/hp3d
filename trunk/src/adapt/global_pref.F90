!--------------------------------------------------------------------
!> Purpose : routine performs a global p-enrichment
!!
!> @date July 2019
!--------------------------------------------------------------------
subroutine global_pref
!
   use data_structure3D, only: NRELES,NRNODS,NODES,MAXNODM, &
                               get_subd,set_subd
   use par_mesh        , only: DISTRIBUTED,get_elem_nodes
   use MPI_param       , only: RANK
!      
   implicit none
!
   integer :: mdle,i,iel,nord,nrnodm,subd
   integer :: nodm(MAXNODM)
!
!..loop over active elements (middle nodes)
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      nord = NODES(mdle)%order
      select case(NODES(mdle)%type)
         case('mdln','mdld'); nord = nord+1
         case('mdlp'); nord = nord+11
         case('mdlb'); nord = nord+111
      end select
!     -------- begin setting subdomain for distributed mesh
      !..should not be necessary anymore with after recent changes in distr_mesh,
      !  because subd values should already be set correctly within subdomain
      !  so that nodmod will work fine.
!      call get_subd(mdle, subd)
!      if ((DISTRIBUTED .eq. 1) .and. (subd .eq. RANK)) then
!         call get_elem_nodes(mdle, nodm,nrnodm)
!         do i=1,nrnodm
!            call set_subd(nodm(i),subd)
!         enddo
!      endif
!     -------- end setting subdomain for distributed mesh
      call nodmod(mdle, nord)
   enddo
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
   use data_structure3D, only: NRELES,NRNODS,NODES
   use par_mesh        , only: DISTRIBUTED
!      
   implicit none
!
   integer :: mdle,iel,nord
!
   if (DISTRIBUTED .eq. 1) then
      write(*,*) 'global_punref: not implemented for distributed mesh. stop.'
      stop
   endif
!
!..loop over active elements
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
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
