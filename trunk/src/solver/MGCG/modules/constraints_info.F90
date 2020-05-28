!----------------------------------------------------------------------
!
!   module name        - constraints_info
!
!----------------------------------------------------------------------
!
!   latest revision    - FEB 18
!
!   purpose            - modules sets up the interface for logic_macro
!
!----------------------------------------------------------------------
!
   module constraints_info

   use mg_data_structure, only: GRID
!
!
   CONTAINS
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
   subroutine allocate_constr(Igrid)
!
   use assembly_sc, only : NRDOF_TOT
   implicit none
!
   integer, intent(in) :: Igrid
   integer :: iel, mdle, ndofb, nrdofb, nreles
   integer :: nrnod_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
!
!---------------------------------------------------------------------------------
!
   nreles = GRID(Igrid)%nreles
   allocate(GRID(Igrid)%constr(nreles))
   allocate(GRID(Igrid)%macro_elem(nreles))
!
   nrdofb = 0
!
!$omp parallel default(shared)  &
!$omp private(iel,mdle,nrnod_macro,ndofb,nrdofH_macro,nrdofE_macro,nrdofV_macro)
!$omp do reduction(+:nrdofb) schedule(dynamic)
   do iel = 1, nreles
      mdle = GRID(Igrid)%mdlel(iel)

      call constraints_compute_bounds(Igrid,mdle, nrnod_macro, nrdofH_macro,  &
                                      nrdofE_macro, nrdofV_macro)
!
      GRID(Igrid)%macro_elem(iel)%nrnod_macro = nrnod_macro
!
      allocate(GRID(Igrid)%macro_elem(iel)%nod_macro(nrnod_macro))
      allocate(GRID(Igrid)%macro_elem(iel)%ndofH_macro(nrnod_macro))
      allocate(GRID(Igrid)%macro_elem(iel)%ndofE_macro(nrnod_macro))
      allocate(GRID(Igrid)%macro_elem(iel)%ndofV_macro(nrnod_macro))
!
      call macro_elem_info(Igrid,mdle,                                      &
                           GRID(Igrid)%macro_elem(iel)%ndof_macro, ndofb ,  &
                           GRID(Igrid)%macro_elem(iel)%nrnod_macro,         &
                           GRID(Igrid)%macro_elem(iel)%nod_macro,           &
                           GRID(Igrid)%macro_elem(iel)%ndofH_macro,         &
                           GRID(Igrid)%macro_elem(iel)%ndofE_macro,         &
                           GRID(Igrid)%macro_elem(iel)%ndofV_macro)
!
      nrdofb = nrdofb + ndofb
!
      GRID(Igrid)%constr(iel)%nrnod_macro = nrnod_macro
      GRID(Igrid)%constr(iel)%nrdofH_macro = nrdofH_macro
      GRID(Igrid)%constr(iel)%nrdofE_macro = nrdofE_macro
      GRID(Igrid)%constr(iel)%nrdofV_macro = nrdofV_macro
      allocate(GRID(Igrid)%constr(iel)%nod_macro(nrnod_macro))
      allocate(GRID(Igrid)%constr(iel)%ndofH_macro(nrnod_macro))
      allocate(GRID(Igrid)%constr(iel)%ndofE_macro(nrnod_macro))
      allocate(GRID(Igrid)%constr(iel)%ndofV_macro(nrnod_macro))
!
      GRID(Igrid)%constr(iel)%nod_macro =  GRID(Igrid)%macro_elem(iel)%nod_macro
!
   enddo
!$omp end do
!$omp end parallel

   NRDOF_TOT = nrdofb
!
   end subroutine allocate_constr
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
   subroutine deallocate_constr(Igrid)
!
   implicit none

   integer, intent(in) :: Igrid
   integer :: iel
   do iel = 1, GRID(Igrid)%nreles
!
      deallocate(GRID(Igrid)%macro_elem(iel)%nod_macro)
      deallocate(GRID(Igrid)%macro_elem(iel)%ndofH_macro)
      deallocate(GRID(Igrid)%macro_elem(iel)%ndofE_macro)
      deallocate(GRID(Igrid)%macro_elem(iel)%ndofV_macro)
!
      deallocate(GRID(Igrid)%constr(iel)%nacH_macro)
      deallocate(GRID(Igrid)%constr(iel)%nacE_macro)
      deallocate(GRID(Igrid)%constr(iel)%nacV_macro)
!
      deallocate(GRID(Igrid)%constr(iel)%nrconH_macro)
      deallocate(GRID(Igrid)%constr(iel)%nrconE_macro)
      deallocate(GRID(Igrid)%constr(iel)%nrconV_macro)
!
      deallocate(GRID(Igrid)%constr(iel)%constrH_macro)
      deallocate(GRID(Igrid)%constr(iel)%constrE_macro)
      deallocate(GRID(Igrid)%constr(iel)%constrV_macro)
!
      deallocate(GRID(Igrid)%constr(iel)%nod_macro)
      deallocate(GRID(Igrid)%constr(iel)%ndofH_macro)
      deallocate(GRID(Igrid)%constr(iel)%ndofE_macro)
      deallocate(GRID(Igrid)%constr(iel)%ndofV_macro)
!
   enddo
   deallocate(GRID(Igrid)%constr)
   deallocate(GRID(Igrid)%macro_elem)
! !
   end subroutine deallocate_constr


   subroutine store_constraints(Igrid)
!
   use parameters, only: NACDIM
   implicit none
!
   integer, intent(in) :: Igrid
   integer :: iel, mdle, nrdofH_macro, nrdofE_macro, nrdofV_macro
!
!---------------------------------------------------------------------------------
!
!$omp parallel default(shared)    &
!$omp private(iel,mdle,nrdofH_macro,nrdofE_macro,nrdofV_macro)
!$omp do schedule(guided)
   do iel = 1, GRID(Igrid)%nreles
      mdle = GRID(Igrid)%mdlel(iel)

      nrdofH_macro = GRID(Igrid)%constr(iel)%nrdofH_macro
      nrdofE_macro = GRID(Igrid)%constr(iel)%nrdofE_macro
      nrdofV_macro = GRID(Igrid)%constr(iel)%nrdofV_macro

      allocate(GRID(Igrid)%constr(iel)%constrH_macro(NACDIM,nrdofH_macro))
      allocate(GRID(Igrid)%constr(iel)%constrE_macro(NACDIM,nrdofE_macro))
      allocate(GRID(Igrid)%constr(iel)%constrV_macro(NACDIM,nrdofV_macro))
!
      allocate(GRID(Igrid)%constr(iel)%nacH_macro(NACDIM,nrdofH_macro))
      allocate(GRID(Igrid)%constr(iel)%nacE_macro(NACDIM,nrdofE_macro))
      allocate(GRID(Igrid)%constr(iel)%nacV_macro(NACDIM,nrdofV_macro))
!
      allocate(GRID(Igrid)%constr(iel)%nrconH_macro(nrdofH_macro))
      allocate(GRID(Igrid)%constr(iel)%nrconE_macro(nrdofE_macro))
      allocate(GRID(Igrid)%constr(iel)%nrconV_macro(nrdofV_macro))
!
      call logic_macro(Igrid, mdle,                                                      &
           GRID(Igrid)%constr(iel)%nod_macro,     GRID(Igrid)%constr(iel)%nrnod_macro,   &
           nrdofH_macro, nrdofE_macro, nrdofV_macro,                                     &
           GRID(Igrid)%constr(iel)%nrnodm,        GRID(Igrid)%constr(iel)%nodm,          &
           GRID(Igrid)%constr(iel)%ndofmH,        GRID(Igrid)%constr(iel)%ndofmE,        &
           GRID(Igrid)%constr(iel)%ndofmV,        GRID(Igrid)%constr(iel)%ndofH_macro,   &
           GRID(Igrid)%constr(iel)%ndofE_macro,   GRID(Igrid)%constr(iel)%ndofV_macro,   &
           GRID(Igrid)%constr(iel)%nrconH_macro,  GRID(Igrid)%constr(iel)%nacH_macro,    &
           GRID(Igrid)%constr(iel)%constrH_macro, GRID(Igrid)%constr(iel)%nrconE_macro,  &
           GRID(Igrid)%constr(iel)%nacE_macro,    GRID(Igrid)%constr(iel)%constrE_macro, &
           GRID(Igrid)%constr(iel)%nrconV_macro,  GRID(Igrid)%constr(iel)%nacV_macro,    &
           GRID(Igrid)%constr(iel)%constrV_macro)
!
   enddo
!$omp end do
!$omp end parallel
!
!
   end subroutine store_constraints
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
! -----------------------------------------------------------------------
!
!    routine name      - constraints_compute_bounds
!
! -----------------------------------------------------------------------
!
!    latest revision   - Mar 2018
!
!    purpose           - returns Nrnod_macro, NrdofH_macro,
!                        NrdofE_macro, NrdofV_macro
!
!   arguments :
!     in:
!             MdleC    - middle node of a coarse element
!     out:
!          Nrnod_macro - number of nodes
!         NrdofH_macro - number of H1      dof for the macro element
!         NrdofE_macro - number of H(curl) dof for the macro element
!         NrdofV_macro - number of H(div)  dof for the macro element
!
! ----------------------------------------------------------------------
!
   subroutine constraints_compute_bounds(Igrid,MdleC, Nrnod_macro, NrdofH_macro, &
                                         NrdofE_macro, NrdofV_macro)
!
   use data_structure3D,  only: NODES,MAXNODM, ndof_nod, NRNODS
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer, intent(in)  :: Igrid,MdleC
   integer, intent(out) :: Nrnod_macro,  NrdofH_macro
   integer, intent(out) :: NrdofE_macro, NrdofV_macro
!
!..locals
!..work space for elem_nodes and logic_nodes
   integer              :: nodesl(27),norientl(27), nodm(MAXNODM)
!
   integer              :: mdle, imdle, i, nod, master
   integer              :: nrnodm, nH, nE, nV, nQ
   integer, allocatable :: nvisit(:)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   allocate(nvisit(NRNODS)); nvisit = 0
!
!..Step 1: determine the first dof offsets for active nodes
   NrdofH_macro  = 0 ; NrdofE_macro  = 0; NrdofV_macro  = 0
   mdle=0 ;   Nrnod_macro = 0 ;   imdle = 1
   call nelcon_macro(MdleC,mdle, mdle)
   if (mdle .ne. 0) then
      imdle = 0
   endif
   do while (mdle.ne.0)
      imdle = imdle +1
      call get_connect_info(mdle, nodesl,norientl)
!
!   ..get nodes of the modified coarse element
      call logic_nodes(mdle, nodesl, nodm, nrnodm )
!
      do i = 1,nrnodm-1
         nod = nodm(i)
         if (nvisit(nod) .eq. 1) cycle
         call find_master(Igrid,nod,master)
         select case(NODES(master)%type)
         case('mdlb','mdln','mdlp','mdld')
         case default
            Nrnod_macro = Nrnod_macro + 1
            call ndof_nod(NODES(nod)%type,NODES(nod)%order,nH,nE,nV,nQ)
            NrdofH_macro=NrdofH_macro+nH
            NrdofE_macro=NrdofE_macro+nE
            NrdofV_macro=NrdofV_macro+nV
            nvisit(nod) = 1
         end select
      enddo
!
      call nelcon_macro(MdleC,mdle, mdle)
   enddo
!
   deallocate(nvisit)
!
   end subroutine constraints_compute_bounds
! !
! !
   end module constraints_info




