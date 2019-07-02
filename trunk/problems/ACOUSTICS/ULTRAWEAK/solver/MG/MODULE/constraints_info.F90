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

   use macro_grid_info
   use assembly 


   type constraints
!..modified coarse grid element
   integer              :: nrnodm, nodm(MAXNODM), nrdof
   integer              :: ndofmH(MAXNODM), ndofmE(MAXNODM), ndofmV(MAXNODM)
!..macro-element
   integer              :: nrnod_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
   integer, allocatable :: nod_macro(:)
   integer, allocatable :: ndofH_macro(:),  ndofE_macro(:),  ndofV_macro(:)
   integer, allocatable :: nrconH_macro(:), nrconE_macro(:), nrconV_macro(:)
   integer, allocatable :: nacH_macro(:,:), nacE_macro(:,:), nacV_macro(:,:)
   real*8,  allocatable :: constrH_macro(:,:),constrE_macro(:,:), constrV_macro(:,:)

   end type constraints

   type (constraints), allocatable :: CONSTR(:)


! 
   CONTAINS
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
   subroutine allocate_constr
!
   use m_assembly, only : NRDOF_CON, NRDOF_TOT
   implicit none
!   
   integer :: iel, mdle, ndofb, nrdofb, nrdof
   integer :: nrnod_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
!
!---------------------------------------------------------------------------------
!
   allocate(CONSTR(NRELES_COARSE))
   nrdof = 0; nrdofb = 0
!
!$omp parallel default(shared)  &
!$omp private(iel,mdle,nrnod_macro,ndofb,nrdofH_macro,nrdofE_macro,nrdofV_macro)
!$omp do reduction(+:nrdof, nrdofb) schedule(dynamic)
   do iel = 1, NRELES_COARSE
      mdle = MDLE_MACRO(iel) 

      call constraints_compute_bounds(mdle, nrnod_macro, nrdofH_macro,  &
                                      nrdofE_macro, nrdofV_macro)
!
      MACRO_ELEM(iel)%nrnod_macro = nrnod_macro
!      
      allocate(MACRO_ELEM(iel)%nod_macro(nrnod_macro))
      allocate(MACRO_ELEM(iel)%ndofH_macro(nrnod_macro))
      allocate(MACRO_ELEM(iel)%ndofE_macro(nrnod_macro))
      allocate(MACRO_ELEM(iel)%ndofV_macro(nrnod_macro))
!      
      call macro_elem_info(mdle,MACRO_ELEM(iel)%ndof_macro, ndofb ,MACRO_ELEM(iel)%nrnod_macro, &
                                MACRO_ELEM(iel)%nod_macro,  MACRO_ELEM(iel)%ndofH_macro, &
                                MACRO_ELEM(iel)%ndofE_macro,MACRO_ELEM(iel)%ndofV_macro)

      nrdofb = nrdofb + ndofb
      nrdof  = nrdof  + MACRO_ELEM(iel)%ndof_macro
!      
      CONSTR(iel)%nrnod_macro = nrnod_macro
      CONSTR(iel)%nrdofH_macro = nrdofH_macro
      CONSTR(iel)%nrdofE_macro = nrdofE_macro
      CONSTR(iel)%nrdofV_macro = nrdofV_macro
      allocate(CONSTR(iel)%nod_macro(nrnod_macro))   
      allocate(CONSTR(iel)%ndofH_macro(nrnod_macro)) 
      allocate(CONSTR(iel)%ndofE_macro(nrnod_macro)) 
      allocate(CONSTR(iel)%ndofV_macro(nrnod_macro)) 
!      
      CONSTR(iel)%nod_macro(1:nrnod_macro) = MACRO_ELEM(iel)%nod_macro(1:nrnod_macro)
!
   enddo   
!$omp end do
!$omp end parallel

   NRDOF_CON = nrdof
   NRDOF_TOT = NRDOF_CON + nrdofb
!
   end subroutine allocate_constr
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!
   subroutine deallocate_constr
!
   implicit none

   integer :: iel
   do iel = 1, NRELES_COARSE
!    
      deallocate(MACRO_ELEM(iel)%nod_macro)
      deallocate(MACRO_ELEM(iel)%ndofH_macro)
      deallocate(MACRO_ELEM(iel)%ndofE_macro)
      deallocate(MACRO_ELEM(iel)%ndofV_macro)
!
      deallocate(CONSTR(iel)%nacH_macro)
      deallocate(CONSTR(iel)%nacE_macro)
      deallocate(CONSTR(iel)%nacV_macro)
!
      deallocate(CONSTR(iel)%nrconH_macro)
      deallocate(CONSTR(iel)%nrconE_macro)
      deallocate(CONSTR(iel)%nrconV_macro)
!
      deallocate(CONSTR(iel)%constrH_macro)
      deallocate(CONSTR(iel)%constrE_macro)
      deallocate(CONSTR(iel)%constrV_macro)
!
      deallocate(CONSTR(iel)%nod_macro)   
      deallocate(CONSTR(iel)%ndofH_macro) 
      deallocate(CONSTR(iel)%ndofE_macro) 
      deallocate(CONSTR(iel)%ndofV_macro) 
!
   enddo   
   deallocate(CONSTR)
!
   end subroutine deallocate_constr


   subroutine store_constraints
!
   implicit none
!
   integer :: iel, mdle, nrdofH_macro, nrdofE_macro, nrdofV_macro
!
!---------------------------------------------------------------------------------
!
!$omp parallel default(shared)    &
!$omp private(iel,mdle,nrdofH_macro,nrdofE_macro,nrdofV_macro)
!$omp do schedule(dynamic)
   do iel = 1, NRELES_COARSE
      mdle = MDLE_MACRO(iel) 

      nrdofH_macro = CONSTR(iel)%nrdofH_macro 
      nrdofE_macro = CONSTR(iel)%nrdofE_macro 
      nrdofV_macro = CONSTR(iel)%nrdofV_macro 

      allocate(CONSTR(iel)%constrH_macro(NACDIM,nrdofH_macro))
      allocate(CONSTR(iel)%constrE_macro(NACDIM,nrdofE_macro))
      allocate(CONSTR(iel)%constrV_macro(NACDIM,nrdofV_macro))
!
      allocate(CONSTR(iel)%nacH_macro(NACDIM,nrdofH_macro))
      allocate(CONSTR(iel)%nacE_macro(NACDIM,nrdofE_macro))
      allocate(CONSTR(iel)%nacV_macro(NACDIM,nrdofV_macro))
!
      allocate(CONSTR(iel)%nrconH_macro(nrdofH_macro))
      allocate(CONSTR(iel)%nrconE_macro(nrdofE_macro))
      allocate(CONSTR(iel)%nrconV_macro(nrdofV_macro))
!        
      call logic_macro(mdle, CONSTR(iel)%nod_macro, CONSTR(iel)%nrnod_macro,        &
                       nrdofH_macro, nrdofE_macro, nrdofV_macro,                    &
                       CONSTR(iel)%nrnodm, CONSTR(iel)%nodm,                        &
                       CONSTR(iel)%ndofmH, CONSTR(iel)%ndofmE, CONSTR(iel)%ndofmV,  &   
                       CONSTR(iel)%ndofH_macro, CONSTR(iel)%ndofE_macro,            & 
                       CONSTR(iel)%ndofV_macro, CONSTR(iel)%nrconH_macro,           &
                       CONSTR(iel)%nacH_macro, CONSTR(iel)%constrH_macro,           &
                       CONSTR(iel)%nrconE_macro, CONSTR(iel)%nacE_macro,            &
                       CONSTR(iel)%constrE_macro, CONSTR(iel)%nrconV_macro,         &
                       CONSTR(iel)%nacV_macro,CONSTR(iel)%constrV_macro)
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
   subroutine constraints_compute_bounds(MdleC, Nrnod_macro, NrdofH_macro, &
                                         NrdofE_macro, NrdofV_macro)
!
   use data_structure3D,  only: NODES,MAXNODM, ndof_nod, NRNODS
!        
   IMPLICIT NONE
!   
!-----------------------------------------------------------------------
!
   integer, intent(in)  :: MdleC
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
         call find_master(nod,master)
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
!
!
   end module constraints_info




