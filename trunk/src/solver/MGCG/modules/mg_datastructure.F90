!----------------------------------------------------------------------
!
!   module name        - mg_data_structure
!
!----------------------------------------------------------------------
!
!   latest revision    - Sept 2018
!
!   purpose            - module sets up the workspace for multigrid 
!                        solver
!
!----------------------------------------------------------------------
!
module mg_data_structure
!      
   use derived_types
   use data_structure3D
   use macro_grid_info
!
!   
!..coarse solver choice
   integer :: COARSE_SOLVER
   integer, parameter :: NO_CSOLVE      = 0
   integer, parameter :: MUMPS_SOLVER   = 1
   integer, parameter :: PARDISO_SOLVER = 2
!
!..choice of storing or not the Schur complements
   integer :: ISTORE
   integer, parameter :: ISTORE_YES = 1
   integer, parameter :: ISTORE_NO  = 0


   integer :: MAXGEN_PR, NRGRIDS
!   
   type(node_mg), allocatable :: NODES_MG(:)
!
!..global structure for multigrid
   type(ssarray), allocatable :: GRID(:)
!
   contains
!
!
!
   subroutine mg_reset_visit
!
   implicit none

   integer :: i
   do i = 1, NRNODS
      NODES_MG(i)%visit = 0   
   enddo
!
   end subroutine mg_reset_visit
!
!
   subroutine mg_init
!   
   implicit none
!
   integer :: i

   ISTORE=ISTORE_YES
!
!..allocate global structure for each grid
   allocate(GRID(NRGRIDS))
   allocate(NRDOF_MACRO(NRGRIDS))

   allocate(NODES_MG(MAXNODS))

   do i = 1, MAXNODS
      allocate(NODES_MG(i)%orderC(NRGRIDS))
      allocate(NODES_MG(i)%master(NRGRIDS))
      allocate(NODES_MG(i)%iel(NRGRIDS))
      NODES_MG(i)%orderC = 0
      NODES_MG(i)%master = 0
      NODES_MG(i)%iel  = 0
      NODES_MG(i)%visit  = 0
   enddo   
!
   end subroutine mg_init
!
!
!
   subroutine mg_finalize

   implicit none

   integer :: i
!
   write(*,*) 'mg_finalize: deallocating NODES_MG array'
   do i = 1, MAXNODS
      deallocate(NODES_MG(i)%orderC)
      deallocate(NODES_MG(i)%master)
      deallocate(NODES_MG(i)%iel)
   enddo   
!
   deallocate(NODES_MG)
!
!..deallocate global structure for each grid
   write(*,*) 'mg_finalize: deallocating GRID array'
   do i = 1, NRGRIDS
      deallocate(GRID(i)%mdlel)
      if (i .lt. NRGRIDS)  then 
         deallocate(GRID(i)%sch)
      endif
   enddo   
   deallocate(GRID, NRDOF_MACRO)   
!
   end subroutine mg_finalize

end module mg_data_structure
   
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
! !..determine number of dof for a higher order node
!    subroutine find_ndofC(Nod, NdofH,NdofE,NdofV,NdofQ)
! !      
!    implicit none  
! !
!    integer, intent(in)  :: Nod
!    integer, intent(out) :: NdofH,NdofE,NdofV,NdofQ
! !
!    call ndof_nod(NODES(Nod)%type,NODES_MG(Nod)%orderC,    &
!                               NdofH,NdofE,NdofV,NdofQ)
! !
!    end subroutine find_ndofC
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
!    subroutine find_orderC(Mdle, Norder)
! !      
!    use element_data
! !   
!    implicit none
! !   
!    integer, intent(in)  :: Mdle
!    integer, intent(out) :: Norder(19)
! !
! !..locals
!    character(len=4)     :: type
!    integer              :: nodesl(27), norientl(27)
!    integer              :: nrv, nre, nrf, i, j, nod
! !
! !--------------------------------------------------------------------
! !
!    call elem_nodes(Mdle, nodesl,norientl)
! !
!    type = NODES(Mdle)%Type
!    nrv  = nvert(type) ; nre = nedge(type); nrf = nface(type)
! !   
!    do i=1,nre+nrf+1
!       j = nrv+i
!       nod = nodesl(j)
!       Norder(i) = NODES_MG(nod)%orderC
!    enddo
! !
! !
!    end subroutine find_orderC
! !
