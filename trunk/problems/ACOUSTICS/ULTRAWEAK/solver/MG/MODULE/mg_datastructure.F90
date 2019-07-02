!----------------------------------------------------------------------
!
!   module name        - mg_data_structure
!
!----------------------------------------------------------------------
!
!   latest revision    - MAR 18
!
!   purpose            - modules sets up the workspace for multigrid 
!                        solver
!
!----------------------------------------------------------------------
!
   module mg_data_structure
!      
   use data_structure3D
!
!   
!  maximum # of generation between coarse and fine mesh
   integer :: MAXGEN_PR
!   
!..node_mg data structure
   type node_mg
!  ...order of a nod on the coarse grid
      integer :: orderC
!  ...master flag      
      integer :: master
!  ...visitation flag      
      integer :: visit
!  ...natural ordering of elements
      integer :: iel
!  ...number of dof supported by the node
      integer, allocatable :: nod_ndof(:)
!            
   end type node_mg

   type(node_mg), allocatable :: NODES_MG(:)
!
   contains
!
!
   subroutine mg_reset_visit
!
   do i = 1, NRNODS
      NODES_MG(i)%visit = 0   
   enddo
!
   end subroutine mg_reset_visit
!
!
   subroutine mg_init
!
   allocate(NODES_MG(MAXNODS))
   NODES_MG(:)%orderC = 0
   NODES_MG(:)%master = 0
   NODES_MG(:)%visit  = 0
!
   end subroutine mg_init
!
!
   subroutine mg_finalize
!
   deallocate(NODES_MG)
!
   end subroutine mg_finalize
!
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
!..determine number of dof for a higher order node
   subroutine find_ndofC(Nod, NdofH,NdofE,NdofV,NdofQ)
!      
   implicit none  
!
   integer, intent(in)  :: Nod
   integer, intent(out) :: NdofH,NdofE,NdofV,NdofQ
!
   call ndof_nod(NODES(Nod)%type,NODES_MG(Nod)%orderC,    &
                              NdofH,NdofE,NdofV,NdofQ)
!
   end subroutine find_ndofC
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
   subroutine find_orderC(Mdle, Norder)
!      
   use element_data
!   
   implicit none
!   
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Norder(19)
!
!..locals
   character(len=4)     :: type
   integer              :: nodesl(27), norientl(27)
   integer              :: nrv, nre, nrf, i, j, nod
!
!--------------------------------------------------------------------
!
   call elem_nodes(Mdle, nodesl,norientl)
!
   type = NODES(Mdle)%Type
   nrv  = nvert(type) ; nre = nedge(type); nrf = nface(type)
!   
   do i=1,nre+nrf+1
      j = nrv+i
      nod = nodesl(j)
      Norder(i) = NODES_MG(nod)%orderC
   enddo
!
!
   end subroutine find_orderC
!
!
   end module mg_data_structure
