!----------------------------------------------------------------------
!
!   module name        - patch_info
!
!----------------------------------------------------------------------
!
!   latest revision    - MAR 18
!
!   purpose            - modules sets up the workspace for patches
!
!----------------------------------------------------------------------
!
   module patch_info

!..number of patches (i.e, number of coarse grid active vertices)
   integer               :: NRPATCH
   integer, allocatable  :: CGRID_VERTICES(:)  
   integer, parameter    :: MAX_NREDGES = 150
   integer, parameter    :: MAX_NRFACES = 150
   integer, parameter    :: MAX_PATCH_MDLE = 100
!   
   type patch1
!   
!  ...mdle number
      integer :: mdle
!  ...number of local dof for each element in the patch
      integer :: ndof
!  ...connectivity arrays for each mdle node in the patch
      integer, allocatable :: lcon(:)
   end type patch1


!..patch information   
   type patch
!  ...number of coarse grid edges and faces in the patch
      integer :: nredges
      integer :: nrfaces
      integer :: nsz    
!  ...edge and face lists (these are coarse grid nodes)
      integer :: nedgel(MAX_NREDGES), nfacel(MAX_NRFACES)    
!  ...the list of nodes for each patch
      integer, allocatable :: nodl(:)
!
!  ...number of mdle nodes contributing to the patch
      integer :: nrmdle
!  ...list of mdle nodes contributing to the patch
      type(patch1) :: mdlel(MAX_PATCH_MDLE)
!  ...total number of dof in the patch      
      integer :: nrdof
!  ...connectivity map form patch to global
      integer, allocatable :: lcon(:)
!  ...Cholesky decomposition of the assembled block in Lapack packed form (dense matrix)
#if C_MODE
      complex*16, allocatable :: zAp(:) 
#else
      real*8,    allocatable :: zAp(:) 
#endif

!  
   end type patch

   type (patch), allocatable :: PTCH(:)

   integer, allocatable, save :: DESCENDANTS(:)
!$omp threadprivate(DESCENDANTS)   

! 
!
   CONTAINS
!
!-----------------------------------------------------------------------
!
!    routine name       - compute_patch_mdle
!
!-----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - routine computes and stores the mdle nodes
!                         contributing to each patch
!
!    arguments          - none
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   subroutine compute_patch_mdle
   use data_structure3D,  only: NODES,MAXNODM
   use mg_data_structure, only: NODES_MG
   use macro_grid_info,   only: NRELES_COARSE, MDLE_MACRO
!
   IMPLICIT NONE
!      
   character*4 ::  type
!..work space for elem_nodes and logic_nodes
   integer     :: nodesl(27),norientl(27), nodm(MAXNODM)
   integer     :: nod, loc, j, i, mdle, iel, ip, nrnodm
   integer     :: iprint
!
!----------------------------------------------------------------------
!

   iprint = 0

   allocate(PTCH(NRPATCH))
!      
!..initialise      
   do j = 1, NRPATCH
      PTCH(j)%nrmdle = 0 
      PTCH(j)%mdlel(:)%mdle = 0 
   enddo   
! 
!$omp parallel default(shared)  &
!$omp private(iel,mdle,nodesl,norientl,nodm,nrnodm,ip,i,nod,type,loc)
!$omp do schedule(dynamic)
   do iel=1,NRELES_COARSE
      mdle = MDLE_MACRO(iel)
!  ...get information from celem
      call get_connect_info(mdle, nodesl,norientl)
!
!   ..get nodes of the modified coarse element 
      call logic_nodes(mdle,nodesl, nodm,nrnodm)
!      
!  ...update the maximum number of local dof
      ip = 0  
!  ...loop through the nodes (no need to go through the last (mdle) node)      
      do i = 1,nrnodm-1
!     ...pick up the node number         
         nod = nodm(i)
!     ...check if the node is a vertex         
         type = NODES(nod)%type
         if (type .eq. 'vert') then
            loc = NODES_MG(nod)%visit
!$omp critical            
            PTCH(loc)%nrmdle=PTCH(loc)%nrmdle+1
!$omp end critical
            ip = PTCH(loc)%nrmdle
            if (ip .gt. MAX_PATCH_MDLE) then
               write(*,*) 'compute_patch_mdle : Increase MAX_PATCH_MDLE'
               write(*,*) 'compute_patch_mdle : Patch = ', loc
               stop 1
            endif
            PTCH(loc)%mdlel(ip)%mdle = mdle  
         endif   
      enddo
   enddo   
!$omp end do
!$omp end parallel   

   if (iprint.eq.1) then 
      do j = 1, NRPATCH
         write(*,999) j, CGRID_VERTICES(j)
  999    format('compute_patch_mdle: patch no, vert = ', 2i7)
         write(*,1000) PTCH(j)%nrmdle
 1000    format('compute_patch_mdle: nrmdle         = ', i7)
         write(*,1001) PTCH(j)%mdlel(1:PTCH(j)%nrmdle)%mdle
 1001    format('compute_patch_mdle: mdlel          = ', 30(5i7,/,37x))        
         write(*,*) ''
      enddo
   endif   

!
!
   end subroutine compute_patch_mdle
!
!
!
! -----------------------------------------------------------------------
!
!    routine name       - patch_nodes
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - routine fills in the patch database
!
!
! ----------------------------------------------------------------------
!
   subroutine patch_nodes
!
   use macro_grid_info,  ONLY: NRELES_COARSE, MDLE_MACRO
!
   IMPLICIT NONE
! 
   integer              :: mdle, iel,ip, nf, ne
   integer              :: maxnlist, max_nrnodes
   integer              :: i, j, N, k, nod, ii, nloc, nlist, nre, nrf
   integer, allocatable, save :: nodl(:), nodla(:)
!$omp threadprivate(nodl, nodla)   
   integer :: iprint
!---------------------------------------------------------------------------------
!
   iprint = 0
!
!..initialise
   PTCH(:)%nredges=0 ; PTCH(:)%nrfaces=0 ; 

! $omp parallel default(shared)  &
! $omp private(iel,mdle)
! $omp do schedule(guided)
   do iel = 1, NRELES_COARSE
      mdle = MDLE_MACRO(iel)
      call patch_contr(mdle)
   enddo
! $omp end do
! $omp end parallel
!


!$omp parallel default(shared)  &
!$omp private(j, maxnlist,max_nrnodes,k,nre,nrf,i,nod,nlist,ii,nloc)
!$omp do schedule(guided)
   do j=1,NRPATCH

      maxnlist = 100
      allocate(DESCENDANTS(maxnlist))

      max_nrnodes = 100
      allocate(PTCH(j)%nodl(max_nrnodes)) ; PTCH(j)%nodl = 0
!
!  ...build final list of patch nodes 
!  ...first unique the vertex node
      PTCH(j)%nodl(1) = CGRID_VERTICES(j)
!  ...now the nodes on the edges and faces     
      k=1
      nre = PTCH(j)%nredges ; nrf = PTCH(j)%nrfaces
!
      allocate(nodl(nre+nrf)) 
      nodl(1:nre) = PTCH(j)%nedgel(1:nre)
      nodl(nre+1:nre+nrf) = PTCH(j)%nfacel(1:nrf)
!      
      do i = 1,nre+nrf
         nod = nodl(i)
         call find_descendants(nod, nlist)
!     ...avoid repetition of nodes
         do ii = 1,nlist
            call locate(DESCENDANTS(ii),PTCH(j)%nodl(1:k),k,nloc)
            if (nloc .eq. 0) then
               if (k .eq. max_nrnodes) then
                  max_nrnodes = 2*k
                  allocate(nodla(max_nrnodes)); nodla(1:k) = PTCH(j)%nodl(1:k)
                  call move_alloc(nodla,PTCH(j)%nodl)
               endif
               k = k + 1
               PTCH(j)%nodl(k) = DESCENDANTS(ii)
            endif
         enddo   
      enddo
      deallocate(nodl)
      deallocate(DESCENDANTS)
      PTCH(j)%nsz = k
   enddo
!$omp end do
!$omp end parallel   
!
   if (iprint .eq. 1) then 
      do j=1,NRPATCH
         write(*,1002) j
 1002    format('patch_nodes: patch no = ', i7)        
         write(*,1003) PTCH(j)%nodl(1:PTCH(j)%nsz)
 1003    format('patch_nodes: nodl     = ', 30(9i7,/,24x))        
         write(*,*) ''
      enddo
   endif   
!
!
   end subroutine patch_nodes
!   
!
!
! -----------------------------------------------------------------------
!
!    routine name       - patch_contr
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - routine computes the list of coarse grid
!                         edge and face nodes for each vertex patch
!                         (both active and inactive)
!
! ----------------------------------------------------------------------
!
!
   subroutine patch_contr(Mdle)
!
   use data_structure3D
   use mg_data_structure
   use constrained_nodes
   use constraints

   IMPLICIT NONE
!
   integer, intent(in)  :: Mdle
!
!..locals
!..element characteristics
   character(len=4)     :: type, ftype, facetype
!   
!..workspace for get_connect_infoC
   integer              :: nodesl(27),norientl(27)
!   
   integer              :: nrv, nre, nrf, ie, iedge(2), if, nrvf, iface(4)  
   integer              :: listv(16), nc, icase, medg, mface
   integer              :: i, j , iv, nodv, jjv, nodvp
   integer              :: jmax, loc, loc1, l
!
!---------------------------------------------------------------------------------
!
!..element type
   type = NODES(Mdle)%type
!
   nrv = nvert(type); nre = nedge(type); nrf = nface(type)
!   
!..get element nodes and build the local database on constraints 
!..for the element
   call get_connect_infoC(Mdle, nodesl,norientl)
!
!..loop through macro-element faces
   do ie = 1, nre
!
!  ...local face nod number
      i = nrv + ie
      medg = nodesl(i)
!
!  ...pick up the edge vertex nodes 
      call edge_to_vert(type,ie, iedge(1),iedge(2))
!      
!  ...establish vertex patches to which `medg' contributes
      j=0
      do iv = 1,2
         nodv = nodesl(iedge(iv))
         if (NODES_MG(nodv)%master.eq.1) then
            j=j+1; listv(j) = nodv
         else
            call decode2(NODES_CONSTR(iedge(iv)), nc,icase)
            select case (icase)
!        ...vertex constrained by an edge               
            case(13,39,49)            
               do jjv=1,2
                  nodvp = NEDG_CONS(jjv,nc)
                  if (NODES_MG(nodvp)%master.ne.1) stop 3
                  j=j+1; listv(j) = nodvp
               enddo
!        ...vertex node constrained by a face
            case(29)
               do jjv=1,4
                  nodvp = NFACE_CONS(4+jjv,nc)
                  if (NODES_MG(nodvp)%master.ne.1) stop 3
                  j=j+1; listv(j) = nodvp
               enddo
            end select   
         endif
      enddo         
!
!  ...add 'medg' to the vertex patches
      jmax=j
      do j=1,jmax
!
!     ...location of the vertex on the coarse grid vertices list
!     ...this visit flag is raised while solving the coarse grid system
         loc = NODES_MG(listv(j))%visit
!
         call locate(medg,PTCH(loc)%nedgel,PTCH(loc)%nredges,loc1)
         if (loc1.eq.0) then
!$omp critical
            PTCH(loc)%nredges = PTCH(loc)%nredges+1
!$omp end critical
            l=PTCH(loc)%nredges
            if (l .gt. MAX_NREDGES) then 
               write(*,*) 'patch_contr: increase MAX_NREDGES' 
               stop 1
            endif
            PTCH(loc)%nedgel(l) = medg
         endif
! 
!  ...end of loop through patches
      enddo
!
!..end loop through macro-element edges
   enddo
!
!..loop through macro-element faces
   do if = 1, nrf
!
!  ...local face nod number
      i = nrv + nre + if
      mface = nodesl(i)
!
!  ...face type
      ftype = facetype(type,if)
!      
!  ...number of vertices on the face (3 or 4)
      nrvf = nvert(ftype)
!
!  ...pick up the face vertex nodes 
      call face_to_vert(type,if, iface(1),iface(2),iface(3),iface(4))
!      
!  ...establish vertex patches to which `mface' contributes
      j=0
      do iv = 1,nrvf
         nodv = nodesl(iface(iv))
         if (NODES_MG(nodv)%master.eq.1) then
               j=j+1; listv(j) = nodv
         else
            call decode2(NODES_CONSTR(iface(iv)), nc,icase)
            select case(icase)
!        ...vertex constrained by an edge               
            case(13,39,49)   
               do jjv=1,2
                  nodvp = NEDG_CONS(jjv,nc)
                  if (NODES_MG(nodvp)%master.ne.1) stop 3
                  j=j+1; listv(j) = nodvp
               enddo
!        ...vertex constrained by a face
            case (29)              
               do jjv=1,nrvf
                  nodvp = NFACE_CONS(4+jjv,nc)
                  if (NODES_MG(nodvp)%master.ne.1) stop 3
                  j=j+1; listv(j) = nodvp
               enddo
            end select   
         endif
      enddo         
! 
!  ...add 'mface' to the vertex patches
      jmax=j
      do j=1,jmax
!
!     ...location of the vertex on the coarse grid vertices list
!     ...this visit flag is raised while solving the coarse grid system
         loc = NODES_MG(listv(j))%visit
!         
         call locate(mface,PTCH(loc)%nfacel,PTCH(loc)%nrfaces,loc1)
         if (loc1.eq.0) then
!$omp critical
            PTCH(loc)%nrfaces = PTCH(loc)%nrfaces+1
!$omp end critical
            l=PTCH(loc)%nrfaces
            if (l .gt. MAX_NRFACES) then 
               write(*,*) 'patch_contr: increase MAX_NRFACES' 
               stop 2
            endif
            PTCH(loc)%nfacel(l) = mface
         endif
!  ...end of loop through patches
      enddo
!
!..end loop through macro-element faces
   enddo
!
!   
   end subroutine patch_contr
!
!
! -----------------------------------------------------------------------
!
!    routine name       - find_descendants
!
! -----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - routine saves the list of active descendants
!                         for an edge of a face of a macro element
!
!            in:        
!                  nod  - edge or face node of a macro element
!               MAXlist - maximum size of the list
!           out:    
!                 List  - list of descendants
!                 Nlist - Size of the list
!
!
! ----------------------------------------------------------------------
!
!
   subroutine find_descendants(Nod, Nlist)
!
   use data_structure3D, only: NODES         
   use refinements_2D,  only : nr_sons
!   
   IMPLICIT NONE
!   
!-----------------------------------------------------------------------
!
   integer, intent(in) :: Nod
!  
   integer, intent(out):: Nlist
!
!..locals
   integer :: maxl, maxnlist
   integer,allocatable :: lista(:), listb(:)
   integer :: lmax, nrsons, nson, nfath,loc, is,l
!
   Nlist=0 ; maxnlist = ubound(DESCENDANTS,1)
!..initial list allocation
   maxl = 50
   allocate(lista(maxl))   
!
!..place 'Nod' on the auxiliary list
   if (NODES(Nod)%ref_kind.eq.0) then
      if (NODES(Nod)%act.eq.1) then
         Nlist=Nlist+1
         DESCENDANTS(Nlist) = Nod
         return
      else
         return
      endif
   else
      lmax=1
      lista(lmax) = Nod
   endif
!
!..loop through (possibly growing) list of descendants
   l=1
   do while (l.le.lmax)
      nfath = lista(l)
      call nr_sons(NODES(nfath)%type,NODES(nfath)%ref_kind, nrsons)
      do is=1,nrsons
         nson = NODES(nfath)%sons(is)
         if (NODES(nson)%ref_kind.eq.0) then
            if (NODES(nson)%act.eq.1) then
               if (Nlist .ge. maxnlist) then 
                  maxnlist = 2*Nlist
                  allocate(listb(maxnlist)); listb(1:Nlist) = DESCENDANTS(1:Nlist)  
                  call move_alloc(listb,DESCENDANTS)
               endif
               Nlist=Nlist+1
               DESCENDANTS(Nlist) = nson
            else
               call locate(nfath,DESCENDANTS,Nlist,loc)
               if (loc .eq. 0) then
                  if (Nlist .ge. maxnlist) then 
                     maxnlist = 2*Nlist ; 
                     allocate(listb(maxnlist)); listb(1:Nlist) = DESCENDANTS(1:Nlist)  
                     call move_alloc(listb,DESCENDANTS)
                  endif
                  Nlist=Nlist+1
                  DESCENDANTS(Nlist) = nfath
               endif
            endif
         else
            if (lmax.eq.maxl) then 
               maxl = 2*lmax
               allocate(listb(maxl)); listb(1:lmax) = lista(1:lmax)
               call move_alloc(listb,lista)
            endif   
            lmax=lmax+1 
            lista(lmax) = nson
         endif
      enddo
      l=l+1
   enddo

   deallocate(lista)
!
   end subroutine find_descendants
!
!


   subroutine deallocate_patches
!
   implicit none
   integer :: i
!   
   do i = 1, NRPATCH
      deallocate(PTCH(i)%lcon)
      deallocate(PTCH(i)%zAp)
      deallocate(PTCH(i)%nodl)
   enddo   

   end subroutine deallocate_patches


   end module patch_info
