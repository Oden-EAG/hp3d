!-----------------------------------------------------------------------
!
!    routine name       - compute_patch_mdleC
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
   subroutine compute_patch_mdleC
   use patch_info
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
! $omp parallel default(shared)  &
! $omp private(iel,mdle,nodesl,norientl,nodm,nrnodm,ip,i,nod,type,loc)
! $omp do schedule(static)
   do iel=1,NRELES_COARSE
      mdle = MDLE_MACRO(iel)
!  ...get information from celem
      call get_connect_infoC(mdle, nodesl,norientl)
!
!   ..get nodes of the modified coarse element 
      call logic_nodesC(mdle,nodesl, nodm,nrnodm)
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
! $omp critical            
            PTCH(loc)%nrmdle=PTCH(loc)%nrmdle+1
! $omp end critical
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
! $omp end do
! $omp end parallel   

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
   end subroutine compute_patch_mdleC