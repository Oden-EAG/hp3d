!-----------------------------------------------------------------------
! 
!    routine name       - mark_masters
! 
!-----------------------------------------------------------------------
! 
!    latest revision    - Jan 2018
! 
!    purpose            - for a given mesh, routine marks all ACTIVE
!                         nodes in the mesh as `master nodes'            
! 
!    arguments          - none
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
   subroutine mark_masters
! 
                               
   use data_structure3D,  ONLY: NODES, NRELES, MAXNODM, MAXNODS, NR_PHYSA 
                               
   use assembly,          ONLY: assembly_begin, assembly_end
   use mg_data_structure, ONLY: NODES_MG
!
   IMPLICIT NONE
!   
!-----------------------------------------------------------------------
!
!..work space for celem
   integer    :: nrdofs(NR_PHYSA)
   integer    :: nodm(MAXNODM),   ndofmH(MAXNODM)
   integer    :: ndofmE(MAXNODM), ndofmV(MAXNODM), ndofmQ(MAXNODM)
   integer    :: nrdofm, nrdofc, nrnodm
   integer    :: nod, mdle, iel, i, idec
#if C_MODE   
   complex*16 :: zvoid
#else
   real*8     :: zvoid
#endif   
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   call assembly_begin 
!
!..clean up and copy order
   do nod=1,MAXNODS
      NODES_MG(nod)%master=0
!  ...save the current order
      NODES_MG(nod)%orderC = NODES(nod)%order     
   enddo
! 
!..loop through active elements
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
! 
!  ...determine nodes of the modified element
      idec=1
      call celem1(mdle,idec, nrdofs,nrdofm,nrdofc,          &
                 nodm,ndofmH,ndofmE,ndofmV,ndofmQ, nrnodm, &
                 zvoid,zvoid)
!
!  ...loop through the nodes of the element and mark them
      do i=1,nrnodm
         nod = nodm(i)
         NODES_MG(nod)%master=1
      enddo
   enddo
!
   call assembly_end
!
!
   end subroutine mark_masters
