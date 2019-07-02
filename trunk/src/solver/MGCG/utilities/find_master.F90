!-----------------------------------------------------------------------
!
!   routine name       - find_master
!
!-----------------------------------------------------------------------
!
!   latest revision    - Jan 2018
!
!   purpose            - for a node in a given mesh, routine returns
!                        the corresponding node ancestor marked as
!                        ``master node''            
!
!
!   arguments:
!   in:
!               Nod    - a node
!   out:
!               Master - the corresponding master node
!
!-----------------------------------------------------------------------
!
   subroutine find_master(Igrid,Nod, Master)
!
   use data_structure3D, ONLY : NODES
   use mg_data_structure
!  
   IMPLICIT NONE   
!-----------------------------------------------------------------------
   integer, intent(in)   :: Igrid
   integer   :: Master, Nod, nfath
!
!...check if not a master node itself
   if (NODES_MG(Nod)%master(Igrid).eq.1) then
      Master = Nod; return
   endif
!
! ...loop through ancestors
   nfath = NODES(Nod)%father
   do while(nfath.gt.0)
      if (NODES_MG(nfath)%master(Igrid).eq.1) then
         Master = nfath
         return
      endif
      nfath = NODES(nfath)%father
   enddo
   write(*,*) 'find_master: INCONSISTENCY ! Nod = ',Nod
   stop 1
!
   end subroutine find_master