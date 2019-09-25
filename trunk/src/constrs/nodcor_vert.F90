!----------------------------------------------------------------------
!
!   routine name       - nodcor_vert
!
!----------------------------------------------------------------------
!
!   latest revision    - July 2019
!
!   purpose            - routine returns vertex coordinates for
!                        a 3D element
!
!   arguments :
!     in:
!           Mdle       - middle node of an element
!     out:
!           Xnod       - the element local geometry dof
!
!----------------------------------------------------------------------
subroutine nodcor_vert(Mdle, Xnod)
!
   use data_structure3D
   use element_data
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real*8,  intent(out) :: Xnod(NDIMEN,8)
!
   integer :: nodesl(27),norientl(27)
   integer :: iv,nrv,nod,iprint
!
   Xnod(1:NDIMEN,1:8) = 0.d0
!
   call elem_nodes(Mdle, nodesl,norientl)
!
   nrv = nvert(NODES(Mdle)%type)
   do iv=1,nrv
     nod = nodesl(iv)
     Xnod(1:NDIMEN,iv) = NODES(nod)%dof%coord(1:NDIMEN,1)
   enddo
!
   iprint=0
   if (iprint.eq.1) then
      write(*,7001) Mdle
 7001 format('nodcor_vert: VERTEX COORDINATES FOR Mdle = ',i6)
      do iv=1,nrv
         write(*,7002) iv,Xnod(1:NDIMEN,iv)
 7002    format('iv = ',i2,2x,3f8.3)
      enddo
   endif
!
end subroutine nodcor_vert
