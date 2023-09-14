!----------------------------------------------------------------------
!
!   routine name       - nodcor_vert
!
!----------------------------------------------------------------------
!
!   latest revision    - July 2019
!
!   purpose            - routine returns unconstrained geometry dof
!                        for the vertices of a 3D element
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
   real(8), intent(out) :: Xnod(3,8)
!
!..modified element nodes and corresponding number of dof
   integer :: nodm  (MAXNODM),ndofmH(MAXNODM), &
              ndofmE(MAXNODM),ndofmV(MAXNODM)
!
   integer :: nrconH(MAXbrickH),nacH(NACDIM,MAXbrickH),  &
              nrconE(MAXbrickE),nacE(NACDIM,MAXbrickE),  &
              nrconV(MAXbrickV),nacV(NACDIM,MAXbrickV)
!
   real(8) :: constrH(NACDIM,MAXbrickH),  &
              constrE(NACDIM,MAXbrickE),  &
              constrV(NACDIM,MAXbrickV)
!
!..modified element dof
   real(8) :: val(3,2*MAXbrickH)
!
   integer :: i,j,k,l,kp,nrv,nrnodm
!
!----------------------------------------------------------------------
!
   Xnod(1:3,1:8) = 0.d0
   val (1:3,1:2*MAXbrickH) = 0.d0
!
!..determine constraints' coefficients
   call logic(Mdle,2,                           &
              nodm,ndofmH,ndofmE,ndofmV,nrnodm, &
              nrconH,nacH,constrH,              &
              nrconE,nacE,constrE,              &
              nrconV,nacV,constrV)
!
!..copy the global dof into the local array
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!..initiate the local dof index
   k = 0
!
!..loop through nodes
   do j=1,nrnodm
#if DEBUG_MODE
      if (.not. associated(NODES(nodm(j))%dof)) then
         write(*,*) 'nodcor_vert: dof not associated.'
         stop
      endif
#endif
      do i=1,ndofmH(j)
         k=k+1
         val(1:3,k) = NODES(nodm(j))%dof%coord(1:3,i)
      enddo
   enddo
!
!..#vertices = #local vertex dofs
   nrv = NVERT(NODES(Mdle)%ntype)
!
!..loop through the local vertex dof
   do k=1,nrv
!  ...accumulate for the values
!     nrconH(k) = #unconstrained dofs for i-th local dof
      do kp=1,nrconH(k)
         l = nacH(kp,k)
         Xnod(1:3,k) = Xnod(1:3,k) + constrH(kp,k)*val(1:3,l)
      enddo
   enddo
!
end subroutine nodcor_vert
