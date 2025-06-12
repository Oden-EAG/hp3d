!-------------------------------------------------------------------------------
!> @brief       Define problem dependent data (multiphysics, BC, approximation)
!!
!> @param[out]  Nelem_order  - order of initial mesh elements
!!
!> @date        July 2023
!-------------------------------------------------------------------------------
subroutine set_initial_mesh(Nelem_order)
!
   use data_structure3D, only: NRELIS
   use physics         , only: NR_PHYSA
   use commonParam
   use mpi_param
!
   implicit none
!
!..polynomial order for initial mesh elements
   integer, intent(out) :: Nelem_order(NRELIS)
!
!..misc
   integer :: attr,bdom,comp,flag,i,nxyz(3)
   integer :: nr_attr,attr_list(NR_PHYSA)
!
!-------------------------------------------------------------------------------
!
!  STEP 1 : set initial element order
!
!..setting uniform isotropic polynomial order "IP"
   call set_order(IP, Nelem_order)

   if (SLAB_GUIDE.eq.1) then 

      do i=1,NRELIS
         ! decode order
         call decod(Nelem_order(i),10,3,nxyz)
         ! replace order for x
         nxyz(1)=1
         ! encode back
         call encod(nxyz,10,3,Nelem_order(i))
      enddo

   endif


!
!  STEP 2 : set up physics
!
!..set up number of physical attributes supported by the element
   nr_attr = NR_PHYSA
   attr_list = [(i, i=1,nr_attr)]
   call set_attr(nr_attr,attr_list)
!
!..BC flag
!  0 - no BC
!  1 - Dirichlet BC
!  2 - impedance BC via penalty term
!  3 - impedance BC via elimination
!
!..physics attributes components
!  attr | comp | index | description
!     1 |  1-2 |   1-2 | Hcurl for Maxwell trace (\hat E,\hat H) (2 components)
!     2 |  1-6 |   3-8 | L2 field for Maxwell (E,H) (6 components)
   attr = 1 ! set BC for Hcurl variable
!
!..boundary domain "0" (Dirichlet BC on E-trace)
   bdom = 0 ! set on all exterior faces with boundary domain "0" (default domain)
   comp = 1 ! E-trace
   flag = 1 ! Dirichlet BC flag
   call set_bcond(bdom,attr,comp,flag)

   if (SLAB_GUIDE.eq.0 .and. IBCFLAG.ge.2) then 
      if (RANK.eq.ROOT) &
      write(*,*) 'set_initial_mesh: IBCFLAG≥2 not implemented for this geometry'
      if (PMLTHUP.eq.0.d0) &
      write(*,*) 'set_initial_mesh: No PML in theta but the mesh will have Dirichlet b.c.'

      call pause
   endif
!
   if (SLAB_GUIDE.eq.1) then
!   ..boundary domain "1"
      bdom=1
!   ..Dirichlet BC on E-trace
      comp = 2 ! H-trace
      flag = 1 ! Dirichlet BC flag
      call set_bcond(bdom,attr,comp,flag)
!    
!   ..boundary domain "2" (Outgoing boundary)
      bdom = 2 ! set on all exterior faces with boundary domain "2"
      ! if (IBCFLAG.eq.2 .or. IBCFLAG.eq.3) then
      if (IBCFLAG.eq.3) then
!     ...impedance BC on H-trace
         comp = 2       ! H-trace
         flag = IBCFLAG ! impedance BC flag
      else
!     ...Dirichlet BC on E-trace
         comp = 2 ! H-trace
         flag = 1 ! Dirichlet BC flag
      endif
      call set_bcond(bdom,attr,comp,flag)
   endif
!
end subroutine set_initial_mesh
