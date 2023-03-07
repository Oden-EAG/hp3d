!-------------------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!-------------------------------------------------------------------------------
!
!     latest revision   - June 2021
!
!     purpose:          - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!            out:       - Nelem_order
!                         (polynomial order of each initial mesh element)
!
!-------------------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use parameters
   use data_structure3D
   use commonParam
   use mpi_param, only: RANK,ROOT
!
   implicit none
!
!..polynomial order for initial mesh elements
   integer, intent(out) :: Nelem_order(NRELIS)
!
!..element type
   integer :: etype
!
   integer :: iel,ndom,i,ifc,neig,max_order
!
!..BC flags; dimension = num_faces * components
   integer :: ibc(6,NRINDEX)
!
!-------------------------------------------------------------------------------
!
   Nelem_order = 0
!
!..check if exceeding maximum order
   max_order = max(ORDER_APPROX_X,ORDER_APPROX_Y,ORDER_APPROX_Z)
   if (max_order.gt.MAXP) then
      write(*,1000) 'set_initial_mesh: ',max_order,MAXP
 1000 format(A,'max_order = ',I3,', MAXP = ', I3, '. stop.')
      stop
   endif
!
!..display fiber-specific geometry number
   if (RANK .eq. ROOT) write(*,1010) GEOM_NO
 1010 format(' GEOM_NO = ',I1,/)
!
!..loop through initial mesh elements
   do iel=1,NRELIS
!
      etype = ELEMS(iel)%etype
!
!  STEP 1 : set up order of approximation (elements may have different
!           orders in each direction)
      if (etype .ne. BRIC .and. etype .ne. PRIS) then
         write(*,*) 'set_initial mesh: element type not equal bric/pris. stop.'
         stop
      endif

!  ...uniform isotropic polynomial order
!      select case(etype)
!         case(BRIC) ; Nelem_order(iel)=ORDER_APPROX*111
!         case(PRIS) ; Nelem_order(iel)=ORDER_APPROX*11
!         case(TETR) ; Nelem_order(iel)=ORDER_APPROX*1
!         case(PYRA) ; Nelem_order(iel)=ORDER_APPROX*1
!      end select

!  ...uniform anisotropic polynomial order (geometry axes: x-y-z)
      select case(etype)
         case(PRIS)
            if (ORDER_APPROX_X .ne. ORDER_APPROX_Y) then
               write(*,*) 'set_initial mesh: pris requires ', &
                          'ORDER_APPROX_X=ORDER_APPROX_Y. stop.'
               stop
            endif
            Nelem_order(iel)= ORDER_APPROX_X*10 + &
                              ORDER_APPROX_Z*1
         case(BRIC)
            Nelem_order(iel)= ORDER_APPROX_X*100 + &
                              ORDER_APPROX_Y*10  + &
                              ORDER_APPROX_Z*1
      end select
!
!  ...non-uniform polynomial order (HEXA)
!  ...set p differently in core and cladding of fiber geometry
!  ...TODO: apply min/max rule afterwards
!      if(etype .eq. BRIC) then
!         call domain_number(iel, ndom)
!         write(*,1020) 'iel = ', iel, ' , ndom = ', ndom
! 1020    format(A,i2,A,i2)
!         if (ndom .eq. 1) then
!            Nelem_order(iel) = 777
!         elseif (ndom .le. 5) then
!            Nelem_order(iel) = 666
!         else
!            Nelem_order(iel) = 555
!         endif
!      endif

!
!  STEP 2 : set up physics
!
!     set up the number of physical attributes supported by the element
      ELEMS(iel)%nrphysics=6
      allocate(ELEMS(iel)%physics(6))
      ELEMS(iel)%physics(1)='tempr'
      ELEMS(iel)%physics(2)='EHtr1'
      ELEMS(iel)%physics(3)='EHtr2'
      ELEMS(iel)%physics(4)='hflux'
      ELEMS(iel)%physics(5)='EHfd1'
      ELEMS(iel)%physics(6)='EHfd2'
!
!  ...initialize BC flags
!     0 - no BC
!     1 - Dirichlet BC everywhere
!     2 - Impedance BC via penalty term
!     3 - Impedance BC via elimination
      ibc(1:6,1:NRINDEX) = 0
!
!  ...physics attributes components (1+2+2+1+6+6 = 18 components)
!     phys| comp | description
!     (1) |  1   | H1 field for heat (1 component)
!     (2) |  2-3 | Hcurl for Maxwell trace (\hat E,\hat H) (signal, 2 comps)
!     (3) |  4-5 | Hcurl for Maxwell trace (\hat E,\hat H) (pump  , 2 comps)
!     (4) |  6   | Hdiv trace for heat (1 component)
!     (5) |  7-12| L2 field for Maxwell (E,H) (signal, 6 components)
!     (6) | 13-18| L2 field for Maxwell (E,H) (pump  , 6 components)
!
      select case(GEOM_NO)
!     ...single cube/brick with Dirichlet
!        perfect electrical conductor (PEC) BC on all faces
         case (1)
            do ifc=1,nface(etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
!                 ...Dirichlet on heat
                     ibc(ifc,1) = 1
!                 ...BCs on EH-traces signal and pump
                     if((IBCFLAG.eq.2..or.IBCFLAG.eq.3).and.(ifc.eq.2)) then
!                    ...Impedance on z=L face
!                       REMARK: if IBCFLAG.eq.3 (i.e., using elimination), then
!                               the routine propagate_flag must be called after
!                               refining the mesh to correctly propagate the
!                               impedance flag from faces to edges and vertices
!                    ...Signal
                        ibc(ifc,3) = IBCFLAG !..sets Impedance flag on H_s trace
!                    ...Pump
                        ibc(ifc,5) = IBCFLAG !..sets Impedance flag on H_p trace
                     else
!                    ...Dirichlet on E-trace
                        ibc(ifc,2) = 1 ! Signal trace \hat E_s
                        ibc(ifc,4) = 1 ! Pump   trace \hat E_p
                     endif
                  case default
!                 ...do nothing
               end select
            enddo
!
!     ...fiber core hexa: fhcor_curv or fhcor_str
         case (4)
            do ifc=1,nface(etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
!                 ...Dirichlet on heat
                     ibc(ifc,1) = 1
!                 ...Dirichlet on E-trace
                     ibc(ifc,2) = 1 ! Signal trace \hat E_s
                     ibc(ifc,4) = 1 ! Pump   trace \hat E_p
               end select
            enddo
!
!     ...full fiber hexa: fhcc_curv or fhcc_str
         case (5)
            do ifc=1,nface(etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
!                 ...Dirichlet on heat
                     ibc(ifc,1) = 1
!                 ...Dirichlet on heat flux (for input/output faces z=0,z=L)
                     !if((ifc.eq.1) .or. (ifc.eq.2)) then
                     !   ibc(ifc,6) = 1
                     !else
                     !   ibc(ifc,1) = 1
                     !endif
!                 ...Dirichlet on E-trace
                     ibc(ifc,2) = 1 ! Signal trace \hat E_s
                     ibc(ifc,4) = 1 ! Pump   trace \hat E_p
               end select
            enddo
!
         case default
            write(*,*) 'set_initial_mesh: invalid GEOM_NO. stop.'
            stop
!
!  ...end select GEOM_NO
      end select
!
!  ...allocate BC flags (one per attribute component)
      allocate(ELEMS(iel)%bcond(NRINDEX))
!
!  ...for each component, encode face BC into a single BC flag
      do i=1,NRINDEX
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
!..end of loop through initial mesh elements
   enddo
!
end subroutine set_initial_mesh
