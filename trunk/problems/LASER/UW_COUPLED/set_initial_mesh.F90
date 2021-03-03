!--------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!--------------------------------------------------------------------
!
!     latest revision   - Apr 2019
!
!     purpose           - TODO
!
!     arguments
!            out:       - Nelem_order
!                       (polynomial order of each initial mesh element)
!
!--------------------------------------------------------------------
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
   character(len=4) :: etype
!
   integer :: iel,ndom,i,ifc,neig,max_order
!
!..dimension = num_faces * num_phys
   integer, dimension(6,6) :: ibc
!
!----------------------------------------------------------------------
!
!..check if have not exceeded the maximum order
   max_order = max(ORDER_APPROX_X,ORDER_APPROX_Y,ORDER_APPROX_Z)
   if (max_order.gt.MAXP) then
      write(*,1000) 'set_initial_mesh: ',max_order,MAXP
 1000 format(A,'max_order = ',I3,', MAXP = ', I3, '. stop.')
      stop
   endif
!
!..
   if (RANK .eq. ROOT) write(*,1010) GEOM_NO
 1010 format(' GEOM_NO = ',I1,/)
!
!..loop through initial mesh elements
   do iel=1,NRELIS
!
      etype = ELEMS(iel)%type
!
!  STEP 1 : set up order of approximation (elements may have different
!           orders in each direction)
      if (etype .ne. 'bric' .and. etype .ne. 'pris') then
         write(*,*) 'set_initial mesh: element type not equal bric/pris. stop'
         stop
      endif

!  ...uniform isotropic polynomial order
!      select case(etype)
!         case('pris') ; Nelem_order(iel)=ORDER_APPROX*11
!         case('bric') ; Nelem_order(iel)=ORDER_APPROX*111
!         case('tetr') ; Nelem_order(iel)=ORDER_APPROX*1
!         case('pyra') ; Nelem_order(iel)=ORDER_APPROX*1
!      end select

!  ...uniform anisotropic polynomial order (geometry axes: x-y-z)
      select case(etype)
         case('pris')
            if (ORDER_APPROX_X .ne. ORDER_APPROX_Y) then
               write(*,*) 'set_initial mesh: pris requires ', &
                          'ORDER_APPROX_X=ORDER_APPROX_Y. stop.'
               stop
            endif
            Nelem_order(iel)= ORDER_APPROX_X*10 + &
                              ORDER_APPROX_Z*1
         case('bric')
            Nelem_order(iel)= ORDER_APPROX_X*100 + &
                              ORDER_APPROX_Y*10  + &
                              ORDER_APPROX_Z*1
      end select
!
!  ...non-uniform polynomial order (HEXA)
!  ...set p differently in core and cladding of fiber geometry
!  ...TODO: apply min/max rule afterwards
!      if(etype .eq. 'bric') then
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
      ibc = 0
      select case(GEOM_NO)
!     ...single cube/brick with Dirichlet
!        perfect electrical conductor (PEC) BC on all faces
         case (1)
            do ifc=1,nface(etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
!                 ...dirichlet on heat
                     ibc(ifc,1) = 1
!                 ...Neumann on heat (for min z /max z faces)
!                    TODO (becomes essential BC for normal trace/flux)
!                 ...BCs on EH-traces signal and pump
                     if((IBCFLAG.eq.3).and.(ifc.eq.2)) then
!                    ...impedance on z=L face
!                    ...signal
                        ibc(ifc,2) = 9 !..sets dirichlet flag on H-trace
!                    ...pump
                        ibc(ifc,3) = 9
                     else
!                    ...dirichlet on E-trace
                        ibc(ifc,2) = 6
                        ibc(ifc,3) = 6
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
!                 ...dirichlet on heat
                     ibc(ifc,1) = 1
!                 ...Neumann on heat (for input/output faces)
!                    TODO (becomes essential BC for normal trace/flux)
!                 ...dirichlet on E-trace
                     ibc(ifc,2) = 6 ! signal
                     ibc(ifc,3) = 6 ! pump
               end select
            enddo
!
!     ...full fiber hexa: fhcc_curv or fhcc_str
         case (5)
            do ifc=1,nface(etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
!                 ...dirichlet on heat
                     ibc(ifc,1) = 1
!                 ...Neumann on heat (for input/output faces)
!                    TODO (becomes essential BC for normal trace/flux)
!                 ...dirichlet on E-trace
                     ibc(ifc,2) = 6 ! signal
                     ibc(ifc,3) = 6 ! pump
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
!  ...allocate BC flags (one per attribute)
      allocate(ELEMS(iel)%bcond(6))
!
!  ...for each physics variable, encode face BC into a single BC flag
      do i=1,NR_PHYSA
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
!..end of loop through initial mesh elements
   enddo
!
end subroutine set_initial_mesh
