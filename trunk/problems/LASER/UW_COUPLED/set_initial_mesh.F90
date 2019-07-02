!--------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!--------------------------------------------------------------------
!
!     latest revision   - Jul 18
!
!     purpose           - TODO
!
!     arguments
!            out:       - Nelem_order
!
!--------------------------------------------------------------------
!
subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use parameters
   use data_structure3D
   use CommonParam
!
   implicit none
!
!..polynomial order for initial mesh elements
   integer, dimension(NRELIS), intent(out) :: Nelem_order
!
   integer                 :: iel,ndom,i,ifc,neig
   integer, dimension(6,6) :: ibc
!
!----------------------------------------------------------------------
!
!..check if have not exceeded the maximum order
   if (ORDER_APPROX.gt.MAXP) then
      write(*,1000) 'set_initial_mesh: ',ORDER_APPROX,MAXP
 1000 format(A,'ORDER_APPROX = ',I3,', MAXP = ', I3, '. stop.')
      stop
   endif
!
!..
   write(*,1010) GEOM_NO
 1010 format(' GEOM_NO = ',I1,/)
!
!..loop through initial mesh elements
   do iel=1,NRELIS
!
!  STEP 1 : set up order of approximation (elements may have different
!           orders in each direction)
      select case(ELEMS(iel)%Type)
         case('pris') ; Nelem_order(iel)=ORDER_APPROX*11
         case('bric') ; Nelem_order(iel)=ORDER_APPROX*111
         case('tetr') ; Nelem_order(iel)=ORDER_APPROX*1
         case('pyra') ; Nelem_order(iel)=ORDER_APPROX*1
      end select
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
!     ...single cube with Dirichlet, PEC BC on all faces
			case(1)
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
                     ibc(ifc,1) = 1
                     if((IBCFLAG.eq.3).and.(ifc.eq.2)) then
                        ibc(ifc,2) = 9
                        ibc(ifc,3) = 9
                     else
                        ibc(ifc,2) = 6
                        ibc(ifc,3) = 6
                     endif
                  case default
!                 ...do nothing
               end select
            enddo
!     ...fiber core prism: fpcor_curv or fpcor_str
			case (2)
            call domain_number(iel, ndom)
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
                     ibc(ifc,1) = 1
                     ibc(ifc,2) = 6
                     ibc(ifc,3) = 6
                  case default
!                 ...do nothing
               end select
            enddo
!
!     ...full fiber prism: fpcc_curv or fpcc_str
         case (3)
            call domain_number(iel, ndom)
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
                     ibc(ifc,1) = 1
                     ibc(ifc,2) = 6
                     ibc(ifc,3) = 6
                  case default
!                 ...do nothing
               end select
          enddo
!
!     ...fiber core only hexa: fhcor_curv or fhcor_str
         case (4)
            call domain_number(iel, ndom)
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
                     ibc(ifc,1) = 1
                     ibc(ifc,2) = 6
                     ibc(ifc,3) = 6
               end select
            enddo
!
!     ...full fiber hexa: fhcc_curv or fhcc_str
         case (5)
            call domain_number(iel, ndom)
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
                     ibc(ifc,1) = 1
                     ibc(ifc,2) = 6
                     ibc(ifc,3) = 6
               end select
            enddo
!
!     ...single prism with Dirichlet, PEC BC on 5 out of 6 faces
!     ...Impedance one one face
         case(6)
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
                  case(0)
                     ibc(ifc,1) = 1
                     ibc(ifc,2) = 6
                     ibc(ifc,3) = 6
               end select
            enddo
!
!  ...end select GEOM_NO
      end select
!
!  ...allocate BC flags (one per attribute)
      allocate(ELEMS(iel)%bcond(6))
!
!  ...for each attribute, encode face BC into a single BC flag
      do i=1,6
         call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
      enddo
!
!..end of loop through initial mesh elements
   enddo
!
end subroutine set_initial_mesh
