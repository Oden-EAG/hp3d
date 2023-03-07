!--------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 17
!
!     purpose:          - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!
!     in:
!     out:
!           Nelem_order - order for initial mesh elements
!
!---------------------------------------------------------------------
!
   subroutine set_initial_mesh(Nelem_order)
!
   use GMP
   use common_prob_data
   use data_structure3D
!
   implicit none
!
!----------------------------------------------------------------------
!
   integer,dimension(NRELIS),intent(out) :: Nelem_order
!..BC flags
   integer, dimension(6,NR_PHYSA) :: ibc
!..miscellaneous
   integer :: iprint,ifc,iel,neig,iat,iDisplacement
   integer, parameter :: adj_elems(1:6) = (/123,165,171,172,178,220/)

!
!------------------------------------------------------------------------------------
!..initialize
   iprint=0
!
!..check if have not exceeded the maximum order
   if (IP.gt.MAXP) then
      write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
      stop
   endif
!
!..set BC
!..loop over initial mesh elements
   do iel=1,NRELIS
!
!  ...set physics
      ELEMS(iel)%nrphysics = NR_PHYSA
      allocate(ELEMS(iel)%physics(NR_PHYSA))
      do iat=1,NR_PHYSA
         ELEMS(iel)%physics(iat) = PHYSA(iat)
      enddo
!
!  ...set order of approximation
      if (IP.gt.0) then
!     ...uniform order of approximation
         select case(ELEMS(iel)%etype)
            case(BRIC); Nelem_order(iel) = 111*IP
            case(PRIS); Nelem_order(iel) = 11*IP
            case(TETR); Nelem_order(iel) = 1*IP
            case(PYRA); Nelem_order(iel) = 1*IP
         end select
      else
!     ...custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
1003     format('ERROR in set_initial_mesh: IP = ',i3)
         stop
      endif
!
!  ...set BC flags: 0 - no BC ; 1 - Dirichlet ; 2 - Neumann ; 3 - Robin ; >3 - Mixed
      ibc(1:6,1:NR_PHYSA) = 0
!
   select case(IBC_PROB)
!
!  ...uniform BC
      case(BC_DIRICHLET)
!     ...if exterior face, set boundary condition to IBC_PROB
         do ifc=1,nface(ELEMS(iel)%etype)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
!              ibcflag      -> physics variable
               ibc(ifc,1) = 6  ! trace (Hcurl)
            end select
         enddo
!
      case default
         do ifc=1,nface(ELEMS(iel)%etype)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
               ibc(ifc,1) = 9  ! fluxv (H(div))
            end select
         enddo
   end select
!
!  ...allocate BC flags (one per attribute)
      allocate(ELEMS(iel)%bcond(NR_PHYSA))
!  ...for each attribute, encode face BC into a single BC flag
      do iat=1,NR_PHYSA
         call encodg(ibc(1:6,iat),10,6, ELEMS(iel)%bcond(iat))
      enddo
!
!  ...print order of approximation
      if (iprint.eq.1 .and. IP.gt.0) then
         write(*,*) '-- uniform order of approximation --'
         write(*,999) NRELIS
         select case(ELEMS(iel)%etype)
            case(BRIC); write(*,1002) IP,IP,IP
            case(PRIS); write(*,1001) IP,IP
            case(TETR); write(*,1000) IP
            case(PYRA); write(*,1000) IP
         end select
         write(*,*)
!
 999     format('Element# : ',i4)
1000     format(' p = ',       i3)
1001     format(' p, q = ',   2i3)
1002     format(' p, q, r = ',3i3)
      endif
!
!..end of loop over initial mesh elements
   enddo
!
!
   end subroutine set_initial_mesh
