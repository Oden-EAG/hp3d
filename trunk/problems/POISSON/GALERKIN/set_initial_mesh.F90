!--------------------------------------------------------------------
!
!     routine name      - set_initial_mesh
!
!--------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - define problem dependent data
!                         (multiphysics, BC, approximation)
!
!     arguments:
!
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
      integer, intent(out) :: Nelem_order(NRELIS)
!
!  ...BC flags
      integer :: ibc(6,NRINDEX)
!
!  ...miscellaneous
      integer :: iprint,ifc,iel,neig
!
!----------------------------------------------------------------------
!
      iprint=0
!
!  ...check if the requested order does not exceed the maximum order
      if (IP.gt.MAXP) then
        write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
        stop 1
      endif
!
!  ...loop over initial mesh elements
      do iel=1,NRELIS
!
!  .....set physics
        ELEMS(iel)%nrphysics = 1
        allocate(ELEMS(iel)%physics(1))
        ELEMS(iel)%physics(1) ='field'
!
!  .....set order of approximation
        if (IP.gt.0) then
!
!  .......uniform order of approximation
          select case(ELEMS(iel)%etype)
            case(TETR); Nelem_order(iel) = 1*IP
            case(PYRA); Nelem_order(iel) = 1*IP
            case(PRIS); Nelem_order(iel) = 11*IP
            case(BRIC); Nelem_order(iel) = 111*IP
          end select
        else
!
!  .......custom order of approximation (NOT IMPLEMENTED)
          write(*,1010) IP
 1010     format('ERROR in set_initial_mesh: IP = ',i3)
          stop 1
        endif
!
!  .....set BC flags: 0 - no BC ; 1 - Dirichlet
        ibc(1:6,1:NRINDEX) = 0
!
        select case(IBC_PROB)
!
!  .....Dirichlet BC
        case(BC_DIRICHLET)
!
!  .......if exterior face, set boundary condition to IBC_PROB
          do ifc=1,nface(ELEMS(iel)%etype)
             neig = ELEMS(iel)%neig(ifc)
             select case(neig)
               case(0); ibc(ifc,1) = 1 ! trace (H1)
             end select
          enddo
!
        case default
          write(*,*) 'set_initial_mesh: problem currently supports Dirichlet BC only'
          stop 1
        end select
!
!  .....allocate BC flags (one per attribute)
        allocate(ELEMS(iel)%bcond(1))
!
!  .....encode face BC into a single BC flag
        call encodg(ibc(1:6,1),10,6, ELEMS(iel)%bcond(1))
!
!  .....print order of approximation
        if (iprint.eq.1 .and. IP.gt.0) then
          write(*,*) '-- uniform order of approximation --'
          write(*,1020) NRELIS
          select case(ELEMS(iel)%etype)
            case(BRIC); write(*,1030) IP,IP,IP
            case(PRIS); write(*,1040) IP,IP
            case(TETR); write(*,1050) IP
            case(PYRA); write(*,1050) IP
          end select
          write(*,*)
!
 1020     format('Element# : ',i4)
 1030     format(' p, q, r = ',3i3)
 1040     format(' p, q = ',   2i3)
 1050     format(' p = ',       i3)
        endif
!
!  ...end of loop over initial mesh elements
      enddo
!
!
      end subroutine set_initial_mesh
