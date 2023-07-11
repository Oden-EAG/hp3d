!----------------------------------------------------------------------
!> @brief       Define problem dependent data (multiphysics, BC, order)
!!
!> @param[out]   Nelem_order     - order of initial mesh elements
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine set_initial_mesh(Nelem_order)
!
      use GMP
      use data_structure3D
      use common_prob_data_UW
!
      implicit none
!
      integer, intent(out) :: Nelem_order(NRELIS)
!
      integer :: ibc(6,NRINDEX)
      integer :: ifc,iel,neig,iat,iDisplacement
!
!  ...for cavity problem; list of elements adjacent to cavity
      integer, parameter :: adj_elems(1:6) = (/123,165,171,172,178,220/)
!
      integer :: iprint=0
!
!------------------------------------------------------------------------------------
!
!  ...check if have not exceeded the maximum order
      if (IP.gt.MAXP) then
         write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
         stop 1
      endif
!
!  ...set BC
!  ...loop over initial mesh elements
      do iel=1,NRELIS
!
!     ...set physics
         ELEMS(iel)%nrphysics = NR_PHYSA
         allocate(ELEMS(iel)%physics(NR_PHYSA))
         do iat=1,NR_PHYSA
            ELEMS(iel)%physics(iat) = PHYSA(iat)
         enddo
!
!     ...set order of approximation
         if (IP.gt.0) then
!        ...uniform order of approximation
            select case(ELEMS(iel)%etype)
               case(TETR); Nelem_order(iel) = 1*IP
               case(PYRA); Nelem_order(iel) = 1*IP
               case(PRIS); Nelem_order(iel) = 11*IP
               case(BRIC); Nelem_order(iel) = 111*IP
            end select
         else
!        ...custom order of approximation (NOT IMPLEMENTED)
            write(*,1003) IP
1003        format('ERROR in set_initial_mesh: IP = ',i3)
            stop 2
         endif
!
!     ...set BC flags: 0 - no BC ; 1 - Dirichlet ; 2 - Neumann ; 3 - Robin ; >3 - Mixed
         ibc(1:6,1:NRINDEX) = 0
!
         select case(IBC_PROB)
!
!     ...uniform BC
         case(BC_DIRICHLET)
!        ...if exterior face, set boundary condition to IBC_PROB
            do ifc=1,nface(ELEMS(iel)%etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
               case(0)
                  ibc(ifc,1) = 1    ! trace (H1)
                  ibc(ifc,2) = 0    ! fluxv (H(div))
                  ibc(ifc,3:6) = 0  ! field (L2)
               end select
            enddo
         case(BC_NEUMANN)
!        ...if exterior face, set boundary condition
            do ifc=1,nface(ELEMS(iel)%etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
               case(0)
                  ibc(ifc,1) = 0    ! trace (H1)
                  ibc(ifc,2) = 1    ! fluxv (H(div))
                  ibc(ifc,3:6) = 0  ! field (L2)
               end select
            enddo
         case(BC_IMPEDANCE)
!        ...impedance: eliminate the H(div) dof on the boundary
            do ifc=1,nface(ELEMS(iel)%etype)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
               case(0)
                  ibc(ifc,1) = 0    ! trace (H1)
                  ibc(ifc,2) = 3    ! fluxv (H(div))
                  ibc(ifc,3:6) = 0  ! field (L2)
               end select
            enddo
!
!     ...cavity problem
         case(BC_CAVITY)
!
            if (PROB_KIND .ne. PROB_CAVITY) then
               write(*,*) 'set_initial_mesh: BC_CAVITY not set up for free space'
               stop 1
            endif
!
!        ...select the 6 elements adjacent to the "cavity (for now just a cube)"
!           Dirichlet flag on the H(Div) variable on the "inner face"
            select case(iel)
!
!        ...bottom hexa
            case(adj_elems(1))
!           ...top face is Dirichlet
               ibc(2,2) = 1
!
!        ...front hexa
            case(adj_elems(2))
!           ...back face is Dirichlet
               ibc(5,2) = 1
!
!        ...left hexa
            case(adj_elems(3))
!           ...right face is Dirichlet
               ibc(4,2) = 1
!
!        ...right hexa
            case(adj_elems(4))
!           ...left face is Dirichlet
               ibc(6,2) = 1
!
!        ...back hexa
            case(adj_elems(5))
!           ...front face is Dirichlet
               ibc(3,2) = 1
!
!        ...top hexa
            case(adj_elems(6))
!           ...bottom face is Dirichlet
               ibc(1,2) = 1
!     
            case default
               do ifc=1,nface(ELEMS(iel)%etype)
                  neig = ELEMS(iel)%neig(ifc)
                  select case(neig)
                  case(0)
                     ibc(ifc,1) = 0    ! trace (H1)
                     ibc(ifc,2) = 3    ! fluxv (H(div))
                     ibc(ifc,3:6) = 0  ! field (L2)
                  end select
               enddo
            end select
         end select
!
!     ...allocate BC flags (one per attribute)
         allocate(ELEMS(iel)%bcond(NRINDEX))
!
!     ...for each attribute, encode face BC into a single BC flag
         do iat=1,NRINDEX
            call encodg(ibc(1:6,iat),10,6, ELEMS(iel)%bcond(iat))
         enddo
!
!     ...print order of approximation
         if (iprint.eq.1 .and. IP.gt.0) then
            write(*,*) '-- uniform order of approximation --'
            write(*,999) NRELIS
            select case(ELEMS(iel)%etype)
               case(TETR); write(*,1000) IP
               case(PYRA); write(*,1000) IP
               case(PRIS); write(*,1001) IP,IP
               case(BRIC); write(*,1002) IP,IP,IP
            end select
            write(*,*)
!
 999        format('Element# : ',i4)
1000        format(' p = ',       i3)
1001        format(' p, q = ',   2i3)
1002        format(' p, q, r = ',3i3)
         endif
!
!  ...end of loop over initial mesh elements
      enddo
!
   end subroutine set_initial_mesh
