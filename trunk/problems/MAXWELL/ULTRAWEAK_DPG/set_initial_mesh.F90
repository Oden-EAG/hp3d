
!-------------------------------------------------------------------------------
!> @brief       Define problem dependent data (multiphysics, BC, approximation)
!!
!> @param[out]  Nelem_order  - order of initial mesh elements
!!
!> @date        July 2023
!-------------------------------------------------------------------------------
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
!  ...polynomial order for initial mesh elements
      integer, intent(out) :: Nelem_order(NRELIS)
!
!  ...BC flags; dimension = num_faces * components
      integer :: ibc(6,NRINDEX)
!
!  ...misc
      integer :: etype, iel, ndom, i, ifc, neig, max_order
!
!-------------------------------------------------------------------------------
!
      Nelem_order = 0
!
!  ...check if exceeding maximum order
      if (IP.gt.MAXP) then
         write(*,1000) 'set_initial_mesh: ',max_order,MAXP
 1000 format(A,'max_order = ',I3,', MAXP = ', I3, '. stop.')
         stop
      endif
!
!  ...loop through initial mesh elements
      do iel=1,NRELIS
!
         etype = ELEMS(iel)%etype
!
!     ...set order of approximation
         if (IP.gt.0) then
!        ...uniform order of approximation
            select case(etype)
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
!
!  STEP 2 : set up physics
!
!     ...set up the number of physical attributes supported by the element
         ELEMS(iel)%nrphysics=2
         allocate(ELEMS(iel)%physics(2))
         ELEMS(iel)%physics(1)='EHtrc'
         ELEMS(iel)%physics(2)='EHfld'
!
!     ...initialize BC flags
!        0 - no BC
!        1 - Dirichlet BC everywhere
!        2 - Impedance BC via penalty term
!        3 - Impedance BC via elimination
         ibc(1:6,1:NRINDEX) = 0
!
!     ...physics attributes components (2+6 = 8 components)
!        phys| comp | description
!        (1) |  1-2 | Hcurl for Maxwell trace (\hat E,\hat H) (2 components)
!        (2) |  3-8 | L2 field for Maxwell (E,H) (6 components)
!
!     ...set BCs on exterior faces
         do ifc=1,nface(etype)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
               case(0)
!              ...BCs on EH-traces
                  if((IBCFLAG.eq.2..or.IBCFLAG.eq.3).and.(ifc.eq.2)) then
!                    REMARK: if IBCFLAG.eq.3 (i.e., using elimination), then
!                            the routine propagate_flag must be called after
!                            refining the mesh to correctly propagate the
!                            impedance flag from faces to edges and vertices
                     ibc(ifc,2) = IBCFLAG !..sets Impedance flag on H_s trace
                  else
!                       ...Dirichlet on E-trace
                     ibc(ifc,1) = 1
                  endif
               case default
                  continue
            end select
         enddo
!
!     ...allocate BC flags (one per attribute component)
         allocate(ELEMS(iel)%bcond(NRINDEX))
!
!     ...for each component, encode face BC into a single BC flag
         do i=1,NRINDEX
            call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
         enddo
!
!  ...end of loop through initial mesh elements
      enddo
!
   end subroutine set_initial_mesh
