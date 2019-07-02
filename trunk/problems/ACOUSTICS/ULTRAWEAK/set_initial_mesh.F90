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
   use data_structure3D
   use common_prob_data
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
!
!------------------------------------------------------------------------------------
!..initialize
   iprint=0
!
!..check if have not exceeded the maximum order
   if (IP.gt.MAXP) then
      write(*,*) 'set_initial_mesh: IP, MAXP = ', IP,MAXP
      stop 1
   endif
!
!..add flag `8' to the list of Dirichlet flags
   call add_dirichlet_to_list(8)
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
         select case(ELEMS(iel)%Type)
            case('tetr'); Nelem_order(iel) = 1*IP
            case('pyra'); Nelem_order(iel) = 1*IP
            case('pris'); Nelem_order(iel) = 11*IP
            case('bric'); Nelem_order(iel) = 111*IP
         end select
      else
!     ...custom order of approximation (NOT IMPLEMENTED)
         write(*,1003) IP
1003     format('ERROR in set_initial_mesh: IP = ',i3)
         stop 2
      endif
!
!  ...set BC flags: 0 - no BC ; 1 - Dirichlet ; 2 - Neumann ; 3 - Robin ; >3 - Mixed
      ibc(1:6,1:NR_PHYSA) = 0
      select case(IBC_PROB)
!
!  ...uniform BC
      case(BC_DIRICHLET)
!     ...if exterior face, set boundary condition to IBC_PROB
         do ifc=1,nface(ELEMS(iel)%Type)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
!              ibcflag      -> physics variable
               ibc(ifc,1) = 1  ! trace (H1)
               ibc(ifc,2) = 0  ! fluxv (H(div))
               ibc(ifc,3) = 0  ! field (L2)
!              Note that L2 variables should not take values at the boundary,
!              so their ibc should always be 0 (no BC).
            end select
         enddo
       case(BC_NEUMANN)
!     ...if exterior face, set boundary condition
         do ifc=1,nface(ELEMS(iel)%Type)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
               ibc(ifc,1) = 0  ! trace (H1)
               ibc(ifc,2) = 1  ! fluxv (H(div))
               ibc(ifc,3) = 0  ! field (L2)
            end select
         enddo
      case(BC_IMPEDANCE)
!     ...impedance: eliminate the H(div) dof on the boundary
         do ifc=1,nface(ELEMS(iel)%Type)
            neig = ELEMS(iel)%neig(ifc)
            select case(neig)
            case(0)
               ibc(ifc,1) = 0  ! trace (H1)
               ibc(ifc,2) = 9  ! fluxv (H(div))
               ibc(ifc,3) = 0  ! field (L2)
            end select
         enddo
!  ...cavity problem
      case(BC_CAVITY)

      if (PROB_KIND .ne. PROB_CAVITY) then
         write(*,*) 'set_initial_mesh: BC_CAVITY NOT SET UP FOR FREE SPACE'
         stop 1
      endif   


!     ...select the 6 elements adjacent to the "cavity (for now just a cube)"
!        Dirichlet flag on the H(Div) variable on the "inner face"
         select case(iel)
!
!     ...bottom hexa             
         case(5)
!        ...top face is Dirichlet,    bottom face is impedance
            ibc(2,2) = 1;             ibc(1,2) = 9  
!
!     ...front hexa             
         case(11)
!        ...back face is Dirichlet,   front face is impedance
            ibc(5,2) = 1;             ibc(3,2) = 9  
!     ...left hexa             
         case(13)
!        ...right face is Dirichlet,  left face is impedance
           ibc(4,2) = 1;              ibc(6,2) = 9  
!     ...right hexa             
         case(14)
!        ...left face is Dirichlet,   right face is impedance
           ibc(6,2) = 1;              ibc(4,2) = 9  
!     ...back hexa             
         case(16)
!        ...front face is Dirichlet,  back face is impedance
           ibc(3,2) = 1;              ibc(5,2) = 9  
!     ...top hexa             
         case(22)
!        ...bottom face is Dirichlet, top face is impedance
           ibc(1,2) = 1;              ibc(2,2) = 9  
!     
         case default
            do ifc=1,nface(ELEMS(iel)%Type)
               neig = ELEMS(iel)%neig(ifc)
               select case(neig)
               case(0)
                  ibc(ifc,1) = 0  ! trace (H1)
                  ibc(ifc,2) = 9  ! fluxv (H(div))
                  ibc(ifc,3) = 0  ! field (L2)
               end select
            enddo
         end select
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
         select case(ELEMS(iel)%Type)
            case('tetr'); write(*,1000) IP
            case('pyra'); write(*,1000) IP
            case('pris'); write(*,1001) IP,IP
            case('bric'); write(*,1002) IP,IP,IP
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
