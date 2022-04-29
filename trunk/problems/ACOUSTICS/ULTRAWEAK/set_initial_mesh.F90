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
      ! integer, parameter :: adj_elems(1:6) = (/5,11,13,14,16,22/)
      ! integer, parameter :: adj_elems(1:6) = (/38,58,62,63,67,87/)
      integer, parameter :: adj_elems(1:6) = (/123,165,171,172,178,220/)
   
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
   !     ...scattering by a thin plate
            if (PROB_KIND .eq. PROB_SCAT_PLANE_PML) then
   !        ...additional dirichlet boundary on interior faces
               select case(iel)
   !
   !        ...4 left hexas             
               ! case(22,26,38,42) 
               case(22) 
   !           ...right face is Dirichlet,
                  ibc(4,2) = 1;     
   !        ...4 right hexas             
               ! case(23,27,39,43)
               case(23)
   !           ...left face is Dirichlet
                  ibc(6,2) = 1;       
               end select
   !
            endif  
   
   
   
   !  ...cavity problem
         case(BC_CAVITY,BC_CAVITY_SCAT)
   
         if (PROB_KIND .ne. PROB_CAVITY .and. PROB_KIND .ne. PROB_SCAT_CAVITY) then
            write(*,*) 'set_initial_mesh: BC_CAVITY NOT SET UP FOR FREE SPACE'
            stop 1
         endif   
   
         ! select case(iel)
         ! case(adj_elems(1),adj_elems(2),adj_elems(3),  &
         !      adj_elems(4),adj_elems(5),adj_elems(6))
   
         ! write(*,*) 'iel = ', iel
         ! write(*,*) 'adj_elems = ', adj_elems(1:6)
   
         ! do ifc = 1, 6
         !    write(*,*) 'ifc  = ', ifc 
         !    write(*,*) 'neig = ', ELEMS(iel)%neig(ifc)
         ! enddo
         ! call pause
         ! end select   
   
   
   !     ...select the 6 elements adjacent to the "cavity (for now just a cube)"
   !        Dirichlet flag on the H(Div) variable on the "inner face"
            select case(iel)
   !
   !     ...bottom hexa             
            case(adj_elems(1))
   !        ...top face is Dirichlet,    bottom face is impedance
               ibc(2,2) = 1;             !ibc(1,2) = 9  
   !
   !     ...front hexa             
            case(adj_elems(2))
   !        ...back face is Dirichlet,   front face is impedance
               ibc(5,2) = 1;             !ibc(3,2) = 9  
   !     ...left hexa             
            case(adj_elems(3))
   !        ...right face is Dirichlet,  left face is impedance
               ibc(4,2) = 1;              !ibc(6,2) = 9  
   !     ...right hexa             
            case(adj_elems(4))
   !        ...left face is Dirichlet,   right face is impedance
               ibc(6,2) = 1;              !ibc(4,2) = 9  
   !     ...back hexa             
            case(adj_elems(5))
   !        ...front face is Dirichlet,  back face is impedance
               ibc(3,2) = 1;              !ibc(5,2) = 9  
   !     ...top hexa             
            case(adj_elems(6))
   !        ...bottom face is Dirichlet, top face is impedance
               ibc(1,2) = 1;              !ibc(2,2) = 9  
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
   
   !  ...scattering from a sphere
         case(BC_SPHERE_SCAT)
   
   
         if (PROB_KIND .ne. PROB_SCAT_SPHERE .and. PROB_KIND .ne. PROB_SPHERE ) then
            write(*,*) 'set_initial_mesh: CHECK BC ' 
            write(*,*) 'PROB_KIND=', PROB_KIND
            stop 1
         endif   
   
   
   !     ...select the 16 elements adjacent to the sphere"
   !        Dirichlet flag on the H(Div) variable on the "inner face"
   
            if (iel .ge. 57 .and. iel .le. 80) then
   !
   !        ...top face is Dirichlet
               ibc(2,2) = 1;        
   !
            else
   !        ...the rest element adjacent to the boundary            
               do ifc=1,nface(ELEMS(iel)%Type)
                  neig = ELEMS(iel)%neig(ifc)
                  select case(neig)
                  case(0)
                     ibc(ifc,1) = 0  ! trace (H1)
                     ibc(ifc,2) = 9  ! fluxv (H(div))
                     ibc(ifc,3) = 0  ! field (L2)
                  end select
               enddo
            endif
   
   
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
   