!---------------------------------------------------------------------
!   latest revision    - June 2021
!
!   purpose            - define index for a node using the nodal case
!                        and boundary condition flags for the node
!
!   arguments
!     in:
!          Nod         - a node
!     out:
!          Indexd      - Vector indicating presence and kind of
!                        particular variables 
!
!          Explanation of index
!          For the i-th component supported by the node:
!
!          Indexd(i) = 0  component does not exist
!                    = 1  H1 component with Dirichlet BC flag
!                    = 2  free H1 component
!                    = 3  H(curl) component with Dirichlet BC flag
!                    = 4  free H(curl) component
!                    = 5  H(div) component with Dirichlet BC flag
!                    = 6  free H(div) component
!                    = 7  L2 component with Dirichlet BC flag
!                    = 8  free L2 component
!---------------------------------------------------------------------
!
subroutine get_index(Nod, Indexd)
!
   use data_structure3D
!
   implicit none
!
!..arguments
   integer, intent(in)  :: Nod
   integer, intent(out) :: Indexd(NRINDEX)
!
!..decimal version of NODES(Nod)%case
   integer :: ncase(NR_PHYSA)
!
!..decimal version of NODES(Nod)%bcond
   integer :: ibcd(NRINDEX)
!
!..misc
   integer :: ic,iphys,ivar
!
#if DEBUG_MODE
   integer :: iprint = 0
#endif
!
!-----------------------------------------------------------------------
!
!..decode the case and the BC flags
   call decod(NODES(Nod)%case,2,NR_PHYSA, ncase)
   call decod(NODES(Nod)%bcond,2,NRINDEX, ibcd )
!
!..initiate index counter
   ic=0
!
!..loop through the physics attributes
   do iphys=1,NR_PHYSA
!
      select case(ncase(iphys))
!
!     ...physical attribute is absent, skip its components
         case(0)
            do ivar=1,NR_COMP(iphys)
               ic=ic+1 ; Indexd(ic)=0
            enddo
!
!     ...physical attribute is present
         case(1)
            select case(D_TYPE(iphys))
!
!           ...H1 variable
               case(CONTIN)
!
!              ...loop through components
                  do ivar=1,NR_COMP(iphys)
                     ic=ic+1
                     select case(ibcd(ic))
!
!                    ...free H1 component
                        case(0); Indexd(ic)=2
!
!                    ...component known from Dirichlet BC
                        case(1); Indexd(ic)=1
                     end select
                  enddo
!
!           ...H(curl) variable
               case(TANGEN)
!
!              ...loop through components
                  do ivar=1,NR_COMP(iphys)
                  ic=ic+1
                     select case(ibcd(ic))
!
!                    ...free H(curl) component
                        case(0); Indexd(ic)=4
!
!                    ...component known from Dirichlet BC
                        case(1); Indexd(ic)=3
                     end select
                  enddo
!
!           ...H(div) variable
               case(NORMAL)
!
!              ...loop through components
                  do ivar=1,NR_COMP(iphys)
                     ic=ic+1
                     select case(ibcd(ic))
!
!                    ...free H(div) component
                        case(0); Indexd(ic)=6
!
!                    ...component known from Dirichlet BC
                        case(1); Indexd(ic)=5
                     end select
                  enddo
!
!           ...L2 variable
               case(DISCON)
!
!              ...loop through components
                  do ivar=1,NR_COMP(iphys)
                     ic=ic+1
                     select case(ibcd(ic))
!
!                    ...free H(div) component
                        case(0); Indexd(ic)=8
!
!                    ...component known from Dirichlet BC
                        case(1); Indexd(ic)=7
                     end select
                  enddo
            end select
      end select
!
!..end of loop through physics attributes
   enddo
!
   if (ic.ne.NRINDEX) then
      write(*,*) 'get_index: INCONSISTENCY.'
      stop
   endif
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7010) nod,ncase(1:NR_PHYSA)
 7010   format('get_index: nod = ',i8,' ncase = ',10i2)
        write(*,7020) ibcd(1:NRINDEX)
 7020   format('           ibcd  = ',30i1)
        write(*,7030) Indexd
 7030   format('           Index = ',30i1)
        call pause
      endif
#endif
!
!
end subroutine get_index
