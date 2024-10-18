!-----------------------------------------------------------------------------
!> @brief Routine performs maps master element of the child (obtained
!         after isotropic refinement of the coarse element) elements to the
!         master element corresponding to the coarse element.
!> @param[in] Iel    - Index of the child (e.g.: 1-8 for isotropic refinement of hex element)
!> @param[in] Xi     - quadrature point correspoding to the child
!> @param[in] Etype  - element type
!> @param[out] Xis   -  quadrature point  corresponding to the coarse element
!> @date June 2024
!-----------------------------------------------------------------------------
subroutine fine_to_coarse_gp_map(Iel,Xi,Xis,Etype)
!
   use element_data
   implicit none
!
   real(8), intent(in) :: Xi(3)
   integer, intent(in) :: Iel
   integer, intent(in) :: Etype
!
   real(8), intent(out) :: Xis(3)
!
   select case(etype)
      case(MDLB)
         select case(iel)
               case(1)
                  xis  = 0.5 * Xi
               case(2)
                  xis(1) = 0.5 + 0.5 * Xi(1)
                  Xis(2) = 0.5 * Xi(2)
                  Xis(3) = 0.5 * Xi(3)
               case(3)
                  Xis(1) = 0.5 + 0.5 * Xi(1)
                  Xis(2) = 0.5 + 0.5 * Xi(2)
                  Xis(3) = 0.5 * Xi(3)
               case(4)
                  Xis(1) = 0.5 * Xi(1)
                  Xis(2) = 0.5 + 0.5 * Xi(2)
                  Xis(3) = 0.5 * Xi(3)
               case(5)
                  Xis(1) = 0.5 * Xi(1)
                  Xis(2) = 0.5 * Xi(2)
                  Xis(3) = 0.5 + 0.5 * Xi(3)
               case(6)
                  Xis(1) = 0.5 + 0.5 * Xi(1)
                  Xis(2) = 0.5 * Xi(2)
                  Xis(3) = 0.5 + 0.5 * Xi(3)
               case(7)
                  Xis(1) = 0.5 + 0.5 * Xi(1)
                  Xis(2) = 0.5 + 0.5 * Xi(2)
                  Xis(3) = 0.5 + 0.5 * Xi(3)
               case(8)
                  Xis(1) = 0.5 * Xi(1)
                  Xis(2) = 0.5 + 0.5 * Xi(2)
                  Xis(3) = 0.5 + 0.5 * Xi(3)
               case default
                  write(*,*) "Wrong son"
         end select
      case default
         write(*,*) "not implemented yet"
   end select
!
end subroutine fine_to_coarse_gp_map
!
!-----------------------------------------------------------------------------
!> @brief Routine performs maps master element of the child (obtained
!         after isotropic refinement of the coarse element) elements to the
!         master element corresponding to the h-ref candidate.
!> @param[in] Etype  - element type
!> @param[in] Kref   - h-ref flag of the h-ref candidate
!> @param[in] Sel    - index of the child element
!> @param[in] Xi     - quadrature point correspoding to the child element
!> @param[out] Xis   - mapped quadrature point  corresponding to the h-ref candidate
!> @date June 2024
!-----------------------------------------------------------------------------
subroutine  fine_to_subson_gp_map(Etype,Kref,Sel,Xi,Xis)
!
   use element_data
   implicit none
   integer,    intent(in)  :: Kref
   integer,    intent(in)  :: Sel
   integer,    intent(in)  :: Etype
!
   real(8), intent(in)  :: Xi(3)
   real(8), intent(out) :: Xis(3)
!
   if(Etype .eq. MDLB) then
      select case(Kref)
         case(100)
            select case(Sel)
               case(1,2)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,0.5d0,0.5d0,Xi,Xis)
               case(4,3)
                  call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,1.d0,0.5d0,0.5d0,Xi,Xis)
               case(5,6)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,1.d0,0.5d0,0.5d0,Xi,Xis)
               case(8,7)
                  call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.5d0,1.d0,0.5d0,0.5d0,Xi,Xis)

               case default
                  write(*,*) "wrong fine son"
            end select
!
         case(010)
            select case(Sel)
               case(1,4)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,1.d0,0.5d0,Xi,Xis)
               case(2,3)
                  call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.d0,0.5d0,1.d0,0.5d0,Xi,Xis)
               case(5,8)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,0.5d0,1.d0,0.5d0,Xi,Xis)
               case(6,7)
                  call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.5d0,0.5d0,1.d0,0.5d0,Xi,Xis)
               case default
                  write(*,*) "wrong fine son"
            end select
!
         case(001)
            select case(Sel)
               case(1,5)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,0.5d0,1.d0,Xi,Xis)
               case(2,6)
                  call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.d0,0.5d0,0.5d0,1.d0,Xi,Xis)
               case(3,7)
                  call singlefineson_to_subson_gp_map(0.5d0,0.5d0,0.d0,0.5d0,0.5d0,1.d0,Xi,Xis)
               case(4,8)
                  call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,0.5d0,0.5d0,1.d0,Xi,Xis)
               case default
                  write(*,*) "wrong fine son"
            end select
!
         case(110)
            select case(Sel)
               case(1,2,3,4)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,1.d0,0.5d0,Xi,Xis)
               case(5,6,7,8)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,1.d0,1.d0,0.5d0,Xi,Xis)
            end select
!
         case(101)
            select case(Sel)
               case(1,2,5,6)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,0.5d0,1.d0,Xi,Xis)
               case(4,3,8,7)
                  call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,1.d0,0.5d0,1.d0,Xi,Xis)
            end select
!
         case(011)
            select case(Sel)
               case(1,4,5,8)
                  call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,Xi,Xis)
               case(2,3,6,7)
                  call singlefineson_to_subson_gp_map(0.5d0,0.0d0,0.d0,0.5d0,1.d0,1.d0,Xi,Xis)
            end select
!
         case(111)
            Xis = Xi
         case default
            write(*,*) "wrong h-ref option"
            stop 
      end select
   endif
!
end subroutine fine_to_subson_gp_map
!
!-----------------------------------------------------------------------------
!> @brief   Routine performs mapping between a master element corresponding 
!           to a fine child element and the master element of a h-ref
!           candidate using offset along x,y,z and scaling factors along x,y,z.
!> @param[in] ofx    - offset along x axis
!> @param[in] ofy    - offset along y axis
!> @param[in] ofz    - offset along z axis
!> @param[in] fx     - scaling along x axis
!> @param[in] fy     - scaling along y axis
!> @param[in] fz     - scaling along z axis
!
!> @param[in] Xi     -  quadrature point correspoding to the child element
!> @param[out] Xis   -  mapped quadrature point  corresponding to the h-ref candidate
!> @date June 2024
!-----------------------------------------------------------------------------
subroutine singlefineson_to_subson_gp_map(Ofx,Ofy,Ofz,Fx,Fy,Fz,Xi,Xis)
!
   implicit none
   real(8), intent(in) :: Ofx
   real(8), intent(in) :: Ofy
   real(8), intent(in) :: Ofz
   real(8), intent(in) :: Fx
   real(8), intent(in) :: Fy
   real(8), intent(in) :: Fz
   real(8),dimension(3),   intent(in)    :: Xi
!
   real(8),dimension(3),   intent(out)   :: Xis

    Xis(1) = Ofx + Fx * Xi(1)
    Xis(2) = Ofy + Fy * Xi(2)
    Xis(3) = Ofz + Fz * Xi(3)
!
end subroutine singlefineson_to_subson_gp_map

!---------------------------------------------------------------------------------------------------------------
!> @brief    Routine interpolates polynomial order for the children elements if the ref applied is different
!>           from the one intented (Generally happens due to 1-irregularity enforcement).
!>  @param[in] Etype : Type of the element
!>  @param[in] Kref_intent    - Intented refinement obtained through projections
!>  @param[in] Kref_appl      - the refinement applied by the code while preserving 1-irreguarity
!>  @param[in] Nr_sons_intent - number of sons for the intented refinement
!>  @param[in] Nr_sons_appl   - number of sons produced due to applied refinement
!>  @param[in] Pref_intent    - p-ref that is compatible with kref_intent
!
!>  @param[out] Pref_appl     - p-ref that should be applied for the applied h-refinement
!> @date Feb 2023
!---------------------------------------------------------------------------------------------------------------
subroutine subson_one_irregularity_map(Etype,Kref_intent,Kref_appl,Nr_sons_intent,Pref_intent,Nr_sons_appl,Pref_appl)
!
   use element_data
   implicit none
!
   integer,    intent(in) :: Etype
   integer,    intent(in) :: Kref_intent
   integer,    intent(in) :: Kref_appl
   integer,    intent(in) :: Nr_sons_intent
   integer,    intent(in) :: Pref_intent(Nr_sons_intent)
   integer,    intent(in) :: Nr_sons_appl
   integer,    intent(out):: Pref_appl(Nr_sons_appl)
!    
   integer :: i
   integer, allocatable :: el_pmap(:)
!
   if(Etype .eq. MDLB) then
!
      allocate(el_pmap(8))
!  ...making an transfer map using an isotropic refinement
      select case(Kref_intent)
!
         case(100)
            do i = 1,8
               select case(i)
                     case(1,4,5,8)
                        el_pmap(i) = Pref_intent(1)
                     case(2,3,6,7)
                        el_pmap(i) = Pref_intent(2)
               end select
            enddo
!
         case(010)
            do i = 1,8
               select case(i)
                  case(1,2,5,6)
                     el_pmap(i) = Pref_intent(1)
                  case(3,4,7,8)
                     el_pmap(i) = Pref_intent(2)
               end select
            enddo
!
         case(001)
            do i = 1,8
               select case(i)
                  case(1,2,3,4)
                     el_pmap(i) = Pref_intent(1)
                  case(5,6,7,8)
                     el_pmap(i) = Pref_intent(2)
               end select
            enddo
!
         case(110)
            do i = 1,8
               select case(i)
                  case(1,5)
                     el_pmap(i) = Pref_intent(1)
                  case(2,6)
                     el_pmap(i) = Pref_intent(2)
                  case(3,7)
                     el_pmap(i) = Pref_intent(3)
                  case(4,8)
                     el_pmap(i) = Pref_intent(4)
               end select
            enddo
!
         case(101)
            do i = 1,8
               select case(i)
                  case(1,4)
                     el_pmap(i) = Pref_intent(1)
                  case(2,3)
                     el_pmap(i) = Pref_intent(2)
                  case(6,7)
                     el_pmap(i) = Pref_intent(3)
                  case(5,8)
                     el_pmap(i) = Pref_intent(4)
               end select
            enddo
!
         case(011)
            do i = 1,8
               select case(i)
                  case(1,2)
                     el_pmap(i) = Pref_intent(1)
                  case(3,4)
                     el_pmap(i) = Pref_intent(2)
                  case(7,8)
                     el_pmap(i) = Pref_intent(3)
                  case(5,6)
                     el_pmap(i) = Pref_intent(4)
               end select
            enddo
!
         case(111)
            do i = 1,8
               el_pmap(i) = Pref_intent(i)
            enddo
!
      end select
!  ...transferring the p-order to kref_close childs
      select case(Kref_appl)
!
         case(100)
            Pref_appl(1) = (el_pmap(1)+el_pmap(4)+el_pmap(5)+el_pmap(8))/4
            Pref_appl(2) = (el_pmap(2)+el_pmap(3)+el_pmap(6)+el_pmap(7))/4
!
         case(010)
            Pref_appl(1) = (el_pmap(1)+el_pmap(2)+el_pmap(5)+el_pmap(6))/4
            Pref_appl(2) = (el_pmap(3)+el_pmap(4)+el_pmap(7)+el_pmap(8))/4
!
         case(001)
            Pref_appl(1) = (el_pmap(1)+el_pmap(2)+el_pmap(3)+el_pmap(4))/4
            Pref_appl(2) = (el_pmap(5)+el_pmap(6)+el_pmap(7)+el_pmap(8))/4
!
         case(110)
            Pref_appl(1) = (el_pmap(1) + el_pmap(5))/2
            Pref_appl(2) = (el_pmap(2) + el_pmap(6))/2
            Pref_appl(3) = (el_pmap(3) + el_pmap(7))/2
            Pref_appl(4) = (el_pmap(4) + el_pmap(8))/2
!
         case(101)
            Pref_appl(1) = (el_pmap(1) + el_pmap(4))/2
            Pref_appl(2) = (el_pmap(2) + el_pmap(3))/2
            Pref_appl(3) = (el_pmap(6) + el_pmap(7))/2
            Pref_appl(4) = (el_pmap(5) + el_pmap(8))/2
!
         case(011)
            Pref_appl(1) = (el_pmap(1) + el_pmap(2))/2
            Pref_appl(2) = (el_pmap(3) + el_pmap(4))/2
            Pref_appl(3) = (el_pmap(7) + el_pmap(8))/2
            Pref_appl(4) = (el_pmap(5) + el_pmap(6))/2
!
         case(111)
            do i = 1,8
               Pref_appl(i) = el_pmap(i)
            enddo
      end select
!
   endif
!
end subroutine subson_one_irregularity_map