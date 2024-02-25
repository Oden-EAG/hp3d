!---------------------------------------------------------------------------------
!> @brief  Given order of an element middle node, determine the implied order for
!!         its edges and faces accounting for the orientation of the nodes
!!
!! @param[in]  Mdle      - element middle node
!! @param[in]  Norientl  - element nodes orientations
!! @param[out] Norder    - order for the element nodes
!!
!! @date Feb 2023
!---------------------------------------------------------------------------------
subroutine element_order(Mdle,Norientl, Norder)
!
      use element_data
      use data_structure3D
!
      implicit none
!
      integer,                intent(in)  :: Mdle
      integer, dimension(27), intent(in)  :: Norientl
      integer, dimension(19), intent(out) :: Norder
!
      integer :: ntype,iprint,nord,nord1,nord2,nord3,norda,i,ie,ifc,nordh,nordv
!
      select case(Mdle)
         case(469)   ; iprint=0
         case default; iprint=0
      end select
!
      ntype = NODES(Mdle)%ntype
      nord = NODES(Mdle)%order
!
      if (iprint.eq.1) then
        write(*,*) 'element_order: Mdle, type, nord = ', Mdle, S_Type(ntype), nord
      endif
!
      select case(ntype)
!  ...tet
      case(MDLN)
        Norder(1:11) = nord
!
!  ...prism
      case(MDLP)
        call decode(nord, nord1,nord2)
!
!  .....horizontal edges
        Norder(1:6) = nord1
!
!  .....vertical edges
        Norder(7:9) = nord2
!
!  .....triangular faces
        Norder(10:11) = nord1
!
!  .....rectangular faces
        do i=1,3
          select case(Norientl(17+i))
          case(0,2,5,7)
            Norder(11+i) = nord1*10+nord2
          case(1,3,4,6)
            Norder(11+i) = nord2*10+nord1
          case default
            write(*,*) 'element_order: i,Norientl(17+i) = ',i,Norientl(17+i)
            stop 1
          end select
        enddo
        Norder(15) = nord
!
!  ...hexa
      case(MDLB)
        call decode(nord, norda,nord3)
        call decode(norda, nord1,nord2)
!
!  .....loop through edges
        do ie=1,12
          select case(ie)
!
!  .......edges parallel to x1 axis
          case(1,3,5,7); Norder(ie) = nord1
!
!  .......edges parallel to x2 axis
          case(2,4,6,8); Norder(ie) = nord2
!
!  .......edges parallel to x3 axis
          case(9,10,11,12); Norder(ie) = nord3
          end select
        enddo
!
!  .....loop through faces
        do ifc=1,6
          select case(ifc)
          case(1,2); nordh=nord1; nordv=nord2
          case(3,5); nordh=nord1; nordv=nord3
          case(4,6); nordh=nord2; nordv=nord3
          end select
          select case(Norientl(20+ifc))
          case(0,2,5,7)
            Norder(12+ifc) = nordh*10+nordv
          case(1,3,4,6)
            Norder(12+ifc) = nordv*10+nordh
          case default
            write(*,*) 'element_order: ifc,Norientl(20+ifc) = ',ifc,Norientl(20+ifc)
            stop 1
          end select
        enddo
        Norder(19) = nord
!
!  ...pyramid
      case(MDLD)
        Norder(1:8) = nord
        Norder(9) = nord*10+nord
        Norder(10:14) = nord
      end select
      if (iprint.eq.1) then
        write(*,*) 'element_order: norder = ',norder
        call pause
      endif
!
end subroutine element_order
!
!---------------------------------------------------------------------------------
!> @brief  Given order of element nodes, use the min rule in reverse to
!!         determine the necessary new order for the middle node
!!
!!
!! @param[in]  Mdle        - element middle node
!! @param[in]  Norientl    - element nodes orientations
!! @param[in]  Norder      - order for the element nodes
!! @param[out] Nordm       - order for the middle node
!!
!! @date Feb 2023
!---------------------------------------------------------------------------------
subroutine element_middle_node_order(Mdle,Norientl,Norder, Nordm)
!
      use element_data
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Norientl(27)
      integer, intent(in)  :: Norder(19)
      integer, intent(out) :: Nordm
!
      integer :: ntype,nord,nord1,nord2,nord3,norda,i,ie,ifc,nordh,nordv
!
      ntype = NODES(Mdle)%ntype
      nord = NODES(Mdle)%order
!
      select case(ntype)
!
!  ...tet
      case(MDLN)
        do i=1,10
          nord = max(nord,Norder(i))
        enddo
        Nordm = nord
!
!  ...prism
      case(MDLP)
        call decode(nord, nord1,nord2)
!
!  .....horizontal edges
        do i=1,6
          nord1 = max(nord1,Norder(i))
        enddo
!
!  .....vertical edges
        do i=7,9
          nord2 = max(nord2,Norder(i))
        enddo
!
!  .....triangular faces
        do i=10,11
          nord1 = max(nord1,Norder(i))
        enddo
!
!  .....rectangular faces
        do i=1,3
          call decode(Norder(11+i), nordh,nordv)
          select case(Norientl(17+i))
          case(0,2,5,7)
            nord1 = max(nord1,nordh)
            nord2 = max(nord2,nordv)
          case(1,3,4,6)
            nord1 = max(nord1,nordv)
            nord2 = max(nord2,nordh)
          case default
            write(*,*) 'element_order: i,Norientl(17+i) = ',i,Norientl(17+i)
            stop 1
          end select
        enddo
        Nordm = nord1*10+nord2
!
!  ...hexa
      case(MDLB)
        call decode(nord, norda,nord3)
        call decode(norda, nord1,nord2)
!
!  .....loop through edges
        do ie=1,12
          select case(ie)
!
!  .......edges parallel to x1 axis
          case(1,3,5,7)
            nord1 = max(nord1,Norder(ie))
!
!  .......edges parallel to x2 axis
          case(2,4,6,8)
            nord2 = max(nord2,Norder(ie))
!
!  .......edges parallel to x3 axis
          case(9,10,11,12)
            nord3 = max(nord3,Norder(ie))
          end select
        enddo
!
!  .....loop through faces
        do ifc=1,6
          call decode(Norder(12+ifc), nordh,nordv)
          select case(ifc)
          case(1,2)
            select case(Norientl(20+ifc))
            case(0,2,5,7)
              nord1 = max(nord1,nordh)
              nord2 = max(nord2,nordv)
            case(1,3,4,6)
              nord1 = max(nord1,nordv)
              nord2 = max(nord2,nordh)
            end select
          case(3,5)
            select case(Norientl(20+ifc))
            case(0,2,5,7)
              nord1 = max(nord1,nordh)
              nord3 = max(nord3,nordv)
            case(1,3,4,6)
              nord1 = max(nord1,nordv)
              nord3 = max(nord3,nordh)
            end select
          case(4,6)
            select case(Norientl(20+ifc))
            case(0,2,5,7)
              nord2 = max(nord2,nordh)
              nord3 = max(nord3,nordv)
            case(1,3,4,6)
              nord2 = max(nord2,nordv)
              nord3 = max(nord3,nordh)
            end select
          end select
        enddo
        Nordm = nord1*100+nord2*10+nord3
!
!  ...pyramid
      case(MDLD)
        do i=1,8
          nord = max(nord,Norder(i))
        enddo
        call decode(Norder(9), nordh,nordv)
        nord = max(nord,nordh)
        nord = max(nord,nordv)
        do i=1,13
          nord = max(nord,Norder(i))
        enddo
        Nordm = nord
      end select
!
end subroutine element_middle_node_order
