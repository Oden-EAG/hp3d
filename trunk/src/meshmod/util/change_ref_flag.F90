!-----------------------------------------------------------------------------------
!> @brief      change face refinement flag wrt orientation
!!
!> @param[in]  How   - 'l2g' : local to global ; 'g2l' - global to local
!> @param[in]  Ntype - face type
!> @param[in]  Kref  - face refinement flag (tria : 1,2,3,4 ; rect : 11,10,01 )
!> @param[in]  Nort  - face orientation
!> @param[out] Krefm - modified face refinement flag
!!
!> @date       Feb 2023
!-----------------------------------------------------------------------------------
subroutine change_ref_flag(How,Ntype,Kref,Nort, Krefm)
!
      use node_types
      implicit none
      character(len=3), intent(in) :: How
!
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: Kref, Nort
      integer, intent(out) :: Krefm
!
! ** Locals for triangle
      integer, parameter, dimension(1:3,0:5) :: loc_to_glob  &
       = reshape( (/1,2,3, 3,1,2, 2,3,1, 1,3,2, 2,1,3, 3,2,1/), &
       (/3,6/) )

      integer, parameter, dimension(1:3,0:5) :: glob_to_loc(1:3,0:5) &
           = reshape( (/1,2,3, 2,3,1, 3,1,2, 1,3,2, 2,1,3, 3,2,1/), &
           (/3,6/) )
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!-----------------------------------------------------------------------------------
!
      select case(How)
!===================================================================================
!  Local -> Global                                                                 |
!===================================================================================
      case('l2g')
        select case(Ntype)
!
        case(MDLT)
          select case(Kref)
          case(1)     ; Krefm=1
          case(2,3,4) ; Krefm=1+loc_to_glob(Kref-1,Nort)
          endselect
!
        case(MDLQ)
          select case(Kref)
          case(11)        ; Krefm=11
          case(10,01)
            select case(Nort)
            case(0,2,5,7) ; Krefm=Kref
            case(1,3,4,6) ; Krefm=11-Kref
            endselect
          endselect
        endselect
!
!===================================================================================
!  Global -> Local                                                                 |
!===================================================================================
      case('g2l')
        select case(Ntype)
!
        case(MDLT)
          select case(Kref)
          case(1)     ; Krefm=1
          case(2,3,4) ; Krefm=1+glob_to_loc(Kref-1,Nort)
          endselect
!
        case(MDLQ)
          select case(Kref)
          case(11)        ; Krefm=11
          case(10,01)
            select case(Nort)
            case(0,2,5,7) ; Krefm=Kref
            case(1,3,4,6) ; Krefm=11-Kref
            endselect
          endselect
        endselect
!
      endselect
!
#if HP3D_DEBUG
!  ...printing
      if (iprint.eq.1) then
         write(*,7001) How,S_Type(Ntype),Kref,Nort,Krefm
7001     format('change_ref_flag: How,Type,Kref,Nort,Krefm = ', &
          a3,2x,a5,2x,i2,2x,i2,5x,i2)
      endif
#endif
!
end subroutine change_ref_flag
