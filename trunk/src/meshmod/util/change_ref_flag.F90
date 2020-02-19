!-----------------------------------------------------------------------------------
!> Purpose : change face refinement flag wrt orientation
!!
!> @param[in]  How   - 'l2g' : local to global ; 'g2l' - global to local
!> @param[in]  Type  - face type
!> @param[in]  Kref  - face refinement flag (tria : 1,2,3,4 ; rect : 11,10,01 )
!> @param[in]  Nort  - face orientation
!> @param[out] Krefm - modified face refinement flag
!
!> rev@Mar 13
!-----------------------------------------------------------------------------------
subroutine change_ref_flag(How,Type,Kref,Nort, Krefm)
!      
      implicit none
      character(len=3), intent(in) :: How
!      
      character(len=4), intent(in) :: Type
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
      integer :: iprint
!-----------------------------------------------------------------------------------
!    
      iprint = 0
    
      select case(How)
!===================================================================================      
!  Local -> Global                                                                 |
!===================================================================================      
      case('l2g')
        select case(Type)
!
        case('mdlt')
          select case(Kref)
          case(1)     ; Krefm=1
          case(2,3,4) ; Krefm=1+loc_to_glob(Kref-1,Nort)
          endselect
!
        case('mdlq')
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
        select case(Type)
!
        case('mdlt')
          select case(Kref)
          case(1)     ; Krefm=1
          case(2,3,4) ; Krefm=1+glob_to_loc(Kref-1,Nort)
          endselect
!
        case('mdlq')
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
!  ...printing
#if DEBUG_FLAG
      if (iprint.eq.1) then
         write(*,7001) How,Type,Kref,Nort,Krefm
7001     format('change_ref_flag: How,Type,Kref,Nort,Krefm = ', &
          a3,2x,a5,2x,i2,2x,i2,5x,i2)
      endif
#endif
!
!
end subroutine change_ref_flag
