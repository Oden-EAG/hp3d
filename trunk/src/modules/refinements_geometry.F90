!----------------------------------------------------------------------
!   latest revision    - Jun 20
!
!   purpose            - module defines geometry of different
!                        h-refinements
!----------------------------------------------------------------------

      module refinements_geometry
!
!  ...1D segment
      real(8), parameter, dimension (2,2) :: XVERT_medg =   &
      reshape(                                              &
      (/ 0.d0,.5d0, .5d0,1.d0/)                             &
      ,(/2,2/) )
!
!  ...quad
      real(8), parameter, dimension (2,4,4) :: XVERT_mdlq11 = &
      reshape(                                                &
      (/ 0.d0,0.d0, .5d0,0.d0, .5d0,.5d0, 0.d0,.5d0,          &
         .5d0,0.d0, 1.d0,0.d0, 1.d0,.5d0, .5d0,.5d0,          &
         .5d0,.5d0, 1.d0,.5d0, 1.d0,1.d0, .5d0,1.d0,          &
         0.d0,.5d0, .5d0,.5d0, .5d0,1.d0, 0.d0,1.d0/)         &
      ,(/2,4,4/) )
!
!  ...brick    
      real(8), parameter, dimension (3,8,8) :: XVERT_mdlb111 =             & 
      reshape(                                                             &
      (/0.d0,0.d0,0.d0, .5d0,0.d0,0.d0, .5d0,.5d0,0.d0, 0.d0,.5d0,0.d0,    &
        0.d0,0.d0,.5d0, .5d0,0.d0,.5d0, .5d0,.5d0,.5d0, 0.d0,.5d0,.5d0,    &
        .5d0,0.d0,0.d0, 1.d0,0.d0,0.d0, 1.d0,.5d0,0.d0, .5d0,.5d0,0.d0,    &
        .5d0,0.d0,.5d0, 1.d0,0.d0,.5d0, 1.d0,.5d0,.5d0, .5d0,.5d0,.5d0,    &
        .5d0,.5d0,0.d0, 1.d0,.5d0,0.d0, 1.d0,1.d0,0.d0, .5d0,1.d0,0.d0,    &
        .5d0,.5d0,.5d0, 1.d0,.5d0,.5d0, 1.d0,1.d0,.5d0, .5d0,1.d0,.5d0,    &
        0.d0,.5d0,0.d0, .5d0,.5d0,0.d0, .5d0,1.d0,0.d0, 0.d0,1.d0,0.d0,    &
        0.d0,.5d0,.5d0, .5d0,.5d0,.5d0, .5d0,1.d0,.5d0, 0.d0,1.d0,.5d0,    &
        0.d0,0.d0,.5d0, .5d0,0.d0,.5d0, .5d0,.5d0,.5d0, 0.d0,.5d0,.5d0,    &
        0.d0,0.d0,1.d0, .5d0,0.d0,1.d0, .5d0,.5d0,1.d0, 0.d0,.5d0,1.d0,    &
        .5d0,0.d0,.5d0, 1.d0,0.d0,.5d0, 1.d0,.5d0,.5d0, .5d0,.5d0,.5d0,    &
        .5d0,0.d0,1.d0, 1.d0,0.d0,1.d0, 1.d0,.5d0,1.d0, .5d0,.5d0,1.d0,    &
        .5d0,.5d0,.5d0, 1.d0,.5d0,.5d0, 1.d0,1.d0,.5d0, .5d0,1.d0,.5d0,    &
        .5d0,.5d0,1.d0, 1.d0,.5d0,1.d0, 1.d0,1.d0,1.d0, .5d0,1.d0,1.d0,    &
        0.d0,.5d0,.5d0, .5d0,.5d0,.5d0, .5d0,1.d0,.5d0, 0.d0,1.d0,.5d0,    &
        0.d0,.5d0,1.d0, .5d0,.5d0,1.d0, .5d0,1.d0,1.d0, 0.d0,1.d0,1.d0/)   &
      ,(/3,8,8/) )
!
      real(8), parameter, dimension (3,8,4) :: XVERT_mdlb011 =             &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .1d1,.0d0,.0d0, .1d1,.5d0,.0d0, .0d0,.5d0,.0d0,    &
        .0d0,.0d0,.5d0, .1d1,.0d0,.5d0, .1d1,.5d0,.5d0, .0d0,.5d0,.5d0,    &
        .0d0,.5d0,.0d0, .1d1,.5d0,.0d0, .1d1,.1d1,.0d0, .0d0,.1d1,.0d0,    &
        .0d0,.5d0,.5d0, .1d1,.5d0,.5d0, .1d1,.1d1,.5d0, .0d0,.1d1,.5d0,    &
        .0d0,.5d0,.5d0, .1d1,.5d0,.5d0, .1d1,.1d1,.5d0, .0d0,.1d1,.5d0,    &
        .0d0,.5d0,.1d1, .1d1,.5d0,.1d1, .1d1,.1d1,.1d1, .0d0,.1d1,.1d1,    &
        .0d0,.0d0,.5d0, .1d1,.0d0,.5d0, .1d1,.5d0,.5d0, .0d0,.5d0,.5d0,    &
        .0d0,.0d0,.1d1, .1d1,.0d0,.1d1, .1d1,.5d0,.1d1, .0d0,.5d0,.1d1/)   &
        ,(/3,8,4/) )
!
      real(8), parameter, dimension (3,8,4) :: XVERT_mdlb101 =             &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .5d0,.1d1,.0d0, .0d0,.1d1,.0d0,    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .5d0,.1d1,.5d0, .0d0,.1d1,.5d0,    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .1d1,.1d1,.0d0, .5d0,.1d1,.0d0,    &
        .5d0,.0d0,.5d0, .1d1,.0d0,.5d0, .1d1,.1d1,.5d0, .5d0,.1d1,.5d0,    &
        .5d0,.0d0,.5d0, .1d1,.0d0,.5d0, .1d1,.1d1,.5d0, .5d0,.1d1,.5d0,    &
        .5d0,.0d0,.1d1, .1d1,.0d0,.1d1, .1d1,.1d1,.1d1, .5d0,.1d1,.1d1,    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .5d0,.1d1,.5d0, .0d0,.1d1,.5d0,    &
        .0d0,.0d0,.1d1, .5d0,.0d0,.1d1, .5d0,.1d1,.1d1, .0d0,.1d1,.1d1/)   &
        ,(/3,8,4/) )
!
      real(8), parameter, dimension (3,8,4) :: XVERT_mdlb110 =             &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .5d0,.5d0,.0d0, .0d0,.5d0,.0d0,    &
        .0d0,.0d0,.1d1, .5d0,.0d0,.1d1, .5d0,.5d0,.1d1, .0d0,.5d0,.1d1,    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .1d1,.5d0,.0d0, .5d0,.5d0,.0d0,    &
        .5d0,.0d0,.1d1, .1d1,.0d0,.1d1, .1d1,.5d0,.1d1, .5d0,.5d0,.1d1,    &
        .5d0,.5d0,.0d0, .1d1,.5d0,.0d0, .1d1,.1d1,.0d0, .5d0,.1d1,.0d0,    &
        .5d0,.5d0,.1d1, .1d1,.5d0,.1d1, .1d1,.1d1,.1d1, .5d0,.1d1,.1d1,    &
        .0d0,.5d0,.0d0, .5d0,.5d0,.0d0, .5d0,.1d1,.0d0, .0d0,.1d1,.0d0,    &
        .0d0,.5d0,.1d1, .5d0,.5d0,.1d1, .5d0,.1d1,.1d1, .0d0,.1d1,.1d1/)   &
        ,(/3,8,4/) )
!
      real(8), parameter, dimension (3,8,2) :: XVERT_mdlb100 =             &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .5d0,.1d1,.0d0, .0d0,.1d1,.0d0,    &
        .0d0,.0d0,.1d1, .5d0,.0d0,.1d1, .5d0,.1d1,.1d1, .0d0,.1d1,.1d1,    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .1d1,.1d1,.0d0, .5d0,.1d1,.0d0,    &
        .5d0,.0d0,.1d1, .1d1,.0d0,.1d1, .1d1,.1d1,.1d1, .5d0,.1d1,.1d1/)   &
        ,(/3,8,2/) )
!
      real(8), parameter, dimension (3,8,2) :: XVERT_mdlb010 =             &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .1d1,.0d0,.0d0, .1d1,.5d0,.0d0, .0d0,.5d0,.0d0,    &
        .0d0,.0d0,.1d1, .1d1,.0d0,.1d1, .1d1,.5d0,.1d1, .0d0,.5d0,.1d1,    &
        .0d0,.5d0,.0d0, .1d1,.5d0,.0d0, .1d1,.1d1,.0d0, .0d0,.1d1,.0d0,    &
        .0d0,.5d0,.1d1, .1d1,.5d0,.1d1, .1d1,.1d1,.1d1, .0d0,.1d1,.1d1/)   &
        ,(/3,8,2/) )
!
      real(8), parameter, dimension (3,8,2) :: XVERT_mdlb001 =             &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .1d1,.0d0,.0d0, .1d1,.1d1,.0d0, .0d0,.1d1,.0d0,    &
        .0d0,.0d0,.5d0, .1d1,.0d0,.5d0, .1d1,.1d1,.5d0, .0d0,.1d1,.5d0,    &
        .0d0,.0d0,.5d0, .1d1,.0d0,.5d0, .1d1,.1d1,.5d0, .0d0,.1d1,.5d0,    &
        .0d0,.0d0,.1d1, .1d1,.0d0,.1d1, .1d1,.1d1,.1d1, .0d0,.1d1,.1d1/)   &
        ,(/3,8,2/) )
!
!  ...tetrahedron
      real(8), parameter, dimension (3,4,8) :: XVERT_mdln11 =              &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .0d0,.5d0,.0d0, .0d0,.0d0,.5d0,    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .5d0,.5d0,.0d0, .5d0,.0d0,.5d0,    &
        .0d0,.5d0,.0d0, .5d0,.5d0,.0d0, .0d0,.1d1,.0d0, .0d0,.5d0,.5d0,    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .0d0,.5d0,.5d0, .0d0,.0d0,.1d1,    &
        .5d0,.0d0,.0d0, .0d0,.5d0,.5d0, .0d0,.0d0,.5d0, .5d0,.0d0,.5d0,    &
        .5d0,.0d0,.0d0, .0d0,.5d0,.5d0, .5d0,.0d0,.5d0, .5d0,.5d0,.0d0,    &
        .5d0,.0d0,.0d0, .0d0,.5d0,.5d0, .5d0,.5d0,.0d0, .0d0,.5d0,.0d0,    &
        .5d0,.0d0,.0d0, .0d0,.5d0,.5d0, .0d0,.5d0,.0d0, .0d0,.0d0,.5d0/)   &
        ,(/3,4,8/) )
!
      real(8), parameter, dimension (3,4,8) :: XVERT_mdln12 =              &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .0d0,.5d0,.0d0, .0d0,.0d0,.5d0,    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .5d0,.5d0,.0d0, .5d0,.0d0,.5d0,    &
        .0d0,.5d0,.0d0, .5d0,.5d0,.0d0, .0d0,.1d1,.0d0, .0d0,.5d0,.5d0,    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .0d0,.5d0,.5d0, .0d0,.0d0,.1d1,    &
        .0d0,.5d0,.0d0, .5d0,.0d0,.5d0, .5d0,.0d0,.0d0, .5d0,.5d0,.0d0,    &
        .0d0,.5d0,.0d0, .5d0,.0d0,.5d0, .5d0,.5d0,.0d0, .0d0,.5d0,.5d0,    &
        .0d0,.5d0,.0d0, .5d0,.0d0,.5d0, .0d0,.5d0,.5d0, .0d0,.0d0,.5d0,    &
        .0d0,.5d0,.0d0, .5d0,.0d0,.5d0, .0d0,.0d0,.5d0, .5d0,.0d0,.0d0/)   &
        ,(/3,4,8/) )
!
      real(8), parameter, dimension (3,4,8) :: XVERT_mdln13 =              &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .0d0,.5d0,.0d0, .0d0,.0d0,.5d0,    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .5d0,.5d0,.0d0, .5d0,.0d0,.5d0,    &
        .0d0,.5d0,.0d0, .5d0,.5d0,.0d0, .0d0,.1d1,.0d0, .0d0,.5d0,.5d0,    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .0d0,.5d0,.5d0, .0d0,.0d0,.1d1,    &
        .0d0,.0d0,.5d0, .5d0,.5d0,.0d0, .5d0,.0d0,.5d0, .5d0,.0d0,.0d0,    &
        .0d0,.0d0,.5d0, .5d0,.5d0,.0d0, .5d0,.0d0,.0d0, .0d0,.5d0,.0d0,    &
        .0d0,.0d0,.5d0, .5d0,.5d0,.0d0, .0d0,.5d0,.0d0, .0d0,.5d0,.5d0,    &
        .0d0,.0d0,.5d0, .5d0,.5d0,.0d0, .0d0,.5d0,.5d0, .5d0,.0d0,.5d0/)   &
        ,(/3,4,8/) )
!
!  ...prism
      real(8), parameter, dimension (3,6,8) :: XVERT_mdlp11 =              &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .0d0,.5d0,.0d0,                    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .0d0,.5d0,.5d0,                    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .5d0,.5d0,.0d0,                    &
        .5d0,.0d0,.5d0, .1d1,.0d0,.5d0, .5d0,.5d0,.5d0,                    &
        .0d0,.5d0,.0d0, .5d0,.5d0,.0d0, .0d0,.1d1,.0d0,                    &
        .0d0,.5d0,.5d0, .5d0,.5d0,.5d0, .0d0,.1d1,.5d0,                    &
        .5d0,.5d0,.0d0, .0d0,.5d0,.0d0, .5d0,.0d0,.0d0,                    &
        .5d0,.5d0,.5d0, .0d0,.5d0,.5d0, .5d0,.0d0,.5d0,                    &
        .0d0,.0d0,.5d0, .5d0,.0d0,.5d0, .0d0,.5d0,.5d0,                    &
        .0d0,.0d0,.1d1, .5d0,.0d0,.1d1, .0d0,.5d0,.1d1,                    &
        .5d0,.0d0,.5d0, .1d1,.0d0,.5d0, .5d0,.5d0,.5d0,                    &
        .5d0,.0d0,.1d1, .1d1,.0d0,.1d1, .5d0,.5d0,.1d1,                    &
        .0d0,.5d0,.5d0, .5d0,.5d0,.5d0, .0d0,.1d1,.5d0,                    &
        .0d0,.5d0,.1d1, .5d0,.5d0,.1d1, .0d0,.1d1,.1d1,                    &
        .5d0,.5d0,.5d0, .0d0,.5d0,.5d0, .5d0,.0d0,.5d0,                    &
        .5d0,.5d0,.1d1, .0d0,.5d0,.1d1, .5d0,.0d0,.1d1/)                   &
        ,(/3,6,8/) )
!
      real(8), parameter, dimension (3,6,4) :: XVERT_mdlp10 =              &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .5d0,.0d0,.0d0, .0d0,.5d0,.0d0,                    &
        .0d0,.0d0,.1d1, .5d0,.0d0,.1d1, .0d0,.5d0,.1d1,                    &
        .5d0,.0d0,.0d0, .1d1,.0d0,.0d0, .5d0,.5d0,.0d0,                    &
        .5d0,.0d0,.1d1, .1d1,.0d0,.1d1, .5d0,.5d0,.1d1,                    &
        .0d0,.5d0,.0d0, .5d0,.5d0,.0d0, .0d0,.1d1,.0d0,                    &
        .0d0,.5d0,.1d1, .5d0,.5d0,.1d1, .0d0,.1d1,.1d1,                    &
        .5d0,.5d0,.0d0, .0d0,.5d0,.0d0, .5d0,.0d0,.0d0,                    &
        .5d0,.5d0,.1d1, .0d0,.5d0,.1d1, .5d0,.0d0,.1d1/)                   &
        ,(/3,6,4/) )
!
      real(8), parameter, dimension (3,6,2) :: XVERT_mdlp01 =              &
      reshape(                                                             &
      (/.0d0,.0d0,.0d0, .1d1,.0d0,.0d0, .0d0,.1d1,.0d0,                    &
        .0d0,.0d0,.5d0, .1d1,.0d0,.5d0, .0d0,.1d1,.5d0,                    &
        .0d0,.0d0,.5d0, .1d1,.0d0,.5d0, .0d0,.1d1,.5d0,                    &
        .0d0,.0d0,.1d1, .1d1,.0d0,.1d1, .0d0,.1d1,.1d1/)                   &
        ,(/3,6,2/) )
!
      contains
!
      subroutine get_son_coord(Type,Kref,Is, Nvert,Xvert)
      use parameters, only: NDIMEN
!
      character(len=4),        intent(in)  :: Type
      integer,                 intent(in)  :: Kref,Is
      integer,                 intent(out) :: Nvert
      real*8, dimension (NDIMEN,8), intent(out) :: Xvert
!
      integer :: iprint
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'get_son_coord: Type,Kref,Is = ',Type,Kref,Is
      endif
!
      select case(Type)
      case('mdlb')
        Nvert=8
        select case(Kref)
        case(111); Xvert(1:3,1:8) = XVERT_mdlb111(1:3,1:8,Is)
        case(011); Xvert(1:3,1:8) = XVERT_mdlb011(1:3,1:8,Is)
        case(101); Xvert(1:3,1:8) = XVERT_mdlb101(1:3,1:8,Is)
        case(110); Xvert(1:3,1:8) = XVERT_mdlb110(1:3,1:8,Is)
        case(100); Xvert(1:3,1:8) = XVERT_mdlb100(1:3,1:8,Is)
        case(010); Xvert(1:3,1:8) = XVERT_mdlb010(1:3,1:8,Is)
        case(001); Xvert(1:3,1:8) = XVERT_mdlb001(1:3,1:8,Is)
        case default
          write(*,7010) Kref; stop 1
 7010     format(' get_son_coord: UNFINISHED, Kref = ', i10)
        end select
      case('mdln')
        Nvert=4
        select case(Kref)
        case(11); Xvert(1:3,1:4) = XVERT_mdln11(1:3,1:4,Is)
        case(12); Xvert(1:3,1:4) = XVERT_mdln12(1:3,1:4,Is)
        case(13); Xvert(1:3,1:4) = XVERT_mdln13(1:3,1:4,Is)
        case default
          write(*,7010) Kref; stop 2
        end select
      case('mdlp')
        Nvert=6
        select case(Kref)
        case(11); Xvert(1:3,1:6) = XVERT_mdlp11(1:3,1:6,Is)
        case(10); Xvert(1:3,1:6) = XVERT_mdlp10(1:3,1:6,Is)
        case(01); Xvert(1:3,1:6) = XVERT_mdlp01(1:3,1:6,Is)
        case default
          write(*,7010) Kref; stop 3
        end select
      case default
        write(*,7020) Type ; stop 4
 7020   format(' get_son_coord: UNFINISHED, Type = ', a10)
      end select
      if (iprint.eq.1) then
        write(*,*) 'get_son_coord: Xvert = ',Xvert
        call pause
      endif
!
      end subroutine get_son_coord
!
      end module refinements_geometry
