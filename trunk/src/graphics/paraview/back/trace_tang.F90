!--------------------------------------------------------------------------------------      
!> Purpose : compute tangential trace of a vector field      
!!
!> @param[in ] Mdle      
!> @param[in ] Xi
!> @param[in ] Dxdxi     
!> @param[in ] Field     
!> @param[out] Field_t     
!!
!> @date Oct 14      
!--------------------------------------------------------------------------------------      
! 
      subroutine trace_tang(Mdle,Xi,Dxdxi,Field, Field_t)
!
      use data_structure3D , only : NODES      
      use element_data     , only : face_param,edge_param,nsign_param
      use control          , only : GEOM_TOL
!
      implicit none      
      integer,                 intent(in ) :: Mdle
      real*8 , dimension(3)  , intent(in ) :: Xi
      real*8 , dimension(3,3), intent(in ) :: Dxdxi
      real*8 , dimension(3)  , intent(in ) :: Field
      real*8 , dimension(3)  , intent(out) :: Field_t
!      
      character(len=4) :: etype
      integer :: ivert,iedge,iface,i,j
      real*8, dimension(2)   :: t
      real*8, dimension(3,2) :: dxidt
      real*8, dimension(3,2) :: dxdt
      real*8, dimension(3)   :: void,field_n,rn
      real*8 :: s
!      
!--------------------------------------------------------------------------------------      
!      
!     initialize
      Field_t(1:3)=0.d0
!
!     element type 
      etype=NODES(Mdle)%type
!
!     find face adjacency      
      call aux1(etype,Xi, ivert,iedge,iface,t)
!
!     VERTEX (trace does not exist)
      if (ivert /= 0)  return
!
!     EDGE
      if (iedge /= 0) then
!
!       compute unit tangent              
        call edge_param(etype,iedge,t(1), void(1:3),dxidt(1:3,1))
        dxdt(1:3,1:2)=0.d0
        do j=1,3
          dxdt(1:3,1) = dxdt(1:3,1) + dxdxi(1:3,j)*dxidt(j,1)
        enddo
        call normalize(dxdt(1:3,1))        
!
!       tangential component        
        call scalar_product(Field,dxdt(1:3,1), s) 
        Field_t = s*dxdt(1:3,1)
!
        return
!        
      endif
!      
!     FACE      
      if (iface /= 0) then
!      
!       compute unit normal              
        call face_param(etype,iface,t, xi,dxidt)
        dxdt(1:3,1:2)=0.d0
        do i=1,2
          do j=1,3
            dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
          enddo
        enddo
        call cross_product(dxdt(1:3,1),dxdt(1:3,2), rn(1:3))
        call normalize(rn(1:3))
        rn(1:3) = rn(1:3)*nsign_param(etype,iface)
!
!       normal component        
        call scalar_product(Field,rn, s) 
        field_n = s*rn
!        
!       tangential component
        Field_t = Field - field_n
!
!       check
        call scalar_product(Field_t,rn, s)
        if (abs(s) > GEOM_TOL) then
          write(*,*)'trace_tang : FAIL! iface = ',iface
        endif
!        
      endif
!
!
      endsubroutine trace_tang

!
!> Purpose : classifies master element point as vertex,edge,or face

      subroutine aux1(Etype,Xi, Ivert,Iedge,Iface,T)
!
      use control , only : GEOM_TOL      
!
      implicit none      
      character(len=4)   , intent(in ) :: Etype
      real*8,dimension(3), intent(in ) :: Xi
      integer            , intent(out) :: Ivert,Iedge,Iface
      real*8,dimension(2), intent(out) :: T
      real*8 :: s

      select case(Etype)
      case('mdlp','pris')
      case('mdlb','bric')
      case('mdln','tetr')
!
!     initialize              
      Ivert=0 ; Iedge=0 ; Iface=0 ; T=0.d0
!
!     1st vertex
      call norm(Xi - (/0.d0,0.d0,0.d0/), s)
      if (s < GEOM_TOL) then
        Ivert=1
        return
      endif
!
!     2nd vertex
      call norm(Xi - (/1.d0,0.d0,0.d0/), s)
      if (s < GEOM_TOL) then
        Ivert=2
        return
      endif
!
!     3rd vertex
      call norm(Xi - (/0.d0,1.d0,0.d0/), s)
      if (s < GEOM_TOL) then
        Ivert=3
        return
      endif
!
!     4th vertex
      call norm(Xi - (/0.d0,0.d0,1.d0/), s)
      if (s < GEOM_TOL) then 
        Ivert=4
        return
      endif
!
!     1st edge      
      if (abs(Xi(2)) < GEOM_TOL .AND. &
          abs(Xi(3)) < GEOM_TOL       ) then
        Iedge=1
        T(1)=Xi(1)
        return
      endif
!
!     2nd edge
      if (abs(Xi(1)+Xi(2)-1.d0) < GEOM_TOL .AND. &
          abs(Xi(3))            < GEOM_TOL       ) then
        Iedge=2
        T(1)=1.d0-Xi(1)
        return
      endif
!
!     3rd edge      
      if (abs(Xi(1)) < GEOM_TOL .AND. &
          abs(Xi(3)) < GEOM_TOL       ) then
        Iedge=3
        T(1)=Xi(2)
        return
      endif
!
!     4th edge      
      if (abs(Xi(1)) < GEOM_TOL .AND. & 
          abs(Xi(2)) < GEOM_TOL       ) then
        Iedge=4
        T(1)=Xi(3)
        return
      endif
! 
!     5th edge
      if (abs(Xi(1)+Xi(3)-1.d0) < GEOM_TOL .AND. &
          abs(Xi(2)           ) < GEOM_TOL       ) then
        Iedge=5
        T(1)=Xi(3)
        return
      endif
!      
!     6th edge
      if (abs(Xi(2)+Xi(3)-1.d0) < GEOM_TOL .AND. &
          abs(Xi(1)           ) < GEOM_TOL       ) then
        Iedge=6
        T(1)=Xi(3)
        return
      endif
!
!     1st face
      if (abs(Xi(3)) < GEOM_TOL) then
        Iface=1
        T(1)=Xi(1)
        T(2)=Xi(2)
        return
      endif
!
!     2nd face
      if (abs(Xi(2)) < GEOM_TOL) then
        Iface=2 
        T(1)=Xi(1)
        T(2)=Xi(3)
        return
      endif
!
!     3rd face
      if (abs(Xi(1)+Xi(2)+Xi(3)-1.d0) < GEOM_TOL) then
        Iface=3 
        T(1)=Xi(2)
        T(2)=Xi(3)
        return
      endif
!
!     4th face
      if (abs(Xi(1)) < GEOM_TOL) then
        Iface=4
        T(1)=Xi(2)
        T(2)=Xi(3)
        return
      endif
!
      case('mdld','pyra')
      endselect
!
!      
      endsubroutine aux1
