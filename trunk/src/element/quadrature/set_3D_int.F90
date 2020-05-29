!----------------------------------------------------------------------
!
!     routine name      - set_3Dint
!
!     DEPRECATED, kept for backward compatibility (Oct 2019)
!                 use set_3D_int instead
!----------------------------------------------------------------------
subroutine set_3Dint(Type,Norder, Nint,Xiloc,Waloc)
!
      use parameters, only : MAXP,MAX_NINT3
      implicit none
!
      character(len=4)               , intent(in)  :: Type
      integer, dimension(19)         , intent(in)  :: Norder
      integer                        , intent(out) :: Nint
      real(8), dimension(3,MAX_NINT3), intent(out) :: Xiloc
      real(8), dimension(  MAX_NINT3), intent(out) :: Waloc
!
      call set_3Dint_aux(Type,Norder,MAXP,MAX_NINT3, Nint,Xiloc,Waloc)
!
end subroutine set_3Dint
!
!----------------------------------------------------------------------
!
!     routine name: set_3D_int
!
!     latest rev  : Oct 2019
!
!     purpose     : routine sets up quadrature data for a standard
!                   3D element, accounting for different element
!                   types and orders of approximation in the element
!                   system of coordinates
!
!     arguments:
!
!     in:
!          Type         - element type
!          Norder       - order of approximation
!          Norient_face - face orientations of the element
!
!     out:
!          Nint         - number of integration points
!          Xiloc        - integration points
!          Waloc        - weights
!
!----------------------------------------------------------------------
subroutine set_3D_int(Type,Norder,Norient_face, Nint,Xiloc,Waloc)
!
      use parameters, only : MAX_NINT3
      implicit none
!
      character(len=4)               , intent(in)  :: Type
      integer, dimension(19)         , intent(in)  :: Norder
      integer, dimension(6)          , intent(in)  :: Norient_face
      integer                        , intent(out) :: Nint
      real(8), dimension(3,MAX_NINT3), intent(out) :: Xiloc
      real(8), dimension(  MAX_NINT3), intent(out) :: Waloc
!
      integer, dimension(19) :: norder_loc
!
      call find_order_loc(Type,Norder,Norient_face, norder_loc)
      call set_3Dint(Type,norder_loc, Nint,Xiloc,Waloc)
!
end subroutine set_3D_int
!
!----------------------------------------------------------------------
!
!     routine name      - set_3Dint_DPG
!
!     DEPRECATED, kept for backward compatibility (Oct 2019)
!                 use set_3D_int instead
!----------------------------------------------------------------------
subroutine set_3Dint_DPG(Type,Norder, Nint,Xiloc,Waloc)
!
      use parametersDPG, only : MAXPP,MAXNINT3ADD
      implicit none
!
      character(len=4)                 , intent(in)  :: Type
      integer, dimension(19)           , intent(in)  :: Norder
      integer                          , intent(out) :: Nint
      real(8), dimension(3,MAXNINT3ADD), intent(out) :: Xiloc
      real(8), dimension(  MAXNINT3ADD), intent(out) :: Waloc
!
      call set_3Dint_aux(Type,Norder,MAXPP,MAXNINT3ADD, Nint,Xiloc,Waloc)
!
end subroutine set_3Dint_DPG
!
!----------------------------------------------------------------------
!
!     routine name: set_3D_int_DPG
!
!     latest rev  : Oct 2019
!
!     purpose     : routine sets up quadrature data for a DPG
!                   3D element, accounting for different element
!                   types and orders of approximation in the element
!                   system of coordinates
!
!     arguments:
!
!     in:
!          Type         - element type
!          Norder       - order of approximation
!          Norient_face - face orientations of the element
!
!     out:
!          Nint         - number of integration points
!          Xiloc        - integration points
!          Waloc        - weights
!
!----------------------------------------------------------------------
subroutine set_3D_int_DPG(Type,Norder,Norient_face, Nint,Xiloc,Waloc)
!
      use parametersDPG, only : MAXNINT3ADD
      implicit none
!
      character(len=4)                 , intent(in)  :: Type
      integer, dimension(19)           , intent(in)  :: Norder
      integer, dimension(6)            , intent(in)  :: Norient_face
      integer                          , intent(out) :: Nint
      real(8), dimension(3,MAXNINT3ADD), intent(out) :: Xiloc
      real(8), dimension(  MAXNINT3ADD), intent(out) :: Waloc
!
      integer, dimension(19) :: norder_loc
!
      call find_order_loc(Type,Norder,Norient_face, norder_loc)
      call set_3Dint_DPG(Type,norder_loc, Nint,Xiloc,Waloc)
!
end subroutine set_3D_int_DPG
!
!----------------------------------------------------------------------
!
!     routine name      - set_3Dint_aux
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine sets up quadrature data for a 3D
!                         element, accouting for different element
!                         types and orders of approximation
!
!     arguments:
!
!     in:
!             Type      - element type
!             Norder    - order of approximation
!             Maxp      - maximum p
!             Max_nint3 - maximum number of integration points
!
!     out:
!             Nint      - number of integration points
!             Xiloc     - integration points
!             Waloc     - weights
!
!----------------------------------------------------------------------
subroutine set_3Dint_aux(Type,Norder,Maxp,Max_nint3, Nint,Xiloc,Waloc)
!
      use parameters,        only : MODORDER
      use control,           only : INTEGRATION
      use gauss_quadrature,  only : INITIALIZED,                &
                                    NSELECT     ,NRGAUPO,       &
                                    NSELECT_TETS,NRGAUPO_TETS,  &
                                    XIGAUSS     ,WAGAUSS,       &
                                    XIGAUS1     ,WAGAUS1,       &
                                    XIGAUSS_TETS,WAGAUSS_TETS
      implicit none
!
      character(len=4)               , intent(in)  :: Type
      integer, dimension(19)         , intent(in)  :: Norder
      integer                        , intent(in)  :: Maxp,Max_nint3
      integer                        , intent(out) :: Nint
      real(8), dimension(3,Max_nint3), intent(out) :: Xiloc
      real(8), dimension(  Max_nint3), intent(out) :: Waloc
!
      integer :: nordhv(2), nordxyz(3)
      integer :: nord,nordh,nordx,nordy,nordz
      integer :: kint,nintx,ninty,nintz
      integer :: i,iprint,l,l1,l2,l3
      real(8) :: factor
!
!----------------------------------------------------------------------
!
!  ...initialize if needed
      if (.not. INITIALIZED) call init_gauss_quadrature
!
      iprint=0
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7001) Type,Norder
 7001   format('set_3Dint_aux: Type, Norder = ',a4,2x,19i4)
      endif
#endif
!
      select case(Type)
!
!======================================================================
!  PRISM                                                              |
!======================================================================
      case('mdlp','pris')
!
!  .....determine order of approximation
        nordh=0 ; nordz=0
        do i=1,15
          select case(i)
!         horizontal edges & faces
          case(1,2,3,4,5,6, 10,11) ; nordh=max(nordh,Norder(i))
!         vertical edges
          case(7,8,9)              ; nordz=max(nordz,Norder(i))
!         vertical faces & interior
          case(12,13,14, 15)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordh=max(nordh,nordhv(1)) ; nordz=max(nordz,nordhv(2))
          endselect
        enddo
!
!  .....account for over-integration
        nordh=min(nordh+INTEGRATION,Maxp)
        nordz=min(nordz+INTEGRATION,Maxp)
!
!  .....compute number of integration points
        kint =NSELECT(nordh)
        nintx=NRGAUPO(kint)
        nintz=nordz+1
        Nint =nintx*nintz
!
!  .....compute integration points and weights
        l=0
        do l2=1,nintz
          do l1=1,nintx
            l=l+1
            Xiloc(1:2,l)=XIGAUSS(2:3,l1,kint)
            Xiloc(  3,l)=XIGAUS1(l2,nintz)
            Waloc(    l)=WAGAUSS(l1,kint)/2.d0*WAGAUS1(l2,nintz)
          enddo
        enddo
!
!======================================================================
!  BRICK                                                              |
!======================================================================
      case('mdlb','bric')
!
!  .....determine order of approximation
        nordx=0 ; nordy=0 ; nordz=0
        do i=1,19
          select case(i)
          case(1,3,5,7)
            nordx = max(nordx,Norder(i))
          case(2,4,6,8)
            nordy = max(nordy,Norder(i))
          case(9,10,11,12)
            nordz = max(nordz,Norder(i))
          case(13,14)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordx = max(nordx,nordhv(1))
            nordy = max(nordy,nordhv(2))
          case(15,17)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordx = max(nordx,nordhv(1))
            nordz = max(nordz,nordhv(2))
          case(16,18)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordy = max(nordy,nordhv(1))
            nordz = max(nordz,nordhv(2))
          case(19)
            call decod(Norder(i),MODORDER,3, nordxyz)
            nordx = max(nordx,nordxyz(1))
            nordy = max(nordy,nordxyz(2))
            nordz = max(nordz,nordxyz(3))
          end select
        enddo
!
!  .....account for overintegration
        nordx=min(nordx+INTEGRATION,Maxp)
        nordy=min(nordy+INTEGRATION,Maxp)
        nordz=min(nordz+INTEGRATION,Maxp)
!
!  .....compute number of integration points
        nintx=nordx+1
        ninty=nordy+1
        nintz=nordz+1
        Nint=nintx*ninty*nintz
!
!  .....compute integration points and weights
        l=0
        do l3=1,nintz
          do l2=1,ninty
            do l1=1,nintx
              l=l+1
              Xiloc(1,l)=XIGAUS1(l1,nintx)
              Xiloc(2,l)=XIGAUS1(l2,ninty)
              Xiloc(3,l)=XIGAUS1(l3,nintz)
              Waloc(  l)=WAGAUS1(l1,nintx)*WAGAUS1(l2,ninty)* &
                         WAGAUS1(l3,nintz)
            enddo
          enddo
        enddo
!
!======================================================================
!  TETRAHEDRON                                                        |
!======================================================================
      case('mdln','tetr')
!
!  .....determine order of approximation
        nord=0
        do i=1,11
          nord=max(nord,Norder(i))
        enddo
!
!  .....account for overintegration
        nord=min(nord+INTEGRATION,Maxp)
!
!  .....retrive number of integration points
        kint=NSELECT_TETS(nord)
        Nint=NRGAUPO_TETS(kint)
!
!  .....retrieve integration points and weights
        Xiloc(1:3,1:Nint)=XIGAUSS_TETS(1:3,1:Nint,kint)
        Waloc(    1:Nint)=WAGAUSS_TETS(    1:Nint,kint)
!
!======================================================================
!  PYRAMID                                                            |
!======================================================================
      case('mdld','pyra')
!
!  .....determine order of approximation
        nordx=0 ; nordy=0 ; nordz=0
        do i=1,14
          select case(i)
          case(1,3)
            nordx = max(nordx,Norder(i))
          case(2,4)
            nordy = max(nordy,Norder(i))
          case(5,6,7,8)
            nordz = max(nordz,Norder(i))
          case(9)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordx = max(nordx,nordhv(1))
            nordy = max(nordy,nordhv(2))
          case(10,12)
            nordx = max(nordx,Norder(i))
            nordz = max(nordz,Norder(i))
          case(11,13)
            nordy = max(nordy,Norder(i))
            nordz = max(nordz,Norder(i))
          case(14)
            nordx = max(nordx,Norder(i))
            nordy = max(nordy,Norder(i))
            nordz = max(nordz,Norder(i))
          end select
        enddo
!
!  .....account for overintegration
        nordx=min(nordx+INTEGRATION,Maxp)
        nordy=min(nordy+INTEGRATION,Maxp)
        nordz=min(nordz+INTEGRATION,Maxp)
!
!  .....compute number of integration points
        nintx=nordx+1
        ninty=nordy+1
        nintz=nordz+2
        Nint=nintx*ninty*nintz
!
!  .....compute integration points and weights
        l=0
        do l3=1,nintz
          do l2=1,ninty
            do l1=1,nintx
              l=l+1
              Xiloc(3,l)=XIGAUS1(l3,nintz)
              factor    =1.d0-Xiloc(3,l)
              Xiloc(1,l)=XIGAUS1(l1,nintx)*factor
              Xiloc(2,l)=XIGAUS1(l2,ninty)*factor
              Waloc(  l)=WAGAUS1(l1,nintx)*WAGAUS1(l2,ninty)  &
                        *WAGAUS1(l3,nintz)*factor**2
            enddo
          enddo
        enddo
!
      case default
        write(*,*) 'set_3Dint: Type = ',Type
        stop
      endselect
!
end subroutine set_3Dint_aux
!
!----------------------------------------------------------------------
!
!     routine name      - test_set_3Dint
!
!----------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine tests quadrature data for a 3D
!                         standard element
!
!     arguments         - none
!
!----------------------------------------------------------------------
subroutine test_set_3Dint
!
      use parameters
      use parametersDPG
      implicit none
!
      character(len=4) :: type
!
      integer :: norder(19)
      real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
      real(8) :: s,sexact,x,y,z,wa
      integer :: iprint,ix,iy,iz,l,nint,np
!
!-----------------------------------------------------------------------
!
      type = 'mdld'
!
      iprint=1
!
!  ...loop through orders of approximation
      do np=1,8
        norder(1:14)=np; norder(9)=np*10+np
        if (iprint.eq.1) then
          write(*,7001) type,norder
 7001     format('test_set_3Dint: type, norder = ',a4,2x,19i4)
        endif
        call set_3Dint(type,norder, nint,xiloc,waloc)
!
!  .....loop through monomials
        do iz=0,2*np
          do iy=0,2*np
            do ix=0,2*np
              if (ix+iy+iz.gt.2*np) cycle
              s=0.d0
              do l=1,nint
                x=xiloc(1,l) ; y=xiloc(2,l) ; z=xiloc(3,l)
                wa=waloc(l)
                s = s + (x**ix*y**iy*(1.d0-z)**iz)*wa
              enddo
              sexact = 1.d0/((ix+1.d0)*(iy+1.d0)*(ix+iy+iz+3.d0))
              if (iprint.eq.1) then
                write(*,7002) ix,iy,iz,s,sexact
 7002           format('test_set_3Dint: ix,iy,iz,s,sexact = ',3i2,2x,2e25.15)
              endif
              if (abs(s-sexact).gt.1.d-14) call pause
            enddo
          enddo
        enddo
      enddo
!
end subroutine test_set_3Dint
