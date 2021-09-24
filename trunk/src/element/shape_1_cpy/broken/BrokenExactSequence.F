c Routines:
c  - shape1HH
c  - shape1QQ
c  - shape2HH
c  - shape2EE
c  - shape2VV
c  - shape2QQ
c  - shape3HH
c  - shape3EE
c  - shape3VV
c  - shape3QQ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c                            1D: H1--->L2
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c                                 1D H1
c----------------------------------------------------------------------
c
c     routine name      - shape1HH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - evaluate 1D BROKEN H1 shape functions
c
c     arguments:
c
c     in:
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofH          - number of the element shape functions
c       ShapH           - values of shape functions
c       GradH           - gradients of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1HH(Xi,NordM, NrdofH,ShapH,GradH)
c
      use parametersDPG , only : MAXPP
c
      implicit none
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(MAXPP+1)
      double precision, intent(out) :: GradH(MAXPP+1)
c
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = NordM
c
      call checkorder('segm','contin',norder,MAXPP, nsize)
      call shape1DHBrokenSeg(Xi,NordM,nsize, NrdofH,ShapH,GradH)
c
      end subroutine shape1HH
c
c----------------------------------------------------------------------
c                                 1D L2
c----------------------------------------------------------------------
c
c     routine name      - shape1QQ
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - evaluate 1D BROKEN L2 shape functions
c
c     arguments:
c
c     in:
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofQ          - number of the element shape functions
c       ShapQ           - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1QQ(Xi,NordM, NrdofQ,ShapQ)
c
      use parametersDPG , only : MAXPP
c
      implicit none
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapQ(MAXPP)
c
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = NordM
c
      call checkorder('segm','discon',norder,MAXPP, nsize)
      call shape1DQBrokenSeg(Xi,NordM,nsize, NrdofQ,ShapQ)
c
      end subroutine shape1QQ
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c                        2D: H1--->Hcurl--->L2
c                            H1--->Hdiv --->L2 (rotated)
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c                                 2D H1
c----------------------------------------------------------------------
c
c     routine name      - shape2HH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN H1 shape
c                         functions for 2D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofH          - number of the element shape functions
c       ShapH           - values of shape functions
c       GradH           - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2HH(Type,Xi,NordM, NrdofH,ShapH,GradH)
c
      use parametersDPG , only : MAXPP,MAXquadHH
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapH(MAXquadHH)
      double precision, intent(out) :: GradH(2,MAXquadHH)
c
      integer :: nsize(2),norder(19)
      norder = 0
c
      select case(Type)
      case('tria','mdlt')
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Type,'contin',norder,MAXPP, nsize)
        call shape2DHBrokenTri(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case('quad','mdlq','rect')
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Type,'contin',norder,MAXPP, nsize)
        call shape2DHBrokenQuad(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape2HH: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2HH
c
c----------------------------------------------------------------------
c                                2D Hcurl
c----------------------------------------------------------------------
c
c     routine name      - shape2EE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN H(curl) shape
c                         functions for 2D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofE          - number of the element shape functions
c       ShapE           - shape functions
c       CurlE           - curls of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2EE(Type,Xi,NordM, NrdofE,ShapE,CurlE)
c
      use parametersDPG , only : MAXPP,MAXquadEE
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofE
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapE(2,MAXquadEE)
      double precision, intent(out) :: CurlE(MAXquadEE)
c
      integer :: nsize(2),norder(19)
      norder = 0
c
      select case(Type)
      case('tria','mdlt')
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Type,'tangen',norder,MAXPP, nsize)
        call shape2DEBrokenTri(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case('quad','mdlq','rect')
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Type,'tangen',norder,MAXPP, nsize)
        call shape2DEBrokenQuad(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape2EE: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2EE
c
c----------------------------------------------------------------------
c                          2D Hdiv (rotated Hcurl)
c----------------------------------------------------------------------
c
c     routine name      - shape2VV
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN H(div) shape
c                         functions for 2D elements of various types
c                         NOTE: only relevant in 2D problems!!
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofV          - number of the element shape functions
c       ShapV           - shape functions
c       DivV            - divergence of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2VV(Type,Xi,NordM, NrdofV,ShapV,DivV)
c
      use parametersDPG , only : MAXPP,MAXquadVV
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofV
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapV(2,MAXquadVV)
      double precision, intent(out) :: DivV(MAXquadVV)
c
      integer :: nsize(2),norder(19)
      norder = 0;
c
      select case(Type)
      case('tria','mdlt')
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Type,'normal',norder,MAXPP, nsize)
        call shape2DVBrokenTri(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case('quad','mdlq','rect')
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Type,'normal',norder,MAXPP, nsize)
        call shape2DVBrokenQuad(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape2VV: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2VV
c
c----------------------------------------------------------------------
c                                 2D L2
c----------------------------------------------------------------------
c
c     routine name      - shape2QQ
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN L2 shape
c                         functions for 2D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofQ          - number of the element shape functions
c       ShapQ           - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2QQ(Type,Xi,NordM, NrdofQ,ShapQ)
c
      use parametersDPG , only : MAXPP,MAXquadQQ
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapQ(MAXquadQQ)
c
      integer :: nsize(2),norder(19)
      norder = 0
c
      select case(Type)
      case('tria','mdlt')
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Type,'discon',norder,MAXPP, nsize)
        call shape2DQBrokenTri(Xi,NordM,nsize, NrdofQ,ShapQ)
      case('quad','mdlq','rect')
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Type,'discon',norder,MAXPP, nsize)
        call shape2DQBrokenQuad(Xi,NordM,nsize, NrdofQ,ShapQ)
      case default
        write(*,*)'shape2QQ: Type = ', Type
        stop 1
      end select
c
      end subroutine shape2QQ
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c                      3D: H1--->Hcurl--->Hdiv--->L2
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c                                 3D H1
c----------------------------------------------------------------------
c
c     routine name      - shape3HH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN H1 shape
c                         functions for 3D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofH          - number of the element shape functions
c       ShapH           - values of shape functions
c       GradH           - gradients of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3HH(Type,Xi,NordM, NrdofH,ShapH,GradH)
c
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickHH
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(MAXbrickHH)
      double precision, intent(out) :: GradH(3,MAXbrickHH)
c
      integer :: nsize(2),norder(19),nordF
      norder = 0
c
      select case(Type)
      case('bric','mdlb')
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Type,'contin',norder,MAXPP, nsize)
        call shape3DHBrokenHexa(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case('tetr','mdln')
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Type,'contin',norder,MAXPP, nsize)
        call shape3DHBrokenTet(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case('pris','mdlp')
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Type,'contin',norder,MAXPP, nsize)
        call shape3DHBrokenPris(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case('pyra','mdld')
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Type,'contin',norder,MAXPP, nsize)
        call shape3DHBrokenPyra(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape3HH: Type = ', Type
        stop
      end select
c
      end subroutine shape3HH
c
c----------------------------------------------------------------------
c                                3D Hcurl
c----------------------------------------------------------------------
c
c     routine name      - shape3EE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN H(curl) shape
c                         functions for 3D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofE          - number of the element shape functions
c       ShapE           - shape functions
c       CurlE           - curls of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3EE(Type,Xi,NordM, NrdofE,ShapE,CurlE)
c
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickEE
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofE
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(3,MAXbrickEE)
      double precision, intent(out) :: CurlE(3,MAXbrickEE)
c
      integer :: nsize(2),norder(19),nordF
      norder = 0
c
      select case(Type)
      case('bric','mdlb')
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Type,'tangen',norder,MAXPP, nsize)
        call shape3DEBrokenHexa(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case('tetr','mdln')
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Type,'tangen',norder,MAXPP, nsize)
        call shape3DEBrokenTet(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case('pris','mdlp')
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Type,'tangen',norder,MAXPP, nsize)
        call shape3DEBrokenPris(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case('pyra','mdld')
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Type,'tangen',norder,MAXPP, nsize)
        call shape3DEBrokenPyra(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape3EE: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3EE
c
c----------------------------------------------------------------------
c                                3D Hdiv
c----------------------------------------------------------------------
c
c     routine name      - shape3VV
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN H(div) shape
c                         functions for 3D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofV          - number of the element shape functions
c       ShapV           - shape functions
c       DivV            - divergence of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3VV(Type,Xi,NordM, NrdofV,ShapV,DivV)
c
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickVV
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofV
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(3,MAXbrickVV)
      double precision, intent(out) :: DivV(MAXbrickVV)
c
      integer :: nsize(2),norder(19),nordF
      norder = 0
c
      select case(Type)
      case('bric','mdlb')
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Type,'normal',norder,MAXPP, nsize)
        call shape3DVBrokenHexa(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case('tetr','mdln')
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Type,'normal',norder,MAXPP, nsize)
        call shape3DVBrokenTet(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case('pris','mdlp')
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Type,'normal',norder,MAXPP, nsize)
        call shape3DVBrokenPris(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case('pyra','mdld')
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Type,'normal',norder,MAXPP, nsize)
        call shape3DVBrokenPyra(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape3VV: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3VV
c
c----------------------------------------------------------------------
c                                 3D L2
c----------------------------------------------------------------------
c
c     routine name      - shape3QQ
c
c----------------------------------------------------------------------
c
c     latest revision:  - Oct 14, Apr 17
c
c     purpose:          - routine evaluates BROKEN L2 shape
c                         functions for 3D elements of various types
c
c     arguments:
c
c     in:
c       Type            - element type
c       Xi              - master element coordinates
c       NordM           - polynomial order for middle node (H1 sense)
c
c     out:
c       NrdofQ          - number of the element shape functions
c       ShapQ           - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3QQ(Type,Xi,NordM, NrdofQ,ShapQ)
c
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickQQ
c
      implicit none
      character(len=4), intent(in)  :: Type
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(MAXbrickQQ)
c
      integer :: nsize(2),norder(19),nordF
      norder = 0
c
      select case(Type)
      case('bric','mdlb')
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Type,'discon',norder,MAXPP, nsize)
        call shape3DQBrokenHexa(Xi,NordM,nsize, NrdofQ,ShapQ)
      case('tetr','mdln')
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Type,'discon',norder,MAXPP, nsize)
        call shape3DQBrokenTet(Xi,NordM,nsize, NrdofQ,ShapQ)
      case('pris','mdlp')
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Type,'discon',norder,MAXPP, nsize)
        call shape3DQBrokenPris(Xi,NordM,nsize, NrdofQ,ShapQ)
      case('pyra','mdld')
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Type,'discon',norder,MAXPP, nsize)
        call shape3DQBrokenPyra(Xi,NordM,nsize, NrdofQ,ShapQ)
      case default
        write(*,*)'shape3QQ: Type = ', Type
        stop 1
      end select
c
      end subroutine shape3QQ
c
