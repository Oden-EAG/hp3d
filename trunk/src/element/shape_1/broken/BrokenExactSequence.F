! Routines:
!  - shape1HH
!  - shape1QQ
!  - shape2HH
!  - shape2EE
!  - shape2VV
!  - shape2QQ
!  - shape3HH
!  - shape3EE
!  - shape3VV
!  - shape3QQ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                            1D: H1--->L2
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                                 1D H1
!----------------------------------------------------------------------
!
!     routine name      - shape1HH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - evaluate 1D BROKEN H1 shape functions
!
!     arguments:
!
!     in:
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofH          - number of the element shape functions
!       ShapH           - values of shape functions
!       GradH           - gradients of the shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape1HH(Xi,NordM, NrdofH,ShapH,GradH)
!
      use parametersDPG , only : MAXPP
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(MAXPP+1)
      double precision, intent(out) :: GradH(MAXPP+1)
!
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = NordM
!
      call checkorder(SEGM,CONTIN,norder,MAXPP, nsize)
      call shape1DHBrokenSeg(Xi,NordM,nsize, NrdofH,ShapH,GradH)
!
      end subroutine shape1HH
!
!----------------------------------------------------------------------
!                                 1D L2
!----------------------------------------------------------------------
!
!     routine name      - shape1QQ
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - evaluate 1D BROKEN L2 shape functions
!
!     arguments:
!
!     in:
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofQ          - number of the element shape functions
!       ShapQ           - values of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape1QQ(Xi,NordM, NrdofQ,ShapQ)
!
      use parametersDPG , only : MAXPP
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapQ(MAXPP)
!
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = NordM
!
      call checkorder(SEGM,DISCON,norder,MAXPP, nsize)
      call shape1DQBrokenSeg(Xi,NordM,nsize, NrdofQ,ShapQ)
!
      end subroutine shape1QQ
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                        2D: H1--->Hcurl--->L2
!                            H1--->Hdiv --->L2 (rotated)
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                                 2D H1
!----------------------------------------------------------------------
!
!     routine name      - shape2HH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN H1 shape
!                         functions for 2D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofH          - number of the element shape functions
!       ShapH           - values of shape functions
!       GradH           - gradients of the shape functions
!
!-----------------------------------------------------------------------
!
      subroutine shape2HH(Ntype,Xi,NordM, NrdofH,ShapH,GradH)
!
      use parametersDPG , only : MAXPP,MAXquadHH
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapH(MAXquadHH)
      double precision, intent(out) :: GradH(2,MAXquadHH)
!
      integer :: nsize(2),norder(19)
      norder = 0
!
      select case(Ntype)
      case(TRIA,MDLT)
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Ntype,CONTIN,norder,MAXPP, nsize)
        call shape2DHBrokenTri(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case(QUAD,MDLQ,RECT)
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Ntype,CONTIN,norder,MAXPP, nsize)
        call shape2DHBrokenQuad(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape2HH: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape2HH
!
!----------------------------------------------------------------------
!                                2D Hcurl
!----------------------------------------------------------------------
!
!     routine name      - shape2EE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN H(curl) shape
!                         functions for 2D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofE          - number of the element shape functions
!       ShapE           - shape functions
!       CurlE           - curls of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape2EE(Ntype,Xi,NordM, NrdofE,ShapE,CurlE)
!
      use parametersDPG , only : MAXPP,MAXquadEE
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofE
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapE(2,MAXquadEE)
      double precision, intent(out) :: CurlE(MAXquadEE)
!
      integer :: nsize(2),norder(19)
      norder = 0
!
      select case(Ntype)
      case(TRIA,MDLT)
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Ntype,TANGEN,norder,MAXPP, nsize)
        call shape2DEBrokenTri(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case(QUAD,MDLQ,RECT)
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Ntype,TANGEN,norder,MAXPP, nsize)
        call shape2DEBrokenQuad(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape2EE: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape2EE
!
!----------------------------------------------------------------------
!                          2D Hdiv (rotated Hcurl)
!----------------------------------------------------------------------
!
!     routine name      - shape2VV
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN H(div) shape
!                         functions for 2D elements of various types
!                         NOTE: only relevant in 2D problems!!
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofV          - number of the element shape functions
!       ShapV           - shape functions
!       DivV            - divergence of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape2VV(Ntype,Xi,NordM, NrdofV,ShapV,DivV)
!
      use parametersDPG , only : MAXPP,MAXquadVV
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofV
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapV(2,MAXquadVV)
      double precision, intent(out) :: DivV(MAXquadVV)
!
      integer :: nsize(2),norder(19)
      norder = 0;
!
      select case(Ntype)
      case(TRIA,MDLT)
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Ntype,NORMAL,norder,MAXPP, nsize)
        call shape2DVBrokenTri(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case(QUAD,MDLQ,RECT)
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Ntype,NORMAL,norder,MAXPP, nsize)
        call shape2DVBrokenQuad(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape2VV: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape2VV
!
!----------------------------------------------------------------------
!                                 2D L2
!----------------------------------------------------------------------
!
!     routine name      - shape2QQ
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN L2 shape
!                         functions for 2D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofQ          - number of the element shape functions
!       ShapQ           - values of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape2QQ(Ntype,Xi,NordM, NrdofQ,ShapQ)
!
      use parametersDPG , only : MAXPP,MAXquadQQ
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: ShapQ(MAXquadQQ)
!
      integer :: nsize(2),norder(19)
      norder = 0
!
      select case(Ntype)
      case(TRIA,MDLT)
        norder(1:3)=1
        norder(4)=NordM
        call checkorder(Ntype,DISCON,norder,MAXPP, nsize)
        call shape2DQBrokenTri(Xi,NordM,nsize, NrdofQ,ShapQ)
      case(QUAD,MDLQ,RECT)
        norder(1:4)=1
        norder(5)=NordM
        call checkorder(Ntype,DISCON,norder,MAXPP, nsize)
        call shape2DQBrokenQuad(Xi,NordM,nsize, NrdofQ,ShapQ)
      case default
        write(*,*)'shape2QQ: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape2QQ
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!                      3D: H1--->Hcurl--->Hdiv--->L2
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!                                 3D H1
!----------------------------------------------------------------------
!
!     routine name      - shape3HH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN H1 shape
!                         functions for 3D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofH          - number of the element shape functions
!       ShapH           - values of shape functions
!       GradH           - gradients of the shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape3HH(Ntype,Xi,NordM, NrdofH,ShapH,GradH)
!
      use parameters    , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickHH
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapH(MAXbrickHH)
      double precision, intent(out) :: GradH(3,MAXbrickHH)
!
      integer :: nsize(2),norder(19),nordF
      norder = 0
!
      select case(Ntype)
      case(BRIC,MDLB)
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Ntype,CONTIN,norder,MAXPP, nsize)
        call shape3DHBrokenHexa(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case(TETR,MDLN)
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Ntype,CONTIN,norder,MAXPP, nsize)
        call shape3DHBrokenTet(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case(PRIS,MDLP)
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Ntype,CONTIN,norder,MAXPP, nsize)
        call shape3DHBrokenPris(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case(PYRA,MDLD)
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Ntype,CONTIN,norder,MAXPP, nsize)
        call shape3DHBrokenPyra(Xi,NordM,nsize, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape3HH: Type = ', S_Type(Ntype)
        stop
      end select
!
      end subroutine shape3HH
!
!----------------------------------------------------------------------
!                                3D Hcurl
!----------------------------------------------------------------------
!
!     routine name      - shape3EE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN H(curl) shape
!                         functions for 3D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofE          - number of the element shape functions
!       ShapE           - shape functions
!       CurlE           - curls of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape3EE(Ntype,Xi,NordM, NrdofE,ShapE,CurlE)
!
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickEE
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofE
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapE(3,MAXbrickEE)
      double precision, intent(out) :: CurlE(3,MAXbrickEE)
!
      integer :: nsize(2),norder(19),nordF
      norder = 0
!
      select case(Ntype)
      case(BRIC,MDLB)
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Ntype,TANGEN,norder,MAXPP, nsize)
        call shape3DEBrokenHexa(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case(TETR,MDLN)
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Ntype,TANGEN,norder,MAXPP, nsize)
        call shape3DEBrokenTet(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case(PRIS,MDLP)
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Ntype,TANGEN,norder,MAXPP, nsize)
        call shape3DEBrokenPris(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case(PYRA,MDLD)
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Ntype,TANGEN,norder,MAXPP, nsize)
        call shape3DEBrokenPyra(Xi,NordM,nsize, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape3EE: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape3EE
!
!----------------------------------------------------------------------
!                                3D Hdiv
!----------------------------------------------------------------------
!
!     routine name      - shape3VV
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN H(div) shape
!                         functions for 3D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofV          - number of the element shape functions
!       ShapV           - shape functions
!       DivV            - divergence of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape3VV(Ntype,Xi,NordM, NrdofV,ShapV,DivV)
!
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickVV
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofV
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapV(3,MAXbrickVV)
      double precision, intent(out) :: DivV(MAXbrickVV)
!
      integer :: nsize(2),norder(19),nordF
      norder = 0
!
      select case(Ntype)
      case(BRIC,MDLB)
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Ntype,NORMAL,norder,MAXPP, nsize)
        call shape3DVBrokenHexa(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case(TETR,MDLN)
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Ntype,NORMAL,norder,MAXPP, nsize)
        call shape3DVBrokenTet(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case(PRIS,MDLP)
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Ntype,NORMAL,norder,MAXPP, nsize)
        call shape3DVBrokenPris(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case(PYRA,MDLD)
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Ntype,NORMAL,norder,MAXPP, nsize)
        call shape3DVBrokenPyra(Xi,NordM,nsize, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape3VV: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape3VV
!
!----------------------------------------------------------------------
!                                 3D L2
!----------------------------------------------------------------------
!
!     routine name      - shape3QQ
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates BROKEN L2 shape
!                         functions for 3D elements of various types
!
!     arguments:
!
!     in:
!       Ntype           - element type
!       Xi              - master element coordinates
!       NordM           - polynomial order for middle node (H1 sense)
!
!     out:
!       NrdofQ          - number of the element shape functions
!       ShapQ           - values of shape functions
!
!----------------------------------------------------------------------
!
      subroutine shape3QQ(Ntype,Xi,NordM, NrdofQ,ShapQ)
!
      use parameters , only : MODORDER
      use parametersDPG , only : MAXPP,MAXbrickQQ
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: NordM
      integer, intent(out) :: NrdofQ
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: ShapQ(MAXbrickQQ)
!
      integer :: nsize(2),norder(19),nordF
      norder = 0
!
      select case(Ntype)
      case(BRIC,MDLB)
        norder(1:12)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(13:18)=nordF
        norder(19)=NordM
        call checkorder(Ntype,DISCON,norder,MAXPP, nsize)
        call shape3DQBrokenHexa(Xi,NordM,nsize, NrdofQ,ShapQ)
      case(TETR,MDLN)
        norder(1:10)=1
        norder(11)=NordM
        call checkorder(Ntype,DISCON,norder,MAXPP, nsize)
        call shape3DQBrokenTet(Xi,NordM,nsize, NrdofQ,ShapQ)
      case(PRIS,MDLP)
        norder(1:11)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(12:14)=nordF
        norder(15)=NordM
        call checkorder(Ntype,DISCON,norder,MAXPP, nsize)
        call shape3DQBrokenPris(Xi,NordM,nsize, NrdofQ,ShapQ)
      case(PYRA,MDLD)
        norder(1:8)=1
        call encod((/1,1/),MODORDER,2, nordF)
        norder(9)=nordF
        norder(10:13)=1
        norder(14)=NordM
        call checkorder(Ntype,DISCON,norder,MAXPP, nsize)
        call shape3DQBrokenPyra(Xi,NordM,nsize, NrdofQ,ShapQ)
      case default
        write(*,*)'shape3QQ: Type = ', S_Type(Ntype)
        stop 1
      end select
!
      end subroutine shape3QQ
!
