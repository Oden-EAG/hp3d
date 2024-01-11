c Routines:
c  - shape1DH
c  - shape1DQ
c  - shape2DH
c  - shape2DE
c  - shape2DV
c  - shape2DQ
c  - shape3DH
c  - shape3DE
c  - shape3DV
c  - shape3DQ
c
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
c     routine name      - shape1DH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - evaluate 1D H1 shape functions
c
c     arguments:
c
c     in:
c       Xi              - master element coordinates
c       Nord            - polynomial order of edge node (H1 sense)
c
c     out:
c       NrdofH          - number of shape functions
c       ShapH           - values of shape functions
c       GradH           - gradients of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1DH(Xi,Nord, NrdofH,ShapH,GradH)
c
      use parameters , only : MAXP
      use node_types
      use physics
c
      implicit none
      integer, intent(in)  :: Nord
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(MAXP+1)
      double precision, intent(out) :: GradH(MAXP+1)
c
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = Nord
c
      call checkorder(SEGM,CONTIN,norder,MAXP, nsize)
      call shape1DHSeg(Xi,Nord,nsize, NrdofH,ShapH,GradH)
c
      end subroutine shape1DH
c
c----------------------------------------------------------------------
c                                 1D L2
c----------------------------------------------------------------------
c
c     routine name      - shape1DQ
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - evaluate 1D L2 shape functions
c
c     arguments:
c
c     in:
c       Xi              - master element coordinates
c       Nord            - polynomial order of edge node (H1 sense)
c
c     out:
c       NrdofQ          - number of shape functions
c       ShapQ           - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1DQ(Xi,Nord, NrdofQ,ShapQ)
c
      use parameters , only : MAXP
      use node_types
      use physics
c
      implicit none
      double precision, intent(in)  :: Xi
      integer         , intent(in)  :: Nord
      integer         , intent(out) :: NrdofQ
      double precision, intent(out) :: ShapQ(MAXP)
c
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = Nord
c
      call checkorder(SEGM,DISCON,norder,MAXP, nsize)
      call shape1DQSeg(Xi,Nord,nsize, NrdofQ,ShapQ)
c
      end subroutine shape1DQ
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
c     routine name      - shape2DH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine returns values of a 2D element
c                         H1 shape functions and their derivatives
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c       NoriE           - edge orientations
c
c     out:
c       NrdofH          - number of dof
c       ShapH           - values of the shape functions at the point
c       GradH           - gradients of the shape functions
c
c-----------------------------------------------------------------------
c
      subroutine shape2DH(NType,Xi,Norder,NoriE, NrdofH,ShapH,GradH)
c
      use parameters , only : MAXP,MAXquadH
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(in)  :: NoriE(4)
      integer         , intent(out) :: NrdofH
      double precision, intent(out) :: ShapH(MAXquadH)
      double precision, intent(out) :: GradH(2,MAXquadH)
c
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
c
      ShapH = 0.d0; GradH = 0.d0
c
      call checkorder(NType,CONTIN,norder3D,MAXP, nsize)
      select case(NType)
      case(TRIA,MDLT)
        call shape2DHTri(Xi,Norder,NoriE,nsize, NrdofH,ShapH,GradH)
      case(QUAD,MDLQ,RECT)
        call shape2DHQuad(Xi,Norder,NoriE,nsize, NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape2DH: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape2DH
c
c----------------------------------------------------------------------
c                                2D Hcurl
c----------------------------------------------------------------------
c
c     routine name      - shape2DE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine returns values of a 2D element
c                         H(curl) shape functions and their curl
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c       NoriE           - edge orientations
c     out:
c       NrdofE          - number of dof
c       ShapE           - values of the shape functions at the point
c       CurlE           - curl of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2DE(NType,Xi,Norder,NoriE, NrdofE,ShapE,CurlE)
c
      use parameters , only : MAXP,MAXquadE
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(in)  :: NoriE(4)
      integer         , intent(out) :: NrdofE
      double precision, intent(out) :: ShapE(2,MAXquadE)
      double precision, intent(out) :: CurlE(MAXquadE)
c
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
c
      ShapE = 0.d0; CurlE = 0.d0
c
      call checkorder(NType,TANGEN,norder3D,MAXP, nsize)
      select case(NType)
      case(TRIA,MDLT)
        call shape2DETri(Xi,Norder,NoriE,nsize, NrdofE,ShapE,CurlE)
      case(QUAD,MDLQ,RECT)
        call shape2DEQuad(Xi,Norder,NoriE,nsize, NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape2DE: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape2DE
c
c----------------------------------------------------------------------
c                          2D Hdiv (rotated Hcurl)
c----------------------------------------------------------------------
c
c     routine name      - shape2DV
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine returns values of a 2D element
c                         H(div) shape functions and their divergence
c                         NOTE: only relevant in 2D problems!!
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c       NoriE           - edge orientations
c
c     out:
c       NrdofV          - number of dof
c       ShapV           - values of the shape functions at the point
c       DivV            - divergences of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape2DV(NType,Xi,Norder,NoriE, NrdofV,ShapV,DivV)
c
      use parameters , only : MAXP,MAXquadV
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(in)  :: NoriE(4)
      integer         , intent(out) :: NrdofV
      double precision, intent(out) :: ShapV(2,MAXquadV)
      double precision, intent(out) :: DivV(MAXquadV)
c
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
c
      ShapV = 0.d0; DivV = 0.d0
c
      call checkorder(NType,NORMAL,norder3D,MAXP, nsize)
      select case(NType)
      case(TRIA,MDLT)
        call shape2DVTri(Xi,Norder,NoriE,nsize, NrdofV,ShapV,DivV)
      case(QUAD,MDLQ,RECT)
        call shape2DVQuad(Xi,Norder,NoriE,nsize, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape2DV: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape2DV
c
c----------------------------------------------------------------------
c                                 2D L2
c----------------------------------------------------------------------
c
c     routine name      - shape2DQ
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine returns values of a 2D element
c                         L2 shape functions
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c
c     out:
c       NrdofQ          - number of dof
c       ShapQ           - values of the shape functions at the point
c
c----------------------------------------------------------------------
c
      subroutine shape2DQ(NType,Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXP,MAXquadQ
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(out) :: NrdofQ
      double precision, intent(out) :: ShapQ(MAXquadQ)
c
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
c
      ShapQ = 0.d0
c
      call checkorder(NType,DISCON,norder3D,MAXP, nsize)
      select case(NType)
      case(TRIA,MDLT)
        call shape2DQTri(Xi,Norder,nsize, NrdofQ,ShapQ)
      case(QUAD,MDLQ,RECT)
        call shape2DQQuad(Xi,Norder,nsize, NrdofQ,ShapQ)
      case default
        write(*,*)'shape2DQ: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape2DQ
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
c     routine name      - shape3DH
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine evaluates H1 shape functions for 3D
c                         elements of various types
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c       NoriE           - edge orientations
c       NoriF           - face orientations
c
c     out:
c       NrdofH          - number of the element shape functions
c       ShapH           - values of shape functions
c       GradH           - gradients of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3DH(NType,Xi,Norder,NoriE,NoriF,
     .                    NrdofH,ShapH,GradH)
c
      use parameters , only : MAXP,MAXbrickH
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(in)  :: NoriE(12)
      integer         , intent(in)  :: NoriF(6)
      integer         , intent(out) :: NrdofH
      double precision, intent(out) :: ShapH(MAXbrickH)
      double precision, intent(out) :: GradH(3,MAXbrickH)
c
      integer :: nsize(2)
c
      ShapH = 0.d0; GradH = 0.d0
c
      call checkorder(NType,CONTIN,Norder,MAXP, nsize)
      select case(NType)
      case(BRIC,MDLB)
        call shape3DHHexa(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofH,ShapH,GradH)
      case(TETR,MDLN)
        call shape3DHTet(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofH,ShapH,GradH)
      case(PRIS,MDLP)
        call shape3DHPris(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofH,ShapH,GradH)
      case(PYRA,MDLD)
        call shape3DHPyra(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape3DH: Type = ', S_Type(Ntype)
        stop
      end select
c
      end subroutine shape3DH
c
c----------------------------------------------------------------------
c                                3D Hcurl
c----------------------------------------------------------------------
c
c     routine name      - shape3DE
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine evaluates H(curl) shape functions for
c                         3D elements of various types
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c       NoriE           - edge orientations
c       NoriF           - face orientations
c
c     out:
c       NrdofE          - number of the element shape functions
c       ShapE           - shape functions
c       CurlE           - curls of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3DE(NType,Xi,Norder,NoriE,NoriF,
     .                    NrdofE,ShapE,CurlE)
c
      use parameters , only : MAXP,MAXbrickE
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(in)  :: NoriE(12)
      integer         , intent(in)  :: NoriF(6)
      integer         , intent(out) :: NrdofE
      double precision, intent(out) :: ShapE(3,MAXbrickE)
      double precision, intent(out) :: CurlE(3,MAXbrickE)
c
      integer :: nsize(2)
c
      ShapE = 0.d0; CurlE = 0.d0
c
      call checkorder(NType,TANGEN,Norder,MAXP, nsize)
      select case(NType)
      case(BRIC,MDLB)
        call shape3DEHexa(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofE,ShapE,CurlE)
      case(TETR,MDLN)
        call shape3DETet(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofE,ShapE,CurlE)
      case(PRIS,MDLP)
        call shape3DEPris(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofE,ShapE,CurlE)
      case(PYRA,MDLD)
        call shape3DEPyra(Xi,Norder,NoriE,NoriF,nsize,
     .                                  NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape3DE: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape3DE
c
c----------------------------------------------------------------------
c                                3D Hdiv
c----------------------------------------------------------------------
c
c     routine name      - shape3DV
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine evaluates H(div) shape functions for
c                         3D elements of various types
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order for the nodes (H1 sense)
c       NoriF           - face orientations
c
c     out:
c       NrdofV          - number of the element shape functions
c       ShapV           - shape functions
c       DivV            - divergence of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3DV(NType,Xi,Norder,NoriF, NrdofV,ShapV,DivV)
c
      use parameters , only : MAXP,MAXbrickV
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(in)  :: NoriF(6)
      integer         , intent(out) :: NrdofV
      double precision, intent(out) :: ShapV(3,MAXbrickV)
      double precision, intent(out) :: DivV(MAXbrickV)
c
      integer :: nsize(2)
c
      ShapV = 0.d0; DivV = 0.d0
c
      call checkorder(NType,NORMAL,Norder,MAXP, nsize)
      select case(NType)
      case(BRIC,MDLB)
        call shape3DVHexa(Xi,Norder,NoriF,nsize, NrdofV,ShapV,DivV)
      case(TETR,MDLN)
        call shape3DVTet(Xi,Norder,NoriF,nsize, NrdofV,ShapV,DivV)
      case(PRIS,MDLP)
        call shape3DVPris(Xi,Norder,NoriF,nsize, NrdofV,ShapV,DivV)
      case(PYRA,MDLD)
        call shape3DVPyra(Xi,Norder,NoriF,nsize, NrdofV,ShapV,DivV)
      case default
        write(*,*)'shape3DV: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape3DV
c
c----------------------------------------------------------------------
c                                 3D L2
c----------------------------------------------------------------------
c
c     routine name      - shape3DQ
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine evaluates L2 shape functions for 3D
c                         elements of various types
c
c     arguments:
c
c     in:
c       NType           - element type
c       Xi              - master element coordinates
c       Norder          - polynomial order of interior node (H1 sense)
c
c     out:
c       NrdofQ          - number of the element shape functions
c       ShapQ           - values of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3DQ(NType,Xi,Norder, NrdofQ,ShapQ)
c
      use parameters , only : MAXP,MAXbrickQ
      use node_types
      use physics
c
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(out) :: NrdofQ
      double precision, intent(out) :: ShapQ(MAXbrickQ)
c
      integer :: nsize(2)
c
      ShapQ = 0.d0
c
      call checkorder(NType,DISCON,Norder,MAXP, nsize)
      select case(NType)
      case(BRIC,MDLB)
        nsize=(/MAXP,MAXbrickQ/)
        call shape3DQHexa(Xi,Norder,nsize, NrdofQ,ShapQ)
      case(TETR,MDLN)
        call shape3DQTet(Xi,Norder,nsize, NrdofQ,ShapQ)
      case(PRIS,MDLP)
        call shape3DQPris(Xi,Norder,nsize, NrdofQ,ShapQ)
      case(PYRA,MDLD)
        call shape3DQPyra(Xi,Norder,nsize, NrdofQ,ShapQ)
      case default
        write(*,*)'shape3DQ: Type = ', S_Type(Ntype)
        stop 1
      end select
c
      end subroutine shape3DQ
c
