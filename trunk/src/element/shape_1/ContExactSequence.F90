! Routines:
!  - shape1DH
!  - shape1DQ
!  - shape2DH
!  - shape2DE
!  - shape2DV
!  - shape2DQ
!  - shape3DH
!  - shape3DE
!  - shape3DV
!  - shape3DQ
!
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
!     routine name      - shape1DH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - evaluate 1D H1 shape functions
!
!     arguments:
!
!     in:
!       Xi              - master element coordinates
!       Nord            - polynomial order of edge node (H1 sense)
!
!     out:
!       NrdofH          - number of shape functions
!       ShapH           - values of shape functions
!       GradH           - gradients of shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape1DH(Xi,Nord, NrdofH,ShapH,GradH)
!
      use parameters , only : MAXP
      use node_types
      use physics
!
      implicit none
      integer, intent(in)  :: Nord
      integer, intent(out) :: NrdofH
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(MAXP+1)
      double precision, intent(out) :: GradH(MAXP+1)
!
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = Nord
!
      call checkorder(SEGM,CONTIN,norder,MAXP, nsize)
      call shape1DHSeg(Xi,Nord,nsize, NrdofH,ShapH,GradH)
!
   end subroutine shape1DH
!
!----------------------------------------------------------------------
!                                 1D L2
!----------------------------------------------------------------------
!
!     routine name      - shape1DQ
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - evaluate 1D L2 shape functions
!
!     arguments:
!
!     in:
!       Xi              - master element coordinates
!       Nord            - polynomial order of edge node (H1 sense)
!
!     out:
!       NrdofQ          - number of shape functions
!       ShapQ           - values of shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape1DQ(Xi,Nord, NrdofQ,ShapQ)
!
      use parameters , only : MAXP
      use node_types
      use physics
!
      implicit none
      double precision, intent(in)  :: Xi
      integer         , intent(in)  :: Nord
      integer         , intent(out) :: NrdofQ
      double precision, intent(out) :: ShapQ(MAXP)
!
      integer :: nsize(2),norder(19)
      norder = 0; norder(1) = Nord
!
      call checkorder(SEGM,DISCON,norder,MAXP, nsize)
      call shape1DQSeg(Xi,Nord,nsize, NrdofQ,ShapQ)
!
   end subroutine shape1DQ
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
!     routine name      - shape2DH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine returns values of a 2D element
!                         H1 shape functions and their derivatives
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!       NoriE           - edge orientations
!
!     out:
!       NrdofH          - number of dof
!       ShapH           - values of the shape functions at the point
!       GradH           - gradients of the shape functions
!
!-----------------------------------------------------------------------
!
   subroutine shape2DH(NType,Xi,Norder,NoriE, NrdofH,ShapH,GradH)
!
      use parameters , only : MAXP,MAXquadH
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(in)  :: NoriE(4)
      integer         , intent(out) :: NrdofH
      double precision, intent(out) :: ShapH(MAXquadH)
      double precision, intent(out) :: GradH(2,MAXquadH)
!
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
!
      ShapH = 0.d0; GradH = 0.d0
!
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
!
   end subroutine shape2DH
!
!----------------------------------------------------------------------
!                                2D Hcurl
!----------------------------------------------------------------------
!
!     routine name      - shape2DE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine returns values of a 2D element
!                         H(curl) shape functions and their curl
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!       NoriE           - edge orientations
!     out:
!       NrdofE          - number of dof
!       ShapE           - values of the shape functions at the point
!       CurlE           - curl of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape2DE(NType,Xi,Norder,NoriE, NrdofE,ShapE,CurlE)
!
      use parameters , only : MAXP,MAXquadE
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(in)  :: NoriE(4)
      integer         , intent(out) :: NrdofE
      double precision, intent(out) :: ShapE(2,MAXquadE)
      double precision, intent(out) :: CurlE(MAXquadE)
!
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
!
      ShapE = 0.d0; CurlE = 0.d0
!
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
!
   end subroutine shape2DE
!
!----------------------------------------------------------------------
!                          2D Hdiv (rotated Hcurl)
!----------------------------------------------------------------------
!
!     routine name      - shape2DV
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine returns values of a 2D element
!                         H(div) shape functions and their divergence
!                         NOTE: only relevant in 2D problems!!
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!       NoriE           - edge orientations
!
!     out:
!       NrdofV          - number of dof
!       ShapV           - values of the shape functions at the point
!       DivV            - divergences of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape2DV(NType,Xi,Norder,NoriE, NrdofV,ShapV,DivV)
!
      use parameters , only : MAXP,MAXquadV
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(in)  :: NoriE(4)
      integer         , intent(out) :: NrdofV
      double precision, intent(out) :: ShapV(2,MAXquadV)
      double precision, intent(out) :: DivV(MAXquadV)
!
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
!
      ShapV = 0.d0; DivV = 0.d0
!
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
!
   end subroutine shape2DV
!
!----------------------------------------------------------------------
!                                 2D L2
!----------------------------------------------------------------------
!
!     routine name      - shape2DQ
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine returns values of a 2D element
!                         L2 shape functions
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!
!     out:
!       NrdofQ          - number of dof
!       ShapQ           - values of the shape functions at the point
!
!----------------------------------------------------------------------
!
   subroutine shape2DQ(NType,Xi,Norder, NrdofQ,ShapQ)
!
      use parameters , only : MAXP,MAXquadQ
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(2)
      integer         , intent(in)  :: Norder(5)
      integer         , intent(out) :: NrdofQ
      double precision, intent(out) :: ShapQ(MAXquadQ)
!
      integer :: nsize(2),norder3D(19)
      norder3D = 0; norder3D(1:5) = Norder
!
      ShapQ = 0.d0
!
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
!
   end subroutine shape2DQ
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
!     routine name      - shape3DH
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates H1 shape functions for 3D
!                         elements of various types
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!       NoriE           - edge orientations
!       NoriF           - face orientations
!
!     out:
!       NrdofH          - number of the element shape functions
!       ShapH           - values of shape functions
!       GradH           - gradients of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape3DH(NType,Xi,Norder,NoriE,NoriF, &
                          NrdofH,ShapH,GradH)
!
      use parameters , only : MAXP,MAXbrickH
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(in)  :: NoriE(12)
      integer         , intent(in)  :: NoriF(6)
      integer         , intent(out) :: NrdofH
      double precision, intent(out) :: ShapH(MAXbrickH)
      double precision, intent(out) :: GradH(3,MAXbrickH)
!
      integer :: nsize(2)
!
      ShapH = 0.d0; GradH = 0.d0
!
      call checkorder(NType,CONTIN,Norder,MAXP, nsize)
      select case(NType)
      case(BRIC,MDLB)
        call shape3DHHexa(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofH,ShapH,GradH)
      case(TETR,MDLN)
        call shape3DHTet(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofH,ShapH,GradH)
      case(PRIS,MDLP)
        call shape3DHPris(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofH,ShapH,GradH)
      case(PYRA,MDLD)
        call shape3DHPyra(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofH,ShapH,GradH)
      case default
        write(*,*)'shape3DH: Type = ', S_Type(Ntype)
        stop
      end select
!
   end subroutine shape3DH
!
!----------------------------------------------------------------------
!                                3D Hcurl
!----------------------------------------------------------------------
!
!     routine name      - shape3DE
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates H(curl) shape functions for
!                         3D elements of various types
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!       NoriE           - edge orientations
!       NoriF           - face orientations
!
!     out:
!       NrdofE          - number of the element shape functions
!       ShapE           - shape functions
!       CurlE           - curls of shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape3DE(NType,Xi,Norder,NoriE,NoriF, &
                          NrdofE,ShapE,CurlE)
!
      use parameters , only : MAXP,MAXbrickE
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(in)  :: NoriE(12)
      integer         , intent(in)  :: NoriF(6)
      integer         , intent(out) :: NrdofE
      double precision, intent(out) :: ShapE(3,MAXbrickE)
      double precision, intent(out) :: CurlE(3,MAXbrickE)
!
      integer :: nsize(2)
!
      ShapE = 0.d0; CurlE = 0.d0
!
      call checkorder(NType,TANGEN,Norder,MAXP, nsize)
      select case(NType)
      case(BRIC,MDLB)
        call shape3DEHexa(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofE,ShapE,CurlE)
      case(TETR,MDLN)
        call shape3DETet(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofE,ShapE,CurlE)
      case(PRIS,MDLP)
        call shape3DEPris(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofE,ShapE,CurlE)
      case(PYRA,MDLD)
        call shape3DEPyra(Xi,Norder,NoriE,NoriF,nsize, &
                                        NrdofE,ShapE,CurlE)
      case default
        write(*,*)'shape3DE: Type = ', S_Type(Ntype)
        stop 1
      end select
!
   end subroutine shape3DE
!
!----------------------------------------------------------------------
!                                3D Hdiv
!----------------------------------------------------------------------
!
!     routine name      - shape3DV
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates H(div) shape functions for
!                         3D elements of various types
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order for the nodes (H1 sense)
!       NoriF           - face orientations
!
!     out:
!       NrdofV          - number of the element shape functions
!       ShapV           - shape functions
!       DivV            - divergence of shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape3DV(NType,Xi,Norder,NoriF, NrdofV,ShapV,DivV)
!
      use parameters , only : MAXP,MAXbrickV
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(in)  :: NoriF(6)
      integer         , intent(out) :: NrdofV
      double precision, intent(out) :: ShapV(3,MAXbrickV)
      double precision, intent(out) :: DivV(MAXbrickV)
!
      integer :: nsize(2)
!
      ShapV = 0.d0; DivV = 0.d0
!
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
!
   end subroutine shape3DV
!
!----------------------------------------------------------------------
!                                 3D L2
!----------------------------------------------------------------------
!
!     routine name      - shape3DQ
!
!----------------------------------------------------------------------
!
!     latest revision:  - Feb 2023
!
!     purpose:          - routine evaluates L2 shape functions for 3D
!                         elements of various types
!
!     arguments:
!
!     in:
!       NType           - element type
!       Xi              - master element coordinates
!       Norder          - polynomial order of interior node (H1 sense)
!
!     out:
!       NrdofQ          - number of the element shape functions
!       ShapQ           - values of shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape3DQ(NType,Xi,Norder, NrdofQ,ShapQ)
!
      use parameters , only : MAXP,MAXbrickQ
      use node_types
      use physics
!
      implicit none
      integer         , intent(in)  :: Ntype
      double precision, intent(in)  :: Xi(3)
      integer         , intent(in)  :: Norder(19)
      integer         , intent(out) :: NrdofQ
      double precision, intent(out) :: ShapQ(MAXbrickQ)
!
      integer :: nsize(2)
!
      ShapQ = 0.d0
!
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
!
   end subroutine shape3DQ
!
