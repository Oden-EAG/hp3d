c Routines:
c  - AffineSegment
c  - AffineQuadrilateral
c  - AffineTriangle
c  - AffineHexahedron
c  - AffineTetrahedron
c  - AffinePrism
c  - AffinePyramid
c----------------------------------------------------------------------
c  Define relevant affine coordinates for each element
c  People who want different master element geometries only need to
c  modify this file.
c----------------------------------------------------------------------
      subroutine AffineSegment(Xi, Mu,DMu)
c
      implicit none
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: Mu(0:1),DMu(0:1)
c
c  ...Define affine coordinates and their gradients
      Mu(0)  = 1.d0-Xi; Mu(1)  = Xi
      DMu(0) = -1.d0;   DMu(1) = 1.d0
c
      end subroutine AffineSegment
c----------------------------------------------------------------------
      subroutine AffineQuadrilateral(Xi, Mu,DMu)
c
      implicit none
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: Mu(1:2,0:1),DMu(1:2,1:2,0:1)
c
c  ...Define affine coordinates
      Mu(1,0) = 1.d0-Xi(1); Mu(1,1) = Xi(1)
      Mu(2,0) = 1.d0-Xi(2); Mu(2,1) = Xi(2)
c  ...and their gradients
      DMu(1:2,1:2,0:1) = 0.d0
      DMu(1,1,0) = -1.d0; DMu(1,1,1) = 1.d0
      DMu(2,2,0) = -1.d0; DMu(2,2,1) = 1.d0
c
      end subroutine AffineQuadrilateral
c----------------------------------------------------------------------
      subroutine AffineTriangle(X, Nu,DNu)
c
      implicit none
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: Nu(0:2),DNu(1:2,0:2)
c
c  ...Define affine coordinates
      Nu(0) = 1.d0-X(1)-X(2); Nu(1) = X(1); Nu(2) = X(2)
c  ...and their gradients
      DNu(1,0) = -1.d0;  DNu(1,1) = 1.d0; DNu(1,2) = 0.d0
      DNu(2,0) = -1.d0;  DNu(2,1) = 0.d0; DNu(2,2) = 1.d0
c
      end subroutine AffineTriangle
c----------------------------------------------------------------------
      subroutine AffineHexahedron(Xi, Mu,DMu)
c
      implicit none
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: Mu(1:3,0:1),DMu(1:3,1:3,0:1)
c
c  ...Define affine coordinates
      Mu(1,0) = 1.d0-Xi(1); Mu(1,1) = Xi(1)
      Mu(2,0) = 1.d0-Xi(2); Mu(2,1) = Xi(2)
      Mu(3,0) = 1.d0-Xi(3); Mu(3,1) = Xi(3)
c  ...and their gradients
      DMu(1:3,1:3,0:1) = 0.d0
      DMu(1,1,0) = -1.d0; DMu(1,1,1) = 1.d0
      DMu(2,2,0) = -1.d0; DMu(2,2,1) = 1.d0
      DMu(3,3,0) = -1.d0; DMu(3,3,1) = 1.d0
c
      end subroutine AffineHexahedron
c----------------------------------------------------------------------
      subroutine AffineTetrahedron(X, Lam,DLam)
c
      implicit none
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: Lam(0:3),DLam(1:3,0:3)
c
c  ...Define affine coordinates
      Lam(0) = 1.d0-X(1)-X(2)-X(3); Lam(1) = X(1)
      Lam(2) = X(2);                Lam(3) = X(3)
c  ...and their gradients
      DLam(1:3,0:3) = 0.d0
      DLam(1,0) = -1.d0;  DLam(1,1) = 1.d0
      DLam(2,0) = -1.d0;  DLam(2,2) = 1.d0
      DLam(3,0) = -1.d0;  DLam(3,3) = 1.d0
c
      end subroutine AffineTetrahedron
c----------------------------------------------------------------------
      subroutine AffinePrism(X, Mu,DMu,Nu,DNu)
c
      implicit none
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(out) :: Nu(0:2),DNu(1:3,0:2)
c
c  ...Define triangle affine coordinates
      Nu(0) = 1.d0-X(1)-X(2); Nu(1) = X(1); Nu(2) = X(2)
c  ...and their gradients
      DNu(1:3,0:2) = 0.d0
      DNu(1,0) = -1.d0;  DNu(1,1) = 1.d0
      DNu(2,0) = -1.d0;  DNu(2,2) = 1.d0
c  ...Define segment affine coordinates
      Mu(0) = 1.d0-X(3); Mu(1) = X(3)
c  ...and their gradients
      DMu(1:3,0:1) = 0.d0
      DMu(3,0) = -1.d0;  DMu(3,1) = 1.d0
c
      end subroutine AffinePrism
c----------------------------------------------------------------------
      subroutine AffinePyramid(Xi, Lam,DLam,Mu,DMu,Nu,DNu,MuZ,DMuZ)
c
      implicit none
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: Lam(1:5),DLam(1:3,1:5)
      double precision, intent(out) :: Mu(1:2,0:1),DMu(1:3,1:2,0:1)
      double precision, intent(out) :: Nu(1:2,0:2),DNu(1:3,1:2,0:2)
      double precision, intent(out) :: MuZ(0:1),DMuZ(1:3,0:1)
      double precision :: zeta, eps
c
      eps = 1.0d-12
c
      zeta = Xi(3)
c  ...First define the two sets of triangle affine coordinates
      Nu(1,0) = 1.d0-Xi(1)-zeta; Nu(1,1) = Xi(1); Nu(1,2) = zeta
      Nu(2,0) = 1.d0-Xi(2)-zeta; Nu(2,1) = Xi(2); Nu(2,2) = zeta
c  ...and their gradients
      DNu(1:3,1:2,0:2) = 0.d0
      DNu(1,1,0) = -1.d0;  DNu(1,1,1) = 1.d0
      DNu(3,1,0) = -1.d0;  DNu(3,1,2) = 1.d0
      DNu(2,2,0) = -1.d0;  DNu(2,2,1) = 1.d0
      DNu(3,2,0) = -1.d0;  DNu(3,2,2) = 1.d0
c
c  ...Define segment affine coordinates over the height
      MuZ(0) = 1.d0-zeta; MuZ(1) = zeta
c  ...Don't divide by zero
      if (abs(MuZ(0)) < eps)  then
        MuZ(0) = 1.d0-eps; MuZ(1) = eps
      endif
c  ...and their gradients
      DMuZ(1:3,0:1) = 0.d0
      DMuZ(3,0) = -1.d0;  DMuZ(3,1) = 1.d0
c
c  ...Next the two sets of scaled segment affine coordinates
      Mu(1,0) = 1.d0-Xi(1)/MuZ(0); Mu(1,1) = Xi(1)/MuZ(0)
      Mu(2,0) = 1.d0-Xi(2)/MuZ(0); Mu(2,1) = Xi(2)/MuZ(0)
c  ...and their gradients
      DMu(1:3,1:2,0:1) = 0.d0
      DMu(1,1,0) = -1.d0/MuZ(0);     DMu(1,1,1) = 1.d0/MuZ(0)
      DMu(3,1,0) = -Xi(1)/MuZ(0)**2; DMu(3,1,1) = Xi(1)/MuZ(0)**2
      DMu(2,2,0) = -1.d0/MuZ(0);     DMu(2,2,1) = 1.d0/MuZ(0)
      DMu(3,2,0) = -Xi(2)/MuZ(0)**2; DMu(3,2,1) = Xi(2)/MuZ(0)**2
c
c  ...Finally the pyramid affine-like coordinates
      Lam(1) = Nu(1,0)*Mu(2,0)
      Lam(2) = Nu(2,0)*Mu(1,1)
      Lam(3) = Nu(1,1)*Mu(2,1)
      Lam(4) = Nu(2,1)*Mu(1,0)
      Lam(5) = zeta
c  ...and their gradients
      DLam(1:3,1) = Nu(1,0)*DMu(1:3,2,0)+DNu(1:3,1,0)*Mu(2,0)
      DLam(1:3,2) = Nu(2,0)*DMu(1:3,1,1)+DNu(1:3,2,0)*Mu(1,1)
      DLam(1:3,3) = Nu(1,1)*DMu(1:3,2,1)+DNu(1:3,1,1)*Mu(2,1)
      DLam(1:3,4) = Nu(2,1)*DMu(1:3,1,0)+DNu(1:3,2,1)*Mu(1,0)
      DLam(1,5) = 0.d0; DLam(2,5) = 0.d0; DLam(3,5) = 1.d0
c
      end subroutine AffinePyramid
c----------------------------------------------------------------------
