! Routines:
!  - AffineSegment
!  - AffineQuadrilateral
!  - AffineTriangle
!  - AffineHexahedron
!  - AffineTetrahedron
!  - AffinePrism
!  - AffinePyramid
!----------------------------------------------------------------------
!  Define relevant affine coordinates for each element
!  People who want different master element geometries only need to
!  modify this file.
!----------------------------------------------------------------------
   subroutine AffineSegment(Xi, Mu,DMu)
!
      implicit none
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: Mu(0:1),DMu(0:1)
!
!  ...Define affine coordinates and their gradients
      Mu(0)  = 1.d0-Xi; Mu(1)  = Xi
      DMu(0) = -1.d0;   DMu(1) = 1.d0
!
   end subroutine AffineSegment
!----------------------------------------------------------------------
   subroutine AffineQuadrilateral(Xi, Mu,DMu)
!
      implicit none
      double precision, intent(in)  :: Xi(2)
      double precision, intent(out) :: Mu(0:1,1:2),DMu(1:2,0:1,1:2)
!
!  ...Define affine coordinates
      Mu(0,1) = 1.d0-Xi(1); Mu(1,1) = Xi(1)
      Mu(0,2) = 1.d0-Xi(2); Mu(1,2) = Xi(2)
!  ...and their gradients
      DMu(1:2,0:1,1:2) = 0.d0
      DMu(1,0,1) = -1.d0; DMu(1,1,1) = 1.d0
      DMu(2,0,2) = -1.d0; DMu(2,1,2) = 1.d0
!
   end subroutine AffineQuadrilateral
!----------------------------------------------------------------------
   subroutine AffineTriangle(X, Nu,DNu)
!
      implicit none
      double precision, intent(in)  :: X(2)
      double precision, intent(out) :: Nu(0:2),DNu(1:2,0:2)
!
!  ...Define affine coordinates
      Nu(0) = 1.d0-X(1)-X(2); Nu(1) = X(1); Nu(2) = X(2)
!  ...and their gradients
      DNu(1,0) = -1.d0;  DNu(1,1) = 1.d0; DNu(1,2) = 0.d0
      DNu(2,0) = -1.d0;  DNu(2,1) = 0.d0; DNu(2,2) = 1.d0
!
   end subroutine AffineTriangle
!----------------------------------------------------------------------
   subroutine AffineHexahedron(Xi, Mu,DMu)
!
      implicit none
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
!
!  ...Define affine coordinates
      Mu(0,1) = 1.d0-Xi(1); Mu(1,1) = Xi(1)
      Mu(0,2) = 1.d0-Xi(2); Mu(1,2) = Xi(2)
      Mu(0,3) = 1.d0-Xi(3); Mu(1,3) = Xi(3)
!  ...and their gradients
      DMu(1:3,0:1,1:3) = 0.d0
      DMu(1,0,1) = -1.d0; DMu(1,1,1) = 1.d0
      DMu(2,0,2) = -1.d0; DMu(2,1,2) = 1.d0
      DMu(3,0,3) = -1.d0; DMu(3,1,3) = 1.d0
!
   end subroutine AffineHexahedron
!----------------------------------------------------------------------
   subroutine AffineTetrahedron(X, Lam,DLam)
!
      implicit none
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: Lam(0:3),DLam(1:3,0:3)
!
!  ...Define affine coordinates
      Lam(0) = 1.d0-X(1)-X(2)-X(3); Lam(1) = X(1)
      Lam(2) = X(2);                Lam(3) = X(3)
!  ...and their gradients
      DLam(1:3,0:3) = 0.d0
      DLam(1,0) = -1.d0;  DLam(1,1) = 1.d0
      DLam(2,0) = -1.d0;  DLam(2,2) = 1.d0
      DLam(3,0) = -1.d0;  DLam(3,3) = 1.d0
!
   end subroutine AffineTetrahedron
!----------------------------------------------------------------------
   subroutine AffinePrism(X, Mu,DMu,Nu,DNu)
!
      implicit none
      double precision, intent(in)  :: X(3)
      double precision, intent(out) :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(out) :: Nu(0:2),DNu(1:3,0:2)
!
!  ...Define triangle affine coordinates
      Nu(0) = 1.d0-X(1)-X(2); Nu(1) = X(1); Nu(2) = X(2)
!  ...and their gradients
      DNu(1:3,0:2) = 0.d0
      DNu(1,0) = -1.d0;  DNu(1,1) = 1.d0
      DNu(2,0) = -1.d0;  DNu(2,2) = 1.d0
!  ...Define segment affine coordinates
      Mu(0) = 1.d0-X(3); Mu(1) = X(3)
!  ...and their gradients
      DMu(1:3,0:1) = 0.d0
      DMu(3,0) = -1.d0;  DMu(3,1) = 1.d0
!
   end subroutine AffinePrism
!----------------------------------------------------------------------
   subroutine AffinePyramid(Xi, Lam,DLam,Mu,DMu,Nu,DNu,MuZ,DMuZ)
!
      implicit none
      double precision, intent(in)  :: Xi(3)
      double precision, intent(out) :: Lam(1:5),DLam(1:3,1:5)
      double precision, intent(out) :: Mu(0:1,1:2),DMu(1:3,0:1,1:2)
      double precision, intent(out) :: Nu(0:2,1:2),DNu(1:3,0:2,1:2)
      double precision, intent(out) :: MuZ(0:1),DMuZ(1:3,0:1)
      double precision :: zeta, eps
!
      eps = 1.0d-12
!
      zeta = Xi(3)
!  ...First define the two sets of triangle affine coordinates
      Nu(0,1) = 1.d0-Xi(1)-zeta; Nu(1,1) = Xi(1); Nu(2,1) = zeta
      Nu(0,2) = 1.d0-Xi(2)-zeta; Nu(1,2) = Xi(2); Nu(2,2) = zeta
!  ...and their gradients
      DNu(1:3,0:2,1:2) = 0.d0
      DNu(1,0,1) = -1.d0;  DNu(1,1,1) = 1.d0
      DNu(3,0,1) = -1.d0;  DNu(3,2,1) = 1.d0
      DNu(2,0,2) = -1.d0;  DNu(2,1,2) = 1.d0
      DNu(3,0,2) = -1.d0;  DNu(3,2,2) = 1.d0
!
!  ...Define segment affine coordinates over the height
      MuZ(0) = 1.d0-zeta; MuZ(1) = zeta
!  ...Don't divide by zero
      if (abs(MuZ(0)) < eps)  then
        MuZ(0) = 1.d0-eps; MuZ(1) = eps
      endif
!  ...and their gradients
      DMuZ(1:3,0:1) = 0.d0
      DMuZ(3,0) = -1.d0;  DMuZ(3,1) = 1.d0
!
!  ...Next the two sets of scaled segment affine coordinates
      Mu(0,1) = 1.d0-Xi(1)/MuZ(0); Mu(1,1) = Xi(1)/MuZ(0)
      Mu(0,2) = 1.d0-Xi(2)/MuZ(0); Mu(1,2) = Xi(2)/MuZ(0)
!  ...and their gradients
      DMu(1:3,0:1,1:2) = 0.d0
      DMu(1,0,1) = -1.d0/MuZ(0);     DMu(1,1,1) = 1.d0/MuZ(0)
      DMu(3,0,1) = -Xi(1)/MuZ(0)**2; DMu(3,1,1) = Xi(1)/MuZ(0)**2
      DMu(2,0,2) = -1.d0/MuZ(0);     DMu(2,1,2) = 1.d0/MuZ(0)
      DMu(3,0,2) = -Xi(2)/MuZ(0)**2; DMu(3,1,2) = Xi(2)/MuZ(0)**2
!
!  ...Finally the pyramid affine-like coordinates
      Lam(1) = Nu(0,1)*Mu(0,2)
      Lam(2) = Nu(0,2)*Mu(1,1)
      Lam(3) = Nu(1,1)*Mu(1,2)
      Lam(4) = Nu(1,2)*Mu(0,1)
      Lam(5) = zeta
!  ...and their gradients
      DLam(1:3,1) = Nu(0,1)*DMu(1:3,0,2)+DNu(1:3,0,1)*Mu(0,2)
      DLam(1:3,2) = Nu(0,2)*DMu(1:3,1,1)+DNu(1:3,0,2)*Mu(1,1)
      DLam(1:3,3) = Nu(1,1)*DMu(1:3,1,2)+DNu(1:3,1,1)*Mu(1,2)
      DLam(1:3,4) = Nu(1,2)*DMu(1:3,0,1)+DNu(1:3,1,2)*Mu(0,1)
      DLam(1,5) = 0.d0; DLam(2,5) = 0.d0; DLam(3,5) = 1.d0
!
   end subroutine AffinePyramid
!----------------------------------------------------------------------
