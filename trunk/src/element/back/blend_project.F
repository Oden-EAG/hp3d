c----------------------------------------------------------------------
c
c   routine name       - blend_project_quad
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 14
c
c   purpose            - routine returns values of edge blending
c                        functions and quad to edge projections
c                        along with their derivatives
c
c   arguments :
c     in:
c        Xi            - master element coordinates of a point
c
c     out:
c        Psi           - blending function values
c        Gradpsi       - gradient of psi
c        Te            - element to edge projections
c        Gradte        - gradient of Te
c
c   required  routines -
c
c----------------------------------------------------------------------

      subroutine blend_project_quad(Xi, Psi,Gradpsi,Te,Gradte)
c
#include "syscom.blk"
c
      dimension Xi(2),Psi(4),Gradpsi(2,4),Te(4),Gradte(2,4)
c
c----------------------------------------------------------------------
c
c  ...blending functions
      Psi(1) = 1.d0-Xi(2); Gradpsi(1,1) =  0.d0; Gradpsi(2,1) = -1.d0
      Psi(2) =      Xi(1); Gradpsi(1,2) =  1.d0; Gradpsi(2,2) =  0.d0
      Psi(3) =      Xi(2); Gradpsi(1,3) =  0.d0; Gradpsi(2,3) =  1.d0
      Psi(4) = 1.d0-Xi(1); Gradpsi(1,4) = -1.d0; Gradpsi(2,4) =  0.d0
c
c  ...projections
      Te(1) = Xi(1); Gradte(1,1) = 1.d0; Gradte(2,1) = 0.d0
      Te(2) = Xi(2); Gradte(1,2) = 0.d0; Gradte(2,2) = 1.d0
      Te(3) = Xi(1); Gradte(1,3) = 1.d0; Gradte(2,3) = 0.d0
      Te(4) = Xi(2); Gradte(1,4) = 0.d0; Gradte(2,4) = 1.d0
c
c
      end subroutine blend_project_quad
c
c----------------------------------------------------------------------
c
c   routine name       - blend_project_edges_bric
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 13
c
c   purpose            - routine returns values of brick edge
c                        blending functions along with brick to edge
c                        projections with their gradients
c
c   arguments :
c     in:
c        Xi            - master element coordinates of a point
c
c     out:
c        Psie          - edge blending function values
c        Gradpsie      - gradient of psie
c        Te            - element to edge projections
c        Gradte        - gradient of Te
c
c----------------------------------------------------------------------

      subroutine blend_project_edges_bric(Xi, Psie,Gradpsie,Te,Gradte)
c
#include "syscom.blk"
c
      dimension Xi(3),
     .          Psie(12),Gradpsie(3,12),Te(12),Gradte(3,12)
      double precision z,z1
c
      iprint=0
c
c----------------------------------------------------------------------
c
      Gradpsie = 0.d0; Gradte = 0.d0
      x = Xi(1); x1 = 1.d0-x
      y = Xi(2); y1 = 1.d0-y
      z = Xi(3); z1 = 1.d0-z
c
c  ...edge blending functions;
c
c  ...edges parallel to x-axis
      Psie(1) = y1*z1; Gradpsie(2,1) = -z1; Gradpsie(3,1) = -y1
      Psie(3) = y *z1; Gradpsie(2,3) =  z1; Gradpsie(3,3) = -y
      Psie(5) = y1*z ; Gradpsie(2,5) = -z ; Gradpsie(3,5) =  y1
      Psie(7) = y *z ; Gradpsie(2,7) =  z ; Gradpsie(3,7) =  y
c
c  ...edges parallel to y-axis
      Psie(2) = x *z1; Gradpsie(1,2) =  z1; Gradpsie(3,2) = -x
      Psie(4) = x1*z1; Gradpsie(1,4) = -z1; Gradpsie(3,4) = -x1
      Psie(6) = x *z ; Gradpsie(1,6) =  z ; Gradpsie(3,6) =  x
      Psie(8) = x1*z ; Gradpsie(1,8) = -z ; Gradpsie(3,8) =  x1
c
c  ...edges parallel to z-axis
      Psie( 9) = x1*y1; Gradpsie(1, 9) = -y1; Gradpsie(2, 9) = -x1
      Psie(10) = x *y1; Gradpsie(1,10) =  y1; Gradpsie(2,10) = -x
      Psie(11) = x *y ; Gradpsie(1,11) =  y ; Gradpsie(2,11) =  x
      Psie(12) = x1*y ; Gradpsie(1,12) = -y ; Gradpsie(2,12) =  x1
c
c  ...edge projections:
c
c  ...edges parallel to the x-axis
      Te( 1) = x; Gradte(1, 1) = 1.d0
      Te( 3) = x; Gradte(1, 3) = 1.d0
      Te( 5) = x; Gradte(1, 5) = 1.d0
      Te( 7) = x; Gradte(1, 7) = 1.d0
c
c  ...edges parallel to the y-axis
      Te( 2) = y; Gradte(2, 2) = 1.d0
      Te( 4) = y; Gradte(2, 4) = 1.d0
      Te( 6) = y; Gradte(2, 6) = 1.d0
      Te( 8) = y; Gradte(2, 8) = 1.d0
c
c  ...edges parallel to the z-axis
      Te( 9) = z; Gradte(3, 9) = 1.d0
      Te(10) = z; Gradte(3,10) = 1.d0
      Te(11) = z; Gradte(3,11) = 1.d0
      Te(12) = z; Gradte(3,12) = 1.d0
c
      if (iprint.eq.1) then
        write(*,7001) Xi(1:3)
 7001   format('blend_project_edges_bric: Xi = ',3f8.3)
        write(*,7002) Te(1:12)
 7002   format('Te = ',12f8.3)
        call pause
      endif
c
c
      end subroutine blend_project_edges_bric


c----------------------------------------------------------------------
c
c   routine name       - blend_project_faces_bric
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 13
c
c   purpose            - routine returns values of brick face blending
c                        functions along with brick to face projections
c                        with their gradients
c
c   arguments :
c     in:
c        Xi            - master element coordinates of a point
c
c     out:
c        Psif          - face blending function values
c        Gradpsif      - gradient of psif
c        Tf            - element to face projections
c        Gradtf        - gradient of Tf
c
c----------------------------------------------------------------------

      subroutine blend_project_faces_bric(Xi, Psif,Gradpsif,Tf,Gradtf)
c
#include "syscom.blk"
c
      dimension Xi(3),
     .          Psif(6),Gradpsif(3,6),Tf(2,6),Gradtf(2,3,6)
      double precision z,z1
c
c----------------------------------------------------------------------
c
      Gradpsif = 0.d0; Gradtf = 0.d0
      x = Xi(1); x1 = 1.d0-x
      y = Xi(2); y1 = 1.d0-y
      z = Xi(3); z1 = 1.d0-z
c
c  ...face blending functions:
      Psif(1) = z1; Gradpsif(3,1) = -1.d0
      Psif(2) = z ; Gradpsif(3,2) =  1.d0
      Psif(3) = y1; Gradpsif(2,3) = -1.d0
      Psif(5) = y ; Gradpsif(2,5) =  1.d0
      Psif(4) = x ; Gradpsif(1,4) =  1.d0
      Psif(6) = x1; Gradpsif(1,6) = -1.d0
c
c  ...face projections:
      Tf(1,1)=x; Tf(2,1)=y; Gradtf(1,1,1) = 1.d0; Gradtf(2,2,1) = 1.d0
      Tf(1,2)=x; Tf(2,2)=y; Gradtf(1,1,2) = 1.d0; Gradtf(2,2,2) = 1.d0
      Tf(1,3)=x; Tf(2,3)=z; Gradtf(1,1,3) = 1.d0; Gradtf(2,3,3) = 1.d0
      Tf(1,5)=x; Tf(2,5)=z; Gradtf(1,1,5) = 1.d0; Gradtf(2,3,5) = 1.d0
      Tf(1,4)=y; Tf(2,4)=z; Gradtf(1,2,4) = 1.d0; Gradtf(2,3,4) = 1.d0
      Tf(1,6)=y; Tf(2,6)=z; Gradtf(1,2,6) = 1.d0; Gradtf(2,3,6) = 1.d0
c
c
      end subroutine blend_project_faces_bric
c
c----------------------------------------------------------------------
c
c   routine name       - blend_pyram
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 14
c
c   purpose            - routine computes a number of auxiliary
c                        quantities common to all pyramid shape
c                        functions routines
c
c   arguments :
c     in:
c        Xip           - master element coordinates of a point
c
c     out:
c        X,GradXmY,GradY,Z,GradZ - infinite element coordinates
c                        with gradients
c        Psi,Gradpsi   - common factors for edge vertex shape
c                        functions with gradients
c        Alpha,GradAlpha - modified edge vertex shape functions
c                        with gradients
c        Vshap,Gradv   - vertex shape functions an dtheir gradients
c
c----------------------------------------------------------------------

      subroutine blend_pyram(Xip,X,GradX,Y,GradY,Z,GradZ,
     .                       Psi,GradPsi,Alpha,GradAlpha,
     .                       Vshap,Gradv)
c
#include "syscom.blk"
c
      dimension Xip(3),
     .          GradX(3),GradY(3),GradZ(3),
     .          Psi(4),Gradpsi(3,4),
     .          Alpha(2,4),Gradalpha(3,2,4),
     .          Vshap(5),Gradv(3,5)
c
      double precision Z,zeta
      real*8,parameter :: eps = 1.d-12
c
      xi = Xip(1); eta = Xip(2); zeta = Xip(3)
      Gradv = 0.d0
c
c  ...avoid singularity at zeta = 1
      if (abs(zeta - 1.d0) < eps)  zeta=1.d0-eps
c
c  ...infinite element coordinates with gradients
      X = xi/(1.d0-zeta)
      GradX = (/1.d0/(1.d0-zeta), 0.d0, xi/(1.d0-zeta)**2/)
      Y = eta/(1.d0-zeta)
      GradY = (/0.d0, 1.d0/(1.d0-zeta), eta/(1.d0-zeta)**2/)
      Z = zeta/(1.d0-zeta)
      GradZ = (/0.d0, 0.d0, 1.d0/(1.d0-zeta)**2/)
c
c  ...common factors for edge vertex shape functions
      Psi(1) = 1.d0-Y; GradPsi(1:3,1) = -gradY(1:3)
      Psi(2) =      X; GradPsi(1:3,2) =  gradX(1:3)
      Psi(3) =      Y; GradPsi(1:3,3) =  gradY(1:3)
      Psi(4) = 1.d0-X; GradPsi(1:3,4) = -gradX(1:3)
c
c  ...modified edge vertex shape functions
      do i=1,4
        select case(i)
        case(1,3)
          Alpha(1,i) = 1.d0-xi-zeta
          GradAlpha(1:3,1,i) = (/-1.d0, 0.d0,-1.d0/)
          Alpha(2,i) =      xi
          GradAlpha(1:3,2,i) = (/ 1.d0, 0.d0, 0.d0/)
        case(2,4)
          Alpha(1,i) = 1.d0-eta-zeta
          GradAlpha(1:3,1,i) = (/ 0.d0,-1.d0,-1.d0/)
          Alpha(2,i) =      eta
          GradAlpha(1:3,2,i) = (/ 0.d0, 1.d0, 0.d0/)
        end select
      enddo
c
c  ...vertex shape functions
      Vshap(1) =     (1.d0-X)*(1.d0-eta-zeta)
      Gradv(1:3,1) =   -GradX*(1.d0-eta-zeta)
     .             + (1.d0-X)*(/0.d0,-1.d0,-1.d0/)
      Vshap(2) =           X *(1.d0-eta-zeta)
      Gradv(1:3,2) =    GradX*(1.d0-eta-zeta)
     .             +        X*(/0.d0,-1.d0,-1.d0/)
      Vshap(3) =           X * eta
      Gradv(1:3,3) =    GradX* eta
     .             +        X*(/0.d0, 1.d0, 0.d0/)
      Vshap(4) =     (1.d0-X)* eta
      Gradv(1:3,4) =   -GradX* eta
     .             + (1.d0-X)*(/0.d0, 1.d0, 0.d0/)
      Vshap(5) = zeta
      Gradv(1:3,5) = (/0.d0, 0.d0, 1.d0/)
c
c
      end subroutine blend_pyram
