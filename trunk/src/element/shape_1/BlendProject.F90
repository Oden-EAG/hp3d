! Routines:
!  - BlendSegV
!  - BlendQuadV
!  - BlendProjectQuadE
!  - BlendTriV
!  - ProjectTriE
!  - BlendHexaV
!  - BlendProjectHexaE
!  - BlendProjectHexaF
!  - BlendTetV
!  - ProjectTetE
!  - ProjectTetF
!  - BlendPrisV
!  - BlendProjectPrisME
!  - BlendProjectPrisQE
!  - BlendProjectPrisTF
!  - ProjectPrisQF
!  - BlendPyraV
!  - BlendProjectPyraME
!  - ProjectPyraTE
!  - ProjectPyraQF
!  - BlendProjectPyraTF
!  - ProjectPyraLamTF
!----------------------------------------------------------------------
!  b=blending p=projecting V=vertex E=edge F=face M=mixed T=tri Q=quad
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  Data consistent with module element_data:
!     numbering of topological entities - vertices, edges, faces
!     local orientations of edges and faces
!----------------------------------------------------------------------
   subroutine BlendSegV(Mu,DMu, MubV,DMubV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Mu(0:1),DMu(0:1)
      double precision, intent(out) :: MubV(1:2),DMubV(1:2)
!
!  ...Info from module element_data - coordinates,connectivities:
!         V=((0),(1))=>(v0,v1)
!
      N=1
!
!  ...2 vertices, each with one blending function
!
!  ...v=1 --> v0=(0)
      v=1
      MubV(v) = Mu(0); DMubV(v) = DMu(0)
!  ...v=2 --> v1=(1)
      v=2
      MubV(v) = Mu(1); DMubV(v) = DMu(1)
!
   end subroutine BlendSegV
!----------------------------------------------------------------------
   subroutine BlendQuadV(Mu,DMu, MubV,DMubV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Mu(0:1,1:2),DMu(1:2,0:1,1:2)
      double precision, intent(out) :: MubV(1:2,1:4),DMubV(1:2,1:2,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!         V=((0,0),(1,0),(1,1),(0,1))=>(v1,v2,v3,v4)
!
      N=2
!
!  ...4 vertices, each with two blending functions
!
!  ...v=1 --> v1=(0,0)
      v=1
      MubV(1,v) = Mu(0,1); DMubV(1:N,1,v) = DMu(1:N,0,1)
      MubV(2,v) = Mu(0,2); DMubV(1:N,2,v) = DMu(1:N,0,2)
!  ...v=2 --> v2=(1,0)
      v=2
      MubV(1,v) = Mu(1,1); DMubV(1:N,1,v) = DMu(1:N,1,1)
      MubV(2,v) = Mu(0,2); DMubV(1:N,2,v) = DMu(1:N,0,2)
!  ...v=3 --> v3=(1,1)
      v=3
      MubV(1,v) = Mu(1,1); DMubV(1:N,1,v) = DMu(1:N,1,1)
      MubV(2,v) = Mu(1,2); DMubV(1:N,2,v) = DMu(1:N,1,2)
!  ...v=4 --> v4=(0,1)
      v=4
      MubV(1,v) = Mu(0,1); DMubV(1:N,1,v) = DMu(1:N,0,1)
      MubV(2,v) = Mu(1,2); DMubV(1:N,2,v) = DMu(1:N,1,2)
!
   end subroutine BlendQuadV
!----------------------------------------------------------------------
   subroutine BlendProjectQuadE(Mu,DMu, MubE,DMubE,MupE,DMupE,IdecE)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecE
      double precision, intent(in)  :: Mu(0:1,1:2),DMu(1:2,0:1,1:2)
      double precision, intent(out) :: MubE(1:4),DMubE(1:2,1:4)
      double precision, intent(out) :: MupE(0:1,1:4),DMupE(1:2,0:1,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!         V=((0,0),(1,0),(1,1),(0,1))=>(v1,v2,v3,v4)
!          E=>((v1->v2),(v2->v3),(v4->v3),(v1->v4))
!
      N=2
!
!  ...4 edges, each with a blending function and a locally oriented
!     pair representing a projection
!
!  ...e=1 --> edge12 with local orientation v1->v2
      e=1
      MubE(e) = Mu(0,2); DMubE(1:N,e) = DMu(1:N,0,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,1); MupE(1,e) = Mu(1,1)
      DMupE(1:N,0,e) = DMu(1:N,0,1); DMupE(1:N,1,e) = DMu(1:N,1,1)
!  ...e=2 --> edge23 with local orientation v2->v3
      e=2
      MubE(e) = Mu(1,1); DMubE(1:N,e) = DMu(1:N,1,1)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,2); MupE(1,e) = Mu(1,2);
      DMupE(1:N,0,e) = DMu(1:N,0,2); DMupE(1:N,1,e) = DMu(1:N,1,2);
!  ...e=3 --> edge34 with local orientation v4->v3
      e=3
      MubE(e) = Mu(1,2); DMubE(1:N,e) = DMu(1:N,1,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,1); MupE(1,e) = Mu(1,1)
      DMupE(1:N,0,e) = DMu(1:N,0,1); DMupE(1:N,1,e) = DMu(1:N,1,1)
!  ...e=4 --> edge41 with local orientation v1->v4
      e=4
      MubE(e) = Mu(0,1); DMubE(1:N,e) = DMu(1:N,0,1)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,2); MupE(1,e) = Mu(1,2)
      DMupE(1:N,0,e) = DMu(1:N,0,2); DMupE(1:N,1,e) = DMu(1:N,1,2)
!
!  ...projected coordinates are Mu, so IdecE=true for all edges
      IdecE = .TRUE.
!
   end subroutine BlendProjectQuadE
!----------------------------------------------------------------------
   subroutine BlendTriV(Nu,DNu, NubV,DNubV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Nu(0:2),DNu(1:2,0:2)
      double precision, intent(out) :: NubV(1:3),DNubV(1:2,1:3)
!
!  ...Info from module element_data - coordinates,connectivities:
!                V=((0,0),(1,0),(0,1))=>(v0,v1,v2)
!
      N=2
!
!  ...3 vertices, each with one blending function
!
!  ...v=1 --> v0=(0,0)
      v=1
      NubV(v) = Nu(0); DNubV(1:N,v) = DNu(1:N,0)
!  ...v=2 --> v1=(1,0)
      v=2
      NubV(v) = Nu(1); DNubV(1:N,v) = DNu(1:N,1)
!  ...v=3 --> v2=(0,1)
      v=3
      NubV(v) = Nu(2); DNubV(1:N,v) = DNu(1:N,2)
!
   end subroutine BlendTriV
!----------------------------------------------------------------------
   subroutine ProjectTriE(Nu,DNu, NupE,DNupE,IdecE)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecE
      double precision, intent(in)  :: Nu(0:2),DNu(1:2,0:2)
      double precision, intent(out) :: NupE(0:1,1:3),DNupE(1:2,0:1,1:3)
!
!  ...Info from module element_data - coordinates,connectivities:
!                V=((0,0),(1,0),(0,1))=>(v0,v1,v2)
!                 E=>((v0->v1),(v1->v2),(v0->v2))
!
      N=2
!
!  ...3 edges, each with a locally oriented pair representing
!     a projection
!
!  ...e=1 --> edge01 with local orientation v0->v1
      e=1
      NupE(0,e) = Nu(0); NupE(1,e) = Nu(1)
      DNupE(1:N,0,e) = DNu(1:N,0); DNupE(1:N,1,e) = DNu(1:N,1)
!  ...e=2 --> edge12 with local orientation v1->v2
      e=2
      NupE(0,e) = Nu(1); NupE(1,e) = Nu(2)
      DNupE(1:N,0,e) = DNu(1:N,1); DNupE(1:N,1,e) = DNu(1:N,2)
!  ...e=3 --> edge20 with local orientation v0->v2
      e=3
      NupE(0,e) = Nu(0); NupE(1,e) = Nu(2)
      DNupE(1:N,0,e) = DNu(1:N,0); DNupE(1:N,1,e) = DNu(1:N,2)
!
!  ...projected coordinates are Nu, so IdecE=false for all edges
      IdecE = .FALSE.
!
   end subroutine ProjectTriE
!----------------------------------------------------------------------
   subroutine BlendHexaV(Mu,DMu, MubV,DMubV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision, intent(out) :: MubV(1:3,1:8),DMubV(1:3,1:3,1:8)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
!           (0,0,1),(1,0,1),(1,1,1),(0,1,1))=>(v1,v2,v3,v4,v5,v6,v7,v8)
!
      N=3
!
!  ...8 vertices, each with three blending functions
!
!  ...v=1 --> v1=(0,0,0)
      v=1
      MubV(1,v) = Mu(0,1); DMubV(1:N,1,v) = DMu(1:N,0,1)
      MubV(2,v) = Mu(0,2); DMubV(1:N,2,v) = DMu(1:N,0,2)
      MubV(3,v) = Mu(0,3); DMubV(1:N,3,v) = DMu(1:N,0,3)
!  ...v=2 --> v2=(1,0,0)
      v=2
      MubV(1,v) = Mu(1,1); DMubV(1:N,1,v) = DMu(1:N,1,1)
      MubV(2,v) = Mu(0,2); DMubV(1:N,2,v) = DMu(1:N,0,2)
      MubV(3,v) = Mu(0,3); DMubV(1:N,3,v) = DMu(1:N,0,3)
!  ...v=3 --> v3=(1,1,0)
      v=3
      MubV(1,v) = Mu(1,1); DMubV(1:N,1,v) = DMu(1:N,1,1)
      MubV(2,v) = Mu(1,2); DMubV(1:N,2,v) = DMu(1:N,1,2)
      MubV(3,v) = Mu(0,3); DMubV(1:N,3,v) = DMu(1:N,0,3)
!  ...v=4 --> v4=(0,1,0)
      v=4
      MubV(1,v) = Mu(0,1); DMubV(1:N,1,v) = DMu(1:N,0,1)
      MubV(2,v) = Mu(1,2); DMubV(1:N,2,v) = DMu(1:N,1,2)
      MubV(3,v) = Mu(0,3); DMubV(1:N,3,v) = DMu(1:N,0,3)
!  ...v=5 --> v5=(0,0,1)
      v=5
      MubV(1,v) = Mu(0,1); DMubV(1:N,1,v) = DMu(1:N,0,1)
      MubV(2,v) = Mu(0,2); DMubV(1:N,2,v) = DMu(1:N,0,2)
      MubV(3,v) = Mu(1,3); DMubV(1:N,3,v) = DMu(1:N,1,3)
!  ...v=6 --> v6=(1,0,1)
      v=6
      MubV(1,v) = Mu(1,1); DMubV(1:N,1,v) = DMu(1:N,1,1)
      MubV(2,v) = Mu(0,2); DMubV(1:N,2,v) = DMu(1:N,0,2)
      MubV(3,v) = Mu(1,3); DMubV(1:N,3,v) = DMu(1:N,1,3)
!  ...v=7 --> v7=(1,1,1)
      v=7
      MubV(1,v) = Mu(1,1); DMubV(1:N,1,v) = DMu(1:N,1,1)
      MubV(2,v) = Mu(1,2); DMubV(1:N,2,v) = DMu(1:N,1,2)
      MubV(3,v) = Mu(1,3); DMubV(1:N,3,v) = DMu(1:N,1,3)
!  ...v=8 --> v8=(0,1,1)
      v=8
      MubV(1,v) = Mu(0,1); DMubV(1:N,1,v) = DMu(1:N,0,1)
      MubV(2,v) = Mu(1,2); DMubV(1:N,2,v) = DMu(1:N,1,2)
      MubV(3,v) = Mu(1,3); DMubV(1:N,3,v) = DMu(1:N,1,3)
!
   end subroutine BlendHexaV
!----------------------------------------------------------------------
   subroutine BlendProjectHexaE(Mu,DMu, MubE,DMubE,MupE,DMupE,IdecE)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecE
      double precision, intent(in)  :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision, intent(out) :: MubE(1:2,1:12), &
                DMubE(1:3,1:2,1:12),MupE(0:1,1:12),DMupE(1:3,0:1,1:12)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
!           (0,0,1),(1,0,1),(1,1,1),(0,1,1))=>(v1,v2,v3,v4,v5,v6,v7,v8)
!      E=>((v1->v2),(v2->v3),(v4->v3),(v1->v4),(v5->v6),(v6->v7),
!            (v8->v7),(v5->v8),(v1->v5),(v2->v6),(v3->v7),(v4->v8))
!
      N=3
!
!  ...12 edges, each with two blending functions
!     and a locally oriented pair representing a projection
!
!  ...e=1 --> edge12 with local orientation v1->v2
      e=1
      MubE(1,e) = Mu(0,2); DMubE(1:N,1,e) = DMu(1:N,0,2)
      MubE(2,e) = Mu(0,3); DMubE(1:N,2,e) = DMu(1:N,0,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,1); MupE(1,e) = Mu(1,1)
      DMupE(1:N,0,e) = DMu(1:N,0,1); DMupE(1:N,1,e) = DMu(1:N,1,1)
!  ...e=2 --> edge23 with local orientation v2->v3
      e=2
      MubE(1,e) = Mu(1,1); DMubE(1:N,1,e) = DMu(1:N,1,1)
      MubE(2,e) = Mu(0,3); DMubE(1:N,2,e) = DMu(1:N,0,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,2); MupE(1,e) = Mu(1,2);
      DMupE(1:N,0,e) = DMu(1:N,0,2); DMupE(1:N,1,e) = DMu(1:N,1,2);
!  ...e=3 --> edge34 with local orientation v4->v3
      e=3
      MubE(1,e) = Mu(1,2); DMubE(1:N,1,e) = DMu(1:N,1,2)
      MubE(2,e) = Mu(0,3); DMubE(1:N,2,e) = DMu(1:N,0,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,1); MupE(1,e) = Mu(1,1)
      DMupE(1:N,0,e) = DMu(1:N,0,1); DMupE(1:N,1,e) = DMu(1:N,1,1)
!  ...e=4 --> edge41 with local orientation v1->v4
      e=4
      MubE(1,e) = Mu(0,1); DMubE(1:N,1,e) = DMu(1:N,0,1)
      MubE(2,e) = Mu(0,3); DMubE(1:N,2,e) = DMu(1:N,0,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,2); MupE(1,e) = Mu(1,2)
      DMupE(1:N,0,e) = DMu(1:N,0,2); DMupE(1:N,1,e) = DMu(1:N,1,2)
!  ...e=5 --> edge56 with local orientation v5->v6
      e=5
      MubE(1,e) = Mu(0,2); DMubE(1:N,1,e) = DMu(1:N,0,2)
      MubE(2,e) = Mu(1,3); DMubE(1:N,2,e) = DMu(1:N,1,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,1); MupE(1,e) = Mu(1,1)
      DMupE(1:N,0,e) = DMu(1:N,0,1); DMupE(1:N,1,e) = DMu(1:N,1,1)
!  ...e=6 --> edge67 with local orientation v6->v7
      e=6
      MubE(1,e) = Mu(1,1); DMubE(1:N,1,e) = DMu(1:N,1,1)
      MubE(2,e) = Mu(1,3); DMubE(1:N,2,e) = DMu(1:N,1,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,2); MupE(1,e) = Mu(1,2);
      DMupE(1:N,0,e) = DMu(1:N,0,2); DMupE(1:N,1,e) = DMu(1:N,1,2);
!  ...e=7 --> edge78 with local orientation v8->v7
      e=7
      MubE(1,e) = Mu(1,2); DMubE(1:N,1,e) = DMu(1:N,1,2)
      MubE(2,e) = Mu(1,3); DMubE(1:N,2,e) = DMu(1:N,1,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,1); MupE(1,e) = Mu(1,1)
      DMupE(1:N,0,e) = DMu(1:N,0,1); DMupE(1:N,1,e) = DMu(1:N,1,1)
!  ...e=8 --> edge85 with local orientation v5->v8
      e=8
      MubE(1,e) = Mu(0,1); DMubE(1:N,1,e) = DMu(1:N,0,1)
      MubE(2,e) = Mu(1,3); DMubE(1:N,2,e) = DMu(1:N,1,3)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,2); MupE(1,e) = Mu(1,2)
      DMupE(1:N,0,e) = DMu(1:N,0,2); DMupE(1:N,1,e) = DMu(1:N,1,2)
!  ...e=9 --> edge15 with local orientation v1->v5
      e=9
      MubE(1,e) = Mu(0,1); DMubE(1:N,1,e) = DMu(1:N,0,1)
      MubE(2,e) = Mu(0,2); DMubE(1:N,2,e) = DMu(1:N,0,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,3); MupE(1,e) = Mu(1,3)
      DMupE(1:N,0,e) = DMu(1:N,0,3); DMupE(1:N,1,e) = DMu(1:N,1,3)
!  ...e=10 --> edge26 with local orientation v2->v6
      e=10
      MubE(1,e) = Mu(1,1); DMubE(1:N,1,e) = DMu(1:N,1,1)
      MubE(2,e) = Mu(0,2); DMubE(1:N,2,e) = DMu(1:N,0,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,3); MupE(1,e) = Mu(1,3)
      DMupE(1:N,0,e) = DMu(1:N,0,3); DMupE(1:N,1,e) = DMu(1:N,1,3)
!  ...e=11 --> edge37 with local orientation v3->v7
      e=11
      MubE(1,e) = Mu(1,1); DMubE(1:N,1,e) = DMu(1:N,1,1)
      MubE(2,e) = Mu(1,2); DMubE(1:N,2,e) = DMu(1:N,1,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,3); MupE(1,e) = Mu(1,3)
      DMupE(1:N,0,e) = DMu(1:N,0,3); DMupE(1:N,1,e) = DMu(1:N,1,3)
!  ...e=12 --> edge48 with local orientation v4->v8
      e=12
      MubE(1,e) = Mu(0,1); DMubE(1:N,1,e) = DMu(1:N,0,1)
      MubE(2,e) = Mu(1,2); DMubE(1:N,2,e) = DMu(1:N,1,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0,3); MupE(1,e) = Mu(1,3)
      DMupE(1:N,0,e) = DMu(1:N,0,3); DMupE(1:N,1,e) = DMu(1:N,1,3)
!
!  ...projected coordinates are Mu, so IdecE=true for all edges
      IdecE = .TRUE.
!
   end subroutine BlendProjectHexaE
!----------------------------------------------------------------------
   subroutine BlendProjectHexaF(Mu,DMu, MubF,DMubF,MupF,DMupF,IdecF)
!
      implicit none
      integer :: N,f
      logical, intent(out) :: IdecF(1:2)
      double precision, intent(in)  :: Mu(0:1,1:3),DMu(1:3,0:1,1:3)
      double precision, intent(out) :: MubF(1:6),DMubF(1:3,1:6)
      double precision, intent(out) :: MupF(0:1,1:2,1:6), &
                                                 DMupF(1:3,0:1,1:2,1:6)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
!           (0,0,1),(1,0,1),(1,1,1),(0,1,1))=>(v1,v2,v3,v4,v5,v6,v7,v8)
!      F=>((v1->v2->v3->v4),(v5->v6->v7->v8),(v1->v2->v6->v5)
!            (v2->v3->v7->v6),(v4->v3->v7->v8),(v1->v4->v8->v5))
!
      N=3
!
!  ...6 faces, each with one blending function
!     and a locally oriented quadruple representing a projection
!
!  ...f=1 --> face1234 with local orientation v1->v2->v3->v4
      f=1
      MubF(f) = Mu(0,3); DMubF(1:N,f) = DMu(1:N,0,3)
!     ...locally oriented quadruple representing projection
      MupF(0,1,f) = Mu(0,1); MupF(1,1,f) = Mu(1,1)
      MupF(0,2,f) = Mu(0,2); MupF(1,2,f) = Mu(1,2)
      DMupF(1:N,0,1,f) = DMu(1:N,0,1); DMupF(1:N,1,1,f) = DMu(1:N,1,1)
      DMupF(1:N,0,2,f) = DMu(1:N,0,2); DMupF(1:N,1,2,f) = DMu(1:N,1,2)
!  ...f=2 --> face5678 with local orientation v5->v6->v7->v8
      f=2
      MubF(f) = Mu(1,3); DMubF(1:N,f) = DMu(1:N,1,3)
!     ...locally oriented quadruple representing projection
      MupF(0,1,f) = Mu(0,1); MupF(1,1,f) = Mu(1,1)
      MupF(0,2,f) = Mu(0,2); MupF(1,2,f) = Mu(1,2)
      DMupF(1:N,0,1,f) = DMu(1:N,0,1); DMupF(1:N,1,1,f) = DMu(1:N,1,1)
      DMupF(1:N,0,2,f) = DMu(1:N,0,2); DMupF(1:N,1,2,f) = DMu(1:N,1,2)
!  ...f=3 --> face1265 with local orientation v1->v2->v6->v5
      f=3
      MubF(f) = Mu(0,2); DMubF(1:N,f) = DMu(1:N,0,2)
!     ...locally oriented quadruple representing projection
      MupF(0,1,f) = Mu(0,1); MupF(1,1,f) = Mu(1,1)
      MupF(0,2,f) = Mu(0,3); MupF(1,2,f) = Mu(1,3)
      DMupF(1:N,0,1,f) = DMu(1:N,0,1); DMupF(1:N,1,1,f) = DMu(1:N,1,1)
      DMupF(1:N,0,2,f) = DMu(1:N,0,3); DMupF(1:N,1,2,f) = DMu(1:N,1,3)
!  ...f=4 --> face2376 with local orientation v2->v3->v7->v6
      f=4
      MubF(f) = Mu(1,1); DMubF(1:N,f) = DMu(1:N,1,1)
!     ...locally oriented quadruple representing projection
      MupF(0,1,f) = Mu(0,2); MupF(1,1,f) = Mu(1,2)
      MupF(0,2,f) = Mu(0,3); MupF(1,2,f) = Mu(1,3)
      DMupF(1:N,0,1,f) = DMu(1:N,0,2); DMupF(1:N,1,1,f) = DMu(1:N,1,2)
      DMupF(1:N,0,2,f) = DMu(1:N,0,3); DMupF(1:N,1,2,f) = DMu(1:N,1,3)
!  ...f=5 --> face4378 with local orientation v4->v3->v7->v8
      f=5
      MubF(f) = Mu(1,2); DMubF(1:N,f) = DMu(1:N,1,2)
!     ...locally oriented quadruple representing projection
      MupF(0,1,f) = Mu(0,1); MupF(1,1,f) = Mu(1,1)
      MupF(0,2,f) = Mu(0,3); MupF(1,2,f) = Mu(1,3)
      DMupF(1:N,0,1,f) = DMu(1:N,0,1); DMupF(1:N,1,1,f) = DMu(1:N,1,1)
      DMupF(1:N,0,2,f) = DMu(1:N,0,3); DMupF(1:N,1,2,f) = DMu(1:N,1,3)
!  ...f=6 --> face1485 with local orientation v1->v4->v8->v5
      f=6
      MubF(f) = Mu(0,1); DMubF(1:N,f) = DMu(1:N,0,1)
!     ...locally oriented quadruple representing projection
      MupF(0,1,f) = Mu(0,2); MupF(1,1,f) = Mu(1,2)
      MupF(0,2,f) = Mu(0,3); MupF(1,2,f) = Mu(1,3)
      DMupF(1:N,0,1,f) = DMu(1:N,0,2); DMupF(1:N,1,1,f) = DMu(1:N,1,2)
      DMupF(1:N,0,2,f) = DMu(1:N,0,3); DMupF(1:N,1,2,f) = DMu(1:N,1,3)
!
!  ...projected coordinates are Mu and Mu, so IdecF=(true,true) for
!     all faces
      IdecF(1) = .TRUE.; IdecF(2) = .TRUE.
!
   end subroutine BlendProjectHexaF
!----------------------------------------------------------------------
   subroutine BlendTetV(Lam,DLam, LambV,DLambV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Lam(0:3),DLam(1:3,0:3)
      double precision, intent(out) :: LambV(1:4),DLambV(1:3,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!           V=((0,0,0),(1,0,0),(0,1,0),(0,0,1))=>(v0,v1,v2,v3)
!
      N=3
!
!  ...4 vertices, each with one blending function
!
!  ...v=1 --> v0=(0,0,0)
      v=1
      LambV(v) = Lam(0); DLambV(1:N,v) = DLam(1:N,0)
!  ...v=2 --> v1=(1,0,0)
      v=2
      LambV(v) = Lam(1); DLambV(1:N,v) = DLam(1:N,1)
!  ...v=3 --> v2=(0,1,0)
      v=3
      LambV(v) = Lam(2); DLambV(1:N,v) = DLam(1:N,2)
!  ...v=4 --> v3=(0,0,1)
      v=4
      LambV(v) = Lam(3); DLambV(1:N,v) = DLam(1:N,3)
!
   end subroutine BlendTetV
!----------------------------------------------------------------------
   subroutine ProjectTetE(Lam,DLam, LampE,DLampE,IdecE)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecE
      double precision, intent(in)  :: Lam(0:3),DLam(1:3,0:3)
      double precision, intent(out) :: LampE(0:1,1:6), &
                                                   DLampE(1:3,0:1,1:6)
!
!  ...Info from module element_data - coordinates,connectivities:
!           V=((0,0,0),(1,0,0),(0,1,0),(0,0,1))=>(v0,v1,v2,v3)
!         E=>((v0->v1),(v1->v2),(v0->v2),(v0->v3),(v1->v3),(v2->v3))
!
      N=3
!
!  ...6 edges, each with a locally oriented pair representing
!     a projection
!
!  ...e=1 --> edge01 with local orientation v0->v1
      e=1
      LampE(0,e) = Lam(0); LampE(1,e) = Lam(1)
      DLampE(1:N,0,e) = DLam(1:N,0); DLampE(1:N,1,e) = DLam(1:N,1)
!  ...e=2 --> edge12 with local orientation v1->v2
      e=2
      LampE(0,e) = Lam(1); LampE(1,e) = Lam(2)
      DLampE(1:N,0,e) = DLam(1:N,1); DLampE(1:N,1,e) = DLam(1:N,2)
!  ...e=3 --> edge20 with local orientation v0->v2
      e=3
      LampE(0,e) = Lam(0); LampE(1,e) = Lam(2)
      DLampE(1:N,0,e) = DLam(1:N,0); DLampE(1:N,1,e) = DLam(1:N,2)
!  ...e=4 --> edge03 with local orientation v0->v3
      e=4
      LampE(0,e) = Lam(0); LampE(1,e) = Lam(3)
      DLampE(1:N,0,e) = DLam(1:N,0); DLampE(1:N,1,e) = DLam(1:N,3)
!  ...e=5 --> edge13 with local orientation v1->v3
      e=5
      LampE(0,e) = Lam(1); LampE(1,e) = Lam(3)
      DLampE(1:N,0,e) = DLam(1:N,1); DLampE(1:N,1,e) = DLam(1:N,3)
!  ...e=6 --> edge23 with local orientation v2->v3
      e=6
      LampE(0,e) = Lam(2); LampE(1,e) = Lam(3)
      DLampE(1:N,0,e) = DLam(1:N,2); DLampE(1:N,1,e) = DLam(1:N,3)
!
!  ...projected coordinates are Lam, so IdecE=false for all edges
      IdecE = .FALSE.
!
   end subroutine ProjectTetE
!----------------------------------------------------------------------
   subroutine ProjectTetF(Lam,DLam, LampF,DLampF,IdecF)
!
      implicit none
      integer :: N,f
      logical, intent(out) :: IdecF
      double precision, intent(in)  :: Lam(0:3),DLam(1:3,0:3)
      double precision, intent(out) :: LampF(0:2,1:4), &
                                                   DLampF(1:3,0:2,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!           V=((0,0,0),(1,0,0),(0,1,0),(0,0,1))=>(v0,v1,v2,v3)
!        F=>((v0->v1->v2),(v0->v1->v3),(v1->v2->v3),(v0->v2->v3))
!
      N=3
!
!  ...4 faces, each with a locally oriented triplet representing
!     a projection
!
!  ...f=1 --> face012 with local orientation v0->v1->v2
      f=1
      LampF(0,f) = Lam(0); LampF(1,f) = Lam(1); LampF(2,f) = Lam(2)
      DLampF(1:N,0,f) = DLam(1:N,0); DLampF(1:N,1,f) = DLam(1:N,1);
                                         DLampF(1:N,2,f) = DLam(1:N,2)
!  ...f=2 --> face013 with local orientation v0->v1->v3
      f=2
      LampF(0,f) = Lam(0); LampF(1,f) = Lam(1); LampF(2,f) = Lam(3)
      DLampF(1:N,0,f) = DLam(1:N,0); DLampF(1:N,1,f) = DLam(1:N,1);
                                         DLampF(1:N,2,f) = DLam(1:N,3)
!  ...f=3 --> face123 with local orientation v1->v2->v3
      f=3
      LampF(0,f) = Lam(1); LampF(1,f) = Lam(2); LampF(2,f) = Lam(3)
      DLampF(1:N,0,f) = DLam(1:N,1); DLampF(1:N,1,f) = DLam(1:N,2);
                                         DLampF(1:N,2,f) = DLam(1:N,3)
!  ...f=4 --> face023 with local orientation v0->v2->v3
      f=4
      LampF(0,f) = Lam(0); LampF(1,f) = Lam(2); LampF(2,f) = Lam(3)
      DLampF(1:N,0,f) = DLam(1:N,0); DLampF(1:N,1,f) = DLam(1:N,2);
                                         DLampF(1:N,2,f) = DLam(1:N,3)
!
!  ...projected coordinates are Lam, so IdecF=false for all faces
      IdecF = .FALSE.
!
   end subroutine ProjectTetF
!----------------------------------------------------------------------
   subroutine BlendPrisV(Mu,DMu,Nu,DNu, MubV,DMubV,NubV,DNubV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(in)  :: Nu(0:2),DNu(1:3,0:2)
      double precision, intent(out) :: MubV(1:6),DMubV(1:3,1:6)
      double precision, intent(out) :: NubV(1:6),DNubV(1:3,1:6)
!
!  ...Info from module element_data - coordinates,connectivities:
!        V=((0,0,0),(1,0,0),(0,1,0),
!                    (0,0,1),(1,0,1),(0,1,1))=>(v0,v1,v2,v3,v4,v5)
!
      N=3
!
!  ...6 vertices, each with two blending functions (one mu, one nu)
!
!  ...v=1 --> v0=(0,0,0)
      v=1
      MubV(v) = Mu(0); DMubV(1:N,v) = DMu(1:N,0)
      NubV(v) = Nu(0); DNubV(1:N,v) = DNu(1:N,0)
!  ...v=2 --> v1=(1,0,0)
      v=2
      MubV(v) = Mu(0); DMubV(1:N,v) = DMu(1:N,0)
      NubV(v) = Nu(1); DNubV(1:N,v) = DNu(1:N,1)
!  ...v=3 --> v2=(0,1,0)
      v=3
      MubV(v) = Mu(0); DMubV(1:N,v) = DMu(1:N,0)
      NubV(v) = Nu(2); DNubV(1:N,v) = DNu(1:N,2)
!  ...v=4 --> v3=(0,0,1)
      v=4
      MubV(v) = Mu(1); DMubV(1:N,v) = DMu(1:N,1)
      NubV(v) = Nu(0); DNubV(1:N,v) = DNu(1:N,0)
!  ...v=5 --> v4=(1,0,1)
      v=5
      MubV(v) = Mu(1); DMubV(1:N,v) = DMu(1:N,1)
      NubV(v) = Nu(1); DNubV(1:N,v) = DNu(1:N,1)
!  ...v=6 --> v5=(0,1,1)
      v=6
      MubV(v) = Mu(1); DMubV(1:N,v) = DMu(1:N,1)
      NubV(v) = Nu(2); DNubV(1:N,v) = DNu(1:N,2)
!
   end subroutine BlendPrisV
!----------------------------------------------------------------------
   subroutine BlendProjectPrisME(Mu,DMu,Nu,DNu, &
                                          MubE,DMubE,NupE,DNupE,IdecME)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecME
      double precision, intent(in)  :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(in)  :: Nu(0:2),DNu(1:3,0:2)
      double precision, intent(out) :: MubE(1:6),DMubE(1:3,1:6)
      double precision, intent(out) :: NupE(0:1,1:6),DNupE(1:3,0:1,1:6)
!
!  ...Info from module element_data - coordinates,connectivities:
!        V=((0,0,0),(1,0,0),(0,1,0),
!                    (0,0,1),(1,0,1),(0,1,1))=>(v0,v1,v2,v3,v4,v5)
!      E_part1=>((v0->v1),(v1->v2),(v0->v2),(v3->v4),(v4->v5),(v3->v5))
!
      N=3
!
!  ...6 mixed edges, each with a blending mu function and a locally
!     oriented pair (of nu) representing a projection
!
!  ...e=1 --> edge01 with local orientation v0->v1
      e=1
      MubE(e) = Mu(0); DMubE(1:N,e) = DMu(1:N,0)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0); NupE(1,e) = Nu(1)
      DNupE(1:N,0,e) = DNu(1:N,0); DNupE(1:N,1,e) = DNu(1:N,1)
!  ...e=2 --> edge12 with local orientation v1->v2
      e=2
      MubE(e) = Mu(0); DMubE(1:N,e) = DMu(1:N,0)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(1); NupE(1,e) = Nu(2)
      DNupE(1:N,0,e) = DNu(1:N,1); DNupE(1:N,1,e) = DNu(1:N,2)
!  ...e=3 --> edge20 with local orientation v0->v2
      e=3
      MubE(e) = Mu(0); DMubE(1:N,e) = DMu(1:N,0)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0); NupE(1,e) = Nu(2)
      DNupE(1:N,0,e) = DNu(1:N,0); DNupE(1:N,1,e) = DNu(1:N,2)
!  ...e=4 --> edge34 with local orientation v3->v4
      e=4
      MubE(e) = Mu(1); DMubE(1:N,e) = DMu(1:N,1)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0); NupE(1,e) = Nu(1)
      DNupE(1:N,0,e) = DNu(1:N,0); DNupE(1:N,1,e) = DNu(1:N,1)
!  ...e=5 --> edge45 with local orientation v4->v5
      e=5
      MubE(e) = Mu(1); DMubE(1:N,e) = DMu(1:N,1)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(1); NupE(1,e) = Nu(2)
      DNupE(1:N,0,e) = DNu(1:N,1); DNupE(1:N,1,e) = DNu(1:N,2)
!  ...e=6 --> edge53 with local orientation v3->v5
      e=6
      MubE(e) = Mu(1); DMubE(1:N,e) = DMu(1:N,1)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0); NupE(1,e) = Nu(2)
      DNupE(1:N,0,e) = DNu(1:N,0); DNupE(1:N,1,e) = DNu(1:N,2)
!
!  ...projected coordinates are Nu, so IdecME=false for all edges
      IdecME = .FALSE.
!
   end subroutine BlendProjectPrisME
!----------------------------------------------------------------------
   subroutine BlendProjectPrisQE(Mu,DMu,Nu,DNu, &
                                          NubE,DNubE,MupE,DMupE,IdecQE)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecQE
      double precision, intent(in)  :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(in)  :: Nu(0:2),DNu(1:3,0:2)
      double precision, intent(out) :: NubE(1:3),DNubE(1:3,1:3)
      double precision, intent(out) :: MupE(0:1,1:3),DMupE(1:3,0:1,1:3)
!
!  ...Info from module element_data - coordinates,connectivities:
!        V=((0,0,0),(1,0,0),(0,1,0),
!                    (0,0,1),(1,0,1),(0,1,1))=>(v0,v1,v2,v3,v4,v5)
!                 E_part2=>((v0->v3),(v1->v4),(v2->v5))
!
      N=3
!
!  ...3 quad edges, each with a blending nu function and a locally
!     oriented pair (of mu) representing a projection
!
!  ...e=1 --> edge03 with local orientation v0->v3
      e=1
      NubE(e) = Nu(0); DNubE(1:N,e) = DNu(1:N,0)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0); MupE(1,e) = Mu(1)
      DMupE(1:N,0,e) = DMu(1:N,0); DMupE(1:N,1,e) = DMu(1:N,1)
!  ...e=2 --> edge14 with local orientation v1->v4
      e=2
      NubE(e) = Nu(1); DNubE(1:N,e) = DNu(1:N,1)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0); MupE(1,e) = Mu(1)
      DMupE(1:N,0,e) = DMu(1:N,0); DMupE(1:N,1,e) = DMu(1:N,1)
!  ...e=3 --> edge25 with local orientation v2->v5
      e=3
      NubE(e) = Nu(2); DNubE(1:N,e) = DNu(1:N,2)
!     ...locally oriented pair representing projection
      MupE(0,e) = Mu(0); MupE(1,e) = Mu(1)
      DMupE(1:N,0,e) = DMu(1:N,0); DMupE(1:N,1,e) = DMu(1:N,1)
!
!  ...projected coordinates are Mu, so IdecQE=true for all edges
      IdecQE = .TRUE.
!
   end subroutine BlendProjectPrisQE
!----------------------------------------------------------------------
   subroutine BlendProjectPrisTF(Mu,DMu,Nu,DNu, &
                                          MubF,DMubF,NupF,DNupF,IdecTF)
!
      implicit none
      integer :: N,f
      logical, intent(out) :: IdecTF
      double precision, intent(in)  :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(in)  :: Nu(0:2),DNu(1:3,0:2)
      double precision, intent(out) :: MubF(1:2),DMubF(1:3,1:2)
      double precision, intent(out) :: NupF(0:2,1:2),DNupF(1:3,0:2,1:2)
!
!  ...Info from module element_data - coordinates,connectivities:
!        V=((0,0,0),(1,0,0),(0,1,0),
!                    (0,0,1),(1,0,1),(0,1,1))=>(v0,v1,v2,v3,v4,v5)
!                F_part1=>((v0->v1->v2),(v3->v4->v5)
!
      N=3
!
!  ...2 triangle faces, each with a blending mu function and a locally
!     oriented triplet (of nu) representing a projection
!
!  ...f=1 --> face012 with local orientation v0->v1->v2
      f=1
      MubF(f) = Mu(0); DMubF(1:N,f) = DMu(1:N,0)
!     ...locally oriented triplet representing projection
      NupF(0,f) = Nu(0); NupF(1,f) = Nu(1); NupF(2,f) = Nu(2)
      DNupF(1:N,0,f) = DNu(1:N,0); DNupF(1:N,1,f) = DNu(1:N,1);
                                         DNupF(1:N,2,f) = DNu(1:N,2)
!  ...f=2 --> face345 with local orientation v3->v4->v5
      f=2
      MubF(f) = Mu(1); DMubF(1:N,f) = DMu(1:N,1)
!     ...locally oriented triplet representing projection
      NupF(0,f) = Nu(0); NupF(1,f) = Nu(1); NupF(2,f) = Nu(2)
      DNupF(1:N,0,f) = DNu(1:N,0); DNupF(1:N,1,f) = DNu(1:N,1);
                                         DNupF(1:N,2,f) = DNu(1:N,2)
!
!  ...projected coordinates are Nu, so IdecTF=true for all faces
      IdecTF = .TRUE.
!
   end subroutine BlendProjectPrisTF
!----------------------------------------------------------------------
   subroutine ProjectPrisQF(Mu,DMu,Nu,DNu, STpF,DSTpF,IdecQF)
!
      implicit none
      integer :: N,f
      logical, intent(out) :: IdecQF(1:2,1:3)
      double precision, intent(in)  :: Mu(0:1),DMu(1:3,0:1)
      double precision, intent(in)  :: Nu(0:2),DNu(1:3,0:2)
      double precision, intent(out) :: STpF(0:1,1:2,1:3), &
                                                 DSTpF(1:3,0:1,1:2,1:3)
!
!  ...Info from module element_data - coordinates,connectivities:
!        V=((0,0,0),(1,0,0),(0,1,0),
!                    (0,0,1),(1,0,1),(0,1,1))=>(v0,v1,v2,v3,v4,v5)
!      F_part2=>((v0->v1->v4->v3),(v1->v2->v5->v4),(v0->v2->v5->v3))
!
      N=3
!
!  ...3 quad faces, each with a locally oriented quadruple representing
!     a projection
!  ...simplification flags depend on the local quad face orientations
!
!  ...f=1 --> face0143 with local orientation v0->v1->v4->v3
      f=1
!     ...locally oriented quadruple representing projection
      STpF(0,1,f) = Nu(0); STpF(1,1,f) = Nu(1)
      STpF(0,2,f) = Mu(0); STpF(1,2,f) = Mu(1)
      DSTpF(1:N,0,1,f) = DNu(1:N,0); DSTpF(1:N,1,1,f) = DNu(1:N,1)
      DSTpF(1:N,0,2,f) = DMu(1:N,0); DSTpF(1:N,1,2,f) = DMu(1:N,1)
!     ...simplification flags: projection (Nu;Mu)=>(false,true)
      IdecQF(1,f) = .FALSE.; IdecQF(2,f) = .TRUE.
!  ...f=2 --> face1254 with local orientation v1->v2->v5->v4
      f=2
!     ...locally oriented quadruple representing projection
      STpF(0,1,f) = Nu(1); STpF(1,1,f) = Nu(2)
      STpF(0,2,f) = Mu(0); STpF(1,2,f) = Mu(1)
      DSTpF(1:N,0,1,f) = DNu(1:N,1); DSTpF(1:N,1,1,f) = DNu(1:N,2)
      DSTpF(1:N,0,2,f) = DMu(1:N,0); DSTpF(1:N,1,2,f) = DMu(1:N,1)
!     ...simplification flags: projection (Nu;Mu)=>(false,true)
      IdecQF(1,f) = .FALSE.; IdecQF(2,f) = .TRUE.
!  ...f=3 --> face0253 with local orientation v0->v2->v5->v3
      f=3
!     ...locally oriented quadruple representing projection
      STpF(0,1,f) = Nu(0); STpF(1,1,f) = Nu(2)
      STpF(0,2,f) = Mu(0); STpF(1,2,f) = Mu(1)
      DSTpF(1:N,0,1,f) = DNu(1:N,0); DSTpF(1:N,1,1,f) = DNu(1:N,2)
      DSTpF(1:N,0,2,f) = DMu(1:N,0); DSTpF(1:N,1,2,f) = DMu(1:N,1)
!     ...simplification flags: projection (Nu;Mu)=>(false,true)
      IdecQF(1,f) = .FALSE.; IdecQF(2,f) = .TRUE.
!
   end subroutine ProjectPrisQF
!----------------------------------------------------------------------
   subroutine BlendPyraV(Lam,DLam, LambV,DLambV)
!
      implicit none
      integer :: N,v
      double precision, intent(in)  :: Lam(1:5),DLam(1:3,1:5)
      double precision, intent(out) :: LambV(1:5),DLambV(1:3,1:5)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1))=>(v1,v2,v3,v4,v5)
!
      N=3
!
!  ...5 vertices, each with one blending function
!
!  ...v=1 --> v1=(0,0,0)
      v=1
      LambV(v) = Lam(1); DLambV(1:N,v) = DLam(1:N,1)
!  ...v=2 --> v2=(1,0,0)
      v=2
      LambV(v) = Lam(2); DLambV(1:N,v) = DLam(1:N,2)
!  ...v=3 --> v3=(1,1,0)
      v=3
      LambV(v) = Lam(3); DLambV(1:N,v) = DLam(1:N,3)
!  ...v=4 --> v4=(0,1,0)
      v=4
      LambV(v) = Lam(4); DLambV(1:N,v) = DLam(1:N,4)
!  ...v=5 --> v5=(0,0,1)
      v=5
      LambV(v) = Lam(5); DLambV(1:N,v) = DLam(1:N,5)
!
   end subroutine BlendPyraV
!----------------------------------------------------------------------
   subroutine BlendProjectPyraME(Mu,DMu,Nu,DNu, &
                                          MubE,DMubE,NupE,DNupE,IdecME)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecME
      double precision, intent(in)  :: Mu(0:1,1:2),DMu(1:3,0:1,1:2)
      double precision, intent(in)  :: Nu(0:2,1:2),DNu(1:3,0:2,1:2)
      double precision, intent(out) :: MubE(1:4),DMubE(1:3,1:4)
      double precision, intent(out) :: NupE(0:1,1:4),DNupE(1:3,0:1,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1))=>(v1,v2,v3,v4,v5)
!             E_part1=>((v1->v2),(v2->v3),(v4->v3),(v1->v4))
!
      N=3
!
!  ...4 edges, each with a blending function (mu) and a locally
!     oriented pair (nu) representing a projection
!
!  ...e=1 --> edge12 with local orientation v1->v2
      e=1
      MubE(e) = Mu(0,2); DMubE(1:N,e) = DMu(1:N,0,2)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0,1); NupE(1,e) = Nu(1,1)
      DNupE(1:N,0,e) = DNu(1:N,0,1); DNupE(1:N,1,e) = DNu(1:N,1,1)
!  ...e=2 --> edge23 with local orientation v2->v3
      e=2
      MubE(e) = Mu(1,1); DMubE(1:N,e) = DMu(1:N,1,1)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0,2); NupE(1,e) = Nu(1,2)
      DNupE(1:N,0,e) = DNu(1:N,0,2); DNupE(1:N,1,e) = DNu(1:N,1,2)
!  ...e=3 --> edge34 with local orientation v4->v3
      e=3
      MubE(e) = Mu(1,2); DMubE(1:N,e) = DMu(1:N,1,2)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0,1); NupE(1,e) = Nu(1,1)
      DNupE(1:N,0,e) = DNu(1:N,0,1); DNupE(1:N,1,e) = DNu(1:N,1,1)
!  ...e=4 --> edge41 with local orientation v1->v4
      e=4
      MubE(e) = Mu(0,1); DMubE(1:N,e) = DMu(1:N,0,1)
!     ...locally oriented pair representing projection
      NupE(0,e) = Nu(0,2); NupE(1,e) = Nu(1,2)
      DNupE(1:N,0,e) = DNu(1:N,0,2); DNupE(1:N,1,e) = DNu(1:N,1,2)
!
!  ...projected coordinates are Nu, so IdecME=false for all edges
      IdecME = .FALSE.
!
   end subroutine BlendProjectPyraME
!----------------------------------------------------------------------
   subroutine ProjectPyraTE(Lam,DLam, LampE,DLampE,IdecTE)
!
      implicit none
      integer :: N,e
      logical, intent(out) :: IdecTE
      double precision, intent(in)  :: Lam(1:5),DLam(1:3,1:5)
      double precision, intent(out) :: LampE(0:1,1:4), &
                                                   DLampE(1:3,0:1,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1))=>(v1,v2,v3,v4,v5)
!             E_part2=>((v1->v5),(v2->v5),(v3->v5),(v4->v5))
!
      N=3
!
!  ...4 edges, each with a locally oriented pair (lam) representing
!     a projection
!
!  ...e=1 --> edge15 with local orientation v1->v5
      e=1
      LampE(0,e) = Lam(1); LampE(1,e) = Lam(5)
      DLampE(1:N,0,e) = DLam(1:N,1); DLampE(1:N,1,e) = DLam(1:N,5)
!  ...e=2 --> edge25 with local orientation v2->v5
      e=2
      LampE(0,e) = Lam(2); LampE(1,e) = Lam(5)
      DLampE(1:N,0,e) = DLam(1:N,2); DLampE(1:N,1,e) = DLam(1:N,5)
!  ...e=3 --> edge35 with local orientation v3->v5
      e=3
      LampE(0,e) = Lam(3); LampE(1,e) = Lam(5)
      DLampE(1:N,0,e) = DLam(1:N,3); DLampE(1:N,1,e) = DLam(1:N,5)
!  ...e=4 --> edge45 with local orientation v4->v5
      e=4
      LampE(0,e) = Lam(4); LampE(1,e) = Lam(5)
      DLampE(1:N,0,e) = DLam(1:N,4); DLampE(1:N,1,e) = DLam(1:N,5)
!
!  ...projected coordinates are Lam, so IdecTE=false for all edges
      IdecTE = .FALSE.
!
   end subroutine ProjectPyraTE
!----------------------------------------------------------------------
   subroutine ProjectPyraQF(Mu,DMu, MupF,DMupF,IdecQF)
!
      implicit none
      integer :: N
      logical, intent(out) :: IdecQF(1:2)
      double precision, intent(in)  :: Mu(0:1,1:2),DMu(1:3,0:1,1:2)
      double precision, intent(out) :: MupF(0:1,1:2),DMupF(1:3,0:1,1:2)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1))=>(v1,v2,v3,v4,v5)
!             F_part1=>((v1->v2->v3->v4))
!
      N=3
!
!  ...1 quadrilateral face with  a locally oriented quadruple
!     representing a projection
!
!  ...face1234 with local orientation v1->v2->v3->v4
!     ...locally oriented quadruple representing projection
      MupF(0,1) = Mu(0,1); MupF(1,1) = Mu(1,1)
      MupF(0,2) = Mu(0,2); MupF(1,2) = Mu(1,2)
      DMupF(1:N,0,1) = DMu(1:N,0,1); DMupF(1:N,1,1) = DMu(1:N,1,1)
      DMupF(1:N,0,2) = DMu(1:N,0,2); DMupF(1:N,1,2) = DMu(1:N,1,2)
!
!  ...projected coordinates are (Mu;Mu), so IdecQF=(true,true) for
!     the face
      IdecQF(1) = .TRUE.; IdecQF(2) = .TRUE.
!
   end subroutine ProjectPyraQF
!----------------------------------------------------------------------
   subroutine BlendProjectPyraTF(Mu,DMu,Nu,DNu, &
                                          MubF,DMubF,NupF,DNupF,IdecTF)
!
      implicit none
      integer :: N,f
      logical, intent(out) :: IdecTF
      double precision, intent(in)  :: Mu(0:1,1:2),DMu(1:3,0:1,1:2)
      double precision, intent(in)  :: Nu(0:2,1:2),DNu(1:3,0:2,1:2)
      double precision, intent(out) :: MubF(1:4),DMubF(1:3,1:4)
      double precision, intent(out) :: NupF(0:2,1:4),DNupF(1:3,0:2,1:4)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1))=>(v1,v2,v3,v4,v5)
!      F_part2=>((v1->v2->v5),(v2->v3->v5),(v4->v3->v5),(v1->v4->v5))
!
      N=3
!
!  ...4 triangle faces, each with a blending function (mu) and a
!      locally oriented triplet (nu) representing a projection
!
!  ...f=1 --> face125 with local orientation v1->v2->v5
      f=1
      MubF(f) = Mu(0,2); DMubF(1:N,f) = DMu(1:N,0,2)
!     ...locally oriented pair representing projection
      NupF(0,f) = Nu(0,1); NupF(1,f) = Nu(1,1); NupF(2,f) = Nu(2,1)
      DNupF(1:N,0,f) = DNu(1:N,0,1); DNupF(1:N,1,f) = DNu(1:N,1,1);
                                          DNupF(1:N,2,f) = DNu(1:N,2,1)
!  ...f=2 --> face235 with local orientation v2->v3->v5
      f=2
      MubF(f) = Mu(1,1); DMubF(1:N,f) = DMu(1:N,1,1)
!     ...locally oriented pair representing projection
      NupF(0,f) = Nu(0,2); NupF(1,f) = Nu(1,2); NupF(2,f) = Nu(2,2)
      DNupF(1:N,0,f) = DNu(1:N,0,2); DNupF(1:N,1,f) = DNu(1:N,1,2);
                                          DNupF(1:N,2,f) = DNu(1:N,2,2)
!  ...f=3 --> face345 with local orientation v4->v3->v5
      f=3
      MubF(f) = Mu(1,2); DMubF(1:N,f) = DMu(1:N,1,2)
!     ...locally oriented pair representing projection
      NupF(0,f) = Nu(0,1); NupF(1,f) = Nu(1,1); NupF(2,f) = Nu(2,1)
      DNupF(1:N,0,f) = DNu(1:N,0,1); DNupF(1:N,1,f) = DNu(1:N,1,1);
                                          DNupF(1:N,2,f) = DNu(1:N,2,1)
!  ...f=4 --> face415 with local orientation v1->v4->v5
      f=4
      MubF(f) = Mu(0,1); DMubF(1:N,f) = DMu(1:N,0,1)
!     ...locally oriented pair representing projection
      NupF(0,f) = Nu(0,2); NupF(1,f) = Nu(1,2); NupF(2,f) = Nu(2,2)
      DNupF(1:N,0,f) = DNu(1:N,0,2); DNupF(1:N,1,f) = DNu(1:N,1,2);
                                          DNupF(1:N,2,f) = DNu(1:N,2,2)
!
!  ...projected coordinates are Nu, so IdecTF=true for all faces
      IdecTF = .TRUE.
!
   end subroutine BlendProjectPyraTF
!----------------------------------------------------------------------

   subroutine ProjectPyraLamTF(Lam,DLam, LampF,DLampF,IdecTF)
!
      implicit none
      integer :: N,f
      logical, intent(out) :: IdecTF
      double precision, intent(in)  :: Lam(5),DLam(1:3,5)
      double precision, intent(out) :: LampF(0:2,4),DLampF(3,0:2,4)
!
!  ...Info from module element_data - coordinates,connectivities:
!      V=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1))=>(v1,v2,v3,v4,v5)
!      F_part2=>((v1->v2->v5),(v2->v3->v5),(v4->v3->v5),(v1->v4->v5))
!
      N=3
!
!  ...f=1 --> face125 with local orientation v1->v2->v5
      f=1
!     ...locally oriented pair representing projection
      LampF(0,f) = Lam(1); LampF(1,f) = Lam(2); LampF(2,f) = Lam(5)
      DLampF(1:N,0,f) = DLam(1:N,1); DLampF(1:N,1,f) = DLam(1:N,2);
                                          DLampF(1:N,2,f) = DLam(1:N,5)
!  ...f=2 --> face235 with local orientation v2->v3->v5
      f=2
!     ...locally oriented pair representing projection
      LampF(0,f) = Lam(2); LampF(1,f) = Lam(3); LampF(2,f) = Lam(5)
      DLampF(1:N,0,f) = DLam(1:N,2); DLampF(1:N,1,f) = DLam(1:N,3);
                                          DLampF(1:N,2,f) = DLam(1:N,5)
!  ...f=3 --> face345 with local orientation v4->v3->v5
      f=3
!     ...locally oriented pair representing projection
      LampF(0,f) = Lam(4); LampF(1,f) = Lam(3); LampF(2,f) = Lam(5)
      DLampF(1:N,0,f) = DLam(1:N,4); DLampF(1:N,1,f) = DLam(1:N,3);
                                          DLampF(1:N,2,f) = DLam(1:N,5)
!  ...f=4 --> face415 with local orientation v1->v4->v5
      f=4
!     ...locally oriented pair representing projection
      LampF(0,f) = Lam(1); LampF(1,f) = Lam(4); LampF(2,f) = Lam(5)
      DLampF(1:N,0,f) = DLam(1:N,1); DLampF(1:N,1,f) = DLam(1:N,4);
                                          DLampF(1:N,2,f) = DLam(1:N,5)
!
!  ...projected coordinates are Lam, so IdecTF=false for all faces
      IdecTF = .FALSE.
!
   end subroutine ProjectPyraLamTF
