c Routines:
c  - OrientE
c  - OrientQuad
c  - OrientTri
c----------------------------------------------------------------------
c Routines representing the local to global transformations of edges,
c triangle faces and quadrilateral faces
c----------------------------------------------------------------------
      subroutine OrientE(S,DS,Nori,N, GS,GDS)
c
      implicit none
      integer, intent(in) :: Nori, N
      integer :: Or(0:1,0:1)
      double precision, intent(in)  :: S(0:1),DS(1:N,0:1)
      double precision, intent(out) :: GS(0:1),GDS(1:N,0:1)
c
c     Or(1) - Is the global axis aligned with the (parallel) local one?
c
c     GS(0:1)=S(Or(Nori,0:1))
c
c  ...Nori=0 => (s0,s1)->(s0,s1)
      Or(0,0) = 0; Or(0,1) = 1
c  ...Nori=1 => (s0,s1)->(s1,s0)
      Or(1,0) = 1; Or(1,1) = 0
c
c  ...Local-to-global transformation
      GS(0) = S(Or(Nori,0));           GS(1) = S(Or(Nori,1))
      GDS(1:N,0) = DS(1:N,Or(Nori,0)); GDS(1:N,1) = DS(1:N,Or(Nori,1))
c
      end subroutine OrientE
c----------------------------------------------------------------------
      subroutine OrientQuad(ST,DST,Nori,Idec,N, GST,GDST,GIdec)
c
      implicit none
      integer, intent(in) :: Nori, N
      integer :: OrPa(0:7,1:2),OrGS(0:7,0:1),OrGT(0:7,0:1)
      logical, intent(in)  :: Idec(2)
      logical, intent(out) :: GIdec(2)
      double precision, intent(in)  :: ST(0:1,1:2),DST(1:N,0:1,1:2)
      double precision, intent(out) :: GST(0:1,1:2),GDST(1:N,0:1,1:2)
c
c     OrPa - Order of the pairs S and T (swapping)
c     OrGS(1) - Is the global S axis aligned with the (parallel) local?
c     OrGT(1) - Is the global T axis aligned with the (parallel) local?
c
c     GST(1,0:1)=GS(0:1)=ST(OrPa(Nori,1),OrGS(Nori,0:1))
c     GST(2,0:1)=GT(0:1)=ST(OrPa(Nori,2),OrGT(Nori,0:1))
c
c  ...Nori=0 => ((s0,s1),(t0,t1))->((s0,s1),(t0,t1))
      OrPa(0,1) = 1; OrPa(0,2) = 2
      OrGS(0,0) = 0; OrGS(0,1) = 1
      OrGT(0,0) = 0; OrGT(0,1) = 1
c  ...Nori=1 => ((s0,s1),(t0,t1))->((t0,t1),(s1,s0))
      OrPa(1,1) = 2; OrPa(1,2) = 1
      OrGS(1,0) = 0; OrGS(1,1) = 1
      OrGT(1,0) = 1; OrGT(1,1) = 0
c  ...Nori=2 => ((s0,s1),(t0,t1))->((s1,s0),(t1,t0))
      OrPa(2,1) = 1; OrPa(2,2) = 2
      OrGS(2,0) = 1; OrGS(2,1) = 0
      OrGT(2,0) = 1; OrGT(2,1) = 0
c  ...Nori=3 => ((s0,s1),(t0,t1))->((t1,t0),(s0,s1))
      OrPa(3,1) = 2; OrPa(3,2) = 1
      OrGS(3,0) = 1; OrGS(3,1) = 0
      OrGT(3,0) = 0; OrGT(3,1) = 1
c  ...Nori=4 => ((s0,s1),(t0,t1))->((t0,t1),(s0,s1))
      OrPa(4,1) = 2; OrPa(4,2) = 1
      OrGS(4,0) = 0; OrGS(4,1) = 1
      OrGT(4,0) = 0; OrGT(4,1) = 1
c  ...Nori=5 => ((s0,s1),(t0,t1))->((s1,s0),(t0,t1))
      OrPa(5,1) = 1; OrPa(5,2) = 2
      OrGS(5,0) = 1; OrGS(5,1) = 0
      OrGT(5,0) = 0; OrGT(5,1) = 1
c  ...Nori=6 => ((s0,s1),(t0,t1))->((t1,t0),(s1,s0))
      OrPa(6,1) = 2; OrPa(6,2) = 1
      OrGS(6,0) = 1; OrGS(6,1) = 0
      OrGT(6,0) = 1; OrGT(6,1) = 0
c  ...Nori=7 => ((s0,s1),(t0,t1))->((s0,s1),(t1,t0))
      OrPa(7,1) = 1; OrPa(7,2) = 2
      OrGS(7,0) = 0; OrGS(7,1) = 1
      OrGT(7,0) = 1; OrGT(7,1) = 0
c
c     GST=[GST(1,0),GST(1,1); GST(2,0),GST(2,1)]
c  ...Local-to-global transformation
      GST(0,1) = ST(OrPa(Nori,1),OrGS(Nori,0))
      GST(1,1) = ST(OrPa(Nori,1),OrGS(Nori,1))
      GST(0,2) = ST(OrPa(Nori,2),OrGT(Nori,0))
      GST(1,2) = ST(OrPa(Nori,2),OrGT(Nori,1))
c
      GDST(1:N,0,1) = DST(1:N,OrPa(Nori,1),OrGS(Nori,0))
      GDST(1:N,1,1) = DST(1:N,OrPa(Nori,1),OrGS(Nori,1))
      GDST(1:N,0,2) = DST(1:N,OrPa(Nori,2),OrGT(Nori,0))
      GDST(1:N,1,2) = DST(1:N,OrPa(Nori,2),OrGT(Nori,1))
c
      GIdec(1) = Idec(OrPa(Nori,1)); GIdec(2) = Idec(OrPa(Nori,2))
c
      end subroutine OrientQuad
c
c----------------------------------------------------------------------
      subroutine OrientTri(S,DS,Nori,N, GS,GDS)
c
      implicit none
      integer, intent(in) :: Nori, N
      integer :: Or(0:5,0:2)
      double precision, intent(in)  :: S(0:2),DS(1:N,0:2)
      double precision, intent(out) :: GS(0:2),GDS(1:N,0:2)
c
c     GS(0:2)=S(Or(Nori,0:2))
c
c  ...Nori=0 => (s0,s1,s2)->(s0,s1,s2)
      Or(0,0) = 0; Or(0,1) = 1; Or(0,2) = 2
c  ...Nori=1 => (s0,s1,s2)->(s1,s2,s0)
      Or(1,0) = 1; Or(1,1) = 2; Or(1,2) = 0
c  ...Nori=2 => (s0,s1,s2)->(s2,s0,s1)
      Or(2,0) = 2; Or(2,1) = 0; Or(2,2) = 1
c  ...Nori=3 => (s0,s1,s2)->(s0,s2,s1)
      Or(3,0) = 0; Or(3,1) = 2; Or(3,2) = 1
c  ...Nori=4 => (s0,s1,s2)->(s1,s0,s2)
      Or(4,0) = 1; Or(4,1) = 0; Or(4,2) = 2
c  ...Nori=5 => (s0,s1,s2)->(s2,s1,s0)
      Or(5,0) = 2; Or(5,1) = 1; Or(5,2) = 0
c
c  ...Local-to-global transformation
      GS(0) = S(Or(Nori,0))
      GS(1) = S(Or(Nori,1))
      GS(2) = S(Or(Nori,2))
c
      GDS(1:N,0) = DS(1:N,Or(Nori,0))
      GDS(1:N,1) = DS(1:N,Or(Nori,1))
      GDS(1:N,2) = DS(1:N,Or(Nori,2))
c
      end subroutine OrientTri
