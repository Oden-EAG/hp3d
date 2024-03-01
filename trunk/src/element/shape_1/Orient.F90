! Routines:
!  - OrientE
!  - OrientQuad
!  - OrientTri
!----------------------------------------------------------------------
! Routines representing the local to global transformations of edges,
! triangle faces and quadrilateral faces
!----------------------------------------------------------------------
   subroutine OrientE(S,DS,Nori,N, GS,GDS)
!
      implicit none
      integer, intent(in) :: Nori, N
      integer :: Or(0:1)
      double precision, intent(in)  :: S(0:1),DS(1:N,0:1)
      double precision, intent(out) :: GS(0:1),GDS(1:N,0:1)
!
!     Or(1) - Is the global axis aligned with the (parallel) local one?
!
!     GS(0:1)=S(Or(Nori,0:1))
!
      select case (Nori)
!     ...Nori=0 => (s0,s1)->(s0,s1)
         case(0)
            Or(0) = 0; Or(1) = 1
!     ...Nori=1 => (s0,s1)->(s1,s0)
         case(1)
            Or(0) = 1; Or(1) = 0
         case default
            write(*,*) "Error: invalid orientation in OrientE"
            stop 1
      end select
!
!  ...Local-to-global transformation
      GS(0) = S(Or(0));           GS(1) = S(Or(1))
      GDS(1:N,0) = DS(1:N,Or(0)); GDS(1:N,1) = DS(1:N,Or(1))
!
   end subroutine OrientE
!----------------------------------------------------------------------
   subroutine OrientQuad(ST,DST,Nori,Idec,N, GST,GDST,GIdec)
!
      implicit none
      integer, intent(in) :: Nori, N
      integer :: OrPa(1:2),OrGS(0:1),OrGT(0:1)
      logical, intent(in)  :: Idec(2)
      logical, intent(out) :: GIdec(2)
      double precision, intent(in)  :: ST(0:1,1:2),DST(1:N,0:1,1:2)
      double precision, intent(out) :: GST(0:1,1:2),GDST(1:N,0:1,1:2)
!
!     OrPa - Order of the pairs S and T (swapping)
!     OrGS(1) - Is the global S axis aligned with the (parallel) local?
!     OrGT(1) - Is the global T axis aligned with the (parallel) local?
!
!     GST(1,0:1)=GS(0:1)=ST(OrPa(Nori,1),OrGS(Nori,0:1))
!     GST(2,0:1)=GT(0:1)=ST(OrPa(Nori,2),OrGT(Nori,0:1))
!
      select case (Nori)
!     ...Nori=0 => ((s0,s1),(t0,t1))->((s0,s1),(t0,t1))
         case (0)
            OrPa(1) = 1; OrPa(2) = 2
            OrGS(0) = 0; OrGS(1) = 1
            OrGT(0) = 0; OrGT(1) = 1
!     ...Nori=1 => ((s0,s1),(t0,t1))->((t0,t1),(s1,s0))
         case (1)
            OrPa(1) = 2; OrPa(2) = 1
            OrGS(0) = 0; OrGS(1) = 1
            OrGT(0) = 1; OrGT(1) = 0
!     ...Nori=2 => ((s0,s1),(t0,t1))->((s1,s0),(t1,t0))
         case(2)
            OrPa(1) = 1; OrPa(2) = 2
            OrGS(0) = 1; OrGS(1) = 0
            OrGT(0) = 1; OrGT(1) = 0
!     ...Nori=3 => ((s0,s1),(t0,t1))->((t1,t0),(s0,s1))
         case(3)
            OrPa(1) = 2; OrPa(2) = 1
            OrGS(0) = 1; OrGS(1) = 0
            OrGT(0) = 0; OrGT(1) = 1
!     ...Nori=4 => ((s0,s1),(t0,t1))->((t0,t1),(s0,s1))
         case(4)
            OrPa(1) = 2; OrPa(2) = 1
            OrGS(0) = 0; OrGS(1) = 1
            OrGT(0) = 0; OrGT(1) = 1
!     ...Nori=5 => ((s0,s1),(t0,t1))->((s1,s0),(t0,t1))
         case(5)
            OrPa(1) = 1; OrPa(2) = 2
            OrGS(0) = 1; OrGS(1) = 0
            OrGT(0) = 0; OrGT(1) = 1
!     ...Nori=6 => ((s0,s1),(t0,t1))->((t1,t0),(s1,s0))
         case(6)
            OrPa(1) = 2; OrPa(2) = 1
            OrGS(0) = 1; OrGS(1) = 0
            OrGT(0) = 1; OrGT(1) = 0
!     ...Nori=7 => ((s0,s1),(t0,t1))->((s0,s1),(t1,t0))
         case(7)
            OrPa(1) = 1; OrPa(2) = 2
            OrGS(0) = 0; OrGS(1) = 1
            OrGT(0) = 1; OrGT(1) = 0
         case default
            write(*,*) "Error: invalid orientation in OrientQuad"
            stop 1
      end select
!
!     GST=[GST(1,0),GST(1,1); GST(2,0),GST(2,1)]
!  ...Local-to-global transformation
      GST(0,1) = ST(OrGS(0),OrPa(1))
      GST(1,1) = ST(OrGS(1),OrPa(1))
      GST(0,2) = ST(OrGT(0),OrPa(2))
      GST(1,2) = ST(OrGT(1),OrPa(2))
!
      GDST(1:N,0,1) = DST(1:N,OrGS(0),OrPa(1))
      GDST(1:N,1,1) = DST(1:N,OrGS(1),OrPa(1))
      GDST(1:N,0,2) = DST(1:N,OrGT(0),OrPa(2))
      GDST(1:N,1,2) = DST(1:N,OrGT(1),OrPa(2))
!
      GIdec(1) = Idec(OrPa(1)); GIdec(2) = Idec(OrPa(2))
!
   end subroutine OrientQuad
!
!----------------------------------------------------------------------
   subroutine OrientTri(S,DS,Nori,N, GS,GDS)
!
      implicit none
      integer, intent(in) :: Nori, N
      integer :: Or(0:2)
      double precision, intent(in)  :: S(0:2),DS(1:N,0:2)
      double precision, intent(out) :: GS(0:2),GDS(1:N,0:2)
!
!     GS(0:2)=S(Or(Nori,0:2))
!
      select case (Nori)
!     ...Nori=0 => (s0,s1,s2)->(s0,s1,s2)
         case(0)
            Or(0) = 0; Or(1) = 1; Or(2) = 2
!     ...Nori=1 => (s0,s1,s2)->(s1,s2,s0)
         case(1)
            Or(0) = 1; Or(1) = 2; Or(2) = 0
!     ...Nori=2 => (s0,s1,s2)->(s2,s0,s1)
         case(2)
            Or(0) = 2; Or(1) = 0; Or(2) = 1
!     ...Nori=3 => (s0,s1,s2)->(s0,s2,s1)
         case(3)
            Or(0) = 0; Or(1) = 2; Or(2) = 1
!     ...Nori=4 => (s0,s1,s2)->(s1,s0,s2)
         case(4)
            Or(0) = 1; Or(1) = 0; Or(2) = 2
!     ...Nori=5 => (s0,s1,s2)->(s2,s1,s0)
         case(5)
            Or(0) = 2; Or(1) = 1; Or(2) = 0
         case default
            write(*,*) "Error: invalid orientation in OrientTri"
            stop 1
      end select
!
!  ...Local-to-global transformation
      GS(0) = S(Or(0))
      GS(1) = S(Or(1))
      GS(2) = S(Or(2))
!
      GDS(1:N,0) = DS(1:N,Or(0))
      GDS(1:N,1) = DS(1:N,Or(1))
      GDS(1:N,2) = DS(1:N,Or(2))
!
   end subroutine OrientTri
