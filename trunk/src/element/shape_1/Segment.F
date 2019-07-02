c Routines:
c  - shape1DHSeg
c  - shape1DQSeg
c
c----------------------------------------------------------------------
c
c     routine name      - shape1DHSeg
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Apr 17
c
c     purpose:          - routine returns values of 1D Segment H1
c                         shape functions
c
c     arguments:
c
c     in:
c       Xi              - master segment coordinate
c       Nord            - polynomial order for the nodes (H1 sense)
c       Nsize           - relevant sizes of local arrays
c
c     out:
c       NrdofH          - number of dof
c       ShapH           - values of the shape functions
c       GradH           - gradient of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1DHSeg(Xi,Nord,Nsize, NrdofH,ShapH,GradH)
c
      implicit none
      integer, intent(in ) :: Nord,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: iprint,N,m,v,ndofE,minI,maxI,i
      logical :: IdecE
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(Nsize(2)),GradH(Nsize(2))
      double precision :: Mu(0:1),DMu(0:1),MubV(2),DMubV(2)
      double precision :: phiE(2:Nsize(1)),DphiE(2:Nsize(1))
c
c  ...debugging flag
      iprint=0
c
c  ...spatial dimensions
      N=1
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffineSegment(Xi, Mu,DMu)
c
c  ...VERTEX SHAPE FUNCTIONS
      call BlendSegV(Mu,DMu, MubV,DMubV)
      do v=1,2
        m=m+1
c
        ShapH(m) = MubV(v)
        GradH(m) = DMubV(v)
      enddo
c
c  ...BUBBLE FUNCTIONS
      ndofE = Nord-1
      if (ndofE.gt.0) then
c    ...local parameters
        minI  = 2
        maxI  = Nord
        IdecE = .TRUE.
c    ...construct the shape functions
        call AncPhiE(Mu,DMu,Nord,IdecE,N,
     .                                phiE(minI:maxI),DphiE(minI:maxI))
        do i=minI,maxI
          m=m+1
c
          ShapH(m) = phiE(i)
          GradH(m) = DphiE(i)
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofH = m
c
c  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord
 7001   format('shape1DHSeg: Xi = ',f8.3,/,
     .         'Norder  = ',i2)
c
        write(*,*) 'VERTEX SHAPE FUNCTIONS = '
        do v=1,2
          m=v
          write(*,7002) m,ShapH(m),GradH(m)
 7002     format('k = ',i3,' ShapH, GradH = ',e12.5,3x,e12.5)
        enddo
        if (ndofE.gt.0) then
          write(*,*) 'BUBBLE FUNCTIONS = '
          do i=1,ndofE
            m=m+1
            write(*,7002) m,ShapH(m),GradH(m)
          enddo
        endif
        call pause
      endif
c
      end subroutine shape1DHSeg
c
c
c----------------------------------------------------------------------
c
c     routine name      - shape1DQSeg
c
c----------------------------------------------------------------------
c
c     latest revision:  - Nov 14, Aor 17
c
c     purpose:          - routine returns values of 1D Segment L2
c                         shape functions
c
c     arguments:
c
c     in:
c       Xi              - master segment coordinate
c       Nord            - polynomial order for the nodes (H1 sense)
c       Nsize           - relevant sizes of local arrays
c
c     out:
c       NrdofQ          - number of dof
c       ShapQ           - values of the shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape1DQSeg(Xi,Nord,Nsize, NrdofQ,ShapQ)
c
      implicit none
      integer, intent(in ) :: Nord,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: iprint,N,m,ndofE,minI,maxI,i
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Mu(0:1),DMu(0:1),homP(0:Nsize(1)-1)
c
c  ...debugging flag
      iprint=0
c
c  ...spatial dimensions
      N=1
c
c  ...initiate counter for shape functions
      m=0
c
c  ...Define affine coordinates and gradients
      call AffineSegment(Xi, Mu,DMu)
c
c  ...EDGE FUNCTIONS
      ndofE = Nord
      if (ndofE.gt.0) then
c    ...local parameters
        minI = 0
        maxI = Nord-1
c    ...construct the shape functions
        call HomLegendre(Mu,maxI, homP(minI:maxI))
        do i=minI,maxI
          m=m+1
c
          ShapQ(m) = homP(i)
        enddo
      endif
c
c  ...give total degrees of freedom
      NrdofQ = m
c
c  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord
 7001   format('shape1DQSeg: Xi = ',f8.3,/,
     .         'Norder  = ',i2)
c
        if (ndofE.gt.0) then
          write(*,*) 'EDGE FUNCTIONS = '
          do m=1,ndofE
            write(*,7002) m,ShapQ(m)
 7002       format('k = ',i3,' ShapQ = ',e12.5)
          enddo
        endif
        call pause
      endif
c
      end subroutine shape1DQSeg
c
