! Routines:
!  - shape1DHSeg
!  - shape1DQSeg
!
!----------------------------------------------------------------------
!
!     routine name      - shape1DHSeg
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14, Apr 17
!
!     purpose:          - routine returns values of 1D Segment H1
!                         shape functions
!
!     arguments:
!
!     in:
!       Xi              - master segment coordinate
!       Nord            - polynomial order for the nodes (H1 sense)
!       Nsize           - relevant sizes of local arrays
!
!     out:
!       NrdofH          - number of dof
!       ShapH           - values of the shape functions
!       GradH           - gradient of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape1DHSeg(Xi,Nord,Nsize, NrdofH,ShapH,GradH)
!
      implicit none
      integer, intent(in ) :: Nord,Nsize(2)
      integer, intent(out) :: NrdofH
      integer :: N,m,v,ndofE,minI,maxI,i
      logical :: IdecE
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapH(Nsize(2)),GradH(Nsize(2))
      double precision :: Mu(0:1),DMu(0:1),MubV(2),DMubV(2)
      double precision :: phiE(2:Nsize(1)),DphiE(2:Nsize(1))
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=1
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineSegment(Xi, Mu,DMu)
!
!  ...VERTEX SHAPE FUNCTIONS
      call BlendSegV(Mu,DMu, MubV,DMubV)
      do v=1,2
        m=m+1
!
        ShapH(m) = MubV(v)
        GradH(m) = DMubV(v)
      enddo
!
!  ...BUBBLE FUNCTIONS
      ndofE = Nord-1
      if (ndofE.gt.0) then
!    ...local parameters
        minI  = 2
        maxI  = Nord
        IdecE = .TRUE.
!    ...construct the shape functions
        call AncPhiE(Mu,DMu,Nord,IdecE,N, &
                                      phiE(minI:maxI),DphiE(minI:maxI))
        do i=minI,maxI
          m=m+1
!
          ShapH(m) = phiE(i)
          GradH(m) = DphiE(i)
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofH = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord
 7001   format('shape1DHSeg: Xi = ',f8.3,/, &
               'Norder  = ',i2)
!
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
#endif
!
   end subroutine shape1DHSeg
!
!
!----------------------------------------------------------------------
!
!     routine name      - shape1DQSeg
!
!----------------------------------------------------------------------
!
!     latest revision:  - Nov 14, Aor 17
!
!     purpose:          - routine returns values of 1D Segment L2
!                         shape functions
!
!     arguments:
!
!     in:
!       Xi              - master segment coordinate
!       Nord            - polynomial order for the nodes (H1 sense)
!       Nsize           - relevant sizes of local arrays
!
!     out:
!       NrdofQ          - number of dof
!       ShapQ           - values of the shape functions
!
!----------------------------------------------------------------------
!
   subroutine shape1DQSeg(Xi,Nord,Nsize, NrdofQ,ShapQ)
!
      implicit none
      integer, intent(in ) :: Nord,Nsize(2)
      integer, intent(out) :: NrdofQ
      integer :: N,m,ndofE,minI,maxI,i
      double precision, intent(in)  :: Xi
      double precision, intent(out) :: ShapQ(Nsize(2))
      double precision :: Mu(0:1),DMu(0:1),homP(0:Nsize(1)-1)
!
#if HP3D_DEBUG
!  ...debugging flag
      integer :: iprint
      iprint=0
#endif
!
!  ...spatial dimensions
      N=1
!
!  ...initiate counter for shape functions
      m=0
!
!  ...Define affine coordinates and gradients
      call AffineSegment(Xi, Mu,DMu)
!
!  ...EDGE FUNCTIONS
      ndofE = Nord
      if (ndofE.gt.0) then
!    ...local parameters
        minI = 0
        maxI = Nord-1
!    ...construct the shape functions
        call HomLegendre(Mu,maxI, homP(minI:maxI))
        do i=minI,maxI
          m=m+1
!
          ShapQ(m) = homP(i)
        enddo
      endif
!
!  ...give total degrees of freedom
      NrdofQ = m
!
#if HP3D_DEBUG
!  ...print this when debugging
      if (iprint.eq.1) then
        write(*,7001) Xi,Nord
 7001   format('shape1DQSeg: Xi = ',f8.3,/, &
               'Norder  = ',i2)
!
        if (ndofE.gt.0) then
          write(*,*) 'EDGE FUNCTIONS = '
          do m=1,ndofE
            write(*,7002) m,ShapQ(m)
 7002       format('k = ',i3,' ShapQ = ',e12.5)
          enddo
        endif
        call pause
      endif
#endif
!
   end subroutine shape1DQSeg
!
