c----------------------------------------------------------------------
c
c   routine name       - shapHke
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate edge kernel H^1 shape functions
c
c   arguments :
c     in:
c        T             - local edge coordinate
c        Nord          - edge order
c        Norient       - edge orientation
c     out:
c        Nrdof         - number of kernel shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHke(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
      use control          , only : SHAPE_FLAG
      use parameters       , only : MAXP
      use gauss_quadrature , only : INITIALIZED,XLEG0
c
      implicit none
      real(8),                  intent(in ) :: T
      integer,                  intent(in ) :: Nord
      integer,                  intent(in ) :: Norient
      integer,                  intent(out) :: Nrdof
      real(8),dimension(MAXP-1),intent(out) :: ShapH,DshapH
c
      real(8) :: s,dsdt,x,x2,dx
      integer :: iprint,k,i
      real(8),dimension(8) :: v,d
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...initialize XLEG0 if needed
      if (.NOT.INITIALIZED) call init_gauss_quadrature
c
      ShapH(1:MAXP-1)=0.d0 ; DshapH(1:MAXP-1)=0.d0
c
      if (Nord.lt.2) then
        Nrdof=0 ; return
      endif
c
c  ...check maximum order of approximation
      if (Nord.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format(' shapHke: Nord,MAXP = ',2(i3,1x))
        stop
      endif
c
c  ...number of dof's associated to edge
      Nrdof=Nord-1
c
c  ...account for orientation
      select case(Norient)
      case(0) ; s=T      ; dsdt=1.d0
      case(1) ; s=1.d0-T ; dsdt=-1.d0
      endselect
c
      select case(SHAPE_FLAG)
c
c======================================================================
c  "KERNEL" FOR PEANO SHAPE FUNCTIONS                                 |
c======================================================================
      case(1)
        ShapH(1)=1.d0 ; DshapH(1)=0.d0
        if (Nord.eq.2) return
        x=s-0.5d0
        do k=2,Nrdof
          ShapH( k)=ShapH( k-1)*x
          DshapH(k)=DshapH(k-1)*x + ShapH(k-1)
        enddo
c
c======================================================================
c  "KERNEL" FOR INTEGRATED LEGENDRE POLYNOMIALS                       |
c======================================================================
      case(2)
        if (Nord-1.gt.8) then
          write(*,9999)Nord
9999      format(' shapHke: max. supported order of approximation
     . for int. Legendre polynomials exceeded! Nord = ',i2)
          stop
        endif
c
        x=2.d0*(s-0.5d0) ; dx=2.d0 ; x2=x*x
c
        goto (80,70,60,50,40,30,20,10), Nord-1
   10   v(8)=(((     XLEG0(8,8)*x2 +      XLEG0(8,6))*x2 +
     .                                    XLEG0(8,4))*x2 +
     .                                    XLEG0(8,2))*x
        d(8)= ((7.d0*XLEG0(8,8)*x2 + 5.d0*XLEG0(8,6))*x2 +
     .                               3.d0*XLEG0(8,4))*x2 + XLEG0(8,2)
   20   v(7)= ((     XLEG0(7,7)*x2 +      XLEG0(7,5))*x2 +
     .                                    XLEG0(7,3))*x2 + XLEG0(7,1)
        d(7)= ((6.d0*XLEG0(7,7)*x2 + 4.d0*XLEG0(7,5))*x2 +
     .                               2.d0*XLEG0(7,3))*x
   30   v(6)= ((     XLEG0(6,6)*x2 +      XLEG0(6,4))*x2 +
     .                                    XLEG0(6,2))*x
        d(6)=  (5.d0*XLEG0(6,6)*x2 + 3.d0*XLEG0(6,4))*x2 + XLEG0(6,2)
   40   v(5)=  (     XLEG0(5,5)*x2 +      XLEG0(5,3))*x2 + XLEG0(5,1)
        d(5)=  (4.d0*XLEG0(5,5)*x2 + 2.d0*XLEG0(5,3))*x
   50   v(4)=  (     XLEG0(4,4)*x2 +      XLEG0(4,2))*x
        d(4)=   3.d0*XLEG0(4,4)*x2 +      XLEG0(4,2)
   60   v(3)=        XLEG0(3,3)*x2 +      XLEG0(3,1)
        d(3)=   2.d0*XLEG0(3,3)*x
   70   v(2)=        XLEG0(2,2)*x
        d(2)=        XLEG0(2,2)
   80   v(1)=        XLEG0(1,1)
        d(1)=   0.d0
c
c  .....store and adjust derivatives
        do i=1,Nrdof
          ShapH(i)=v(i) ; DshapH(i)=d(i)*dx
        enddo
c
      endselect
c
c  ...account for orientation
      do i=1,Nrdof ; DshapH(i)=DshapH(i)*dsdt ; enddo
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7002) T,Nord,Norient
 7002   format(' shapHke: T,Nord,Norient = ',f8.3,2x,i2,2x,i1,
     .         ' Shap,Dvshap = ')
        do k=1,Nrdof
          write(*,7003) k,ShapH(k),DshapH(k)
 7003     format('         k = ',i2,2f8.3)
        enddo
        call pause
      endif
c
c
      end subroutine shapHke
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapHbe_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate edge bubble H^1 shape functions
c
c   arguments :
c     in:
c        T             - local edge coordinate
c        Nord          - edge order
c        Norient       - edge orientation
c     out:
c        Nrdof         - number of shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHbe_0(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
      use parameters , only : MAXP
c
      implicit none
      real(8),                  intent(in ) :: T
      integer,                  intent(in ) :: Nord
      integer,                  intent(in ) :: Norient
      integer,                  intent(out) :: Nrdof
      real(8),dimension(MAXP-1),intent(out) :: ShapH,DshapH
c
      real(8) :: v,dv
      integer :: iprint,k
      logical :: JASON
c
c----------------------------------------------------------------------
c
      iprint=0 ; JASON=.TRUE.
c
      ShapH(1:MAXP-1)=0.d0 ; DshapH(1:MAXP-1)=0.d0
c
      if (Nord.lt.2) then
        Nrdof=0 ; return
      endif
c
      if (Nord.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format(' shapHbe: Nord,MAXP = ',2i2)
        stop
      endif
c
c  ...JASON implementation
      if (JASON) then
        call shapHbe_Jason(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
c  ...SOLIN implementation
      else
        call shapHke(      T,Nord,Norient, Nrdof,ShapH,DshapH)
        v=(1.d0-T)*T ; dv=1.d0-2.d0*T
c
        do k=1,Nrdof
          DshapH(k) = dv*ShapH(k)+v*DshapH(k)
          ShapH( k) =  v*ShapH(k)
        enddo
      endif
c
      if (iprint.eq.1) then
        write(*,7002) T,Nord,Norient
 7002   format(' shapHbe: T,Nord,Norient = ',f8.3,2x,i2,2x,i1,
     .         ' Shap,Dvshap = ')
        do k=1,Nrdof
          write(*,7003) k,ShapH(k),DshapH(k)
 7003     format('         k = ',i2,2f8.3)
        enddo
        call pause
      endif
c
c
      endsubroutine shapHbe_0
c
c
c
c----------------------------------------------------------------------
      subroutine shapHbe_Jason(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
      use parameters , only : MAXP
      implicit none
c
      real(8),                  intent(in ) :: T
      integer,                  intent(in ) :: Nord
      integer,                  intent(in ) :: Norient
      integer,                  intent(out) :: Nrdof
      real(8),dimension(MAXP-1),intent(out) :: ShapH,DshapH
c
      real(8),dimension(MAXP+1) :: v,d
      integer :: nvoid,i
      real(8) :: s,dsdt
      real(8),parameter :: factor=0.5d0
c----------------------------------------------------------------------
c
c  ...number of dof's associated to edge
      Nrdof=Nord-1
c
c  ...account for orientation
      select case(Norient)
      case(0) ; s=T      ; dsdt=1.d0
      case(1) ; s=1.d0-T ; dsdt=-1.d0
      endselect
c
      call shape1_Jason(s,Nord, nvoid,v,d)
      do i=1,Nrdof
        ShapH(i)=factor*v(i+2) ; DshapH(i)=factor*d(i+2)*dsdt
      enddo
c
c
      end subroutine shapHbe_Jason
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - shapHkt
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate triangle kernel H^1 shape functions
c
c   arguments :
c     in:
c        T             - local triangle coordinates
c        Nord          - triangle middle node order
c        Norient       - triangle (mddle node) orientation
c     out:
c        Nrdof         - number of kernel shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHkt(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
      use control
      use parameters
      use element_data
      implicit none
c
      integer :: Nord,Norient,Nrdof
      real(8) :: T(2),ShapH(MAXmdltH),DshapH(2,MAXmdltH)
c
c  ...global triangle coordinates
      real(8) :: s(2),dsdt(2,2)
c
c  ...work space
      real(8) :: dshapHds(2)
c
      integer :: i,is,j,k,l,lsave
c----------------------------------------------------------------------
c
      integer :: iprint=0
      if (iprint.eq.1) then
        write(*,7005) T(1:2),Nord,Norient
 7005   format('shapHkt: T,Nord,Norient = ',2f8.3,2x,2i3)
      endif
c
      if (Nord.lt.3) then
        Nrdof=0; return
      endif
      if (Nord.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format('shapHkt: Nord,MAXP = ',2i2)
        stop 1
      endif
c
c  ...local to global coordinates
      call local2global(TRIA,T,Norient, s,dsdt)
ccc      select case(Norient)
ccc      case (0); s(1:2) = T(1:2)
ccc      case (1); s(1) = T(2); s(2) = 1.d0 - T(1) - T(2)
ccc      case (2); s(1) = 1.d0 - T(1) - T(2); s(2) = T(1)
ccc      case (3); s(1) = T(2); s(2) = T(1)
ccc      case (4); s(1) = 1.d0 - T(1) - T(2); s(2) = T(2)
ccc      case (5); s(1) = T(1); s(2) = 1.d0 - T(1) - T(2)
ccc      case default
ccc        write(*,7002) Norient
ccc 7002   format('shapHkt: UNKNOWN ORIENTATION', i3)
ccc        stop 1
ccc      end select
c
c  ...evaluate kernel shape functions in global coordinate

      select case(SHAPE_FLAG)
c
c  ...Peano and Legendre shape functions (temporary)
      case(1,2)
        k=1
        ShapH(k)    = 1.d0
        DshapH(1:2,k) = 0.d0
c
c  .....save current index
        lsave=k

c  .....loop through orders higher than 3:
        do i=4,Nord
          l=lsave
c
c  .......first function in row of Pascal triangle:
          k=k+1
          ShapH(k) = ShapH(l)*s(1)
          DshapH(1,k) = DshapH(1,l)*s(1) + ShapH(l)
          DshapH(2,k) = DshapH(2,l)*s(1)
c
c  .......save current index - of the first function in row:
          lsave=k
c
c  .......loop through remaining functions in the row:
          do j=0,i-4
            k=k+1
            ShapH(k) = ShapH(l+j)*s(2)
            DshapH(1,k) = DshapH(1,l+j)*s(2)
            DshapH(2,k) = DshapH(2,l+j)*s(2) + ShapH(l+j)
          enddo
c
c  .....end of loop through degrees of polynomials
        enddo
c
      end select
c
c  ...save the number of dof
      Nrdof=k
c
c  ...account for the orientation
      do k=1,Nrdof
        dshapHds(1:2) = DshapH(1:2,k)
        do is=1,2
          DshapH(is,k) = dshapHds(1)*dsdt(1,is)
     .                 + dshapHds(2)*dsdt(2,is)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,7003)
 7003   format('shapHkt: Shap,Dvshap = ')
        do k=1,Nrdof
          write(*,7004) k,ShapH(k),DshapH(1:2,k)
 7004     format('         k = ',i2,3f8.3)
        enddo
        call pause
      endif
c
c
      end subroutine shapHkt
c
c----------------------------------------------------------------------
c
c   routine name       - shapHbt
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate triangle bubble H^1 shape functions
c
c   arguments :
c     in:
c        T             - local triangle coordinates
c        Nord          - triangle middle node order
c        Norient       - triangle orientation
c     out:
c        Nrdof         - number of shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHbt(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
      use parameters
      implicit none
c
      integer :: Nord,Norient,Nrdof
      real(8) :: T(2),ShapH(MAXmdltH),DshapH(2,MAXmdltH)
c
      integer :: k
c
c  ...bubble
      real(8) :: v,dv(2)
c
      if (Nord.lt.3) then
        Nrdof=0; return
      endif
      if (Nord.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format('shapHbt: Nord,MAXP = ',2i2)
        stop 1
      endif
c
c  ...evaluate kernel shape functions
      call shapHkt(T,Nord,Norient, Nrdof,ShapH,DshapH)
c
c  ...evaluate bubble shape functions
      v = (1.d0 - T(1) - T(2))*T(1)*T(2)
      dv(1) = -T(1)*T(2) + (1.d0 - T(1) - T(2))*T(2)
      dv(2) = -T(1)*T(2) + (1.d0 - T(1) - T(2))*T(1)
      do k=1,Nrdof
        DshapH(1:2,k) = dv(1:2)*ShapH(k) + v*DshapH(1:2,k)
        ShapH(k) = v*ShapH(k)
      enddo
c
c
      end subroutine shapHbt
c
c----------------------------------------------------------------------
c
c   routine name       - shapHb2_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate 2D bubble H^1 shape functions
c
c   arguments :
c     in:
c        Type          - middle node type
c        T             - local triangle coordinates
c        Nord          - element middle node order
c        Norient       - middle node orientation
c     out:
c        Nrdof         - number of shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHb2_0(Type,T,Nord,Norient, Nrdof,ShapH,DshapH)
c
      use element_data , only : NFAXES
      use parameters
      implicit none
c
      integer :: Nord,Norient,Nrdof
      real(8) :: T(2),ShapH(MAXmdlqH),DshapH(2,MAXmdlqH)
c
c  ...1D bubble shape functions
      real(8) :: shapH1(MAXP-1),dshapH1(MAXP-1),
     .           shapH2(MAXP-1),dshapH2(MAXP-1)
c
      character (len=4) :: Type
c
      integer :: i,j,k,nordh,nordv,nrdofH1,nrdofH2
c
c  ...bubble
c      dimension dv(2)
c
      integer :: iprint
c
c----------------------------------------------------------------------
c
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Type,Nord,Norient
 7001   format('shapHb2: Type,Nord,Norient = ',a5,2x,2i3)
        call pause
      endif

      select case(Type)
      case('tria','mdlt')
        call shapHbt(T,Nord,Norient, Nrdof,ShapH,DshapH)
      case('rect','mdlq')
        k=0
        call decode(Nord, nordh,nordv)
        select case(NFAXES(3,Norient))
c
c  .....the axes have NOT been reversed
        case(0)
          call shapHbe(T(1),nordh,NFAXES(1,Norient),
     .                 nrdofH1,shapH1,dshapH1)
          call shapHbe(T(2),nordv,NFAXES(2,Norient),
     .                 nrdofH2,shapH2,dshapH2)
          do j=1,nrdofH2
            do i=1,nrdofH1
              k=k+1
              ShapH(k) = shapH1(i)*shapH2(j)
              DshapH(1,k) = dshapH1(i)*shapH2(j)
              DshapH(2,k) = shapH1(i)*dshapH2(j)
            enddo
          enddo
c
c  .....the axes HAVE been reversed
        case(1)
          call shapHbe(T(1),nordv,NFAXES(1,Norient),
     .                 nrdofH1,shapH1,dshapH1)
          call shapHbe(T(2),nordh,NFAXES(2,Norient),
     .                 nrdofH2,shapH2,dshapH2)
c
c  .......the order of loops gets reversed...
          do i=1,nrdofH1
            do j=1,nrdofH2
              k=k+1
              ShapH(k) = shapH1(i)*shapH2(j)
              DshapH(1,k) = dshapH1(i)*shapH2(j)
              DshapH(2,k) = shapH1(i)*dshapH2(j)
            enddo
          enddo
        end select
        Nrdof=k
      end select
c
      if (iprint.eq.1) then
        write(*,7002) Nrdof
 7002   format('shapHb2: Nrdof = ',i3)
        do k=1,Nrdof
          write(*,7003) k,ShapH(k),DshapH(1:2,k)
 7003     format('k = ',i3,2x,' ShapH = ',e12.5,' DshapH = ',2e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine shapHb2_0
c
c----------------------------------------------------------------------
c
c   routine name       - shapHkn
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate tetrahedroN kernel H^1 shape functions
c
c   arguments :
c     in:
c        T             - master element coordinates
c        Nord          - middle node order
c     out:
c        Nrdof         - number of kernel shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHkn(T,Nord, Nrdof,ShapH,DshapH)
c
      use control
      use parameters
      implicit none
c
      integer :: Nord,Nrdof
      real(8) :: T(3),ShapH(MAXmdlnH),DshapH(3,MAXmdlnH)
c
      integer :: i,j,k,l,lsave
c
      integer :: iprint=1
c
      if (Nord.lt.4) then
        Nrdof=0; return
      endif
      if (Nord.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format('shapHkn: Nord,MAXP = ',2i2)
        stop 1
      endif
c
      select case(SHAPE_FLAG)
c
c  ...Peano and Legendre shape functions (temporary)
      case(1,2)
        k=1
        ShapH(k)    = 1.d0
        DshapH(1:3,k) = 0.d0
c
c  .....save current index
        lsave=k

c  .....loop through orders higher than 4:
        do i=5,Nord
          l=lsave
c
c  .......first function in bottom row
          k=k+1
          ShapH(k) = ShapH(l)*T(1)
          DshapH(1,k) = DshapH(1,l)*T(1) + ShapH(l)
          DshapH(2,k) = DshapH(2,l)*T(1)
          DshapH(3,k) = DshapH(3,l)*T(1)
c
c  .......save current index - of the first function in bottom row:
          lsave=k
c
c  .......loop through remaining functions in the bottom row:
          do j=0,i-5
            k=k+1
            ShapH(k) = ShapH(l+j)*T(2)
            DshapH(1,k) = DshapH(1,l+j)*T(2)
            DshapH(2,k) = DshapH(2,l+j)*T(2) + ShapH(l+j)
            DshapH(3,k) = DshapH(3,l+j)*T(2)
          enddo
c
c  .......loop through remaining functions (upper triangle)
          do j=0,(i-4)*(i-3)/2-1
            k=k+1
            ShapH(k) = ShapH(l+j)*T(3)
            DshapH(1,k) = DshapH(1,l+j)*T(3)
            DshapH(2,k) = DshapH(2,l+j)*T(3)
            DshapH(3,k) = DshapH(3,l+j)*T(3) + ShapH(l+j)
          enddo
c
c  .....end of loop through degrees of polynomials
        enddo
c
      end select
c
c  ...save the number of dof
      Nrdof=k
c
c
      end subroutine shapHkn
c
c----------------------------------------------------------------------
c
c   routine name       - shapHbn
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - evaluate tetrahedroN bubble H^1 shape functions
c
c   arguments :
c     in:
c        T             - master element coordinates
c        Nord          - middle node order
c     out:
c        Nrdof         - number of kernel shape functions
c        ShapH         - values of shape functions
c        DshapH        - values of derivatives of the shape functions
c                        wrt to master element coordinates
c
c----------------------------------------------------------------------
c
      subroutine shapHbn(T,Nord, Nrdof,ShapH,DshapH)
c
      use control
      use parameters
      implicit none
c
      integer :: Nord,Nrdof
      real(8) :: T(3),ShapH(MAXmdlnH),DshapH(3,MAXmdlnH)
c
c  ...bubble
      real(8) :: v,dv(3)
c
      integer :: k
c
      integer :: iprint=1
c
      if (Nord.lt.4) then
        Nrdof=0; return
      endif
      if (Nord.gt.MAXP) then
        write(*,7001) Nord,MAXP
 7001   format('shapHkn: Nord,MAXP = ',2i2)
        stop 1
      endif
c
c  ...determine kernel shape functions
      call shapHkn(T,Nord, Nrdof,ShapH,DshapH)
c
c  ...evaluate bubble shape functions
      v     = (1.d0-T(1)-T(2)-T(3))*T(1)*T(2)*T(3)
      dv(1) = (1.d0-2.d0*T(1)-T(2)-T(3))*T(2)*T(3)
      dv(2) = (1.d0-T(1)-2.d0*T(2)-T(3))*T(1)*T(3)
      dv(3) = (1.d0-T(1)-T(2)-2.d0*T(3))*T(1)*T(2)
      do k=1,Nrdof
        DshapH(1:3,k) = dv(1:3)*ShapH(k) + v*DshapH(1:3,k)
        ShapH(k) = v*ShapH(k)
      enddo
c
c
      end subroutine shapHbn
