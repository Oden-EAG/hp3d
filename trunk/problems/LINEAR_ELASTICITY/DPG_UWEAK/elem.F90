!--------------------------------------------------------------------------
!> Purpose : return stiffness matrix and load vector for element
!!
!! @param[in]  Mdle   - an element (middle node) number
!! @param[out] Nrdof  - number of dof for a single component
!! @param[out] Itest  - index for assembly
!! @param[out] Itrial - index for assembly
!--------------------------------------------------------------------------
!
subroutine elem(Mdle, Itest,Itrial)
  ! use parameters, only : ZERO
  use physics   , only : NR_PHYSA
  ! use assembly  , only : ALOC,BLOC,NR_RHS
  use data_structure3D
!--------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!--------------------------------------------------------------------------

  Itest (1:NR_PHYSA) = 0
  Itrial(1:NR_PHYSA) = 0

  ! ALOC(1,1)%array = ZERO; ALOC(1,2)%array = ZERO; ALOC(1,3)%array = ZERO; ALOC(1,4)%array = ZERO; ALOC(1,5)%array = ZERO; BLOC(1)%array = ZERO
  ! ALOC(2,1)%array = ZERO; ALOC(2,2)%array = ZERO; ALOC(2,3)%array = ZERO; ALOC(2,4)%array = ZERO; ALOC(2,5)%array = ZERO; BLOC(2)%array = ZERO
  ! ALOC(3,1)%array = ZERO; ALOC(3,2)%array = ZERO; ALOC(3,3)%array = ZERO; ALOC(3,4)%array = ZERO; ALOC(3,5)%array = ZERO; BLOC(3)%array = ZERO
  ! ALOC(4,1)%array = ZERO; ALOC(4,2)%array = ZERO; ALOC(4,3)%array = ZERO; ALOC(4,4)%array = ZERO; ALOC(4,5)%array = ZERO; BLOC(4)%array = ZERO
  ! ALOC(5,1)%array = ZERO; ALOC(5,2)%array = ZERO; ALOC(5,3)%array = ZERO; ALOC(5,4)%array = ZERO; ALOC(5,5)%array = ZERO; BLOC(5)%array = ZERO

  select case(NODES(Mdle)%case)
  !  we wish to support both the displacement field variables and the Cauchy stress variables at each point
  !  (all physical attributes are supported when NODES(Mdle)%case == 2**NR_PHYSA-1)
  case(31)

    Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
    call elem_DPG_UWEAK(Mdle)
    ! call elem_DPG_UWEAK(  &
    !             Mdle,BLOC(1)%nrow,BLOC(2)%nrow,BLOC(3)%nrow,BLOC(4)%nrow,BLOC(5)%nrow  &
    !             ALOC(1,1)%array,ALOC(1,2)%array,ALOC(1,3)%array,ALOC(1,4)%array,ALOC(1,5)%array,BLOC(1)%array,  &
    !             ALOC(2,1)%array,ALOC(2,2)%array,ALOC(2,3)%array,ALOC(2,4)%array,ALOC(2,5)%array,BLOC(2)%array,  &
    !             ALOC(3,1)%array,ALOC(3,2)%array,ALOC(3,3)%array,ALOC(3,4)%array,ALOC(3,5)%array,BLOC(3)%array,  &
    !             ALOC(4,1)%array,ALOC(4,2)%array,ALOC(4,3)%array,ALOC(4,4)%array,ALOC(4,5)%array,BLOC(4)%array,  &
    !             ALOC(5,1)%array,ALOC(5,2)%array,ALOC(5,3)%array,ALOC(5,4)%array,ALOC(5,5)%array,BLOC(5)%array)
  case default
    write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
               Mdle,NODES(Mdle)%case
    call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
  end select
!
end subroutine elem

!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for primal DPG elasticity problem
!! @param[in]  Mdle      - middle node number
!! @param[in]  Nrow_B1   - number of rows of component 1  of load vector
!! @param[in]  Nrow_B2   - number of rows of component 2 of load vector
!! @param[in]  Nrow_B3   - number of rows of component 3 of load vector
!! @param[in]  Nrow_B4   - number of rows of component 4 of load vector
!! @param[in]  Nrow_B5   - number of rows of component 5 of load vector
!!
!! @param[out] Aloc11,Aloc12,Aloc13,Aloc14,Aloc15,
!!             Aloc21,Aloc22,Aloc23,Aloc24,Aloc25,
!!             Aloc31,Aloc32,Aloc33,Aloc34,Aloc35,   - elem stiffness matrix
!!             Aloc41,Aloc42,Aloc43,Aloc44,Aloc45,
!!             Aloc51,Aloc52,Aloc53,Aloc54,Aloc55
!! @param[out] Bloc1,Bloc2,Bloc3,Bloc4,Bloc5         - elem load vector(s)
!------------------------------------------------------------------------------------------
!
!                |   \tau \in H(div)^3   |       v \in (H1)^3     |
!
!                   - <\hat u,(\tau n)>  +            0
!  +                        0            + - <\hat (\sigma n),v>           
!  + \int_\Omega [   u \cdot div(\tau)   +            0           ]
!  + \int_\Omega [    A \sigma : \tau    +    \sigma : grad(v)    ]
!  + \int_\Omega [     \omega : \tau     +            0           ]
!  = \int_\Omega [          0            +        f \cdot v       ]
!
!------------------------------------------------------------------------------------------
!
subroutine elem_DPG_UWEAK(Mdle)
! subroutine elem_DPG_UWEAK(Mdle,Nrow_Bloc1,Nrow_Bloc2,Nrow_Bloc3,Nrow_Bloc4,Nrow_Bloc5,  &
!                           Aloc11,Aloc12,Aloc13,Aloc14,Aloc15,Bloc1,  &
!                           Aloc21,Aloc22,Aloc23,Aloc24,Aloc25,Bloc2,  &
!                           Aloc31,Aloc32,Aloc33,Aloc34,Aloc35,Bloc3,  &
!                           Aloc41,Aloc42,Aloc43,Aloc44,Aloc45,Bloc4,  &
!                           Aloc51,Aloc52,Aloc53,Aloc54,Aloc55,Bloc5)
      ! use uweak_module
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use isotropic_elast_material
      use physics   , only : NR_PHYSA      
      use assembly  , only : ALOC,BLOC,NR_RHS!, super_array
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
!------------------------------------------------------------------------------------------
      implicit none
      integer,                                  intent(in)  :: Mdle

!  ...Necessary constant since there is no dynamic allocation
      integer, parameter :: MAXNRHS_MOD = 1
!
!  ...MATRICES
!     stiffnes and load matrices for the enriched test space
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickH) :: EnrTraceDispl!,EnrTraceDisplc
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickV) :: EnrTraceStress!,EnrTraceStressc
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickQ) :: EnrFieldDispl!,EnrFieldDisplc
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,6*MAXbrickQ) :: EnrFieldStress!,EnrFieldStressc
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*MAXbrickQ) :: EnrFieldOmega!,EnrFieldOmegac
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,MAXNRHS_MOD) :: EnrLoad
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+12*MAXbrickQ) :: EnrStiffness
      real*8, dimension(3*MAXbrickVV+3*MAXbrickHH,3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+12*MAXbrickQ+MAXNRHS_MOD) :: EnrEverything
      ! real*8, dimension(3*MAXbrickH+3*MAXbrickV+12*MAXbrickQ,3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+12*MAXbrickQ+MAXNRHS_MOD) :: FullDPG
!     Gram matrix for the local Riesz matrix in LAPACK format
      real*8, dimension((3*MAXbrickVV+3*MAXbrickHH)*(3*MAXbrickVV+3*MAXbrickHH+1)/2) :: Gram

      real*8, allocatable :: FullDPG(:,:)

      ! integer,                                  intent(in)  :: Nrow_Bloc1
      ! integer,                                  intent(in)  :: Nrow_Bloc2
      ! integer,                                  intent(in)  :: Nrow_Bloc3
      ! integer,                                  intent(in)  :: Nrow_Bloc4
      ! integer,                                  intent(in)  :: Nrow_Bloc5
      ! real*8, dimension(Nrow_Bloc1,Nrow_Bloc1), intent(out) :: Aloc11
      ! real*8, dimension(Nrow_Bloc1,Nrow_Bloc2), intent(out) :: Aloc12
      ! real*8, dimension(Nrow_Bloc1,Nrow_Bloc3), intent(out) :: Aloc13
      ! real*8, dimension(Nrow_Bloc1,Nrow_Bloc4), intent(out) :: Aloc14
      ! real*8, dimension(Nrow_Bloc1,Nrow_Bloc5), intent(out) :: Aloc15
      ! real*8, dimension(Nrow_Bloc2,Nrow_Bloc1), intent(out) :: Aloc21
      ! real*8, dimension(Nrow_Bloc2,Nrow_Bloc2), intent(out) :: Aloc22
      ! real*8, dimension(Nrow_Bloc2,Nrow_Bloc3), intent(out) :: Aloc23
      ! real*8, dimension(Nrow_Bloc2,Nrow_Bloc4), intent(out) :: Aloc24
      ! real*8, dimension(Nrow_Bloc2,Nrow_Bloc5), intent(out) :: Aloc25
      ! real*8, dimension(Nrow_Bloc3,Nrow_Bloc1), intent(out) :: Aloc31
      ! real*8, dimension(Nrow_Bloc3,Nrow_Bloc2), intent(out) :: Aloc32
      ! real*8, dimension(Nrow_Bloc3,Nrow_Bloc3), intent(out) :: Aloc33
      ! real*8, dimension(Nrow_Bloc3,Nrow_Bloc4), intent(out) :: Aloc34
      ! real*8, dimension(Nrow_Bloc3,Nrow_Bloc5), intent(out) :: Aloc35
      ! real*8, dimension(Nrow_Bloc4,Nrow_Bloc1), intent(out) :: Aloc41
      ! real*8, dimension(Nrow_Bloc4,Nrow_Bloc2), intent(out) :: Aloc42
      ! real*8, dimension(Nrow_Bloc4,Nrow_Bloc3), intent(out) :: Aloc43
      ! real*8, dimension(Nrow_Bloc4,Nrow_Bloc4), intent(out) :: Aloc44
      ! real*8, dimension(Nrow_Bloc4,Nrow_Bloc5), intent(out) :: Aloc45
      ! real*8, dimension(Nrow_Bloc5,Nrow_Bloc1), intent(out) :: Aloc51
      ! real*8, dimension(Nrow_Bloc5,Nrow_Bloc2), intent(out) :: Aloc52
      ! real*8, dimension(Nrow_Bloc5,Nrow_Bloc3), intent(out) :: Aloc53
      ! real*8, dimension(Nrow_Bloc5,Nrow_Bloc4), intent(out) :: Aloc54
      ! real*8, dimension(Nrow_Bloc5,Nrow_Bloc5), intent(out) :: Aloc55
      ! real*8, dimension(Nrow_Bloc1,NR_RHS),     intent(out) :: Bloc1
      ! real*8, dimension(Nrow_Bloc2,NR_RHS),     intent(out) :: Bloc2
      ! real*8, dimension(Nrow_Bloc3,NR_RHS),     intent(out) :: Bloc3
      ! real*8, dimension(Nrow_Bloc4,NR_RHS),     intent(out) :: Bloc4
      ! real*8, dimension(Nrow_Bloc5,NR_RHS),     intent(out) :: Bloc5
!------------------------------------------------------------------------------------------
!
!  ...element and face type
      character(len=4) :: etype,ftype
!
!  ...number of topological entities (vertices,edges,faces)
      integer :: nrv,nre,nrf
!
!  ...element and face order, enriched order
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
      integer                :: nordP
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...SHAPE FUNCTIONS
!     H1  (geometry and trial)
      real*8, dimension(  MAXbrickH)  :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
      integer                         :: nrdofH
!     H(div)  (trial)
      real*8, dimension(3,MAXbrickV)  :: shapV
      real*8, dimension(  MAXbrickV)  :: divV
      real*8, dimension(  MAXbrickV)  :: shapV_n
      integer                         :: nrdofV
!     L2  (trial)
      real*8, dimension(  MAXbrickQ)  :: shapQ
      integer                         :: nrdofQ
!     H1   (test)
      real*8, dimension(  MAXbrickHH) :: shapHH
      real*8, dimension(3,MAXbrickHH) :: gradHH
      integer                         :: nrdofHH
!     H(div)  (test)
      real*8, dimension(3,MAXbrickVV) :: shapVV
      real*8, dimension(  MAXbrickVV) :: divVV
      real*8, dimension(  MAXbrickVV) :: shapVV_n
      integer                         :: nrdofVV
!
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt,rt
      integer                        :: nsign
!
!  ...tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: A, AA
!
!  ...source term (don't need Neumann term)
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, dimension(3,MAXNINT3ADD) :: xiloc
      real*8, dimension(MAXNINT3ADD)   :: wxi
!
!  ...2D quadrature data for boundary terms
      real*8, dimension(2,MAXNINT2ADD) :: tloc
      real*8, dimension(MAXNINT2ADD)   :: wt
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,k1,k2,k3,k4,k5,m1,m2,m3,m4,m5,n1,n2,n3,ipt,ifc,  &
                 icomp,jcomp,kcomp,lcomp,nint,iprint,iflag,info,info1,info2,info3,  &
                 kH,kV,kQ,lH,lV,lQ,kmin,kmax,enrdof,                                      &
                 nrHbub,nrEbub,nrVbub,nrQbub,       gdump,bdump6,bdump7
      integer, dimension(NR_PHYSA) :: ndofphysics
      real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax,alpha
!
!  ...LAPACK stuff
      character uplo,transa,transb
! NOTE: nk is a "statement function"
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
      iprint=0

      alpha = 1.d0
!
!  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...set the enriched order of appoximation
      select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
      end select
!
! !  ...initialize the local element matrices and load vectors
!       Aloc11 = ZERO; Aloc12 = ZERO; Aloc13 = ZERO; Aloc14 = ZERO; Aloc15 = ZERO; Bloc1 = ZERO
!       Aloc21 = ZERO; Aloc22 = ZERO; Aloc23 = ZERO; Aloc24 = ZERO; Aloc25 = ZERO; Bloc2 = ZERO
!       Aloc31 = ZERO; Aloc32 = ZERO; Aloc33 = ZERO; Aloc34 = ZERO; Aloc35 = ZERO; Bloc3 = ZERO
!       Aloc41 = ZERO; Aloc42 = ZERO; Aloc43 = ZERO; Aloc44 = ZERO; Aloc45 = ZERO; Bloc4 = ZERO
!       Aloc51 = ZERO; Aloc52 = ZERO; Aloc53 = ZERO; Aloc54 = ZERO; Aloc55 = ZERO; Bloc5 = ZERO
!  ...initialize the enriched local element stiffness matrices and load vectors
      EnrTraceDispl=ZERO; EnrTraceStress=ZERO; EnrFieldDispl=ZERO; EnrFieldStress=ZERO; EnrFieldOmega=ZERO
      EnrLoad=ZERO
! !  ...initialize the Gram matrix
      Gram=ZERO
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = NORD_ADD
      call set_3dint_DPG(etype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....Compute shape functions needed for test/trial field variables and geometry
!       H1 (geometry)
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!       L2 (field trial)
        call shape3DQ(etype,xi,norder, nrdofQ,shapQ)
!       H1 (test)
        call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!       H(div) (test)
        call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
                    x,dxdxi,dxidx,rjac,iflag)
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif
!
!  .....Change coordinates so the shape functions are on the physical element
!       L2 (trial)
        shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac
!       H1 (test)
        do k=1,nrdofHH
          gradHH(1:3,k) = gradHH(1,k)*dxidx(1,1:3)  &
                        + gradHH(2,k)*dxidx(2,1:3)  &
                        + gradHH(3,k)*dxidx(3,1:3)
        enddo
!       H(div) (test)
        do k=1,nrdofVV
          shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
                        + dxdxi(1:3,2)*shapVV(2,k)  &
                        + dxdxi(1:3,3)*shapVV(3,k)
        enddo
        shapVV(1:3,1:nrdofVV) = shapVV(1:3,1:nrdofVV)/rjac
        divVV(1:nrdofVV) = divVV(1:nrdofVV)/rjac
!
!  .....integration weight
        weight = wa*rjac
!
!  .....compute the compliance tensor
        call getA(x, A)
!
!  .....need this for the adjoint graph norm
        call getAA(X, AA)
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!  NOTE: Because nrddofVV \neq nrdofQ, there is no convenient way of organizing
!        the matrices testing all of the n-th basis functions in turn before
!        moving on to the (n+1)-th basis functions. Therefore, and which will lead
!        to a sparser system as well, we test all the basis functions representing
!        \tau and then all those representing v.
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!  .....FIRST OUTER loop through enriched H(div) dofs
        do k1=1,nrdofVV
!  .......OUTER loop through components
          do jcomp=1,3
            m1 = (k1-1)*3+jcomp
!
!           E N R I C H E D   L O A D   V E C T O R
!
!   0
!
!           G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (A:\tau_2,A:\tau)+(div(\tau_2),div(tau))+(\tau_2,\tau)+alpha*(\tau_2,\tau)
!
            case(1)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-1)/3)+1
                do kcomp=1,3
                  do lcomp=1,3
                    if ((icomp.eq.jcomp.and.kcomp.eq.lcomp).or. &
                        (icomp.eq.kcomp.and.jcomp.eq.lcomp).or. &
                        (icomp.eq.lcomp.and.jcomp.eq.kcomp))    &
                      Gram(k) = Gram(k)  &
                              + ( AA(icomp,kcomp,jcomp,lcomp)                 &
                                 *shapVV(lcomp,k1)*shapVV(kcomp,k2) )*weight
                  enddo
                enddo
                if (icomp.eq.jcomp) then
                  Gram(k) = Gram(k)  + divVV(k1)*divVV(k2)*weight
                  Gram(k) = Gram(k)  &
                          + (1.d0+alpha)*(               &
                            + shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*weight
                endif
              enddo
!
!   (A:\tau_2, \grad v)
!
              do m2=3*nrdofVV+1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                k = nk(m1,m2)
                k2 = int((m2-3*nrdofVV-1)/3)+1
                do kcomp=1,3
                  do lcomp=1,3
                      if ((icomp.eq.jcomp.and.kcomp.eq.lcomp).or. &
                        (icomp.eq.kcomp.and.jcomp.eq.lcomp).or. &
                        (icomp.eq.lcomp.and.jcomp.eq.kcomp))    &
                        Gram(k) = Gram(k)  &
                                + ( A(icomp,kcomp,jcomp,lcomp)                 &
                                   *shapVV(lcomp,k1)*gradHH(kcomp,k2) )*weight
                  enddo
                enddo
              enddo
!
!   MATHEMATICIAN'S
!   (\tau_2,\tau)+(div(\tau_2),div(tau))
!
            case(2)
              do m2=m1,3*nrdofVV
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( divVV(k1)*divVV(k2)  &
                            + shapVV(1,k1)*shapVV(1,k2)  &
                            + shapVV(2,k1)*shapVV(2,k2)  &
                            + shapVV(3,k1)*shapVV(3,k2) )*weight
                endif
              enddo
            end select
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   A \sigma : \tau
!
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
!
!  .........INNER loop through trial dofs for Cauchy stress
            do k3=1,nrdofQ
              do icomp=1,6
                m3 = (k3-1)*6+icomp
                if (icomp.eq.1) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(1,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.2) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(2,2,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.3) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(3,3,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.4) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(1,2,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                        + A(2,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.5) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(1,3,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                        + A(3,1,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                elseif (icomp.eq.6) then
                  tmp = 0.d0
                  do n=1,3
                    tmp = tmp  &
                        + A(2,3,jcomp,n)*shapQ(k3)*shapVV(n,k1) &
                        + A(3,2,jcomp,n)*shapQ(k3)*shapVV(n,k1)
                  enddo
                endif
                EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                      + tmp*weight
              enddo
            enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   u \cdot div(\tau)
!
!  .........INNER loop through trial dofs for displacement
            do k4=1,nrdofQ
              do icomp=1,3
                m4 = (k4-1)*3+icomp
                if (icomp.eq.jcomp) then
                  EnrFieldDispl(m1,m4) = EnrFieldDispl(m1,m4)  &
                                       + shapQ(k4)*divVV(k1)*weight
                endif
              enddo
            enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                   (   L A G R A N G E    M U L T I P L I E R   )
!
!   + omega : \tau
!
!          (   0     omega3  -omega2 )    1
!  omega = (-omega3    0      omega1 ) = --- ( grad(u)-grad(u)^T)
!          ( omega2  -omega1    0    )    2
!
!  .........INNER loop through trial dofs for omega (antisymmetric part of displacement gradient)
            do k5=1,nrdofQ
              do icomp=1,3
                m5 = (k5-1)*3+icomp
                if (icomp.eq.1) then
                  if (jcomp.eq.2) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(2,k1)*weight
                  endif
                elseif (icomp.eq.2) then
                  if (jcomp.eq.1) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(1,k1)*weight
                  endif
                elseif (icomp.eq.3) then
                  if (jcomp.eq.1) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         + shapQ(k5)*shapVV(2,k1)*weight
                  elseif (jcomp.eq.2) then
                    EnrFieldOmega(m1,m5) = EnrFieldOmega(m1,m5)  &
                                         - shapQ(k5)*shapVV(1,k1)*weight
                  endif
                endif
              enddo
            enddo
!  .....END OUTER LOOP through test stresses
          enddo
        enddo
!
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows) 
!
!  .....SECOND OUTER loop through enriched H1 dofs
        do k1=1,nrdofHH
!  .......OUTER loop through components
          do jcomp=1,3
            ! counter of row
            m1 = 3*nrdofVV+(k1-1)*3+jcomp
!
!
!            E N R I C H E D   L O A D   V E C T O R
!
!   f \cdot v
!
            EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                                 + fval(jcomp,1:NR_RHS)*shapHH(k1)*weight
!
!            G R A M   M A T R I X
!
            select case(TEST_NORM)
!
!   ADJOINT GRAPH
!   (\grad(v_1),A:\tau_2+\grad(v_2))+alpha*(v_1,v_2)
!
            case(1)

              do m2=m1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( alpha*shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
!
!   MATHEMATICIAN'S
!   (v_2,v) + (grad(v_2),grad(v))
!
            case(2)
              do m2=m1,3*nrdofVV+3*nrdofHH
                icomp = mod(m2-1,3)+1
                if (icomp.eq.jcomp) then
                  k = nk(m1,m2)
                  k2 = int((m2-3*nrdofVV-1)/3)+1
                  Gram(k) = Gram(k)  &
                          + ( shapHH(k1)*shapHH(k2)  &
                            + gradHH(1,k1)*gradHH(1,k2)  &
                            + gradHH(2,k1)*gradHH(2,k2)  &
                            + gradHH(3,k1)*gradHH(3,k2) )*weight
                endif
              enddo
            end select
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   + \sigma : grad(v)
!
!          ( sigma1   sigma4  sigma5 )
!  sigma = ( sigma4   sigma2  sigma6 )
!          ( sigma5   sigma6  sigma3 )
!
!  .........INNER loop through trial dofs for Cauchy stress
            do k3=1,nrdofQ
              do icomp=1,6
                m3 = (k3-1)*6+icomp
                if (icomp.eq.1) then
                  if (jcomp.eq.1) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(1,k1)*weight
                  endif
                elseif (icomp.eq.2) then
                  if (jcomp.eq.2) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(2,k1)*weight
                  endif
                elseif (icomp.eq.3) then
                  if (jcomp.eq.3) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(3,k1)*weight
                  endif
                elseif (icomp.eq.4) then
                  if (jcomp.eq.1) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(2,k1)*weight
                  elseif (jcomp.eq.2) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(1,k1)*weight
                  endif
                elseif (icomp.eq.5) then
                  if (jcomp.eq.1) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(1,k1)*weight
                  endif
                elseif (icomp.eq.6) then
                  if (jcomp.eq.2) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(3,k1)*weight
                  elseif (jcomp.eq.3) then
                    EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
                                          + shapQ(k3)*gradHH(2,k1)*weight
                  endif
                endif
              enddo
            enddo
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   0
!
!
!           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
!                   (   L A G R A N G E    M U L T I P L I E R   )
!
!   0
!
!
!  .......END OUTER LOOP
          enddo
        enddo
!
!  ...end of loop through integration points
      enddo
!
!-----------------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L                                           |
!-----------------------------------------------------------------------------------
!
!  ...get number of H1 and Hdiv bubbles
    call find_ndof(mdle,nrHbub,nrEbub,nrVbub,nrQbub)
!  ...loop through element faces
      do ifc=1,nrf
!
!  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,ifc)
!
!  .....face type
        ftype = face_type(etype,ifc)
!
!  .....face order of approximation
        call face_order(etype,ifc,norder, nordf)
!
!  .....set up the face quadrature
        INTEGRATION = NORD_ADD
        call set_2dint_DPG(ftype,nordf, nint,tloc,wt)
        INTEGRATION = 0
!
!  .....loop through face integration points
        do ipt=1,nint
!
!  .......face coordinates
          t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
!
!  .......master element coordinates using face parameterization
          call face_param(etype,ifc,t, xi,dxidt)
!
!  .......Compute shape functions needed for test/trial field variables and geometry
!         H1 (geometry and bdry trial)
          call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                       nrdofH,shapH,gradH)
!         H(div) (bdry trial)
          call shape3DV(etype,xi,norder,nface_orient,  &
                       nrdofV,shapV,divV)
!         H1 (test)
          call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!         H(div) (test)
          call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!  .......geometry map
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
                       x,dxdxi,dxidx,rjac,dxdt,rn,brjac)
          weight = wa*brjac
!
!  .......Change coordinates (only of what is needed in this case) so the shape functions are on the physical element 
!         H(div) (trial)
          do k=1,nrdofV-nrVbub
            shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                         + dxdxi(1:3,2)*shapV(2,k)  &
                         + dxdxi(1:3,3)*shapV(3,k)
            shapV_n(k)  = shapV(1,k)*rn(1)  &
                        + shapV(2,k)*rn(2)  &
                        + shapV(3,k)*rn(3)
          enddo
          shapV_n(1:nrdofV-nrVbub) = shapV_n(1:nrdofV-nrVbub)/rjac
!
!         H(div) (test)
          do k=1,nrdofVV
            shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
                          + dxdxi(1:3,2)*shapVV(2,k)  &
                          + dxdxi(1:3,3)*shapVV(3,k)
            shapVV_n(k) = shapVV(1,k)*rn(1)  &
                        + shapVV(2,k)*rn(2)  &
                        + shapVV(3,k)*rn(3)
          enddo
          shapVV_n(1:nrdofVV) = shapVV_n(1:nrdofVV)/rjac
!
!
!    P A R T  1 : go through \tau\in H(div)^3 test space (this fills the first set of rows)
!
!
!  .......OUTER loop through enriched H(div) test functions
          do k1=1,nrdofVV
            do jcomp=1,3
              m1 = (k1-1)*3+jcomp
!
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   - <\hat u,(\tau n)>
!
!  ...........INNER loop through H1 bdry trial dofs
              do k2=1,nrdofH-nrHbub
                do icomp=1,3
                  m2 = (k2-1)*3+icomp
                  if (icomp.eq.jcomp) then
                    EnrTraceDispl(m1,m2) = EnrTraceDispl(m1,m2)  &
                                         - shapH(k2)*shapVV_n(k1)*weight
                  endif
                enddo
              enddo
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   C A U C H Y    S T R E S S   )
!
!   0
!
!
!  .........END OUTER LOOP
            enddo
          enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space (this fills the second set of rows) 
!
!  .......SECOND OUTER loop through enriched H1 test function
          do k1=1,nrdofHH
!  .........OUTER loop through components
            do jcomp=1,3
              ! counter of row
              m1 = 3*nrdofVV+(k1-1)*3+jcomp
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   0
!
!
!
!           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
!                          (   D I S P L A C E M E N T   )
!
!   - <\hat \sigma,v>
!
!  ...........INNER loop through H(div) bdry trial dofs
              do k3=1,nrdofV-nrVbub
                do icomp=1,3
                  m3 = (k3-1)*3+icomp
                  if (icomp.eq.jcomp) then
                    EnrTraceStress(m1,m3) = EnrTraceStress(m1,m3)  &
                                          - shapV_n(k3)*shapHH(k1)*weight
                  endif
                enddo
              enddo
!
!  .........END OUTER LOOP
            enddo
          enddo          
!  .....end of loop over integration points
        enddo
!
!  ...end of loop over faces
      enddo
!
      if (iprint.eq.2) then
        write(*,*) 'Gram = '
        do i=1,25
          write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6000     format('i = ',i3,'  ',25e12.5)
        enddo

        call pause

        write(*,*) 'EnrFieldDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldDispl(i,1:3*nrdofQ)
 6001     format('i = ',i4,'  ',19(/,10e12.5))
        enddo

        call pause

        write(*,*) 'EnrFieldStress = '
        do i=1,3*nrdofVV
          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
        enddo

        call pause

        do i=1+3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrFieldStress(i,1:6*nrdofQ)
        enddo

        call pause

        write(*,*) 'EnrFieldOmega = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrFieldOmega(i,1:3*nrdofQ)
        enddo

        call pause

        write(*,*) 'EnrTraceDispl = '
        do i=1,3*nrdofVV+1
          write(*,6001) i,EnrTraceDispl(i,1:3*(nrdofH-nrHbub))
        enddo

        call pause

        write(*,*) 'EnrTraceStress = '
        do i=3*nrdofVV,3*nrdofVV+3*nrdofHH
          write(*,6001) i,EnrTraceStress(i,1:3*(nrdofV-nrVbub))
        enddo

        call pause

      endif

!
!-----------------------------------------------------------------------------------
!     D P G    L O C A L    A S S E M B L Y                                        |
!-----------------------------------------------------------------------------------
!
!  ...Compact enriched number of rows (total enriched test dof)
      enrdof = 3*nrdofVV+3*nrdofHH
! 
 !            gdump=75
 !        open(unit=gdump,file='output/gram', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! 5999     format(189(e12.5,","))
 !        enddo
 !        close(gdump)
!
!  ...factor the Gram matrix
      uplo = 'U'
      call DPPTRF(uplo, enrdof, Gram, info)
      if (info.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info = ',info ; stop
      endif
! 
 !            gdump=76
 !        open(unit=gdump,file='output/gram_fact', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(gdump,5999)Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 ! ! 5999     format(189(e12.5,","))
 !        enddo
 !        close(gdump)
!
!  ...Create vector of indices with dof of each physics variable (in the same order of physics file)
      ndofphysics = (/3*(nrdofH-nrHbub),3*(nrdofV-nrVbub),3*nrdofQ,6*nrdofQ,3*nrdofQ/)
!
!  ...Construct EnrEverything by packing all Enriched Matrices in the order H1,Hcurl,Hdiv,L2 
!     (in the same order of physics file) and load
      EnrEverything=ZERO
!  ...EnrTraceDispl (H1)
      kmin=0
      kmax=kmin+ndofphysics(1)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrTraceStress (H(div))
      kmin=kmax
      kmax=kmin+ndofphysics(2)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrTraceStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldDispl (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(3)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldDispl(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldStress (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(4)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldStress(1:enrdof,1:(kmax-kmin))
!  ...EnrFieldOmega (L2)
      kmin=kmax
      kmax=kmin+ndofphysics(5)
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrFieldOmega(1:enrdof,1:(kmax-kmin))
!  ...EnrLoad
      kmin=kmax
      kmax=kmin+NR_RHS
      EnrEverything(1:enrdof,kmin+1:kmax)=EnrLoad(1:enrdof,1:(kmax-kmin))
!
      if (iprint.eq.1) then
        write(*,*) 'EnrEverything = '
        do i=1,enrdof
          write(*,6001) i,EnrEverything(i,1:kmax)
        enddo
      endif
 !            bdump7=82
 !      open(unit=bdump7,file='output/EnrEverything', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump7,5995) EnrEverything(i,1:kmax)
 ! ! 5995     format(66(e12.5,","))        
 !        enddo
 !      close(bdump7)
!  ...Save copy of EnrStiffness which implies deleting the load part from EnrEverything
      EnrStiffness=ZERO
      EnrStiffness(1:enrdof,1:kmin)=EnrEverything(1:enrdof,1:kmin)
!
! !  ...save copies of enriched stiffness matrices
!       EnrTraceDisplc=EnrTraceDispl; EnrTraceStressc=EnrTraceStress
!       EnrFieldDisplc=EnrFieldDispl; EnrFieldStressc=EnrFieldStress
!       EnrFieldOmegac=EnrFieldOmega
!
!  ...G^-1 * EnrEverything
      call DPPTRS(uplo,enrdof,kmax,Gram,EnrEverything(:,1:kmax),3*MAXbrickVV+3*MAXbrickHH,info1)
      if (info1.ne.0) then
        write(*,*) 'elem_DPG_UWEAK: info1 = ',info1 ; stop
      endif
! 
 !            bdump7=88
 !      open(unit=bdump7,file='output/EnrEverything_post', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,enrdof
 !          write(bdump7,5995) EnrEverything(i,1:kmax)
 ! ! 5995     format(66(e12.5,","))        
 !        enddo
 !      close(bdump7)
!
!  ...Build full DPG matrix (stiffness + load) in one go
      allocate(FullDPG(kmin,kmax))
      FullDPG = ZERO
      transa = 'T'
      transb = 'N'
      m = kmin
      n = kmax
      k = enrdof
      l = 3*MAXbrickVV+3*MAXbrickHH
      call DGEMM(transa,transb,m,n,k,1.d0,EnrStiffness,l,EnrEverything,l,0.d0,FullDPG,m)!3*(MAXbrickH-MAXmdlbH)+3*(MAXbrickV-MAXmdlbV)+12*MAXbrickQ)
!
!     ULTIMATE DPG LOAD VECTORS AND STIFFNESS MATRICES THROUGH STATIC CONDENSATION
!
!  ...Populate the assembly matrices accordingly
      n1=0
      do i=1,NR_PHYSA
        n=ndofphysics(i)
        n2=n1+n
        m1=0
        do j=1,NR_PHYSA
          m=ndofphysics(j)
          m2=m1+m
!  .......First initialize
          ALOC(i,j)%array = ZERO
!  .......Then populate
          ALOC(i,j)%array(1:n,1:m) = FullDPG(n1+1:n2,m1+1:m2)
          m1=m2  
        enddo
        m2=m1+NR_RHS
!  .....First initialize
        BLOC(i)%array = ZERO
!  .....Then populate
        BLOC(i)%array(1:n,1:NR_RHS) = FullDPG(n1+1:n2,m1+1:m2)
        n1=n2
      enddo
! 
 !        bdump6=81      
 !      open(unit=bdump6,file='output/FullDPG_e1', &
 !          form='formatted',access='sequential',status='unknown')
 !        do i=1,kmin
 !          write(bdump6,5995) FullDPG(i,1:kmax)
 ! 5995     format(55(e22.15,","))        
 !        enddo
 !        close(bdump6)
! 
      if (iprint.eq.1) then
        write(*,*) 'EnrStiffness = '
        do i=1,enrdof
          write(*,6001) i,EnrStiffness(i,1:kmin)
        enddo
        write(*,*) 'Gram^-1*EnrEverything = '
        do i=1,enrdof
          write(*,6001) i,EnrEverything(i,1:kmax)
        enddo
        write(*,*) 'FullDPG = '
        do i=1,kmin
          write(*,6001) i,FullDPG(i,1:kmax)
        enddo
      endif
! 
      deallocate(FullDPG)
! !
!  ...print statments (not complete)
      if (iprint.ge.2) then
        write(*,*)    ''
        write(*,*) 'ndofphysics=',ndofphysics
        write(*,*)    ''
        do i=1,NR_PHYSA
          write(*,*)    'i,BLOC(i)%nrow=',i,BLOC(i)%nrow
          write(*,7011) BLOC(i)%array(1:BLOC(i)%nrow,1)
          write(*,*)    ''
          do j=i,NR_PHYSA
            write(*,*) 'i,ALOC(i,j)%nrow=',i,ALOC(i,j)%nrow
            write(*,*) 'j,ALOC(i,j)%ncol=',j,ALOC(i,j)%ncol
            do k=1,ALOC(i,j)%nrow
              write(*,7013) k,ALOC(i,j)%array(k,1:ALOC(i,j)%ncol)
            enddo
            write(*,*)    ''
          enddo
        enddo
! 
 7011   format(10(/,18e11.4))
 7013   format('i = ',i3,10(/,18e11.4))
! 
      endif
!
!
end subroutine elem_DPG_UWEAK
