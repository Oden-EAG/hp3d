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
  use parameters, only : ZERO
  use physics   , only : NR_PHYSA
  use assembly  , only : ALOC,BLOC,NR_RHS
  use data_structure3D
!--------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!--------------------------------------------------------------------------

  Itest (1:NR_PHYSA) = 0
  Itrial(1:NR_PHYSA) = 0

  !  we wish to support both the displacement field variables and the Cauchy stress variables at each point
  !  (all physical attributes are supported when NODES(Mdle)%case == 2**NR_PHYSA-1)
  if (NODES(Mdle)%case==(2**NR_PHYSA-1)) then
    Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
    call elem_DPG_UWEAK(Mdle)
  else
    write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
               Mdle,NODES(Mdle)%case
    call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
  endif
!
end subroutine elem

!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for primal DPG elasticity problem
!! @param[in]  Mdle      - middle node number (THE most important index of all!!!)
!! @param[in]  Nrow_B1,Nrow_B2,Nrow_B3,Nrow_B4,Nrow_B5   - number of rows of each physics variable (see physics file)
!!
!! @param[out] Aloc11,Aloc12,Aloc13,Aloc14,Aloc15
!!             Aloc21,Aloc22,Aloc23,Aloc24,Aloc25
!!             Aloc31,Aloc32,Aloc33,Aloc34,Aloc35   - elem stiffness matrix
!!             Aloc41,Aloc42,Aloc43,Aloc44,Aloc45
!!             Aloc51,Aloc52,Aloc53,Aloc54,Aloc55
!! @param[out] Bloc1,Bloc2,Bloc3,Bloc4,Bloc5        - elem load vector(s)
!------------------------------------------------------------------------------------------
!
!   9 equations.
!
!                |   \tau \in H(div)^3   |       v \in (L2)^3     |  q \in L2_skew |
!
!    \int_\Omega [    A \sigma : \tau    + - div(\sigma) \cdot v  +   \sigma : q   ]
!  + \int_\Omega [   u \cdot div(\tau)   +            0           +        0       ]
!  + \int_\Omega [      - p : \tau       +            0           +        0       ]
!  +                - <\hat u,(\tau n)>  +            0           +        0
!  = \int_\Omega [          0            +        f \cdot v       +        0       ]
!
!------------------------------------------------------------------------------------------
!
subroutine elem_DPG_UWEAK(Mdle)
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use isotropic_elast_material
      use physics   , only : NR_PHYSA
      use assembly  , only : ALOC,BLOC,NR_RHS
      use common_prob_data, only: SYMMETRY_TOL, TEST_NORM
!------------------------------------------------------------------------------------------
      implicit none
      integer,                                  intent(in)  :: Mdle
!       integer,                                  intent(in)  :: Nrow_Bloc1
!       integer,                                  intent(in)  :: Nrow_Bloc2
!       integer,                                  intent(in)  :: Nrow_Bloc3
!       integer,                                  intent(in)  :: Nrow_Bloc4
!       integer,                                  intent(in)  :: Nrow_Bloc5
!       real*8, dimension(Nrow_Bloc1,Nrow_Bloc1), intent(out) :: Aloc11
!       real*8, dimension(Nrow_Bloc1,Nrow_Bloc2), intent(out) :: Aloc12
!       real*8, dimension(Nrow_Bloc1,Nrow_Bloc3), intent(out) :: Aloc13
!       real*8, dimension(Nrow_Bloc1,Nrow_Bloc4), intent(out) :: Aloc14
!       real*8, dimension(Nrow_Bloc1,Nrow_Bloc5), intent(out) :: Aloc15
!       real*8, dimension(Nrow_Bloc2,Nrow_Bloc1), intent(out) :: Aloc21
!       real*8, dimension(Nrow_Bloc2,Nrow_Bloc2), intent(out) :: Aloc22
!       real*8, dimension(Nrow_Bloc2,Nrow_Bloc3), intent(out) :: Aloc23
!       real*8, dimension(Nrow_Bloc2,Nrow_Bloc4), intent(out) :: Aloc24
!       real*8, dimension(Nrow_Bloc2,Nrow_Bloc5), intent(out) :: Aloc25
!       real*8, dimension(Nrow_Bloc3,Nrow_Bloc1), intent(out) :: Aloc31
!       real*8, dimension(Nrow_Bloc3,Nrow_Bloc2), intent(out) :: Aloc32
!       real*8, dimension(Nrow_Bloc3,Nrow_Bloc3), intent(out) :: Aloc33
!       real*8, dimension(Nrow_Bloc3,Nrow_Bloc4), intent(out) :: Aloc34
!       real*8, dimension(Nrow_Bloc3,Nrow_Bloc5), intent(out) :: Aloc35
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc1), intent(out) :: Aloc41
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc2), intent(out) :: Aloc42
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc3), intent(out) :: Aloc43
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc4), intent(out) :: Aloc44
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc5), intent(out) :: Aloc45
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc1), intent(out) :: Aloc51
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc2), intent(out) :: Aloc52
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc3), intent(out) :: Aloc53
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc4), intent(out) :: Aloc54
!       real*8, dimension(Nrow_Bloc4,Nrow_Bloc5), intent(out) :: Aloc55
!       real*8, dimension(Nrow_Bloc1,NR_RHS),     intent(out) :: Bloc1
!       real*8, dimension(Nrow_Bloc2,NR_RHS),     intent(out) :: Bloc2
!       real*8, dimension(Nrow_Bloc3,NR_RHS),     intent(out) :: Bloc3
!       real*8, dimension(Nrow_Bloc4,NR_RHS),     intent(out) :: Bloc4
!       real*8, dimension(Nrow_Bloc5,NR_RHS),     intent(out) :: Bloc5
! !------------------------------------------------------------------------------------------
! !
! !  ...element and face type
!       character(len=4) :: etype,ftype
! !
! !  ...number of topological entities (vertices,edges,faces)
!       integer :: nrv,nre,nrf
! !
! !  ...element and face order, enriched order
!       integer, dimension(19) :: norder
!       integer, dimension(5)  :: nordf
!       integer                :: nordP
! !
! !  ...node edge and face orientations
!       integer, dimension(12) :: nedge_orient
!       integer, dimension(6)  :: nface_orient
! !
! !  ...SHAPE FUNCTIONS
! !     H1  (geometry and trial)
!       real*8, dimension(  MAXbrickH)  :: shapH
!       real*8, dimension(3,MAXbrickH)  :: gradH
!       integer                         :: nrdofH
! !     H(div)  (field)
!       real*8, dimension(3,MAXbrickV)  :: shapV
!       real*8, dimension(  MAXbrickV)  :: divV
!       integer                         :: nrdofV
! !     L2  (field)
!       real*8, dimension(  MAXbrickQ)  :: shapQ
!       integer                         :: nrdofQ
! !     H(div)  (test)
!       real*8, dimension(3,MAXbrickVV) :: shapVV
!       real*8, dimension(  MAXbrickVV) :: divVV
!       real*8, dimension(  MAXbrickVV) :: shapVV_n
!       integer                         :: nrdofVV
! !     L2  (test)
!       real*8, dimension(  MAXbrickQQ) :: shapQQ
!       integer                         :: nrdofQQ
! !
! !  ...MATRICES
! !     stiffnes matrices (and copies) for the enriched test space
!       real*8, dimension(3*MAXbrickVV+6*MAXbrickQQ,3*MAXbrickH) :: EnrTraceDispl,EnrTraceDisplc
!       real*8, dimension(3*MAXbrickVV+6*MAXbrickQQ,3*MAXbrickV) :: EnrFieldStress,EnrFieldStressc
!       real*8, dimension(3*MAXbrickVV+6*MAXbrickQQ,3*MAXbrickQ) :: EnrFieldDispl,EnrFieldDisplc
!       real*8, dimension(3*MAXbrickVV+6*MAXbrickQQ,3*MAXbrickQ) :: EnrFieldLagra,EnrFieldLagrac
! !     load vector for the enriched space
!       real*8, dimension(3*MAXbrickVV+6*MAXbrickQQ,NR_RHS) :: EnrLoad
! !     Gram matrix for the local Riesz matrix in LAPACK format
!       ! real*8, dimension((3*MAXbrickVV)*(3*MAXbrickVV+1)/2) :: GramStress
!       ! real*8, dimension((3*MAXbrickQQ)*(3*MAXbrickQQ+1)/2) :: GramDispl
!       ! real*8, dimension((3*MAXbrickQQ)*(3*MAXbrickQQ+1)/2) :: GramLagra
!       real*8, dimension((3*MAXbrickVV+6*MAXbrickQQ)*(3*MAXbrickVV+6*MAXbrickQQ+1)/2) :: Gram
! !
! !  ...geometry
!       real*8, dimension(3,MAXbrickH) :: xnod
!       real*8, dimension(3)           :: xi,x,rn
!       real*8, dimension(3,3)         :: dxdxi,dxidx
!       real*8, dimension(2)           :: t
!       real*8, dimension(3,2)         :: dxidt,dxdt,rt
!       integer                        :: nsign
! !
! !  ...tensors in physical coordinates
!       real*8, dimension(3,3,3,3) :: A
! !
! !  ...source term (don't need Neumann term)
!       real*8, dimension(3,MAXNRHS) :: fval
! !
! !  ...3D quadrature data
!       real*8, dimension(3,MAXNINT3ADD) :: xiloc
!       real*8, dimension(MAXNINT3ADD)   :: wxi
! !
! !  ...2D quadrature data for boundary terms
!       real*8, dimension(2,MAXNINT2ADD) :: tloc
!       real*8, dimension(MAXNINT2ADD)   :: wt
! !
! !  ...miscellaneous
!       integer :: i,j,k,l,m,n,k1,k2,k3,k4,k5,m1,m2,m3,m4,m5,n1,n2,n3,ipt,ifc,  &
!                  icomp,jcomp,nint,iprint,iflag,info,info1,info2,info3,  &
!                  kH,kV,kQ,lH,lV,lQ
!       real*8  :: weight,wa,rjac,brjac,tmp,diffmax,dmax
! !
! !  ...LAPACK stuff
!       character uplo
! ! NOTE: nk is a "statement function"
!       integer :: nk
!       nk(k1,k2) = (k2-1)*k2/2+k1
! !
! !-----------------------------------------------------------------------------------
! !      I N I T I A L I Z A T I O N                                                 |
! !-----------------------------------------------------------------------------------
! !
!       iprint=0
! !
! !  ...element type
!       etype = NODES(Mdle)%type
!       nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
! !
! !  ...order of approximation, orientations, geometry dof's (don't need bc flags)
!       call find_order (Mdle, norder)
!       call find_orient(Mdle, nedge_orient,nface_orient)
!       call nodcor     (Mdle, xnod)
! !
! !  ...set the enriched order of appoximation
!       select case(etype)
!       case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
!       case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
!       case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
!       end select
! !
! !  ...initialize the local element matrices and load vectors
! 	  ALOC(1,1)%array = ZERO; ALOC(1,2)%array = ZERO; ALOC(1,3)%array = ZERO; ALOC(1,4)%array = ZERO; ALOC(1,5)%array = ZERO; BLOC(1)%array = ZERO
! 	  ALOC(2,1)%array = ZERO; ALOC(2,2)%array = ZERO; ALOC(2,3)%array = ZERO; ALOC(2,4)%array = ZERO; ALOC(2,5)%array = ZERO; BLOC(2)%array = ZERO
! 	  ALOC(3,1)%array = ZERO; ALOC(3,2)%array = ZERO; ALOC(3,3)%array = ZERO; ALOC(3,4)%array = ZERO; ALOC(3,5)%array = ZERO; BLOC(3)%array = ZERO
! 	  ALOC(4,1)%array = ZERO; ALOC(4,2)%array = ZERO; ALOC(4,3)%array = ZERO; ALOC(4,4)%array = ZERO; ALOC(4,5)%array = ZERO; BLOC(4)%array = ZERO
! 	  ALOC(5,1)%array = ZERO; ALOC(5,2)%array = ZERO; ALOC(5,3)%array = ZERO; ALOC(5,4)%array = ZERO; ALOC(5,5)%array = ZERO; BLOC(5)%array = ZERO
! !  ...initialize the enriched local element stiffness matrices and load vectors
!       EnrTraceDispl=ZERO; EnrFieldStress=ZERO; EnrFieldDispl=ZERO; EnrFieldLagra=ZERO
!       EnrLoad=ZERO
! ! !  ...initialize the Gram matrix
!       Gram=ZERO
! !
! !-----------------------------------------------------------------------------------
! !      E L E M E N T    I N T E G R A L                                            |
! !-----------------------------------------------------------------------------------
! !
! !  ...set up the element quadrature
!       INTEGRATION = NORD_ADD
!       call set_3Dint(etype,norder, nint,xiloc,wxi)
!       INTEGRATION = 0
! !
! !  ...loop through integration points
!       do ipt=1,nint
!         xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
! !
! !  .....Compute shape functions needed for test/trial field variables and geometry
! !       H1 (geometry)
!         call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
!                      nrdofH,shapH,gradH)
! !       H(div) (trial)
!         call shape3V(etype,xi,norder,nface_orient,  &
!                      nrdofV,shapV,divV)
! !       L2 (trial)
!         call shape3Q(etype,xi,norder, nrdofQ,shapQ)
! !       H(div) (test)
!         call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
! !       L2 (test)
!         call shape3QQ(etype,xi,nordP, nrdofQQ,shapQQ)
! !
! !  .....geometry map
!         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,  &
!                     x,dxdxi,dxidx,rjac,iflag)
!         if (iflag.ne.0) then
!           write(*,1000) Mdle,rjac
!  1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
!           stop
!         endif
! !
! !  .....Change coordinates so the shape functions are on the physical element
! !       H(div) (trial)
!         do k=1,nrdofV
!           shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
!                        + dxdxi(1:3,2)*shapV(2,k)  &
!                        + dxdxi(1:3,3)*shapV(3,k)
!         enddo
!         shapV(1:3,1:nrdofV) = shapV(1:3,1:nrdofV)/rjac
!         divV(1:nrdofV) = divV(1:nrdofV)/rjac
! !       L2 (trial)
!         shapQ(1:nrdofQ) = shapQ(1:nrdofQ)/rjac
! !       H(div) (test)
!         do k=1,nrdofVV
!           shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
!                         + dxdxi(1:3,2)*shapVV(2,k)  &
!                         + dxdxi(1:3,3)*shapVV(3,k)
!         enddo
!         shapVV(1:3,1:nrdofVV) = shapVV(1:3,1:nrdofVV)/rjac
!         divVV(1:nrdofVV) = divVV(1:nrdofVV)/rjac
! !       L2 (test)
!         shapQQ(1:nrdofQQ) = shapQQ(1:nrdofQQ)/rjac
! !
! !  .....integration weight
!         weight = wa*rjac
! !
! !  .....compute the compliance tensor
!         call getA(x, A)
! !
! !  .....get the source term
!         call getf(Mdle,x, fval)
! !
! !  NOTE: Because nrddofVV \neq nrdofQ, there is no convenient way of organizing
! !        the matrices testing all of the n-th basis functions in turn before
! !        moving on to the (n+1)-th basis functions. Therefore, and which will lead
! !        to a sparser system as well, we test all the basis functions representing
! !        \tau and then all those representing v and then all those representing q.
! !
! !
! !    P A R T  1 : go through H(div)^3 test space (this fills the first submatrices)
! !
! !
! !  .....FIRST OUTER loop through enriched H(div) dofs
!         do k1=1,nrdofVV
! !  .......OUTER loop through components
!           do jcomp=1,3
!             m1 = (k1-1)*3+jcomp
! !
! !           E N R I C H E D   L O A D   V E C T O R
! !
! !   0
! !
! !           G R A M   M A T R I X
! !
!             select case(TEST_NORM)
! !
! !   (\sigma,\tau)+(div(\sigma),div(tau))
! !
!             case(2)
!               do m2=m1,3*nrdofVV
!                 icomp = mod(m2-1,3)+1
!                 if (icomp.eq.jcomp) then
!                   k = nk(m1,m2)
!                   k2 = int((m2-1)/3)+1
!                   Gram(k) = Gram(k)  &
!                           + ( divVV(k1)*divVV(k2)  &
!                             + shapVV(1,k1)*shapVV(1,k2)  &
!                             + shapVV(2,k1)*shapVV(2,k2)  &
!                             + shapVV(3,k1)*shapVV(3,k2) )*weight
!                 endif
!               enddo
!             end select
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                          (   C A U C H Y    S T R E S S   )
! !
! !   A \sigma : \tau
! !
! !  .........INNER loop through trial dofs for Cauchy stress
!             do k3=1,nrdofV
!               do icomp=1,3
!                 m3 = (k3-1)*3+icomp
! !
! ! USING REGULAR A
!                 ! tmp = 0.d0
!                 ! do m=1,3; do n=1,3
!                 !   tmp = tmp  &
!                 !       + A(icomp,m,jcomp,n)*shapV(m,k3)*shapVV(n,k1)
!                 ! enddo; enddo
! !
! ! USING EXTENDED A
! !               A \sigma = 1/(2*MU)*(\sigma - LAMBDA/(2*MU+3*LAMBDA)*tr(\sigma)*I)
!                 tmp = -LAMBDA/(2*MU+3*LAMBDA)*shapV(icomp,k3)*shapVV(jcomp,k1)
!                 if (icomp.eq.jcomp) then
!                   do m=1,3
!                     tmp = tmp  &
!                         + shapV(m,k3)*shapVV(m,k1)
!                   enddo
!                 endif
!                 tmp = tmp/(2*MU)
! !
!                 if (abs(tmp).gt.1.0d-15) then
!                   EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
!                                         + tmp*weight
!                 endif
!               enddo
!             enddo
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                          (   D I S P L A C E M E N T   )
! !
! !   u \cdot div(\tau)
! !
! !  .........INNER loop through trial dofs for displacement
!             do k4=1,nrdofQ
!               !  contribute only when components match (populate only the diagonal)
!               m4 = (k4-1)*3+jcomp
!               EnrFieldDispl(m1,m4) = EnrFieldDispl(m1,m4)  &
!                                    + shapQ(k4)*divVV(k1)*weight
!             enddo
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                   (   L A G R A N G E    M U L T I P L I E R   )
! !
! !   - p : \tau
! !
! !      ( 0   p3  -p2 )
! !  p = (-p3   0   p1 )
! !      ( p2 -p1    0 )
! !
! !  .........INNER loop through trial dofs for Lagrange multiplier
!             do k5=1,nrdofQ
!               do icomp=1,3
!                 m5 = (k5-1)*3+icomp
!                 if (icomp.eq.1) then
!                   if (jcomp.eq.2) then
!                     EnrFieldLagra(m1,m5) = EnrFieldLagra(m1,m5)  &
!                                          - shapQ(k5)*shapVV(3,k1)*weight
!                   elseif (jcomp.eq.3) then
!                     EnrFieldLagra(m1,m5) = EnrFieldLagra(m1,m5)  &
!                                          + shapQ(k5)*shapVV(2,k1)*weight
!                   endif
!                 elseif (icomp.eq.2) then
!                   if (jcomp.eq.1) then
!                     EnrFieldLagra(m1,m5) = EnrFieldLagra(m1,m5)  &
!                                          + shapQ(k5)*shapVV(3,k1)*weight
!                   elseif (jcomp.eq.3) then
!                     EnrFieldLagra(m1,m5) = EnrFieldLagra(m1,m5)  &
!                                          - shapQ(k5)*shapVV(1,k1)*weight
!                   endif
!                 elseif (icomp.eq.3) then
!                   if (jcomp.eq.1) then
!                     EnrFieldLagra(m1,m5) = EnrFieldLagra(m1,m5)  &
!                                          - shapQ(k5)*shapVV(2,k1)*weight
!                   elseif (jcomp.eq.2) then
!                     EnrFieldLagra(m1,m5) = EnrFieldLagra(m1,m5)  &
!                                          + shapQ(k5)*shapVV(1,k1)*weight
!                   endif
!                 endif
!               enddo
!             enddo
! !  .....END OUTER LOOP through test stresses
!           enddo
!         enddo
! !
! !  .....SECOND OUTER loop through enriched L2 dofs
!         do k1=1,nrdofQQ
! !  .......OUTER loop through components
!           do jcomp=1,3
!             ! counter for part 2
!             m1 = 3*nrdofVV+(k1-1)*3+jcomp
!             ! counter for part 3
!             n1 = m1+3*nrdofQQ
! !
! !
! !    P A R T  2 : go through (L2)^3 test space (this fills the second submatrices)
! !
! !
! !
! !            E N R I C H E D   L O A D   V E C T O R
! !
! !   f \cdot v
! !
!             EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
!                                  + fval(jcomp,1:NR_RHS)*shapQQ(k1)*weight
! !
! !            G R A M   M A T R I X
! !
!             select case(TEST_NORM)
! !
! !   (u,v)
! !
!             case(2)
!               do m2=m1,3*nrdofVV+3*nrdofQQ
!                 icomp = mod(m2-1,3)+1
!                 if (icomp.eq.jcomp) then
!                   k = nk(m1,m2)
!                   k2 = int((m2-3*nrdofVV-1)/3)+1
!                   Gram(k) = Gram(k) + shapQQ(k2)*shapQQ(k1)*weight
!                 endif
!               enddo
!             end select
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                          (   C A U C H Y    S T R E S S   )
! !
! !   - div(\sigma) \cdot v
! !
! !  .........INNER loop through trial dofs for Cauchy stress
!             do k3=1,nrdofV
!               !  contribute only when components match (populate only the diagonal)
!               m3 = (k3-1)*3+jcomp
!               EnrFieldStress(m1,m3) = EnrFieldStress(m1,m3)  &
!                                     - divV(k3)*shapQQ(k1)*weight
!             enddo
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                          (   D I S P L A C E M E N T   )
! !
! !   0
! !
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                   (   L A G R A N G E    M U L T I P L I E R   )
! !
! !   0
! !
! !
! !
! !    P A R T  3 : go through L2_skew test space (this fills the last submatrices)
! !
! !
! !
! !            E N R I C H E D   L O A D   V E C T O R
! !
! !   0
! !
! !
! !            G R A M   M A T R I X
! !
!             select case(TEST_NORM)
! !
! !   (p,q)
! !
!             case(2)
!               do n2=n1,3*nrdofVV+6*nrdofQQ
!                 icomp = mod(n2-1,3)+1
!                 if (icomp.eq.jcomp) then
!                     k = nk(n1,n2)
!                     k2 = int((n2-3*nrdofVV-3*nrdofQQ-1)/3)+1
!                     Gram(k) = Gram(k) + 2*shapQQ(k2)*shapQQ(k1)*weight
!                 endif
!               enddo
!             end select
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                          (   C A U C H Y    S T R E S S   )
! !
! !   \sigma : q
! !
! !      ( 0   q3  -q2 )
! !  q = (-q3   0   q1 )
! !      ( q2 -q1    0 )
! !
! !  .........INNER loop through trial dofs for Cauchy stress tensor
!             do k3=1,nrdofV
!               do icomp=1,3
!                 n3 = (k3-1)*3+icomp
!                 if (icomp.eq.1) then
!                   if (jcomp.eq.2) then
!                     EnrFieldStress(n1,n3) = EnrFieldStress(n1,n3)  &
!                                           - shapV(3,k3)*shapQQ(k1)*weight
!                   elseif (jcomp.eq.3) then
!                     EnrFieldStress(n1,n3) = EnrFieldStress(n1,n3)  &
!                                           + shapV(2,k3)*shapQQ(k1)*weight
!                   endif
!                 elseif (icomp.eq.2) then
!                   if (jcomp.eq.1) then
!                     EnrFieldStress(n1,n3) = EnrFieldStress(n1,n3)  &
!                                           + shapV(3,k3)*shapQQ(k1)*weight
!                   elseif (jcomp.eq.3) then
!                     EnrFieldStress(n1,n3) = EnrFieldStress(n1,n3)  &
!                                           - shapV(1,k3)*shapQQ(k1)*weight
!                   endif
!                 elseif (icomp.eq.3) then
!                   if (jcomp.eq.1) then
!                     EnrFieldStress(n1,n3) = EnrFieldStress(n1,n3)  &
!                                           - shapV(2,k3)*shapQQ(k1)*weight
!                   elseif (jcomp.eq.2) then
!                     EnrFieldStress(n1,n3) = EnrFieldStress(n1,n3)  &
!                                           + shapV(1,k3)*shapQQ(k1)*weight
!                   endif
!                 endif
!               enddo
!             enddo
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                          (   D I S P L A C E M E N T   )
! !
! !   0
! !
! !           E N R I C H E D   F I E L D   S T I F F N E S S   M A T R I X
! !                   (   L A G R A N G E    M U L T I P L I E R   )
! !
! !   0
! !
! !  ......END OUTER LOOP
!           enddo
!         enddo
! !
! !  ...end of loop through integration points
!       enddo
! !
! !-----------------------------------------------------------------------------------
! !     B O U N D A R Y    I N T E G R A L                                           |
! !-----------------------------------------------------------------------------------
! !
! !  ...loop through element faces
!       do ifc=1,nrf
! !
! !  .....sign factor to determine the OUTWARD normal unit vector
!         nsign = nsign_param(etype,ifc)
! !
! !  .....face type
!         ftype = face_type(etype,ifc)
! !
! !  .....face order of approximation
!         call face_order(etype,ifc,norder, nordf)
! !
! !  .....set up the face quadrature
!         INTEGRATION = NORD_ADD
!         call set_2Dint(ftype,nordf, nint,tloc,wt)
!         INTEGRATION = 0
! !
! !  .....loop through face integration points
!         do ipt=1,nint
! !
! !  .......face coordinates
!           t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
! !
! !  .......master element coordinates using face parameterization
!           call face_param(etype,ifc,t, xi,dxidt)
! !
! !  .......Compute shape functions needed for test/trial field variables and geometry
! !         H1 (geometry and trial)
!           call shape3H(etype,xi,norder,nedge_orient,nface_orient,  &
!                      nrdofH,shapH,gradH)
! !
! !         H(div) (test)
!           call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
! !
! !  .......geometry map
!           call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,  &
!                        x,dxdxi,dxidx,rjac,dxdt,rn,brjac)
!           weight = wa*brjac
! !
! !  .......Change coordinates so the shape functions are on the physical element
! !         H(div) (test)
!           do k=1,nrdofVV
!             shapVV(1:3,k) = dxdxi(1:3,1)*shapVV(1,k)  &
!                           + dxdxi(1:3,2)*shapVV(2,k)  &
!                           + dxdxi(1:3,3)*shapVV(3,k)
!             shapVV_n(k) = shapVV(1,k)*rn(1)  &
!                         + shapVV(2,k)*rn(2)  &
!                         + shapVV(3,k)*rn(3)
!           enddo
!           shapVV_n(1:nrdofVV) = shapVV_n(1:nrdofVV)/rjac
! !
! !
! !    P A R T  1 : go through H(div)^3 test space (this fills the first submatrices)
! !
! !
! !  .......OUTER loop through enriched H1 test functions
!           do k1=1,nrdofVV
!             do jcomp=1,3
!               m1 = (k1-1)*3+jcomp
! !
! !  ...........INNER loop through enriched dofs
!               do k2=1,nrdofH
!                 !  contribute only when components match (populate only the diagonal)
!                 m2 = (k2-1)*3+jcomp
! !
! !
! !           E N R I C H E D   T R A C E   S T I F F N E S S   M A T R I X
! !
! !   - <\hat u,(\tau n)>
! !
!                 EnrTraceDispl(m1,m2) = EnrTraceDispl(m1,m2)  &
!                                      - shapH(k2)*shapVV_n(k1)*weight
!               enddo
! !
! !  .......OUTER loop
!             enddo
!           enddo
! !
! !
! !    P A R T  2 and 3 : no contributions from the remaining test spaces
! !
! !
! !  .....end of loop over integration points
!         enddo
! !
! !  ...end of loop over faces
!       enddo

!       if (iprint.eq.1) then
!         write(*,*) 'Gram = '
!         do i=1,25
!           write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
!  6000     format('i = ',i3,'  ',25e12.5)
!         enddo

!         call pause

!         write(*,*) 'EnrFieldDispl = '
!         do i=1,3*nrdofVV+1
!           write(*,6001) i,EnrFieldDispl(i,1:3*NrdofQ)
!  6001     format('i = ',i4,'  ',15(/,10e12.5))
!         enddo

!         call pause

!         write(*,*) 'EnrFieldStress = '
!         do i=1,3*nrdofVV
!           write(*,6001) i,EnrFieldStress(i,1:3*NrdofV)
!         enddo

!         call pause

!         do i=1+3*nrdofVV,3*nrdofVV+3*nrdofQQ
!           write(*,6001) i,EnrFieldStress(i,1:3*NrdofV)
!         enddo

!         call pause

!         do i=1+3*nrdofVV+3*nrdofQQ,3*nrdofVV+3*nrdofQQ+3*nrdofQQ+1
!           write(*,6001) i,EnrFieldStress(i,1:3*NrdofV)
!         enddo

!         call pause

!         write(*,*) 'EnrFieldLagra = '
!         do i=1,3*nrdofVV+1
!           write(*,6001) i,EnrFieldLagra(i,1:3*NrdofQ)
!         enddo

!         call pause

!         write(*,*) 'EnrTraceDispl = '
!         do i=1,3*nrdofVV+1
!           write(*,6001) i,EnrTraceDispl(i,1:3*NrdofH)
!         enddo

!         call pause

!       endif

! !
! !-----------------------------------------------------------------------------------
! !     D P G    L O C A L    A S S E M B L Y                                        |
! !-----------------------------------------------------------------------------------
! !
! !  ...factor the Gram matrix
!       uplo = 'U'
!       call DPPTRF(uplo, 3*nrdofVV+6*nrdofQQ, Gram, info)
!       if (info.ne.0) then
!         write(*,*) 'elem_DPG_MIXED: info = ',info ; stop
!       endif
! !
! !  ...save copies of enriched stiffness matrices
!       EnrTraceDisplc=EnrTraceDispl; EnrFieldStressc=EnrFieldStress
!       EnrFieldDisplc=EnrFieldDispl; EnrFieldLagrac=EnrFieldLagra
! !
! !  ...G^-1 * Load
!       call DPPTRS(uplo,3*nrdofVV+6*nrdofQQ,NR_RHS,Gram,EnrLoad,3*MAXbrickVV+6*MAXbrickQQ,info1)
!       if (info1.ne.0) then
!         write(*,*) 'elem_DPG_MIXED: info1 = ',info1 ; stop
!       endif
! !
! !  ...G^-1 * EnrTraceDispl
!       call DPPTRS(uplo,3*nrdofVV+6*nrdofQQ,3*nrdofH,Gram,EnrTraceDisplc,3*MAXbrickVV+6*MAXbrickQQ,info2)
!       if (info2.ne.0) then
!         write(*,*) 'elem_DPG_MIXED: info2 = ',info2 ; stop
!       endif
! !
! !  ...G^-1 * EnrFieldStress
!       call DPPTRS(uplo,3*nrdofVV+6*nrdofQQ,3*nrdofV,Gram,EnrFieldStressc,3*MAXbrickVV+6*MAXbrickQQ,info3)
!       if (info3.ne.0) then
!         write(*,*) 'elem_DPG_MIXED: info3 = ',info3 ; stop
!       endif
! !
! !  ...G^-1 * EnrFieldDispl
!       call DPPTRS(uplo,3*nrdofVV+6*nrdofQQ,3*nrdofQ,Gram,EnrFieldDisplc,3*MAXbrickVV+6*MAXbrickQQ,info2)
!       if (info2.ne.0) then
!         write(*,*) 'elem_DPG_MIXED: info2 = ',info2 ; stop
!       endif
! !
! !  ...G^-1 * EnrFieldLagra
!       call DPPTRS(uplo,3*nrdofVV+6*nrdofQQ,3*nrdofQ,Gram,EnrFieldLagrac,3*MAXbrickVV+6*MAXbrickQQ,info2)
!       if (info2.ne.0) then
!         write(*,*) 'elem_DPG_MIXED: info2 = ',info2 ; stop
!       endif
! !
! !
! !     ULTIMATE DPG LOAD VECTORS AND STIFFNESS MATRICES THROUGH STATIC CONDENSATION
! !
! !     ROW 1
!       do kH=1,3*nrdofH
! !
! !       EnrTraceDispl^T * (G^-1 * EnrLoad)
!         do m=1,3*nrdofVV+6*nrdofQQ
!           Bloc1(kH,1:NR_RHS) = Bloc1(kH,1:NR_RHS)  &
!                              + EnrTraceDispl(m,kH)*EnrLoad(m,1:NR_RHS)
!         enddo
! !
! !       EnrTraceDispl^T * (G^-1 * EnrTraceDispl)
!         do lH=1,3*nrdofH
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc11(kH,lH) = Aloc11(kH,lH)  &
!                           + EnrTraceDispl(m,kH)*EnrTraceDisplc(m,lH)
!           enddo
!         enddo
! !
! !       EnrTraceDispl^T * (G^-1 * EnrFieldStress)
!         do lV=1,3*nrdofV
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc12(kH,lV) = Aloc12(kH,lV)  &
!                           + EnrTraceDispl(m,kH)*EnrFieldStressc(m,lV)
!           enddo
!         enddo
! !
!         do lQ=1,3*nrdofQ
!           do m=1,3*nrdofVV+6*nrdofQQ
! !
! !           EnrTraceDispl^T * (G^-1 * EnrFieldDispl)
!             Aloc13(kH,lQ) = Aloc13(kH,lQ)  &
!                           + EnrTraceDispl(m,kH)*EnrFieldDisplc(m,lQ)
! !
! !           EnrTraceDispl^T * (G^-1 * EnrFieldLagra)
!             Aloc14(kH,lQ) = Aloc14(kH,lQ)  &
!                           + EnrTraceDispl(m,kH)*EnrFieldLagrac(m,lQ)
!           enddo
!         enddo
!       enddo
! !
! !     ROW 2
!       do kV=1,3*nrdofV
! !
! !       EnrFieldStress^T * (G^-1 * EnrLoad)
!         do m=1,3*nrdofVV+6*nrdofQQ
!           Bloc2(kV,1:NR_RHS) = Bloc2(kV,1:NR_RHS)  &
!                              + EnrFieldStress(m,kV)*EnrLoad(m,1:NR_RHS)
!         enddo
! !
! !       EnrFieldStress^T * (G^-1 * EnrTraceDispl)
!         do lH=1,3*nrdofH
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc21(kV,lH) = Aloc21(kV,lH)  &
!                           + EnrFieldStress(m,kV)*EnrTraceDisplc(m,lH)
!           enddo
!         enddo
! !
! !       EnrFieldStress^T * (G^-1 * EnrFieldStress)
!         do lV=1,3*nrdofV
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc22(kV,lV) = Aloc22(kV,lV)  &
!                           + EnrFieldStress(m,kV)*EnrFieldStressc(m,lV)
!           enddo
!         enddo
! !
!         do lQ=1,3*nrdofQ
!           do m=1,3*nrdofVV+6*nrdofQQ
! !
! !           EnrFieldStress^T * (G^-1 * EnrFieldDispl)
!             Aloc23(kV,lQ) = Aloc23(kV,lQ)  &
!                           + EnrFieldStress(m,kV)*EnrFieldDisplc(m,lQ)
! !
! !           EnrFieldStress^T * (G^-1 * EnrFieldLagra)
!             Aloc24(kV,lQ) = Aloc24(kV,lQ)  &
!                           + EnrFieldStress(m,kV)*EnrFieldLagrac(m,lQ)
!           enddo
!         enddo
!       enddo
! !
! !     ROW 3
!       do kQ=1,3*nrdofQ
! !
! !       EnrFieldDispl^T * (G^-1 * EnrLoad)
!         do m=1,3*nrdofVV+6*nrdofQQ
!           Bloc3(kQ,1:NR_RHS) = Bloc3(kQ,1:NR_RHS)  &
!                              + EnrFieldDispl(m,kQ)*EnrLoad(m,1:NR_RHS)
!         enddo
! !
! !       EnrFieldDispl^T * (G^-1 * EnrTraceDisp)
!         do lH=1,3*nrdofH
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc31(kQ,lH) = Aloc31(kQ,lH)  &
!                           + EnrFieldDispl(m,kQ)*EnrTraceDisplc(m,lH)
!           enddo
!         enddo
! !
! !       EnrFieldDispl^T * (G^-1 * EnrFieldStress)
!         do lV=1,3*nrdofV
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc32(kQ,lV) = Aloc32(kQ,lV)  &
!                           + EnrFieldDispl(m,kQ)*EnrFieldStressc(m,lV)
!           enddo
!         enddo
! !
!         do lQ=1,3*nrdofQ
!           do m=1,3*nrdofVV+6*nrdofQQ
! !
! !           EnrFieldDispl^T * (G^-1 * EnrFieldDispl)
!             Aloc33(kQ,lQ) = Aloc33(kQ,lQ)  &
!                           + EnrFieldDispl(m,kQ)*EnrFieldDisplc(m,lQ)
! !
! !           EnrFieldDispl^T * (G^-1 * EnrFieldLagra)
!             Aloc34(kQ,lQ) = Aloc34(kQ,lQ)  &
!                           + EnrFieldDispl(m,kQ)*EnrFieldLagrac(m,lQ)
!           enddo
!         enddo
!       enddo
! !
! !     ROW 4
!       do kQ=1,3*nrdofQ
! !
! !       EnrFieldLagra^T * (G^-1 * EnrLoad)
!         do m=1,3*nrdofVV+6*nrdofQQ
!           Bloc4(kQ,1:NR_RHS) = Bloc4(kQ,1:NR_RHS)  &
!                              + EnrFieldLagra(m,kQ)*EnrLoad(m,1:NR_RHS)
!         enddo
! !
! !       EnrFieldLagra^T * (G^-1 * EnrTraceDisp)
!         do lH=1,3*nrdofH
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc41(kQ,lH) = Aloc41(kQ,lH)  &
!                           + EnrFieldLagra(m,kQ)*EnrTraceDisplc(m,lH)
!           enddo
!         enddo
! !
! !       EnrFieldLagra^T * (G^-1 * EnrFieldStress)
!         do lV=1,3*nrdofV
!           do m=1,3*nrdofVV+6*nrdofQQ
!             Aloc42(kQ,lV) = Aloc42(kQ,lV)  &
!                           + EnrFieldLagra(m,kQ)*EnrFieldStressc(m,lV)
!           enddo
!         enddo
! !
!         do lQ=1,3*nrdofQ
!           do m=1,3*nrdofVV+6*nrdofQQ
! !
! !           EnrFieldLagra^T * (G^-1 * EnrFieldDispl)
!             Aloc43(kQ,lQ) = Aloc43(kQ,lQ)  &
!                           + EnrFieldLagra(m,kQ)*EnrFieldDisplc(m,lQ)
! !
! !           EnrFieldLagra^T * (G^-1 * EnrFieldLagra)
!             Aloc44(kQ,lQ) = Aloc44(kQ,lQ)  &
!                           + EnrFieldLagra(m,kQ)*EnrFieldLagrac(m,lQ)
!           enddo
!         enddo
!       enddo
! !
! !-----------------------------------------------------------------------------------
! !     T E S T S    A N D    P R I N T    S T A T E M E N T S                       |
! !-----------------------------------------------------------------------------------
! !
! !  ...check symmetry (not complete)
!       diffmax = ZERO; dmax = ZERO
!       do kH=1,3*nrdofH
!         do lH=kH,3*nrdofH
!           diffmax = max(diffmax,abs(Aloc11(kH,lH)-Aloc11(lH,kH)))
!           dmax = max(dmax,abs(Aloc11(kH,lH)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7021) diffmax, dmax
!  7021   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc11 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kV=1,3*nrdofV
!         do lV=kV,3*nrdofV
!           diffmax = max(diffmax,abs(Aloc22(kV,lV)-Aloc22(lV,kV)))
!           dmax = max(dmax,abs(Aloc22(kV,lV)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7022) diffmax, dmax
!  7022   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc22 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kQ=1,3*nrdofQ
!         do lQ=kQ,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc33(kQ,lQ)-Aloc33(lQ,kQ)))
!           dmax = max(dmax,abs(Aloc33(kQ,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7023) diffmax, dmax
!  7023   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc33 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kQ=1,3*nrdofQ
!         do lQ=kQ,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc44(kQ,lQ)-Aloc44(lQ,kQ)))
!           dmax = max(dmax,abs(Aloc44(kQ,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7024) diffmax, dmax
!  7024   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc44 = ',2e12.5)
!         call pause
!       endif
! !
!       diffmax = ZERO; dmax = ZERO
!       do kH=1,3*nrdofH
!         do lV=1,3*nrdofV
!           diffmax = max(diffmax,abs(Aloc12(kH,lV)-Aloc21(lV,kH)))
!           dmax = max(dmax,abs(Aloc12(kH,lV)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7025) diffmax, dmax
!  7025   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc12 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kH=1,3*nrdofH
!         do lQ=1,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc13(kH,lQ)-Aloc31(lQ,kH)))
!           dmax = max(dmax,abs(Aloc13(kH,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7026) diffmax, dmax
!  7026   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc13 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kH=1,3*nrdofH
!         do lQ=1,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc14(kH,lQ)-Aloc41(lQ,kH)))
!           dmax = max(dmax,abs(Aloc14(kH,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7027) diffmax, dmax
!  7027   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc14 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kV=1,3*nrdofV
!         do lQ=1,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc23(kV,lQ)-Aloc32(lQ,kV)))
!           dmax = max(dmax,abs(Aloc23(kV,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7028) diffmax, dmax
!  7028   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc23 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kV=1,3*nrdofV
!         do lQ=1,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc24(kV,lQ)-Aloc42(lQ,kV)))
!           dmax = max(dmax,abs(Aloc24(kV,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7029) diffmax, dmax
!  7029   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc24 = ',2e12.5)
!         call pause
!       endif
!       diffmax = ZERO; dmax = ZERO
!       do kQ=1,3*nrdofQ
!         do lQ=1,3*nrdofQ
!           diffmax = max(diffmax,abs(Aloc34(kQ,lQ)-Aloc43(lQ,kQ)))
!           dmax = max(dmax,abs(Aloc34(kQ,lQ)))
!         enddo
!       enddo
!       if (diffmax/dmax.gt.SYMMETRY_TOL) then
!         write(*,7029) diffmax, dmax
!  7030   format('elem_DPG_MIXED: diffmax,dmax FOR Aloc34 = ',2e12.5)
!         call pause
!       endif
! !
! !  ...print statments (not complete)
!       if (iprint.ge.1) then
!         write(*,7010)
!  7010   format('elem_DPG_MIXED: Bloc1,Bloc2 = ')
!         write(*,7011) Bloc1(1:3*NrdofH,1)
!         write(*,7011) Bloc2(1:3*NrdofV,1)
!  7011   format(10e12.5)
!         write(*,7012)
!  7012   format('elem_DPG_MIXED: Aloc11 = ')
!         do i=1,3*NrdofH
!           write(*,7013) i,Aloc11(i,1:3*NrdofH)
!  7013     format('i = ',i3,10(/,10e12.5))
!         enddo
!         write(*,7014)
!  7014   format('elem_DPG_MIXED: Aloc12 = ')
!         do i=1,3*NrdofH
!           write(*,7013) i,Aloc12(i,1:3*NrdofV)
!         enddo
!         write(*,7015)
!  7015   format('elem_DPG_MIXED: Aloc22 = ')
!         do i=1,3*NrdofV
!           write(*,7013) i,Aloc22(i,1:3*NrdofV)
!         enddo
!       endif
! !
! !
end subroutine elem_DPG_UWEAK
