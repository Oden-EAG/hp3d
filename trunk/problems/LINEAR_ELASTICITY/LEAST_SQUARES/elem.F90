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

  ALOC(1,1)%array = ZERO; ALOC(1,2)%array = ZERO; BLOC(1)%array = ZERO
  ALOC(2,1)%array = ZERO; ALOC(2,2)%array = ZERO; BLOC(2)%array = ZERO

  select case(NODES(Mdle)%case)
  !  we wish to support both the displacement field variables and the Cauchy stress variables at each point
  !  (all physical attributes are supported when NODES(Mdle)%case == 2**NR_PHYSA-1)
  case(3)

    Itest(1:2)=1; Itrial(1:2)=1
    call elem_LEAST_SQUARES(  &
                Mdle,BLOC(1)%nrow,BLOC(2)%nrow,  &
                ALOC(1,1)%array,ALOC(1,2)%array,BLOC(1)%array,  &
                ALOC(2,1)%array,ALOC(2,2)%array,BLOC(2)%array)
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
!! @param[in]  Nrow_B1   - number of rows of 1-component of load vector
!! @param[in]  Nrow_B2   - number of rows of 2-component of load vector
!!
!! @param[out] Aloc11,Aloc12,Aloc11,Aloc12    - elem stiffness matrix
!! @param[out] Bloc1,Bloc2                    - elem load vector(s)
!------------------------------------------------------------------------------------------
!
subroutine elem_LEAST_SQUARES(Mdle,Nrow_Bloc1,Nrow_Bloc2,  &
                            Aloc11,Aloc12,Bloc1,Aloc21,Aloc22,Bloc2)
      use control, only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use isotropic_elast_material
      use assembly, only: NR_RHS
      use common_prob_data, only: SYMMETRY_TOL
!------------------------------------------------------------------------------------------
      implicit none
      integer,                                  intent(in)  :: Mdle
      integer,                                  intent(in)  :: Nrow_Bloc1
      integer,                                  intent(in)  :: Nrow_Bloc2
      real*8, dimension(Nrow_Bloc1,Nrow_Bloc1), intent(out) :: Aloc11
      real*8, dimension(Nrow_Bloc1,Nrow_Bloc2), intent(out) :: Aloc12
      real*8, dimension(Nrow_Bloc2,Nrow_Bloc1), intent(out) :: Aloc21
      real*8, dimension(Nrow_Bloc2,Nrow_Bloc2), intent(out) :: Aloc22
      real*8, dimension(Nrow_Bloc1,NR_RHS),     intent(out) :: Bloc1
      real*8, dimension(Nrow_Bloc2,NR_RHS),     intent(out) :: Bloc2
!------------------------------------------------------------------------------------------
!
!  ...element and face type
      character(len=4) :: etype
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
!     H1
      real*8, dimension(  MAXbrickH) :: shapH
      real*8, dimension(3,MAXbrickH) :: gradH
      integer                        :: nrdofH
!     H(div)
      real*8, dimension(3,MAXbrickV) :: shapV
      real*8, dimension(  MAXbrickV) :: divV
      integer                        :: nrdofV
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x,rn
      real*8, dimension(3,3)         :: dxdxi,dxidx
!
!  ...stiffness tensors in physical coordinates
      real*8, dimension(3,3,3,3) :: C
!
!  ...source term (don't need Neumann term)
      real*8, dimension(3,MAXNRHS) :: fval
!
!  ...3D quadrature data
      real*8, dimension(3,MAXNINT3ADD) :: xiloc
      real*8, dimension(MAXNINT3ADD)   :: wxi
!
!  ...miscellaneous
      integer :: i,j,k,l,m,n,k1,k2,k_H,l_H,k_V,l_V,kk_H,ll_H,kk_V,ll_V,  &
                 ipt,icomp,kcomp,nint,iprint,iflag
      real*8  :: weight,wa,rjac,tmp,diffmax,dmax
!
!  ...LAPACK stuff
      character uplo
! NOTE: nk is a "statement function"
      integer :: nk,ij
      nk(k1,k2) = (k2-1)*k2/2+k1
      ij(i,j)   = (i-1)*3+j
!
!-----------------------------------------------------------------------------------
!      I N I T I A L I Z A T I O N                                                 |
!-----------------------------------------------------------------------------------
!
      iprint=0
!
!  ...element type
      etype = NODES(Mdle)%type
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...initialize the local element matrices and load vectors
      Aloc11 = ZERO; Aloc12 = ZERO; Bloc1 = ZERO
      Aloc21 = ZERO; Aloc22 = ZERO; Bloc2 = ZERO
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L                                            |
!-----------------------------------------------------------------------------------
!
!  ...set up the element quadrature
      call set_3Dint(etype,norder, nint,xiloc,wxi)
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....Compute shape functions needed for test/trial field variables
!       H1 (trial/geometry)
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,  &
                     nrdofH,shapH,gradH)
!       H(div) (trial)
        call shape3DV(etype,xi,norder,nface_orient,  &
                     nrdofV,shapV,divV)
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
!       H1 (trial)
        do k=1,nrdofH
          gradH(1:3,k) = gradH(1,k)*dxidx(1,1:3)  &
                       + gradH(2,k)*dxidx(2,1:3)  &
                       + gradH(3,k)*dxidx(3,1:3)
        enddo
!       H(div) (trial)
        do k=1,nrdofV
          shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                       + dxdxi(1:3,2)*shapV(2,k)  &
                       + dxdxi(1:3,3)*shapV(3,k)
        enddo
        shapV(1:3,1:nrdofV) = shapV(1:3,1:nrdofV)/rjac
        divV(1:nrdofV) = divV(1:nrdofV)/rjac
!
!  .....integration weight
        weight = wa*rjac
!
!  .....compute the stiffness tensor
        call getC(x, C)
!
!  .....get the source term
        call getf(Mdle,x, fval)
!
!
!    P A R T  1 : go through (H^1)^3 test space (this fills the top half of the stiffness matrix)
!
!
!  .....OUTER loops through test displacement
        do k_H=1,nrdofH
          do kcomp=1,3
            kk_H = (k_H-1)*3+kcomp
!
!           L O A D   V E C T O R
!
!   0
!
!           S T I F F N E S S   M A T R I X
!   (   D I S P L A C E M E N T   -   D I S P L A C E M E N T  )
!
!   ( C:grad(u) , C:grad(v) )
!
!  .........INNER loops through trial displacement
            do l_H=1,nrdofH
              do icomp=1,3
                ll_H = (l_H-1)*3+icomp
                tmp = ZERO
                ! ! ( grad(u) , grad(v) )
                ! if (kcomp.eq.icomp) then
                !     tmp = tmp  &
                !         + 0.5d0 *dot_product(gradH(:,l_H)*gradH(:,k_H)) &
                !         - 0.25d0*gradH(icomp,lH)*gradH(kcomp,k_H)
                ! else 
                !     tmp = tmp &
                !         + 0.5d0*gradH(icomp,lH)*gradH(kcomp,k_H)
                ! endif
                !
                do j=1,3; do l=1,3
                  do m=1,3; do n=1,3
                    tmp = tmp  &
                        + C(m,n,icomp,j)*gradH(j,l_H)  &
                         *C(m,n,kcomp,l)*gradH(l,k_H)
                  enddo; enddo
                enddo; enddo
                Aloc11(kk_H,ll_H) = Aloc11(kk_H,ll_H) + tmp*weight
              enddo
            enddo
!
!           S T I F F N E S S   M A T R I X
!   (   C A U C H Y    S T R E S S   -   D I S P L A C E M E N T  )
!
!   - ( sigma , C:grad(v) )
!
!  .........INNER loops through trial Cauchy stress
            do l_V=1,nrdofV
              do icomp=1,3
                ll_V = (l_V-1)*3+icomp
                tmp = ZERO
                do l=1,3; do j=1,3
                  tmp = tmp  &
                      - shapV(j,l_V)*C(icomp,j,kcomp,l)*gradH(l,k_H)
                enddo; enddo
                Aloc12(kk_H,ll_V) = Aloc12(kk_H,ll_V) + tmp*weight
              enddo
            enddo
!
!  .....END OUTER loops
          enddo
        enddo
!
!
!    P A R T  2 : go through (H(div))^3 test space (this fills the bottom half of the stiffness matrix)
!
!
!  .....OUTER loops through test Cauchy stress
        do k_V=1,nrdofV
          do kcomp=1,3
            kk_V = (k_V-1)*3+kcomp
!
!            L O A D   V E C T O R
!
!   - ( f , div(tau) )
!
            Bloc2(kk_V,1:NR_RHS) = Bloc2(kk_V,1:NR_RHS) &
                                 - fval(kcomp,1:NR_RHS)*divV(k_V)*weight
!
!           F I E L D   S T I F F N E S S   M A T R I X
!   (   D I S P L A C E M E N T   -   C A U C H Y    S T R E S S  )
!
!   - ( C:grad(u) , tau )
!
!  .........INNER loops through trial displacement
            do l_H=1,nrdofH
              do icomp=1,3
                ll_H = (l_H-1)*3+icomp
                tmp = ZERO
                do l=1,3; do j=1,3
                  tmp = tmp  &
                      - C(kcomp,l,icomp,j)*gradH(j,l_H)*shapV(l,k_V)
                enddo; enddo
                Aloc21(kk_V,ll_H) = Aloc21(kk_V,ll_H) + tmp*weight
              enddo
            enddo
!
!           F I E L D   S T I F F N E S S   M A T R I X
!   (   C A U C H Y    S T R E S S   -   C A U C H Y    S T R E S S  )
!
!   ( sigma , tau) + ( div(sigma) , div(tau) )
!
!  .........INNER loops through trial Cauchy stress
            do l_V=1,nrdofV
              !  contribute only when components match
              ll_V = (l_V-1)*3+kcomp
              Aloc22(kk_V,ll_V) = Aloc22(kk_V,ll_V)  &
                                + ( shapV(1,l_V)*shapV(1,k_V)  &
                                   +shapV(2,l_V)*shapV(2,k_V)  &
                                   +shapV(3,l_V)*shapV(3,k_V)  &
                                   +divV(l_V)*divV(k_V) )*weight
            enddo
!
!  .....END OUTER loops
          enddo
        enddo

!
!  ...end of loop through integration points
     enddo
!
!-----------------------------------------------------------------------------------
!     T E S T S    A N D    P R I N T    S T A T E M E N T S                       |
!-----------------------------------------------------------------------------------
!
      if (iprint.ge.2) then
  !  ...check symmetry
        diffmax = ZERO; dmax = ZERO
        do k1=1,3*nrdofH
          do k2=k1,3*nrdofH
            diffmax = max(diffmax,abs(Aloc11(k1,k2)-Aloc11(k2,k1)))
            dmax = max(dmax,abs(Aloc11(k1,k2)))
          enddo
        enddo
        if (diffmax/dmax.gt.SYMMETRY_TOL) then
          write(*,7021) diffmax, dmax
   7021   format('elem_LEAST_SQUARES: diffmax,dmax FOR Aloc11 = ',2e12.5)
          call pause
        endif
        diffmax = ZERO; dmax = ZERO
        do k1=1,nrdofV
          do k2=k1,nrdofV
            diffmax = max(diffmax,abs(Aloc22(k1,k2)-Aloc22(k2,k1)))
            dmax = max(dmax,abs(Aloc22(k1,k2)))
          enddo
        enddo
        if (diffmax/dmax.gt.SYMMETRY_TOL) then
          write(*,7022) diffmax, dmax
   7022   format('elem_LEAST_SQUARES: diffmax,dmax FOR Aloc22 = ',2e12.5)
          call pause
        endif
        diffmax = ZERO; dmax = ZERO
        do k1=1,nrdofH
          do k2=1,nrdofV
            diffmax = max(diffmax,abs(Aloc12(k1,k2)-Aloc21(k2,k1)))
            dmax = max(dmax,abs(Aloc12(k1,k2)))
          enddo
        enddo
        if (diffmax/dmax.gt.SYMMETRY_TOL) then
          write(*,7023) diffmax, dmax
   7023   format('elem_LEAST_SQUARES: diffmax,dmax FOR Aloc12 = ',2e12.5)
          call pause
        endif
      endif
!
!  ...print statments
      if (iprint.ge.1) then
        write(*,7010)
 7010   format('elem_LEAST_SQUARES: Bloc1,Bloc2 = ')
        write(*,7011) Bloc1(1:3*NrdofH,1)
        write(*,7011) Bloc2(1:3*NrdofV,1)
 7011   format(10e12.5)
        write(*,7012)
 7012   format('elem_LEAST_SQUARES: Aloc11 = ')
        do i=1,3*NrdofH
          write(*,7013) i,Aloc11(i,1:3*NrdofH)
 7013     format('i = ',i3,10(/,10e12.5))
        enddo
        write(*,7014)
 7014   format('elem_LEAST_SQUARES: Aloc12 = ')
        do i=1,3*NrdofH
          write(*,7013) i,Aloc12(i,1:3*NrdofV)
        enddo
        write(*,7015)
 7015   format('elem_LEAST_SQUARES: Aloc22 = ')
        do i=1,3*NrdofV
          write(*,7013) i,Aloc22(i,1:3*NrdofV)
        enddo
      endif
!
!
end subroutine elem_LEAST_SQUARES
