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
  use control
  use physics
  use assembly
  use error
  use data_structure3D
!--------------------------------------------------------------------------
  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!--------------------------------------------------------------------------
  !
  Itest(1)=1; Itrial(1)=1
  !
  ALOC(1,1)%array = ZERO ; BLOC(1)%array = ZERO
  !
  select case(NODES(Mdle)%case)
  ! elasticity
  case(1)
     call elem_elast(Mdle,BLOC(  1)%array,BLOC(  1)%nrow,  &
                          ALOC(1,1)%array,ALOC(1,1)%ncol)
  case default
    write(*,*) 'mdle,case=',mdle,NODES(Mdle)%case
     call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
  end select
!
end subroutine elem

!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for elasticity problem
!! @param[in]  Mdle      - middle node number
!! @param[out] Bloc      - elem load vector(s)
!! @param[in]  Nrow_Bloc - number of rows of load vector
!! @param[out] Aloc      - elem stiffness matrix
!! @param[in]  Ncol_Aloc - number of columns of stiffness matrix
!------------------------------------------------------------------------------------------
!
!  Since integrals are computed on the master element, it is convenient to transform the
!  coefficients of the differential operator from physical to master coordinates:
!
!    2nd order term (Einstein summation convention)
!
!      EE_ijkl dv_i/dx_j du_k/dx_l =
!                                                                (chain rule)
!      EE_ijkl (dv_i/dxi_m dxi_m/dx_j) (du_k/dxi_n dxi_n/dx_l) =
!                                                                (rearrange)
!     (EE_ijkl dxi_m/dx_j dxi_n/dx_l) dv_i/dxi_m du_k/dxi_n
!                                                                (define)
!     \tilde{EE}_imkn dv_i/dxi_m du_k/dxi_n
!
!------------------------------------------------------------------------------------------
!
subroutine elem_elast(Mdle,Bloc,Nrow_Bloc,Aloc,Ncol_Aloc)
      use control , only: EXGEOM
      use data_structure3D
      use element_data
! #if TRANSVERSE_MODE
!       use transverse_isotropic_elast_material
! #else
!       use isotropic_elast_material
! #endif
      use assembly , only: NR_RHS
      use mpi_param, only: RANK
      use hyperelasticity
      use nl_solver_module, only: LINESEARCH_FACTOR

!------------------------------------------------------------------------------------------
      implicit none
      integer,                               intent(in)  :: Mdle
      integer,                               intent(in)  :: Nrow_Bloc,Ncol_Aloc
      real*8, dimension(Nrow_Bloc,NR_RHS),   intent(out) :: Bloc
      real*8, dimension(Nrow_Bloc,Ncol_Aloc),intent(out) :: Aloc
!------------------------------------------------------------------------------------------
!
!  ...element and face type
      character(len=4) :: etype,ftype
!
!  ...element and face order
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
!
!  ...node edge and face orientations
      integer, dimension(12) :: nedge_orient
      integer, dimension(6)  :: nface_orient
!
!  ...shape functions and their derivatives wrt master coordinates
      real*8, dimension(  MAXbrickH) :: vshapH
      real*8, dimension(3,MAXbrickH) :: dvshapH
!
!  ...approximate solution
      real*8, dimension(MAXEQNH,MAXbrickH) :: dofH
      real*8, dimension(MAXEQNE,MAXbrickE) :: dofE
      real*8, dimension(MAXEQNV,MAXbrickV) :: dofV
      real*8, dimension(MAXEQNQ,MAXbrickQ) :: dofQ
      real*8, dimension(  MAXEQNH  )       :: solH
      real*8, dimension(  MAXEQNH,3)       :: dsolH
      real*8, dimension(3,MAXEQNE  )       :: solE
      real*8, dimension(3,MAXEQNE  )       :: curlE
      real*8, dimension(3,MAXEQNV  )       :: solV
      real*8, dimension(  MAXEQNV  )       :: divV
      real*8, dimension(  MAXEQNQ  )       :: solQ
      real*8, dimension(  NRHVAR  )       :: solH_up
      real*8, dimension(  NRHVAR,3)       :: dsolH_up
      ! real*8, dimension(3,NREVAR  )       :: solE_up
      ! real*8, dimension(3,NREVAR  )       :: curlE_up
      ! real*8, dimension(3,NRVVAR  )       :: solV_up
      ! real*8, dimension(  NRVVAR  )       :: divV_up
      ! real*8, dimension(  NRQVAR  )       :: solQ_up
!
!  ...geometry
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt
!
!  ...elasticities in master and physical coordinates,
      real*8, dimension(3,3,3,3) :: EE,za
      real*8, dimension(3,3) :: c
! ....hyperelasticity 
      real*8 :: W,dWdF(3,3),d2WdF(3,3,3,3), Ftensor(3,3)
      integer:: imat
!
!  ...source term, Neumann term
      real*8, dimension(3,MAXNRHS) :: zfval,zgval
!
!  ...3D quadrature data
      real*8, dimension(3,MAX_NINT3) :: xiloc
      real*8, dimension(MAX_NINT3)   :: wxi
!
!  ...2D quadrature data for boundary terms
      real*8, dimension(2,MAXquadH) :: tloc
      real*8, dimension(MAXquadH)   :: wt
!
!  ...external unit normal vector
      real*8, dimension(3) :: rn
!
!  ...BC's flags
      integer, dimension(6,NR_PHYSA) :: ibc
!
!     miscellaneous
      integer :: i,j,k,l,m,n,k1,k2,ipt,icomp,kcomp,nint,ifig,nord,m1,m2, &
                 nsign,iprint,nrdofH,iflag,ie,ndof
      integer :: ibeg,iend,jbeg,jend,istep,nordh,nordv,jstep
      real*8  :: weight,wa,rjac,zs
!
!-----------------------------------------------------------------------------------
      iprint=0
!
!  ...order of approximation, orientations, geometry dof's, BC flags
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
      call find_bc    (Mdle, ibc)
!
!  ...clear spaces for the element matrices
      Bloc = ZERO
      Aloc = ZERO
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L S                                          |
!-----------------------------------------------------------------------------------
!

!
!  ...determine solution dof
      call solelm(Mdle, dofH,dofE,dofV,dofQ)



!  ...set up the element quadrature
      etype = NODES(Mdle)%type
      call set_3Dint(etype,norder, nint,xiloc,wxi)
!
!  ...loop through integration points
      do ipt=1,nint
        xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!  .....derivatives and values of the shape functions at the point
        call shape3DH(etype,xi,norder,nedge_orient,nface_orient,nrdofH, &
                     vshapH,dvshapH)
!
!  .....geometry map
        select case(EXGEOM)
        case(0)
          x(1:3) = 0.d0 ; dxdxi(1:3,1:3) = 0.d0
          do k=1,nrdofH
            x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
            do i=1,3
              dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dvshapH(i,k)
            enddo
          enddo
        case(1)
          call exact_geom(Mdle,xi, x,dxdxi)
        end select
!
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
        if (iflag.ne.0) then
          write(*,1000) Mdle,rjac
 1000     format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
          stop
        endif


!
!  .....compute the approximate solution on the PHYSICAL element
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,dofH,dofE,dofV,dofQ,1, &
                     x,dxdxi,solH,dsolH,solE,curlE,solV,divV,solQ)


        ! solH_up (1:NRHVAR  ) = solH (NRHVAR+1:2*NRHVAR  ) + LINESEARCH_FACTOR*solH (1:NRHVAR  )
        dsolH_up(1:NRHVAR,:) = dsolH(NRHVAR+1:2*NRHVAR,:) + LINESEARCH_FACTOR*dsolH(1:NRHVAR,:)

        Ftensor = DEL + dsolH_up

        call find_material(Mdle,imat)

        if (MATERIALS(Imat)%FLAG_INCOM.and.MATERIALS(Imat)%CONSTIT.eq.0) then
          write(*,*) 'elem: primal formulation does not support an incompressible linear elastic material!'
          stop
        endif

        call eval_strain_energy_W_F(imat,x,Ftensor,W,dWdF,d2WdF)
!
!  .....Change coordinates so the shape functions are on the physical element
        do k=1,nrdofH
          dvshapH(1:3,k) = dvshapH(1,k)*dxidx(1,1:3)  &
                         + dvshapH(2,k)*dxidx(2,1:3)  &
                         + dvshapH(3,k)*dxidx(3,1:3)
        enddo
!
!  ...compute the elastictiy tensor
      ! call getC(X, EE)
! !
! !  .....transform coefficients to master coordinates
!         za = ZERO ; c = 0.d0
! !
! !  .....loop through displacement components
!         do m=1,3; do n=1,3
! !  .......loop through coordinates
!           do i=1,3; do k=1,3
!             do j=1,3; do l=1,3
!               za(i,m,k,n) = za(i,m,k,n)  &
!                           + dxidx(m,j)*dxidx(n,l)*EE(i,j,k,l)
!             enddo; enddo
!             ! if (i.eq.k) c(i,k) = - RHO*OMEG**2
!           enddo; enddo
!         enddo; enddo
!
!  .....total weight
        weight = wa*rjac
!
!  .....get the source term
        call getf(Mdle,x, zfval)
!
!  .....OUTER loop through dofs
        do k1=1,nrdofH
          do icomp=1,3
            m1 = (k1-1)*3+icomp
!
!           L O A D   V E C T O R
            Bloc(m1,1:NR_RHS) = Bloc(m1,1:NR_RHS) +                        &
                                ( zfval(icomp,1:NR_RHS)*vshapH(k1)         &
                                 -dWdF(icomp,1)*dvshapH(1,k1)              &
                                 -dWdF(icomp,2)*dvshapH(2,k1)              &
                                 -dWdF(icomp,3)*dvshapH(3,k1)      )*weight
!
!  .........INNER loop through dofs
            do k2=1,nrdofH
              do kcomp=1,3
                m2 = (k2-1)*3+kcomp
!
!               S T I F F N E S S   M A T R I X
                zs = ZERO
                do m=1,3; do n=1,3
                  zs = zs +  &
                       d2WdF(icomp,m,kcomp,n)*dvshapH(m,k1)*dvshapH(n,k2)
                enddo; enddo
                ! zs = zs + c(icomp,kcomp)*vshapH(k1)*vshapH(k2)
                Aloc(m1,m2) = Aloc(m1,m2) + zs*weight
!
!  .........INNER loop
              enddo
            enddo
!  .....OUTER loop
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
!  ...loop through element faces
      do ifig=1,nface(etype)
!
!  .....skip if not a Neumann interface
        if ((ibc(ifig,1).ne.2) .and. (ibc(ifig,1).ne.8)) cycle
!
!  .....determine order for the face
        call face_order(etype,ifig, norder,nordf)
!
!  .....set up the face quadrature
        ftype = face_type(etype,ifig)
        call set_2Dint(ftype,nordf, nint,tloc,wt)
        nsign = Nsign_param(etype,ifig)
!
!  .....loop through face integration points
        do ipt=1,nint
          t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
!
!  .......determine master element coordinates using face parameterization
          call face_param(etype,ifig,t, xi,dxidt)
!
!  .......derivatives and values of the shape functions at the point
          call shape3DH(etype,xi,norder,nedge_orient,nface_orient,nrdofH,  &
                       vshapH,dvshapH)
!
!  .......geometry map
          select case(EXGEOM)
          case(0)
            x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0
            do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
              do i=1,3
                dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dvshapH(i,k)
              enddo
            enddo
          case(1)
            call exact_geom(Mdle,xi,  x,dxdxi)
          end select
!
!  .......determine normal unit vector and jacobian
          dxdt(1:3,1:2) = 0.d0
          do i=1,2
            do j=1,3
              dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
            enddo
          enddo
          call cross_product(dxdt(1:3,1),dxdt(1:3,2), rn)
          call norm(rn, rjac)
          rn(1:3) = rn(1:3)*nsign/rjac
          weight = wa*rjac
!
!  .......get the Neumann data
          call getg(Mdle,ibc(ifig,1),x,rn, zgval)
          if (iprint.eq.1) then
            write(*,4005) ipt,zgval(1:3,1)
 4005       format('elem_A: ipt,zgval = ',i3,2x,3(2e12.5,2x))
          endif
!
!  .......NEUMANN
          if (ibc(ifig,1).eq.2) then
!
!  .........loop through test functions
            do k1=1,nrdofH
              do icomp=1,3
                m1 = (k1-1)*3+icomp
!
!  .............accumulate for the load vector
                  Bloc(m1,1:NR_RHS) = Bloc(m1,1:NR_RHS)   &
                            + zgval(icomp,1:NR_RHS)*vshapH(k1)*weight
              enddo
            enddo
!
!  .......MIXED (Neumann on first 2 components, Dirichlet on last)
          elseif (ibc(ifig,1).eq.8) then
!
!  .........loop through test functions
            do k1=1,nrdofH
              do icomp=1,2
                m1 = (k1-1)*3+icomp
!
!  .............accumulate for the load vector
                  Bloc(m1,1:NR_RHS) = Bloc(m1,1:NR_RHS)   &
                            + zgval(icomp,1:NR_RHS)*vshapH(k1)*weight
              enddo
            enddo
          endif
!
!  .......loop through test functions
          do k1=1,nrdofH
            do icomp=1,3
              m1 = (k1-1)*3+icomp
!
!  ...........accumulate for the load vector
                Bloc(m1,1:NR_RHS) = Bloc(m1,1:NR_RHS)   &
                          + zgval(icomp,1:NR_RHS)*vshapH(k1)*weight
            enddo
          enddo
!
!  .....end of loop over integration points
        enddo
!
!  ...end of loop over faces
      enddo
!
!
!
        ! write(*,*) 'elem_E: Bloc = '
        ! do i=1,Nr_RHS
        !   write(*,7004) Bloc(1:3*nrdofH,i)
        ! enddo
        ! write(*,*) 'elem_E: Aloc = '
        ! do j=1,3*nrdofH,3
        !   write(*,7004) (Aloc(j,i),i=1,3*nrdofH,3)
        ! enddo

!
      if (iprint.ge.1) then
        select case(etype)
        case('mdlp')
          write(*,7011) 1,18
 7011     format('elem_E: VERTEX SHAPE FUNCTIONS = ',i2,'...',i2)
          k=18
          do ie=1,9
            ndof = (norder(ie)-1)*3
            if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
 7012       format('        EDGE ',i2,' SHAPE FUNCTIONS = ',i2,'...',i2)
            k = k+ndof
          enddo
          do ifig=1,5
            nord = norder(9+ifig)
            select case(ifig)
            case(1,2); ndof = (nord-2)*(nord-1)/2*3
            case(3,4,5); call decode(nord, nordh,nordv)
              ndof=(nordh-1)*(nordv-1)*3
            end select
            if (k+ndof.ge.k+1) write(*,7013) ifig,k+1,k+ndof
 7013       format('        FACE ',i2,' SHAPE FUNCTIONS = ',i2,'...',i2)
            k=k+ndof
          enddo
          if (nrdofH*3.ge.k+1) write(*,7014) k+1,nrdofH*3
 7014     format('        MIDDLE NODE SHAPE FUNCTIONS = ',i2,'...',i2)
        case('mdln')
          write(*,7011) 1,12
          k=12
          do ie=1,6
            ndof = (norder(ie)-1)*3
            if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
            k = k+ndof
          enddo
          do ifig=1,4
            nord = norder(6+ifig)
            ndof = (nord-2)*(nord-1)/2*3
            if (k+ndof.ge.k+1) write(*,7013) ifig,k+1,k+ndof
            k=k+ndof
          enddo
          if (nrdofH*3.ge.k+1) write(*,7014) k+1,nrdofH*3
        end select
!
        write(*,*) 'elem_E: Bloc = '
        do i=1,Nr_RHS
          write(*,7004) Bloc(1:3*nrdofH,i)
        enddo
#if C_MODE
 7004   format(5(2e11.4,2x))
#else
 7004   format(10e12.5)
#endif
 50     write(*,*) 'elem_E: Aloc = '
        write(*,*) 'elem_E: SET jbeg,jend,jstep,ibeg,iend,istep'
        read(*,*) jbeg,jend,jstep, ibeg,iend,istep
!!!        if (jbeg.eq.0) go to 60
        do j=jbeg,jend,jstep
          write(*,7005) j,(Aloc(j,i),i=ibeg,iend,istep)
#if C_MODE
 7005     format(i2,2x,5(2e11.4,2x))
#else
 7005     format(i2,2x,10e12.5)
#endif
        enddo
        go to 50
      endif
!
!
end subroutine elem_elast
