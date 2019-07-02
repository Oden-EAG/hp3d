!------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for elasticity problem 
!! @param[in]  Mdle      - middle node number
!! @param[out] Bloc      - elem load vector(s)
!! @param[in]  Nrow_Bloc - number of rows of load vector
!! @param[out] Aloc      - elem stiffness matrix
!! @param[in]  Ncol_Aloc - number of columns of stiffness matrix
!------------------------------------------------------------------------------------------
!
subroutine elem_elast(Mdle,Bloc,Nrow_Bloc,Aloc,Ncol_Aloc)
      use control , only: EXGEOM
      use data_structure3D
      use element_data
      use dome
      use assembly , only: NR_RHS
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
!  ...geometry 
      real*8, dimension(3,MAXbrickH) :: xnod
      real*8, dimension(3)           :: xi,x
      real*8, dimension(3,3)         :: dxdxi,dxidx
      real*8, dimension(2)           :: t
      real*8, dimension(3,2)         :: dxidt,dxdt
!
!     Kronecker's delta
      real*8, parameter, dimension(3,3) :: del = &
        reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0 /),(/3,3/))
!  ...elasticities in master and physical coordinates,
      real*8, dimension(3,3,3,3) :: zaa,za
      real*8, dimension(3,3) :: c
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
!  ...external unit vector
      real*8, dimension(3) :: rn
!
!  ...BC's flags
      integer, dimension(6,NR_PHYSA) :: ibc
!
      integer :: i,ii,j,jj,k,kk,k1,k2,l,ll,nint,ivar1,ivar2,ifig,nord,m1,m2, &
                 nsign,iprint,nrdofH,iflag,ie,ndof
      integer :: ibeg,iend,jbeg,jend,istep,nordh,nordv,jstep
      real*8 :: weight,wa,rjac,zs
!
!-----------------------------------------------------------------------------------
!
!  ...order of approximation, orientations, geometry dof's, BC flags
      call find_order     (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
      call find_bc    (Mdle, ibc)
!
!  ...define the tensor of elestictiy
      do ii=1,3; do jj=1,3; do i=1,3; do j=1,3
        zaa(ii,jj,i,j) = MU*(del(ii,jj)*del(i,j) + del(ii,j)*del(i,jj)) +  &
                         LAMBDA*del(ii,i)*del(jj,j)
      enddo; enddo; enddo; enddo
!
!  ...clear spaces for the element matrices                    
      Bloc = ZERO
      Aloc = ZERO
!
!-----------------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L S                                          |
!-----------------------------------------------------------------------------------
!
!  ...set up the element quadrature
      etype = NODES(Mdle)%type
      call set_3Dint(etype,norder, nint,xiloc,wxi)
!
!  ...loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l); wa = wxi(l)
!
!  .....derivatives and values of the shape functions at the point
        call shape3H(etype,xi,norder,nedge_orient,nface_orient,nrdofH, &
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
          write(*,*) 'elem: NEGATIVE JACOBIAN for Mdle = ',Mdle
          stop
        endif
!
!  .....transform coefficients to master coordinates
        za = ZERO ; c = 0.d0
!
!  .....loop through displacement components
        do ii=1,3; do jj=1,3
!  .......loop through coordinates
          do i=1,3; do j=1,3
            do kk=1,3; do ll=1,3
              za(ii,jj,i,j) = za(ii,jj,i,j) +  &
                              dxidx(i,kk)*dxidx(j,ll)*zaa(ii,jj,kk,ll)
            enddo; enddo
          enddo; enddo
          if (ii.eq.jj) c(ii,jj) = - RHO*OMEGA**2
        enddo; enddo
!
!  .....total weight 
        weight = wa*rjac
!
!  .....get the source term
        call getf(Mdle,xi,x, zfval)
!
!  .....OUTER loop through dofs
        do k1=1,nrdofH
          do ivar1=1,3
            m1 = (k1-1)*3+ivar1
!
!           L O A D   V E C T O R  
            Bloc(m1,1:NR_RHS) = Bloc(m1,1:NR_RHS) + & 
                                  zfval(ivar1,1:NR_RHS)*vshapH(k1)*weight
!
!  .........INNER loop through dofs
            do k2=1,nrdofH
              do ivar2=1,3
                m2 = (k2-1)*3+ivar2
!
!               S T I F F N E S S   M A T R I X
                zs = ZERO
                do i=1,3; do j=1,3
                  zs = zs +  &
                       za(ivar1,ivar2,i,j)*dvshapH(i,k1)*dvshapH(j,k2)
                enddo; enddo
                zs = zs + c(ivar1,ivar2)*vshapH(k1)*vshapH(k2)
                Aloc(m1,m2) = Aloc(m1,m2) + zs*weight
!                
              enddo
            enddo
!  .......INNER loop            
          enddo
!  .....OUTER loop            
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
        if ((ibc(ifig,1).ne.2)) cycle
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
        do l=1,nint
          t(1:2) = tloc(1:2,l); wa = wt(l)
!
!  .......determine master element coordinates using face parameterization
          call face_param(etype,ifig,t, xi,dxidt)
!
!  .......derivatives and values of the shape functions at the point
          call shape3H(etype,xi,norder,nedge_orient,nface_orient,nrdofH,  &
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
          call getg(Mdle,ibc(ifig,1),xi,x,rn, zgval)
          if (iprint.eq.1) then
            write(*,4005) l,zgval(1:3,1)
 4005       format('elem_A: l,zgval = ',i3,2x,3(2e12.5,2x))
          endif
!
!  .......loop through test functions
          do k1=1,nrdofH
            do ivar1=1,3
              m1 = (k1-1)*3+ivar1
!
!  ...........accumulate for the load vector
                Bloc(m1,1:NR_RHS) = Bloc(m1,1:NR_RHS)   &
                          + zgval(ivar1,1:NR_RHS)*vshapH(k1)*weight
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
      iprint=0
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
