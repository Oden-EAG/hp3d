!>  Purpose:          - routine returns unconstrained (ordinary) 
!!                      stiffness matrix and load vector  
!!                      for various projection problems
!!                      (for testing shape functions routines)
!! @param[in]   Mdle  - middle node
!! @param[out]  Itest - index for assembly
!! @param[out]  Itrial- index for assembly
!-----------------------------------------------------------------------
subroutine elem(Mdle, Itest,Itrial)
  !
  use data_structure3D , only : NODES
  use parameters       , only : ZERO
  use physics          , only : NR_PHYSA
  use assembly         , only : ALOC, BLOC

  implicit none
  integer,                    intent(in)  :: Mdle
  integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial

  Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
  !  ...thermoelasticity
  call elem_ET( Mdle, &
       BLOC(1)%nrow,BLOC(2)%nrow, &
       BLOC(1)%array,BLOC(2)%array, &
       ALOC(1,1)%array,ALOC(1,2)%array, &
       ALOC(2,1)%array,ALOC(2,2)%array)

end subroutine elem
!
!> Purpose:            - routine returns unconstrained (ordinary) 
!!                       stiffness matrix and load vector  
!!                       for the thermoelasticity problem
!! @param[in]  Mdle    - middle node
!! @param[in]  NrE,NrT - dimensions of load and stiffness arrays
!! @param[out] Zbloc1, Zbloc2 - load vector
!! @param[out] Zaloc11, Zaloc12, Zaloc21, Zaloc22 - stiffness matrix
!-----------------------------------------------------------------------
subroutine elem_ET( &
     Mdle,NrE,NrT, &
     Zbloc1,Zbloc2, &
     Zaloc11,Zaloc12,zaloc21,Zaloc22)
  !
  use control
  use data_structure3D
  use element_data
  use assembly , only: NR_RHS
  use GMP

  use thermo_elast

  !-----------------------------------------------------------------------
  implicit none
  integer, intent(in)  :: Mdle
  integer, intent(in)  :: NrE, NrT

  real*8, intent(out) :: &
       Zbloc1(NrE,NR_RHS),Zbloc2(NrT,NR_RHS)
  real*8, intent(out) :: &
       Zaloc11(NrE,NrE),Zaloc12(NrE,NrT), &
       Zaloc21(NrT,NrE),Zaloc22(NrT,NrT)
  !-----------------------------------------------------------------------
  character(len=4) :: type
  !
  !  ...element, face order, geometry dof
  integer :: norder(19),nordf(5)
  real*8 :: xnod(3,MAXbrickH)
  !
  !  ...node orientations
  integer :: nedge_orient(12), nface_orient(6)
  !
  !  ...shape functions and their derivatives wrt master coordinates
  real*8 :: &
       shapH(MAXbrickH),  dshapH(3,MAXbrickH), dshapHx(3,MAXbrickH), &
       shapE(3,MAXbrickE),curlE(3,MAXbrickE), &
       shapV(3,MAXbrickV),diveV(MAXbrickV) 
  !
  !  ...geometry 
  real*8 :: xi(3),x(3),dxdxi(3,3),dxidx(3,3),t(2),dxidt(3,2),dxdt(3,2)
  !
  !  ...3D quadrature data
  real*8 :: xiloc(3,MAXbrickH),wxi(MAXbrickH)
  !
  !  ...2D quadrature data for boundary terms
  real*8 :: tloc(2,MAXquadH),wt(MAXquadH)
  !
  !  ...external unit vector
  real*8 :: rn(3)
  !
  !  ...boundary condition flags
  integer :: ibc(6,NR_PHYSA),ibcE(6),ibcT(6)
  !
  !  ...body force
  real*8 :: zf(3,MAXNRHS),zq(MAXNRHS)

  !  ...traction
  real*8 :: zgval(3,MAXNRHS),zhflx(MAXNRHS)
  !
  !  ...elasticities and other const. in master and physical coordinates
  real*8 :: zaa(3,3,3,3), za(3,3,3,3), zka(3,3), zca(3,3)
  !
  real*8 :: zalpha, zk, zla, zlam, zmu, zs

  !     Kronecker's delta
  real*8 :: del(3,3) = &
       reshape( (/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/),(/3,3/) )
  !
  integer :: &
       i,j,k,l,k1,k2,n1,n2, ii,jj,kk,ll, ie,if, &
       l1,l2, idir1,idir2, ivar1,ivar2, &
       ndof,nint,nord,nordh,nordv,nrdofH,nrdofE,nrdofV,nrdofQ, nsign, &
       ibeg,iend,jbeg,jend,iflag,iprint 

!  integer :: i,j,k,l,k1,k2,l1,l2,ii,jj,ifig,ie,if,idirl,idir2, nint,nrdofH,iflag,ixi,nsign,ibeg,iend,jbeg,jend,iprint

  real*8 ::  rjac,wa,weight,sign
  !---------------------------------------------------------------------
  !
  select case(Mdle)
  case(1)
     iprint=0
  case default
     iprint=0
  end select
  !
  if (iprint.eq.1) then
     write(*,7001) Mdle,NODES(Mdle)%type
7001 format(' elem_ET: Mdle,type = ',i10,2x,a5)
  endif
  !
  !  ...determine order of approximation
  call find_order(Mdle, norder)
  !
  !  ...determine the node orientation
  call find_orient(Mdle, nedge_orient,nface_orient)
  if (iprint.eq.1) then
     write(*,7003) Mdle,nedge_orient(1:6),nface_orient(1:4)
7003 format(' elem_ET: Mdle, nedge_orient,nface_orient = ', i3,6i2,2x,4i2)
     call pause
  endif
  !                                                                     
  !  ...determine nodes coordinates 
  call nodcor(Mdle, xnod)
  if (iprint.eq.1) then
     call celndof(NODES(Mdle)%type,Norder, nrdofH,nrdofE,nrdofV,nrdofQ)
     do k=1,nrdofH
        write(*,7002) k,xnod(1:3,k)
7002    format(' elem_ET: k,xnod(1:3,k) = ',i3,3f8.3)
     enddo
  endif
  !
  !  ...get the element boundary conditions flags
  call find_bc(Mdle, ibc)
  ibcE(1:6)=ibc(1:6,1)
  ibcT(1:6)=ibc(1:6,2)
  if (iprint.ge.1) then
     write(*,8001) 'E',Mdle,ibcE(1:nface(NODES(Mdle)%type))
     write(*,8001) 'T',Mdle,ibcT(1:nface(NODES(Mdle)%type))
8001 format('elem_ET: Mdle,ibc',a1,' = ',i8,6i2)
  endif
  !
  !  ...clear spaces for the element matrices                    
  Zbloc1 = ZERO; Zbloc2 = ZERO
  Zaloc11 = ZERO; Zaloc12 = ZERO; Zaloc21 = ZERO; Zaloc22 = ZERO 
  !
  !-----------------------------------------------------------------------
  !
  !  ...Step 1: compute element integrals
  !
  !  ...set up the element quadrature
  call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
  !
  !  ...loop through integration points
  do l=1,nint
     !      
     xi(1:3) = xiloc(1:3,l); wa = wxi(l)
     if (iprint.eq.1) then
        write(*,7004) l,xi(1:3)
7004    format(' elem_ET: l,xi = ',i3,2x,(3e12.5,2x))
     endif
     !
     !  .....compute appropriate shape functions at integration point and
     !       set up appropriate number of dof's        
     !
     !  .....evaluate H^1 shape functions at the integration point
     call shape3H( &
          NODES(Mdle)%type,xi,norder, &
          nedge_orient,nface_orient, &
          nrdofH,shapH,dshapH)
     !
     !  .....determine physical coordinates and the derivatives of 
     !       the physical coordinates wrt master element coordinates
     x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
     do k=1,nrdofH
        x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
        do i=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
        enddo
     enddo
     !
     !  .....evaluate the inverse derivatives and jacobian
     iflag = 0
     call geom(dxdxi, dxidx,rjac,iflag) 
     if (iflag.ne.0) then
        write(*,*) 'elem_ET: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
        write(*,*) '        rjac = ',rjac
        call geom(dxdxi, dxidx,rjac,iflag)
        call result
        call print_GMP
     endif
     !
     !  .....total weight
     weight = wa*rjac
     !
     !  .....get source terms
     call getf(Mdle,x, zf,zq)
     !
     !  .....get material data...........................................
     call getmat(Mdle,x, zmu,zla,zk,zalpha)
     !
     !  .....define the elasticities
     do ii=1,3; do jj=1,3; do i=1,3; do j=1,3
        zaa(ii,jj,i,j) =  &
             zmu*(del(ii,jj)*del(i,j) + del(ii,j)*del(i,jj)) + zla*del(ii,i)*del(jj,j)
     enddo; enddo; enddo; enddo
     !
     !  .....transform elasticity coefficients to master coordinates
     za = ZERO
     !
     !  .....loop through displacement components
     do ii=1,3; do jj=1,3
        !
        !  .......loop through coordinates
        do i=1,3; do j=1,3
           do kk=1,3; do ll=1,3
              za(ii,jj,i,j) = za(ii,jj,i,j) + dxidx(i,kk)*dxidx(j,ll)*zaa(ii,jj,kk,ll)
           enddo; enddo
        enddo; enddo
     enddo; enddo
     !
     !  .....transform thermal conductivity to master coordinates....
     zka = ZERO
     !  .....loop through coordinates
     do kk=1,3; do ll=1,3
        !  ........summation loop
        do ii=1,3
           zka(kk,ll)=zka(kk,ll)+zk*dxidx(kk,ii)*dxidx(ll,ii)
        enddo
     enddo; enddo
     !
     !  .....transform "c=alpha*(2*mu+3*lam)" to master coordinates...
     zca = ZERO
     !  .....loop through coordinates
     do kk=1,3; do ii=1,3           
        zca(kk,ii)=zalpha*(2*zmu+3*zlam)*dxidx(kk,ii)
     enddo; enddo           

     !
     !---------------------------------------------------------------- 11 - 12
     !
     !  .....elasticity blocks
     !
     !  .....loop through test functions
     do k1=1,nrdofH
        !
        !  .......loop through displacement components
        do ivar1=1,3
           !
           !  .........elasticty test function #
           n1 = (k1-1)*3+ivar1
           !
           !  .........accumulate for the load vector
           Zbloc1(n1,1:Nr_RHS) = Zbloc1(n1,1:Nr_RHS) &
                + zf(ivar1,1:Nr_RHS)*shapH(k1)*weight
           !
           !  .........loop through displacement trial functions
           do k2=1,nrdofH
              !
              !  ...........loop through components
              do ivar2=1,3
                 !
                 !  .............elasticity trial function #
                 n2 = (k2-1)*3+ivar2
                 !
                 !  .............accumulate for matrix Zaloc11
                 zs = ZERO
                 do i=1,3; do j=1,3
                    zs = zs + za(ivar1,ivar2,i,j)*dshapH(i,k1)*dshapH(j,k2)
                 enddo; enddo
                 Zaloc11(n1,n2) = Zaloc11(n1,n2)+zs*weight
              enddo
           enddo
           !
           !  .........loop through temperature trial functions
           do k2=1,nrdofH
              !
              !  ............accumulate for matrix Zaloc12
              do kk=1,3
                 Zaloc12(n1,k2) = Zaloc12(n1,k2) &
                      - zca(kk,ivar1)*dshapH(kk,k1)*shapH(k2)*weight
              enddo

           enddo
        enddo
     enddo
     !
     !------------------------------------------------------------------ 21- 22
     !
     !  .....heat blocks
     !
     !  .....loop through temperature test functions
     do k1=1,nrdofH
        !
        !  .......accumulate for the load vector
        Zbloc2(k1,1:Nr_RHS) = Zbloc2(k1,1:Nr_RHS) + zq(1:Nr_RHS)*shapH(k1)*weight
        !
        !  .......loop through displacement trial functions
        do k2=1,nrdofH
           !
           !  .........loop through components
           do ivar2=1,3
              !
              !  ...........elasticity trial function #
              n2 = (k2-1)*3+ivar2
              !
              !  ...........accumulate for matrix Zaloc21
              Zaloc21(k1,n2) = Zaloc21(k1,n2) + 0.d0
           enddo
        enddo
        !
        !  .......loop through temperature trial functions
        do k2=1,nrdofH
           !
           !  .........accumulate for matrix Zaloc22
           do idir1=1,3; do idir2=1,3
              Zaloc22(k1,k2) = Zaloc22(k1,k2) &
                   + zka(idir1,idir2)*dshapH(idir1,k1)*dshapH(idir2,k2)*weight
           enddo; enddo

        enddo
     enddo
     !  ...end of loop through integration points
  enddo
  !
  !======================================================================

  !
  !  ...Step 2: compute boundary integrals
  !
  !  ...loop through element faces
  do if=1,nface(NODES(Mdle)%type)
     !
     if (ibcE(if).ne.2 .OR. ibcT(if).ne.2) cycle
     !     
     !  .....determine order for the face
     call face_order(NODES(Mdle)%type,if,norder, nordf)
     !
     !  .....set up the face quadrature
     type = face_type(NODES(Mdle)%type,if)
     call set_2Dint(type,nordf, nint,tloc,wt)
     if (iprint.eq.1) then
        write(*,4001) if,ibc(if,1),nint,nordf
4001    format('elem_ET: if,ibc(if,1),nint,nordf = ', i2,2x,i2,2x,i3,2x,4i2,i3)
     endif
     nsign = Nsign_param(NODES(Mdle)%type,If)
     !
     !  .....loop through face integration points
     do l=1,nint
        t(1:2) = tloc(1:2,l); wa = wt(l)
        !
        !  .......use the face parametrization to determine the master
        !         element coordinates
        call face_param(NODES(Mdle)%type,if,t, xi,dxidt)
        !
        !  .......evaluate derivatives and values of the shape functions
        !         at the point
        call shape3H(NODES(Mdle)%type,xi,norder, &
             nedge_orient,nface_orient, &
             nrdofH,shapH,dshapH)
        !
        !  .......determine physical coordinates and the derivatives of 
        !         the physical coordinates wrt master element coordinates
        x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
        do k=1,nrdofH
           x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
           do i=1,3
              dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
           enddo
        enddo
        if (iprint.eq.1) then
           write(*,4002) if,(dxdxi(i,1:3),i=1,3)
4002       format('elem_ET: if,dxdxi = ',i2,3(/,3f8.3))
           write(*,4003) (dxidt(i,1:2),i=1,3)
4003       format('elem_ET: dxidt = ',3(/,2f8.3))
        endif
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
        if (iprint.eq.2) then
           write(*,4004) if,sign,rn
4004       format('elem_ET: if,nsign,rn = ',i2,2x,i2,2x,3f8.3)
           call pause
        endif
        !
        !  .......get the Neumann data
        call getg(Mdle,ibc(if,1),xi,x,rn, zgval,zhflx)
        if (iprint.eq.1) then
           write(*,4005) l,zgval(1:3,1)
4005       format('elem_ET: l,zgval = ',i3,2x,3(2e12.5,2x))
        endif
        !
        ! .......loop through test functions
        do k1=1,nrdofH
           do ivar1=1,3
              n1 = (k1-1)*3+ivar1
              !
              ! ............accumulate for the load vectors
              Zbloc1(n1,1:Nr_RHS) = Zbloc1(n1,1:Nr_RHS) &
                   + zgval(ivar1,1:Nr_RHS)*shapH(k1)*weight
           enddo
           Zbloc2(k1,1:Nr_RHS) = Zbloc2(k1,1:Nr_RHS) &
                + zhflx(1:Nr_RHS)*shapH(k1)*weight

        enddo

        !
        !  .....end of loop through integration points
     enddo
     !
     !  ...end of loop through faces
  enddo

  !
  !
  !=====================================================================
  !
  if (iprint.eq.1) then
     !
     select case(NODES(Mdle)%type)
     case('mdlp')
        write(*,7011) 1,2,3,4,5,6
7011    format('elem_ET: VERTEX SHAPE FUNCTIONS = ',8i2)
        k=6
        do ie=1,9
           ndof = norder(ie)-1
           if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
7012       format('        EDGE ',i2,' SHAPE FUNCTIONS = ',i3,'...',i3)
           k = k+ndof
        enddo
        do if=1,5
           nord = norder(9+if)
           select case(if)
           case(1,2); ndof = (nord-2)*(nord-1)/2
           case(3,4,5); call decode(nord, nordh,nordv)
              ndof=(nordh-1)*(nordv-1)
           end select
           if (k+ndof.ge.k+1) write(*,7013) if,k+1,k+ndof
7013       format('        FACE ',i2,' SHAPE FUNCTIONS = ',i3,'...',i3)
           k=k+ndof
        enddo
        if (nrdofH.ge.k+1) write(*,7014) k+1,nrdofH
7014    format('        MIDDLE NODE SHAPE FUNCTIONS = ',i3,'...',i3)
     case('mdln')
        write(*,7011) 1,2,3,4
        k=4
        do ie=1,6
           ndof = norder(ie)-1
           if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
           k = k+ndof
        enddo
        do if=1,4
           nord = norder(6+if)
           ndof = (nord-2)*(nord-1)/2
           if (k+ndof.ge.k+1) write(*,7013) if,k+1,k+ndof
           k=k+ndof
        enddo
        if (nrdofH.ge.k+1) write(*,7014) k+1,nrdofH
     case('mdld')
        write(*,7011) 1,2,3,4,5
        k=5
        do ie=1,8
           ndof = norder(ie)-1
           if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
           k = k+ndof
        enddo
        call decode(norder(9), nordh,nordv)
        ndof = (nordh-1)*(nordv-1)
        if (k+ndof.ge.k+1) write(*,7013) 1,k+1,k+ndof
        k=k+ndof
        do if=2,5
           nord = norder(8+if)
           ndof = (nord-2)*(nord-1)/2
           if (k+ndof.ge.k+1) write(*,7013) if,k+1,k+ndof
           k=k+ndof
        enddo
        if (nrdofH.ge.k+1) write(*,7014) k+1,nrdofH
     case('mdlb')
        write(*,7011) 1,2,3,4,5,6,7,8
        k=8
        do ie=1,12
           ndof = norder(ie)-1
           if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
           k = k+ndof
        enddo
        do if=1,6
           nord = norder(12+if)
           call decode(nord, nordh,nordv)
           ndof = (nordh-1)*(nordv-1)
           if (k+ndof.ge.k+1) write(*,7013) if,k+1,k+ndof
           k=k+ndof
        enddo
        if (nrdofH.ge.k+1) write(*,7014) k+1,nrdofH
     end select
     !
50   write(*,*) 'elem_ET: Zbloc1 = '
     write(*,*) 'elem_ET: SET jbeg,jend,l1,l2 (jbeg=0 to exit)'
     read(*,*) jbeg,jend,l1,l2
     if (jbeg.eq.0) go to 55
     do j=jbeg,jend
        write(*,7005) j,Zbloc1(j,l1:l2)
     enddo
     go to 50
55   write(*,*) 'elem_ET: Zaloc = '
     write(*,*) 'elem_ET: SET jbeg,jend,ibeg,iend (jbeg=0 to exit)'
     read(*,*) jbeg,jend,ibeg,iend
     if (jbeg.eq.0) go to 60
     do j=jbeg,jend
        write(*,7005) j,Zaloc11(j,ibeg:iend)
#if C_MODE
7005    format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
#else
7005    format(i3,2x,10(e12.5,2x),10(/,5x,10e12.5))
#endif
     enddo
     go to 55
60   continue
  endif
  !
  !
end subroutine elem_ET

