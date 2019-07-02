!-------------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for EM problem
!! @param[in]  Mdle      - element (middle node) number 
!! @param[out] Bloc      - element load vector(s)
!! @param[in]  Nrow_Bloc - number of rows of load vector
!! @param[out] Aloc      - element stiffness matrix
!! @param[in]  Ncol_Aloc - number of columns of stiffness matrix
!-------------------------------------------------------------------------------------------
!
subroutine elem_em(Mdle,Bloc,Nrow_Bloc,Aloc,Ncol_Aloc) 
  use control  , only : EXGEOM
  use assembly , only : NR_RHS
  use data_structure3D
  use em
!-------------------------------------------------------------------------------------------
  implicit none
  integer,                                   intent(in)  :: Mdle
  integer,                                   intent(in)  :: Nrow_Bloc,Ncol_Aloc
  complex*16,dimension(Nrow_Bloc,NR_RHS),    intent(out) :: Bloc
  complex*16,dimension(Nrow_Bloc,Ncol_Aloc), intent(out) :: Aloc
!-------------------------------------------------------------------------------------------
! element and face type 
  character(len=4) :: etype,ftype
!     
! element order, face order, edge and face orientations
  integer,dimension(19) :: norder
  integer,dimension(5)  :: nordf
  integer,dimension(12) :: nedge_orient
  integer,dimension(6)  :: nface_orient
!
! BC flagx
  integer,dimension(6,NR_PHYSA) :: ibc
!
! shape functions and their derivatives wrt master coordinates
  real*8,dimension(3,MAXbrickE) ::  shapE,  curlE
  real*8,dimension(3,MAXbrickE) ::  shapEx, curlEx
  real*8,dimension(  MAXbrickE) ::  shapEr, curlEr
  real*8,dimension(  MAXbrickH) ::  shapH
  real*8,dimension(3,MAXbrickH) :: dshapH
!
! geometry 
  real*8,dimension(3,MAXbrickH) :: xnod
  real*8,dimension(3)           :: xi,x
  real*8,dimension(3,3)         :: dxdxi,dxidx
  real*8,dimension(2)           :: t
  real*8,dimension(3,2)         :: dxidt,dxdt
!
! 2D and 3D quadrature
  real*8,dimension(3,MAXquadH ) :: tloc
  real*8,dimension(  MAXquadH ) :: wt
  real*8,dimension(3,MAXbrickH) :: xiloc
  real*8,dimension(  MAXbrickH) :: wxi
!
! impressed current and impressed surface current
  complex*16,dimension(3,NR_RHS) :: j_imp,j_imp_s
!
! external unit vector
  real*8,dimension(3) :: rn
  !
  real*8                       :: wa, weight, rjac
  !
  integer :: i,j,k,l,n,k1,k2,nint,iprint,iflag,ifig,nsign
  integer :: nrdofE,nrdofH
  !
  complex*16                       :: zkappa 
  complex*16                       :: z_tmp_a, z_tmp_b
  !
  integer :: l1, l2, ibeg, iend, jbeg, jend
!-------------------------------------------------------------------------------------------
!
  iprint=0
!
! determine order of approximation, orientations, nodes coordinates
  call find_order     (Mdle, norder)
  call find_orient(Mdle, nedge_orient,nface_orient)
  call nodcor     (Mdle, xnod)
  call find_bc    (Mdle,ibc) 
!
! clear spaces for stiffness matrix and load vector
  Aloc = ZERO ; Bloc = ZERO
!
!-------------------------------------------------------------------------------------------
! E L E M E N T    I N T E G R A L S                                                       !
!-------------------------------------------------------------------------------------------
!
! collect integration points
  etype = NODES(Mdle)%type
  call set_3Dint(etype,norder, nint,xiloc,wxi)
!
! loop over integration points
  do l=1,nint
     xi(1:3) = xiloc(1:3,l) ; wa = wxi(l)
!
!  ..calculate shape function for Hcurl and H1
     call shape3E(etype,xi,norder,nedge_orient,nface_orient, nrdofE,shapE,curlE)
     call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,dshapH)
!
!  ..geometry map
     select case(EXGEOM)
     case(0)
       x(1:3)= 0.d0 ; dxdxi(1:3,1:3) = 0.d0
       do k=1,nrdofH
          x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
          do i=1,3
             dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
          enddo
       enddo
     case(1)
       call exact_geom(Mdle,Xi, x,dxdxi)
     endselect
!
!  ..evaluate the inverse derivatives and jacobian
     call geom(dxdxi, dxidx,rjac,iflag)
     if (iflag.ne.0) then
        write(*,*) 'elem_em: NEGATIVE JACOBIAN, Mdle, rjac= ',Mdle, rjac
        stop 1
     endif
!     
!  ..total weight
     weight = wa*rjac
!
     zkappa = EPSIL*(omega**2) - Z_i*OMEGA*SIGMA
!
!  ..evaluate the shapE and curlE wrt physical coordinate (Piola transform)
     shapEx = 0.d0 ; curlEx = 0.d0
     do k=1,nrdofE
       do i=1,3
         do j=1,3
           shapEx(i,k) = shapEx(i,k) + (dxidx(j,i)*shapE(j,k))
           curlEx(i,k) = curlEx(i,k) + (dxdxi(i,j)*curlE(j,k))/rjac
         enddo
       enddo 
     enddo 
!
!  ..get impressed current
     call getf(Mdle,Xi,x, j_imp)

!  ..loop over TEST functions
     do k1=1,nrdofE
!  .....loop over equations     
        do n=1,NR_RHS
           z_tmp_a = Z_0
           do i=1,3
              z_tmp_a = z_tmp_a + Z_i*OMEGA*j_imp(i,n)*shapEx(i,k1)
           enddo
           
!          L O A D    V E C T O R           
           Bloc(k1,n) = Bloc(k1,n) - (z_tmp_a*weight)           
!
!  .....loop over equations
        enddo              
!        
!  .....loop over TRIAL functions
        do k2=1,nrdofE
           z_tmp_a = Z_0
           z_tmp_b = Z_0
           do i=1,3
              z_tmp_a = z_tmp_a + curlEx(i,k1)*curlEx(i,k2)
              z_tmp_b = z_tmp_b + shapEx(i,k1)*shapEx(i,k2)
           enddo
!           
!          S T I F F N E S S    M A T R I X                                    
           Aloc(k1,k2) = Aloc(k1,k2) + &                                       
                ( 1.d0/MU*z_tmp_a - zkappa*z_tmp_b )*weight                    
!                
!  .....loop over TRIAL functions
        enddo
!        
!  ..loop over TEST functions        
     enddo
!     
! loop over integration points     
  enddo 
!
!-------------------------------------------------------------------------------------------
! B O U N D A R Y    I N T E G R A L                                                       !
!-------------------------------------------------------------------------------------------
!
! loop over element faces
  do ifig=1,nface(etype)
!
! ..skip of not a Neumann interface
    if (ibc(ifig,1).ne.2) cycle
!
! ..determine face order
    call face_order(etype,ifig, norder,nordf)  
!
! ..set up face quadrature
    ftype = face_type(etype,ifig)
    call set_2Dint(ftype,nordf, nint,tloc,wt)
    nsign = nsign_param(etype,ifig)
!
! ..loop over integration points
    do l=1,nint
!    
       t(1:2) = xiloc(1:2,l) ; wa = wxi(l)
       call face_param(etype,ifig,t, xi,dxidt)
!       
!  ....calculate shape function for Hcurl and H1
       call shape3E(etype,xi,norder,nedge_orient,nface_orient, nrdofE,shapE,curlE)
       call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,dshapH)
!
!  ....geometry map
       select case(EXGEOM)
       case(0)
         x(1:3) = 0.d0 ; dxdxi(1:3,1:3) = 0.d0
         do k=1,nrdofH
           x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
           do i=1,3
             dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
           enddo
         enddo
       case(1)
         call exact_geom(Mdle,xi, x,dxdxi)
       end select
!       
!  ....normal unit vector and jacobian
       dxdt(1:3,1:2) = 0.d0
       do i=1,2
         do j=1,3
           dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
         enddo
       enddo
       call cross_product(dxdt(1:3,1),dxdt(1:3,2), rn)
       call norm(rn, rjac)
       rn(1:3) = rn(1:3)*nsign/rjac
!
!  ....total weight
       weight = wa*rjac
!
       zkappa = EPSIL*(OMEGA**2) - Z_i*OMEGA*SIGMA
!
!  ....evaluate the shapE and curlE wrt physical coordinate (Piola transform)
       shapEx = 0.d0 ; curlEx = 0.d0
       do k=1,nrdofE
         do i=1,3
           do j=1,3
             shapEx(i,k) = shapEx(i,k) + (dxidx(j,i)*shapE(j,k))
             curlEx(i,k) = curlEx(i,k) + (dxdxi(i,j)*curlE(j,k))/rjac
           enddo
         enddo
       enddo
!       
!  ....get impressed surface current
       call getg(Mdle,ibc(ifig,1),Xi,x,rn, j_imp_s)

!  ....loop OVER test functions
       do k1=1,nrdofE
         do n=1,NR_RHS
           z_tmp_a=Z_0
           do i=1,3
             z_tmp_a = z_tmp_a + Z_i*OMEGA*j_imp_s(i,n)*shapE(i,k1)
           enddo
!
!          L O A D    V E C T O R
           Bloc(k1,n) = Bloc(k1,n) - (z_tmp_a*weight)           
!           
         enddo
!  ....loop over TEST functions        
       enddo
!
!  ..loop over integration points
     enddo 
!
! loop over faces
  enddo







  !
  iprint = 0
  !
  if (iprint.eq.-1) then
     write(*,*) 'elem_em: Mdle = ', Mdle
50   write(*,*) 'elem_em: Zbloc, nrdofE = ', nrdofE
     write(*,*) 'elem_em: SET jbeg,jend,l1,l2(jbeg=0 to exit)'
     read(*,*) jbeg,jend,l1,l2
     if (jbeg.eq.0) go to 55
     do j=jbeg,jend
        write(*,7005) j,Bloc(j,l1:l2)
     enddo
     go to 50
55   write(*,*) 'elem_em: Zaloc, nrdofE = ', nrdofE
     write(*,*) 'elem_em: SET jbeg,jend,ibeg,iend(jbeg=0 to exit)'
     read(*,*) jbeg,jend,ibeg,iend
     if (jbeg.eq.0) go to 60
     do j=jbeg,jend
        write(*,7005) j,Aloc(j,ibeg:iend)
7005    format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
     enddo
     go to 55
60   continue
  endif
  !
end subroutine elem_em
