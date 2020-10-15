subroutine generate_visual_poly_elem_hni
  use data_structure3D_poly
  use parametersDPG
  use connectivity_poly
  use physics
  use common_prob_data
  use geometry_polydpg, only: ALLC_DOM
!
  implicit none
  
  integer :: Mdle
  real*8,dimension(12)  :: errorQ
  
!
!     node case (decimal form)
  integer,dimension(NR_PHYSA) :: icased
!
!     element, face order, geometry dof
  integer,dimension(19)          :: norder
  real*8 ,dimension(3,MAXbrickH) :: xnod
  integer,dimension(12)          :: nedge_orient
  integer,dimension(6)           :: nface_orient
!
!     geometry
  real*8,dimension(3)   :: xi,x
  real*8,dimension(2)   :: t
  real*8,dimension(3,3) :: dxidx,dxdxi
  real*8                :: rjac
!
!     3D quadrature data
  real*8, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!     approximate solution
  real*8, dimension(  MAXEQNH  ) ::  zsolH
  real*8, dimension(  MAXEQNH,3) :: zdsolH
  real*8, dimension(3,MAXEQNE  ) ::  zsolE
  real*8, dimension(3,MAXEQNE  ) :: zcurlE
  real*8, dimension(3,MAXEQNV  ) ::  zsolV
  real*8, dimension(  MAXEQNV  ) ::  zdivV
  real*8, dimension(  MAXEQNQ  ) ::  zsolQ
!
!     exact solution
  real*8, dimension(  MAXEQNH    ) ::   zvalH
  real*8, dimension(  MAXEQNH,3  ) ::  zdvalH
  real*8, dimension(  MAXEQNH,3,3) :: zd2valH
  real*8, dimension(3,MAXEQNE    ) ::   zvalE
  real*8, dimension(3,MAXEQNE,3  ) ::  zdvalE
  real*8, dimension(3,MAXEQNE,3,3) :: zd2valE
  real*8, dimension(3,MAXEQNE    ) ::  zcurlE_ex
  real*8, dimension(3,MAXEQNV    ) ::   zvalV
  real*8, dimension(3,MAXEQNV,3  ) ::  zdvalV
  real*8, dimension(3,MAXEQNV,3,3) :: zd2valV
  real*8, dimension(  MAXEQNV    ) ::   zdivV_ex
  real*8, dimension(  MAXEQNU    ) ::   zvalU
  real*8, dimension(  MAXEQNF    ) ::   zvalF
  real*8, dimension(  MAXEQNQ    ) ::   zvalQ
  real*8, dimension(  MAXEQNQ,3  ) ::  zdvalQ
  real*8, dimension(  MAXEQNQ,3,3) :: zd2valQ
  real*8,dimension(  MAXbrickQ) :: shapQ


  real*8, dimension(  MAXtetraHH) :: shapHH
  real*8, dimension(3,MAXtetraHH) :: gradHH
  real*8, dimension(3,MAXtetraVV) :: shapVV
  real*8, dimension(  MAXtetraVV) :: divVV

!
      integer, dimension(MAXtetraHH) ::degH
      integer, dimension(MAXtetraVV) ::degV
      integer, dimension(MAXtetraQ ) ::degQ
!
!     miscellanea
  integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag,k,nordp, idom
  
  real*8  :: weight,wa
!
!     printing flag
  integer :: iprint
  character(len=4) :: etype,ftype
  integer :: ndofE,ndofH,ndofV,ndofQ,nrf,nrv,nre,nrdofQ,nrv_f,ntri,ntri_total, &
              npt, npt_total,jf,mdlf,dump_shape,dump_exa,dump_num,dump_err,dump_e3d,jv, &
              nrdofHH,nrdofVV
! 
  integer, allocatable :: nfaces(:),Norientf(:),Nedges(:),Nverts(:),nverts_f(:),trian(:,:)
  real*8, allocatable :: xverts(:,:),affverts_f(:,:),xverts_f(:,:),verts(:,:),xiloc(:,:),xloc(:,:)
!
  real*8  :: R_e,r_f,Area_f
  real*8, dimension(3) :: Xg_e,xg_f,fn,a0,a1,a2,a3,b0,b1,b2,aux
  real*8, dimension(3,4) :: a_aff
  real*8, dimension(3,3) :: b_aff
  real*8, dimension(3,2) :: dxi_eta
!
!---------------------------------------------------------------------------------------
!
  iprint=0
!   
  dump_exa=155
    open(unit=dump_exa,file='output/elem_exa', &
      form='formatted',access='sequential',status='unknown')
  dump_num=156
    open(unit=dump_num,file='output/elem_num', &
      form='formatted',access='sequential',status='unknown')
  dump_err=157
    open(unit=dump_err,file='output/elem_err', &
      form='formatted',access='sequential',status='unknown')
  dump_shape=158
    open(unit=dump_shape,file='output/elem_shape', &
      form='formatted',access='sequential',status='unknown')
  dump_e3d=159
      open(unit=dump_e3d,file='output/elem_3d', &
        form='formatted',access='sequential',status='unknown')
!
  ntri_total = 0
  npt_total = 0
  do mdle = 1,NRELES
!
  errorQ=0.d0 ; aux = 0.d0
!
  call find_ndof_poly(mdle, ndofH,ndofE,ndofV,ndofQ)
! 
  zdofQ = 0.d0
  zdofQ(1:NRQVAR,1:ndofQ) = NODES(mdle)%zdofQ(1:NRQVAR,1:ndofQ)
! 
  etype = 'tetr'
! 
!     retrieve element node data
  call elem_nodes_poly(Mdle,Nfaces,Norientf,Nrf,Nedges,Nre,Nverts,Nrv)
! !     allocate array for vertices 
  allocate(xverts(3,Nrv))
!     get vertices coordinates
  do jv = 1, nrv
    xverts(1:3,jv) = NODES(nverts(jv))%coord(1:3,1)
  enddo
!
!     get element affine coordinates points a0,a1,a2,a3
  call poly_affine_3d(xverts,nrv,a0,a1,a2,a3)
!     construct array a_aff with a0,a1,a2,a3 as columns
  a_aff(:,1) = a0; a_aff(:,2) = a1; a_aff(:,3) = a2; a_aff(:,4) = a3
! ! computational of jacobian determinant
! ! get normal of face a1a2a3
!       call face_normal(a_aff(:,2:4),fn_a)
! ! get tetrahedron volume
!       call tetra_centr_vol(a0,a1,a2,a3,fn_a,x_aux,a_vol)
! ! correct
!       rjac_e = 6*a_vol
! construct linear transformation associated to affine map. Interpreted as classic jacobian
  dxdxi(:,1) = a1 - a0
  dxdxi(:,2) = a2 - a0
  dxdxi(:,3) = a3 - a0
! get inverse and determinant of linear transformation
  call geom(dxdxi,dxidx,rjac,iflag)
! 
  if (iflag.ne.0) then
    write(*,8999) Mdle,rjac
 8999   format('Negative Jacobian!Mdle,rjac=',i8,2x,e12.5)
    stop
  endif
! 
!     get order, in the norder structure of a uniform tetrahedral element
  norder = 0
  norder(1:15) = NODES(Mdle)%order   
! get enriched test functions nominal order
  nordP = NODES(Mdle)%order + NORD_ADD
!     set up the element quadrature

  idom = ALLC_DOM(mdle)

  do jf=1,Nrf

    mdlf = Nfaces(jf)

    call face_vert_list(mdlf,nrv_f,nverts_f)
    allocate(xverts_f(3,nrv_f))
    do jv=1,Nrv_f
      xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)
    enddo
  ! get face affine coordinates points b0,b1,b2
    call poly_affine_2d(xverts_f,nrv_f,b0,b1,b2)
    b_aff(:,1) = b0
    b_aff(:,2) = b1
    b_aff(:,3) = b2
!       get face normal in physical coordinates
    call face_normal(xverts_f(:,1:3),fn)


    allocate(affverts_f(2,Nrv_f))
    affverts_f = 0.d0
    call poly_x_to_affine_2d(xverts_f,nrv_f,b_aff,affverts_f)

!   allocate memory for evaluation points on current face
    allocate(verts(2,(nrv_f-2)*12+nrv_f),trian(3,16*(nrv_f-2)))

    call break_poly2d(mdlf,NODES(Mdlf)%order - 1,nrv_f,affverts_f,ntri,npt,trian,verts)

    do l=1,ntri   
      write(dump_e3d,*) mdle,",",npt_total+trian(1,l),",", &
                                 npt_total+trian(2,l),",", &
                                 npt_total+trian(3,l),","
    enddo 


    allocate(xiloc(3,npt))
!   transform evaluation points' face affine coordinates to element affine coordinates
    xiloc = 0.d0
    dxi_eta = 0.d0
    call poly_affine_2d_to_3d(verts,npt,b_aff,a_aff,xiloc,dxi_eta)    
!
!   get physical points from face affine coordinates
    allocate(xloc(3,npt))
    xloc = 0.d0
    call poly_affine_2d_to_x(verts,npt,b_aff,xloc)
!
!         loop through evaluation points
    do l=1,npt
!     retrieve evaluation point in face affine coordinates t (eta)
      t(1:2)=verts(1:2,l)
!     retrieve it in element affine coordinates xi
      xi(1:3) = xiloc(1:3,l)
!     retrieve it in physical coordinates x 
      x(1:3) = xloc(1:3,l)
!
!                  -- APPROXIMATE SOLUTION --
!----------------------COPIED FROM SOLEVAL---------------------------
!     L2 shape functions
      ! call shape3Q(etype,Xi,Norder, nrdofQ,shapQ)

      ! call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)

      ! call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!        .....Compute shape functions needed for test/trial field variables and geometry
!             L2 (field trial)
      call shap3Q_mono(xi,NODES(Mdle)%order, nrdofQ,shapQ,degQ)              
!             H1 (test)
      call shap3HH_mono(xi,nordP, nrdofHH,shapHH,gradHH,degH)              
!             H(div) (test)
      call shap3VV_mono(xi,nordP, nrdofVV,shapVV,divVV,degV)
!      
!     Piola transform
      shapQ(1:nrdofQ)=shapQ(1:nrdofQ)/rjac
!      
!     evaluate the approximate solution
      ZsolQ(1:MAXEQNQ)=ZERO
      do k=1,nrdofQ
        do ivar=1,MAXEQNQ
          ZsolQ(ivar) = ZsolQ(ivar) + zdofQ(ivar,k)*shapQ(k)
        enddo
      enddo
!--------------------------------------------------------------------
!
!           -- EXACT SOLUTION --
      call exact(x,aux,0, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV,zvalU, zvalF, zvalQ,zdvalQ,zd2valQ)
!--------------------------------------------------------------------
!       POINTWISE ERROR        
      ErrorQ(1:NRQVAR) = (zvalQ(1:NRQVAR) - zsolQ(1:NRQVAR))

!--------------------------------------------------------------------

      write(dump_exa,9877) mdle,npt_total+l,x(1:3),zvalQ(1:12)
      write(dump_num,9877) mdle,npt_total+l,x(1:3),zsolQ(1:12)
      write(dump_err,9877) mdle,npt_total+l,x(1:3),errorQ(1:12)
 9877   format(2(i4.2,","),14(e22.15,","),e22.15)


      write(dump_shape,9875) shapQ(1),shapQ(2),shapQ(8), &
                             shapHH(4),shapHH(6),shapHH(11), &
                             shapVV(1:3,1),shapVV(1:3,14)
 9875 format(11(e22.15,","),e22.15)
!
!         loop over evaluaion points
    enddo

    deallocate(affverts_f,xverts_f,trian,verts,xiloc,xloc)

    ntri_total = ntri_total + ntri
    npt_total  = npt_total  + npt      

!   end loop over faces
  enddo

  deallocate(xverts,nfaces,norientf,nedges,nverts)
! end loop over elements
  enddo

  write(*,*) "***Element visualization: Dumped data for ",ntri_total," triangles and ",npt_total," points"

  close(dump_num)
  close(dump_exa)
  close(dump_err)
  close(dump_shape)
  close(dump_e3d)
end subroutine generate_visual_poly_elem_hni

subroutine generate_visual_poly_face_ct

    use control          , only : INTEGRATION
    use data_structure3D_poly
    use geometry_polydpg
    use parametersDPG
    use connectivity_poly
    use environment      , only : L2PROJ
    use physics
!
    implicit none
    integer, dimension(NR_PHYSA) :: Flag
    integer                      :: Mdlf
    real*8, dimension(NRFVAR)    :: errorF
    real*8, dimension(NRUVAR)    :: erroru
    real*8                       :: rnormF,rnormU
!
!     node case (decimal form)
    integer,dimension(NR_PHYSA) :: icased
!
!     element, face order, geometry dof
    integer,dimension(5)           :: nordf
    real*8 ,dimension(3,MAXbrickH) :: xnod
    integer,dimension(12)          :: nedge_orient
    integer,dimension(6)           :: nface_orient
!
!     geometry
    real*8,dimension(3)   :: x
    real*8,dimension(2)   :: t
    ! real*8,dimension(3,3) :: dxidx,dxdxi
    real*8                :: rjac_f
!
!     3D quadrature data
    real*8,allocatable :: tloc(:,:)
    real*8,allocatable :: wt(:)
!
!     approximate solution dof's
    real*8, dimension(MAXEQNU,MAXtriaH) :: zdofU
    real*8, dimension(MAXEQNF,MAXtriaQ) :: zdofF
!
!     approximate solution
    real*8, dimension(  MAXEQNF  )   ::  zsolF
    real*8, dimension(  MAXEQNU  )   ::  zsolU
!
!     exact solution
    real*8, dimension(  MAXEQNU    ) ::   zvalU
    real*8, dimension(  MAXEQNF    ) ::   zvalF
    !     exact solution
    real*8, dimension(  MAXEQNH    ) ::   zvalH
    real*8, dimension(  MAXEQNH,3  ) ::  zdvalH
    real*8, dimension(  MAXEQNH,3,3) :: zd2valH
    real*8, dimension(3,MAXEQNE    ) ::   zvalE
    real*8, dimension(3,MAXEQNE,3  ) ::  zdvalE
    real*8, dimension(3,MAXEQNE,3,3) :: zd2valE
    real*8, dimension(3,MAXEQNE    ) ::  zcurlE_ex
    real*8, dimension(3,MAXEQNV    ) ::   zvalV
    real*8, dimension(3,MAXEQNV,3  ) ::  zdvalV
    real*8, dimension(3,MAXEQNV,3,3) :: zd2valV
    real*8, dimension(  MAXEQNV    ) ::   zdivV_ex
    real*8, dimension(  MAXEQNQ    ) ::   zvalQ
    real*8, dimension(  MAXEQNQ,3  ) ::  zdvalQ
    real*8, dimension(  MAXEQNQ,3,3) :: zd2valQ
    real*8, dimension(  MAXtriaQ)    :: shapQ_f
    real*8,dimension(  MAXtriaH)    :: shapH_f
    real*8,dimension(2,MAXtriaH)    :: gradH_f
    real*8,dimension(6) :: sigma
!
!     miscellanea
    integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag,k
    
    real*8  :: weight,wa
!
!     printing flag
    integer :: iprint

    integer :: ndofE,ndofH,ndofV,ndofQ,nrv_f,nrdofQ_f , &
     ndofF, nrdofF,jv,jf,dump_exa,dump_num,dump_err,dump_f3d,npt,npt_total,ntri,ntri_total,original, &
     nord_mdlf,ndofU,nrdofH_f, nf_act, idom, iblock, mdle, ires
    integer, dimension(4)  :: norie
    integer, allocatable :: Nverts_f(:)
    real*8, allocatable :: Xverts_f(:,:),affverts_f(:,:),xloc(:,:)
    real*8, allocatable :: verts(:,:)
    integer,allocatable :: trian(:,:)
    real*8  :: b_area
    real*8, dimension(3) :: Xg_f,fn,x_aux,b0,b1,b2
    real*8, dimension(3,3) :: b_aff

    character(len=4) :: etype,ftype
!
!---------------------------------------------------------------------------------------
!
    iprint=0
    dump_exa=55
      open(unit=dump_exa,file='output/face_exa', &
        form='formatted',access='sequential',status='unknown')
    dump_num=56
      open(unit=dump_num,file='output/face_num', &
        form='formatted',access='sequential',status='unknown')
    dump_err=57
      open(unit=dump_err,file='output/face_err', &
        form='formatted',access='sequential',status='unknown')
    dump_f3d=58
      open(unit=dump_f3d,file='output/face_3d', &
        form='formatted',access='sequential',status='unknown')

    ntri_total = 0
    npt_total = 0
    nf_act = 0
    do jf=1,NUMF

      mdlf = NUMC+NUMV+NUME+jf

      if (NODES(mdlf)%act.ne.1) cycle
      nf_act = nf_act + 1

      iblock = TRIANGLES(jf)%BlockNo(1)

      call decode(iblock,mdle,ires)
      idom = ALLC_DOM(mdle)
!
!     initialize global quantities
      errorU=0.d0 ; rnormU=0.d0; errorF=0.d0 ; rnormF=0.d0
!
      call find_ndof_poly(mdlf, ndofH,ndofE,ndofV,ndofQ)
      nord_mdlf= NODES(mdlf)%order
      ndofU = (nord_mdlf+1)*(nord_mdlf+2)/2      
      ndofF=ndofV
      

      if (NRFVAR.gt.0) then
      zdofF = 0.d0
      zdofF(1:NRFVAR,1:ndofF) = NODES(mdlf)%zdofF(1:NRFVAR,1:ndofF)
 !      write(*,*)'zdofF ='
 !      do k=1,NRFVAR
 !        write(*,1235) zdofF(k,1:3)
 ! 1235   format(3e12.5)
      ! enddo
      endif
      if (NRUVAR.gt.0) then
      zdofU = 0.d0
      zdofU(1:NRUVAR,1:ndofU) = NODES(mdlf)%zdofU(1:NRUVAR,1:ndofU)
      ! write(*,*)'zdofU ='
      ! do k=1,NRUVAR
        ! write(*,1235) zdofU(k,1:3)
 ! 1235   format(3e12.5)
      ! enddo
      endif

      ftype = 'tria'

      call face_vert_list(mdlf,nrv_f,nverts_f)
      allocate(xverts_f(3,nrv_f))
      do jv=1,Nrv_f
        xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)
      enddo

!
! get face affine coordinates points b0,b1,b2
      call poly_affine_2d(xverts_f,nrv_f,b0,b1,b2)
      b_aff(:,1) = b0
      b_aff(:,2) = b1
      b_aff(:,3) = b2
! 
!       set up face order in usual structure
      nordf = 0
      nordf(1:4) = NODES(Mdlf)%order
      norie = 0
      if (nrv_f.eq.3) then
          norie(1:3) = NODES(Mdlf)%face_orient(1:3)
          norie(3) = 1 - norie(3)
      endif
!       transform vertices to element coordinate system
      
!       get face normal in physical coordinates
      call face_normal(xverts_f(:,1:3),fn)

      allocate(affverts_f(2,Nrv_f))
      affverts_f = 0.d0
      call poly_x_to_affine_2d(xverts_f,nrv_f,b_aff,affverts_f)

  !   allocate memory for evaluation points on current face
      allocate(verts(2,(nrv_f-2)*12+nrv_f),trian(3,16*(nrv_f-2)))

      call break_poly2d(mdlf,NODES(Mdlf)%order,nrv_f,affverts_f,ntri,npt,trian,verts)
!
  !   get physical points from face affine coordinates
      allocate(xloc(3,npt))
      xloc = 0.d0
      call poly_affine_2d_to_x(verts,npt,b_aff,xloc)
      call trian_centr_area(b0,b1,b2,x_aux,b_area)
      rjac_f = 2*b_area
      do l=1,ntri   
        write(dump_f3d,*) mdlf,",",npt_total+trian(1,l),",", &
                                   npt_total+trian(2,l),",", &
                                   npt_total+trian(3,l),","
      enddo

      ! write(*,*)'face_error: mdlf, fn=',mdlf,fn


!    .....loop through face integration points
      do l=1,npt
!
!    .......quadrature point in face coordinates
        t(1:2) = verts(1:2,l); x(1:3) = xloc(1:3,l)
!           evaluate trial (2d face) and test (3d element) functions
        call shape2DQ(ftype,t,nordf, nrdofQ_f,ShapQ_f)
!
        shapQ_f(1:nrdofQ_f) = shapQ_f(1:nrdofQ_f)/rjac_f
!           evaluate trial (2d face) and test (3d element) functions
        call shape2DH(ftype,t,nordf, norie, nrdofH_f,ShapH_f,gradH_f)
!             
!
!     evaluate the approximate solution
        ZsolU(1:MAXEQNU)=ZERO
        do k=1,nrdofH_f
          do ivar=1,MAXEQNU
            ZsolU(ivar) = ZsolU(ivar) + zdofU(ivar,k)*shapH_f(k)
          enddo
        enddo
!
!     evaluate the approximate solution
        ZsolF(1:MAXEQNF)=ZERO
        do k=1,nrdofQ_f
          do ivar=1,MAXEQNF
            ZsolF(ivar) = ZsolF(ivar) + zdofF(ivar,k)*shapQ_f(k)
          enddo
        enddo

!           -- EXACT SOLUTION --

        

        call exact(x,fn,idom, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV,zvalU, zvalF, zvalQ,zdvalQ,zd2valQ)
!
!
        ErrorU(1:NRUVAR) = zvalU(1:NRUVAR) - zsolU(1:NRUVAR)
        ErrorF(1:NRFVAR) = zvalF(1:NRFVAR) - zsolF(1:NRFVAR)
        if (l.le.nrv_f) then
          original = 1
        else
          original = 0
        endif
        write(dump_exa,9876) nf_act,npt_total+l,original,x(1:3),zvalU( 1:3),zvalF( 1:3)
        write(dump_num,9876) nf_act,npt_total+l,original,x(1:3),zsolU( 1:3),zsolF( 1:3)
        write(dump_err,9876) nf_act,npt_total+l,original,x(1:3),errorU(1:3),errorF(1:3)
 9876   format(3(i4.2,","),8(e22.15,","),e22.15)
!         loop over integration points
      enddo

      ! write(*,*) 'mdlf,ErrorF,rnormF =',mdlf,ErrorF,rnormF

      deallocate(affverts_f,xverts_f,trian,verts,xloc)

      ntri_total = ntri_total + ntri
      npt_total  = npt_total  + npt      

    enddo

    write(*,*) "***Face visualization: Dumped data for ",ntri_total," triangles and ",npt_total," points"
    
    close(dump_num)
    close(dump_exa)
    close(dump_err)
    close(dump_f3d)

end subroutine generate_visual_poly_face_ct



! subroutine break_poly3d(mdle,Nrf,nord,ntri,npt,trian,verts)

!   use parametersDPG, only: MAXNINT3ADD,MAX_NRVF
!   use data_structure3D_poly
!   use connectivity_poly

!   implicit none

!   integer,              intent(in ) :: mdle, nord,nrf
!   integer,              intent(out) :: nint
!   real*8, dimension(3,(Nrf-3)*(MAX_NRVF-2)*MAXNINT3ADD),  intent(out) :: Xiloc
!   real*8, dimension((Nrf-3)*(MAX_NRVF-2)*MAXNINT3ADD  ),  intent(out) :: Wxi

!   integer, dimension(19)            :: norder
!   integer                           :: nrfc,nrv,nre,nint_tet,jf,jv,kv,ipt,nrv_f,mdlf,loc,iflag
!   integer, allocatable              :: nverts(:),Nfaces(:),Nedges(:),Norientf(:),nverts_f(:)
!   real*8, allocatable               :: Xvertl(:,:),xverts_f(:,:)
!   real*8                            :: r_e,maptetdet
!   real*8, dimension(3)              :: Xg_e,vn1,vn2,vn3,vn4
!   real*8, dimension(3,3)            :: qrot_e,maptet,maptetinv
!   real*8, dimension(3,MAXNINT3ADD)  :: xiq_tet
!   real*8, dimension(MAXNINT3ADD)    :: wq_tet
!   real*8, dimension(3,6)            :: newverts_master,newverts

!   norder = 0
!   norder(1:15) = nord 



!   newverts_master (1,:) = (/0.5d0,0.d0 ,0.d0 /)
!   newverts_master (2,:) = (/0.d0 ,0.5d0,0.d0 /)
!   newverts_master (3,:) = (/0.d0 ,0.d0 ,0.5d0/)
!   newverts_master (4,:) = (/0.5d0,0.5d0,0.d0 /)
!   newverts_master (5,:) = (/0.d0 ,0.5d0,0.5d0/)
!   newverts_master (6,:) = (/0.5d0,0.5d0,0.d0 /)
! ! retrieve element node data
!   call elem_nodes_poly(Mdle,Nfaces,Norientf,Nrfc,Nedges,Nre,Nverts,Nrv)
! ! allocate array for vertices
!   allocate(Xvertl(3,Nrv))

!   xiloc= 0.d0
!   wxi = 0.d0

! ! compute centroid, radius, rotation (unitary) matrix and vertices in local coord
!   call elem_vertex_coordl(Mdle,Nfaces(1:Nrf),Norientf(1:Nrf),Nrf,Nverts(1:Nrv),Nrv, &
!                               Xg_e,R_e,Qrot_e,Xvertl)
! !
! ! save the first vertex of the polyhedron as vertex 1 of all subtetrahedra
!   vn1(1:3) = Xvertl(1:3,1)
!   call set_3Dint('tetr',norder, nint_tet,xiq_tet,wq_tet)
!   ntet = 0
!   npt = 0
! !     loop over faces
!   do jf = 1,nrf
!     mdlf = Nfaces(jf)
!     call face_vert_list(mdlf,nrv_f,nverts_f)
! !       Check if element vertex 1 is on the face vertex list.
! !       If not, the face is a valid tetrahedron base.
!     call locate(Nverts(1),nverts_f,nrv_f,loc)
!     if (loc.eq.0) then
!       allocate(xverts_f(3,nrv_f))
!       do jv=1,Nrv_f
!         xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)-Xg_e(1:3)
!       enddo
! !       transform vertices to element coordinate system
!       xverts_f = matmul(transpose(Qrot_e),xverts_f)
!       xverts_f = xverts_f / R_e
! !     complete subtetrahedron vertex list vn2,3,4
!       vn2(1:3) = xverts_f(1:3,1)
! !         loop on subtetrahedra
!       do kv = 1,Nrv_f-2
! !         check if face normal locally goes outward (Norientf = 0)
!         if(Norientf(jf).eq.0) then
!           vn3(1:3) = xverts_f(1:3,kv+1)
!           vn4(1:3) = xverts_f(1:3,kv+2)
!         else
! !           if normal goes inward swap order of vertices 3 and 4            
!           vn3(1:3) = xverts_f(1:3,kv+2)
!           vn4(1:3) = xverts_f(1:3,kv+1)
!         endif

!         call tetra_affine_map(vn1,vn2,vn3,vn4,maptet)
!         newverts = 0.d0
!         if (nord.le.2) then
!           subtets = 1
!         else 
!           subtets = 8
          
!         endif


!         endif

!         call geom(maptet,maptetinv,maptetdet,iflag)
! !  
! !        ...loop through integration points
!         do ipt=1,nint_tet

!           npt = npt + 1
! !         copy quadrature point
!           xiloc(1:3,nint) = xiq_tet(1:3,ipt)
! !         transform point to the local coordinates
!           xiloc(1:3,nint) = matmul(maptet,xiloc(1:3,nint)) + vn1
! !         copy weight and multiply by tetrahedron jacobian
!           wxi(nint) = wq_tet(ipt)*maptetdet

!         enddo
!       enddo

!       deallocate(xverts_f,nverts_f)

!     endif

!   enddo
!   deallocate(Nfaces,Nedges,Norientf,Nverts)


! end subroutine break_poly3d




subroutine break_poly2d(mdlf,nord,nrv_f,rotverts_f,ntri,npt,trian,verts)
  use parametersDPG, only : MAXNINT2ADD
  use control,       only : INTEGRATION

  implicit none

  integer,                               intent(in ) :: mdlf, nord,nrv_f
  real*8, dimension(2,nrv_f)            ,intent(in ) :: rotverts_f
  integer,                               intent(out) :: ntri,npt
  integer,dimension(3,16*(nrv_f-2))      ,intent(out) :: trian
  real*8, dimension(2,(nrv_f-2)*12+nrv_f),intent(out) :: verts

  integer                                    :: kv,ipt
  real*8, dimension(2)                       :: vt1,vt2,vt3
  real*8, dimension(2,2)                     :: maptri
  real*8, dimension(2,12)                    :: newverts,newverts_master  
  integer,dimension(15)                      :: it


  newverts_master(:,1) = (/0.5d0, 0.d0/)
  newverts_master(:,2) = (/0.d0 ,0.5d0/)
  newverts_master(:,3) = (/0.5d0,0.5d0/)
  newverts_master(:,4) = (/0.25d0,0.d0/)
  newverts_master(:,5) = (/0.25d0,0.25d0/)
  newverts_master(:,6) = (/0.d0,0.25d0/)
  newverts_master(:,7) = (/0.75d0,0.d0/)
  newverts_master(:,8) = (/0.75d0,0.25d0/)
  newverts_master(:,9) = (/0.5d0,0.25d0/)
  newverts_master(:,10) = (/0.25d0,0.5d0/)
  newverts_master(:,11) = (/0.25d0,0.75d0/)
  newverts_master(:,12) = (/0.d0,0.75d0/)

! define first vertex of every subtriangle always as vertex 1 of current face
  verts = 0.d0
  verts(1:2,1:nrv_f) = Rotverts_f(1:2,1:nrv_f)
  vt1(1:2) = verts(1:2,1)
!
! loop on subtriangles
  ntri = 0
  npt  = nrv_f
  do kv = 1, nrv_f-2
    
    if (nord.le.1) then      
      trian(1:3,ntri+1)=(/1,kv+1,kv+2/)
      ntri = ntri +1
    else
      vt2(1:2) = verts(1:2,kv+1)
      vt3(1:2) = verts(1:2,kv+2)
      call trian_affine_map(vt1(1:2),vt2(1:2),vt3(1:2),maptri)
      newverts = matmul(maptri,newverts_master)
      do ipt=1,12
        newverts(:,ipt) = newverts(:,ipt) + vt1(:)
      enddo

      verts(:,npt+1:npt+12)=newverts(:,1:12)
      
      it(1) = 1
      it(2) = kv + 1
      it(3) = kv + 2
      do ipt=1,12
        it(3+ipt) = npt+ ipt
      enddo

      npt = npt + 12

      trian(:,ntri + 1) = (/it(1 ),it(7 ),it(9 ) /)
      trian(:,ntri + 2) = (/it(7 ),it(8 ),it(9 ) /)
      trian(:,ntri + 3) = (/it(7 ),it(4 ),it(8 ) /)
      trian(:,ntri + 4) = (/it(9 ),it(8 ),it(5 ) /)
      trian(:,ntri + 5) = (/it(4 ),it(12),it(8 ) /)
      trian(:,ntri + 6) = (/it(4 ),it(10),it(12) /)
      trian(:,ntri + 7) = (/it(8 ),it(12),it(13) /)
      trian(:,ntri + 8) = (/it(8 ),it(13),it(5 ) /)
      trian(:,ntri + 9) = (/it(5 ),it(13),it(15) /)
      trian(:,ntri + 10)= (/it(10),it(2 ),it(11) /)
      trian(:,ntri + 11)= (/it(10),it(11),it(12) /)
      trian(:,ntri + 12)= (/it(12),it(11),it(6 ) /)
      trian(:,ntri + 13)= (/it(12),it(6 ),it(13) /)
      trian(:,ntri + 14)= (/it(13),it(6 ),it(14) /)
      trian(:,ntri + 15)= (/it(13),it(14),it(15) /)
      trian(:,ntri + 16)= (/it(15),it(14),it(3 ) /)

      ntri = ntri +16

    endif
  enddo

end subroutine break_poly2d