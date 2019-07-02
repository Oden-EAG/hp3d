subroutine check_jacobian

  !-------------------------------------------------------------------
  use data_structure3D
  use control , only : EXGEOM
  implicit none
  !-------------------------------------------------------------------
  ! order of approximation
  integer, dimension(19) :: norder
  integer, dimension(5)  :: nordf
  ! reference and physical coordinates
  real*8, dimension(3)   :: xi,x
  real*8, dimension(2)   :: t
  real*8, dimension(3,3) :: dxdxi,dxidx
  real*8, dimension(3,2) :: dxdt ,dxidt
  ! Gauss points and weights 
  real*8, dimension(3,MAX_NINT3) :: xiloc
  real*8, dimension(  MAX_NINT3) :: wxi
  real*8, dimension(2,MAXquadH ) :: tloc
  real*8, dimension(  MAXquadH ) :: wt
  ! miscellanea
  real*8 :: rjac
  integer :: mdle, i, nint, l, iflag, ndom, nrdofH, k, j, ifig, &
       iprint,nv,icheck
  character(len=4) :: ftype

  ! shape function
  real*8, dimension(  MAXbrickH) :: vshapH 
  real*8, dimension(3,MAXbrickH) :: dvshapH
  integer, dimension(12)         :: nedge_orient
  integer, dimension(6)          :: nface_orient
  real*8, dimension(3,MAXbrickH) :: xnod
  !-------------------------------------------------------------------

  iprint=0
  write(*,*)'checking elements jacobians...'

  !  ...loop over active elements
  mdle=0
  do i=1,NRELES
     call nelcon(mdle, mdle) 

     if (iprint.eq.1) then
        write(*,9999)mdle,NODES(mdle)%type
9999    format(' mdle,type = ',i7,2x,a4)        
     endif

     !  ...order, orientations, gdof's        
     call find_order( mdle, norder)
     call find_orient(mdle, nedge_orient,nface_orient)
     if (EXGEOM.eq.0) call nodcor(mdle, xnod)

     !-------------------------------------------------------------------------        
     !   E L E M E N T   I N T E R I O R                                      | 
     !-------------------------------------------------------------------------        

     !  ...integration points
     call set_3Dint(NODES(mdle)%type,norder, nint,xiloc,wxi)

     !  ...number of element vertices, neeed to perfrom incremental check
     nv=nvert(NODES(mdle)%type)

     !  ...loop over integration points
     do l=1,nint

        !  ...integration point
        xi(1:3) = xiloc(1:3,l)

        !  ...flag             
        icheck=0
1000    continue                        

        !  ...geometry map 
        select case(EXGEOM)
           !  ...parametric element
        case(0)
           !  ...shape functions
           call shape3H(NODES(Mdle)%type,xi,norder,nedge_orient, &
                nface_orient, nrdofH,vshapH,dvshapH)
           !  ...accumulate
           x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0 
           do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
              do j=1,3
                 dxdxi(1:3,j) = dxdxi(1:3,j) + xnod(1:3,k)*dvshapH(j,k)
              enddo

              !  ...incremental check of godf's and jacobian
              if (icheck.ne.0) then
                 rjac=0.d0 ; if (k.ge.nv) call geom(dxdxi, dxidx,rjac,iflag)
                 write(*,6999)k,xnod(1:3,k),rjac
6999             format(' k,xnod(1:3,k),rjac = ',i3,2x,3(e12.5,2x),2x,e12.5)
              endif

           enddo
           !  ...exact geometry element  
        case(1)
           call exact_geom(mdle,xi, x,dxdxi)
        endselect

        !  ...return if check was performed
        if (icheck.ne.0) return

        !  ...check jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
        if (iflag.ne.0) then
           icheck=1
           call find_domain(mdle, ndom) 
           write(*,*)'Interior point FAIL'
           write(*,7000) mdle,NODES(mdle)%type,ndom,xi,rjac
7000       format(' mdle,type,ndom,xi,rjac = ',i8,2x,a4,2x,i2,2x, &
                3(e12.5,1x),2x,e12.5)            
           goto 1000
        endif

        !  ...loop over integration points
     enddo

     !-------------------------------------------------------------------------        
     !   E L E M E N T   F A C E S                                            | 
     !-------------------------------------------------------------------------        

     !  ...loop over element faces
     do ifig=1,nface(NODES(mdle)%type)

        !  ...face integration points
        ftype=face_type(NODES(mdle)%type,ifig)
        call face_order(NODES(mdle)%type,ifig,norder, nordf)
        call set_2Dint(ftype,nordf, nint,tloc,wt)

        !  ...loop through face integration points
        do l=1,nint

           !  ...Gauss point
           t(1:2)=tloc(1:2,l)

           !  ...determine the master element coordinates
           call face_param(NODES(mdle)%type,ifig,t, xi,dxidt)

           select case(EXGEOM)
              !  ...parametric element            
           case(0)
              !  ...derivatives and values of the shape functions
              call shape3H(NODES(mdle)%type,xi,norder,nedge_orient, &
                   nface_orient, nrdofH,vshapH,dvshapH)
              !  ...accumulate
              x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0 
              do k=1,nrdofH
                 x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
                 do j=1,3
                    dxdxi(1:3,j) = dxdxi(1:3,j) + xnod(1:3,k)*dvshapH(j,k)
                 enddo
              enddo
              !  ...exact geometry element              
           case(1)
              call exact_geom(mdle,xi, x,dxdxi)
           endselect

           !  ...check jacobian
           call geom(dxdxi, dxidx,rjac,iflag)
           if (iflag.ne.0) then
              call find_domain(mdle, ndom) 
              write(*,7001)ifig
7001          format(' Face point FAIL, ifig = ',i1)
              write(*,7000) mdle,NODES(mdle)%type,ndom,xi,rjac
           endif

           !  ...end of loop over integration points
        enddo

        !  ...end of loop over faces
     enddo

     !  ...loop over active elements
  enddo

  write(*,*)'elements jacobians checked'

endsubroutine check_jacobian
