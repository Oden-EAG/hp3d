#if DEBUG_MODE

!------------------------------------------------------------------------
!> Purpose : calculate H1- or L2-norm of the PB interpolation error
!            for GMP parametrizations in reference coordinates
!!
!! @param[out] Err   - norm of the error
!! @param[out] Rnorm - norm of the exact geometry map
!!
!! @date Dec 14
!------------------------------------------------------------------------
!  REMARK
!
!  Recall that:
!
!     || u - u_h ||_1 =< C h^min{p,r} || u ||_r+1
!
!  The #dof's is estimated as:
!
!     #dof's = p^3 #elem's
!
!  Furthermore:
!
!     Fully isotropic ref.      -->  #elem's = meas{Omega}/h^3
!     Isotropic ref. in 2 dir.  -->  #elem's = ...        /h^2
!     Isotropic ref. in 1 dir.  -->  #elem's = ...        /h^1
!
!  Thus:
!
!     #dof's = p^3 meas{Omega}/h^{3,2,1}
!
!  In terms of h:
!
!     h = C #dof's^{-1/3,-1/2,-1}
!
!  The final estimate in terms of #dof's is:
!
!     || u - u_h ||_1 =< C #dof's^min{p,r}/{3,2,1} || u ||_r+1
!
!  Therefore, for a smooth solution (r = \infty) and fully isotropic
!  refinements:
!
!     p = 2  -->  rate = -2/3
!     p = 3  -->  rate = -3/3
!     p = 4  -->  rate = -4/3
!     ...
!
!  For a smooth solution and isotropic refinements in 2 directions :
!
!     p = 2  -->  rate = -2/2
!     p = 3  -->  rate = -3/2
!     p = 4  -->  rate = -4/3
!     ...
!------------------------------------------------------------------------
!
subroutine geometry_error(Err,Rnorm)
!
      use data_structure3D , only : NRELES,ELEM_ORDER
      use environment      , only : QUIET_MODE,L2GEOM
!
      implicit none
      real(8), intent(out) :: Err, Rnorm
!
      real(8) :: derr, dnorm, err_rate
      integer :: iprint, mdle, iel, i, nrgdof, nvoid , ic
!
      integer,parameter :: nin=13
      integer,parameter :: maxvis=2000
!
      integer, save :: ivis = 0
      integer, save :: nrgdof_save
      real(8), save ::    err_save
      real(8), dimension(maxvis,4), save :: rwork
      integer, dimension(maxvis,1), save :: iwork
!-------------------------------------------------------------------------
!
      iprint=0
!
!     initialize global quantities
      Err=0.d0 ; Rnorm=0.d0
!
!     loop over active elements
      do iel=1,NRELES
        mdle = ELEM_ORDER(iel)
        call geometry_error_elem(mdle, derr,dnorm)
!
!       accumulate
        Err=Err+derr ; Rnorm=Rnorm+dnorm
!
!       printing
        if (iprint.eq.1) then
          write(*,7004)mdle,derr,dnorm
 7004     format(' geometry_error: mdle,err^2,rnorm^2 = ',i7,2x,2(e12.5,2x))
        endif
!
!     end of loop over active elements
      enddo
!
!     the much neglected square root!
      Err=sqrt(Err) ; Rnorm=sqrt(Rnorm)
!
!     number of geometry dof, namely number of H1 dof for a single component
      call find_nrdof(nrgdof,nvoid,nvoid,nvoid)
!
      err_rate=0.d0
!
!     if not 1st visit, compute rate
      if (ivis.gt.0) then
        if (nrgdof.gt.nrgdof_save) then
          if (Err.gt.0.d0) then
            err_rate = log(err_save/Err)/log(float(nrgdof_save)/nrgdof)
      endif ; endif ; endif
!
!     save error and number of gdofs
      err_save=Err ; nrgdof_save=nrgdof
!
!     raise visitation flag
      ivis=ivis+1
!
if (.not. QUIET_MODE) then
!
!     check
      if (ivis > maxvis) then
        write(*,*) 'geometry_error: increase maxvis!'
        stop
      endif
!
!     store
      rwork(ivis,1)=Err
      rwork(ivis,2)=Rnorm
      rwork(ivis,3)=Err/Rnorm
      rwork(ivis,4)=err_rate
      iwork(ivis,1)=nrgdof
!
endif
!
!     printing
!
!     -- 1st visit --
      if (ivis == 1) then
!
!       open file for printing
        open(unit   = nin                    , &
             file   = './files/dump_geo_err' , &
             form   = 'formatted'            , &
             access = 'sequential'           , &
             status = 'unknown'              , &
             iostat = ic)
        if (ic /= 0) then
          write(*,*)'error_geom: COULD NOT OPEN FILE! [0]'
          stop
        endif
!
!       print header
if (.not. L2GEOM) then
        write(nin,*)'-- Geometry Error Report --'
else
        write(nin,*)'-- Geometry Error Report (L2 only) --'
endif
        write(nin,1000)
 1000   format('          Gdofs // ' , &
                 '        Error // ' , &
                 '         Norm // ' , &
                 '   Rel. Error //'  , &
                   '       Rate '        )
!
!     -- subsequent visits --
      else
!
!       append to file
        open(unit     = nin                    , &
             file     = './files/dump_geo_err' , &
             form     = 'formatted'            , &
             access   = 'sequential'           , &
             status   = 'old'                  , &
             position = 'append'               , &
             iostat   = ic)
        if (ic /= 0) then
          write(*,*)'error_geom: COULD NOT OPEN FILE! [1]'
          stop
        endif
      endif
!
!     print to file
      write(nin,9998) ivis,nrgdof,Err,Rnorm,Err/Rnorm,err_rate
 9998 format(1x,i3,' ; ',2x,i6,' ; ',2x,3(e12.5,' ; ',2x),f9.6)
!
!     print to screen
if (.not. QUIET_MODE) then ; write(*,*)''
  if (.not. L2GEOM) then   ; write(*,*)'-- Geometry Error Report --'
  else                     ; write(*,*)'-- Geometry Error Report (L2 only) --'
  endif
                             write(*,1000)
        do i=1,ivis        ; write(*,9998)i,iwork(i,1),rwork(i,1:4)
        enddo
                             write(*,*)''
endif
!
!     close file
      close(unit=nin,iostat=ic)
      if (ic /= 0) then
        write(*,*)'error_geom: COULD NOT CLOSE FILE!'
        stop
      endif
!
!
end subroutine geometry_error
!
!
!
!-------------------------------------------------------------------------
!> Purpose : evaluate element PB interpolation error (squared) for
!            the exact geometry map (in reference coordinates)
!!
!> @param[in   ] Mdle  - middle node number
!> @param[inout] Derr  - error
!> @param[inout] Dnorm - norm of the GMP map
!!
!> @date Dec 14
!-------------------------------------------------------------------------
subroutine geometry_error_elem(Mdle, Derr,Dnorm)
!
      use element_data
      use data_structure3D
      use control          , only : INTEGRATION
      use environment      , only : L2GEOM
!
      implicit none
      integer, intent(in   ) :: Mdle
      real(8), intent(inout) :: Dnorm, Derr
!
!     order of approx., gdof's, orientations
      integer, dimension(19)          :: norder
      real(8), dimension(3,MAXbrickH) :: xnod
      integer :: nedge_orient(12), nface_orient(6)
!
!     reference coordinates of the element vertices
      integer :: no,iflag
      real(8) :: etav(3,8)
!
!     shape functions
      real(8) :: shapH(MAXbrickH),gradH(3,MAXbrickH)
!
!     reference geometry
      real(8) :: eta(3),detadxi(3,3),dxideta(3,3),rjac
      integer :: error_flag
!
!     exact and approximate geometry
      real(8), dimension(3)   :: xi,xex,xhp
      real(8), dimension(3,3) :: dxhpdxi,dxhpdeta,dxexdeta
!
!     quadrature
      real(8) :: xiloc(3,MAX_NINT3),wxi(MAX_NINT3),wa,weight
!
      integer :: ntype
7001  format(' geometry_error_elem: Mdle,type = ',i10,2x,a4)
!
      integer :: nrv,i,j,k,l,nint,nrdofH,iprint
!
!---------------------------------------------------------------------
!
      iprint=0
!
      if (iprint.eq.1) then
         write(*,7001) Mdle,S_Type(NODES(Mdle)%ntype)
      endif
!
!     order of approximation, orientations, geometry dofs
      call find_order(Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor(Mdle, xnod)
!
!     determine the corresponding GMP block and reference coordinates
!     of the element vertices
      call refel(Mdle, iflag,no,etav)
!
!     initialize
      Derr=0.d0 ; Dnorm=0.d0

!     set up the element quadrature (use overintegration for the exact geometry map)
      ntype=NODES(Mdle)%ntype
      INTEGRATION=1
      call set_3Dint(ntype,norder, nint,xiloc,wxi)
      INTEGRATION=0
!
!     loop over integration points
      do l=1,nint
        xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!       evaluate appropriate shape functions at the point
        call shape3DH(ntype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
!
!       ISOPARAMETRIC MAP : x_hp = x_hp(xi)
        xhp(1:3)=0.d0 ; dxhpdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
           xhp(1:3)=xhp(1:3)+xnod(1:3,k)*shapH(k)
           do i=1,3
              dxhpdxi(1:3,i)=dxhpdxi(1:3,i) + xnod(1:3,k)*gradH(i,k)
           enddo
        enddo
!
!       evaluate the reference geometry
        nrv = nvert(ntype)
        call refgeom3D(Mdle,xi,etav,shapH,gradH,nrv, &
                       eta,detadxi,dxideta,rjac,error_flag)
!
!       use chain formula to compute derivatives wrt reference coordinates
        do i=1,3
          dxhpdeta(1:3,i) = dxhpdxi(1:3,1)*dxideta(1,i) &
                          + dxhpdxi(1:3,2)*dxideta(2,i) &
                          + dxhpdxi(1:3,3)*dxideta(3,i)
        enddo
!
!       compute the exact geometry map and its derivatives wrt reference coordinates
        select case(iflag)
        case(5) ; call prism(no,eta, xex,dxexdeta)
        case(6) ; call  hexa(no,eta, xex,dxexdeta)
        case(7) ; call tetra(no,eta, xex,dxexdeta)
        case(8) ; call pyram(no,eta, xex,dxexdeta)
        case default
          write(*,*) 'geometry_error_elem: Mdle,type,iflag = ',Mdle,ntype,iflag
          call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
        endselect
!
!       total weight
        weight=wa*rjac

        if (iprint.eq.1) then
           write(*,7010) l,xi(1:3)
7010       format(' l,xi = ',i3,2x,3f8.3)
           write(*,7011) xex(1:3)-xhp(1:3)
7011       format(' xex - xhp = ',3(e12.5,2x),3x)
        endif
!
!       L2 contribution
        do i=1,3
           Dnorm=Dnorm+(xex(i)       )**2*weight
           Derr =Derr +(xex(i)-xhp(i))**2*weight
        enddo
!
!       H1 seminorm contribution
        if (.not. L2GEOM) then
          do i=1,3
            do j=1,3
              Dnorm=Dnorm+(dxexdeta(i,j)              )**2*weight
              Derr =Derr +(dxexdeta(i,j)-dxhpdeta(i,j))**2*weight
            enddo
          enddo
        endif
!
        if (iprint.eq.1) then
           write(*,9999) Mdle,Derr,Dnorm
9999       format('geometry_error_elem: Mdle,error,norm = ',i7,2x,2(e12.5,2x))
        endif
!
!     end of loop over integration points
      enddo
!
!     store error in data structure
      NODES(Mdle)%error(:,0)=0.d0 ; NODES(Mdle)%error(0,0)=Derr
!
!
end subroutine geometry_error_elem

#endif
