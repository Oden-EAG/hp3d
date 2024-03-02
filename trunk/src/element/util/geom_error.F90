#if HP3D_DEBUG

!------------------------------------------------------------------------
!> @brief calculate geometry error. Refer to hp book, vol. 2,
!!           page 103, formula 5.4.
!!
!> @param[out] Err   - norm of the error
!> @param[out] Rnorm - norm of the exact geometry map
!!
!> @date Dec 14
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
   subroutine geom_error(Err,Rnorm)
!
      use data_structure3D , only : NRELES,ELEM_ORDER
      use environment      , only : QUIET_MODE,L2GEOM
!
      implicit none
      real(8), intent(out) ::  Err, Rnorm
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
        call geom_error_elem(mdle, derr,dnorm)
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
!     check that sons error is less than father error
      call check_geom_error
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
        write(*,*) 'geom_error: increase maxvis!'
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
end subroutine geom_error
!
!
!
!-------------------------------------------------------------------------
!> @brief geometry error square on element Mdle
!!
!> @param[in   ] Mdle  - middle node number
!> @param[in   ] Ierr  - 0 : L2 error only; 1 : full H1 error
!> @param[inout] Derr  - error
!> @param[inout] Dnorm - norm of solution
!!
!> @date Feb 2023
!-------------------------------------------------------------------------
subroutine geom_error_elem(Mdle, Derr,Dnorm)
!
      use data_structure3D , only : NODES,MAXbrickH,MAX_NINT3
      use control          , only : INTEGRATION
      use environment      , only : L2GEOM
      use node_types
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
!     shape functions
      real(8) :: shapH(MAXbrickH),dshapH(3,MAXbrickH)
!
!     geometry
      real(8), dimension(3)   :: xi,x_ex,x_hp
      real(8), dimension(3,3) :: dx_hpdxi,dxidx_ex,dx_exdxi,dx_hpdx_ex
!
!     quadrature
      real(8) :: xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
!
      integer :: ntype
7001  format(' geom_error_elem: Mdle,type = ',i10,2x,a4)
!
      real(8),dimension(3,3),parameter :: del = &
      reshape((/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/),(/3,3/))
!
      integer :: i,j,k,l,nint,nrdofH,iprint,iflag
      real(8) :: wa,weight,rjac,s
!---------------------------------------------------------------------
!
      iprint=0
!
      if (iprint.eq.1) then
         write(*,7001) Mdle,S_Type(NODES(Mdle)%ntype)
      endif
!
!     order of approximation, orientations, geometry dofs
      call find_order( Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor(     Mdle, xnod)
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
        call shape3DH(ntype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,dshapH)
!
!       ISOPARAMETRIC MAP : x_hp = x_hp(xi)
        x_hp(1:3)=0.d0 ; dx_hpdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
           x_hp(1:3)=x_hp(1:3)+xnod(1:3,k)*shapH(k)
           do i=1,3
              dx_hpdxi(1:3,i)=dx_hpdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
           enddo
        enddo
!
!       EXACT GEOMETRY MAP : x_ex = x_ex(xi)
        call exact_geom(Mdle,xi, x_ex,dx_exdxi)
        call geom(dx_exdxi, dxidx_ex,rjac,iflag)
!
!       check Jacobian
        if (iflag.ne.0) then
           write(*,*) 'geometry_error_element: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
           write(*,*) '        rjac = ',rjac
           call pause
        endif
!
!       total weight
        weight=wa*rjac

        if (iprint.eq.1) then
           write(*,7010) l,xi(1:3)
7010       format(' l,xi = ',i3,2x,3f8.3)
           write(*,7011) x_ex(1:3)-x_hp(1:3)
7011       format(' x_ex - x_hp = ',3(e12.5,2x),3x)
        endif
!
!       L2 contribution
        do i=1,3
           Dnorm=Dnorm+(x_ex(i)        )**2*weight
           Derr =Derr +(x_ex(i)-x_hp(i))**2*weight
        enddo
!
!       H1 contribution
        if (.not. L2GEOM) then
!
!          dx_hp / dx_ex = dx_hp / dxi * dxi / dx_ex
           do i=1,3
              do j=1,3
                 s=0.d0
                 do k=1,3
                    s=s+dx_hpdxi(i,k)*dxidx_ex(k,j)
                 enddo
                 dx_hpdx_ex(i,j)=s
              enddo
           enddo
!
!          accumulate
!          notice : dx_ex / dx_ex = delta
           do i=1,3
              do j=1,3
                 Dnorm=Dnorm+(del(i,j)                )**2*weight
                 Derr =Derr +(del(i,j)-dx_hpdx_ex(i,j))**2*weight
              enddo
           enddo

        endif
!
        if (iprint.eq.1) then
           write(*,9999)Mdle,Derr,Dnorm
9999       format(' Mdle,error,norm = ',i7,2x,2(e12.5,2x))
        endif
!
!     end of loop over integration points
      enddo
!
!     store error in data structure
      NODES(Mdle)%error(:,0)=0.d0 ; NODES(Mdle)%error(0,0)=Derr
!
end subroutine geom_error_elem
!
!
!------------------------------------------------------------------------
subroutine check_geom_error
    use data_structure3D
    implicit none
    integer :: mdle,i,j,nfath,ns,nrsons,nfail
    real(8) :: err_fath,err_sons

!  ...lower visitation flag
      call reset_visit
!
!  ...loop over active elements
      nfail=0
      do i=1,NRELES
        mdle = ELEM_ORDER(i)
!
!  .....if no father, skip
        nfath=NODES(mdle)%father
        if (nfath.lt.0) cycle
!
!  .....if father has already been visited, skip
        if (NODES(nfath)%visit.ne.0) cycle
!
        err_fath = NODES(nfath)%error(0,0)
        call find_nsons(nfath, nrsons)
!
!  .....loop over children
        err_sons=0.d0
        do j=1,nrsons
!          ns=NODES(nfath)%sons(j)
          ns=Son(nfath,j)
!
          select case(NODES(ns)%ntype)
          case(MDLB,MDLP,MDLN,MDLD)
            err_sons=err_sons+NODES(ns)%error(0,0)
          endselect
        enddo
!
!  .....C H E C K
!  .....skip if within machine precision
        if (abs(err_sons-err_fath).lt.1.d-14) cycle
        if (err_sons.gt.err_fath) then
          nfail=nfail+1
          write(*,9000)nfail,nfath,err_fath
9000      format(' nfail,nfath,err_fath = ',2(i7,2x),e12.5)
!  .......loop over children
          err_sons=0.d0
          do j=1,nrsons
!            ns=NODES(nfath)%sons(j)
            ns=Son(nfath,j)
!
            select case(NODES(ns)%ntype)
            case(MDLB,MDLP,MDLN,MDLD)
9001          format('   ns,type,err,err_sons = ',i7,2x,a4,2x,2(e12.5,2x))
              err_sons=err_sons+NODES(ns)%error(0,0)
              write(*,9001)ns,S_Type(NODES(ns)%ntype),NODES(ns)%error(0,0),err_sons
            endselect
          enddo
        endif
!
!  .....raise visitation flag
        NODES(nfath)%visit=1
!
!  ...end of loop over active elements
      enddo
!
!  ...lower visitation flag
      call reset_visit
!
end subroutine check_geom_error
!
!------------------------------------------------------------------------
subroutine display_geom_error(Nfath)
      use data_structure3D
      implicit none
      integer,intent(in) :: Nfath
      real(8) :: err_sons
      integer :: j,nrsons,ns
!
      write(*,5000)Nfath,NODES(Nfath)%error(0,0)
5000  format(' nfath,err_fath = ',i7,2x,e12.5)
!
      call find_nsons(nfath, nrsons)
!
!  ...loop over children
      err_sons=0.d0
      do j=1,nrsons
!        ns=NODES(nfath)%sons(j)
        ns=Son(nfath,j)
!
        select case(NODES(ns)%ntype)
        case(MDLB,MDLP,MDLN,MDLD)
          err_sons=err_sons+NODES(ns)%error(0,0)
          write(*,5001)ns,NODES(ns)%error(0,0),err_sons
5001      format('   ns,err,err_sons = ',i7,2x,2(e12.5))
        endselect
      enddo
      if (err_sons.gt.NODES(Nfath)%error(0,0)) then
        write(*,*)'  F A I L'
      endif
!
!
   end subroutine display_geom_error

#endif
