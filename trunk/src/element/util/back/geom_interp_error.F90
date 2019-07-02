!------------------------------------------------------------------------
!> Purpose : routine evaluates the geometry H_0^1 interpolation error. 
!!           Refer to hp book, vol. 2, page 103, formula 5.4.
!!
!! @param[out] Err  
!!
!! @revision Nov 12
!------------------------------------------------------------------------
subroutine geom_interp_error(Err)
!
      use data_structure3D , only : NRELES, NODES 
      use parameters       , only : MAXbrickH
!
      implicit none
      real*8, intent(out) :: Err
!
      real*8  :: derr, dnorm, err_rate
      integer :: iprint, mdle, iel, ierr, i,nb,nrgdof,nvoid
      real*8,dimension(3,8)         :: xsub
      real*8,dimension(3,MAXbrickH) :: xnod
      integer,dimension(19) :: norder
      integer,dimension(12) :: nedge_orient
      integer,dimension( 6) :: nface_orient
!
!  ...variables saved between subsequent calls
      integer, save :: ivis        = 0
      integer, save :: nrgdof_save
      real*8 , save ::    err_save
      real*8 , dimension(10,4), save :: rwork
      integer, dimension(10,1), save :: iwork
!-----------------------------------------------------------------------
!
      iprint=0
!
!  ...initialize global quantities  
      Err=0.d0
!      
!  ...loop over active elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call geom_interp_error_elem(mdle, derr)
        Err=Err+derr
!
!  .....printing
        if (iprint.eq.1) then
          write(*,7004) mdle,derr
7004      format(' geom_interp_error: mdle,err^2 = ',i7,2x,e12.5,2x)
        endif
      enddo
      Err=sqrt(Err)
!  
!  ...number of geometry dof, namely number of H1 dof for a single component
      call find_nrdof(nrgdof,nvoid,nvoid,nvoid)
!  
      err_rate=0.d0 
!
!  ...if not 1st visit, compute rate
      if (ivis.gt.0) then
        if (nrgdof.gt.nrgdof_save) then
          if (Err.gt.0.d0) then
            err_rate = log(err_save/Err)/log(float(nrgdof_save)/nrgdof)
      endif ; endif ; endif      
!
!  ...save error and number of geometry dofs
      err_save=Err ; nrgdof_save=nrgdof
!      
!  ...raise visitation flag
      ivis=ivis+1
!     
!  ...store
      rwork(ivis,1)=Err
      rwork(ivis,2)=err_rate
      iwork(ivis,1)=nrgdof
!
!  ...printing      
      write(*,*)'-- Geometry Interp. Error Report --'
      do i=1,ivis     
        write(*,9999)i,iwork(i,1),rwork(i,1:2)
9999    format(' i,gdofs,err,rate = ',i2,2x,i6,2x,e12.5,2x,f9.6)
      enddo
      write(*,*)''
!
!  ...additional printing
      if (iprint.eq.1) then
        write(*,7002)nrgdof_save,nrgdof
7002    format(' geometry_error: nrgdof old,new = ',2(i7,2x))
      endif
!
!
end subroutine geom_interp_error
!
!
!----------------------------------------------------------------------
!> Purpose : element H01-interpolation error^2    
!!
!! @param[in ] Mdle - an element middle node number 
!! @param[out] Derr - the square of the error
!!
!! @revision Nov 12     
!----------------------------------------------------------------------
subroutine  geom_interp_error_elem(Mdle, Derr)
!
      use data_structure3D , only : NODES
      use parameters       , only : MAXbrickH, MAX_NINT3
!
      implicit none
      integer,intent(in ) :: Mdle
      real*8, intent(out) :: Derr
!      
      character(len=4) :: etype
      integer,dimension(12) :: nedge_orient
      integer,dimension( 6) :: nface_orient
      integer,dimension(19) :: norder
      real*8,dimension(3,MAXbrickH) :: xnod
      real*8,dimension(3) :: x_ex,x_hp,xi
      real*8,dimension(  MAXbrickH) :: shapH
      real*8,dimension(3,MAXbrickH) :: dshapH
      real*8,dimension(3,3) :: dx_hpdxi,dx_exdxi,dxidx_ex,dx_hpdx_ex
      real*8,dimension(3,MAX_NINT3) :: xiloc
      real*8,dimension(  MAX_NINT3) :: wxi
      real*8,dimension(3,3),parameter :: del = &
      reshape((/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/),(/3,3/))
      real*8  :: s,rjac,weight,wa
      integer :: nrdof,nint,i,j,k,l,iflag
      integer :: iprint
!----------------------------------------------------------------------
!
      iprint=0
!     
      etype=NODES(Mdle)%type
!
!  ...element order of approximation, orientations, gdof's
      call find_order( Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor(     Mdle, xnod)
!     
!  ...Gauss integration points and weights
      call set_3Dint(etype,norder, nint,xiloc,wxi)
!
!  ...initialize
      Derr=0.d0
!      
!  ...loop through integration points
      do l=1,nint
!
!  .....integration point and weight
        xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!        
!  .....shape functions
        call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdof,shapH,dshapH)
!
!  .....ISOPARAMETRIC MAP : x_hp = x_hp(xi)
        x_hp(1:3)=0.d0 ; dx_hpdxi(1:3,1:3)=0.d0 
        do k=1,nrdof
          x_hp(1:3) = x_hp(1:3) + xnod(1:3,k)*shapH(k)
          do i=1,3
            dx_hpdxi(1:3,i) = dx_hpdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
          enddo
        enddo
!
!  .....EXACT GEOMETRY MAP : x_ex = x_ex(xi)
        call exact_geom(mdle,xi, x_ex,dx_exdxi) 
        call geom(dx_exdxi, dxidx_ex,rjac,iflag)
!
!  .....dx_hp / dx_ex = dx_hp / dxi * dxi / dx_ex
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
!  .....accumulate
!       notice : dx_ex / dx_ex = delta
        weight=wa*rjac
        do i=1,3
          do j=1,3
            Derr = Derr + (del(i,j)-dx_hpdx_ex(i,j))**2*weight
          enddo
        enddo
!     
!  ...end of loop through integration points
      enddo
!      
!  ...printing
      if (iprint.eq.1) then
        write(*,7000)Mdle,etype
 7000   format(' Mdle,etype,Derr = ',i8,2x,a4,2x,e12.5)       
      endif
!
!
endsubroutine geom_interp_error_elem
