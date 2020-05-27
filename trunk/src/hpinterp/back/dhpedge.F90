!> Purpose : update the dirichle BC on the edge
!!           H1 projection for ZnodH, L2 projection for ZnodE
!! @param[in]  Iflag        - a flag specifying which of the objects
!!                            5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Type         - element (middle node) type
!! @param[in]  Icase        - the mid-edge node case
!! @param[in]  Nedge_orient - edge orientation (not needed really)
!! @param[in]  Nface_orient - face orientation (not needed really)
!! @param[in]  Norder       - element order
!! @param[in]  Iedg         - edge number
!! @param[in]  ZdofH        - H1 dof for the element (vertex values really)
!!
!! @param[out] ZnodH        - H1 dof for the edge
!! @param[out] ZnodE        - H(curl) dof for the edge

#include "implicit_none.h"
subroutine dhpedge( &
     Iflag,No,Etav, Type,Icase, &
     Nedge_orient,Nface_orient,Norder,Iedg,  &
     ZdofH, &
     ZnodH, ZnodE)
  use parameters
  use physics
  use element_data
  implicit none
  ! ** Arguments
  !--------------------------------------------------------------------
  integer,                                    intent(in)  :: Iflag,No
  integer,                                    intent(in)  :: Icase,Iedg
  real*8,  dimension(3,8),                    intent(in)  :: Etav
  character(len=4),                           intent(in)  :: Type
  integer, dimension(12),                     intent(in)  :: Nedge_orient
  integer, dimension(6),                      intent(in)  :: Nface_orient
  integer, dimension(19),                     intent(in)  :: Norder
  VTYPE,   dimension(MAXEQNH,MAXbrickH),      intent(in)  :: ZdofH
  VTYPE,   dimension(NRCOMS*NREQNH(Icase),*), intent(out) :: ZnodH
  VTYPE,   dimension(NRCOMS*NREQNE(Icase),*), intent(out) :: ZnodE
  !
  ! ** Locals
  !--------------------------------------------------------------------
  ! for vertices
  integer, dimension(19)             :: norder_1
  real*8,  dimension(MAXbrickH)      ::  shapH
  real*8,  dimension(3,MAXbrickH)    :: dshapH
  !
  ! edge parameterization and shape function
  real*8,  dimension(3)              :: xi, eta, x, dxdt, dxidt
  real*8,  dimension(3,3)            :: dxdxi,dxdeta,detadxi
  real*8,  dimension(MAXP-1)         :: shapH_edge,dshapH_edge
  real*8,  dimension(MAXP)           :: shapE_edge

  ! geometry for solelm
  real*8 :: xnod(3,MAXbrickH)

  VTYPE :: &
       zdofH_c(MAXEQNH,MAXbrickH),zdofE_c(MAXEQNE,MAXbrickE),&
       zdofV_c(MAXEQNV,MAXbrickV),zdofQ_c(MAXEQNQ,MAXbrickQ)

  VTYPE :: &
       ZsolH(MAXEQNH),ZgradH(MAXEQNH,3), &
       ZsolE(3,MAXEQNE),ZcurlE(3,MAXEQNE), &
       ZsolV(3,MAXEQNV),ZdivV(MAXEQNV), &
       ZsolQ(MAXEQNQ),zEt

  integer :: mdle, nrp, nflag
  real*8 :: dt

  !
  ! Dirichlet BC data at a point and its Piola transformations
  VTYPE :: &
       zvalH(MAXEQNH), &
       zdvalH(MAXEQNH,3), zdvalHdxi(MAXEQNH,3),zdvalHdt(MAXEQNH), &
       zvalE(3,MAXEQNE), zvalExi(3,MAXEQNE), zvalEt(MAXEQNE), &
       zdvalE(3,MAXEQNE,3), zdvalExi(3,MAXEQNE,3), &
       zvalV(3,MAXEQNV), zdvalV(3,MAXEQNV,3)

  integer, dimension(NR_PHYSA)       :: ncase

  ! linear systems for H1, Hcurl
  real*8,  dimension(MAXP-1,MAXP-1)  :: aaH
  real*8,  dimension(MAXP, MAXP)     :: aaE

  ! load vector and solution
  VTYPE,   dimension(MAXP-1,MAXEQNH) :: zbH,zuH
  VTYPE,   dimension(MAXP  ,MAXEQNE) :: zbE,zuE

  integer :: &
       i,j,k,l,ii,ivar,ivarH,ivarE,nvarH,nvarE,naH,naE,iprint, &
       nint,ndofH_edge,ndofE_edge,nrdofH, nord, nord_old
  real*8,  dimension(MAXP+1)         :: xi_list, wa_list
  real*8                             :: t, wa

  ! to reuse precalculated one
  save nord_old, aaH, aaE
  data nord_old /0/
  !----------------------------------------------------------------------
  iprint=0 !iprint_parameter

  ! decide whether it reuse saved one
  nord = Norder(Iedg)
  if (nord.ne.nord_old) then
     aaH = 0.d0; aaE = 0.d0
  endif
  zbH = ZERO; zbE = ZERO

  call initiate_order(Type, norder_1)
  call set_1Dint(nord, nint, xi_list, wa_list)

  ! loop over integration points
  do l=1,nint
     t  = xi_list(l)
     wa = wa_list(l)

     call edge_param(Type,Iedg,t, xi,dxidt)
     if (iprint.eq.1) then
        write(*,8000) l,t
8000    format('dphmedge: INTEGRATION POINT            = ',i3,e12.5)
        write(*,8001) xi(1:3),dxidt(1:3)
8001    format('dhpmedge: EDGE PARAMETERATION xi,dxidt = ',3f8.3,2x,3f8.3)
     endif

     ! compute element shape functions
     call shape3H( &
          Type,xi, &
          norder_1,Nedge_orient,Nface_orient, &
          nrdofH,shapH,dshapH)

     ! evaluate reference coordinates of the point
     eta(1:3)         = 0.d0
     detadxi(1:3,1:3) = 0.d0
     do k=1,nrdofH
        eta(1:3) = eta(1:3) + etav(1:3,k)*shapH(k)
        do i=1,3
           detadxi(1:3,i) = detadxi(1:3,i) + etav(1:3,k)*dshapH(i,k)
        enddo
     enddo

     ! remove the contribution from vertices
     zdvalHdxi(1:MAXEQNH,1:3) = ZERO
     zvalExi(1:3,1:MAXEQNE) = ZERO
     do k=1,nrdofH
        do i=1,3
           zdvalHdxi(1:MAXEQNH,i) = &
                zdvalHdxi(1:MAXEQNH,i) - ZdofH(1:MAXEQNH,k)*dshapH(i,k)
        enddo
     enddo

     ! evaluate derivatives of the physical coordinates wrt the reference
     ! coordinates
     select case(Iflag)
     case(5);        call prism(No,eta, x,dxdeta)
     case(6);        call  hexa(No,eta, x,dxdeta)
     case(7);        call tetra(No,eta, x,dxdeta)
     case(8);        call pyram(No,eta, x,dxdeta)
     case default
        write(*,*) 'dhpedge: Type = ', Type
        call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
     end select

     ! evaluate the Dirichlet data
     call dirichlet( &
          x,Icase, &
          zvalH,zdvalH, &
          zvalE,zdvalE, &
          zvalV,zdvalV)

     if (iprint.eq.1) then
        write(*,7000) x(1:3)
7000    format('dhpedge: x = ',3f8.3)
        do ivar=1,MAXEQNH
           write(*,7001) ivar,zvalH(ivar),zdvalH(ivar,1:3)
7001       format(' ivar,zvalH(ivar),zdvalH(ivar,1:3) = ', &
                i2,2x,2e12.5,3x,3(2e12.5,2x))
        enddo
        do ivar=1,MAXEQNE
           write(*,7002) ivar,zvalE(1:3,ivar)
#if C_MODE
7002       format(' ivar,zvalE(1:3,ivar) = ', &
                i2,2x,3(2e12.5,2x))
#else
7002       format(' ivar,zvalE(1:3,ivar) = ', &
                i2,2x,3e12.5)
#endif
        enddo
     endif

     !---------------------------------------------
     ! ** parameterization
     ! evaluate dxdxi (using the exact geometry)
     dxdxi(1:3,1:3) = 0.d0
     do i=1,3
        do j=1,3
           dxdxi(1:3,i) = dxdxi(1:3,i) + dxdeta(1:3,j)*detadxi(j,i)
        enddo
     enddo
     if (iprint.eq.1) then
       write(*,7010) dxdxi(1:3,1:3)
7010   format('dhpedge: dxdxi = ',3(3f8.3,2x))
     endif

     ! add derivatives of Dirichlet data wrt master coordinates
     do i=1,3
        do j=1,3
           zdvalHdxi(1:MAXEQNH,i) = &
                zdvalHdxi(1:MAXEQNH,i) + zdvalH(1:MAXEQNH,j)*dxdxi(j,i)
           zvalExi(i,1:MAXEQNE) = &
                zvalExi(i,1:MAXEQNE) + zvalE(j,1:MAXEQNE)*dxdxi(j,i)
        enddo
     enddo
     if (iprint.eq.1) then
       do ivarE=1,MAXEQNE
         write(*,7020) ivarE,zvalExi(1:3,ivarE)
7020     format('dhpedge: ivarE = ',i2,' zvalExi = ',3(2e12.5,2x))
       enddo
     endif

     ! evalute derivative wrt edge coordinate
     zdvalHdt(1:MAXEQNH) = ZERO
     zvalEt(1:MAXEQNE) = ZERO
     do i=1,3
        zdvalHdt(1:MAXEQNH) = zdvalHdt(1:MAXEQNH) &
             + zdvalHdxi(1:MAXEQNH,i)*dxidt(i)
        zvalEt(1:MAXEQNE) = zvalEt(1:MAXEQNE) &
             + zvalExi(i,1:MAXEQNE)*dxidt(i)
     enddo
     if (iprint.ne.0) then
       write(*,7035) l,zvalEt(1)
7035   format('dhpedge: l,zvalEt = ',i3,2e12.5)
     endif
     !---------------------------------------------
     ! ** 1D shape function
     call shapHbe(t,nord,Nedge_orient(Iedg), &
          ndofH_edge,shapH_edge,dshapH_edge)
     call shapEbe(t,nord,Nedge_orient(Iedg), &
          ndofE_edge,shapE_edge)

     !---------------------------------------------
     ! ** construct matrix for projections
     do j=1,ndofH_edge
        zbH(j,1:MAXEQNH) = zbH(j,1:MAXEQNH) &
             + zdvalHdt(1:MAXEQNH)*dshapH_edge(j)*wa
        ! accumulate for the stiffness matrix
        if (nord.ne.nord_old) then
           do i=1,ndofH_edge
              aaH(i,j) = aaH(i,j) + dshapH_edge(i)*dshapH_edge(j)*wa
           enddo
        endif
     enddo

     do j=1,ndofE_Edge
        zbE(j,1:MAXEQNE) = zbE(j,1:MAXEQNE) &
            + zvalEt(1:MAXEQNE)*shapE_edge(j)*wa
        ! accumulate for the stiffness matrix
        if (nord.ne.nord_old) then
           do i=1,ndofE_Edge
              aaE(i,j) = aaE(i,j) + shapE_edge(i)*shapE_edge(j)*wa
           enddo
        endif
     enddo
     !---------------------------------------------
  enddo

  iprint = 1

  if (iprint.eq.1) then
     write(*,*) 'dhpedge: LOAD VECTOR AND STIFFNESS MATRIX (H1) = '
     do ii=1,ndofH_edge
        write(*,7036) zbH(ii,1:MAXEQNH)
        write(*,7005) aaH(ii,1:ndofH_edge)
     enddo

     write(*,*) 'dhpedge: LOAD VECTOR AND STIFFNESS MATRIX (Hcurl) = '
     do ii=1,ndofE_edge
        write(*,7036) zbE(ii,1:MAXEQNE)
        write(*,7005) aaE(ii,1:ndofE_edge)
     enddo

7036 format(4(2e12.5,2x))
7005 format(1x,10f8.3)

  endif

  !-----------------------------------
  ! solve the linear system
  naH = MAXP-1
  naE = MAXP

  ! if it is not same order, factorization
  if (nord.ne.nord_old) then
     call tri(aaH,naH,ndofH_edge)
     call tri(aaE,naE,ndofE_edge)
     nord_old = nord
  endif

  ! solve the matrix
  if (ndofH_edge.gt.0) then
     do ivar=1,MAXEQNH
        call zrhsub(aaH,zuH(1:ndofH_edge,ivar),zbH(1:ndofH_edge,ivar), &
             naH,ndofH_edge)
     enddo
  endif

  if (ndofE_edge.gt.0) then
     do ivar=1,MAXEQNE
        call zrhsub(aaE,zuE(1:ndofE_edge,ivar),zbE(1:ndofE_edge,ivar), &
             naE,ndofE_edge)
     enddo
  endif
  !---------------------------------------------
  if (iprint.eq.1) then
     write(*,*) 'dhpedge: H1, k,zu(k) = '
     do k=1,ndofH_edge
        write(*,*) k,zuH(k,1:MAXEQNH)
     enddo

     write(*,*) 'dhpedge: Hcurl, k,zu(k) = '
     do k=1,ndofE_edge
        write(*,*) k,zuE(k,1:MAXEQNE)
     enddo
     call pause
  endif
  !---------------------------------------------
  ! save the dof, skipping irrelevant entries
  call decod(Icase,2,NR_PHYSA, ncase)
  if (iprint.eq.1) then
     write(*,*) 'dhpedge: ncase = ',ncase
  endif
  !
  !---------------------------------------------
  ivarH=0; nvarH=0;  ivarE=0; nvarE=0

  ! loop through multiple copies of variables
  do j=1,NRCOMS

     ! loop through physical attributes
     do i=1,NR_PHYSA

        ! loop through components
        do k=1,NR_COMP(i)

           select case(DTYPE(i))
           case('contin')
              ivarH = ivarH + 1
              if (ncase(i).eq.1) then
                 nvarH = nvarH + 1
                 ZnodH(nvarH,1:ndofH_edge) = zuH(1:ndofH_edge,ivarH)
              endif
           case('tangen')
              ivarE = ivarE + 1
              if (ncase(i).eq.1) then
                 nvarE = nvarE + 1
                 ZnodE(nvarE,1:ndofE_edge) = zuE(1:ndofE_edge,ivarE)
              endif
           end select

        enddo
     enddo
  enddo

  if (iprint.ne.0) then
    mdle=iprint
    nrp=10; dt=1.d0/nrp
    call nodcor(mdle,xnod)
    call solelm(mdle, zdofH_c,zdofE_c,zdofV_c,zdofQ_c)
    ! loop through points along the edge
    Nflag = 1
    do i=0,nrp
      t = i*dt
      call edge_param(Type,Iedg,t, xi,dxidt)
      call shapEbe(t,nord,Nedge_orient(Iedg), &
                   ndofE_edge,shapE_edge)
      zEt = ZERO
      do j=1,ndofE_edge
        write(*,6048) j,zuE(j,1),shapE_edge(j)
6048    format('dhpedge: j,zuE(j,1),shapE_edge(j) = ',i2,2e12.5,2x,e12.5)
        zEt = zEt + zuE(j,1)*shapE_edge(j)
      enddo
      call soleval(mdle,xi, &
            Nedge_orient,Nface_orient,Norder,&
            Xnod,ZdofH_c,ZdofE_c,ZdofV_c,ZdofQ_c,Nflag, &
            x,dxdxi,&
            ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ)
      dxdt(1:3) = dxdxi(1:3,1)*dxidt(1) + dxdxi(1:3,2)*dxidt(2) &
                + dxdxi(1:3,3)*dxidt(3)

      call dirichlet( &
           x,Icase, &
           zvalH,zdvalH, &
           zvalE,zdvalE, &
           zvalV,zdvalV)

      write(*,9002) iprint,i,t,xi(1:3),x(1:3)
9002  format('mdle = ',i5,' i = ',i2,' t = ',f8.3,' xi = ',3f8.3,' x = ',3f8.3)
      write(*,9001) zvalE(1,1)*dxdt(1)+zvalE(2,1)*dxdt(2)+zvalE(3,1)*dxdt(3),&
                    zsolE(1,1)*dxdt(1)+zsolE(2,1)*dxdt(2)+zsolE(3,1)*dxdt(3), zEt
9001  format('zvalE = ',(2e12.5,2x),' zsolE = ',(2e12.5,2x),' zEt = ',2e12.5)
    enddo
    call pause
  endif

end subroutine dhpedge

