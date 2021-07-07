!---------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for H(div) projection
!!
!! @param[in]  Mdle      - element (middle node) number 
!! @param[out] Bloc      - element load vector(s)
!! @param[in]  Nrow_Bloc - number of rows of load vector
!! @param[out] Aloc      - element stiffness matrix
!! @param[in]  Ncol_Aloc - number of columns of stiffness matrix
!!
!! @date Nov 14
!---------------------------------------------------------------------------
!
subroutine elem_hdiv(Mdle, Bloc,Nrow_Bloc,Nrhs,Aloc,Ncol_Aloc) 
!                           
      use control          , only : EXGEOM
      use data_structure3D , only : NODES
      use physics          , only : NR_COMP
      use parameters 
!
#include"typedefs.h"
!-------------------------------------------------------------------------
      implicit none
      integer,                              intent(in)  :: Mdle
      integer,                              intent(in)  :: Nrow_Bloc
      integer,                              intent(in)  :: Ncol_Aloc
      integer,                              intent(in)  :: Nrhs
      VTYPE,dimension(Nrow_Bloc,Nrhs),      intent(out) :: Bloc
      VTYPE,dimension(Nrow_Bloc,Ncol_Aloc), intent(out) :: Aloc
!-------------------------------------------------------------------------
!     element and face type 
      character(len=4) :: etype
!    
!     element order, face order, edge and face orientations
      integer,dimension(19) :: norder
      integer,dimension(12) :: nedge_orient
      integer,dimension(6)  :: nface_orient
!
!     shape functions and their derivatives 
      real*8,dimension(  MAXbrickH) :: shapH
      real*8,dimension(3,MAXbrickH) :: gradH
      real*8,dimension(3,MAXbrickE) :: shapE, shapEx, curlE, curlEx
      real*8,dimension(3,MAXbrickV) :: shapV, shapVx
      real*8,dimension(  MAXbrickV) :: divV, divVx
!
!     geometry 
      real*8,dimension(3,MAXbrickH) :: xnod
      real*8,dimension(3)           :: xi,x
      real*8,dimension(3,3)         :: dxdxi,dxidx
!
!     2D and 3D quadrature
      real*8,dimension(3,MAX_NINT3) :: xiloc
      real*8,dimension(  MAX_NINT3) :: wxi
!
!     exact solution routine
      VTYPE, dimension(  MAXEQNH    ) ::   zvalH 
      VTYPE, dimension(  MAXEQNH,3  ) ::  zdvalH
      VTYPE, dimension(  MAXEQNH,3,3) :: zd2valH
      VTYPE, dimension(3,MAXEQNE    ) ::   zvalE
      VTYPE, dimension(3,MAXEQNE,3  ) ::  zdvalE
      VTYPE, dimension(3,MAXEQNE,3,3) :: zd2valE
      VTYPE, dimension(3,MAXEQNV    ) ::   zvalV
      VTYPE, dimension(3,MAXEQNV,3  ) ::  zdvalV
      VTYPE, dimension(3,MAXEQNV,3,3) :: zd2valV
      VTYPE, dimension(  MAXEQNQ    ) ::   zvalQ
      VTYPE, dimension(  MAXEQNQ,3  ) ::  zdvalQ
      VTYPE, dimension(  MAXEQNQ,3,3) :: zd2valQ
!      
!     divergence of H(div) exact solution
      VTYPE, dimension(  MAXEQNV    ) :: zdiv_valV
!
!     misc
      real*8  :: wa, weight, rjac
      integer :: i,j,k,l,k1,k2,nint,iflag,ivoid
      integer :: iprint, iprint_octave, iprint_eig
      integer :: l1, l2, ibeg, iend, jbeg, jend
      integer :: nrdofV,nrdofE,nrdofH
      integer :: ncomp,irhs,ii,ivar1,ivar2,m1,m2
!      
!     to use eigenvalue lapack routine
      character(len=1) :: jobz = 'N'
      character(len=1) :: uplo = 'U'
      integer :: lda, lwork_lapack, lapack_info, nnulleig, nnegeig
      real*8 :: eigenvalues_lapack(Nrow_Bloc), & 
           upper_Aloc(Nrow_Bloc,Nrow_Bloc)
      real*8 :: work_lapack((3*Nrow_Bloc-1)*2)
      real*8 :: eps = 1e-14
!
      integer :: icheck
!
!-------------------------------------------------------------------------
!
!     check exact sequence (0 - No ; 1 - Yes)
      icheck=0
!      
!     printing flag
      iprint=0
!
!     number of components of physical attribute
      ncomp=NR_COMP(3)
!
!     order of approximation, orientations, nodes coordinates
      call find_elem_nodes(Mdle, norder,nedge_orient,nface_orient)
      call nodcor(         Mdle, xnod)
!    
!     clear spaces for stiffness matrix and load vector
      Aloc=ZERO ; Bloc=ZERO
!
!-------------------------------------------------------------------------
! E L E M E N T    I N T E G R A L S                                      
!-------------------------------------------------------------------------
!
!     collect integration points
      etype = NODES(Mdle)%type
      call set_3Dint(etype,norder, nint,xiloc,wxi)
!
!     loop over integration points
      do l=1,nint
!      
        xi(1:3) = xiloc(1:3,l) ; wa = wxi(l)
!
!       calculate shape function for H1, [ H(curl) ], H(div)
        call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
IF (icheck /= 0) THEN
        call shape3E(etype,xi,norder,nedge_orient,nface_orient, nrdofE,shapE,curlE)
ENDIF        
        call shape3V(etype,xi,norder,             nface_orient, nrdofV,shapV,divV )
!
!       geometry map
        select case(EXGEOM)
        case(0)
           x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0
           do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
              do i=1,3
                 dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*gradH(i,k)
              enddo
           enddo
        case(1)
           call exact_geom(Mdle,Xi, x,dxdxi)
        endselect
!
!       evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag)
        iflag=0
        if (iflag.ne.0) then
           write(*,*) 'elem_hdiv: NEGATIVE JACOBIAN, Mdle, rjac= ',Mdle, rjac
           stop
        endif
!     
!       total weight
        weight = wa*rjac
!
!       evaluate the shapV and divV wrt physical coordinates (Piola transform)
IF (icheck /= 0) THEN
        curlEx=0.d0
ENDIF        
        shapVx=0.d0 ; divVx=0.d0
        do k=1,nrdofV
           do i=1,3
              do j=1,3
IF (icheck /= 0) THEN
                 curlEx(i,k) = curlEx(i,k) + dxdxi(i,j)*curlE(j,k)/rjac
ENDIF        
                 shapVx(i,k) = shapVx(i,k) + dxdxi(i,j)*shapV(j,k)/rjac
              enddo
           enddo
           divVx(k) = divVx(k) + divV(k)/rjac
        enddo
!
!       get rhs terms
        ivoid=1
        call exact(x,ivoid, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ )
!             
!       divergence of the exact solution
        zdiv_valV(1:MAXEQNV) = zdvalV(1,1:MAXEQNV,1) + zdvalV(2,1:MAXEQNV,2) + zdvalV(3,1:MAXEQNV,3)
!
        if (iprint.eq.1) then
          write(*,7011) l,zvalV(1:3,1),zdiv_valV(1)
 7011     format('elem: l,zvalV(1:3,1) = ',i3,2x,3e12.5,' zdiv = ',e12.5)
          call pause
        endif
!
!       P.G., Nov 14 : not too sure about this...
!       shortcut to verify the exact sequence property
IF (icheck /= 0) THEN
!       V = curl E                 ; div V = 0
        zvalV(1:3,1)=curlEx(1:3,1) ; zdiv_valV(1)=0.d0
ENDIF
!
!       loop over TEST functions
        do k1=1,nrdofV
!     
!         1st loop over components
          do ivar1=1,ncomp
!
            m1=(k1-1)*ncomp+ivar1
!
!           loop over RHS
            do irhs=1,Nrhs
!
              ii = ncomp*(irhs-1) + ivar1
!
!             seminorm
              Bloc(m1,irhs) = & 
              Bloc(m1,irhs) + zdiv_valV(ii)*divVx(k1)*weight
!
!             L2 norm   
              do i=1,3
              Bloc(m1,irhs) = &
              Bloc(m1,irhs) + zvalV(i,ii)*shapVx(i,k1)*weight
              enddo
            enddo
!          
!           loop over TRIAL functions
            do k2=1,nrdofV
!
!             2nd loop over components
              do ivar2=1,ncomp
!
                m2=(k2-1)*ncomp+ivar2
!                
!               seminorm
                Aloc(m1,m2) = Aloc(m1,m2) + divVx(k1)*divVx(k2)*weight 
!
!               L2 norm
                do i=1,3
                  Aloc(m1,m2) = Aloc(m1,m2) + shapVx(i,k1)*shapVx(i,k2)*weight   
                enddo
!               
              enddo
            enddo
          enddo
        enddo
!
!     loop over integration points
      enddo
!
!!!  ! ----------------------------------------------------------
!!!  ! FOR DEBUGGING PURPOSES
!!!  iprint_octave=0
!!!  if (iprint_octave == 1) then
!!!     numelem=numelem+1
!!!     call int2str(numelem,aux_charnum)
!!!     
!!!     if (IERROR_PROB.eq.IERROR_DIV) then
!!!        aux_string='ALOC_DIV_octave_'//trim(aux_charnum)//'.txt'
!!!     else
!!!        aux_string='ALOC_L2_octave_'//trim(aux_charnum)//'.txt'
!!!     endif
!!!     call matrix2octave(Aloc,Nrow_Bloc,Ncol_Aloc,trim(aux_string))
!!!
!!!     if (IERROR_PROB.eq.IERROR_DIV) then
!!!        aux_string='BLOC_DIV_octave_'//trim(aux_charnum)//'.txt'
!!!     else
!!!        aux_string='BLOC_L2_octave_'//trim(aux_charnum)//'.txt'
!!!     endif
!!!     call matrix2octave(Bloc,Nrow_Bloc,1,trim(aux_string))
!!!
!!!  endif
!!!  ! ----------------------------------------------------------
!!!  
!!!  ! ----------------------------------------------------------
!!!  ! FOR DEBUGGING PURPOSES
!!!  iprint_eig=0
!!!  if (iprint_eig == 1) then
!!!     lda=Nrow_Bloc
!!!     upper_Aloc=ALOC
!!!     lwork_lapack=size(work_lapack)
!!!     call DSYEV( jobz, uplo, Nrow_Bloc, upper_Aloc, lda, eigenvalues_lapack, &
!!!          work_lapack, lwork_lapack, lapack_info )
!!!     if (lapack_info /= 0) then
!!!        write(*,*) 'lapack INFO: ', lapack_info
!!!        call pause
!!!     endif
!!!     write(*,*) 'Eigenvalues of ALOC:'
!!!!!$     do l=1,size(eigenvalues_lapack)
!!!!!$        write(*,*) eigenvalues_lapack(l)
!!!!!$     enddo
!!!     ! calculate number of numerically zero eigenvalues
!!!     nnulleig=COUNT(abs(eigenvalues_lapack) < eps)
!!!     write(*,*) 'Number of null eigenvalues of ALOC: ', nnulleig
!!!     ! calculate number of negative eigenvalues
!!!     nnegeig=COUNT(eigenvalues_lapack < 0.d0)
!!!     write(*,*) 'Number of negative eigenvalues of ALOC: ', nnegeig
!!!!!$     pause
!!!  endif

  ! ----------------------------------------------------------
  !
  iprint=0
  if (iprint.eq.-1) then
     write(*,*) 'elem_hdiv: Mdle = ', Mdle
50   write(*,*) 'elem_hdiv: Zbloc, nrdofV = ', nrdofV
     write(*,*) 'elem_hdiv: SET jbeg,jend,l1,l2(jbeg=0 to exit)'
     read(*,*) jbeg,jend,l1,l2
     if (jbeg.eq.0) go to 55
     do j=jbeg,jend
        write(*,7005) j,Bloc(j,l1:l2)
     enddo
     go to 50
55   write(*,*) 'elem_hdiv: Zaloc, nrdofV = ', nrdofV
     write(*,*) 'elem_hdiv: SET jbeg,jend,ibeg,iend(jbeg=0 to exit)'
     read(*,*) jbeg,jend,ibeg,iend
     if (jbeg.eq.0) go to 60
     do j=jbeg,jend
        write(*,7005) j,Aloc(j,ibeg:iend)
7005    format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
     enddo
     go to 55
60   continue
  endif


endsubroutine elem_hdiv
