!---------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for H(curl) projection
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
subroutine elem_hcurl(Mdle,Bloc,Nrow_Bloc,Nrhs,Aloc,Ncol_Aloc) 
!
      use control          , only : EXGEOM
      use data_structure3D , only : NODES
      use physics          , only : NR_COMP
      use parameters
!
#include"typedefs.h"
!
!-------------------------------------------------------------------------
      implicit none
      integer,                               intent(in)  :: Mdle
      integer,                               intent(in)  :: Nrow_Bloc
      integer,                               intent(in)  :: Ncol_Aloc
      integer,                               intent(in)  :: Nrhs
      VTYPE, dimension(Nrow_Bloc,Nrhs   ),   intent(out) :: Bloc
      VTYPE, dimension(Nrow_Bloc,Ncol_Aloc), intent(out) :: Aloc
!-------------------------------------------------------------------------
!     element and face type 
      character(len=4) :: etype
!    
!     element order, face order, edge and face orientations
      integer,dimension(19) :: norder
      integer,dimension(12) :: nedge_orient
      integer,dimension(6)  :: nface_orient
!
!     shape functions and their derivatives wrt master coordinates
      real*8,dimension(3,MAXbrickE) :: shapE, shapEx, curlE, curlEx
      real*8,dimension(  MAXbrickH) :: shapH
      real*8,dimension(3,MAXbrickH) :: gradH
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
!     exact solution
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
!     curl of H(curl) exact solution
      VTYPE, dimension(3,MAXEQNE    ) :: zcvalE
!
!     misc
      real*8  :: wa, weight, rjac
      integer :: i,j,k,l,k1,k2,nint,iflag,ivoid
      integer :: iprint, iprint_octave, iprint_eig
      integer :: l1, l2, ibeg, iend, jbeg, jend
      integer :: nrdofE,nrdofH
      integer :: m1,m2,ivar1,ivar2,irhs,ncomp,ii
      character(len=30) :: aux_string,aux_charnum
    !  integer :: stat_alloc
      !
      ! to use eigenvalue lapack routine
      character(len=1) :: jobz = 'N'
      character(len=1) :: uplo = 'U'
      integer :: lda, lwork_lapack, lapack_info, nnulleig, nnegeig
      double precision :: eigenvalues_lapack(Nrow_Bloc), & 
           upper_Aloc(Nrow_Bloc,Nrow_Bloc)
      double precision :: work_lapack((3*Nrow_Bloc-1)*2)
      double precision :: eps = 1e-14
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
      ncomp=NR_COMP(2)
!
!     order of approximation, orientations, nodes coordinates
      call find_elem_nodes(Mdle, norder,nedge_orient,nface_orient)
      call nodcor(         Mdle, xnod)
    
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
         xi(1:3) = xiloc(1:3,l) ; wa = wxi(l)
!
!        calculate shape function for Hcurl and H1
         call shape3E(etype,xi,norder,nedge_orient,nface_orient, nrdofE,shapE,curlE)
         call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
!
!        geometry map
         select case(EXGEOM)
         case(0)
            x(1:3)= 0.d0 ; dxdxi(1:3,1:3) = 0.d0
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
         !  ..evaluate the inverse derivatives and jacobian
         call geom(dxdxi, dxidx,rjac,iflag)
         iflag=0
         if (iflag.ne.0) then
            write(*,*) 'elem_hcurl: NEGATIVE JACOBIAN, Mdle, rjac= ',Mdle, rjac
            stop 1
         endif
         !     
         !  ..total weight
         weight = wa*rjac
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
!        get rhs terms
         ivoid=1
         call exact(x,ivoid, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                             zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ )
!
!        compute curl
         zcvalE(1,1:MAXEQNE) = zdvalE(3,1:MAXEQNE,2) - zdvalE(2,1:MAXEQNE,3)
         zcvalE(2,1:MAXEQNE) = zdvalE(1,1:MAXEQNE,3) - zdvalE(3,1:MAXEQNE,1)
         zcvalE(3,1:MAXEQNE) = zdvalE(2,1:MAXEQNE,1) - zdvalE(1,1:MAXEQNE,2)
!      
!        P.G., Nov 14 : not too sure about this...
!        shortcut to verify the exact sequence property
IF (icheck /= 0) THEN
!        E = grad H                 ; curl E = 0        
         zvalE(1:3,1)=gradH(1:3,14) ; zcvalE(1:3,1)=0.d0
ENDIF
    
         
!  .....loop over TEST functions
        do k1=1,nrdofE
!
!  .......1st loop over components
          do ivar1=1,ncomp
!
            m1=(k1-1)*ncomp+ivar1
!
!  .........loop over RHS
            do irhs=1,Nrhs
!
              ii = ncomp*(irhs-1) + ivar1
!
              do i=1,3
                Bloc(m1,irhs) = Bloc(m1,irhs) + zcvalE(i,ii)*curlEx(i,k1)*weight &
                                              +  zvalE(i,ii)*shapEx(i,k1)*weight
              enddo
            enddo
!        
!  .........loop over TRIAL functions
            do k2=1,nrdofE
!            
!  ...........2nd loop over components
              do ivar2=1,ncomp
!
                m2=(k2-1)*ncomp+ivar2
!
                do i=1,3
                  Aloc(m1,m2) = Aloc(m1,m2) + curlEx(i,k1)*curlEx(i,k2)*weight &
                                            + shapEx(i,k1)*shapEx(i,k2)*weight 
                enddo
              enddo
            enddo
          enddo
        enddo
!
!  ...end of loop over integration points
      enddo
  !
!!!  ! ----------------------------------------------------------
!!!  ! FOR DEBUGGING PURPOSES
!!!  iprint_octave=0
!!!  if (iprint_octave == 1) then
!!!     numelem=numelem+1
!!!     call int2str(numelem,aux_charnum)
!!!     
!!!     if (IERROR_PROB.eq.IERROR_CURL) then
!!!        aux_string='ALOC_CURL_octave_'//trim(aux_charnum)//'.txt'
!!!     else
!!!        aux_string='ALOC_L2_octave_'//trim(aux_charnum)//'.txt'
!!!     endif
!!!     call matrix2octave(Aloc,Nrow_Bloc,Ncol_Aloc,trim(aux_string))
!!!
!!!     if (IERROR_PROB.eq.IERROR_CURL) then
!!!        aux_string='BLOC_CURL_octave_'//trim(aux_charnum)//'.txt'
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
!!!
!!!  ! ----------------------------------------------------------
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


endsubroutine elem_hcurl
