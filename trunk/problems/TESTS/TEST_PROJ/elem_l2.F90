!------------------------------------------------------------------------------------
!> Purpose : element stiffness matrix and load vector for elasticity problem 
!! @param[in]  Mdle      - middle node number
!! @param[out] Bloc      - elem load vector(s)
!! @param[in]  Nrow_Bloc - number of rows of load vector
!! @param[out] Aloc      - elem stiffness matrix
!! @param[in]  Ncol_Aloc - number of columns of stiffness matrix
!------------------------------------------------------------------------------------
!
subroutine elem_l2(Mdle,Bloc,Nrow_Bloc,Nr_RHS,Aloc,Ncol_Aloc)
!
      use data_structure3D
      use element_data
      use control  , only : EXGEOM
      use physics  , only : NR_COMP
!
#include"typedefs.h"
!
!------------------------------------------------------------------------------------
      implicit none
      integer,                              intent(in)    :: Mdle
      integer,                              intent(in)    :: Nrow_Bloc,Ncol_Aloc
      integer,                              intent(in)    :: Nr_RHS
      VTYPE, dimension(Nrow_Bloc,Nr_RHS),   intent(inout) :: Bloc
      VTYPE, dimension(Nrow_Bloc,Ncol_Aloc),intent(inout) :: Aloc
!------------------------------------------------------------------------------------
!
!  ...element type
      character(len=4) :: etype
!
!  ...element order, face order, edge and face orientations
      integer,dimension(19) :: norder
      integer,dimension(12) :: nedge_orient
      integer,dimension(6)  :: nface_orient
!
!  ...shape functions and their derivatives 
      real*8,dimension(  MAXbrickH) :: shapH
      real*8,dimension(3,MAXbrickH) :: gradH
      real*8,dimension(3,MAXbrickV) :: shapV
      real*8,dimension(  MAXbrickV) :: divV,divVx
      real*8,dimension(  MAXbrickQ) :: shapQ
!
!  ...geometry 
      real*8,dimension(3,MAXbrickH) :: xnod
      real*8,dimension(3)           :: xi,x
      real*8,dimension(3,3)         :: dxdxi,dxidx
!
!  ...3D quadrature data
      real*8,dimension(3,MAX_NINT3) :: xiloc
      real*8,dimension(  MAX_NINT3) :: wxi
!
!  ...exact solution routine
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
!  ...miscellaneous
      integer :: icase,iflag 
      integer :: i,j,k,l, k1,k2, nint, iprint
      integer :: ibeg,iend,jbeg,jend,l1,l2
      integer :: nrdofH,nrdofV,nrdofQ
      real*8 ::  rjac,wa,weight
      integer :: ncomp,ivar1,ivar2,m1,m2,irhs,ii
!
!------------------------------------------------------------------------------------
!     
      iprint = 0
!
!  ...number of components of attribute
      ncomp=NR_COMP(4)
!
!  ...order of approximation, element nodes, orientations
      call find_elem_nodes(Mdle, norder, nedge_orient,nface_orient)
      call nodcor(Mdle, xnod)
!
!  ...clear spaces for the element matrices                    
      Bloc=ZERO ; Aloc=ZERO
!
!------------------------------------------------------------------------------------
!     E L E M E N T   I N T E G R A L S
!------------------------------------------------------------------------------------
!
!  ...set up the element quadrature
!     REMARK : We use quadrature for order p (although L2 is exactly of order p-1), 
!     in order to use same integration points as geometry mapping H1.
      etype=NODES(Mdle)%type
      call set_3Dint(etype,norder, nint,xiloc,wxi)
!
!  ...loop through integration points
      do l=1,nint
        xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!  .....evaluate appropriate shape functions at the point
        call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
        call shape3V(etype,xi,norder,             nface_orient, nrdofV,shapV,divV )
        call shape3Q(etype,xi,norder,                           nrdofQ,shapQ      )
!
!  .....geometry map
        select case(EXGEOM)
        case(0)
           x(1:3) = 0.d0 ; dxdxi(1:3,1:3) = 0.d0
           do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
              do i=1,3
                 dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*gradH(i,k)
              enddo
           enddo
        case(1)
           call exact_geom(Mdle,xi, x,dxdxi)
        endselect
!
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac, iflag) 
        if (iflag.ne.0) then
           write(*,*) 'elem_proj: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
           write(*,*) '        rjac = ',rjac
           stop
        endif
!
!
!  .....piola transform for L2 is simply to divide by "rjac"
        divVx(1:nrdofV) = divV(1:nrdofV)/rjac
        shapQ(1:nrdofQ)=shapQ(1:nrdofQ)/rjac
!
!  .....total weight
        weight = wa*rjac
!
!  .....get rhs terms
        icase=1
        call exact(x,icase, &
             zvalH,zdvalH,zd2valH, &
             zvalE,zdvalE,zd2valE, &
             zvalV,zdvalV,zd2valV, &
             zvalQ,zdvalQ,zd2valQ)
!
!  ...shortcut to verify the exact sequence property
!!!    zvalQ(1) = divVx(1)
!      
!  .....loop over test functions   
        do k1=1,nrdofQ
!
!  .......1st loop over components
          do ivar1=1,ncomp
!
            m1=(k1-1)*ncomp+ivar1
!            
!  .........loop over RHS
            do irhs=1,Nr_RHS
!
              ii = ncomp*(irhs-1) + ivar1
!              
              Bloc(m1,irhs) = &
              Bloc(m1,irhs) + zvalQ(ii)*shapQ(k1)*weight
            enddo
!
!  .........loop over trial functions
            do k2=1,nrdofQ
!
!  ...........2nd loop over components
              do ivar2=1,ncomp
!
                m2=(k2-1)*ncomp+ivar2
                Aloc(m1,m2) = Aloc(m1,m2) + shapQ(k1)*shapQ(k2)*weight
              enddo 
            enddo 
          enddo
        enddo
!   
!  ...loop over integration points
      enddo 
     !

!!!  ! ----------------------------------------------------------
!!!  ! FOR DEBUGGING PURPOSES
!!!  numelem=numelem+1
!!!  call int2str(numelem,aux_charnum)
!!!
!!!  aux_string='ALOC_octave_'//trim(aux_charnum)//'.txt'
!!!  call matrix2octave(Aloc,Nrow_Bloc,Ncol_Aloc,trim(aux_string))
!!!  aux_string='BLOC_octave_'//trim(aux_charnum)//'.txt'
!!!  call matrix2octave(Bloc,Nrow_Bloc,1,trim(aux_string))
!!!  ! ----------------------------------------------------------
 
  iprint=0
  if (iprint.eq.-1) then
     !
     !
50   write(*,*) 'elem_l2: Bloc, nrdofQ = ',nrdofQ
     write(*,*) 'elem_l2: SET jbeg,jend,l1,l2'
     read(*,*) jbeg,jend,l1,l2
     if (jbeg.eq.0) go to 55
     do j=jbeg,jend
        write(*,7005) j,Bloc(j,l1:l2)
     enddo
     go to 50
55   write(*,*) 'elem_l2: Aloc = '
     write(*,*) 'elem_l2: SET jbeg,jend,ibeg,iend'
     read(*,*) jbeg,jend,ibeg,iend
     if (jbeg.eq.0) go to 60
     do j=jbeg,jend
        write(*,7005) j,Aloc(j,ibeg:iend)
#if C_MODE
7005    format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
#else
7005    format(i3,2x,10e12.5,10(/,5x,10e12.5))
#endif
     enddo
     go to 55
60   continue
  endif
  !
  !
end subroutine elem_l2
