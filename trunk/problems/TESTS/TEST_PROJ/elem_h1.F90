!--------------------------------------------------------------------
!                                                                     
!     routine name      - elem_h1
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - May 10
!                                                                     
!     purpose:          - routine returns unconstrained (ordinary) 
!                         stiffness matrix and load vector  
!                         for the L2 or H1 projection problem
!     arguments:                                                     
!                                                                     
!     in:              
!             Mdle      - an element middle node number, identified
!                         with the element
!             Nrow_Zbloc,Ncol_Zaloc - dimensions
!                         of element arrays
!             Nrhs    - number of right-hand sides (loads)
!        
!     out:     
!             Zbloc     - element load vector
!             Zaloc     - element stiffness matrices
!
!-----------------------------------------------------------------------
!
subroutine elem_h1(Mdle, Zbloc,Nrow_Zbloc,Nrhs,Zaloc,Ncol_Zaloc)
!
      use parameters
      use control
      use element_data
      use data_structure3D
      use GMP
      use physics          , only : NR_COMP
!      
#include "syscom.blk"
!
      dimension Zbloc(Nrow_Zbloc,Nrhs)
      dimension Zaloc(Nrow_Zbloc,Ncol_Zaloc)
!-----------------------------------------------------------------------      
!
      character(len=4) :: type
!
!  ...element order, geometry dof
      dimension norder(19),xnod(3,MAXbrickH)
!
!  ...node orientations
      dimension nedge_orient(12), nface_orient(6)
!
!  ...shape functions and their derivatives wrt master coordinates
      dimension shapH(MAXbrickH),dshapH(3,MAXbrickH), &
                dshapHx(3,MAXbrickH)
!
!  ...geometry 
      dimension xi(3),x(3),dxdxi(3,3),dxidx(3,3)
!
!  ...3D quadrature data
      dimension xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
!
!  ...exact solution routine
      dimension                        &
                zvalH(    MAXEQNH    ),&    
                zdvalH(   MAXEQNH,3  ),&  
                zd2valH(  MAXEQNH,3,3),&
                zvalE(  3,MAXEQNE    ),&
                zdvalE( 3,MAXEQNE,3  ),&  
                zd2valE(3,MAXEQNE,3,3),&
                zvalV(  3,MAXEQNV    ),&
                zdvalV( 3,MAXEQNV,3  ),&
                zd2valV(3,MAXEQNV,3,3),&
                zvalQ(    MAXEQNQ    ),&
                zdvalQ(   MAXEQNQ,3  ),&
                zd2valQ(  MAXEQNQ,3,3)
!      
!---------------------------------------------------------------------
!
      iprint=0 
!
      ncomp=NR_COMP(1)
!
      if (iprint.ge.1) then
        write(*,7001) Mdle,NODES(Mdle)%type
 7001   format('elem_h1: Mdle,type = ',i10,2x,a5)
      endif
!      
      type = NODES(Mdle)%type
!     
!  ...determine order of approximation
      call find_order(Mdle, norder)
!
!  ...determine the element nodes and their orientation
      call find_orient(Mdle, nedge_orient,nface_orient)
!      
!  ...determine nodes coordinates 
      call nodcor(Mdle, xnod)
!
      if (iprint.eq.2) then
        call celndof(NODES(Mdle)%type,Norder, nrdofH,nrdofE,nrdofV,nrdofQ)
        do k=1,nrdofH
          write(*,7002) k,xnod(1:3,k)
 7002     format('elem_h1: k,xnod(1:3,k) = ',i3,3f8.3)
        enddo
      endif
!
!  ...clear spaces for the element matrices                    
      Zbloc = ZERO
      Zaloc = ZERO
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L
!-----------------------------------------------------------------------
!
!  ...set up the element quadrature
      call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
      if (iprint.eq.1) write(*,*) 'elem_h1: nint = ',nint 
!
!  ...loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l); wa = wxi(l)
!
!  .....evaluate shape functions
!       at the point
        call shape3H(NODES(Mdle)%type,xi,norder,nedge_orient,nface_orient, &
                     nrdofH,shapH,dshapH)
!
!  .....geometry map
        select case(EXGEOM)
        case(0)
          x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0 
          do k=1,nrdofH
            x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
            do i=1,3
              dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
            enddo
          enddo
        case(1)
          call exact_geom(Mdle,xi, x,dxdxi)
        endselect
!
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag) 
        if (iflag /= 0) then
          write(*,*) 'elem_h1: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
          write(*,*) '        rjac = ',rjac
          stop 
        endif
!
!  .....compute derivatives wrt physical coordinates
        do k=1,nrdofH
          dshapHx(1:3,k) = 0.d0
          do ixi=1,3
            dshapHx(1:3,k) = dshapHx(1:3,k)  &
                           + dshapH(ixi,k)*dxidx(ixi,1:3)
          enddo
        enddo
!
!  .....total weight
        weight = wa*rjac
!
!  .....compute the exact solution 
        icase=1
        call exact(x,icase, zvalH,zdvalH,zd2valH, &
                            zvalE,zdvalE,zd2valE, &
                            zvalV,zdvalV,zd2valV, &
                            zvalQ,zdvalQ,zd2valQ)
!
!  .....loop through test functions
        do k1=1,nrdofH
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
!  ...........accumulate for the load vector
              Zbloc(m1,irhs) = Zbloc(m1,irhs)  &
                             + zvalH(ii)*shapH(k1)*weight
              do i=1,3
                Zbloc(m1,irhs) = Zbloc(m1,irhs) &
                               + zdvalH(ii,i)*dshapHx(i,k1)*weight
              enddo
            enddo
!
!  .........loop through trial functions
            do k2=1,nrdofH
!
!  ...........2nd loop over components
              do ivar2=1,ncomp
!
                m2=(k2-1)*ncomp+ivar2
!
!  .............accumulate for the stiffness matrix....
                zs = shapH(k1)*shapH(k2)
                zs = zs &
                   + dshapHx(1,k1)*dshapHx(1,k2)  &
                   + dshapHx(2,k1)*dshapHx(2,k2)  &
                   + dshapHx(3,k1)*dshapHx(3,k2)
!
                Zaloc(m1,m2) = Zaloc(m1,m2) + zs*weight
!                
              enddo
            enddo
          enddo
        enddo
!
!  ...end of loop through integration points
      enddo
!
!-----------------------------------------------------------------------

      iprint=0
      if (iprint.ge.1) then
!              
        select case(type)
!  .....P R I S M        
        case('mdlp')
          write(*,7011) 1,18
 7011     format('elem_h1: VERTEX SHAPE FUNCTIONS = ',i2,'...',i2)
          k=18
          do ie=1,9
            ndof = (norder(ie)-1)*3
            if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
 7012       format('        EDGE ',i2,' SHAPE FUNCTIONS = ',i2,'...',i2)
            k = k+ndof
          enddo 
          do ifig=1,5
            nord = norder(9+ifig)
            select case(ifig)
            case(1,2); ndof = (nord-2)*(nord-1)/2*3
            case(3,4,5); call decode(nord, nordh,nordv)
              ndof=(nordh-1)*(nordv-1)*3
            end select
            if (k+ndof.ge.k+1) write(*,7013) ifig,k+1,k+ndof
 7013       format('        FACE ',i2,' SHAPE FUNCTIONS = ',i2,'...',i2)
            k=k+ndof
          enddo
          if (nrdofH*3.ge.k+1) write(*,7014) k+1,nrdofH*3
 7014     format('        MIDDLE NODE SHAPE FUNCTIONS = ',i2,'...',i2)
!  .....H E X A          
        case('mdlb')
!  .....T E T R A              
        case('mdln')
          write(*,7011) 1,12
          k=12
          do ie=1,6
            ndof = (norder(ie)-1)*3
            if (k+ndof.ge.k+1) write(*,7012) ie, k+1, k+ndof
            k = k+ndof
          enddo 
          do ifig=1,4
            nord = norder(6+ifig)
            ndof = (nord-2)*(nord-1)/2*3
            if (k+ndof.ge.k+1) write(*,7013) ifig,k+1,k+ndof
            k=k+ndof
          enddo
          if (nrdofH*3.ge.k+1) write(*,7014) k+1,nrdofH*3
!  .....P Y R A M I D                
        case('mdld')
        endselect
!
 50     write(*,*) 'elem_h1: Zbloc = '
        write(*,*) 'nrdofH = ',nrdofH
        write(*,*) 'elem_h1: SET jbeg,jend,l1,l2 (jbeg=0 to exit)'
        read(*,*) jbeg,jend,l1,l2
        if (jbeg.eq.0) go to 55
        do j=jbeg,jend
          write(*,7005) j,Zbloc(j,l1:l2)
        enddo
        go to 50
 55     write(*,*) 'elem_h1: Zaloc = '
        write(*,*) 'elem_h1: SET jbeg,jend,ibeg,iend (jbeg=0 to exit)'
        read(*,*) jbeg,jend,ibeg,iend
        if (jbeg.eq.0) go to 60
        do j=jbeg,jend
          write(*,7005) j,Zaloc(j,ibeg:iend)
#if C_MODE
 7005     format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
#else
 7005     format(i3,2x,10e12.5,10(/,5x,10e12.5))
#endif
        enddo
        go to 55
 60     continue
      endif
!
!
endsubroutine elem_h1
