!--------------------------------------------------------------------
!     
!     routine name      - elem_hcurl
!     
!--------------------------------------------------------------------
!     
!     latest revision:  - Jan 13
!     
!     purpose:          - routine returns unconstrained (ordinary) 
!                         stiffness matrix and load vector  
!                         for H(curl) projection
!     arguments:                                                     
!     
!     in:              
!           Mdle       - an element middle node number, identified
!                        with the element
!           Nrow_Zbloc,
!           Ncol_Zaloc - dimensions of element arrays
!           Nr_RHS     - number of right-hand sides (loads)
!     
!     out:     
!           Zbloc      - element load vectors
!           Zaloc      - element stiffness matrices
!     
!-----------------------------------------------------------------------
!     
subroutine elem_hcurl_back(Mdle, Zbloc,Nrow_Zbloc,Nr_RHS, &
                            Zaloc,Ncol_Zaloc          )
!     
      use control
      use element_data
      use data_structure3D
      use GMP
      use PROJ
#include "syscom.blk"
      dimension Zbloc(Nrow_Zbloc,Nr_RHS)
      dimension Zaloc(Nrow_Zbloc,Ncol_Zaloc)
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
      dimension shapE(3,MAXbrickE),curlE(3,MAXbrickE), &
                shapEx(3,MAXbrickE),curlEx(3,MAXbrickE)
      dimension shapH(MAXbrickH),dshapH(3,MAXbrickH), &
                dshapHx(3,MAXbrickH)
!     
!  ...geometry 
      dimension xi(3),x(3),dxdxi(3,3),dxidx(3,3)
!     
!  ...3D quadrature data
      dimension xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
!     
!  ...metrics, rhs
      dimension aa(3,3),bb(3,3),zc(3),zd(3)
!     
!  ...exact solution routine
      dimension zvalH(  MAXEQNH),zdvalH(  MAXEQNH,3),zd2valH(  MAXEQNH,3,3),  &
                zvalE(3,MAXEQNE),zdvalE(3,MAXEQNE,3),zd2valE(3,MAXEQNE,3,3),  &
                zvalV(3,MAXEQNV),zdvalV(3,MAXEQNV,3),zd2valV(3,MAXEQNV,3,3),  &
                zvalQ(  MAXEQNQ),zdvalQ(  MAXEQNQ,3),zd2valQ(  MAXEQNQ,3,3)
!
!  ...curl of exact solution
      dimension zcvalE(3,MAXEQNE)
!     
!  ...for activating L2 projection only
      logical :: L2PROJ=.FALSE.
!
!---------------------------------------------------------------------
!     
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  ...HACK FOR PYRAMID
      if ( (NODES(Mdle)%type.eq.'pyra').or.         &
           (NODES(Mdle)%type.eq.'mdld')     ) then
        Zbloc=ZERO ; Zaloc=ZERO
        do i=1,Ncol_Zaloc ; Zaloc(i,i)=ZONE ; enddo
        return
      endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      iprint=0 
!     
      if (iprint.eq.1) then
        write(*,7001) Mdle,NODES(Mdle)%type
 7001   format(' elem_hcurl: Mdle,type = ',i10,2x,a5)
      endif
!     
!  ...determine order of approximation
      call find_order(Mdle, norder)
!     
!  ...determine the node orientation
      call find_orient(Mdle, nedge_orient,nface_orient)
!     
!  ...determine nodes coordinates 
      call nodcor(Mdle, xnod)

      if (iprint.eq.2) then
        call celndof(NODES(Mdle)%type,Norder, &
          nrdofH,nrdofE,nrdofV,nrdofQ)
        do k=1,nrdofH
          write(*,7002) k,xnod(1:3,k)
 7002     format(' elem_hcurl: k,xnod(1:3,k) = ',i3,3f8.3)
        enddo
      endif
!     
!     ...clear spaces for the element matrices                    
      Zbloc = ZERO
      Zaloc = ZERO
!     
!-----------------------------------------------------------------------
!     
!  ...compute element integrals
!     
!  ...set up the element quadrature
      call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
      if (iprint.eq.1) write(*,*) 'elem_hcurl: nint = ',nint 
!     
!  ...loop through integration points
      do l=1,nint
        xi(1:3) = xiloc(1:3,l); wa = wxi(l)
!     
!  .....compute appropriate shape functions at integration point and
!       set up appropriate number of dof's        
        call shape3E(NODES(Mdle)%type,xi,norder, &
          nedge_orient,nface_orient, &
          nrdofE,shapE,curlE)

!!!          do i=1,nrdofE
!!!            write(*,9999)i,shapE(1:3,i),curlE(1:3,i)
!!! 9999       format(' i,shapE,curlE = ',i3,2x,3(e12.5,2x),2x,3(e12.5,2x))           
!!!          enddo

!     
!  .....evaluate H^1 shape functions at integration point
        call shape3H(NODES(Mdle)%type,xi,norder, &
          nedge_orient,nface_orient, &
          nrdofH,shapH,dshapH)
!     
!  .....determine physical coordinates and the derivatives of 
!       the physical coordinates wrt master element coordinates
        x(1:3) = 0.d0; dxdxi(1:3,1:3) = 0.d0 
        do k=1,nrdofH
          x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dshapH(i,k)
          enddo
        enddo
!     
!  .....evaluate the inverse derivatives and jacobian
        call geom(dxdxi, dxidx,rjac,iflag) 
        if (iflag.ne.0) then
          write(*,*) 'elem_A: NEGATIVE JACOBIAN FOR Mdle = ',Mdle
          write(*,*) '        rjac = ',rjac
          stop 1
        endif
!
!  .....evaluate the shapE and curlE wrt physical coordinate
        shapEx = ZERO; curlEx = ZERO;
        do k=1,nrdofE
          do i=1,3
            do j=1,3
              shapEx(i,k) = shapEx(i,k)+ &
                            (dxidx(j,i)*shapE(j,k))
              curlEx(i,k) = curlEx(i,k)+ &
                            (dxdxi(i,j)*curlE(j,k))/rjac
            enddo
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
        do k1=1,nrdofE
!     
!  .......curl of exact solution
          zcvalE(1,1:MAXEQNE)=zdvalE(3,1:MAXEQNE,2)-zdvalE(2,1:MAXEQNE,3)
          zcvalE(2,1:MAXEQNE)=zdvalE(1,1:MAXEQNE,3)-zdvalE(3,1:MAXEQNE,1)
          zcvalE(3,1:MAXEQNE)=zdvalE(2,1:MAXEQNE,1)-zdvalE(1,1:MAXEQNE,2)
!     
!  .......loop over equations
          do n=1,MAXEQNE
!
!  .........accumulate              
            do i=1,3       
              Zbloc(k1,n) = Zbloc(k1,n) &
                          + zvalE(i,n)*shapEx(i,k1)*weight
              if (.NOT. L2PROJ) then
                Zbloc(k1,n) = Zbloc(k1,n) &
                            + zcvalE(i,n)*curlEx(i,k1)*weight
              endif
            enddo
!     
!  .......end of loop over equations            
          enddo
!     
!  .......loop through trial functions
          do k2=1,nrdofE
!     
!  .........loop through components
            do i=1,3
              Zaloc(k1,k2) = Zaloc(k1,k2) &
                           + shapEx(i,k1)*shapEx(i,k2)*weight
              if (.NOT. L2PROJ) then
                Zaloc(k1,k2) = Zaloc(k1,k2) &
                             + curlEx(i,k1)*curlEx(i,k2)*weight
              endif

                if ((k1.eq.30).and.(k2.eq.1)) then
!!!                  write(*,*)'Zaloc = ',zaloc(k1,k2)
                endif

            enddo
!     
!  .......end of loop through trial functions
          enddo
!     
!  .....end of loop through test functions
        enddo
!     
!  ...end of loop through integration points
      enddo
!     
!-----------------------------------------------------------------------
!     
!     
      iprint=0
      if (iprint.eq.1) then
!     
!     
 50     write(*,*) 'elem_hcurl: Zbloc = '
        write(*,*) 'elem_hcurl: SET jbeg,jend,l1,l2(jbeg=0 to exit)'
        read(*,*) jbeg,jend,l1,l2
        if (jbeg.eq.0) go to 55
        do j=jbeg,jend
          write(*,7005) j,Zbloc(j,l1:l2)
        enddo
        go to 50
 55     write(*,*) 'elem_hcurl: Zaloc = '
        write(*,*) 'elem_hcurl: SET jbeg,jend,ibeg,iend(jbeg=0 to exit)'
        read(*,*) jbeg,jend,ibeg,iend
        if (jbeg.eq.0) go to 60
        do j=jbeg,jend
          write(*,7005) j,Zaloc(j,ibeg:iend)
#if C_MODE
 7005     format(i3,2x,5(2e11.4,2x),10(/,5x,5(2e11.4,2x)))
#else
 7005     format(i3,2x,10(e12.5,2x),10(/,5x,10e12.5))
#endif
        enddo
        go to 55
 60     continue
      endif
!     
!     
endsubroutine elem_hcurl_back
