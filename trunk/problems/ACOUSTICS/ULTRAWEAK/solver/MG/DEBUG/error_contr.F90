!--------------------------------------------------------------------
!                                                                     
!     routine name      - error_contr
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Feb 2018
!                                                                     
!     purpose:            - routine evaluates L^2 norm of the difference
!                           between the coarse and fine-grid solutions
!                           over the fine grid element boundary
!                         
!     arguments:                                                     
!                                                                     
!     in:     
!           MdleC         - coarse grid element         
!           Mdle          - fine grid element descendant of MdleC
!           Nflag         - vector of length NR_PHYSA indicating for which
!                           attribute the error should be computed.
!           NorderC       - element order of approximation 
!           Nedge_orientC - orientation of element edges
!           Nface_orientC - orientation of element face
!           XnodC         - element geometry dof
!           ZdofH_C       - coarse element H1      dof
!           ZdofE_C       - coarse element H(curl) dof
!           ZdofV_C       - coarse element H(div)  dof
!     out:              
!           DerrorH       - H1      contribution to the element error (squared)
!           DerrorE       - H(curl) contribution to the element error (squared)
!           DerrorV       - H(div)  contribution to the element error (squared)
!
!---------------------------------------------------------------------
!
   subroutine error_contr(MdleC,Mdle,Nflag,NorderC,Nedge_orientC,Nface_orientC,  &
                          XnodC,ZdofH_C,ZdofE_C,ZdofV_C,ZdofQ_C,                 &
                          ErrorH,ErrorE,ErrorV)
!
   use control,          only : INTEGRATION
   use data_structure3D
   use parametersDPG
   use physics
!
   IMPLICIT NONE
!
   integer,    intent(in ) :: MdleC, mdle
   integer,    intent(in ) :: Nflag(NR_PHYSA)
   integer,    intent(in ) :: norderC(19),Nedge_orientC(12), Nface_orientC(6)
   real*8,     intent(in ) :: XnodC(3,MAXbrickH)
#if C_MODE
   complex*16, intent(in ) :: ZdofH_C(MAXEQNH,MAXbrickH),    &
                              ZdofE_C(MAXEQNE,MAXbrickE),    &
                              ZdofV_C(MAXEQNV,MAXbrickV),    &
                              ZdofQ_C(MAXEQNQ,MAXbrickQ)

#else       
   real*8,     intent(in ) :: ZdofH_C(MAXEQNH,MAXbrickH),    &
                              ZdofE_C(MAXEQNE,MAXbrickE),    &
                              ZdofV_C(MAXEQNV,MAXbrickV),    &   
                              ZdofQ_C(MAXEQNQ,MAXbrickQ)   
#endif           
   real*8,     intent(out) :: ErrorH(MAXEQNH),ErrorE(MAXEQNE),ErrorV(MAXEQNV)
!
!..locals

!
!..boundary adjacency flags for faces
   integer  :: nadj(6), iflag
!
!..element order and orientation
   integer  :: norder(19),nedge_orient(12), nface_orient(6), norderf(5)
!
!..relative vertex coordinates wrt MdleC
   real*8   :: xsub(3,8)
!
!..geometry dof
   real*8   :: xnod(3,MAXbrickH) 

!..fine grid approximate solution dof
#if C_MODE
   complex*16 :: zdofH(MAXEQNH,MAXbrickH), zdofE(MAXEQNE,MAXbrickE),   &
                 zdofV(MAXEQNV,MAXbrickV), zdofQ(MAXEQNQ,MAXbrickQ)
#else
   real*8     :: zdofH(MAXEQNH,MAXbrickH), zdofE(MAXEQNE,MAXbrickE),   &
                 zdofV(MAXEQNV,MAXbrickV), zdofQ(MAXEQNQ,MAXbrickQ)
#endif                 
!       
!..geometry
   real*8 :: xi(3),x(3),dxdxi(3,3),dxidx(3,3), rjac
!
!..coarse grid coordinates
   real*8 :: xiC(3)
!
!
#if C_MODE
!..approximate fine grid solution 
   complex*16 ::  zsolH(MAXEQNH),   zdsolH(MAXEQNH,3),   &
                  zsolE(MAXEQNE,3), zcurlE(MAXEQNE),     &
                  zsolV(MAXEQNV,3), zdivV(MAXEQNV),      &
                  zsolQ(MAXEQNQ)
!..approximate coarse grid solution 
   complex*16 ::  zsolH_C(MAXEQNH),   zdsolH_C(MAXEQNH,3),           &
                  zsolE_C(MAXEQNE,3), zcurlE_C(MAXEQNE), zsolEn(3),  &
                  zsolV_C(MAXEQNV,3), zdivV_C(MAXEQNV),  zsolVn,     &
                  zsolQ_C(MAXEQNQ)
#else
!..approximate fine grid solution 
   real*8     ::  zsolH(MAXEQNH),   zdsolH(MAXEQNH,3),           &
                  zsolE(MAXEQNE,3), zcurlE(MAXEQNE), zsolEn(3),  &
                  zsolV(MAXEQNV,3), zdivV(MAXEQNV),  zsolVn,     &
                  zsolQ(MAXEQNQ)
!..approximate coarse grid solution 
   real*8     ::  zsolH_C(MAXEQNH),   zdsolH_C(MAXEQNH,3),   &
                  zsolE_C(MAXEQNE,3), zcurlE_C(MAXEQNE),     &
                  zsolV_C(MAXEQNV,3), zdivV_C(MAXEQNV),      &
                  zsolQ_C(MAXEQNQ)
#endif
!
   character(len=4) :: type, ftype
! 
   integer :: norderl(19) 
   integer :: norientl(27) 
!
   integer :: nrv, nre, nrf
!..node physics attributes flags
   integer :: ncase(NR_PHYSA), iprint,j, icase
   integer :: nvarH, nvarE, nvarV,nord,nsign,if, nint
!..2D quadrature data
   real*8  :: tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!..integration loop
   real*8  :: shapH(MAXbrickH), gradH(3,MAXbrickH)
   real*8  :: rn(3), t(2), dxidt(3,2), dxdt(3,2), bjac, weight
   integer :: l,icomp, iattr, ivarH, ivarE, ivarV, k, iload, nrdofH
!
!---------------------------------------------------------------------
!
   type = NODES(mdle)%type
   xi = 0.0d0; xnod = 0.0d0
   select case(mdle)
   case(118)
      iprint = 0
   case default
      iprint = 0
   end select
   if (iprint.eq.1) then
      write(*,7001) MdleC,Mdle,Nflag(1:NR_PHYSA)
 7001 format('error_contr: MdleC,Mdle = ',2i5,' Nflag = ',10i2)
   endif
!
!..initialize global quantities
   ErrorH=0.d0; ErrorE=0.d0; ErrorV=0.d0
!
!..determine relative vertex coordinates of the element
   call locate_coordC(MdleC,Mdle, xsub)
   if (iprint.eq.2) then
      write(*,7002) (xsub(1:3,j),j=1,8)
 7002 format('xsub = ', 3f8.3)
      call pause
   endif
!
!..check if faces adjacent to the boundary (hexas only for now)
   nadj=1
   select case(type)
   case('mdlb')
!  ...bottom face      
      if (xsub(3,1) .gt. 0.d0 .or. xsub(3,2) .gt. 0.d0 .or.     &
          xsub(3,3) .gt. 0.d0 .or. xsub(3,4) .gt. 0.d0  )  nadj(1)=0
!  ...top face   
      if (xsub(3,5) .lt. 1.d0 .or. xsub(3,6) .lt. 1.d0 .or.     &
          xsub(3,7) .lt. 1.d0 .or. xsub(3,8) .lt. 1.d0  )  nadj(2)=0
!  ...left face
      if (xsub(2,1) .gt. 0.d0 .or. xsub(2,2) .gt. 0.d0 .or.     &
          xsub(2,5) .gt. 0.d0 .or. xsub(2,6) .gt. 0.d0  )  nadj(3)=0
!  ...front face
      if (xsub(1,2) .lt. 1.d0 .or. xsub(1,3) .lt. 1.d0 .or.     &
          xsub(1,6) .lt. 1.d0 .or. xsub(1,7) .lt. 1.d0  )  nadj(4)=0
!  ...right face
      if (xsub(2,3) .lt. 1.d0 .or. xsub(2,4) .lt. 1.d0 .or.     &
          xsub(2,7) .lt. 1.d0 .or. xsub(2,8) .lt. 1.d0  )  nadj(5)=0
!  ...back face
      if (xsub(1,1) .gt. 0.d0 .or. xsub(1,4) .gt. 0.d0 .or.     &
          xsub(1,5) .gt. 0.d0 .or. xsub(1,8) .gt. 0.d0  )  nadj(6)=0   
   case default
      write(*,*) 'error_contr: NOT FINISHED. type = ', type 
      stop 1      
   end select
!..quit if the element is not adjacent to the boundary
   if (sum(nadj(1:6)).eq.0) return
!
!..determine order and orientation for nodes
   nrv = nvert(type); nre = nedge(type); nrf = nface(type)
   call find_order( Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..determine element dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..decode the physical attributes of the node
   call decod(NODES(Mdle)%case,2,NR_PHYSA, ncase)

   norderl(1:19) = (/1,1,1,1,1,1,1,1,1,1,1,1, 11,11,11,11,11,11, 111/)
!
!..identify number of components
   icase = NODES(Mdle)%case
   nvarH = NREQNH(icase)
   nvarE = NREQNE(icase)
   nvarV = NREQNV(icase)
!
   ! write(*,*) 'mdle = ', mdle
!..loop through the element faces
   do if=1,nrf
      if (nadj(if).eq.0) cycle
!
!  ...normal sign corresponding to face parametrization
      nsign = nsign_param(type,if)
!
!  ...face type
      ftype = face_type(type,if)
!
!  ...face order of approximation
      call face_order(type,if,norder, norderf)
!
!  ...set 2D quadrature
      INTEGRATION = 2
      call set_2Dint_DPG(ftype,norderf, nint,tloc,wtloc)
      INTEGRATION = 0
!
!  ...loop through integration points
      do l=1,nint
!
!     ...face coordinates
         t(1:2) = tloc(1:2,l)
!         
!     ...face parametrization
         call face_param(type,if,t, xi,dxidt)
!
!     ...coarse grid element coordinates
         xiC(1:3) = 0.d0

         call shape3H(type,xi,norderl,nedge_orient,nface_orient,nrdofH,shapH,gradH)
         do k=1,nrdofH
            xiC(1:3) = xiC(1:3) + xsub(1:3,k)*shapH(k)
         enddo
!
!     ...evaluate the fine approximate solution at the point 
         iflag = 1


 !        call shape3H(NODES(Mdle)%type,Xi,Norder,Nedge_orient,Nface_orient, &
 !                     nrdofH,shapH,gradH)

 !        write(*,*) 'error_contr:shapH, nrdofH = ', nrdofH
 !        write(*,1000) shapH(1:nrdofH)
 ! 1000 format(40(4(e10.2,4x),/))     
 !        zsolH = zero
 !        do k=1,nrdofH
 !           ZsolH(1) = ZsolH(1) + ZdofH(1,k)*shapH(k)
 !           write(*,*) 'error_contr:k, ZsolH(1) = ', k, zsolH(1)
 !        enddo

 !        call shape3H(NODES(Mdle)%type,Xi,NorderC,Nedge_orientC,Nface_orientC, &
 !                     nrdofH,shapH,gradH)

 !        write(*,*) 'error_contr:shapH_C, nrdofH = ', nrdofH
 !        write(*,1000) shapH(1:nrdofH)

 !        zsolH_C = zero
 !        do k=1,nrdofH
 !           ZsolH_C(1) = ZsolH_C(1) + ZdofH_C(1,k)*shapH(k)

 !           write(*,*) 'error_contr:k, ZsolH_C(1) =', k, zsolH_C(1)
 !        enddo

 !        write(*,*) 'error_contr:zsolH = ', zsolH(1)
 !        write(*,*) 'error_contr:zsolH_C = ', zsolH_C(1)
 !        call pause

         call soleval(Mdle,xi,nedge_orient,nface_orient,norder,xnod,  &
                      zdofH,zdofE,zdofV,zdofQ,iflag,x,dxdxi, &
                      zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ) 

!     ...geometry
         call bgeom3D(mdle,xi,xnod,shapH,gradH,nrdofH,    &
                      dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
! 
!     ...weight        
         weight = bjac*wtloc(l)
!
!     ...evaluate the coarse approximate solution at the point 
         call soleval(MdleC,xiC,Nedge_orientC,Nface_orientC,norderC,XnodC,  &
                      zdofH_C,zdofE_C,zdofV_C,zdofQ_C,iflag,x,dxdxi,        &
                      zsolH_C,zdsolH_C,zsolE_C,zcurlE_C,zsolV_C,zdivV_C,zsolQ_C) 
!
!     ...accumulate for the errors
! 
!     ...loop through physical attributes of the node
         do iattr=1,NR_PHYSA
!
!        ...if the error not needed, skip
            if (Nflag(iattr) .eq. 0) cycle
!         
!        ...if attribute is absent, skip
            if (ncase(iattr) .eq. 0) cycle
!
!        ...loop over number of rhs
            do iload=1,NRCOMS
!               
!           ...loop through components of physical attribute
               do icomp=1,NR_COMP(iattr)
!                  
                  select case(DTYPE(iattr))
!
!              ...H1                   
                  case('contin')
                     ivarH = (iload-1)*NRHVAR+ADRES(iattr)+icomp
!
                     ErrorH(ivarH) = ErrorH(ivarH)          &
                                   + abs(zsolH(ivarH)-zsolH_C(ivarH))**2*weight
!
                  case('tangen')
                     ivarE = (iload-1)*NREVAR+ADRES(iattr)+icomp

                     call zcross_product(rn,zsolE-zsolE_C,zsolEn)

                     ErrorE(ivarE) = ErrorE(ivarE)             & 
                                   + (abs(zsolEn(1))**2        &
                                   +  abs(zsolEn(2))**2        &
                                   +  abs(zsolEn(3))**2)*weight
!
                  case('normal')
                     ivarV = (iload-1)*NRVVAR+ADRES(iattr)+icomp
!

                     call zdot_product(rn,zsolV-zsolV_C,zsolVn)
!                  
                     ErrorV(ivarV) = ErrorV(ivarV)          &
                                   + abs(zsolVn)**2*weight
                  end select                    
!
!           ...end of loop though components
               enddo
!  
!        ...end of loop though rhs's
            enddo
!
!     ...end of loop though physical attributes
         enddo
!
!  ...end of loop through integration points
      enddo

      ! write(*,*) 'if, ErrorH = ', if, ErrorH 
!
!..end of loop through faces
   enddo
!
!
   end subroutine error_contr
!
!
!--------------------------------------------------------------------
!                                                                     
!     routine name      - zcross_product
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Feb 2018
!                                                                     
!     purpose:          - compute cross product of real and
!                         (possibly) complex valued vectors
!                                                                    
!---------------------------------------------------------------------
!
   subroutine zcross_product(Rn, Za, Zcross)
!
   IMPLICIT NONE

#if C_MODE
   complex*16 :: Za(3), Zcross(3)
#else
   real*8     :: Za(3), Zcross(3)   
#endif
   real*8     :: Rn(3)
!                                                                    
!---------------------------------------------------------------------
!
   Zcross(1) =   Rn(2)*Za(3) - Rn(3)*Za(2)
   Zcross(2) = - Rn(1)*Za(3) + Rn(3)*Za(1)
   Zcross(3) =   Rn(1)*Za(2) - Rn(2)*Za(1)
!
!
   end subroutine zcross_product
!                                                                    
!                                                                    
!--------------------------------------------------------------------
!                                                                     
!     routine name      - zdot_product
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Feb 2018
!                                                                     
!     purpose:          - compute dot product of real and
!                         (possibly) complex valued vector
!                                                                    
!---------------------------------------------------------------------
!
   subroutine zdot_product(Rn, Za, Zprod)
!
   IMPLICIT NONE

#if C_MODE
   complex*16 :: Za(3), Zprod
#else
   real*8     :: Za(3), Zprod
#endif
   real*8     :: Rn(3)
!
!---------------------------------------------------------------------
!
   Zprod = Rn(1)*Za(1)+Rn(2)*Za(2)+Rn(3)*Za(3)
!
!
   end subroutine zdot_product