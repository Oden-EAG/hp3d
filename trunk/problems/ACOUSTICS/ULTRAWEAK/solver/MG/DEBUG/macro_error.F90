!--------------------------------------------------------------------
!
!     routine name      - prolong_error
!
!--------------------------------------------------------------------
!
!     latest revision:  - Feb 2018
!
!     purpose:          - routine evaluates L2 norm of the difference
!                         between the coarse and fine grid solutions
!                         on the boundary of a coarse grid element
!
!     arguments:
!
!---------------------------------------------------------------------
!
   subroutine prolong_error
!
   use macro_grid_info 
!
   IMPLICIT NONE
!
   integer :: iel, mdle
   real*8  :: errH, errE, errV, errorH, errorE, errorV
!
   errorH = 0.0d0; errorE = 0.0d0; errorV = 0.d0
!
   do iel = 1, NRELES_COARSE
      mdle = MDLE_MACRO(iel)
      call macro_error(iel, mdle, errH, errE, errV)

      ! write(*,*) 'macro_error: mdle = ', mdle
      ! write(*,*) 'macro_error:mdle, errH = ',mdle, errH
      ! write(*,*) 'macro_error: errV = ', errV


      errorH = errorH + errH 
      errorE = errorE + errE 
      errorV = errorV + errV
   enddo    
!
   write(*,*) ''
   write(*,1001) ErrorH
   write(*,1002) ErrorE
   write(*,1003) ErrorV
   write(*,*) ''
 1001 format(' macro_error: ErrorH = ',12x, es13.4)    
 1002 format(' macro_error: ErrorE = ',12x, es13.4)    
 1003 format(' macro_error: ErrorV = ',12x, es13.4)    
!
!
!
   end subroutine prolong_error
!
!--------------------------------------------------------------------
!
!     routine name      - macro_error
!
!--------------------------------------------------------------------
!
!     latest revision:  - Feb 2018
!
!     purpose:          - routine evaluates L2 norm of the difference
!                         between the coarse and fine grid solutions
!                         on the boundary of a coarse grid element
!
!     arguments:
!
!     in:
!              MdleC    - coarse grid element
!
!---------------------------------------------------------------------
!
   subroutine macro_error(Iel,MdleC,ErrH,ErrE,ErrV)
!
   use data_structure3D
   use macro_grid_info, ONLY: ZSOL_C
!   
   IMPLICIT NONE
!   
!
   integer, intent(in)  :: Iel, MdleC
   real*8,  intent(out) :: errH, errE, errV
!
!..locals
!..flag indicating for which variables the error should be computed
   integer    :: nflag(NR_PHYSA), mdle
!
!..dof for the coarse grid
#if C_MODE   
   complex*16 :: zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE)
   complex*16 :: zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
#else
   real*8     :: zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE)
   real*8     :: zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
#endif
!
!..workspace for error_contr
   integer    :: nedge_orient(12), nface_orient(6) 
   integer    :: norder(19)
   real*8     :: xnod(NDIMEN,MAXbrickH)
!..the total error for the element
   real*8     :: errorH( MAXEQNH), errorE( MAXEQNE), errorV( MAXEQNV)
   real*8     :: derrorH(MAXEQNH), derrorE(MAXEQNE), derrorV(MAXEQNV)   
   integer    :: iprint
!
!---------------------------------------------------------------------
!
   select case(MdleC)
   case(236)
      iprint = 1
   case default
      iprint = 0
   end select    
!     
!..compute the error on the traces and fluxes 
   nflag=0; nflag(1:2)=1
!
!..coarse grid order 
   norder = ZSOL_C(Iel)%norder
!..coarse grid nod coordinates
   xnod = ZSOL_C(Iel)%xnod
!   
   nedge_orient = ZSOL_C(Iel)%nedge_orient
   nface_orient = ZSOL_C(Iel)%nface_orient
!..determine element dof
   zdofH = ZSOL_C(Iel)%zdofH
   zdofE = ZSOL_C(Iel)%zdofE
   zdofV = ZSOL_C(Iel)%zdofV
   zdofQ = ZERO
!
!..initiate errors and norms
   errorH = 0.d0
   errorE = 0.d0
   errorV = 0.d0
!
!..loop through fine grid sub-elements
   mdle=0
   do
      call nelcon_macro(MdleC,mdle, mdle)
!
      if (mdle.eq.0) exit
!      
      call error_contr(MdleC,mdle,nflag,norder,nedge_orient, nface_orient, xnod,   &
                       zdofH,zdofE,zdofV,zdofQ, derrorH, derrorE, derrorV)
!  ...accumulate
      errorH = errorH + derrorH 
      errorE = errorE + derrorE 
      errorV = errorV + derrorV
   enddo
!
   ErrH = errorH(1)
   ErrE = errorE(1)
   ErrV = errorV(1)
!
!
   end subroutine macro_error
!
!
!--------------------------------------------------------------------
!
!     routine name      - get_sol_dof
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2018
!
!     purpose:          - returns solution degrees of freedom for an 
!                         element including Dirichlet dof
!
!     arguments:
!
!     in:
!              iel      - natural oder of a coarse element
!     out: 
!            zdofH      - H1      solution degrees of freedom
!            zdofE      - H(curl) solution degrees of freedom
!            zdofV      - H(div)  solution degrees of freedom
!
!---------------------------------------------------------------------
!  
   subroutine get_sol_dof(iel, zdofH, zdofE, zdofV)
!      
   use data_structure3D, only: NODES, get_index, ndof_nod, NACDIM
   use physics,          only: NRHVAR, NREVAR, NRVVAR, NRINDEX, NR_PHYSA     
   use parameters,       only: NRCOMS, MAXEQNH, MAXEQNE, MAXEQNV, MAXEQNQ,   &
                               MAXbrickH, MAXbrickE, MAXbrickV, MAXbrickQ, ZERO            
   use assembly,         only: NR_RHS, MAXDOFM, MAXbrickH,MAXbrickE,       &
                               MAXbrickV, MAXbrickQ, NRHVAR, NREVAR,       &
                               NRVVAR, NRQVAR, MAXDOFS, MAXDOFC, NEXTRACT, &
                               IDBC, ZDOFD, ZERO, BLOC, AAUX,ALOC, ZBMOD,  &
                               ZAMOD, NR_PHYSA, MAXNODM
   use macro_grid_info,  only: ZSOL_C, NRELES_COARSE, MDLE_MACRO
!   
   IMPLICIT NONE
!
   integer,    intent(in)  :: Iel

#if C_MODE
   complex*16, intent(out) :: zdofH(MAXEQNH,MAXbrickH),     &
                              zdofE(MAXEQNE,MAXbrickE),     &
                              zdofV(MAXEQNV,MAXbrickV)
#else
   real*8,     intent(out) :: zdofH(MAXEQNH,MAXbrickH),     &
                              zdofE(MAXEQNE,MAXbrickE),     &
                              zdofV(MAXEQNV,MAXbrickV)
#endif
!
#if C_MODE
   complex*16              :: zvoid(1)
#else
   real*8                  :: zvoid(1)     
#endif

!
#if C_MODE
   complex*16 :: zdofHm(MAXEQNH,2*MAXbrickH),     &
                 zdofEm(MAXEQNE,2*MAXbrickE),     &
                 zdofVm(MAXEQNV,2*MAXbrickV)
#else
   real*8     :: zdofHm(MAXEQNH,2*MAXbrickH),     &
                 zdofEm(MAXEQNE,2*MAXbrickE),     &
                 zdofVm(MAXEQNV,2*MAXbrickV)
#endif                 
!..element order of approximation
   integer    :: norder(19)
!
!..modified element nodes and corresponding number of dof
   integer    :: nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),   &
                 ndofmV(MAXNODM),ndofmQ(MAXNODM)
!
   integer    :: nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV),   &
                 nacH(NACDIM,MAXbrickH),nacE(NACDIM,MAXbrickE),           &
                 nacV(NACDIM,MAXbrickV)
   real*8     :: constrH(NACDIM,MAXbrickH),constrE(NACDIM,MAXbrickE),  &   
                 constrV(NACDIM,MAXbrickV)               
!
!..component counters for the nodes (used in case of multiple loads)
   integer, allocatable    :: mvarH(:),mvarE(:), mvarV(:)
!
   integer                 :: iH(MAXEQNH), iE(MAXEQNE), iV(MAXEQNV)  
   integer                 :: nrnodm, nvarH, nvarE, nvarV, ivar
   integer                 :: naH, naE, naV,i,k,load,nn,nod,j
   integer                 :: kH, kE, kV, kp,lH, lE,lV
   integer                 :: nrdofm,nrdofc, mdle, nvoid
   integer                 :: nrdoflH, nrdoflE, nrdoflV, nrdoflQ
   integer                 :: l, kQ
!   
!..decoded index for a node
   integer                 :: index(NRINDEX)
!
!
!-----------------------------------------------------------------------------
!
   mdle = MDLE_MACRO(iel)
   call find_order(mdle, norder)

!
!..determine number of local dof
   call celndof(NODES(Mdle)%type,norder, nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
!
!..determine constraints' coefficients
   call logic(Mdle,2, nodm,ndofmH,ndofmE,ndofmV,nrnodm,nrconH,nacH,constrH,   &
              nrconE,nacE,constrE,nrconV,nacV,constrV)
!
!---------------------------------------------------------------------
!
!
!
!..remove the mdle node
   nrnodm = nrnodm - 1

!..initiate dof counter
   nn = 0; 
!..initialize the number of H1,H(curl),H(div),L2
   ivar=0

!..initialize component counters
   allocate(mvarH(nrnodm)); mvarH = 0
   allocate(mvarE(nrnodm)); mvarE = 0
   allocate(mvarV(nrnodm)); mvarV = 0

!..initialize dof matrices
   zdofHm = ZERO; zdofEm = ZERO; zdofVm = ZERO 
!
!..loop through loads
   do load=1,NR_RHS
!      
!  ...H1 dof 
      if (NRHVAR.eq.0) goto 100
!      
!  ...initialize H1 column indices
      iH = 0
!
!  ...loop through coarse element nodes (expept the mdle node)
      do i = 1, nrnodm
!         
!     ...pick up the node number from the list   
         nod = nodm(i)
!         
!     ...compute the number of active H1 variables for the node         
         call get_index(nod, index)
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,ndofmH(i)
!        ...loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
               select case(index(k))
               case(1)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  zdofHm(ivar,iH(ivar)) = NODES(nod)%zdofH(ivar,j)
               case(2)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  nn = nn+1
                  zdofHm(ivar,iH(ivar)) = ZSOL_C(iel)%coarse(nn)
               end select
!
!        ...end of loop through NRINDEX               
            enddo   
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarH(i) = ivar   

!  ...end of loop through coarse element nodes            
      enddo   
!
  100 continue
!      
!  ...H(curl) dof 
      if (NREVAR.eq.0) goto 200
!      
!  ...initialize H(curl) column indices
      iE = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnodm
!     ...pick up the node number from the list   
         nod = nodm(i)
!         
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         if (nvarE.eq.0) cycle
!         
         do j=1,ndofmE(i)
!
!        ...loop through the components
            ivar=mvarE(i)
            do k=1,NRINDEX
               select case(index(k))
               case(3)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  zdofEm(ivar,iE(ivar)) = NODES(nod)%zdofE(ivar,j)
               case(4)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  nn=nn+1
                  zdofEm(ivar,iE(ivar)) = ZSOL_C(iel)%coarse(nn)
               end select
!
!        ...end of loop through NRINDEX               
            enddo
!            
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarE(i) = ivar

!  ...end of loop through coarse element nodes            
      enddo  
!
  200 continue
!
!  ...H(div) dof 
      if (NRVVAR.eq.0) goto 300
!      
!  ...initialize H(div) column indices
      iV = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnodm
!     ...pick up the node number from the list   
         nod = nodm(i)
!         
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,ndofmV(i)!/nvarV
!
!        ...loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
               select case(index(k))
               case(5)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  zdofVm(ivar,iV(ivar)) = NODES(nod)%zdofV(ivar,j)
               case(6)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  nn=nn+1
                  zdofVm(ivar,iV(ivar)) = ZSOL_C(iel)%coarse(nn)
               end select
!
!        ...end of loop through NRINDEX               
            enddo
!            
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarV(i) = ivar

!  ...end of loop through coarse element nodes            
      enddo  

  300 continue

!..end of loop through loads
   enddo  
!
   deallocate(mvarH, mvarE, mvarV)
!
!..initiate the dof
   ZdofH = ZERO
   ZdofE = ZERO
   ZdofV = ZERO
!      
!
!..loop through the local dof
   do k=1,nrdoflH
!      
!  ...accumulate for the values
      do kp=1,nrconH(k)
         l = nacH(kp,k)
         do ivar=1,NRHVAR*NRCOMS
            ZdofH(ivar,k) = ZdofH(ivar,k) + constrH(kp,k)*zdofHm(ivar,l)
         enddo
      enddo
!..loop through local dof
   enddo
!
!..loop through the local dof
   do kE=1,nrdoflE
!
!  ...accumulate for the values
      do kp=1,nrconE(kE)
         l = nacE(kp,kE)
         do ivar=1,NREVAR*NRCOMS
            ZdofE(ivar,kE) = ZdofE(ivar,kE) + constrE(kp,kE)*zdofEm(ivar,l)
         enddo
      enddo
!
!..loop through local H(curl) dof
   enddo

!..loop through the local dof
   do kV=1,nrdoflV
!
!  ...accumulate for the values
      do kp=1,nrconV(kV)
         l = nacV(kp,kV)
         do ivar=1,NRVVAR*NRCOMS
            ZdofV(ivar,kV) = ZdofV(ivar,kV) + constrV(kp,kV)*zdofVm(ivar,l)
         enddo
      enddo
!
!..loop through local H(div) dof
   enddo
!
!
   end subroutine get_sol_dof