!
!
!--------------------------------------------------------------------
!
!     routine name      - get_sol_dofMG
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
   subroutine get_sol_dofMG(iel, zdofH, zdofE, zdofV)
!      
   use mg_data_structure
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
   call find_orderC(mdle, norder)

!
!..determine number of local dof
   call celndof(NODES(Mdle)%type,norder, nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
!
!..determine constraints' coefficients
   call logicC(Mdle,2, nodm,ndofmH,ndofmE,ndofmV,nrnodm,nrconH,nacH,constrH,   &
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
   end subroutine get_sol_dofMG







!-----------------------------------------------------------------------
!
!    routine name       - store_elem_solMG
!
!-----------------------------------------------------------------------
!
!    latest revision    - Mar 2018
!
!    purpose            - store element modified solution vector, 
!                         order and geometry dof
!
!----------------------------------------------------------------------
!
   subroutine store_elem_solMG(Iel, Zelem, Ndof)
! 
   use mg_data_structure
   use data_structure3D
   use macro_grid_info
!
   IMPLICIT NONE
!
   integer, intent(in) :: Iel, Ndof
#if C_MODE
   complex*16, intent(in) :: Zelem(Ndof)
#else
   real*8,     intent(in) :: Zelem(Ndof)
#endif   

   ZSOL_C(iel)%ndof_coarse = Ndof
   allocate(ZSOL_C(iel)%coarse(Ndof))
   ZSOL_C(iel)%coarse = Zelem

   call find_orderC(MDLE_MACRO(iel), ZSOL_C(iel)%norder)
   call find_orient(MDLE_MACRO(iel),ZSOL_C(iel)%nedge_orient,ZSOL_C(iel)%nface_orient)
   call nodcor(MDLE_MACRO(iel),ZSOL_C(iel)%xnod)

!..determine element dof
   allocate(ZSOL_C(iel)%zdofH(MAXEQNH,MAXbrickH))
   allocate(ZSOL_C(iel)%zdofE(MAXEQNE,MAXbrickE))
   allocate(ZSOL_C(iel)%zdofV(MAXEQNV,MAXbrickV))
   allocate(ZSOL_C(iel)%zdofQ(MAXEQNQ,MAXbrickQ))
!
   call get_sol_dofMG(Iel, ZSOL_C(iel)%zdofH, ZSOL_C(iel)%zdofE, ZSOL_C(iel)%zdofV) 
!
!
   end subroutine store_elem_solMG