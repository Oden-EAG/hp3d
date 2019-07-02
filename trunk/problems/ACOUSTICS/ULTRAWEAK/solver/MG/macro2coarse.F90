!-----------------------------------------------------------------------
!
!   routine name       - macro2coarse
!
!-----------------------------------------------------------------------
!
!   latest revision    - Mar 2018
!
!   purpose            - routine applies the transpose of the prolongation 
!                        operator to macro-grid residual vector using the 
!                        coefficients from logic_macro
!                        (including Dirichlet dof)
!
!   arguments :
!     in:
!             Iel      - natural order of elements
!             Zu_m     - residual vector on the modified macro-element
!             nrdof_m  - size of Zu_m
!     out:
!             Zu_c     - residual vector on the modified coarse element
!             nrdof_c  - size of Zu_c
!
!----------------------------------------------------------------------
!
   subroutine macro2coarse(Iel,Zu_m,Nrdof_m,Zu_c,Nrdof_c)
!
   use data_structure3D, only: NODES, get_index
   use physics,          only: NRHVAR, NREVAR, NRVVAR, NRINDEX     
   use parameters,       only: NRCOMS, MAXEQNH, MAXEQNE, MAXEQNV,      &
                               MAXbrickH, MAXbrickE, MAXbrickV, ZERO                   
   use constraints_info, only: CONSTR
   use assembly,         only: NR_RHS 
!
   IMPLICIT NONE
!
!
   integer,    intent(in)  :: Iel, Nrdof_m, Nrdof_c
#if C_MODE
   complex*16, intent(in)  :: Zu_m(nrdof_m) 
   complex*16, intent(out) :: Zu_c(nrdof_c)
#else
   real*8,     intent(in)  :: Zu_m(nrdof_m) 
   real*8,     intent(out) :: Zu_c(nrdof_c)
#endif
!
!..locals
#if C_MODE
   complex*16           :: zdofH_c(MAXEQNH,MAXbrickH),     &
                           zdofE_c(MAXEQNE,MAXbrickE),     &
                           zdofV_c(MAXEQNV,MAXbrickV) 
   complex*16           :: zdofH_m(MAXEQNH,CONSTR(Iel)%nrdofH_macro),   &
                           zdofE_m(MAXEQNE,CONSTR(Iel)%nrdofE_macro),   &
                           zdofV_m(MAXEQNV,CONSTR(Iel)%nrdofV_macro) 
#else
   real*8               :: zdofH_c(MAXEQNH,MAXbrickH),     &
                           zdofE_c(MAXEQNE,MAXbrickE),     &
                           zdofV_c(MAXEQNV,MAXbrickV)    
   real*8               :: zdofH_m(MAXEQNH,CONSTR(Iel)%nrdofH_macro),   &
                           zdofE_m(MAXEQNE,CONSTR(Iel)%nrdofE_macro),   &
                           zdofV_m(MAXEQNV,CONSTR(Iel)%nrdofV_macro)                                  
#endif
!
!..component counters for the nodes (used in case of multiple loads)
   integer, allocatable :: mvarH(:),mvarE(:), mvarV(:)
!
   integer              :: iH(MAXEQNH), iE(MAXEQNE), iV(MAXEQNV)  
   integer              :: nrnod_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
   integer              :: nrnodm, nvarH, nvarE, nvarV, ivar
   integer              :: naH, naE, naV,i,k,load,nn,nod,j
   integer              :: kH, kE, kV, kp,lH, lE,lV
!   
!..decoded index for a node
   integer              :: index(NRINDEX)
!   
!--------------------------------------------------------------------------
!
!..macro-element number of nodes
   nrnod_macro = CONSTR(Iel)%nrnod_macro
!   
!..modified coarse element number of nodes
   nrnodm   = CONSTR(Iel)%nrnodm
!
!--------------------------------------------------------------------------
!
!..compute total number of dof for the macro-element
   naH = 0; naE = 0; naV = 0
   do i=1,nrnod_macro
      naH = naH + CONSTR(Iel)%ndofH_macro(i)
      naE = naE + CONSTR(Iel)%ndofE_macro(i)
      naV = naV + CONSTR(Iel)%ndofV_macro(i)
   enddo

   nrdofH_macro = naH; nrdofE_macro = naE; nrdofV_macro = naV
!
!..initiate dof counter
   nn = 0; 
!..initialize the number of H1,H(curl),H(div)
   ivar=0

!..initialize the macro dof matrices
   ZdofH_m = ZERO; ZdofE_m = ZERO; ZdofV_m = ZERO
!
!..initiate component counters
   allocate(mvarH(nrnod_macro)); mvarH = 0
   allocate(mvarE(nrnod_macro)); mvarE = 0
   allocate(mvarV(nrnod_macro)); mvarV = 0
!
!..loop through loads
   do load=1,NR_RHS
!      
!  ...H1 dof 
      if (NRHVAR.eq.0) goto 400
!
!  ...initialize H1 column indices
      iH = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnod_macro
!         
!     ...pick up the node number from the list   
         nod = CONSTR(Iel)%nod_macro(i)
!         
!     ...compute the number of active H1 variables for the node         
         call get_index(nod, index)
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,CONSTR(Iel)%ndofH_macro(i)
!        ...loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
               select case(index(k))
               case(1)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  ! zdofH_m(ivar,iH(ivar)) = NODES(nod)%zdofH(ivar,j)
               case(2)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  nn = nn+1
                  zdofH_m(ivar,iH(ivar)) = Zu_m(nn) 
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
  400 continue
!      
!  ...H(curl) dof 
      if (NREVAR.eq.0) goto 500
!      
!  ...initialize H(curl) column indices
      iE = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnod_macro
!     ...pick up the node number from the list   
         nod = CONSTR(Iel)%nod_macro(i)
!         
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,CONSTR(Iel)%ndofE_macro(i)!/nvarE
!
!        ...loop through the components
            ivar=mvarE(i)
            do k=1,NRINDEX
               select case(index(k))
               case(3)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  ! zdofE_m(ivar,iE(ivar)) = NODES(nod)%zdofE(ivar,j) 
               case(4)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  nn=nn+1
                  zdofE_m(ivar,iE(ivar)) = Zu_m(nn)
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

  500 continue
!
!  ...H(div) dof 
      if (NRVVAR.eq.0) goto 600
!      
!  ...initialize H(div) column indices
      iV = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnod_macro

!     ...pick up the node number from the list   
         nod = CONSTR(Iel)%nod_macro(i)
!         
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,CONSTR(Iel)%ndofV_macro(i)!/nvarV
!
!        ...loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
               select case(index(k))
               case(5)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  ! zdofV_m(ivar,iV(ivar)) = NODES(nod)%zdofV(ivar,j) 
               case(6)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  nn=nn+1
                  zdofV_m(ivar,iV(ivar)) = Zu_m(nn)  
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
!
  600 continue
!
!..end of loop through loads
   enddo  
!
   deallocate(mvarH, mvarE, mvarV)
!
!...calculate the coarse element dof
!
!..initialize the coarse dof matrices
   ZdofH_c = ZERO; ZdofE_c = ZERO; ZdofV_c = ZERO
!      
!..loop through the macro-element H1 dof
   do kH=1,nrdofH_macro
!
!  ...accumulate for the values
      do kp=1,CONSTR(Iel)%nrconH_macro(kH)
         lH = CONSTR(Iel)%nacH_macro(kp,kH)
         do ivar=1,NRHVAR*NRCOMS
            ZdofH_c(ivar,lH) = ZdofH_c(ivar,lH)   &
                             + CONSTR(Iel)%constrH_macro(kp,kH)*zdofH_m(ivar,kH)
         enddo
      enddo
!...end of loop through macro-element H1 dof
   enddo
!      
!..loop through the macro-element H(curl) dof
   do kE=1,nrdofE_macro
!
!  ...accumulate for the values
      do kp=1,CONSTR(Iel)%nrconE_macro(kE)
         lE = CONSTR(Iel)%nacE_macro(kp,kE)
         do ivar=1,NREVAR*NRCOMS
            zdofE_c(ivar,lE) = zdofE_c(ivar,lE)    &
                             + CONSTR(Iel)%ConstrE_macro(kp,kE)*ZdofE_m(ivar,kE)
         enddo
      enddo
!...end of loop through macro-element H(curl) dof
   enddo
!
!..loop through the macro-element H(curl) dof
   do kV=1,nrdofV_macro
!
!  ...accumulate for the values
      do kp=1,CONSTR(Iel)%nrconV_macro(kV)
         lV = CONSTR(Iel)%nacV_macro(kp,kV)
         do ivar=1,NRVVAR*NRCOMS
            zdofV_c(ivar,lV) = zdofV_c(ivar,lV) &
                             + CONSTR(Iel)%constrV_macro(kp,kV)*ZdofV_m(ivar,kV)
         enddo
      enddo
!...end of loop through macro-element H(curl) dof
   enddo
!
!..reconstruct the coarse element residual vector
!
!..initiate dof counter
   nn = 0; 
!..initialize the number of H1,H(curl),H(div),L2 copied so far
   ivar=0

!..initialize component counters
   allocate(mvarH(nrnodm)); mvarH = 0
   allocate(mvarE(nrnodm)); mvarE = 0
   allocate(mvarV(nrnodm)); mvarV = 0
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
!  ...loop through coarse element nodes
      do i = 1, nrnodm
!         
!     ...pick up the node number from the list   
         nod = CONSTR(Iel)%nodm(i)
!         
!     ...compute the number of active H1 variables for the node         
         call get_index(nod, index)
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,CONSTR(Iel)%ndofmH(i)
!        ...loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
               select case(index(k))
               case(1)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  ! NODES(nod)%zdofH(ivar,j) = zdofH_c(ivar,iH(ivar)) 
               case(2)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  nn = nn+1
                  Zu_c(nn) = zdofH_c(ivar,iH(ivar)) 
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
         nod = CONSTR(Iel)%nodm(i)
!         
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         if (nvarE.eq.0) cycle
!         
         do j=1,CONSTR(Iel)%ndofmE(i)
!
!        ...loop through the components
            ivar=mvarE(i)
            do k=1,NRINDEX
               select case(index(k))
               case(3)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  ! NODES(nod)%zdofE(ivar,j) = zdofE_c(ivar,iH(ivar))                   
               case(4)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  nn=nn+1
                  Zu_c(nn) = zdofE_c(ivar,iE(ivar))
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
         nod = CONSTR(Iel)%nodm(i)
!         
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,CONSTR(Iel)%ndofmV(i)!/nvarV
!
!        ...loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
               select case(index(k))
               case(5)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  ! NODES(nod)%zdofV(ivar,j) = zdofV_c(ivar,iH(ivar))                   
               case(6)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  nn=nn+1
                  Zu_c(nn) = zdofV_c(ivar,iV(ivar)) 
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
!
   end subroutine macro2coarse