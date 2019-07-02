!-----------------------------------------------------------------------
!
!   routine name       - coarse2macro
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2018
!
!   purpose            - routine calculates solution vector 
!                        for a macro element (including Dirichlet dof)
!
!   arguments :
!     in:
!             Iel      - natural order of elements
!             Zu_c     - solution vector for the modified coarse element
!             nrdof_c  - length of Zu_c
!     out:
!             Zu_m     - solution vector for the modified macro-element
!             nrdof_m  - length of Zu_m
!
!----------------------------------------------------------------------

   subroutine coarse2macro_test(Iel,Zu_c,Nrdof_c,Zu_m,Nrdof_m)
!
   use data_structure3D
   use constraints_info 
!
   IMPLICIT NONE
!
!
   integer,    intent(in)  :: Iel,  Nrdof_c, Nrdof_m
#if C_MODE
   complex*16, intent(in)  :: Zu_c(nrdof_c) 
   complex*16, intent(out) :: Zu_m(nrdof_m)
#else
   real*8,     intent(in)  :: Zu_c(nrdof_c) 
   real*8,     intent(out) :: Zu_m(nrdof_m)
#endif
!
!..locals
#if C_MODE
   complex*16              :: zdofH_c(MAXEQNH,MAXbrickH),     &
                              zdofE_c(MAXEQNE,MAXbrickE),     &
                              zdofV_c(MAXEQNV,MAXbrickV) 
   complex*16              :: zdofH_m(MAXEQNH,MAX_MACRO_H),   &
                              zdofE_m(MAXEQNE,MAX_MACRO_E),   &
                              zdofV_m(MAXEQNV,MAX_MACRO_V) 
#else
   real*8                  :: zdofH_c(MAXEQNH,MAXbrickH),     &
                              zdofE_c(MAXEQNE,MAXbrickE),     &
                              zdofV_c(MAXEQNV,MAXbrickV)    
   real*8                  :: zdofH_m(MAXEQNH,MAX_MACRO_H),   &
                              zdofE_m(MAXEQNE,MAX_MACRO_E),   &
                              zdofV_m(MAXEQNV,MAX_MACRO_V)                                  
#endif


!..work space for constraints
   integer :: nr_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
   integer :: nrnodm, nvarH, nvarE, nvarV, ivar
   integer :: naH, naE, naV,i,k,load,nn,nod,j
!..decoded index for a node
   integer :: index(NRINDEX)
!   
!..component counters for the nodes (use in case of multiple loads)
   integer, allocatable :: mvarH(:),mvarE(:), mvarV(:)
   integer              :: kH, kE, kV, kp,lH, lE,lV, iH, iE, iV
!               
!
!--------------------------------------------------------------------------
!
!..initialize 
   zdofH_c = ZERO; zdofE_c = ZERO; zdofV_c = ZERO 
!

!..macro-element
   nr_macro = CONSTR(Iel)%nr_macro
!..modified coarse element
   nrnodm   = CONSTR(Iel)%nrnodm
!
!--------------------------------------------------------------------------
!
!..compute total number of dof for the macro-element
   naH = 0; naV = 0
   do i=1,nr_macro
      naH = naH + CONSTR(Iel)%ndofH_macro(i)
      naV = naV + CONSTR(Iel)%ndofV_macro(i)
   enddo
!
   nrdofH_macro = naH; nrdofV_macro = naV
!
!..initiate dof counter
   nn = 0 
!   
!..initiate H1 dof counter
   iH = 0
!..loop through coarse element nodes
   do i = 1, nrnodm
!         
!  ...pick up the node number from the list   
      nod = CONSTR(Iel)%nodm(i)
!         
!  ...compute the number of active H1 variables for the node         
      call get_index(nod, index)

!  ...loop through the nodal dof (potentially NONE)
      do k=1,NRINDEX
         select case(index(k))
         case(1)
            do j=1,CONSTR(Iel)%ndofmH(i)
               iH = iH + 1
               zdofH_c(1,iH) = NODES(nod)%zdofH(1,j)
            enddo   
         case(2)
            do j=1,CONSTR(Iel)%ndofmH(i)
               iH = iH + 1
               nn = nn+1
               zdofH_c(1,iH) = Zu_c(nn)
            enddo   
         end select
!
!  ...end of loop through NRINDEX               
      enddo   
!
!..end of loop through coarse element nodes            
   enddo   
!
!..initiate H(div) dof counter
   iV = 0   
!   
!..loop through coarse element nodes
   do i = 1, nrnodm
!      
!  ...pick up the node number from the list   
      nod = CONSTR(Iel)%nodm(i)
!         
!  ...compute the number of active H(curl) variables for the node
      call get_index(nod, index)
!
      do k=1,NRINDEX
         select case(index(k))
         case(5)
            do j=1,CONSTR(Iel)%ndofmV(i)
               iV = iV + 1
               zdofV_c(1,iV) = NODES(nod)%zdofV(1,j)
            enddo   
         case(6)
            do j=1,CONSTR(Iel)%ndofmV(i)
               iV = iV + 1
               nn=nn+1
               zdofV_c(1,iV) = Zu_c(nn)
            enddo   
         end select
!
!  ...end of loop through NRINDEX               
      enddo
!            
!..end of loop through the nodal dof
   enddo
!
!
!--------------------------------------------------------------------
!
!...calculate the macro element dof
!
!..initiate the dof
   ZdofH_m = ZERO; ZdofE_m = ZERO; ZdofV_m = ZERO
!      
!..loop through the macro-element H1 dof
   do kH=1,nrdofH_macro
!
!  ...accumulate for the values
      do kp=1,CONSTR(Iel)%nrconH_macro(kH)
         lH = CONSTR(Iel)%nacH_macro(kp,kH)
         ZdofH_m(1,kH) = ZdofH_m(1,kH)     &
                       + CONSTR(Iel)%constrH_macro(kp,kH)*zdofH_c(1,lH)
      enddo
!...end of loop through macro-element H1 dof
   enddo
!      
!..loop through the macro-element H(curl) dof
   do kV=1,nrdofV_macro
!
!  ...accumulate for the values
      do kp=1,CONSTR(Iel)%nrconV_macro(kV)
         lV = CONSTR(Iel)%nacV_macro(kp,kV)
         ZdofV_m(1,kV) = ZdofV_m(1,kV)      &
                       + CONSTR(Iel)%constrV_macro(kp,kV)*zdofV_c(1,lV)
      enddo
!...end of loop through macro-element H(curl) dof
   enddo

!   
!--------------------------------------------------------------------
!
!..reconstruct the macro element solution vector
!
!..initiate dof counter
   nn = 0; 
!
!..initiate H1 dof counter
   iH = 0
!..loop through coarse element nodes
   do i = 1, nr_macro
!         
!  ...pick up the node number from the list   
      nod = CONSTR(Iel)%nod_macro(i)
!         
!  ...compute the number of active H1 variables for the node         
      call get_index(nod, index)

      do k=1,NRINDEX
         select case(index(k))
         case(1)
            do j=1,CONSTR(Iel)%ndofH_macro(i)
               iH = iH + 1
            enddo   
         case(2)
            do j=1,CONSTR(Iel)%ndofH_macro(i)
               iH = iH + 1
               nn = nn+1
               Zu_m(nn) = zdofH_m(1,iH)
            enddo
         end select
!
!  ...end of loop through NRINDEX               
      enddo   
!..end of loop through the nodal dof
   enddo
!
!      
!..initiate H1 dof counter
   iV = 0
!..loop through coarse element nodes
   do i = 1, nr_macro
!         
!  ...pick up the node number from the list   
      nod = CONSTR(Iel)%nod_macro(i)
!         
!  ...compute the number of active H1 variables for the node         
      call get_index(nod, index)

      do k=1,NRINDEX
         select case(index(k))
         case(5)
            do j=1,CONSTR(Iel)%ndofV_macro(i)
               iV = iV + 1
            enddo   
         case(6)
            do j=1,CONSTR(Iel)%ndofV_macro(i)
               iV = iV + 1
               nn = nn+1
               Zu_m(nn) = zdofV_m(1,iV)
            enddo
         end select
!
!  ...end of loop through NRINDEX               
      enddo   
!..end of loop through the nodal dof
   enddo
!
   end subroutine coarse2macro_test