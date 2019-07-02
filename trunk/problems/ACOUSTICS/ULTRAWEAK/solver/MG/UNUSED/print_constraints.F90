


   subroutine print_constraints
   use data_structure3D, ONLY: NODES
!
   implicit none
!
!..workspace for coarse element
   integer              :: nrnodm
   integer, allocatable :: nodm(:), ndofmH(:), ndofmE(:), ndofmV(:)
!
!..workspace for macro element
   integer              :: nrnod_macro
   integer, allocatable :: nod_macro(:),    &
                           ndofH_macro(:), ndofE_macro(:), ndofV_macro(:)
!
!..workspace for constraints
   integer              :: nrconH_macro(MAX_MACRO_H), nacH_macro(NACDIM,MAX_MACRO_H) 
   integer              :: nrconE_macro(MAX_MACRO_E), nacE_macro(NACDIM,MAX_MACRO_E)
   integer              :: nrconV_macro(MAX_MACRO_V), nacV_macro(NACDIM,MAX_MACRO_V)
   real*8               :: constrH_macro(NACDIM,MAX_MACRO_H), &
                           constrE_macro(NACDIM,MAX_MACRO_E), &
                           constrV_macro(NACDIM,MAX_MACRO_V)
!..locals
   integer              :: iel, mdle,i
   integer              :: nrdofH_macro, nrdofE_macro, nrdofV_macro
   integer              :: kH, kE
!
!---------------------------------------------------------------------------------
!
!..loop through all coarse elements
   do iel = 1, NRELES_COARSE
      mdle = MDLE_MACRO(iel) 
!  ...copy to locals
!
!  ...modified coarse element
      nrnodm = CONSTR(iel)%nrnodm

      allocate(nodm(nrnodm));   nodm   = CONSTR(iel)%nodm
      allocate(ndofmH(nrnodm)); ndofmH = CONSTR(iel)%ndofmH
      allocate(ndofmE(nrnodm)); ndofmE = CONSTR(iel)%ndofmE
      allocate(ndofmV(nrnodm)); ndofmV = CONSTR(iel)%ndofmV

!
!  ...macro element 
      nrnod_macro  = CONSTR(iel)%nrnod_macro
      allocate(nod_macro(nrnod_macro));   nod_macro = CONSTR(iel)%nod_macro
      allocate(ndofH_macro(nrnod_macro)); ndofH_macro = CONSTR(iel)%ndofH_macro
      allocate(ndofE_macro(nrnod_macro)); ndofE_macro = CONSTR(iel)%ndofE_macro
      allocate(ndofV_macro(nrnod_macro)); ndofV_macro = CONSTR(iel)%ndofV_macro
!
!  ...constraints
      nrconH_macro = CONSTR(iel)%nrconH_macro
      nrconE_macro = CONSTR(iel)%nrconE_macro
      nrconV_macro = CONSTR(iel)%nrconV_macro
!
      nacH_macro = CONSTR(iel)%nacH_macro
      nacE_macro = CONSTR(iel)%nacE_macro
      nacV_macro = CONSTR(iel)%nacV_macro
!
      constrH_macro = CONSTR(iel)%constrH_macro
      constrE_macro = CONSTR(iel)%constrE_macro
      constrV_macro = CONSTR(iel)%constrV_macro
!
!  ...printing
      write(*,1000) mdle
 1000 format(' mdleC       =', i7) 
      write(*,1001) nrnodm
 1001 format(' nrnodm      =', i7) 
      write(*,1002) nodm(1:nrnodm)
 1002 format(' nodm        =',/, 10(9x,7i7,/)) 
      ! write(*,1003) ndofmH(1:nrnodm)
 ! 1003 format(' ndofmH      =',/, 10(9x,7i7,/)) 
      write(*,1004) ndofmE(1:nrnodm)
 1004 format(' ndofmE      =',/, 10(9x,7i7,/)) 
      ! write(*,1005) ndofmV(1:nrnodm)
 ! 1005 format(' ndofmV      =',/, 10(9x,7i7,/)) 
      write(*,1006) nrnod_macro
 1006 format(' nrnod_macro    =', i7)  
      write(*,1007) nod_macro(1:nrnod_macro)
 1007 format(' nod_macro   =',/, 50(9x,7i7,/))  
      ! write(*,1008) ndofH_macro(1:nrnod_macro)
 ! 1008 format(' ndofH_macro =',/, 50(9x,7i7,/))  
      write(*,1009) ndofE_macro(1:nrnod_macro)
 1009 format(' ndofE_macro =',/, 50(9x,7i7,/)) 
      ! write(*,1010) ndofV_macro(1:nrnod_macro)
 ! 1010 format(' ndofV_macro =',/, 50(9x,7i7,/)) 

!  ...compute total number of dof for the macro-element
      nrdofH_macro = sum(NdofH_macro(1:nrnod_macro))
      nrdofE_macro = sum(NdofE_macro(1:nrnod_macro))
      nrdofV_macro = sum(NdofV_macro(1:nrnod_macro))

!  ...H1 first
!  ...loop through macro element H1 dof
!       do kH = 1, nrdofH_macro
!          write(*,1011) kH , nrconH_macro(kH)
!  1011    format('print_constraints: kH,nrconH_macro(kH)            =' ,2(i6,3x))
!          write(*,1012) nacH_macro(1:nrconH_macro(kH),kH)
!  1012    format('print_constraints: nacH(1:nrconH_macro(kH),kH)    =', i6)
!          write(*,1013) constrH_macro(1:nrconH_macro(kH),kH)
!  1013    format('print_constraints: constrH_macro(1:nrconH(kH),kH) =',5(4e13.5,/))
!          call pause
! !  ...end of loop through macro element H1 dof
!       enddo
      do kE = 1, nrdofE_macro
         write(*,1014) kE , nrconE_macro(kE)
 1014    format('print_constraints: kE,nrconE_macro(kE)            =' ,2(i6,3x))
         write(*,1015) nacE_macro(1:nrconE_macro(kE),kE)
 1015    format('print_constraints: nacE(1:nrconE_macro(kE),kE)    =', i6)
         write(*,1016) constrE_macro(1:nrconE_macro(kE),kE)
 1016    format('print_constraints: constrE_macro(1:nrconE(kE),kE) =',5(4e13.5,/))
         call pause
!  ...end of loop through macro element H(curl) dof
      enddo
!       do kV = 1, nrdofV_macro
!          write(*,1017) kV, nrconV_macro(kV)
!  1017    format('print_constraints: kV,nrconV_macro(kV)            =' ,2(i6,3x))
!          write(*,1018) nacV_macro(1:nrconV_macro(kV),kV)
!  1018    format('print_constraints: nacV(1:nrconH_macro(kV),kV)    =', i6)
!          write(*,1019) constrV_macro(1:nrconV_macro(kV),kV)
!  1019    format('print_constraints: constrV_macro(1:nrconV(kV),kV) =',5(4e13.5,/))
!          call pause
! !  ...end of loop through macro element H(div) dof
!       enddo

      call pause
!
      deallocate(nod_macro, ndofH_macro, ndofE_macro, ndofV_macro)
      deallocate(nodm, ndofmH, ndofmE, ndofmV)
!    
!..end of loop through coarse elements    
   enddo
!
!
   end subroutine print_constraints

