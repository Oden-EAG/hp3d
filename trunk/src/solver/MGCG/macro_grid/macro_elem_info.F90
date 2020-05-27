!
! -----------------------------------------------------------------------
!
!    routine name       - macro elem_info
!
! -----------------------------------------------------------------------
!
!    latest revision    - Feb 2018
!
!    purpose            - returns the macro-element information
!
!   arguments :
!     in:
!             MdleC    - middle node of a coarse element
!     out:     
!             nrdof    - number of condensed element dof 
!         nrnod_macro  - number of nodes
!           nod_macro  - node list
!         ndofmH_macro - the corresponding number of H1 dof
!         ndofmE_macro - the corresponding number of H(curl) dof
!         ndofmV_macro - the corresponding number of H(div) dof
!
! ----------------------------------------------------------------------
!
   subroutine macro_elem_info(Igrid, MdleC, Nrdof,Nrdofb, Nrnod_macro, Nod_macro,  &
                              NdofmH_macro, NdofmE_macro, NdofmV_macro)
!
   use data_structure3D, ONLY: NODES,MAXNODM, NRNODS,NR_PHYSA
   use assembly,         ONLY: assembly_begin_par, assembly_end_par
   use stc,              ONLY: stc_get_nrdof
!
   implicit none
!   
!-----------------------------------------------------------------------
!
   integer, intent(in)   :: Igrid, MdleC, Nrnod_macro
   integer, intent(out)  :: Nrdof, Nrdofb,                                         &
                            nod_macro(Nrnod_macro),    ndofmH_macro(Nrnod_macro),  &
                            ndofmE_macro(Nrnod_macro), ndofmV_macro(Nrnod_macro)
!
!..number of local element dof for each physics variable
   integer, dimension(NR_PHYSA) :: nrdof_i,nrdof_b
!
!..locals
!..work space for celem
   integer              :: nodm(MAXNODM),  ndofmH(MAXNODM),  &
                           ndofmE(MAXNODM),ndofmV(MAXNODM), ndofmQ(MAXNODM)
   character(len=4)     :: type
   integer, allocatable :: nvisitH(:),nvisitE(:), nvisitV(:), nvisit(:)
!   
   integer              :: nrdofs(NR_PHYSA)
   integer              :: nrdof_H,  nrdof_E,  nrdof_V
   integer              :: nrdof_Hb, nrdof_Eb, nrdof_Vb 
   integer              :: mdle, imdle, nrdofc
   integer              :: nrnodm, nrdofm, nrdof_mdl
   complex*16           :: zvoid(1)
   integer              :: i, nod, master, inod
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!..allocate required variables for celem
   call assembly_begin_par
!
   allocate(nvisitH(NRNODS)); nvisitH = -1
   allocate(nvisitE(NRNODS)); nvisitE = -1
   allocate(nvisitV(NRNODS)); nvisitV = -1
   allocate(nvisit(NRNODS)) ; nvisit  = -1

!
!..Step 1: determine the first dof offsets for active nodes
   nrdof_H  = 0 ; nrdof_E  = 0; nrdof_V  = 0
   nrdof_Hb = 0 ; nrdof_Eb = 0; nrdof_Vb = 0
   nrdof_mdl = 0
   mdle=0 ;   inod = 0 ;   imdle = 1
   call nelcon_macro(MdleC,mdle, mdle)
   if (mdle .ne. 0) then
      imdle = 0
   endif
   do while (mdle.ne.0)
      imdle = imdle +1
      call celem_mg(-1,-1,mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                    ndofmE,ndofmV,ndofmQ,nrnodm, zvoid,zvoid)
!
!  ...fill in the lists for macro element nodes      
      do i = 1,nrnodm ! no need to check the mdle node
         nod = nodm(i)
         if (nvisit(nod) .gt. 0) cycle
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check if the the type of the master node 
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
         case default
!        ...update the macro nod counter            
            inod = inod + 1 
!        ...fill in the lists            
            nod_macro(inod) = nod
            ndofmH_macro(inod) = ndofmH(i)
            ndofmE_macro(inod) = ndofmE(i)
            ndofmV_macro(inod) = ndofmV(i)
!            
!        ...raise the visitation flag            
            nvisit(nod) = 1
         end select
      enddo
!  ...Step 1: Compute offsets
!  ...H1 dof
      do i = 1,nrnodm
         nod = nodm(i)
         if (nvisitH(nod).ge.0) cycle

!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check if the the type of the master node 
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!        ...store the first dof offset
            nvisitH(nod) = nrdof_Hb
!        ...update the H1 dof counter
            nrdof_Hb = nrdof_Hb + ndofmH(i)
         case default
!        ...store the first dof offset
            nvisitH(nod) = nrdof_H
!        ...update the H1 dof counter
            nrdof_H = nrdof_H + ndofmH(i)
         end select
      enddo
!
!  ...H(curl) dof
      do i = 1,nrnodm
         nod = nodm(i)
         if (nvisitE(nod).ge.0) cycle
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check if the the type of the master node 
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!        ...store the first dof offset
            nvisitE(nod) = nrdof_Eb
!        ...update the H(curl) dof counter
            nrdof_Eb = nrdof_Eb + ndofmE(i)
         case default
!        ...store the first dof offset
            nvisitE(nod) = nrdof_E
!        ...update the H(curl) dof counter
            nrdof_E = nrdof_E + ndofmE(i)
         end select
      enddo
! 
!  ...H(div) dof
      do i = 1,nrnodm
         nod = nodm(i)
         if (nvisitV(nod).ge.0) cycle
!     ...find its coarse master
         call find_master(Igrid,nod, master)
!     ...check if the the type of the master node 
         type = NODES(master)%type
         select case(type)
         case('mdlb','mdln','mdlp','mdld')
!        ...store the first dof offset
            nvisitV(nod) = nrdof_Vb
!        ...update the H(div) dof counter
            nrdof_Vb = nrdof_Vb + ndofmV(i)
         case default
!        ...store the first dof offset
            nvisitV(nod) = nrdof_V
!        ...update the H(div) dof counter
            nrdof_V = nrdof_V + ndofmV(i)
         end select
      enddo
!
      call stc_get_nrdof(mdle, nrdof_i,nrdof_b)
      nrdof_mdl = nrdof_mdl + sum(nrdof_b)
!
!....end of loop through elements
      call nelcon_macro(MdleC,mdle, mdle)
   enddo
!
   Nrdof = nrdof_H + nrdof_E + nrdof_V
   Nrdofb = nrdof_Hb + nrdof_Eb + nrdof_Vb + nrdof_mdl
!   
   call assembly_end_par
   deallocate(nvisitH,nvisitE,nvisitV, nvisit)
! 
! 
   end subroutine macro_elem_info


