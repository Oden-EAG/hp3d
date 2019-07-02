

   subroutine test_logic(Mdle,idec)

   use data_structure3D

   IMPLICIT NONE
!
!----------------------------------------------------------------------------
!
   integer :: nodm(MAXNODM)
   integer :: ndofmHl(MAXNODM),ndofmEl(MAXNODM),ndofmVl(MAXNODM)
   integer :: nrnodm, idec, mdle
   integer :: nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV), &
              nacH(NACDIM,MAXbrickH), nacE(NACDIM,MAXbrickE),     &
              nacV(NACDIM,MAXbrickV)

   real*8  :: constrH(NACDIM,MAXbrickH), constrE(NACDIM,MAXbrickE),  &
              constrV(NACDIM,MAXbrickV)
   integer :: iprint
!
!----------------------------------------------------------------------------
!
   iprint = 1
!
!   
   select case(idec)
   case(1)
      call logicC(Mdle,idec,                             &
                  nodm,ndofmHl,ndofmEl,ndofmVl,nrnodm,   &
                  nrconH,nacH,constrH,                   &
                  nrconE,nacE,constrE,                   &
                  nrconV,nacV,constrV)
   case default
      call logic(Mdle,idec,                             &
                 nodm,ndofmHl,ndofmEl,ndofmVl,nrnodm,   &
                 nrconH,nacH,constrH,                   &
                 nrconE,nacE,constrE,                   &
                 nrconV,nacV,constrV)
   end select   

   if (iprint .eq. 1) then 
      select case(idec)
      case(1)
         write(*,*) 'TESTING logicC'
      case default   
         write(*,*) 'TESTING logic'
      end select

      write(*,2000) Mdle
 2000 format('test_logic: MdleC = ', i7)
      write(*,2001) nodm(1:nrnodm)
 2001 format('test_logic: nodm = ',/, 10(8i7,/))
      write(*,2002) ndofmHl(1:nrnodm)
 2002 format('test_logicC: ndofmHl = ',/, 10(8i7,/))
      write(*,2003) ndofmEl(1:nrnodm)
 2003 format('test_logicC: ndofmEl = ',/, 10(8i7,/))
      write(*,2004) ndofmHl(1:nrnodm)
 2004 format('test_logicC: ndofmVl = ',/, 10(8i7,/))
   call pause
   endif



   end subroutine test_logic