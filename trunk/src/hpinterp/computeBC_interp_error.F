c--------------------------------------------------------------------
c                                                                     
c     routine name      - compute_BC_interp_error
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Mar 17
c                                                                     
c     purpose:          - routine returns global BC interpolation
c                         error for traces of BC data 
c                         (a verification for interpolating 
c                          and exact solution routines)
c                                                                    
c---------------------------------------------------------------------
c
      subroutine computeBC_interp_error
c
      use data_structure3D
      use environment      , only : QUIET_MODE
#include "syscom.blk"
c
c  ...visitation flag
      integer, save :: ivis = 0
c
c  ...errors for the old mesh
      real*8,  save :: errorHold, errorEold, errorVold
c
c  ...residuals and rates to display
      real*8 , dimension(6,10), save :: rwork
c
      iprint=0
c
c  ...compute total errors and number of dof
      errorH =  0.d0; errorE = 0.d0; errorV = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elemBCIerror(mdle, derrorH,derrorE,derrorV)
        if ((derrorH.eq.0.d0).and.(derrorE.eq.0.d0).and.
     .      (derrorV.eq.0.d0)) cycle
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, derrorH,derrorE,derrorV
 7010     format('compute_interp_error: iel, mdle = ',2i5,
     .           ' element traces errors = ',3e12.5)
        endif
        errorH = errorH + derrorH
        errorE = errorE + derrorE
        errorV = errorV + derrorV
      enddo
      errorH = sqrt(errorH) 
      errorE = sqrt(errorE)
      errorV = sqrt(errorV)
c
c  ...compute rate
      rateH = 0.d0; rateE = 0.d0; rateV = 0.d0 
      if (ivis.ne.0) then
        if (errorH.ne.0.d0.and.errorHold.ne.0.d0) then
          rateH = log(errorHold/errorH)/log(2.d0)
        endif
        if (errorE.ne.0.d0.and.errorEold.ne.0.d0) then
          rateE = log(errorEold/errorE)/log(2.d0)
        endif
        if (errorV.ne.0.d0.and.errorVold.ne.0.d0) then
          rateV = log(errorVold/errorV)/log(2.d0)
        endif
      endif
c
c  ...save current data
      ivis = ivis+1
      errorHold = errorH
      errorEold = errorE
      errorVold = errorV
c
c  ...store data to display
      rwork(1,ivis) = errorH; rwork(2,ivis) = rateH
      rwork(3,ivis) = errorE; rwork(4,ivis) = rateE
      rwork(5,ivis) = errorV; rwork(5,ivis) = rateV
c
c  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Error Report --'
        write(*,7100)
 7100   format(' Mesh,      errorH,  rateH,'
     .         '       errorE,  rateE ',
     .         '       errorV,  rateV ')
        do i=1,ivis     
          write(*,7110)i,rwork(1:6,i)
 7110     format(i3,3x,3(e12.5,f8.2,2x))
        enddo
        write(*,*)''
      endif
c
c
      end subroutine computeBC_interp_error





