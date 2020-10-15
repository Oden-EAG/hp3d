!---------------------------------------------------------------------------------------
!> Purpose : routine computes error for a physical attribute, in the  appropriate energy
!            space.
!
!> param[in] Flag - vector of length NR_PHYSA indicating for which
!                   attribute the error should be computed. The error
!                   is accumulated over components and rhs's.
!
!> param[in] Itag - tag to identify problem being tested
!
!> rev@Jun 2020 - ADDED MPI SUPPORT
!---------------------------------------------------------------------------------------
!  REMARK : miracles do not happen! Routine "exact" providing the exact
!           solution, should use the following ordering:
!
!  H1 :
!
!    comp1|comp2|... comp1|comp2|... ... |comp1|comp2|... comp1|comp2|... ...
!    attr1          |attr2           ... |attr1          |attr2           ...
!    rhs1                                |rhs2                            ...
!
!  H(curl), H(div), L2 : same as above
!---------------------------------------------------------------------------------------
!  REMARK : do NOT change field delimiter ";" in dump file! It is used by
!           testing script!
!---------------------------------------------------------------------------------------
#include "typedefs.h"
!
subroutine compute_error_modif(Flag,Itag)
!
      use control          , only : NEXACT
      use data_structure3D
      use environment      , only : QUIET_MODE,L2PROJ,FILE_ERR
      use physics
      use assembly_sc, only: NRDOF_TOT,NRDOF_CON
      use par_mesh   , only: DISTRIBUTED,HOST_MESH
      use mpi_param  , only: ROOT,RANK
      use MPI        , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8
!
      implicit none
      integer, dimension(NR_PHYSA),intent(in) :: Flag
      integer,                     intent(in) :: Itag
!
      real(8) :: errorH,errorE,errorV,errorQ,errorHEVQ,derrorH,derrorE,derrorV,derrorQ
      real(8) :: rnormH,rnormE,rnormV,rnormQ,rnormHEVQ,drnormH,drnormE,drnormV,drnormQ
      real(8) :: errorH_subd,errorE_subd,errorV_subd,errorQ_subd, &
                 rnormH_subd,rnormE_subd,rnormV_subd,rnormQ_subd
      real(8) :: errorH_rel,errorE_rel,errorV_rel,errorQ_rel,errorHEVQ_rel
      real(8) :: rateH,rateE,rateV,rateQ,rateHEVQ
!
      integer, parameter :: nin = 13
      integer, parameter :: maxvis =2000
!
!     for computing error rate
      integer, save :: ivis = 0
      integer, save :: nrdof_tot_save
      real(8), save :: errorH_save,errorE_save,errorV_save,errorQ_save,errorHEVQ_save
      real(8), dimension(maxvis,10), save :: rwork
      integer, dimension(maxvis,10), save :: iwork
!
!     miscellanea
      integer :: iel,mdle,i,iattr,ic,count,ierr
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
!     check that exact solution is indeed known
      if (NEXACT == 0) then
        write(*,*) 'compute_error: UNKNOWN exact solution!'
        return
      endif
!
!     initialize global quantities
      errorH=0.d0 ; rnormH=0.d0
      errorE=0.d0 ; rnormE=0.d0
      errorV=0.d0 ; rnormV=0.d0
      errorQ=0.d0 ; rnormQ=0.d0
!     initialize subdomain quantities
      errorH_subd=0.d0 ; rnormH_subd=0.d0
      errorE_subd=0.d0 ; rnormE_subd=0.d0
      errorV_subd=0.d0 ; rnormV_subd=0.d0
      errorQ_subd=0.d0 ; rnormQ_subd=0.d0
!
!     loop over active elements
!$OMP PARALLEL DO                                &
!$OMP PRIVATE(derrorH,derrorE,derrorV,derrorQ,   &
!$OMP         drnormH,drnormE,drnormV,drnormQ)   &
!$OMP REDUCTION(+:errorH_subd,rnormH_subd,errorE_subd,rnormE_subd, &
!$OMP             errorV_subd,rnormV_subd,errorQ_subd,rnormQ_subd)
   do iel=1,NRELES_SUBD
      call element_error_modif(ELEM_SUBD(iel),Flag,           &
                         derrorH,derrorE,derrorV,derrorQ,    &
                         drnormH,drnormE,drnormV,drnormQ)
      errorH_subd = errorH_subd + derrorH
      rnormH_subd = rnormH_subd + drnormH
      errorE_subd = errorE_subd + derrorE
      rnormE_subd = rnormE_subd + drnormE
      errorV_subd = errorV_subd + derrorV
      rnormV_subd = rnormV_subd + drnormV
      errorQ_subd = errorQ_subd + derrorQ
      rnormQ_subd = rnormQ_subd + drnormQ
   enddo
!$OMP END PARALLEL DO
! 
      if (DISTRIBUTED .and. (.not. HOST_MESH)) then
        count = 1
        call MPI_REDUCE(errorH_subd,errorH,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rnormH_subd,rnormH,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(errorE_subd,errorE,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rnormE_subd,rnormE,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(errorV_subd,errorV,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rnormV_subd,rnormV,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(errorQ_subd,errorQ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(rnormQ_subd,rnormQ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      else
        errorH = errorH_subd
        rnormH = rnormH_subd
        errorE = errorE_subd
        rnormE = rnormE_subd
        errorV = errorV_subd
        rnormV = rnormV_subd
        errorQ = errorQ_subd
        rnormQ = rnormQ_subd
      endif

!
!     compute total error
      errorHEVQ=sqrt(errorH+errorE+errorV+errorQ)
      errorH   =sqrt(errorH)
      errorE   =sqrt(errorE)
      errorV   =sqrt(errorV)
      errorQ   =sqrt(errorQ)
!
!     compute norm
      rnormHEVQ=sqrt(rnormH+rnormE+rnormV+rnormQ)
      rnormH   =sqrt(rnormH)
      rnormE   =sqrt(rnormE)
      rnormV   =sqrt(rnormV)
      rnormQ   =sqrt(rnormQ)
!
!     compute relative error
      errorHEVQ_rel=0.d0 ; errorH_rel=0.d0 ; errorE_rel=0.d0 ; errorV_rel=0.d0 ; errorQ_rel=0.d0
      if (rnormHEVQ > 0.d0) errorHEVQ_rel=errorHEVQ/rnormHEVQ
      if (rnormH    > 0.d0) errorH_rel   =errorH   /rnormH
      if (rnormE    > 0.d0) errorE_rel   =errorE   /rnormE
      if (rnormV    > 0.d0) errorV_rel   =errorV   /rnormV
      if (rnormQ    > 0.d0) errorQ_rel   =errorQ   /rnormQ
!
!     compute (a rough) total number of dof
      NRDOF_TOT=0
!
!     compute rate
      rateHEVQ=0.d0 ; rateH=0.d0 ; rateE=0.d0 ; rateV=0.d0 ; rateQ=0.d0
      if (ivis /= 0) then
        if (NRDOF_TOT > nrdof_tot_save) then
          rateHEVQ = (log(errorHEVQ_save/errorHEVQ))/log(float(nrdof_tot_save)/float(NRDOF_TOT))
      endif ; endif
!
!     save quantities
      errorHEVQ_save=errorHEVQ
      errorH_save=errorH ; errorE_save=errorE ; errorV_save=errorV ; errorQ_save=errorQ
      nrdof_tot_save=NRDOF_TOT
!
!
      if (RANK .eq. ROOT) then
!     raise visitation flag
        ivis=ivis+1
!
        IF (.NOT. QUIET_MODE) THEN
!
!     check
          if (ivis > maxvis) then
            write(*,*) 'compute_error: increase maxvis!'
            stop
          endif
!
!     store
          rwork(ivis, 1)=errorH
          rwork(ivis, 2)=errorE
          rwork(ivis, 3)=errorV
          rwork(ivis, 4)=errorQ
          rwork(ivis, 5)=errorHEVQ
          rwork(ivis, 6)=rateHEVQ
          iwork(ivis, 1)=Itag
!
        ENDIF
!
!     printing
!
!     -- 1st visit --
        if (ivis == 1) then
!
!       open file
          open(unit=nin,file=trim(FILE_ERR),form='formatted',access='sequential',status='unknown',iostat=ic)
          if (ic /= 0) then
            write(*,*)'compute_error: COULD NOT OPEN FILE! [0]'
            stop
          endif

!       print header
          IF (.NOT. L2PROJ) THEN
            write(nin,*)'-- Error Report --'
          ELSE
            write(nin,*)'-- Error Report (L2 only)--'
          ENDIF
          write(nin,9998)
 9998     format('             H1            //', &
                             ' H(curl)       //', &
                             ' H(div)        //', &
                             ' L2            //', &
                             ' Total         //', &
                             ' Rate       //',    &
                             ' Case tag')
!
!     -- subsequent visits --
        else
!
!       append to file
          open(unit=nin,file=trim(FILE_ERR),form='formatted',access='sequential',status='old',position='append',iostat=ic)
          if (ic /= 0) then
            write(*,*)'compute_error: COULD NOT OPEN FILE! [1]'
            stop
          endif
        endif
!
!     print to file
        write(nin,9999)ivis,errorH,errorE,errorV,errorQ,errorHEVQ,rateHEVQ,Itag
 9999   format(1x,i6,' ; ',2x,5(e12.5,' ; ',2x),f9.6,' ; ',3x,i8)
!
!     print to screen
        IF (.NOT.QUIET_MODE) THEN ; write(*,*)''
        IF (.NOT. L2PROJ   ) THEN ; write(*,*)'-- Error Report --'
        ELSE                      ; write(*,*)'-- Error Report (L2 only)--'
        ENDIF
                              write(*,9998)
        do i=1,ivis         ; write(*,9999)i,rwork(i,1:6),iwork(i,1)
        enddo
                              write(*,*)''
      ENDIF
!
!     close file
      close(unit=nin,iostat=ic)
      if (ic /= 0) then
        write(*,*)'compute_error: COULD NOT CLOSE FILE!'
        stop
      endif
!
    endif
    call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
end subroutine compute_error_modif
!
!
!
!---------------------------------------------------------------------------------------
!>  Purpose : compute element contributions to the total error
!!
!>  @date : Nov 2014
!---------------------------------------------------------------------------------------
!
subroutine element_error_modif(Mdle,Flag, errorH,errorE,errorV,errorQ, rnormH,rnormE,rnormV,rnormQ)
!
      use control          , only : INTEGRATION
      use data_structure3D
      use environment      , only : L2PROJ
      use physics
!
      implicit none
      integer, dimension(NR_PHYSA),intent(in ) :: Flag
      integer,                     intent(in ) :: Mdle
      real(8),                      intent(out) :: errorH,errorE,errorV,errorQ
      real(8),                      intent(out) :: rnormH,rnormE,rnormV,rnormQ
!
!     node case (decimal form)
      integer,dimension(NR_PHYSA) :: icased
!
!     element, face order, geometry dof
      integer,dimension(19)          :: norder
      real(8),dimension(3,MAXbrickH) :: xnod
      integer,dimension(12)          :: nedge_orient
      integer,dimension(6)           :: nface_orient
!
!     geometry
      real(8),dimension(3)   :: xi,x
      real(8),dimension(3,3) :: dxidx,dxdxi
      real(8)                :: rjac
!
!     3D quadrature data
      real(8),dimension(3,MAX_NINT3) :: xiloc
      real(8),dimension(  MAX_NINT3) :: wxi
!
!     approximate solution dof's
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!     approximate solution
      VTYPE, dimension(  MAXEQNH  ) ::  zsolH
      VTYPE, dimension(  MAXEQNH,3) :: zdsolH
      VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
      VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
      VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
      VTYPE, dimension(  MAXEQNV  ) ::  zdivV
      VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!     exact solution
      VTYPE, dimension(  MAXEQNH    ) ::   zvalH
      VTYPE, dimension(  MAXEQNH,3  ) ::  zdvalH
      VTYPE, dimension(  MAXEQNH,3,3) :: zd2valH
      VTYPE, dimension(3,MAXEQNE    ) ::   zvalE
      VTYPE, dimension(3,MAXEQNE,3  ) ::  zdvalE
      VTYPE, dimension(3,MAXEQNE,3,3) :: zd2valE
      VTYPE, dimension(3,MAXEQNE    ) ::  zcurlE_ex
      VTYPE, dimension(3,MAXEQNV    ) ::   zvalV
      VTYPE, dimension(3,MAXEQNV,3  ) ::  zdvalV
      VTYPE, dimension(3,MAXEQNV,3,3) :: zd2valV
      VTYPE, dimension(  MAXEQNV    ) ::   zdivV_ex
      VTYPE, dimension(  MAXEQNQ    ) ::   zvalQ
      VTYPE, dimension(  MAXEQNQ,3  ) ::  zdvalQ
      VTYPE, dimension(  MAXEQNQ,3,3) :: zd2valQ
!
!     miscellanea
      integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag
      real(8) :: weight,wa
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
!     initialize global quantities
      errorH=0.d0 ; rnormH=0.d0
      errorE=0.d0 ; rnormE=0.d0
      errorV=0.d0 ; rnormV=0.d0
      errorQ=0.d0 ; rnormQ=0.d0
!
!     order of approx, orientations, geometry dof's, solution dof's
      call find_order( mdle, norder)
      call find_orient(mdle, nedge_orient,nface_orient)
      call nodcor(     mdle, xnod)
      call solelm(     mdle, zdofH,zdofE,zdofV,zdofQ)
!
!     set up the element quadrature
      INTEGRATION=2
      call set_3Dint(NODES(mdle)%type,norder, nint,xiloc,wxi)
      INTEGRATION=0
!
!     supported physical attributes
      icase=NODES(mdle)%case
      call decod(icase,2,NR_PHYSA, icased)
!
!     loop over physical attributes
      do iattr=1,NR_PHYSA
!
!       if the error not needed, skip
        if (Flag(iattr) == 0) cycle
!
!       if attribute is absent, skip
        if (icased(iattr) == 0) cycle
!
!       address of the 1st component for the attribute
        ibeg=ADRES(iattr)
!
        select case(DTYPE(iattr))
!
!===================================================================================
!  H1 ATTRIBUTE                                                                    |
!===================================================================================
          case('contin')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           -- EXACT SOLUTION --
            call exact(x,Mdle,icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                                zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
9997          format(' element_error: mdle,ndom,rjac = ',i8,2x,i2,2x,e12.5)
            endif
!
!           total weight
            weight=wa*rjac
!
!           loop over rhs's
            do iload=1,NRCOMS
!
!             loop over components of physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(iload-1)*NRHVAR+ibeg+icomp
!
!               accumulate H1 seminorm
                IF (.NOT. L2PROJ) THEN
                do j=1,3
                  errorH = errorH + abs(zdvalH(i,j) - zdsolH(i,j))**2 * weight
                  rnormH = rnormH + abs(zdvalH(i,j)              )**2 * weight
                enddo
                ENDIF
!
!               accumulate L2 norm
                errorH = errorH + abs(zvalH(i) - zsolH(i))**2 * weight
                rnormH = rnormH + abs(zvalH(i)           )**2 * weight
!
!             loop over components
              enddo
!
!           loop over rhs's
            enddo
!
!         loop over integration points
          enddo
!
!===================================================================================
!  H(curl) ATTRIBUTE                                                               |
!===================================================================================
          case('tangen')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           -- EXACT SOLUTION --
            call exact(x,Mdle,icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                                zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
!
!           curl of exact solution
            zcurlE_ex(1,1:MAXEQNE) = zdvalE(3,1:MAXEQNE,2) - zdvalE(2,1:MAXEQNE,3)
            zcurlE_ex(2,1:MAXEQNE) = zdvalE(1,1:MAXEQNE,3) - zdvalE(3,1:MAXEQNE,1)
            zcurlE_ex(3,1:MAXEQNE) = zdvalE(2,1:MAXEQNE,1) - zdvalE(1,1:MAXEQNE,2)
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
!
!           total weight
            weight=wa*rjac
!
!           loop over rhs's
            do iload=1,NRCOMS
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(iload-1)*NREVAR+ibeg+icomp
!
!               accumulate for the error and norm
                do ivar=1,3
                  ErrorE = ErrorE + abs(zsolE(ivar,i) - zvalE(ivar,i))**2 * weight
                  RnormE = RnormE + abs(                zvalE(ivar,i))**2 * weight
                  IF (.NOT. L2PROJ) THEN
                  ErrorE = ErrorE + abs(zcurlE(ivar,i) - zcurlE_ex(ivar,i))**2 * weight
                  RnormE = RnormE + abs(                 zcurlE_ex(ivar,i))**2 * weight
                  ENDIF
                enddo
!
!             loop over components
              enddo
!
!           loop over rhs's
            enddo
!
!         loop over integration points
          enddo
!
!===================================================================================
!  H(div) ATTRIBUTE                                                                |
!===================================================================================
          case('normal')
!
!         loop over integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           -- EXACT SOLUTION --
            call exact(x,Mdle,icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                                zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
!
!           divergence of exact solution
            zdivV_ex(1:MAXEQNV) = zdvalV(1,1:MAXEQNV,1) + &
                                  zdvalV(2,1:MAXEQNV,2) + &
                                  zdvalV(3,1:MAXEQNV,3)
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
!
!           total weight
            weight=wa*rjac
!
!           loop over rhs's
            do iload=1,NRCOMS
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(iload-1)*NRVVAR+ibeg+icomp
!
!               accumulate H(div) seminorm
                IF (.NOT. L2PROJ) THEN
                ErrorV = ErrorV + abs(zdivV_ex(i) - zdivV(i))**2 * weight
                RnormV = RnormV + abs(zdivV_ex(i)           )**2 * weight
                ENDIF
!
!               accumulate L2 norm
                do ivar=1,3
                  ErrorV = ErrorV + abs(zvalV(ivar,i) - zsolV(ivar,i))**2 * weight
                  RnormV = RnormV + abs(zvalV(ivar,i)                )**2 * weight
                enddo
!
              enddo
            enddo
!
!         loop over integration points
          enddo
!
!===================================================================================
!  L2 ATTRIBUTE                                                                    |
!===================================================================================
          case('discon')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           -- EXACT SOLUTION --
            call exact(x,Mdle,icase, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, &
                                zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif

!           total weight
            weight=wa*rjac
!
!           loop over rhs's
            do iload=1,NRCOMS
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(iload-1)*NRQVAR+ibeg+icomp
!
!               accumulate L2 norm
                ErrorQ = ErrorQ + abs(zvalQ(i) - zsolQ(i))**2 * weight
                RnormQ = RnormQ + abs(zvalQ(i)           )**2 * weight
!
              enddo
            enddo
!
!         loop over integration points
          enddo
!
          case default
          write(*,*)'element_error: UNKNOWN PHYSICAL ATTRIBUTE TYPE!'
          stop
!
          endselect
!
!       loop over physical attributes
        enddo
!
!
end subroutine element_error_modif
