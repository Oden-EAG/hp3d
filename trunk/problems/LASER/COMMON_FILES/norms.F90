!
#include "typedefs.h"
!
!-------------------------------------------------------------------
!  subroutine: get_L2NormCOMS
!
!  last modified : Jan 2019
!
!  purpose: - evaluates the L2 norm of the difference of components
!             No1 and No2 (as specified in NRCOMS) of L2 variable
!             [used to evaluate
!
!  input:   - Flag
!           - No1,No2
!
!  output:  - L2NormDiff
!
!-------------------------------------------------------------------
subroutine get_L2NormCOMS(Flag,No1,No2, L2NormDiff)
!
   use control         , only: NEXACT
   use data_structure3D
   use environment     , only: QUIET_MODE,L2PROJ,FILE_ERR
   use physics
   use par_mesh        , only: DISTRIBUTED,HOST_MESH
   use mpi_param       , only: ROOT,RANK
   use MPI             , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8,MPI_Wtime
!
   implicit none
!
   integer, intent(in)  :: Flag(NR_PHYSA)
   integer, intent(in)  :: No1, No2
   real(8), intent(out) :: L2NormDiff
!
   real(8)              :: currL2NormDiff
!
   integer, parameter :: nin = 13
   integer, parameter :: maxvis =2000
!
!..miscellanea
   integer :: mdle,nint,iattr,nrdof_tot,ic,iel
!
!..auxiliary
   integer :: count,ierr
   real(8) :: L2NormDiff_subd
!
   real(8) :: start_time,end_time
!
!-------------------------------------------------------------------
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..fetch active elements
   if (.not. DISTRIBUTED .or. HOST_MESH) then
      if (RANK .ne. ROOT) goto 90
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..initialize global quantities
   L2NormDiff_subd = 0.d0
   currL2NormDiff = 0.d0
!
!..iterate over elements
!
!$OMP PARALLEL DO                    &
!$OMP PRIVATE(mdle,currL2NormDiff)   &
!$OMP REDUCTION(+:L2normDiff_subd)   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call get_elem_L2NormCOMS(mdle,Flag,No1,No2, currL2NormDiff)
!  ...accumulate L2NormDiff
      L2NormDiff_subd = L2NormDiff_subd + currL2NormDiff
   enddo
!$OMP END PARALLEL DO
!
   L2NormDiff = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_ALLREDUCE(L2NormDiff_subd,L2NormDiff,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   else
      L2NormDiff = L2NormDiff_subd
   endif
!
 90 continue
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (RANK.eq.ROOT .and. .not. QUIET_MODE) then
      write(*,6010) end_time-start_time
 6010 format('  get_L2NormCOMS : ',f12.5,'  seconds',/)
   endif
!
!..compute sqrt of total difference
   L2NormDiff = sqrt(L2NormDiff)
!
end subroutine get_L2NormCOMS
!
!
!-------------------------------------------------------------------
!
!  routine: get_elem_L2NormCOMS, TODO: VERIFY CORRECTNESS
!
!  purpose: to evaluate the L2 norm of the difference between components
!           No1 and No2 with components (set using NRCOMS) of L2 variables
!
!  last modified: Jan 2019
!
!  input:   - mdle,Flag,No1,No2
!
!  output:  - FieldNormQ: L2 norm of difference
!
!-------------------------------------------------------------------
subroutine get_elem_L2NormCOMS(mdle,Flag,No1,No2, FieldNormQ)
!
   use laserParam
   use commonParam
!
   use control          , only : INTEGRATION
   use data_structure3D
   use environment      , only : L2PROJ
   use physics
!
   implicit none
!
   integer, dimension(NR_PHYSA),intent(in ) :: Flag
   integer,                     intent(in ) :: mdle,No1,No2
   real(8),                     intent(out) :: FieldNormQ
!
!..node case (decimal form)
   integer,dimension(NR_PHYSA) :: icased
!
!..element, face order, geometry dof
   integer,dimension(19)           :: norder
   real(8) ,dimension(3,MAXbrickH) :: xnod
   integer,dimension(12)           :: nedge_orient
   integer,dimension(6)            :: nface_orient
!
!..geometry
   real(8),dimension(3)   :: xi,x
   real(8),dimension(3,3) :: dxidx,dxdxi
   real(8)                :: rjac
!
!..3D quadrature data
   real(8),dimension(3,MAX_NINT3) :: xiloc
   real(8),dimension(  MAX_NINT3) :: wxi
!
!..approximate solution dof's
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!..approximate solution
   VTYPE, dimension(  MAXEQNH  ) ::  zsolH
   VTYPE, dimension(  MAXEQNH,3) :: zdsolH
   VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
   VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
   VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
   VTYPE, dimension(  MAXEQNV  ) ::  zdivV
   VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!..miscellanea
   integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag
   real(8) :: weight,wa
!
!---------------------------------------------------------------------------------------
!
   nflag=1
!
!..initialize global quantities
   FieldNormQ = 0.d0
!
!..order of approx, orientations, geometry dof's, solution dof's
   call find_order( mdle, norder)
   call find_orient(mdle, nedge_orient,nface_orient)
   call nodcor(     mdle, xnod)
   call solelm(     mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..set up the element quadrature
   INTEGRATION=0
   call set_3D_int_DPG(NODES(mdle)%ntype,norder,nface_orient, nint,xiloc,wxi)
   INTEGRATION=0
!
!..supported physical attributes
   icase=NODES(mdle)%case
   call decod(icase,2,NR_PHYSA, icased)
!
!..loop over physical attributes
   do iattr=1,NR_PHYSA
!
!  ...if the error not needed, skip
      if (Flag(iattr) == 0) cycle
!
!  ...if attribute is absent, skip
      if (icased(iattr) == 0) cycle
!
!  ...address of the 1st component for the attribute
      ibeg=ADRES(iattr)
!
      select case(D_TYPE(iattr))
!
!===================================================================================
!  L2 ATTRIBUTE                                                                    |
!===================================================================================
      case(DISCON)
!
!     ...loop through integration points
         do l=1,nint
!
!        ...Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                         zdofH,zdofE,zdofV,zdofQ,nflag,                 &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!        ...Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
#if DEBUG_MODE
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,*)' get_elem_L2NormCOMS: mdle,ndom,rjac = ', mdle,ndom,rjac
            endif
#endif

!        ...total weight
            weight=wa*rjac
!
!        ...loop over components of the physical attribute
            do icomp=1,NR_COMP(iattr)
!
               i=(No1-1)*NRQVAR+ibeg+icomp
               j=(No2-1)*NRQVAR+ibeg+icomp
!
!           ...accumulate L2 norm
               FieldNormQ = FieldNormQ+(abs(zsolQ(i)-zsolQ(j))**2)*weight
!
            enddo

!     ...end loop over integration points
         enddo
!
         case default
         write(*,*)' get_elem_L2NormCOMS: PHYSICAL ATTRIBUTE MUST BE L2!. stop.'
         stop
!
      end select
!
!..loop over physical attributes
   enddo
!
end subroutine get_elem_L2NormCOMS
!
!
!-------------------------------------------------------------------
!
!  routine: get_Norm, TODO: VERIFY CORRECTNESS
!
!  purpose: evaluates the norm of specified physical attributes
!           and specified component by integrating the current
!           solution over the entire domain
!
!  last modified: Jan 2019
!
!  input:   - Field, No
!
!  output:  - FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!             [(norm of field variable) = \int_{\Omega} |X|^2]
!
!-------------------------------------------------------------------
subroutine get_Norm(Flag,No, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
!
   use control          , only: NEXACT
   use data_structure3D
   use environment      , only: QUIET_MODE,L2PROJ,FILE_ERR
   use physics
   use par_mesh         , only: DISTRIBUTED,HOST_MESH
   use mpi_param        , only: ROOT,RANK
   use MPI              , only: MPI_SUM,MPI_COMM_WORLD,MPI_REAL8,MPI_Wtime
!
   implicit none
!
   integer, intent(in)  :: Flag(NR_PHYSA)
   integer, intent(in)  :: No
!
   real(8), intent(out) :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
   real(8) :: currFieldNormH,currFieldNormE,currFieldNormV,currFieldNormQ
!
   integer, parameter :: nin = 13
   integer, parameter :: maxvis =2000
!
!..miscellanea
   integer :: mdle,nint,iattr,nrdof_tot,ic,iel
!
!..auxiliary
   integer :: count, ierr
   real(8) :: fieldNormH_subd, fieldNormE_subd,  &
              fieldNormV_subd, fieldNormQ_subd
!
!..timer
   real(8) :: start_time,end_time
!
!-------------------------------------------------------------------
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..fetch active elements
   if (.not. DISTRIBUTED .or. HOST_MESH) then
      if (RANK .ne. ROOT) goto 90
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..initialize residual
   fieldNormH_subd = 0.d0
   fieldNormE_subd = 0.d0
   fieldNormV_subd = 0.d0
   fieldNormQ_subd = 0.d0
!
!..loop over active elements
!
!$OMP PARALLEL DO                                                                   &
!$OMP PRIVATE(mdle,currFieldNormH,currFieldNormE,currFieldNormV,currFieldNormQ)     &
!$OMP REDUCTION(+:FieldNormH_subd,FieldNormE_subd,FieldNormV_subd,FieldNormQ_subd)  &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call get_elem_Norm(mdle,Flag,No, &
         currFieldNormH,currFieldNormE,currFieldNormV,currFieldNormQ)
!  ...accumulate
      fieldNormH_subd = fieldNormH_subd + currFieldNormH
      fieldNormE_subd = fieldNormE_subd + currFieldNormE
      fieldNormV_subd = fieldNormV_subd + currFieldNormV
      fieldNormQ_subd = fieldNormQ_subd + currFieldNormQ
   enddo
!$OMP END PARALLEL DO
!
   FieldNormH = 0.d0
   FieldNormE = 0.d0
   FieldNormV = 0.d0
   FieldNormQ = 0.d0
!
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_ALLREDUCE(fieldNormH_subd,FieldNormH,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(fieldNormE_subd,FieldNormE,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(fieldNormV_subd,FieldNormV,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(fieldNormQ_subd,FieldNormQ,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   else
      FieldNormH = fieldNormH_subd
      FieldNormE = fieldNormE_subd
      FieldNormV = fieldNormV_subd
      FieldNormQ = fieldNormQ_subd
   endif
!
 90 continue
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (RANK.eq.ROOT .and. .not. QUIET_MODE) then
      write(*,9010) end_time-start_time
 9010 format('  get_Norm       : ',f12.5,'  seconds',/)
   endif
!
!..compute total error
   FieldNormH = sqrt(FieldNormH)
   FieldNormE = sqrt(FieldNormE)
   FieldNormV = sqrt(FieldNormV)
   FieldNormQ = sqrt(FieldNormQ)
!
end subroutine get_Norm
!
!
!-------------------------------------------------------------------
!
!  routine: get_elem_Norm, TODO: VERIFY CORRECTNESS
!
!  purpose: evaluates the magnitude of a physical attribute and
!           component by integrating the current solution over a
!           given middle node
!
!  last modified: Jan 2019
!
!  input:   Mdle,Flag,No
!
!  output:  FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!           (FieldNorm (norm of field variable))
!
!-------------------------------------------------------------------
subroutine get_elem_Norm(Mdle,Flag,No, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
!
   use laserParam
   use commonParam
!
   use control          , only : INTEGRATION
   use data_structure3D
   use environment      , only : L2PROJ
   use physics
!
   implicit none
!
   integer, intent(in)  :: Flag(NR_PHYSA)
   integer, intent(in)  :: Mdle,No
   real(8), intent(out) :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
!..node case (decimal form)
   integer :: icased(NR_PHYSA)
!
!..element, face order, geometry dof
   integer :: norder(19)
   real(8) :: xnod(3,MAXbrickH)
   integer :: nedge_orient(12)
   integer :: nface_orient(6)
!
!..geometry
   real(8), dimension(3)   :: xi,x
   real(8), dimension(3,3) :: dxidx,dxdxi
   real(8)                 :: rjac
!
!..3D quadrature data
   real(8) :: xiloc(3,MAX_NINT3)
   real(8) :: wxi(MAX_NINT3)
!
!..approximate solution dof's
   VTYPE :: zdofH(MAXEQNH,MAXbrickH)
   VTYPE :: zdofE(MAXEQNE,MAXbrickE)
   VTYPE :: zdofV(MAXEQNV,MAXbrickV)
   VTYPE :: zdofQ(MAXEQNQ,MAXbrickQ)
!
!..approximate solution
   VTYPE, dimension(  MAXEQNH  ) ::  zsolH
   VTYPE, dimension(  MAXEQNH,3) :: zdsolH
   VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
   VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
   VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
   VTYPE, dimension(  MAXEQNV  ) ::  zdivV
   VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!..miscellanea
   integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag
   real(8) :: weight,wa
!
!---------------------------------------------------------------------------------------
!
   nflag=1
!
!..initialize global quantities
   FieldNormH=0.d0
   FieldNormE=0.d0
   FieldNormV=0.d0
   FieldNormQ=0.d0
!
!..order of approx, orientations, geometry dof's, solution dof's
   call find_order( mdle, norder)
   call find_orient(mdle, nedge_orient,nface_orient)
   call nodcor(     mdle, xnod)
   call solelm(     mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..set up the element quadrature
   INTEGRATION=0
   call set_3D_int_DPG(NODES(mdle)%ntype,norder,nface_orient, nint,xiloc,wxi)
   INTEGRATION=0
!
!..supported physical attributes
   icase=NODES(mdle)%case
   call decod(icase,2,NR_PHYSA, icased)
!
!..loop over physical attributes
   do iattr=1,NR_PHYSA
!
!  ...if the error not needed, skip
      if (Flag(iattr) == 0) cycle
!
!  ...if attribute is absent, skip
      if (icased(iattr) == 0) cycle
!
!  ...address of the 1st component for the attribute
      ibeg=ADRES(iattr)
!
      select case(D_TYPE(iattr))
!
!===================================================================================
!  H1 ATTRIBUTE                                                                    |
!===================================================================================
      case(CONTIN)
!
!     ...loop through integration points
         do l=1,nint
!
!        ...Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                         zdofH,zdofE,zdofV,zdofQ,nflag,                 &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!        ...Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
#if DEBUG_MODE
            if (iflag /= 0) then
               call find_domain(mdle, ndom)
               write(*,9997) mdle,ndom,rjac
            endif
 9997       format(' get_elem_Norm: mdle,ndom,rjac = ',i8,2x,i2,2x,e12.5)
#endif
!
!        ...total weight
            weight = wa*rjac
!
!        ...loop over components of physical attribute
            do icomp=1,NR_COMP(iattr)
!
               i=(No-1)*NRHVAR+ibeg+icomp
!
!           ...accumulate H1 seminorm
               if (.not. L2PROJ) then
                  do j=1,3
                     FieldNormH = FieldNormH + abs(zdsolH(i,j))**2 * weight
                  enddo
               endif
!
!           ...accumulate L2 norm
               FieldNormH = FieldNormH + abs(zsolH(i))**2 * weight
!
!        ...end loop over components
            enddo
!
!     ...end loop over integration points
         enddo
!
!===================================================================================
!  H(curl) ATTRIBUTE                                                               |
!===================================================================================
      case(TANGEN)
!
!     ...loop through integration points
         do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                         zdofH,zdofE,zdofV,zdofQ,nflag,                 &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
#if DEBUG_MODE
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
#endif
!
!           total weight
            weight=wa*rjac
!
!             loop over the components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NREVAR+ibeg+icomp
!
!               accumulate for the error and norm
                do ivar=1,3
                  FieldNormE = FieldNormE + abs(zsolE(ivar,i))**2 * weight
                  IF (.NOT. L2PROJ) THEN
                  FieldNormE = FieldNormE + abs(zcurlE(ivar,i))**2 * weight
                  ENDIF
                enddo
!
!        ...end loop over components
            enddo
!
!     ...end loop over integration points
         enddo
!
!===================================================================================
!  H(div) ATTRIBUTE                                                                |
!===================================================================================
      case(NORMAL)
!
!     ...loop over integration points
         do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                         zdofH,zdofE,zdofV,zdofQ,nflag,                 &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
#if DEBUG_MODE
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
#endif
!
!           total weight
            weight=wa*rjac
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NRVVAR+ibeg+icomp
!
!               accumulate H(div) seminorm
                IF (.NOT. L2PROJ) THEN
                FieldNormV = FieldNormV + abs(zdivV(i))**2 * weight
                ENDIF
!
!               accumulate L2 norm
                do ivar=1,3
                  FieldNormV = FieldNormV + abs(zsolV(ivar,i))**2 * weight
                enddo
!
!        ...end loop over components
            enddo
!
!     ...end loop over integration points
         enddo
!
!===================================================================================
!  L2 ATTRIBUTE                                                                    |
!===================================================================================
      case(DISCON)
!
!     ...loop through integration points
         do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                         zdofH,zdofE,zdofV,zdofQ,nflag,                 &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
#if DEBUG_MODE
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
#endif
!
!           total weight
            weight=wa*rjac
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NRQVAR+ibeg+icomp
!
!               accumulate L2 norm
                FieldNormQ = FieldNormQ + abs(zsolQ(i))**2 * weight
!
!        ...end loop over components
            enddo
!
!     ...end loop over integration points
         enddo
!
         case default
         write(*,*)  'get_elem_Norm: UNKNOWN PHYSICAL ATTRIBUTE TYPE. stop.'
         stop
!
!  ...end select data type
      end select
!
!..end loop over physical attributes
   enddo
!
end subroutine get_elem_Norm
!
!
! routines below are currently not used
!
!!-------------------------------------------------------------------
!!Routine to evaluate the electric field norm of UW-Maxwell
!!by integrating H(curl) trace solution on a face of middle node.
!!the input variable fld_flag specifies signal or pump
!!fld_flag = 1 - signal
!!fld_flag = 0 - pump
!!zpoint - value of z at cross section face
!!Last modified : Feb 18
!!input: Mdle,facenumber,fld_flag
!!output: faceNorm,face_diff_norm
!!-------------------------------------------------------------------
!!face norms were used before 'get power' integration (probably not used anymore)
!!
!subroutine compute_element_faceNorm(mdle,facenumber,fld_flag,zpoint, faceNorm,face_diff_norm)
!!
!      use control
!      use data_structure3D
!      use environment      , only : L2PROJ
!      use physics
!      use parametersDPG
!      use commonParam
!      implicit none
!      integer, intent(in)                         :: mdle
!      integer, intent(in)                         :: fld_flag
!      integer, intent(in)                         :: facenumber
!      real(8),  intent(in)                        :: zpoint
!      real(8),  intent(out)                       :: faceNorm,face_diff_norm
!!
!!     element, face order, geometry dof
!      integer,dimension(19)           :: norder
!      real(8) ,dimension(3,MAXbrickH) :: xnod
!      integer,dimension(12)           :: nedge_orient
!      integer,dimension(6)            :: nface_orient
!! ...face order
!      integer, dimension(5)     :: norderf
!! .... number of vertices,edge,faces per element type
!      integer :: nrv, nre, nrf
!!.......declare edge/face type varibles
!      integer :: etype,ftype
!!
!!     geometry
!! ... variables for geometry
!      real(8), dimension(3)      :: xi,x,rn,x_new
!      real(8), dimension(3,2)    :: dxidt,dxdt,rt
!      real(8), dimension(3,3)    :: dxdxi,dxidx
!      real(8), dimension(2)      :: t
!      real(8)                    :: rjac,bjac
!!
!!     2D quadrature data
!      real(8), dimension(2,MAXNINT2ADD)  :: tloc
!      real(8), dimension(MAXNINT2ADD)    :: wtloc
!!
!!     approximate solution dof's
!      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
!      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
!      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
!      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!! ...H1 shape functions
!      integer                         :: nrdofH
!      real(8), dimension(MAXbrickH)   :: shapH
!      real(8), dimension(3,MAXbrickH) :: gradH
!!
!!     approximate solution
!      VTYPE, dimension(  MAXEQNH  )   ::  zsolH
!      VTYPE, dimension(  MAXEQNH,3)   :: zdsolH
!      VTYPE, dimension(3,MAXEQNE  )   ::  zsolE
!      VTYPE, dimension(3,MAXEQNE  )   :: zcurlE
!      VTYPE, dimension(3,MAXEQNV  )   ::  zsolV
!      VTYPE, dimension(  MAXEQNV  )   ::  zdivV
!      VTYPE, dimension(  MAXEQNQ  )   ::  zsolQ
!
!!     exact solution
!      VTYPE,dimension(  MAXEQNH    )  ::   ValH
!      VTYPE,dimension(  MAXEQNH,3  )  ::  DvalH
!      VTYPE,dimension(  MAXEQNH,3,3)  :: d2valH
!      VTYPE,dimension(3,MAXEQNE    )  ::   ValE
!      VTYPE,dimension(3,MAXEQNE,3  )  ::  DvalE
!      VTYPE,dimension(3,MAXEQNE,3,3)  :: d2valE
!      VTYPE,dimension(3,MAXEQNV    )  ::   ValV
!      VTYPE,dimension(3,MAXEQNV,3  )  ::  DvalV
!!
!!     exact solution (UNUSED)
!      VTYPE,dimension(3,MAXEQNV,3,3)  :: d2valV
!      VTYPE,dimension(  MAXEQNQ    )  ::   valQ
!      VTYPE,dimension(  MAXEQNQ,3  )  ::  dvalQ
!      VTYPE,dimension(  MAXEQNQ,3,3)  :: d2valQ
!!
!!     miscellanea
!      integer :: nint,icase,iattr,l,i,j
!      real(8) :: weight,wa
!      integer :: iel,nsign
!      integer :: nflag,iload
!!
!!---------------------------------------------------------------------------------------
!!
!      faceNorm = 0.d0
!      face_diff_norm = 0.0d0
!      nflag = 1
!! ...element type
!      etype = NODES(mdle)%ntype
!      nrv = NVERT(etype); nre = NEDGE(etype); nrf = NFACE(etype)
!      call find_order(mdle, norder)
!      call find_orient(mdle, nedge_orient,nface_orient)
!      call nodcor(mdle, xnod)
!      call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!! .....sign factor to determine the OUTWARD normal unit vector
!      nsign = nsign_param(etype,facenumber)
!!
!! .....face type
!      ftype = face_type(etype,facenumber)
!!
!! .....face order of approximation
!      call face_order(etype,facenumber,norder, norderf)
!!
!! .....set 2D quadrature
!      INTEGRATION = NORD_ADD
!      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
!      INTEGRATION = 0
!
!      do l=1,nint
!!
!! .......face coordinates
!        t(1:2) = tloc(1:2,l)
!!
!! .......face parametrization
!        call face_param(etype,facenumber,t, xi,dxidt)
!!
!! .......determine element H1 shape functions (for geometry)
!        call shape3DH(etype,xi,norder,nedge_orient,nface_orient, &
!                      nrdofH,shapH,gradH)
!!
!! .......geometry
!        call bgeom3D(mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
!                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
!        weight = bjac*wtloc(l)
!
!        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
!                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!        if(NEXACT.eq.1) then
!          call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE,  &
!                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!        endif
!!
!!       accumulate Norm of E - trace variable for signal (fld_flag=1) or pump (fld_flag=0)
!!       i.e., integrate (abs(E(1))**2+abs(E(2))**2+abs(E(3))**2) with:
!!                     E corresponding to signal if fld_flag = 1
!!                     E corresponding to pump if fld_flag = 0
!! .... first check for signal, i.e, if fld_flag = 1
!          if(fld_flag.eq.1) then
!            faceNorm = faceNorm+ &
!            sqrt(abs(zsolE(1,1))**2+abs(zsolE(2,1))**2+abs(zsolE(3,1))**2)*weight
!!   ..... if we have an exact
!            if(NEXACT.eq.1) then
!              face_diff_norm = face_diff_norm   &
!                           + sqrt(abs(valE(1,1)-zsolE(1,1))**2+ &
!                           + abs(valE(2,1)-zsolE(2,1))**2+abs(valE(3,1)-zsolE(3,1))**2)*weight
!            endif
!! .... next check for pump, i.e, if fld_flag = 0
!          else if(fld_flag.eq.0) then
!            faceNorm = faceNorm+ &
!            sqrt(abs(zsolE(1,3))**2+abs(zsolE(2,3))**2+abs(zsolE(3,3))**2)*weight
!!   ..... if we have an exact
!            if(NEXACT.eq.1) then
!              face_diff_norm = face_diff_norm   &
!                           + sqrt(abs(valE(1,3)-zsolE(1,3))**2+ &
!                           + abs(valE(2,3)-zsolE(2,3))**2+abs(valE(3,3)-zsolE(3,3))**2)*weight
!            endif
!          else
!            write(*,*) 'from compute_element_faceNorm: fld_flag must be 0 or 1'
!            stop
!          endif
!
!!       end loop over integration points
!      enddo
!
!end subroutine compute_element_faceNorm
!!
!!-------------------------------------------------------------------
!!Routine to evaluate the electric field norm of UW-Maxwell
!!by along the cross sections specified by the vector of z-values
!!zValues of length Num_zpts
!!the input variables:
!!input: zValues,Num_zpts,fld_flag
!!       zValues: vector of z-values along z-axis to evaluate
!!       Num_zpts: length of zValues vector
!!       fld_flag:
!!          fld_flag = 1 - signal
!!          fld_flag = 0 - pump
!!Last modified : Jan 18
!!output: Norm,diff_norm
!!-------------------------------------------------------------------
!!
!subroutine compute_FaceNorm(zValues,Num_zpts,fld_flag, norm,diff_norm)
!!
!      use data_structure3D
!      implicit none
!
!      integer, intent(in)                       :: Num_zpts
!      real(8), dimension(Num_zpts), intent(in)  :: zValues
!      integer, intent(in)                       :: fld_flag
!      real(8), dimension(Num_zpts), intent(out) :: norm,diff_norm
!      real(8) :: faceNorm, face_diff_norm, elemNorm
!! .... Mdle number
!      integer :: mdle
!! .... element, face order, geometry dof
!      real(8), dimension(3,MAXbrickH) :: xnod
!! .... number of vertices,edge,faces per element type
!      integer :: nrv, nre, nrf
!!
!!     miscellanea
!      integer :: iel,i
!!
!!---------------------------------------------------------------------------------------
!!
!! ... initialize outputs (vector of powers for all z-points)
!      norm = 0.d0
!      diff_norm = 0.d0
!
!! ... initialize running powers computed (elements per z-point)
!      faceNorm = 0.d0
!      face_diff_norm  = 0.d0
!
!      mdle=0
!!
!      do iel=1,NRELES
!        call nelcon(mdle, mdle)
!        call nodcor(mdle, xnod)
!        do i=1,Num_zpts
!          if((zValues(i).lt.xnod(3,8)).and.(zValues(i).gt.xnod(3,1))) then
!            call compute_element_faceNorm(mdle,2,fld_flag,zValues(i), faceNorm,face_diff_norm)
!            norm(i) = norm(i)+faceNorm
!            diff_norm(i)  = diff_norm(i) +face_diff_norm
!          endif
!        enddo
!      enddo
!end subroutine compute_FaceNorm
!!
!!-------------------------------------------------------------------
!!Driver routine to specify values of z-points in order to
!! evaluate the face norm of E-trace variables at those points corresponding
!!to linear problem
!!zValues of length Num_zpts
!!the input variable:  NONE
!!            output: NONE
!!Last modified : Feb 18
!!
!!-------------------------------------------------------------------
!subroutine get_faceNorm
!!
!  use commonParam
!  use laserParam
!  implicit none
!  integer                            :: fld_flag
!  real(8), allocatable, dimension(:) :: zValues
!  real(8), allocatable, dimension(:) :: face_norm
!  real(8), allocatable, dimension(:) :: diff_face_norm
!  integer :: zlength,numPts
!  real(8) :: a,b
!  integer :: i
!  a = 0.05d0 ! Initial value
!  b = 0.005d0
!  numPts = 200
!  if(ZL.gt.1.d0) then
!    zlength = numPts*int(ZL)
!  else
!    zlength = numPts
!  endif
!
!   allocate(zValues(zlength), face_norm(zlength), &
!            diff_face_norm(zlength))
!! ... get z-values
!   zValues = (/(((i-1)*b+a),i=1,zlength)/)
!! ... make z-values irrational numbers so that they are not on element
!! ... boundaries
!   zValues = zValues*PI*(7.d0/22.d0)
!   write(*,*) ' from get_faceNorm: zValues '
!   do i =1,zlength
!    write(*,*) zValues(i)
!   enddo
!! ... get face norm and difference of face norm for linear problem
!   fld_flag = 1
!   call compute_FaceNorm(zValues,zlength,fld_flag, face_norm,diff_face_norm)
!   write(*,*) ' from get_faceNorm: face_norm'
!   do i =1,zlength
!    write(*,*) face_norm(i)
!   enddo
!   write(*,*) ' from get_faceNorm: diff_face_norm'
!   do i =1,zlength
!    write(*,*) diff_face_norm(i)
!   enddo
!   deallocate(zValues,face_norm,diff_face_norm)
!end subroutine get_faceNorm
!!
!
