!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   routine name    - pump_ode_solve
!
!----------------------------------------------------------------------
!
!   latest revision - Mar 2022
!
!   purpose         - Driver routine for computing pump power with
!                     pump ODE model, assumes pump is a plane wave
!                     over the entire core/cladding region
!
!   arguments
!                   - NumPts
!
!----------------------------------------------------------------------
!
subroutine pump_ode_solve
!
   use commonParam
   use laserParam
   use mpi_wrapper
!
   implicit none
!
!..transverse-averaged irradiance values and pump gain along the fiber
   real(8), allocatable :: pump_irr(:), sign_irr(:), gain_p(:)
!..auxiliary arrays
   real(8), allocatable :: zValues(:), dummy(:)
!
!..auxiliary variables
   integer :: numPts, i, j, max_it, ierr, fld
   real(8) :: a, eta, sum1, sum2, Is, Ip, dz, g0, l2norm, l2diff, eps
!
!..pump gain
!   = 1: computing pump gain from transverse-averaged signal irradiance
!   = 2: computing pump gain by integrating over transverse slices
   integer :: pump_gain
!
!----------------------------------------------------------------------
!
   pump_gain = 2
!
   if (.not. allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode: PUMP_VAL has not been initiated yet. stop.'
      stop
   endif
!
!..compute auxiliary variables:
!  - numPts is expected to be the number of elements in z-direction
   numPts = size(PUMP_VAL)
!  - dz is then the (average) element size in z-direction
   dz = ZL/numPts
   a = dz/2.d0
!
!..allocate arrays:
   allocate(zValues(numPts), pump_irr(numPts), gain_p(numPts))
   select case(pump_gain)
      case(1) ; allocate(sign_irr(numPts), dummy(numPts))
      case(2) ;
      case default ; write(*,*) 'pump_ode: invalid pump_gain param.'
   end select
!
!  - zValues (one z-coordinate per element)
   do i=1,numPts
      zValues(i) = (i-1)*dz+a
   enddo
!  - irrationalize z values to avoid points on element interfaces
   zValues = zValues*PI*(113.d0/355.d0)
!
!..compute auxiliary array (signal irradiance)
   if (pump_gain .eq. 1) then
!  ...i) compute signal power in fiber core
!        (note: only sign_irr is filled with valid entries in compute_power)
      fld = 1 ! signal field index
      call compute_power(zValues,numPts,fld, dummy,pump_irr,sign_irr,gain_p)
!     ii) collect signal power values (in fiber core) on each MPI proc
      call MPI_ALLREDUCE(MPI_IN_PLACE,sign_irr,numPts,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, ierr)
!     iii) compute transverse-averaged signal irradiance over the core area
!          from signal optical power in the core area
      sign_irr(:) = sign_irr(:) / (PI*R_CORE*R_CORE)
   endif
!
!..compute pump irradiance
   pump_irr(:) = PUMP_VAL(:) / (PI*R_CLAD*R_CLAD)
!
!..update pump irradiance
!  iterate a few times (fixed-point iteration) since the pump gain "gain_p(i)"
!  on the RHS depends on the pump irradiance "pump_irr(i)"
   max_it = 10; eps = 1.d-6
   do j=1,max_it
!
!  ...compute l2 norm of pump irradiance
      l2norm = NORM2(pump_irr)
!
!  ...compute pump gain
      select case(pump_gain)
      case(1)
         do i=1,numPts
            Is = sign_irr(i)
            Ip = pump_irr(i)
            sum1 = (SIGMA_S_ABS/OMEGA_SIGNAL)*Is+(SIGMA_P_ABS/OMEGA_PUMP)*Ip
            sum2 = ((SIGMA_S_ABS+SIGMA_S_EMS)/OMEGA_SIGNAL)*Is + &
                   ((SIGMA_P_ABS+SIGMA_P_EMS)/OMEGA_PUMP)*Ip
            eta = sum1/(TAU_0+sum2)
            gain_p(i) = -SIGMA_P_ABS + (SIGMA_P_ABS+SIGMA_P_EMS)*eta
         enddo
      case(2)
         fld = 0 ! pump field index
         call compute_gain(zValues,numPts,fld, gain_p)
         call MPI_ALLREDUCE(MPI_IN_PLACE,gain_p,numPts,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, ierr)
      case default ;
      end select
!
!  ...non-dimensional scaling factor for gain function
      g0 = ACTIVE_GAIN*L_0*SIGMA_0*NU_0*N_TOTAL
      gain_p(:) = g0 * gain_p(:) / (PI*R_CORE*R_CORE)
!
!  ...solve the pump ODE by explicit stepping in z-direction
!     (pos. z-direction: co-pumped; neg. z-direction: counter-pumped)
 2002 format(F10.3)
      l2diff = 0.d0
      if (COPUMP.eq.1) then
         do i=1,numPts-1
            !if (RANK.eq.ROOT) write(*,2002) gain_p(i)
            Ip = pump_irr(i+1)
            pump_irr(i+1) = pump_irr(i) + (R_CORE*R_CORE/(R_CLAD*R_CLAD)) * &
                                        dz * gain_p(i) * pump_irr(i)
            l2diff = l2diff + (Ip-pump_irr(i+1))**2.d0
         enddo
      elseif (COPUMP.eq.0) then
         do i=numPts,2,-1
            Ip = pump_irr(i-1)
            pump_irr(i-1) = pump_irr(i) + (R_CORE*R_CORE/(R_CLAD*R_CLAD)) * &
                                        dz * gain_p(i) * pump_irr(i)
            l2diff = l2diff + (Ip-pump_irr(i-1))**2.d0
         enddo
      else
         write(*,*) ' pump_ode_solve: COPUMP must be 1 or 0. stop.'
         stop
      endif
      l2diff = sqrt(l2diff)
!
!  ...update global pump power array based on irradiance update
!     (needed in each iteration since compute_gain uses PUMP_VAL array)
      PUMP_VAL(:) = pump_irr(:) * (PI*R_CLAD*R_CLAD)
!
!  ...check convergence of fixed-point iteration
      if (l2diff / l2norm < eps) then
         if (RANK.eq.ROOT) then
            write(*,*) 'pump_ode_solve: fixed-point iteration converged, #iter = ', j
         endif
         exit
      endif
      if ((j.eq.max_it) .and. (RANK.eq.ROOT)) then
         write(*,*) 'pump_ode_solve: fixed-point iteration not converged, max_it = ', max_it
         write(*,*) 'l2diff, l2norm = ', l2diff, l2norm
      endif
   enddo
!
!..deallocate auxiliary arrays
   deallocate(zValues, pump_irr, gain_p)
   select case(pump_gain)
      case(1) ; deallocate(sign_irr, dummy)
      case default ;
   end select
!
end subroutine pump_ode_solve
!
!----------------------------------------------------------------------
!
subroutine pump_ode_alloc(NumPts)
!
   use laserParam
!
   implicit none
!
   integer, intent(in) :: NumPts
!
   if (allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode: PUMP_VAL has already been allocated.'
   else
      allocate(PUMP_VAL(NumPts))
   endif
!
!..Initiate pump power values
   PUMP_VAL(1:NumPts) = PLANE_PUMP_POWER
!
end subroutine pump_ode_alloc
!
!----------------------------------------------------------------------
!
subroutine pump_ode_dealloc
!
   use laserParam
!
   implicit none
!
   if (.not. allocated(PUMP_VAL)) then
      write(*,*) 'pump_ode_dealloc: PUMP_VAL had not been allocated.'
   else
      deallocate(PUMP_VAL)
   endif
!
end subroutine pump_ode_dealloc
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!   routine name       - compute_gain
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2022
!
!   purpose            - Evaluates the active gain function g
!                        along the cross sections specified by
!                        the vector of zValues in the input
!
!   arguments
!        in:
!                      - ZValues     : sample points in z-direction
!                      - Num_zpts    : number of sample points
!                      - Fld         : 1 (signal) or 0 (pump)
!       out:
!                      - Gain        : Accumulated gain per face
!
!----------------------------------------------------------------------
!
subroutine compute_gain(ZValues,Num_zpts,Fld, Gain)
!
   use commonParam
   use laserParam
   use data_structure3D
   use environment, only : QUIET_MODE
   use mpi_wrapper
   use par_mesh   , only : DISTRIBUTED
!
   implicit none
!
   integer, intent(in)  :: Num_zpts
   real(8), intent(in)  :: ZValues(Num_zpts)
   integer, intent(in)  :: Fld
   real(8), intent(out) :: Gain(Num_zpts)
!
!..auxiliary variables
   real(8)    :: faceGain
!
!..mdle number
   integer :: mdle
!
!..element, face order, geometry dof
   real*8 :: xnod (3,8)
   real*8 :: maxz,minz
!
!..miscellanea
   integer :: iel, i, ndom
!
!..element type
   integer :: etype
!
!..face number over which power is computed
!  (in brick and prism, face 2 is face normal to xi3, at xi3=1)
   integer, parameter :: faceNum = 2
!
!..timer
   real(8) :: start_time,end_time
   integer :: ierr
!
!---------------------------------------------------------------------------------------
!
!..initialize outputs (vector of Gain for all z-points)
   Gain = 0.d0
!
!..initialize running gain values computed (elements per z-point)
   faceGain = 0.d0
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..iterate over elements
!
!$OMP PARALLEL DO                         &
!$OMP PRIVATE(mdle,etype,xnod,maxz,minz,  &
!$OMP         i,ndom,faceGain)            &
!$OMP REDUCTION(+:Gain)                   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      if (GEOM_NO.eq.5) then
!     ...skip elements outside the fiber core (gain region)
         call find_domain(mdle, ndom)
         if (ndom.ne.1 .and. ndom.ne.2) cycle
      endif
      call nodcor_vert(mdle, xnod)
      etype = NODES(Mdle)%ntype
      select case(etype)
         case(MDLB)
            maxz = maxval(xnod(3,1:8))
            minz = minval(xnod(3,1:8))
         case(MDLP)
            maxz = maxval(xnod(3,1:6))
            minz = minval(xnod(3,1:6))
         case default
            write(*,*) 'compute_gain: invalid etype param. stop.'
            stop
      end select
      do i=1,Num_zpts
         if((ZValues(i).le.maxz).and.(ZValues(i).gt.minz)) then
            call compute_faceGain(mdle,faceNum,Fld, faceGain)
            Gain(i) = Gain(i) + faceGain
         endif
      enddo
   enddo
!$OMP END PARALLEL DO
!
   90 continue
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,3010) end_time-start_time
 3010 format('  compute_gain : ',f12.5,'  seconds')
   endif
!
end subroutine compute_gain
!
!
!----------------------------------------------------------------------
!
!   routine name       - compute_face_gain
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2022
!
!   purpose            - Evaluates the gain function g
!                        by integrating H(curl) trace solution
!                        on a face of a middle node.
!
!   arguments
!        in:
!                      - Mdle       : middle element node
!                      - Facenumber : element face used for integration
!                      - Fld        : 1 (signal) or 0 (pump)
!       out:
!                      - FaceGain   : gain function accumulated on face
!
!----------------------------------------------------------------------
!
subroutine compute_faceGain(Mdle,Facenumber,Fld, FaceGain)
!
   use control
   use data_structure3D
   use physics
   use parametersDPG
   use commonParam
   use laserParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Fld
   integer, intent(in)  :: Facenumber
   real(8), intent(out) :: FaceGain
!
!..element, face order, geometry dof
   integer, dimension(19)          :: norder
   real(8), dimension(3,MAXbrickH) :: xnod
   integer, dimension(12)          :: nedge_orient
   integer, dimension(6)           :: nface_orient
!
!..face order
   integer :: norderf(5)
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!
!..declare edge/face type varibles
   integer :: etype,ftype
!
!..variables for geometry
   real(8), dimension(3)   :: xi,x,rn
   real(8), dimension(3,2) :: dxidt,dxdt
   real(8), dimension(3,3) :: dxdxi,dxidx
   real(8), dimension(2)   :: t
   real(8)                 :: rjac,bjac
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD) :: tloc
   real(8), dimension(MAXNINT2ADD)   :: wtloc
!
!..approximate solution dof's
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..H1 shape functions
   integer                         :: nrdofH
   real(8), dimension(MAXbrickH)   :: shapH
   real(8), dimension(3,MAXbrickH) :: gradH
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
!..for Poynting vector
   VTYPE :: EtimesH(3)
!
!..miscellanea
   integer :: nint,l,jz
   real(8) :: weight
   integer :: nsign
   integer :: nflag,numPts
!
!..
   real(8) :: sigma_abs, sigma_ems, Is, Ip
   real(8) :: eta, sum1, sum2, dz, minz, maxz, coord_z
!
!---------------------------------------------------------------------------------------
!
   FaceGain = 0.d0
!
   nflag = 1
!..element type
   etype = NODES(Mdle)%ntype
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
   call find_order(Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
   call nodcor(mdle, xnod)
   !..determine z-coordinate inside the element
   select case(etype)
      case(MDLB)
         maxz = maxval(xnod(3,1:8))
         minz = minval(xnod(3,1:8))
      case(MDLP)
         maxz = maxval(xnod(3,1:6))
         minz = minval(xnod(3,1:6))
      case default
         write(*,*) 'elem_maxwell: unexpected etype=',etype,'. stop.'
         stop
   end select
   coord_z = (minz + maxz) / 2.d0
   call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..face type
   ftype = face_type(etype,Facenumber)
!
!..face order of approximation
   call face_order(etype,Facenumber,norder, norderf)
!
!..set 2D quadrature
   INTEGRATION = NORD_ADD ! why ?
   call set_2D_int(ftype,norderf,nface_orient(Facenumber), nint,tloc,wtloc)
   INTEGRATION = 0
!
!..Fld=1: compute signal gain
   if(Fld.eq.1) then
      sigma_abs = SIGMA_S_ABS
      sigma_ems = SIGMA_S_EMS
   endif
!..Fld=0: compute pump gain
   if(Fld.eq.0) then
      sigma_abs = SIGMA_P_ABS
      sigma_ems = SIGMA_P_EMS
   endif
!
!..loop over integration points
   do l=1,nint
!
!  ...face coordinates
      t(1:2) = tloc(1:2,l)
!
!  ...face parametrization
      call face_param(etype,Facenumber,t, xi,dxidt)
!
!  ...determine element H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,gradH)
!
!  ...geometry
      call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                   x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
      weight = bjac*wtloc(l)
!
      call soleval(Mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                   zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
                   zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!  ...compute signal irradiance
      call zz_cross_product(zsolE(1:3,1),conjg(zsolE(1:3,2)), EtimesH)
      Is = sqrt(real(EtimesH(1))**2+real(EtimesH(2))**2+real(EtimesH(3))**2)
!
!  ...compute pump irradiance
      if (PLANE_PUMP .eq. 1) then
!     ...assume pump is a plane wave in fiber cladding (cladding-pumped)
         Ip = PLANE_PUMP_POWER / (PI*R_CLAD*R_CLAD) ! calculate non-dimensional irradiance
      elseif (PLANE_PUMP .eq. 2) then
         numPts = size(PUMP_VAL)
         dz = ZL / numPts
         jz = min(INT(coord_z/dz), numPts-1) + 1
         Ip = PUMP_VAL(jz) / (PI*R_CLAD*R_CLAD) ! calculate non-dimensional irradiance
      else
         call zz_cross_product(zsolE(1:3,3),conjg(zsolE(1:3,4)), EtimesH)
         Ip = sqrt(real(EtimesH(1))**2+real(EtimesH(2))**2+real(EtimesH(3))**2)
      endif
!
!     accumulate gain function for signal (Fld=1) or pump (Fld=0),
      sum1 = (SIGMA_S_ABS/OMEGA_SIGNAL)*Is+(SIGMA_P_ABS/OMEGA_PUMP)*Ip
      sum2 = ((SIGMA_S_ABS+SIGMA_S_EMS)/OMEGA_SIGNAL)*Is + &
             ((SIGMA_P_ABS+SIGMA_P_EMS)/OMEGA_PUMP)*Ip
      eta = sum1/(TAU_0+sum2)
      FaceGain = FaceGain + (-sigma_abs + (sigma_abs+sigma_ems)*eta)*weight
!..end loop over integration points
   enddo
!
end subroutine compute_faceGain
!
!
