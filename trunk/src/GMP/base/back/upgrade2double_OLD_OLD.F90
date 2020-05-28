!==============================================================================
!  subroutine - upgrade2double
!==============================================================================
!  LATEST REVISION: Jul 09
!
!  PURPOSE: routine upgrades to double precision geometry using Newton
!           method iterations
!
!  REMARKS: new general routine!
!==============================================================================
!
subroutine upgrade2double
!
!------------------------------------------------------------------------------
! MODULES
  use kinds
  use control
  use GMP
!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
! VARIABLES
  real(DP)                 :: geom_tol_orig, prod
  real(DP), dimension(3)   :: x_0, void_1, grad_1, grad_2
  real(DP), dimension(4)   :: void_2
  real(DP), dimension(3,3) :: GRAD
  integer                  :: nt, iv, np, ns_1, ns_2, nc, np_s, i
  integer, dimension(3)    :: surfs
! FUNCTIONS
  integer                  :: mod3
!------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 1
!
! ..set GEOM_TOL to a fairly big value
    geom_tol_orig = GEOM_TOL;  GEOM_TOL = 0.1D-3
#if I_PRINT >= 1
    write(*,*)'upgrade2double: upgrading to double precision geometry...'
#endif
#if I_PRINT >= 1
    write(*,*)'upgrade2double: setting GEOM_TOL = ',GEOM_TOL
    call pause
#endif
! ..initialize void_1, void_2
    void_1 = 0.d0;  void_2 = 0.d0
! ..update points type
    call update_point_type
#if I_PRINT >= 2
    write(*,*)'upgrade2double: updating points type.'
    write(*,*)'************************************************************'
#endif
!
! ..loop through triangles
    do nt = 1, NRTRIAN
! ....cycle if not on a surface
      if (TRIANGLES(nt)%Type .ne. 'PTITri') cycle
#if I_PRINT >= 2
      write(*,*)'upgrade2double: nt, type =',nt,TRIANGLES(nt)%Type
#endif
! ....loop through vertices
      do iv = 1, 3
! ......get point number
        np = TRIANGLES(nt)%VertNo(iv)
#if I_PRINT >= 2
        write(*,*)'upgrade2double: np =',np
#endif
! ......cycle if point has already been visited
        if (POINTS(np)%Type .eq. 'Regular') cycle
! ......get intial guess for Newton-Rapson method
        x_0(1:3) = POINTS(np)%Rdata(1:3)
! ......select number of surfaces
        select case (POINTS(np)%Idata(1))
! ........point lies on ONE surface
          case (1)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 1_surf point.'
#endif
! ..........create two extra planes
            NRSURFS = NRSURFS + 1;  ns_1 = NRSURFS
            NRSURFS = NRSURFS + 1;  ns_2 = NRSURFS
! ..........get curve number
            nc = abs(TRIANGLES(nt)%EdgeNo(iv))
#if I_PRINT >= 2
  write(*,*)'upgrade2double: nc_1 =',nc
#endif
! ..........get curve start point
            np_s = CURVES(nc)%EndPoNo(1)
! ..........compute gradient
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,grad_1)
            else
              call curve(nc,1.d0, void_1,grad_1)
            endif
! ..........normalize gradient
            call normalize(grad_1)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: grad_1 =',grad_1
#endif
! ..........get curve
            nc = abs(TRIANGLES(nt)%EdgeNo(mod3(iv + 2)))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: nc_2 =',nc
#endif
! ..........compute gradient
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,grad_2)
            else
              call curve(nc,1.d0, void_1,grad_2)
            endif
! ..........normalize gradient
            call normalize(grad_2)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: grad_2 =',grad_2
#endif
! ..........create first plane
            SURFACES(ns_1)%Type = 'VecPt'
            allocate (SURFACES(ns_1)%Rdata(6))
            SURFACES(ns_1)%Rdata(1:3) = POINTS(np)%Rdata(1:3)
            SURFACES(ns_1)%Rdata(4:6) = grad_1
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 1st plane defined.'
#endif
! ..........create second plane
            SURFACES(ns_2)%Type = 'VecPt'
            allocate (SURFACES(ns_2)%Rdata(6))
            SURFACES(ns_2)%Rdata(1:3) = POINTS(np)%Rdata(1:3)
            SURFACES(ns_2)%Rdata(4:6) = grad_2
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 2nd plane defined.'
#endif
! ..........create array of surface numbers
            surfs(1) = POINTS(np)%Idata(2);  surfs(2) = ns_1;  surfs(3) = ns_2
! ..........apply Newton-Rapson method
            void_1 = 0.d0
            call mnewt(1,surfs,void_1,void_2,x_0,void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: Newton-Rapson applied.'
#endif
! ..........delete extra surfaces
            SURFACES(ns_1)%Type = 'Void';     SURFACES(ns_2)%Type = 'Void'
            deallocate(SURFACES(ns_1)%Rdata); deallocate(SURFACES(ns_2)%Rdata)
            NRSURFS = NRSURFS - 2
! ........poiint lies on TWO surfaces
          case (2)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 2_surf point.'
#endif
! ..........create one extra plane
10          NRSURFS = NRSURFS + 1;  ns_1 = NRSURFS
! ..........get curve number
            nc = abs(TRIANGLES(nt)%EdgeNo(iv))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: nc_1 = ',nc
#endif
! ..........get curve start point
            np_s = CURVES(nc)%EndPoNo(1)
! ..........compute gradient
            if (np .eq. np_s) then
              call curve(nc,0.d0, void_1,grad_1)
            else
              call curve(nc,1.d0, void_1,grad_1)
            endif
! ..........normalize gradient
            call normalize(grad_1)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: grad_1 = ',grad_1
#endif
! ..........create plane
            SURFACES(ns_1)%Type = 'VecPt'
            allocate (SURFACES(ns_1)%Rdata(6))
            SURFACES(ns_1)%Rdata(1:3) = POINTS(np)%Rdata(1:3)
            SURFACES(ns_1)%Rdata(4:6) = grad_1
#if I_PRINT >= 2
            write(*,*)'upgrade2double: extra plane defined.'
#endif
! ..........create array of surface numbers
            surfs(1:2) = POINTS(np)%Idata(2:3);  surfs(3) = ns_1
! ..........apply Newton-Rapson method
            void_1 = 0.d0
            call mnewt(1,surfs,void_1,void_2,x_0,void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: Newton method applied.'
#endif
! ..........delete extra surface
            SURFACES(ns_1)%Type = 'Void'
            deallocate(SURFACES(ns_1)%Rdata)
            NRSURFS = NRSURFS - 1
! ........point lies on THREE surfaces
          case (3)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: 3_surf point.'
#endif
! ..........compute mixed product of gradients
            do i = 1, 3
              call surf(POINTS(np)%Idata(i + 1),POINTS(np)%Rdata(1:3), void_1,GRAD(1:3,i))
              call normalize(GRAD(1:3,i))
            enddo
            call mixed_product(GRAD(1:3,1),GRAD(1:3,2),GRAD(1:3,3), prod)
#if I_PRINT >= 2
            write(*,*)'upgrade2double: mixed_product = ',prod
#endif
! ..........if a surface is redundant, treat it as a 2_surf point
            if (abs(prod) .lt. GEOM_TOL) goto 10
! ..........apply Newton method
            call mnewt(1,POINTS(np)%Idata(2:4),void_1,void_2,x_0,void_2, POINTS(np)%Rdata(1:3))
#if I_PRINT >= 2
            write(*,*)'upgrade2double: Newton method applied.'
#endif
! ........point lies on FOUR OR MORE surfaces
          case default
            write(*,*)'upgrade2double: case not supported yet!'
            stop
! ......end select number of sufaces
        end select
! ......reset point type
        POINTS(np)%Type = 'Regular'
        deallocate(POINTS(np)%Idata)
#if I_PRINT >= 2
        write(*,*)'upgrade2double: resetting point type.'
#endif
! ....end of loop through vertices
      enddo
#if I_PRINT >= 2
      call pause
#endif
! ..end of loop through triangles
    enddo
!
! ..reset GEOM_TOL to original value
    GEOM_TOL = geom_tol_orig
#if I_PRINT >= 1
    write(*,*)'upgrade2double: done!'
#endif
!
end subroutine upgrade2double
