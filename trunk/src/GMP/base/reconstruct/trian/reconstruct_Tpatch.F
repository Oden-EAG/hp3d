c----------------------------------------------------------------------
c
c   routine name       - reconstruct_patch
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine solves a patch problem determining
c                        second derivatives and derivatives of the
c                        normal; this is a new version based on
c                        reconstructing Hesjans at the points
c
c   arguments
c     in:
c      Np              - point on a reconstructed surface
c      List_curves     - a list of the curves forming the patch
c                        organized in the couterclockwise fashion
c      Nrcurves        - number of curves in the patch
c      Ipatch          = 0  for 360 degrees patch
c                      = 1  for < 360 degress patch
c      Rnorm           - normal vector for the patch
c                        for i-th curve on the list:
c      List_trian(1,i) = a triangle on the surface containing
c                        the curve
c      List_trian(2,i) = local vertex number of Np for the triangle
c      List_trian(3,i) = orientation of the rectangle wrt to the
c                        surface normal,
c                         = 0 if the orientations are consistent
c                         = 1 if the orientations are opposite
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_patch(Np,List_curves,Nrcurves,Ipatch,
     .                             Rnorm,List_trian)
c
      use GMP
      use control
c
      implicit none
c
      integer :: Np,Nrcurves,Ipatch
      integer :: List_curves(Nrcurves),List_trian(1:3,Nrcurves)
      real(8) :: Rnorm(3)
c
c  ...list of tangent vectors and their derivatives
      integer, parameter :: maxc=10
      real(8) :: tangent(1:3,maxc),dtangent(1:3,maxc)
c
c  ...hesjan of the reconstructed surface at the point
      real(8) :: acoeff(3,3)
c
c  ...work space
      real(8) :: temp(3)
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
      integer :: i,i1,io,iv,k,l,nc,nc1,nr,nr_trian
      real(8) :: dd,s
c
      integer :: iprint
c************************************************************************
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,*)'----------------------------------------------------'
        write(*,7001) Np,Rnorm,Ipatch
 7001  format(' reconstruct_patch: Np = ',i6,'; Rnorm = ',3(e12.5,2x),'
     .; Ipatch = ',i1)
      endif
c
c------------------------------------------------------------------------
c
c  STEP 0: reconstruct the Hessian at the point
c
      call reconstruct_curvature(Np,List_curves,Nrcurves,0,
     .                           Rnorm,acoeff)
      if (iprint.eq.1) then
        do i = 1,3
          write(*,1000)i,acoeff(i,1:3)
 1000     format(' Hessian(',i1,',1:3) = ',3(e12.5,2x))
        enddo
      endif
c
c------------------------------------------------------------------------
c
c  STEP 1: for each curve in the patch determine velocity vector
c
c  ...loop over patch curves
      do i = 1,Nrcurves
        nc = List_curves(i)
c  .....collect velocity vector
        call collect_tangent(nc,Np,tangent(1:3,i),dtangent(1:3,i))
        if (iprint.eq.1) then
          write(*,8999)nc,tangent(1:3,i)
 8999     format(' nc = ',i4,' ;   vel  = ',3(e12.5,2x))
        endif
c
c
c  .....correct the derivative of the tangent
ccc        call scalar_product(dtangent(1:3,i),Rnorm, s)
ccc        dtangent(1:3,i) = dtangent(1:3,i) - s*Rnorm(1:3)
ccc        s = 0.d0
ccc        do k=1,3
ccc        do l=1,3
ccc          s = s - acoeff(k,l)*tangent(k,i)*tangent(l,i)
ccc        enddo
ccc        enddo
ccc        dtangent(1:3,i) = dtangent(1:3,i) + s*Rnorm(1:3)
cccc
cccc  .....store the corrected derivative of the tangent vector
ccc        do ip=1,2
ccc          if (CURVES(nc)%EndPoNo(ip).eq.Np) then
ccc            ii = 6+(ip-1)*3
ccc            CURVES(nc)%Rdata(ii+1:ii+3) = dtangent(1:3,i)
ccc          endif
ccc        enddo
c
c
c  .....check consistency
        call scalar_product(Rnorm,tangent(1:3,i), s)
        if (abs(s).gt.GEOM_TOL) then
          write(*,*)'reconstruct_patch: inconsistency b/w normal and
     . tangent'
          write(*,7478) Np,i,nc
 7478     format(' Np = ',i8,'; i = ',i3,'; nc = ',i8)
          write(*,9000)Rnorm,tangent(1:3,i),s
 9000     format(' normal = ',3(e12.5,2x),'; tangent = ',3(e12.5,2x),
     . 'prod = ',e12.5)
          stop
        endif
c  ...end of loop over patch curves
      enddo
c
c------------------------------------------------------------------------
c
c  STEP 2: for each triangle in the patch determine normal component
c     of mixed derivative and store it
c
c  ...determine # of triangles according to patch type
      select case(Ipatch)
c  ...360 degrees patch
      case(0)
        nr_trian = Nrcurves
c  ...less than 360 degrees (open) patch
      case(1)
        nr_trian = Nrcurves - 1
      end select
c
c  ...loop over patch triangles
      do i = 1,nr_trian
c  .....1st curve
        nc  = List_curves(i)
        i1  = imod(i+1,nr_trian)
c  .....2nd curve
        nc1 = List_curves(i1)
        nr  = List_trian(1,i)
        iv  = List_trian(2,i)
        io  = List_trian(3,i)
c
c  .....compute mixed derivative using eq. (2.29)
        dd = 0.d0
        do k = 1,3
          do l = 1,3
             dd = dd - acoeff(k,l)*tangent(k,i)*tangent(l,i1)
          enddo
        enddo
        temp = dd*Rnorm
        if (iprint.eq.1) then
          write(*,1001)temp
 1001     format(' normal component of mixed derivative: ',3(e12.5,2x))
        endif
c
c  .....store in data structure according to orientation
        select case(iv)
        case(1)
          TRIANGLES(nr)%Rdata(27:29) = temp
        case(2)
          TRIANGLES(nr)%Rdata(39:41) = temp
        case(3)
          TRIANGLES(nr)%Rdata(93:95) = temp
        case default
          write(*,*)'reconstruct_patch: unknown vertex!'
          stop
        endselect
c
c  ...end of loop over patch triangles
      enddo
c
c
      end subroutine reconstruct_patch
