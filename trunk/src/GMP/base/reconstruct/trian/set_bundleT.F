c----------------------------------------------------------------------
c
c   routine name       - set_bundle
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 10
c
c   purpose            - routine sets tangent and normal directions
c                        for curves connected to a point on a surface
c
c   arguments
c     in:
c      Np              - a point
c      Nrcurves        - number of curves on the list
c      Rn              - normal at the point
c      Newlist_curves  - list of the curves connected to the point
c                        organized in the couterclockwise fashion;
c                        for i-th curve on the list:
c      List_trian(1,i) = a triangle on the surface containing
c                        the curve
c      List_trian(2,i) = local vertex number of Np for the triangle
c      List_trian(3,i) = orientation of the triangle wrt to the
c                        surface normal,
c                         = 0 if the orientations are consistent
c                         = 1 if the orientations are opposite
c     out:
c                        for each curve in the bundle: unit (!) tangent
c                        vectors and curve principal normal multiplied
c                        by curve curvature stored in data structure
c
c----------------------------------------------------------------------
c
      subroutine set_bundle(Np,Nrcurves,Rn,
     .                      Newlist_curves,List_trian)
c
      use control
      use GMP
      use element_data
c
#include "syscom.blk"
c
      dimension Newlist_curves(Nrcurves),
     .          List_trian(3,Nrcurves)
      dimension Rn(3)
c
c  ...normal components for hexahedra edges "orthogonal" to the
c     plane and adjacent to the point
      dimension edge_comp(10)
c
c  ...transformation matrix to patch coordinates
      dimension transf(3,3)
c
c  ...tangents for a curve
      dimension tangent(3)
c
c  ...work space
      dimension temp(3),temp1(3)
c
c  ...patch coordinates
      dimension xy(2,10),r(10),theta(10),dtheta(10),theta_new(10)
c
c  ...second derivatives of the implicit reconstruction of surface
      dimension acoeff(3,3)
c
c  ...points for a curvature reconstruction
      dimension list_points(20)
c
      data pi /3.14159265358979312d0/
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
      twopi = 2.d0*pi
ccc      theta_min = 1.d0/9.d0*pi
      theta_min = 1.d0/10.d0*pi
      theta_max = pi - theta_min
c
c-----------------------------------------------------------------------
c
      iprint=0
c
 5    continue
c
c  ...there must be at least three curves...
      if (Nrcurves.le.2) then
        write(*,*)'set_bundle: Nrcurves = ',Nrcurves
        stop
      endif
      if (iprint.eq.1) then
        write(*,*) '--------------------------------------------------'
        write(*,7001) Np, Rn
 7001   format(' set_bundle: Np, Rn = ',i4,2x,3f8.3)
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 1: for each tetrahedron adjacent to the point, determine the
c          normal component of the unit vector for the edge connected
c          to the point "Np" that does NOT lie on the reconstructed
c          surface
c
c  ...loop over curves
      do ir = 1,Nrcurves
c  .....get corresponding triangle
        nr = List_trian(1,ir)
c  .....find attached tet
        call decode(TRIANGLES(nr)%BlockNo(1), nh,lab)
c  .....loop over the edges
        do ie = 1,6
        nc = iabs(TETRAS(nh)%EdgeNo(ie))
c  .......loop over endpoints
          do ive = 1,2
            iv = TETRA_EDGE_TO_VERT(ive,ie)
c  .........if endpoint coincides with Np
            if (TETRAS(nh)%VertNo(iv).eq.Np) then
              call locate(nc,Newlist_curves,Nrcurves, ii)
c  ...........if curve is not on list, it is the wanted one
              if (ii.eq.0) then
                ive1 = imod(ive+1,2)
                iv1 = TETRA_EDGE_TO_VERT(ive1,ie)
                np1 = TETRAS(nh)%VertNo(iv1)
                temp(1:3) = POINTS(Np)%Rdata(1:3)-POINTS(np1)%Rdata(1:3)
                call normalize(temp)
ccc                call scalar_product(temp,Rn(1:3), edge_comp(ir))
                call cross_product(temp,Rn(1:3), temp1)
                call norm(temp1, edge_comp(ir))
              endif
            endif
c  .......end of loop over endpoints
          enddo
c  .....end of loop over edges
        enddo
c  ...end of loop over curves
      enddo
c
c  ...check the orientation
      icount=0
      do i = 1,Nrcurves
        if (edge_comp(i).lt.0.d0) icount = icount + 1
      enddo
      if (icount.eq.Nrcurves) then
        Rn(1:3) = -Rn(1:3)
        edge_comp(1:Nrcurves) = -edge_comp(1:Nrcurves)
      elseif (icount.ne.0) then
        write(*,*)'INCONSISTENT NORMAL'
        ! call pause
      endif
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7003) edge_comp(1:Nrcurves)
 7003   format(' set_bundle: edge_comp = ',10f8.3)
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 2: determine the angles for the tangent bundle
c
c  ...compute the transformation matrix to the patch coordinates
      nc = Newlist_curves(1)
c  ...loop over curve endpoints
      do ivc = 1,2
        np1 = CURVES(nc)%EndPoNo(ivc)
        if (np1.ne.Np) then
c  .......determine the secant vector
          tangent(1:3) = POINTS(np1)%Rdata(1:3)-POINTS(Np)%Rdata(1:3)
c  .......project on the tangent plane
          call scalar_product(tangent,Rn(1:3), s)
          tangent(1:3) = tangent(1:3) - s*Rn(1:3)
        endif
c  ...end of loop over curve endpoints
      enddo
      call normalize(tangent)
      transf(1,1:3) = tangent(1:3)
      call cross_product(Rn(1:3),tangent, temp)
      transf(2,1:3) = temp(1:3)
      transf(3,1:3) = Rn(1:3)
c
c  ...compute the patch coordinates...
      do i=1,Nrcurves
        nc = Newlist_curves(i)
        do ivc=1,2
          np1 = CURVES(nc)%EndPoNo(ivc)
          if (np1.ne.np) then
c
c  .........determine the secant vector
            tangent(1:3) = POINTS(np1)%Rdata(1:3)-POINTS(Np)%Rdata(1:3)
c
c  .........project on the tangent plane
            call scalar_product(tangent,Rn(1:3), s)
            tangent(1:3) = tangent(1:3) - s*Rn(1:3)
          endif
        enddo
        do j=1,2
          xy(j,i) = 0.d0
          do k=1,3
            xy(j,i) = xy(j,i) + transf(j,k)*tangent(k)
          enddo
        enddo
        call cart_to_polar(xy(1:2,i), r(i),theta(i))
c
c  .....enforce the range from 0 to twopi
        if (theta(i).lt.0.d0) theta(i) = theta(i)+twopi
      enddo
c
c  ...consistency check
      inconsistent=0
      do i = 2,Nrcurves
        dt = theta(i)-theta(i-1)
        if ((dt.lt.theta_min).or.(dt.gt.theta_max)) inconsistent = 1
      enddo
      if (inconsistent.eq.1) then
ccc        iprint=1
        write(*,*)'set_bundle: BAD ANGLES FOR NP = ',Np
        write(*,7004) Np,POINTS(Np)%Rdata(1:3)
 7004   format(' POINT = ',i6,' ; ',3e12.5)
        call pause
      endif
      if (iprint.eq.1) then
        do ic = 1,Nrcurves
          write(*,7006) ic,xy(1:2,ic), r(ic),theta(ic)
 7006     format(' set_bundle: i,xy,r,theta = ',i1,2x,4f8.3)
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 3: determine new angles by "promoting" the constant angle
c          condition
c
c  ...take equal angles
      dtheta(1:Nrcurves) = twopi/Nrcurves
c
c  ...determine the global rotation
      theta_new(1) = 0.d0
      do i = 2,Nrcurves
        theta_new(i) = theta_new(i-1) + dtheta(i-1)
      enddo
      theta0 = 0.d0
      theta(Nrcurves+1) = twopi
      do i=1,Nrcurves
        theta0 = theta0 + theta(i) - theta_new(i)
      enddo
      theta0 = theta0/Nrcurves
      theta_new(1) = theta0
      do i=2,Nrcurves
        theta_new(i) = theta_new(i-1) + dtheta(i)
      enddo
      theta_new(Nrcurves+1) = twopi + theta_new(1)
c
c  ...interpolate between the new and the old values
c     for the angles
      alpha = 1.d0
      do i=1,Nrcurves
        dtheta_new = theta_new(i+1) - theta_new(i)
        dtheta_old = theta(i+1) - theta(i)
        dt = dtheta_new-dtheta_old
        if (dt.gt.0) then
          alpha = min(alpha,(dtheta_new - theta_min)/dt)
        elseif (dt.lt.0) then
          alpha = min(alpha,(dtheta_new - theta_max)/dt)
        endif
      enddo
      do i=1,Nrcurves
        theta_new(i) = (1.d0-alpha)*theta_new(i) + alpha*theta(i)
      enddo
c
      if (iprint.eq.1) then
        write(*,7008) theta(1:Nrcurves)
 7008   format(' set_bundle: OLD thetas = ',10f8.3)
        write(*,7009) theta_new(1:Nrcurves)
 7009   format(' set_bundle: NEW thetas = ',10f8.3)
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 4: generate the unit (!) tangent vectors and store them along
c     with the curve (!) principal normal rescaled by the curve (!)
c     curvature
c
c  ...collect neighboring points
      ip=0
      do ic = 1,Nrcurves
        nr = List_trian(1,ic)
        do iv = 1,3
          npv = TRIANGLES(nr)%VertNo(iv)
          if (npv.ne.Np) then
            call locate(npv,list_points,ip, ii)
            if (ii.eq.0) then
              ip = ip + 1
              list_points(ip)=npv
            endif
          endif
        enddo
      enddo
c
c  ...reconstruct Hessian of tangent surface
      call reconstruct_curvature(Np,Newlist_curves,Nrcurves,
     .                           0,Rn,acoeff)
c
c  ...loop over connected curves
      do ic = 1,Nrcurves
        nc = Newlist_curves(ic)
        xy(1,ic) = cos(theta_new(ic))
        xy(2,ic) = sin(theta_new(ic))
        do i=1,3
          temp(i)=0.d0
          do j=1,2
            temp(i) = temp(i) + transf(j,i)*xy(j,ic)
          enddo
        enddo
c
c  .....check orthogonality
        call scalar_product(temp,Rn(1:3), s)
        if (abs(s).gt.GEOM_TOL) then
          write(*,7021) s
 7021     format(' set_bundle: s = ',e12.5)
          stop
        endif
c
c  .....compute the curve (!) curvature using eq (2.27)
        s = 0.d0
        do i = 1,3
          do j = 1,3
            s = s + acoeff(i,j)*temp(i)*temp(j)
          enddo
        enddo
c  .....store (curve principal normal) * curvature
        call locate(Np,CURVES(nc)%EndPoNo(1:2),2, i)
        select case(i)
        case(1)
          CURVES(nc)%Rdata(1:3) = temp(1:3)
          CURVES(nc)%Rdata(7:9) = -s*Rn(1:3)
        case(2)
          CURVES(nc)%Rdata(4:6) = -temp(1:3)
          CURVES(nc)%Rdata(10:12) = -s*Rn(1:3)
        end select
c
cc        iprint = 1
c  .....printing
        if (iprint.eq.1) then
          write(*,7020) nc,temp(1:3)
 7020     format(' set_bundle: STORING TANGENT FOR
     . CURVE     ',i5,': ',3(e12.5,2x))
          write(*,7022) nc,-s*Rn(1:3)
 7022     format('   STORING NORMAL RESCALED BY
     . CURVATURE FOR CURVE ',i5,': ',3(e12.5,2x))
        endif
c
c  ...end of loop over conneceted curves
      enddo
c
c  ...printing
      if (iprint.eq.1) then
        write(*,*) 'set_bundle: EXITING'
      endif
c
c
      end subroutine set_bundle
