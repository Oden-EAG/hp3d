c----------------------------------------------------------------------
c
c   routine name       - set_sharp_edge_bundle
c
c----------------------------------------------------------------------
c
c   latest revision    - Sep 07
c
c   purpose            - routine sets tangent and normal dircetions
c                        for curves connected to a point on a sharp
c                        edge
c
c   arguments
c     in:
c      Np              - a point
c      Nrcurves        - number of curves on the list
c      Newlist_curves  - list of the curves connected to the point
c                        organized in the couterclockwise fashion;
c                        for i-th curve on the list:
c      List_recta(1,i) = a rectangle on the surface containing
c                        the curve
c      List_recta(2,i) = local vertex number of Np for the rectangle
c      List_recta(3,i) = orientation of the rectangle wrt to the
c                        surface normal,
c                         = 0 if the orientations are consistent
c                         = 1 if the orientations are opposite
c      Ifirst, Isecond = the edge curve numbers
c      Rn1,Rn2         - normals for the two surfaces
c     out:
c                        unit tangent vectors for the curves and
c                        normal components of the second derivatives
c                        stored in data structure
c
c----------------------------------------------------------------------
c
      subroutine set_sharp_edge_bundle(Np,Nrcurves,
     .                                 Newlist_curves,List_recta,
     .                                 Ifirst,Isecond,Rn1,Rn2)
c
      use control
      use GMP
c
#include "syscom.blk"
c
      dimension Newlist_curves(Nrcurves),
     .          List_recta(3,Nrcurves)
      dimension Rn1(3),Rn2(3)
c
c  ...edge vector, dual basis to Rn1, Rn2 in the plane orthogonal
c     to edge vector
      dimension evect(3),dual1(3),dual2(3)
c
c  ...work space
      dimension temp(3),list(10),Rn1_save(3),Rn2_save(3)
c
c  ...second derivatives of the implicit reconstruction of surfaces
      dimension acoeff1(3,3),acoeff2(3,3)
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
c-----------------------------------------------------------------------
c
ccc      if (Np.eq.428) then
        iprint=0
ccc      else
ccc        iprint=0
ccc      endif
c
c  ...there must be at least three curves in each smooth part of the
c     bundle
      nrcurves1 = Isecond+1-Ifirst
      nrcurves2 = Ifirst+1+Nrcurves-Isecond
      if ((nrcurves1.le.2).or.(nrcurves2.le.2)) then
        write(*,7001) nrcurves1,nrcurves2
 7001   format('set_sharp_edge_bundle: nrcurves1,nrcurves2 = ',2i3)
        stop
      endif
      if (iprint.eq.1) then
        write(*,*) '----------------------------------------------'
        write(*,7002) Np,Rn1,Rn2
 7002   format(' set_sharp_edge_bundle: DEBUGGING FOR Np,Rn1,Rn2 = ',
     .          i4,2x,2(3(e12.5,2x),';'))
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 1: collect the edge vector and determine the basis dual to
c          Rn1,Rn2 in the plane orthogonal to the edge evector
c
      nc = Newlist_curves(Ifirst)
      call locate(Np,CURVES(nc)%EndPoNo(1:2),2, ii)
      select case(ii)
      case(1)
        evect(1:3) = CURVES(nc)%Rdata(1:3)
      case(2)
        evect(1:3) = -CURVES(nc)%Rdata(4:6)
      end select
      call normalize(evect)
c
c  ...check consistency, the edge tangents should be opposite for the
c     two curves on the sharp edge
      nc = Newlist_curves(Isecond)
      call locate(Np,CURVES(nc)%EndPoNo(1:2),2, ii)
      select case(ii)
      case(1)
        temp(1:3) = CURVES(nc)%Rdata(1:3)
      case(2)
        temp(1:3) = -CURVES(nc)%Rdata(4:6)
      end select
      call normalize(temp)
      do ivar = 1,3
        s = abs(evect(ivar)+temp(ivar))
        if (s.gt.GEOM_TOL) then
          write(*,7003) Np
 7003     format('set_sharp_edge_bundle: INCONSISTENT DATA FOR Np = ',
     .           i5)
          stop
        endif
      enddo
c
c  ...set the dual basis
      call cross_product(evect,Rn1, dual2)
      call scalar_product(Rn2,dual2, s)
      dual2(1:3) = dual2(1:3)/s
      call cross_product(evect,Rn2, dual1)
      call scalar_product(Rn1,dual1, s)
      dual1(1:3) = dual1(1:3)/s
c
c-----------------------------------------------------------------------
c
c  STEP 2: generate the Hessians and store them
c
      call reconstruct_curvature(Np,Newlist_curves(1:Isecond),
     .                            nrcurves1,0,Rn1, acoeff1)
      do i = 1,nrcurves2
        ic = imod(Isecond-1+i,Nrcurves)
        list(i) = Newlist_curves(ic)
      enddo
      call reconstruct_curvature(Np,list,nrcurves2,0,Rn2, acoeff2)
c
      do ic = 1,Nrcurves
        nc = Newlist_curves(ic)
c
c       1  s  t    S  U  B  P  A  T  C  H
        if ((ic.eq.Ifirst).or.(ic.eq.Isecond)) then
          s1 = 0.d0; s2 = 0.d0
          do i = 1,3
            do j = 1,3
              s1 = s1 + acoeff1(i,j)*evect(i)*evect(j)
              s2 = s2 + acoeff2(i,j)*evect(i)*evect(j)
            enddo
          enddo
c
c  .......store the normals premultiplied with the curvatures
          call locate(Np,CURVES(nc)%EndPoNo(1:2),2, i)
          select case(i)
          case(1)
            CURVES(nc)%Rdata(7:9) = -s1*dual1(1:3)-s2*dual2(1:3)
          case(2)
            CURVES(nc)%Rdata(10:12) = -s1*dual1(1:3)-s2*dual2(1:3)
          end select
          if (iprint.eq.1) then
            write(*,7004) ic,nc,-s1*dual1(1:3)-s2*dual2(1:3)
 7004       format(' set_sharp_edge_bundle: STORING ',i2,'TH NORMAL
     . FOR CURVE ',i5,2x,3e12.5)
          endif
c
c       2  n  d    S  U  B  P  A  T  C  H
        else
c  .......copy the unit vector
          call locate(Np,CURVES(nc)%EndPoNo(1:2),2, i)
          select case(i)
          case(1)
            temp(1:3) = CURVES(nc)%Rdata(1:3)
          case(2)
            temp(1:3) = -CURVES(nc)%Rdata(4:6)
          end select
c
          if ((ic.gt.Ifirst).and.(ic.lt.Isecond)) then
c
c  .........check orthogonality
            call scalar_product(temp,Rn1, s)
            if (abs(s).gt.GEOM_TOL) then
              write(*,7005) s
 7005         format(' set_sharp_edge_bundle: s = ',e12.5)
              stop
            endif
c
            s = 0.d0
            do i=1,3
              do j=1,3
                s = s + acoeff1(i,j)*temp(i)*temp(j)
              enddo
            enddo
c
c  .........store the normal premultiplied with the curvature
            call locate(Np,CURVES(nc)%EndPoNo(1:2),2, i)
            select case(i)
            case(1)
              CURVES(nc)%Rdata(7:9) = -s*Rn1(1:3)
            case(2)
              CURVES(nc)%Rdata(10:12) = -s*Rn1(1:3)
            end select
            if (iprint.eq.1) then
              write(*,7004) ic,nc,-s*Rn1(1:3)
            endif
          else
c
c  .........check orthogonality
            call scalar_product(temp,Rn2, s)
            if (abs(s).gt.GEOM_TOL) then
              write(*,7005) s
              stop
            endif
c
            s = 0.d0
            do i=1,3
            do j=1,3
              s = s + acoeff2(i,j)*temp(i)*temp(j)
            enddo
            enddo
c
c  .........store the normal premultiplied with the curvature
            call locate(Np,CURVES(nc)%EndPoNo(1:2),2, i)
            select case(i)
            case(1)
              CURVES(nc)%Rdata(7:9) = -s*Rn2(1:3)
            case(2)
              CURVES(nc)%Rdata(10:12) = -s*Rn2(1:3)
            end select
            if (iprint.eq.1) then
              write(*,7004) ic,nc,-s*Rn2(1:3)
            endif
          endif
        endif
      enddo
      if (iprint.eq.1) then
        write(*,*)'exiting set_sharp_edge_bundle'
        call pause
      endif
c
c
      end subroutine set_sharp_edge_bundle
