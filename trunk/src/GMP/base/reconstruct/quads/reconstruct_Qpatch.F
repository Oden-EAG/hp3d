c----------------------------------------------------------------------
c
c   routine name       - reconstruct_patch
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 07
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
c      List_recta(1,i) = a rectangle on the surface containing
c                        the curve
c      List_recta(2,i) = local vertex number of Np for the rectangle
c      List_recta(3,i) = orientation of the rectangle wrt to the
c                        surface normal,
c                         = 0 if the orientations are consistent
c                         = 1 if the orientations are opposite
c
c----------------------------------------------------------------------
c
      subroutine reconstruct_patch(Np,List_curves,Nrcurves,Ipatch,
     .                             Rnorm,List_recta)
c
      use GMP
      use control
c
#include "syscom.blk"
      common /crhsub/ iprint_rhsub
      common /ctri/ iprint_tri
c
      dimension List_curves(Nrcurves),List_recta(1:3,Nrcurves),
     .          Rnorm(3)
c
c  ...list of tangent vectors and their derivatives
      parameter(maxc=10)
      dimension tangent(1:3,maxc), dtangent(1:3,maxc)
c
c  ...hesjan of the reconstructed surface at the point
      dimension acoeff(3,3)
c
c  ...derivatives of normal vectors
      dimension dnorm1(1:3),dnorm2(1:3)
c
c  ...for orientation "io" and vertex number "iv", function
c     nvar1(iv,io) = local axis number * orientation for the first edge
c                    adjacent to the vertex
c     nvar2(iv,io) = local axis number * orientation for the second edge
c                    adjacent to the vertex
      dimension nvar1(4,0:1), nvar2(4,0:1)
c
c  ...orientation sign
      dimension rsign(0:1)
c
c  ...work space
      dimension temp(3)
c
      data nvar1/1, 2,-1,-2,  2,-1,-2, 1/
      data nvar2/2,-1,-2, 1,  1, 2,-1,-2/
      data rsign/1.d0,-1.d0/
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
      ione=1
c
c**********************************************************************
c
      iprint=3
      if (iprint.ge.1) then
        write(*,7001) Np,Ipatch
 7001   format('reconstruct_patch: Np,Ipatch = ',i6,i2)
      endif
c
c  ...reconstruct the curvature at the point
      call reconstruct_curvature(Np,List_curves,Nrcurves,Rnorm,
     .                           acoeff)
c
c  ...collect tangent vectors emanating from point and their derivatives
      do i=1,Nrcurves
        nc = List_curves(i)
        call collect_tangent(nc,Np,tangent(1:3,i),dtangent(1:3,i))
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
c  .....check consistency
        call scalar_product(Rnorm,tangent(1:3,i), s)
        if (abs(s).gt.GEOM_TOL) then
          write(*,7478) Np,i,nc,s
 7478     format('reconstruct_patch: Np,i,nc,Rnorm o tangent = ',
     .            i7,i2,i7,2x,e12.5)
          write(*,7499) Rnorm
 7499     format('                   Rnorm = ',3f8.3)
          write(*,7479) tangent(1:3,i)
 7479     format('                   tangent = ',3f8.3)
          stop 1
        endif
      enddo
c
      select case(Ipatch)
c
c  ...360 degrees patch
      case(0)
        nr_rectas = Nrcurves
c
c  ...less than 360 degrees (open) patch
      case(1)
        nr_rectas = Nrcurves-1
      end select
c
c  ...save the mixed derivatives for the rectangles
      do i=1,nr_rectas
        nc = List_curves(i)
ccc        i1=imod(i+1,nr_rectas)
        i1=imod(i+1,Nrcurves)
        nr = List_recta(1,i)
        iv = List_recta(2,i)
        io = List_recta(3,i)
        if (iprint.eq.3) then
          write(*,8047) nr,iv,io
 8047     format('reconstruct_patch: SAVING DATA FOR RECTANGLE ',
     .           'nr,iv,io = ',i4,i3,i2)
        endif
c
c  .....compute the normal component of the  mixed derivative using
c       the local curvature, eq (2.29)
        dd = 0.d0
        do k=1,3
        do l=1,3
          dd = dd - acoeff(k,l)*tangent(k,i)*tangent(l,i1)
        enddo
        enddo
        if (nc.eq.3214) then
          temp(1:3) = tangent(1:3,i)
          call normalize(temp)
          s = 0.d0
          do k=1,3
          do l=1,3
            s = s + acoeff(k,l)*temp(k)*temp(l)
          enddo
          enddo
          write(*,*)'************',s,'***************'
        endif
c
c  .....determine index for Rdata array
        jj = (iv-1)*3
c  .....determine sign factor
        nsign = (-1)**(iv-1)
c  .....store normal component of mixed derivative
        RECTANGLES(nr)%Rdata(jj+1:jj+3) = nsign*dd*Rnorm(1:3)
c
c  REMARK: SECOND DERIVATIVES ARE INDEPENDENT OF THE ORIENTATION
c          OF THE RECTANGLE...
c
c
c  .....compute the derivatives of the normals at the point
        do ivar=1,3
          dnorm1(ivar)=0.d0
          dnorm2(ivar)=0.d0
          do k=1,3
            dnorm1(ivar) = dnorm1(ivar) + acoeff(ivar,k)*tangent(k,i)
            dnorm2(ivar) = dnorm2(ivar) + acoeff(ivar,k)*tangent(k,i1)
            do l=1,3
              dnorm1(ivar) = dnorm1(ivar) - acoeff(l,k)
     .                      *Rnorm(ivar)*Rnorm(l)*tangent(k,i)
              dnorm2(ivar) = dnorm2(ivar) - acoeff(l,k)
     .                      *Rnorm(ivar)*Rnorm(l)*tangent(k,i1)
            enddo
          enddo
        enddo
c
c  .....verify
        call scalar_product(Rnorm,dnorm1, s1)
        call scalar_product(Rnorm,dnorm2, s2)
        s = max(abs(s1),abs(s2))
        if (s.gt.GEOM_TOL) then
          write(*,*) 'reconstruct_patch: s1,s2 = ',s1,s2
          stop 1
        endif
        call scalar_product(Rnorm,dtangent(1:3,i), s1)
        call scalar_product(dnorm1,tangent(1:3,i), s2)
        s = abs(s1+s2)
        if (s.gt.GEOM_TOL) then
          write(*,8048) 'i ',List_curves(i),io,iv,s1,s2
 8048     format('reconstruct_patch: i,nc,io,iv,s1,s2 = ',
     .            a2,i6,2i3,2x,2e12.5)
          write(*,8039) tangent(1:3,i),dtangent(1:3,i)
 8039     format('                   tangent,dtangent = ',
     .            2(3f8.3,2x))
          write(*,8040) Rnorm(1:3), dnorm1(1:3)
 8040     format('                   Rnorm, dnorm1 = ',
     .            2(3f8.3,2x))
          stop 1
        endif
        call scalar_product(Rnorm,dtangent(1:3,i1), s1)
        call scalar_product(dnorm2,tangent(1:3,i1), s2)
        s = abs(s1+s2)
        if (s.gt.GEOM_TOL) then
          write(*,8048) 'i1',List_curves(i1),io,iv,s1,s2
          stop 1
        endif
c
c  .....store the derivatives of the normal for the vertex
        dnorm1(1:3) = dnorm1(1:3)*sign(ione,nvar1(iv,io))
        dnorm2(1:3) = dnorm2(1:3)*sign(ione,nvar2(iv,io))
        jj=12+(iv-1)*3
        RECTANGLES(nr)%Rdata(jj+1:jj+3) = Rnorm(1:3)*rsign(io)
        jj=24+(abs(nvar1(iv,io))-1)*12+(iv-1)*3
        RECTANGLES(nr)%Rdata(jj+1:jj+3) =
     .    dnorm1(1:3)*rsign(io)
        jj=24+(abs(nvar2(iv,io))-1)*12+(iv-1)*3
        RECTANGLES(nr)%Rdata(jj+1:jj+3) =
     .    dnorm2(1:3)*rsign(io)
c
c  ...end of loop through rectangles
      enddo
c
c
      end
