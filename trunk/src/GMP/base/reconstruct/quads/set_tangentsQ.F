c----------------------------------------------------------------------
c
c   routine name       - set_tangents
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 07
c
c   purpose            - routine defines tangent versors for curves
c                        connected to a point on a reconstructed
c                        surface
c
c   arguments
c     in:
c      Np              - a point
c      List_curves     - list of curves connected to the point
c      Nrcurves        - number of curves on the list
c      Idec            = 1 generate the tangent bundles
c                      = 0 collect only the data on the bundle
c     out:
c      Newlist_curves  - a new list of the curves organized in
c                        the couterclockwise fashion;
c                        for i-th curve on the list:
c      List_recta(1,i) = a rectangle on the surface containing
c                        the curve
c      List_recta(2,i) = local vertex number of Np for the rectangle
c      List_recta(3,i) = orientation of the rectangle wrt to the
c                        surface normal,
c                         = 0 if the orientations are consistent
c                         = 1 if the orientations are opposite
c      Ifirst, Isecond = 0 if no sharp edge has been encountered
c                        the edge curve numbers if an edge has been found
c      Rn1,Rn2         - reconstructed normal(s)
c
c     stored in data structure:
c                        unit versors for curves forming the bundle(s)
c
c----------------------------------------------------------------------
c
      subroutine set_tangents(Np,List_curves,Nrcurves,Idec,
     .                        Newlist_curves,List_recta,
     .                        Ifirst,Isecond,Rn1,Rn2)
c
      use control
      use GMP
c
#include "syscom.blk"
c
      dimension List_curves(Nrcurves),Newlist_curves(Nrcurves),
     .          List_recta(3,Nrcurves)
      dimension Rn1(3),Rn2(3)
c
c  ...work lists containing adjacent surface rectangles data
      dimension list(10),list_rloc(3,10),list_hexa(20)
c
c  ...work space
      dimension temp(3)
c
c  ...secant vectors components
      dimension xyz(3,10)
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
c-----------------------------------------------------------------------
c
      if (Np.eq.11) then
        iprint=0
      else
        iprint=0
      endif
 5    continue
      if (iprint.eq.1) then
        write(*,7046) Np
 7046   format('set_tangents: Np = ',i6)
        write(*,7047) List_curves(1:Nrcurves)
 7047   format('set_tangents: ORIGINAL CURVE LIST = ',20i5)
      endif
c
c  ...there must be at least three curves...
      if (Nrcurves.le.2) then
        write(*,*) 'set_tangents: Nrcurves = ',Nrcurves
        stop 1
      endif
      Ifirst=0; Isecond=0
c
c-----------------------------------------------------------------------
c
c  STEP 1: collect the patch data
c
c  ...determine surface rectangles adjacent to the point
      nrhex=0
      nr_recta = 0
      do ic=1,Nrcurves
        nc = List_curves(ic)
c
c  .....loop through connected rectangles
        do ir=1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(ir)), nr,lab)
          do kk=1,2
            call decode(iabs(RECTANGLES(nr)%BlockNo(kk)), nhex,labh)
            if (labh.eq.2) then
              call locate(nhex, list_hexa,nrhex, loc)
              if (loc.eq.0) then
                nrhex=nrhex+1
                list_hexa(nrhex) = nhex
              endif
            endif
          enddo
c
c  .......skip if already on the list
          call locate(nr,list,nr_recta, ii)
          if (ii.ne.0) cycle
c
c  .......skip if not on the surface
          if (RECTANGLES(nr)%Type.ne.'HermRec') cycle
c
c  .......look for two edge curves on the original list
          do ie1=1,4
            nc1 = iabs(RECTANGLES(nr)%EdgeNo(ie1))
            if (iprint.eq.1) then
              write(*,8045) nc,ir,nr,ie1,nc1
 8045         format('nc,ir,nr,ie1,nc1 = ',i5,i2,i5,i2,i5)
            endif
            if (nc1.eq.nc) cycle
            call locate(nc1,List_curves,Nrcurves, ii)
            if (ii.ne.0) then
c
c  ...........add the rectangle to the list
              nr_recta = nr_recta + 1
              list(nr_recta) = nr
              if (iprint.eq.1) then
                write(*,7038) nc,nr
 7038           format('set_tangents: ADDING nc,nr = ',2i5)
              endif
ccc              exit
            endif
          enddo
        enddo
      enddo
      if (nr_recta.ne.Nrcurves) then
        write(*,7001) Nrcurves,nr_recta
 7001   format('set_tangents: Nrcurves,nr_recta = ',2i3)
        write(*,7049) list_hexa(1:nrhex)
 7049   format('              ADJACENT HEXAS = ',20i6)
        iprint=1
        call pause
        go to 5
        stop 1
      endif
c
c  ...reorganize the list of curves in the counterclokwise fashion
      nr = list(1); List_recta(1,1) = nr
c
c  ...find the adjacent hexahedron
      call decode(RECTANGLES(nr)%BlockNo(1), nh,lab)
      if (lab.ne.2.or.RECTANGLES(nr)%BlockNo(2).ne.0) then
        write(*,7002) Np,nc,nr,RECTANGLES(nr)%BlockNo(1:2)
 7002   format('set_tangents: Np,nc,nr,BLOCKS = ',3i5,2x,2i5)
        stop 1
      endif
c
c  ...find the local face number and the corresponding orientation
      do if=1,6
        call decode(HEXAS(nh)%FigNo(if), nr1,no1)
        if (nr1.eq.nr) then
          select case(if)
          case(2,3,4)
            if (no1.le.3) then
              norient=0
            elseif (no1.le.7) then
              norient=1
            else
              write(*,*) 'set_tangents: no1 = '; stop 1
            endif
          case(1,5,6)
            if (no1.le.3) then
              norient=1
            elseif (no1.le.7) then
              norient=0
            else
              write(*,*) 'set_tangents: no1 = '; stop 1
            endif
          end select
          exit
        endif
      enddo
      call locate(Np, RECTANGLES(nr)%VertNo(1:4),4, iv)
      iv3 = imod(iv+3,4)
      List_recta(2,1) = iv
      List_recta(3,1) = norient
      select case(norient)
c
c  .....orientation of the rectangle consistent with the normal
        case(0)
          Newlist_curves(1) = iabs(RECTANGLES(nr)%EdgeNo(iv))
          Newlist_curves(2) = iabs(RECTANGLES(nr)%EdgeNo(iv3))
c
c  .....orientation of the rectangle opposite to the normal
        case(1)
          Newlist_curves(1) = iabs(RECTANGLES(nr)%EdgeNo(iv3))
          Newlist_curves(2) = iabs(RECTANGLES(nr)%EdgeNo(iv))
      end select
c
c  ...complete the list
      do ic=2,Nrcurves
        do ir=1,Nrcurves
          nr = list(ir)
          call locate(Newlist_curves(ic),
     .                iabs(RECTANGLES(nr)%EdgeNo(1:4)),4, ic1)
          if ((ic1.ne.0).and.(nr.ne.List_recta(1,ic-1))) then
c
c  .........build up the list
            List_recta(1,ic) = nr
            call locate(Np, RECTANGLES(nr)%VertNo(1:4),4, iv)
            iv3 = imod(iv+3,4)
            List_recta(2,ic) = iv
            if (ic1.eq.iv) then
              List_recta(3,ic) = 0
              if (ic.lt.Nrcurves) then
                Newlist_curves(ic+1) = iabs(RECTANGLES(nr)%EdgeNo(iv3))
              endif
            elseif (ic1.eq.iv3) then
              List_recta(3,ic) = 1
              if (ic.lt.Nrcurves) then
                Newlist_curves(ic+1) = iabs(RECTANGLES(nr)%EdgeNo(iv))
              endif
            else
              write(*,*) 'set_tangents: ic1,iv = ',ic1,iv; stop 1
            endif
          endif
        enddo
      enddo
      if (iprint.eq.1) then
        do ic=1,Nrcurves
          write(*,7003) Newlist_curves(ic),List_recta(1:3,ic)
 7003     format('set_tangents: curve,rectangle,vert_no,orient = ',
     .           2i5,2i2)
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...check for the case of a point on a sharp edge
      ierror=0
      ile_edge_curves=0
      do ic=1,Nrcurves
        nc = Newlist_curves(ic)
        call check_edge(nc, id)
        if (id.eq.2) ile_edge_curves = ile_edge_curves+1
      enddo
      select case(ile_edge_curves)
      case(0)
c
c  ...reconstruct the normal
      temp(1:3) = POINTS(Np)%Rdata(4:6)
      call reconstruct_normal(Np,Newlist_curves,Nrcurves, temp)
      POINTS(Np)%Rdata(4:6) = temp(1:3)
c
      if (Idec.eq.0) return
      call set_bundle(Np,Nrcurves,POINTS(Np)%Rdata(4:6),
     .                Newlist_curves,List_recta)
      case(2)
c
c  ...set the tangent vectors corresponding to the edge and two tangent
c     planes...
      if (iprint.eq.1) then
        write(*,8038) Np
 8038   format('set_tangents: SETTING EDGE PATCHES FOR Np = ',i5)
      endif
c
c  ...collect the secant vectors
      Ifirst=0
      do ic=1,Nrcurves
        nc = Newlist_curves(ic)
        do j=1,2
          np1 = CURVES(nc)%EndPoNo(j)
          if (np1.eq.Np) cycle
          xyz(1:3,ic) = POINTS(np1)%Rdata(1:3) - POINTS(Np)%Rdata(1:3)
        enddo
        call check_edge(nc, id)
        if (id.eq.2) then
          if (Ifirst.eq.0) then
            Ifirst=ic
          else
            Isecond=ic
          endif
        endif
      enddo
      if ((Isecond-Ifirst.eq.1).or.(Ifirst+Nrcurves-Isecond.eq.1)) then
        write(*,*) 'set_tangents: Ifirst,Isecond,Nrcurves = ',
     .              Ifirst,Isecond,Nrcurves
        do ic=1,Nrcurves
          nr = List_recta(1,ic)
          call decode(RECTANGLES(nr)%BlockNo(1), nh,lab)
          write(*,7004) Newlist_curves(ic),nr,nh
 7004     format('set_tangents: curve,rectangle,hexa = ',3i5)
        enddo
        if (Isecond-Ifirst.eq.1) then
          nr = List_recta(1,Ifirst)
          call decode(RECTANGLES(nr)%BlockNo(1), nh,lab)
          do if=1,6
            call decode(HEXAS(nh)%FigNo(if), nr1,nor)
            if (nr.eq.nr1) exit
          enddo
          write(*,7005) nh,if
 7005     format('set_tangents: SPLIT HEXA nh,if = ',i5,i2)
          stop 1
        endif
         if (Ifirst+Nrcurves-Isecond.eq.1) then
          nr = List_recta(1,Isecond)
          call decode(RECTANGLES(nr)%BlockNo(1), nh,lab)
          do if=1,6
            call decode(HEXAS(nh)%FigNo(if), nr1,nor)
            if (nr.eq.nr1) exit
          enddo
          write(*,7005) nh,if
          stop 1
        endif
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...set the tangent vectors for the two curves on the sharp edge
      call norm(xyz(1:3,Ifirst), s1)
      call norm(xyz(1:3,Isecond), s2)
      temp(1:3) = xyz(1:3,Ifirst)/s1*s2 - xyz(1:3,Isecond)/s2*s1
      call normalize(temp)
      call scalar_product(xyz(1:3,Ifirst),temp, s)
      xyz(1:3,Ifirst) = s*temp(1:3)
      call normalize(xyz(1:3,Ifirst))
      call scalar_product(xyz(1:3,Isecond),temp, s)
      xyz(1:3,Isecond) = s*temp(1:3)
      call normalize(xyz(1:3,Isecond))
c
c  ...double check that the two edge vectors are parallel to each other
      call cross_product(xyz(1:3,Ifirst),xyz(1:3,Isecond),temp)
      call norm(temp, s)
      if (s.gt.GEOM_TOL) then
        write(*,*) 'set_tangents: INCONSISTENCY'
        stop 1
      endif
c
c-----------------------------------------------------------------------
c
c  ...set up the tangent vectors for the first plane
      call reconstruct_normal(Np,Newlist_curves(Ifirst),
     .                        Isecond+1-Ifirst, Rn1)
c
c  ...enforce orthogonality to the edge tangent
      call scalar_product(Rn1,xyz(1:3,Ifirst), s)
      Rn1(1:3) = Rn1(1:3) - s*xyz(1:3,Ifirst)
      call normalize(Rn1)
c
      next = Ifirst+1
      do while(next.ne.Isecond)
        call scalar_product(xyz(1:3,next),Rn1, s)
        xyz(1:3,next) = xyz(1:3,next) - s*Rn1(1:3)
        call normalize(xyz(1:3,next))
        next=next+1
      enddo
c
c  ...check the orientation of the normal vector
      call mixed_product(xyz(1:3,Ifirst),xyz(1:3,Ifirst+1),Rn1, s)
      if (s.lt.0.d0) then
        Rn1(1:3) = - Rn1(1:3)
      endif
      if (iprint.eq.1) write(*,7006) Rn1
 7006 format('set_tangents: FIRST NORMAL = ',3f8.3)
c
c-----------------------------------------------------------------------
c
c  ...set up the tangent vectors for the second plane
      j=0
      do i=Isecond,Ifirst+Nrcurves
        j=j+1
        ii = imod(i,Nrcurves)
        list(j)=Newlist_curves(ii)
        list_rloc(1:3,j) = List_recta(1:3,ii)
      enddo
      nr_curv_2nd_bundle = j
      call reconstruct_normal(Np,list,nr_curv_2nd_bundle, Rn2)
c
c  ...enforce orthogonality to the edge tangent
      call scalar_product(Rn2,xyz(1:3,Isecond), s)
      Rn2(1:3) = Rn2(1:3) - s*xyz(1:3,Isecond)
      call normalize(Rn2)
c
      next = imod(Isecond+1+Nrcurves, Nrcurves)
      do while(next.ne.Ifirst)
        call scalar_product(xyz(1:3,next),Rn2, s)
        xyz(1:3,next) = xyz(1:3,next) - s*Rn2(1:3)
        call normalize(xyz(1:3,next))
        next = imod(next+1,Nrcurves)
      enddo
c
c  ...check the orientation of the normal vector
      next = imod(Isecond+1+Nrcurves, Nrcurves)
      if (s.lt.0.d0) then
        Rn2(1:3) = - Rn2(1:3)
      endif
      if (iprint.eq.1) write(*,7007) Rn2
 7007 format('set_tangents: SECOND NORMAL = ',3f8.3)
      if (Idec.eq.0) return
c
c-----------------------------------------------------------------------
c
c  ...store the temporary unit vectors for the curves
      do ic=1,Nrcurves
        nc = Newlist_curves(ic)
        call locate(Np,CURVES(nc)%EndPoNo(1:2),2, i)
        select case(i)
        case(1)
          CURVES(nc)%Rdata(1:3) = xyz(1:3,ic)
        case(2)
          CURVES(nc)%Rdata(4:6) = -xyz(1:3,ic)
        end select
      enddo
c
c  ...set up the bundle
      call set_sharp_edge_bundle(Np,Nrcurves,
     .                           Newlist_curves,List_recta,
     .                           Ifirst,Isecond,Rn1,Rn2)
c
      case default
        write(*,*) 'set_tangents: ile_edge_curves = ',ile_edge_curves
        stop 1
      end select
c
c
      end subroutine set_tangents
