c----------------------------------------------------------------------
c
c   routine name       - set_tangents_Tsurf
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 10
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
c      List_trian(1,i) = a triangle on the surface containing
c                        the curve
c      List_trian(2,i) = local vertex number of Np for the triangle
c      List_trian(3,i) = orientation of the triangle wrt to the
c                        surface normal,
c                         = 0 if the orientations are consistent
c                         = 1 if the orientations are opposite
c      Ifirst, Isecond = 0 if no sharp edge has been encountered,
c                        otherwise edge curve numbers if a sharp edge
c                        has been found
c      Rn1,Rn2         - reconstructed normal(s) for sharp edge(s)
c
c     stored in data structure:
c                      - unit versors for curves forming the bundle(s)
c                      - normal for Np
c
c----------------------------------------------------------------------
c
      subroutine set_tangents(Np,List_curves,Nrcurves,Idec,
     .                        Newlist_curves,List_trian,
     .                        Ifirst,Isecond,Rn1,Rn2)
c
      use control
      use GMP
c
#include "syscom.blk"
c
      dimension List_curves(Nrcurves),Newlist_curves(Nrcurves),
     .          List_trian(3,Nrcurves)
      dimension Rn1(3),Rn2(3)
c
c  ...work lists containing adjacent surface triangles data
      dimension list(10),list_rloc(3,10),list_tets(20)
c
c  ...work space
      dimension temp(3),void(3,3),list_aux(Nrcurves),
     .                            list_aux1(3,Nrcurves)
c
c  ...secant vectors components
      dimension xyz(3,10)
      dimension vec(3,Nrcurves)
c
      integer, external :: imod
ccc      imod(j,mod) = j-(j-1)/mod*mod
c
c-----------------------------------------------------------------------
c
      ipass=0
      iprint=0

      if (Np.eq.3) then
        iprint=0
      else
        iprint=0
      endif

c
 5    continue
c  ...printing statement
      if (iprint.ge.1) then
        write(*,*)'----------------------------------------------------'
        write(*,7046) Np
 7046   format(' set_tangents: Np = ',i6)
        write(*,7047) List_curves(1:Nrcurves)
 7047   format(' set_tangents: ORIGINAL CURVE LIST = ',10(i8,2x)/30x
     .                                                ,10(i8,2x))
      endif
c
c  ...there must be at least three curves...
      if (Nrcurves.le.2) then
        write(*,*)'set_tangents: Nrcurves = ',Nrcurves
        stop
      endif
c
c  ...initialize sharp edges to 0
      Ifirst=0; Isecond=0
c
c-----------------------------------------------------------------------
c
c  STEP 1: collect the patch data
c
c  ...determine surface triangles adjacent to the point
      nrtet = 0
      nr_trian = 0
c  ...loop over connected curves
      do ic = 1,Nrcurves
        nc = List_curves(ic)
c  .....loop over triangles connected to curve
        do ir = 1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(ir)), nr,lab)
c  .......loop over blocks adjacent to triangle
          call decode(iabs(TRIANGLES(nr)%BlockNo(1)), ntet,labt)

          if (TRIANGLES(nr)%BlockNo(2).eq.0) then
c   .......if block is a tet
            select case (labt)
c   ........TET
            case(3)
              call locate(ntet,list_tets,nrtet, loc)
c  ...........if tet was not found add to list
              if (loc.eq.0) then
                nrtet = nrtet + 1
                list_tets(nrtet) = ntet
              endif
c  ...........we don't smooth for PRIS, HEXA, PYRA
            case default
              write(*,*)'set_tangents: invalid block type =',labt
              stop
            end select
          endif
c
c  .......cycle if triangle is already on the list or not on surface
          call locate(nr,list,nr_trian, ii)
          if (ii.ne.0) cycle
          if (TRIANGLES(nr)%Type.ne.'G1RecTri') cycle
c  .......loop over triangle edges
          do ie1 = 1,3
            nc1 = iabs(TRIANGLES(nr)%EdgeNo(ie1))
cc            if (iprint.ge.1) then
cc              write(*,8045) nc,ir,nr,ie1,nc1
cc 8045         format(' nc,ir,nr,ie1,nc1 = ',i5,i2,i5,i2,i5)
cc            endif
c  .........if edge coincides with current curve, then cycle
            if (nc1.eq.nc) cycle
            call locate(nc1,List_curves,Nrcurves, ii)
c  .........if edge was found on the curve list, add triangle to list
            if (ii.ne.0) then
              nr_trian = nr_trian + 1
              list(nr_trian) = nr
              if (iprint.ge.1) then
cc                write(*,7038) nc,nr
cc 7038           format(' set_tangents: ADDING nc,nr = ',2i5)
              endif
ccc              exit
            endif
c  .......end of loop over edges
          enddo
c  .....end of loop over triangles
        enddo
c  ...end of loop over curves
      enddo
c
c  ...check
      if (nr_trian.ne.Nrcurves) then
        write(*,7001) Nrcurves,nr_trian
 7001   format(' set_tangents: Nrcurves,nr_trian = ',2i3)
        write(*,7049) list_tets(1:nrtet)
 7049   format('              ADJACENT TETS = ',20i6)
        iprint=1
        call pause
        go to 5
        stop
      endif
c
      if (iprint.ge.1) then
        write(*,9003)list(1:nr_trian)
 9003   format(' set_tangents: triangles list = ',20(i5))
        write(*,9004)list_tets(1:nr_trian)
 9004   format(' set_tangents: tets list      = ',20(i5))
      endif
c
c-----------------------------------------------------------------------
c
c  STEP 2: reorganize the list of curves in the counterclokwise fashion
c
      if (iprint.ge.1) then
        write(*,*)'set_tangents: reorganizing curves for Np = ',Np
      endif
c  ...store 1st triangle in the patch
      nr = list(1);  List_trian(1,1) = nr
c  ...find tet adjacent to 1st triangle
      call decode(TRIANGLES(nr)%BlockNo(1), nh,lab)
c  ...check
      if ((lab.ne.3).or.(TRIANGLES(nr)%BlockNo(2).ne.0)) then
        write(*,7002) Np,nc,nr,TRIANGLES(nr)%BlockNo(1:2)
 7002   format('set_tangents: Np,nc,nr,BLOCKS = ',3i5,2x,2i5)
        stop
      endif
c
c
c  STEP 2a: find local face number and orientation of 1st triangle
c
c  ...loop over faces
      do ifig = 1,4
        call decode(TETRAS(nh)%FigNo(ifig), nr1,no1)
c  .....if face is found
        if (nr1.eq.nr) then
          if (iprint.ge.1) then
            write(*,*)'set_tangents: nr1,no1 = ',nr1,no1
          endif
c  .......select local face number
          select case(ifig)
c  .......clockwise faces (from element_data module!!!)
ccc          case(2,4)
          case(1,4)
            if (no1.le.2) then
              norient = 1
            elseif (no1.le.5) then
              norient = 0
            else
              write(*,*) 'set_tangents: no1 = '
              stop
            endif
c  .......counterclockwise faces (from element_data module!!!)
ccc          case(1,3)
          case(2,3)
            if (no1.le.2) then
              norient = 0
            elseif (no1.le.5) then
              norient = 1
            else
              write(*,*) 'set_tangents: no1 = '
              stop
            endif
          end select
          exit
        endif
c  ...end of loop over faces
      enddo
      if (iprint.ge.1) then
        write(*,*)'set_tangents: norient = ',norient
      endif
c
c
c  STEP 2b: reorganize curves in counterclockwise fashion starting
c           from 1st triangle
c
      call locate(Np,TRIANGLES(nr)%VertNo(1:3),3, iv)
      iv3 = imod(iv+2,3)
c  ...complete List_trian for 1st triangle
      List_trian(2,1) = iv
      List_trian(3,1) = norient
c  ...select orientation
      select case(norient)
c  ...orientation of the triangle consistent with the normal
      case(0)
        Newlist_curves(1) = iabs(TRIANGLES(nr)%EdgeNo(iv))
        Newlist_curves(2) = iabs(TRIANGLES(nr)%EdgeNo(iv3))
c  ...orientation of the triangle opposite to the normal
      case(1)
        Newlist_curves(1) = iabs(TRIANGLES(nr)%EdgeNo(iv3))
        Newlist_curves(2) = iabs(TRIANGLES(nr)%EdgeNo(iv))
      end select
c
c  ...loop over remaining curves
      do ic = 2,Nrcurves
c  .....loop over triangles
        do ir = 1,Nrcurves
          nr = list(ir)
          call locate(Newlist_curves(ic),
     .                iabs(TRIANGLES(nr)%EdgeNo(1:3)),3, ic1)
c  .......if curve is found
          if ((ic1.ne.0).and.(nr.ne.List_trian(1,ic-1))) then
c  .........build up the list
            List_trian(1,ic) = nr
            call locate(Np,TRIANGLES(nr)%VertNo(1:3),3, iv)
            iv3 = imod(iv+2,3)
            List_trian(2,ic) = iv
c  .........if consistent orientation
            if (ic1.eq.iv) then
              List_trian(3,ic) = 0
              if (ic.lt.Nrcurves) then
                Newlist_curves(ic+1) = iabs(TRIANGLES(nr)%EdgeNo(iv3))
              endif
c  .........else inconsistent orientation
            elseif (ic1.eq.iv3) then
              List_trian(3,ic) = 1
              if (ic.lt.Nrcurves) then
                Newlist_curves(ic+1) = iabs(TRIANGLES(nr)%EdgeNo(iv))
              endif
            else
              write(*,*)'set_tangents: ic1,iv = ',ic1,iv; stop 1
            endif
c  .......end if curve found
          endif
c  .....end of loop over triangles
        enddo
c  ...end of loop over remaining curves
      enddo
c
      if (iprint.ge.1) then
        write(*,9000)Newlist_curves(1:Nrcurves)
 9000   format(' set_tangents: curves reorganized: ',10(i8,2x)/30x
     .                                              ,10(i8,2x))
      endif
c  ...check
      if (iprint.ge.1) then
        do ic = 1,Nrcurves
          write(*,7003) Newlist_curves(ic),List_trian(1:3,ic)
 7003    format(' set_tangents: curve,triangle,vert_no,orient = ',
     .           2i8,2x,2i2)
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
cccc
cccc  STEP 3: reconstruct normal
cccc
cccc  ...collect secant vectors
ccc      do i = 1,Nrcurves
ccc        nc = Newlist_curves(i)
ccc        np1 = CURVES(nc)%EndPoNo(1)
ccc        np2 = CURVES(nc)%EndPoNo(2)
ccc        if (Np.eq.np1) then
ccc          vec(1:3,i) = POINTS(np2)%Rdata(1:3) - POINTS(np1)%Rdata(1:3)
ccc        else
ccc          vec(1:3,i) = POINTS(np1)%Rdata(1:3) - POINTS(np2)%Rdata(1:3)
ccc        endif
ccc      enddo
cccc  ...set normal to 0
ccc      POINTS(Np)%Rdata(4:6) = 0.d0
cccc  ...determine initial guess
ccc      do i = 1,Nrcurves
ccc        call cross_product(vec(1:3,i),vec(1:3,imod(i+1,Nrcurves)),
ccc     .                     temp(1:3))
ccc        call normalize(temp(1:3))
cccc  .....accumulate
ccc        POINTS(Np)%Rdata(4:6) = POINTS(Np)%Rdata(4:6) + temp(1:3)
ccc      enddo
ccc      call normalize(POINTS(Np)%Rdata(4:6))
cccc
cccc  ...compute normal
ccc      call reconstruct_curvature(Np,Newlist_curves,Nrcurves,
ccc     .                           POINTS(Np)%Rdata(4:6),void)
cccc
c
c
cccc
cccc  ...reconstruct normal
ccc      iprint=1
ccc      if (iprint.ge.1) then
ccc        write(*,9020)Np
ccc 9020 format(' set_tangents: reconstructing normal for Np = ',i3)
ccc      endif
ccc
cccc     NEED TO CHECK reconstruc_normal!!!
ccc
ccc      temp(1:3) = POINTS(Np)%Rdata(4:6)
ccc      write(*,9001)temp
ccc 9001 format(' set_tangents: before rec. = ',3(e12.5,2x))
cccccc      call reconstruct_normal_Tpatch(Np,Newlist_curves,Nrcurves, temp)
ccc      write(*,9002)temp
ccc 9002 format(' set_tangents: after rec.  = ',3(e12.5,2x))
ccc      POINTS(Np)%Rdata(4:6) = temp(1:3)
cccc
ccc      if (iprint.ge.1) then
ccc        write(*,*)'set_tangents: normal reconstructed.'
ccc      endif
cccc
cccc  ...generate tangent bundle if required
ccc      if (Idec.eq.0) return
ccc        call set_bundle(Np,Nrcurves,POINTS(Np)%Rdata(4:6),
ccc     .                  Newlist_curves,List_trian)
ccc
ccc
c
c
c
c-----------------------------------------------------------------------
c
c  STEP 3: set tagent bundle
c
c  ...check for the case of a point on a sharp edge
      ierror = 0
      ile_edge_curves = 0
c  ...loop over curves
      do ic = 1,Nrcurves
        nc = Newlist_curves(ic)
        call check_edge(nc, id)
        if (iprint.ge.1) then
          write(*,*)'set_tangents: nc,id = ',nc,id
        endif
c  .....if curve is on sharp edge, increase counter
        if (id.eq.2)  ile_edge_curves = ile_edge_curves + 1
c  ...end of loop over curves
      enddo
c
c  ...select # of curves in the patch that are sharp edges
      select case(ile_edge_curves)
c
c     N  O    S  H  A  R  P    E  D  G  E
      case(0)
c  STEP 3a: reconstruct normal
c
c  .....collect secant vectors
        do i = 1,Nrcurves
          nc = Newlist_curves(i)
          np1 = CURVES(nc)%EndPoNo(1)
          np2 = CURVES(nc)%EndPoNo(2)
          if (Np.eq.np1) then
            vec(1:3,i) = POINTS(np2)%Rdata(1:3) - POINTS(np1)%Rdata(1:3)
          else
            vec(1:3,i) = POINTS(np1)%Rdata(1:3) - POINTS(np2)%Rdata(1:3)
          endif
        enddo
c
c  .....determine initial guess for normal
        POINTS(Np)%Rdata(4:6) = 0.d0
        do i = 1,Nrcurves
          call cross_product(vec(1:3,i),vec(1:3,imod(i+1,Nrcurves)),
     .                       temp(1:3))
          call normalize(temp(1:3))
          POINTS(Np)%Rdata(4:6) = POINTS(Np)%Rdata(4:6) + temp(1:3)
        enddo
        call normalize(POINTS(Np)%Rdata(4:6))
c
c  .....compute actual normal and store it in data structure
        call reconstruct_curvature(Np,Newlist_curves,Nrcurves,
     .                             1,POINTS(Np)%Rdata(4:6),void)
c
c
c=======================================================================
c  ...return if only bundle info is needed                             |
      if (Idec.eq.0) return                                            !
c=======================================================================
c
c  STEP 3b: set tangent bundle
        call set_bundle(Np,Nrcurves,POINTS(Np)%Rdata(4:6),
     .                  Newlist_curves,List_trian)
c
c     2    S  H  A  R  P    E  D  G  E  S
      case(2)
c  .....update point type if needed
        if (POINTS(Np)%Type.ne.'SharpPt') then
          temp(1:3) = POINTS(Np)%Rdata(1:3)
          deallocate(POINTS(Np)%Rdata, STAT=is)
          if (is.ne.0) then
            write(*,*)'set_tangents: Rdata not deallocated for
     . Np = ',Np
            stop
          endif
          allocate(POINTS(Np)%Rdata(9), STAT=is)
          if (is.ne.0) then
            write(*,*)'set_tangents: Rdata not allocated for
     . Np = ',Np
            stop
          endif
          if (associated(POINTS(Np)%Idata)) then
            deallocate(POINTS(Np)%Idata, STAT=is)
            if (is.ne.0) then
              write(*,*)'set_tangents: Idata not deallocated for
     . Np = ',Np
              stop
            endif
          endif
          allocate(POINTS(Np)%Idata(Nrcurves), STAT=is)
          if (is.ne.0) then
            write(*,*)'set_tangents: Idata not allocated for
     . Np = ',Np
            stop
          endif
          POINTS(Np)%Type       = 'SharpPt'
          POINTS(Np)%Rdata(1:3) = temp(1:3)
        endif
c
c  STEP 3c: rearrange bundle curves so that Ifirst = 1
ccc        iprint=1
        if (iprint.ge.1) then
          write(*,8038)Np
 8038   format(' set_tangents: SETTING EDGE PATCHES FOR Np = ',i5)
          write(*,8039)Newlist_curves(1:Nrcurves)
 8039     format('                     CURVES: ',15(i8))
        endif
c
c
cccc  STEP 3d: collect the secant vectors
ccc
ccc        Ifirst = 0
cccc  .....loop over patch curves
ccc        do ic = 1,Nrcurves
ccc          nc = Newlist_curves(ic)
cccc  .......loop over curve endpoints
ccc          do j = 1,2
ccc            np1 = CURVES(nc)%EndPoNo(j)
ccc            if (np1.eq.Np) cycle
ccc            xyz(1:3,ic) = POINTS(np1)%Rdata(1:3) - POINTS(Np)%Rdata(1:3)
cccc  .......end of loop over endpoints
ccc          enddo
cccc  .......check if a sharp edge
ccc          call check_edge_Tsurf(nc, id)
cccc  .......store local patch # of sharp edge
ccc          if (id.eq.2) then
ccc            if (Ifirst.eq.0) then
ccc              Ifirst = ic
ccc            else
ccc              Isecond = ic
ccc            endif
ccc          endif
cccc  .....end of loop over patch curves
ccc        enddo
c
c
c  .....loop over patch curves
        do ic = 1,Nrcurves
          nc = Newlist_curves(ic)
c  .......check if a sharp edge
          call check_edge(nc, id)
c  .......store local patch # of sharp edge
          if (id.eq.2) then
            Ifirst = ic
            exit
          endif
c  .....end of loop over patch curves
        enddo
c  .....printing
        if (iprint.ge.1) then
          write(*,8040)Ifirst
 8040     format(' set_tangents: Ifirst = ',i2)
        endif
c  .....loop over patch curves
        do i = 1,Nrcurves
          j = imod(Ifirst+i-1,Nrcurves)
          list_aux(i)      = Newlist_curves(j)
          list_aux1(1:3,i) = List_trian(1:3,j)
          nc = list_aux(i)
          if (i.eq.1) cycle
c  .......skip 1st sharp edge
          call check_edge(nc, id)
          if (id.eq.2) Isecond = i
c  .....end of loop over patch curves
        enddo
c  .....reset Ifirst, Newlist_curves, List_trian
        Ifirst = 1
        Newlist_curves(1:Nrcurves) = list_aux(1:Nrcurves)
        List_trian(1:3,1:Nrcurves) = list_aux1(1:3,1:Nrcurves)
c  .....printing
        if (iprint.ge.1) then
          write(*,7555)Isecond,Newlist_curves(1:Nrcurves)
 7555     format(' set_tangents: Isecond = ',i2,';
     . Newlist_curves = ',15(i8))
          call pause
        endif
c
c  STEP 3d: collect secant vectors
c  .....loop over patch curves
        do i = 1,Nrcurves
          nc = Newlist_curves(i)
c  .......loop over curve endpoints
          do j = 1,2
            np1 = CURVES(nc)%EndPoNo(j)
            if (np1.eq.Np) cycle
            xyz(1:3,i) = POINTS(np1)%Rdata(1:3) - POINTS(Np)%Rdata(1:3)
c  .......end of loop over endpoints
          enddo
c  .....end of loop over patch curves
        enddo
c
c**********************************************************************
c       I F   S U B S E Q U E N T   S H A R P   E D G E S
        if (((Isecond-Ifirst)         .eq.1).or.
     .      ((Ifirst+Nrcurves-Isecond).eq.1))       then
         write(*,*) 'set_tangents: Ifirst,Isecond,Nrcurves = ',
     .              Ifirst,Isecond,Nrcurves
        do ic=1,Nrcurves
          nr = List_trian(1,ic)
          call decode(TRIANGLES(nr)%BlockNo(1), nh,lab)
          write(*,7004) Newlist_curves(ic),nr,nh
 7004     format('set_tangents: curve,triangle,tet = ',3(i8,2x))
        enddo
        if (Isecond-Ifirst.eq.1) then
          nr = List_trian(1,Ifirst)
          call decode(TRIANGLES(nr)%BlockNo(1), nh,lab)
          do if=1,4
            call decode(TETRAS(nh)%FigNo(if), nr1,nor)
            if (nr.eq.nr1) exit
          enddo
          write(*,7005) nh,if
 7005     format('set_tangents: SPLIT TETS nh,if = ',i5,i2)
          call print_GMP
          stop
        endif
         if (Ifirst+Nrcurves-Isecond.eq.1) then
          nr = List_trian(1,Isecond)
          call decode(TRIANGLES(nr)%BlockNo(1), nh,lab)
          do if=1,4
            call decode(TETRAS(nh)%FigNo(if), nr1,nor)
            if (nr.eq.nr1) exit
          enddo
          write(*,7005) nh,if
          stop
       endif
       call pause
      endif
c**********************************************************************
c
c----------------------------------------------------------------------
c
c  ...set the tangent vectors for the two curves on the sharp edges
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
        stop
      endif
c
c-----------------------------------------------------------------------
c
c  ...compute 1st guess for normal
      POINTS(Np)%Rdata(4:6) = 0.d0
      do i = 1,(Isecond-1)
        call cross_product(xyz(1:3,i),xyz(1:3,i+1), temp(1:3))
        call normalize(temp(1:3))
        POINTS(Np)%Rdata(4:6) = POINTS(Np)%Rdata(4:6) + temp(1:3)
      enddo
      call normalize(POINTS(Np)%Rdata(4:6))
      Rn1(1:3) = POINTS(Np)%Rdata(4:6)
c
c  ...compute 2nd guess for normal
      call reconstruct_curvature(Np,Newlist_curves(1:Isecond),
     .                           Isecond,1,Rn1(1:3), void)
c
c
c  ...compute final normal enforcing orthogonality to the edge tangent
      call scalar_product(Rn1,xyz(1:3,Ifirst), s)
      Rn1(1:3) = Rn1(1:3) - s*xyz(1:3,Ifirst)
      call normalize(Rn1)
c
      next = Ifirst + 1
      do while(next.ne.Isecond)
        call scalar_product(xyz(1:3,next),Rn1, s)
        xyz(1:3,next) = xyz(1:3,next) - s*Rn1(1:3)
        call normalize(xyz(1:3,next))
        next = next + 1
      enddo
c
c  ...check the orientation of the normal vector (UNNECESSARY!!!)
      call mixed_product(xyz(1:3,Ifirst),xyz(1:3,Ifirst+1),Rn1, s)
      if (s.lt.0.d0) then
        Rn1(1:3) = - Rn1(1:3)
      endif
c  ...store final normal in data structure
      POINTS(Np)%Rdata(4:6) = Rn1(1:3)
c  ...printing
      if (iprint.ge.1) write(*,7006) Rn1
 7006 format(' set_tangents: FIRST NORMAL  = ',3f8.3)
c
c----------------------------------------------------------------------
c
cccc  ...set up the tangent vectors for the second plane
ccc      j = 0
ccc      do i = Isecond,(Ifirst + Nrcurves)
ccc        j = j + 1
ccc        ii = imod(i,Nrcurves)
ccc        list(j)=Newlist_curves(ii)
ccc        list_rloc(1:3,j) = List_trian(1:3,ii)
ccc      enddo
ccc      nr_curv_2nd_bundle = j
ccc      call reconstruct_normal(Np,list,nr_curv_2nd_bundle, Rn2)
c
c
c
c  ...compute 1st guess for normal
      POINTS(Np)%Rdata(7:9) = 0.d0
      do i = Isecond,Nrcurves
        call cross_product(xyz(1:3,i),xyz(1:3,imod(i+1,Nrcurves)),
     .                     temp(1:3))
        call normalize(temp(1:3))
        POINTS(Np)%Rdata(7:9) = POINTS(Np)%Rdata(7:9) + temp(1:3)
      enddo
      call normalize(POINTS(Np)%Rdata(7:9))
      Rn2(1:3) = POINTS(Np)%Rdata(7:9)
c
c  ...compute 2nd guess for normal
      do i = 0,(Nrcurves-Isecond+1)
        next = imod(i+Isecond,Nrcurves)
        list_aux(i+1) = Newlist_curves(next)
      enddo
      call reconstruct_curvature(Np,
     .          list_aux(1:Nrcurves-Isecond+2),(Nrcurves-Isecond+2),
     .          1, Rn2(1:3),void)
c
c
c  ...compute final normal by enforcing orthogonality to the edge tangent
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
c
c
cccc  ...check the orientation of the normal vector
ccc      next = imod(Isecond+1+Nrcurves, Nrcurves)
ccc      if (s.lt.0.d0) then
ccc        Rn2(1:3) = - Rn2(1:3)
ccc        write(*,*)'here'
ccc        call pause
ccc      endif
c
c
c
c  ...store final normal in data structure
      POINTS(Np)%Rdata(7:9) = Rn2(1:3)
c  ...printing
      if (iprint.ge.1) write(*,7007) Rn2
 7007 format(' set_tangents: SECOND NORMAL = ',3f8.3)
c
c  ...store adjacent triangles
      do i = 1,Nrcurves
c  .....1st subpatch
        if (i.lt.Isecond) then
          POINTS(Np)%Idata(i) = List_trian(1,i)*10 + 1
c  .....2nd subpatch
        else
          POINTS(Np)%Idata(i) = List_trian(1,i)*10 + 2
        endif
      enddo
c
c=======================================================================
c  ...return if only bundle info is needed                             |
      if (Idec.eq.0) return                                            !
c=======================================================================
c
c-----------------------------------------------------------------------

c  ...store the temporary unit vectors for the curves
      do ic = 1,Nrcurves
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
     .                           Newlist_curves,List_trian,
     .                           Ifirst,Isecond,Rn1,Rn2)
c  ...inconsistent case!
      case default
        write(*,*)'set_tangents: ile_edge_curves = ',
     .             ile_edge_curves
        if (ipass.eq.1) stop
        write(*,*)'RERUNNIG ROUTINE...'
        iprint = 1
        ipass = 1
        goto 5
        stop
c
c  ...end select # of sharp edges
      end select
c
      iprint=0
      if (iprint.ge.1) then
        write(*,7777)Np,POINTS(Np)%Rdata(4:6)
 7777   format(' reconstruct_tangents: Np,Rn = ',i5,2x,3(e12.5,2x))
        call pause
      endif
c
c
      end subroutine set_tangents
