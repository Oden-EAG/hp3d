c-----------------------------------------------------------------------
c
c   routine name       - customize_geometry
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 09
c
c   purpose            - routine customizes geometry for a particular
c                        problem, this one is for the AF problem,
c                        it adds layers (including PML) and generates
c                        three membranes (realhead8mesh)
c
c-----------------------------------------------------------------------
c
      subroutine customize_geometry
c
      use control
      use GMP
      use AF
#include "syscom.blk"
c
c  ...layer thickness
      dimension thickness(6),no_dom(6),xp(3),nosrf(3)
c
c  ...split, bounding surfaces and surfaces to conform to
      dimension nsplit(2),ns_bound(10),ns_confm(10)


c      return
c
c  ...add skull, two layers of air and three layers of PML air
      nts=9
      nr_layers=6      !1-skull, 2,3-air,     4,5,6-PML
      thickness(1:6) = (/0.5d0, 2.d0, 2.d0, .75d0, .75d0, 1.5d0/)
c
c  ...add PML domain
      NRDOMAIN = NRDOMAIN+1
      write(*,*) 'customize_geometry: PML domain = ',NRDOMAIN
      no_dom(nr_layers-2:nr_layers)=NRDOMAIN
      no_dom(1) = 9; no_dom(2:3)=7
c
c  ...check PML parameters
      r1 = SURFACES(nts)%Rdata(4)
      do i=1,3
        r1 = r1 + thickness(i)
      enddo
      r2=r1
      do i=4,6
        r2 = r2 + thickness(i)
      enddo
      if ((abs(RPML_MIN-r1).gt.GEOM_TOL).or.
     .    (abs(RPML_MAX-r2).gt.GEOM_TOL)) then
        write(*,*) 'customize_geometry: INCONSISTENT PML DATA'
        write(*,*) 'RPML_MIN,RPML_MAX = ',RPML_MIN,RPML_MAX
        write(*,*) 'SPHERE RADIUSES = ',r1,r2
        stop 1
      endif
      call add_prism_layer(nts,nr_layers,thickness,no_dom)
c
c  ...create a temporary surface
      NRSURFS=NRSURFS+1
      SURFACES(NRSURFS)%Type = 'ThrPt'
      allocate(SURFACES(NRSURFS)%Rdata(9))
      SURFACES(NRSURFS)%Rdata(1:3) = (/0.d0,0.d0,0.d0/)
      SURFACES(NRSURFS)%Rdata(4:6) = (/1.d0,0.d0,0.d0/)
      nosrf(1:3) = (/1,13,NRSURFS/)
c
c  ...reset the domain number for prism within the ear canal
      do npri=1,NRPRISM
c
c  .....skip if not in the skull
        if (PRISMS(npri)%Domain.ne.9) cycle
c
c  .....check if within the cylinder bounding the air canal
        idec=1
        do iv=1,3
          np = PRISMS(npri)%VertNo(iv)
          xp(1:3) = POINTS(np)%Rdata(1:3)
          s = xp(2)**2 + xp(3)**2
          s = sqrt(s)
          if (s.gt..35d0+GEOM_TOL) idec=0
        enddo

        if (idec.eq.1) then
          write(*,*) 'customize_geometry: REDEFINING DOMAIN FOR PRISM',
     .                npri
          PRISMS(npri)%Domain = 7
c
c  .......reproject the points on the inner side of the ear canal
          do iv=1,3
            np = PRISMS(npri)%VertNo(iv)
            xp(1:3) = POINTS(np)%Rdata(1:3)
            s = xp(2)**2 + xp(3)**2
            s = sqrt(s)
            if (abs(s-.35d0).lt.GEOM_TOL) then
              SURFACES(NRSURFS)%Rdata(7:9) = xp(1:3)
              np = PRISMS(npri)%VertNo(iv+3)
              xp(1:3) = POINTS(np)%Rdata(1:3)
              call mnewt(1,nosrf,void1,void2,xp,void3,
     .                   POINTS(np)%Rdata(1:3))
            endif
          enddo
        endif
      enddo
c
c  ...delete the temporary surface
      deallocate(SURFACES(NRSURFS)%Rdata)
      SURFACES(NRSURFS)%Type = 'void'
      NRSURFS=NRSURFS-1
c
c  ...attach rectangles within the ear canal to the cylindrical surface
      do nr=1,NRRECTA
        idec=1
        do iv=1,4
          np = RECTANGLES(nr)%VertNo(iv)
          xp(1:3) = POINTS(np)%Rdata(1:3)
          s = xp(2)**2 + xp(3)**2
          s = sqrt(s)
          if (abs(s-.35d0).gt.GEOM_TOL) idec=0
        enddo
        if (idec.eq.1) then
          write(*,*)'customize geometry: ATTACHING RECTANGLE ', nr,
     .              ' TO THE EAR CANAL'
          RECTANGLES(nr)%Type='PTIRec'
          allocate(RECTANGLES(nr)%Idata(1))
          RECTANGLES(nr)%Idata(1)=1
        endif
      enddo
c      return
c
c
c  ...generate the ear drum...
      nsplit=10
      nr_bound=1; ns_bound(1)=8
      nr_confm=1; ns_confm(1)=1
      dh1 = 0.05d0; dh2 = 0.05d0
      call split_surface(nsplit,nr_bound,ns_bound,
     .                          nr_confm,ns_confm,dh1,dh2)
      write(*,*) 'customize_geometry: EAR DRUM GENERATED'
c
c  ...generate the oval window...
      nsplit=2
      nr_bound=1; ns_bound(1)=3
      nr_confm=2; ns_confm(1)=4; ns_confm(2)=6
      dh1 = 0.d0; dh2 = 0.03d0
      call split_surface(nsplit,nr_bound,ns_bound,
     .                          nr_confm,ns_confm,dh1,dh2)
      write(*,*) 'customize_geometry: COCHLEA OVAL WINDOW GENERATED'
c
c  ...generate the cochlea membrane
      nsplit=6
      nr_bound=3; ns_bound(1)=7; ns_bound(2)=3; ns_bound(3)=12
      nr_confm=3; ns_confm(1)=4; ns_confm(2)=5; ns_confm(3)=2
      dh1 = 0.005d0; dh2 = 0.005d0
      call split_surface(nsplit,nr_bound,ns_bound,
     .                          nr_confm,ns_confm,dh1,dh2)
      write(*,*) 'customize_geometry: COCHLEA MEMBRANE GENERATED'
c
      end
