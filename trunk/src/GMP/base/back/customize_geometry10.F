c-----------------------------------------------------------------------
c   routine name       - customize_geometry
c-----------------------------------------------------------------------
c   latest revision    - Jul 09
c
c   purpose            - routine customizes geometry for a particular
c                        problem, this one is for the AF problem,
c                        it adds layers (including PML) and generates
c                        three membranes (head_10)
c
c   REMARKS  geometry constructed by routine:
c                        1  ossicles
c                        2  small air
c                        3  upper cochlea
c                        4  lower cochlea
c                        5  upper shell containing cochlea
c                        6  lower shell containing cochlea
c                        7  air
c                        8  small skull
c                        9  big skull
c                       10  brain
c                       11  PML
c                       12  ear drum
c                       13  oval window
c                       14  basilar membrane
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
c
c  ...printing flag
#define I_PRINT
c
#ifdef I_PRINT
      write(*,*)'customize_geometry: head_10 model considered.'
#endif
c
c---------------------------------------------------------------------
c  STEP 1: layers specifications
c
c  ...number of truncating sphere, according to head_10
      nts = 9
c  ...number of layers to add: 1-skull; 2,3-air; 4,5,6-PML
      nr_layers = 6
c  ...layers thickness
      thickness(1:6) = (/0.5d0, 2.d0, 2.d0, .75d0, .75d0, 1.5d0/)
c
c---------------------------------------------------------------------
c  STEP 2: set layers domain number
c
      write(*,*)'this custom geo.'
c  ...add PML domain (11)
      NRDOMAIN = NRDOMAIN + 1
#ifdef I_PRINT
      write(*,*) 'customize_geometry: PML domain = ',NRDOMAIN
#endif
c  ...set 3 most external layers to PML domain
      no_dom(nr_layers - 2:nr_layers) = NRDOMAIN
c  ...set 1st layer to BIG SKULL domain (9)
      no_dom(1) = 9
c  ...set 2nd and 3rd layer to AIR domain (7)
      no_dom(2:3) = 7
c
c  ...check PML parameters
      r1 = SURFACES(nts)%Rdata(4)
      do i = 1, 3
        r1 = r1 + thickness(i)
      enddo
      r2 = r1
      do i = 4, 6
        r2 = r2 + thickness(i)
      enddo
      if ((abs(RPML_MIN-r1).gt.GEOM_TOL).or.
     .                           (abs(RPML_MAX-r2).gt.GEOM_TOL)) then
        write(*,*) 'customize_geometry: INCONSISTENT PML DATA'
        write(*,*) 'RPML_MIN,RPML_MAX = ',RPML_MIN,RPML_MAX
        write(*,*) 'SPHERE RADIUSES = ',r1,r2
        stop
      endif
c
c  ...add prism layers
      call add_prism_layer(nts,nr_layers,thickness,no_dom)
#ifdef I_PRINT
      write(*,*)'customize_geometry: prismatic layers added.'
      write(*,*)'customize_geometry: NRSURFS = ',NRSURFS
c     REMARK:  NRSURFS = 13 + 6 = 19
#endif
c
c----------------------------------------------------------------------
c  STEP 3: "drill" a hole in the skull in order to open ear canal,
c    i.e. reset the domain number to 7 (BIG AIR) for prisms on the
c    skull and within the ear canal. Also reproject points on the
c    brim of the ear opening to the ear canal inner cylinder.
c
c  ...create a temporary plane through THREE POINTS
      NRSURFS = NRSURFS + 1
#ifdef I_PRINT
      write(*,*)'customize_geometry: temporary plane added, NSURFS = ',
     .                                                          NRSURFS
c     REMARK:  NRSURFS = 19 + 1 = 20
#endif
      SURFACES(NRSURFS)%Type = 'ThrPt'
      allocate(SURFACES(NRSURFS)%Rdata(9))
c  ...1st POINT
      SURFACES(NRSURFS)%Rdata(1:3) = (/0.d0, 0.d0, 0.d0/)
c  ...2nd POINT
      SURFACES(NRSURFS)%Rdata(4:6) = (/1.d0, 0.d0, 0.d0/)
c  ...surfaces for reprojecting on inner cylinder of ear canal:
c                1  - inner cylinder of ear canal
c               14  - outer sphere of skull
c     NRSURFS = 20  - temporary plane
      nosrf(1:3) = (/1, 14, NRSURFS/)
c
c  ...loop through prisms
      do npri = 1, NRPRISM
c  .....skip if not on the SKULL
        if (PRISMS(npri)%Domain .ne. 9) cycle
c  .....check if within the inner cylinder of ear canal
        idec = 1
c  .....loop through vertices of prism bottom face
        do iv = 1, 3
          np = PRISMS(npri)%VertNo(iv)
          xp(1:3) = POINTS(np)%Rdata(1:3)
          s = xp(2)**2 + xp(3)**2
          s = sqrt(s)
          if (s .gt. (0.35d0 + GEOM_TOL)) idec = 0
        enddo
c  .....if within the ear canal
        if (idec .eq. 1) then
c  .......set domain to BIG AIR (7)
          PRISMS(npri)%Domain = 7
#ifdef I_PRINT
          write(*,*)'customize_geometry: redefining domain for prism = '
     .,npri
#endif
c  .......reproject the points on the inner side of the ear canal
          do iv = 1, 3
            np = PRISMS(npri)%VertNo(iv)
            xp(1:3) = POINTS(np)%Rdata(1:3)
            s = xp(2)**2 + xp(3)**2
            s = sqrt(s)
            if (abs(s - 0.35d0) .lt. GEOM_TOL) then
c  ...........3rd POINT of extra plane
              SURFACES(NRSURFS)%Rdata(7:9) = xp(1:3)
              np = PRISMS(npri)%VertNo(iv + 3)
c  ...........initial guess for Newton method
              xp(1:3) = POINTS(np)%Rdata(1:3)
c  ...........use Newton method
              call mnewt(1,nosrf,void1,void2,xp,void3,
     .                                       POINTS(np)%Rdata(1:3))
            endif
          enddo
        endif
c  ...end of loop through prisms
      enddo
#ifdef I_PRINT
      write(*,*)'customize_geometry: hole drilled in skull!'
#endif
c
c  ...delete temporary plane
      deallocate(SURFACES(NRSURFS)%Rdata)
      SURFACES(NRSURFS)%Type = 'void'
      NRSURFS = NRSURFS - 1
c
c  ...attach all rectangles within the ear canal to cylindrical
c     surface 1 in order to avoid duplication of surfaces
      do nr = 1, NRRECTA
        idec = 1
        do iv = 1, 4
          np = RECTANGLES(nr)%VertNo(iv)
          xp(1:3) = POINTS(np)%Rdata(1:3)
          s = xp(2)**2 + xp(3)**2
          s = sqrt(s)
          if (abs(s - 0.35d0) .gt. GEOM_TOL) idec = 0
        enddo
c  .....if triangle is on the surface
        if (idec .eq. 1) then
#ifdef I_PRINT
          write(*,*)'customize_geometry: attaching rectangle = ',nr
#endif
c  .......update rectangle type
          RECTANGLES(nr)%Type = 'PTIRec'
          allocate(RECTANGLES(nr)%Idata(1))
c  .......set surface number to inner cylinder of ear canal (1)
          RECTANGLES(nr)%Idata(1) = 1
        endif
      enddo
#ifdef I_PRINT
      write(*,*)'customize_geometry: rectangles attached to cylinder.'
#endif
c





#ifdef I_PRINT
      write(*,*)'customize_geometry: surfaces'
      do i = 1, NRSURFS
        write(*,*)'ns, type = ',i,SURFACES(i)%Type
      enddo
      call pause
#endif

c---------------------------------------------------------------------
c  STEP 4: add membranes
c
c  ...generate the ear drum
      nsplit = 11
      nr_bound = 1;  ns_bound(1) = 10
      nr_confm = 2;  ns_confm(1) = 1; ns_confm(2) = 8
      dh1 = 0.05d0;  dh2 = 0.05d0
cc      call split_surface(nsplit,nr_bound,ns_bound,
cc     .                          nr_confm,ns_confm,dh1,dh2)
#ifdef I_PRINT
      write(*,*)'customize_geometry: surfaces'
      do i = 1, NRSURFS
        write(*,*)'ns, type = ',i,SURFACES(i)%Type
      enddo
      call pause
#endif
c
#ifdef I_PRINT
      write(*,*) 'customize_geometry: EAR DRUM added.'
#endif
c  ...generate the oval window
      nsplit = 2
      nr_bound = 1;  ns_bound(1) = 3
      nr_confm = 3;  ns_confm(1) = 4; ns_confm(2) = 6; ns_confm(3) = 8
      dh1 = 0.d0; dh2 = 0.02d0
cc      call split_surface(nsplit,nr_bound,ns_bound,
cc     .                          nr_confm,ns_confm,dh1,dh2)
#ifdef I_PRINT
      write(*,*)'customize_geometry: surfaces'
      do i = 1, NRSURFS
        write(*,*)'ns, type = ',i,SURFACES(i)%Type
      enddo
      call pause
#endif
c
#ifdef I_PRINT
      write(*,*) 'customize_geometry: OVAL WINDOW added.'
#endif
c  ...generate the basilar membrane
      nsplit = 6
      nr_bound = 3;  ns_bound(1) = 7; ns_bound(2) = 3; ns_bound(3) = 13
      nr_confm = 3;  ns_confm(1) = 4; ns_confm(2) = 5; ns_confm(3) = 2
      dh1 = 0.005d0; dh2 = 0.005d0
cc      call split_surface(nsplit,nr_bound,ns_bound,
cc     .                          nr_confm,ns_confm,dh1,dh2)
#ifdef I_PRINT
      write(*,*) 'customize_geometry: BASILAR MEMBRANE added.'
      write(*,*) 'customize_geometry: done!'
#endif
c
c---------------------------------------------------------------------------
#ifdef I_PRINT
      write(*,*)'customize_geometry: surfaces'
      do i = 1, NRSURFS
        write(*,*)'ns, type = ',i,SURFACES(i)%Type
      enddo
      call pause
#endif
c
      end subroutine customize_geometry
