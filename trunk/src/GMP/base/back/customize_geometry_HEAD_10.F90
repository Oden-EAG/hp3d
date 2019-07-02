!-------------------------------------------------------------------------------
!   routine name       - customize_geometry
!-------------------------------------------------------------------------------
!   latest revision    - Jul 09
!
!   purpose            - routine customizes geometry for a particular
!                        problem, this one is for the AF problem,
!                        it adds layers (including PML) and generates 
!                        three membranes (head_10)
!
!
!   INPUT geometry from head_10:
!
!                        NRDOMAIN = 10
!                        1  ossicles
!                        2  air (inside middle ear cavity)
!                        3  upper cochlea
!                        4  lower cochlea
!                        5  upper shell
!                        6  lower shell
!                        7  air (inside ear canal)
!                        8  small skull
!                        9  big skull
!                       10  brain
!
!    OUTPUT geometry:
!
!                        NRDOMAIN = 14
!                        1  ossicles                              OK
!                        2  air (inside middle ear cavity)        OK
!                        3  upper cochlea                         OK
!                        4  lower cochlea                         OK
!                        5  upper shell containing cochlea        OK
!                        6  lower shell containing cochlea        OK
!                        7  air (outside middle ear cavity)       OK
!                        8  small skull                           OK
!                        9  big skull                             OK
!                       10  brain                                 OK
!                       11  PML                                   OK 
!                       12  ear drum                              OK
!                       13  oval window                           OK
!                       14  basilar membrane                      KO!!
!-------------------------------------------------------------------------------
subroutine customize_geometry_head
!-------------------------------------------------------------------------------
! MODULES
  use kinds
  use U2D
  use control
  use AF
!-------------------------------------------------------------------------------
! VARIABLES
! layers thickness
  real(DP), dimension(6) :: thickness
! layers domain number  
  integer, dimension(6)  :: no_dom
  real(DP), dimension(3) :: xp
! bonding and conforming surfaces for membrane generation
  integer, dimension(10) :: ns_bound,ns_confm
! split surface number  
  integer                :: nsplit
  integer                :: nts,nr_layers,i,npri,idec,iv,np,nr
  integer                :: status
  real(DP)               :: r1,r2,s
!-------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 1
!
#if I_PRINT >= 1
    write(*,*)'customize_geometry_NEW: customizing head_10...'
#endif      
#if I_PRINT >= 2
    write(*,*)'customize_geometry_NEW: original number of subdomains =',NRDOMAIN
!   REMARK: NRDOMAIN = 10    
#endif
!-------------------------------------------------------------------------------
! STEP 1: set layers specifications
!      
! ..number of truncating sphere, according to head_10      
    nts = 9
! ..number of layers to add: 1-skull; 2,3-air; 4,5,6-PML      
    nr_layers = 6
! ..layers thickness      
    thickness(1:6) = (/0.5d0, 2.d0, 2.d0, .75d0, .75d0, 1.5d0/)
!---------------------------------------------------------------------
!  STEP 2: set layers domain number
!
! ..add PML domain
    NRDOMAIN = NRDOMAIN + 1
#if I_PRINT  >= 2    
    write(*,*) 'customize_geometry_NEW: PML domain = ',NRDOMAIN
!   REMARK: NRDOMAIN = 11    
#endif
! ..set 3 most external layers to PML domain      
    no_dom(nr_layers - 2:nr_layers) = NRDOMAIN
! ..set 1st layer to BIG SKULL domain (9)      
    no_dom(1) = 9
! ..set 2nd and 3rd layer to AIR domain (7)
    no_dom(2:3) = 7
#if I_PRINT >= 2
    write(*,1)no_dom
1   format(' customize_geometry_NEW: layers domains = ',6(I2,2X))
!   REMARK: no_dom = (9;7,7;11,11,11)    
#endif    
! ..check PML parameters
    r1 = SURFACES(nts)%Rdata(4)
    do i = 1, 3
      r1 = r1 + thickness(i)
    enddo
    r2 = r1
    do i = 4, 6
      r2 = r2 + thickness(i)
    enddo
    if ((abs(RPML_MIN - r1) .gt. GEOM_TOL) .or. (abs(RPML_MAX - r2) .gt. GEOM_TOL)) then
      write(*,*) 'customize_geometry_NEW: INCONSISTENT PML DATA'
      write(*,*) 'RPML_MIN,RPML_MAX = ',RPML_MIN,RPML_MAX
      write(*,*) 'SPHERE RADIUSES = ',r1,r2
      stop
    endif
!
! ..add prism layers      
    call add_prism_layer(nts,nr_layers,thickness,no_dom)
#if I_PRINT >= 2
     write(*,*)'customize_geometry_NEW: prismatic layers added.'
     write(*,*)'customize_geometry_NEW: NRSURFS = ',NRSURFS
!    REMARK:  NRSURFS = 13 + 6 = 19      
#endif      
!----------------------------------------------------------------------
!  STEP 3: "drill" a hole in the skull in order to open ear canal, 
!    i.e. reset the domain number to 7 (BIG AIR) for prisms on the 
!    skull and within the ear canal. Also reproject points on the
!    brim of the ear opening to the ear canal inner cylinder.
!      

! ..allocate variables for point reprojection
    call allocate_EXTRA_PLANES
! ..loop through prisms
    do npri = 1, NRPRISM
! ....skip if not on the SKULL
      if (PRISMS(npri)%Domain .ne. 9) cycle
! ....check if within the inner cylinder of ear canal
      idec = 1
! ....loop through vertices of prism bottom face        
      do iv = 1, 3
        np = PRISMS(npri)%VertNo(iv)
        xp(1:3) = POINTS(np)%Rdata(1:3)
        s = xp(2)**2 + xp(3)**2
        s = sqrt(s)
        if (s .gt. (0.35d0 + GEOM_TOL))  idec = 0
      enddo
! ....if within the ear canal     
      if (idec .eq. 1) then
! ......set domain to BIG AIR (7)                
        PRISMS(npri)%Domain = 7
#if I_PRINT >= 2        
        write(*,*)'customize_geometry_NEW: redefining domain for prism = ',npri
        write(*,*)'customize_geometry_NEW: reprojecting vertices of top face.'
#endif          
! ......reproject the points on the inner side of the ear canal
        do iv = 1, 3
          np = PRISMS(npri)%VertNo(iv)
          xp(1:3) = POINTS(np)%Rdata(1:3)
          s = xp(2)**2 + xp(3)**2
          s = sqrt(s)
! ........if vertex of BOTTOM face is inside cylinder            
          if (abs(s - 0.35d0) .lt. GEOM_TOL) then
! ..........get corresponding vertex of TOP face and reproject                    
            np = PRISMS(npri)%VertNo(iv + 3)
#if I_PRINT >= 2 
            write(*,*)'customize_geometry_NEW: reprojecting np = ',np
            write(*,2) POINTS(np)%Rdata(1:3)
2           format(' *****  original coord. =',3(E12.5,2X))
#endif            
            call reproject_point(np,2,(/1,14/))
#if I_PRINT >= 2 
            write(*,3) POINTS(np)%Rdata(1:3)
3           format(' **  reprojected coord. =',3(E12.5,2X))
#endif            
          endif
        enddo
#if I_PRINT >= 2
      call pause
#endif      
      endif
! ..end of loop through prisms       
    enddo
! ..deallocate planes for reprojection
    call deallocate_EXTRA_PLANES
#if I_PRINT >= 2
    write(*,*)'customize_geometry_NEW: hole drilled in skull!'
#endif      
!
! ..attach all rectangles within the ear canal to cylindrical surface 1 in order to avoid duplication of surfaces
! ..loop through rectangles
    do nr = 1, NRRECTA
      idec = 1 
      do iv = 1, 4
        np = RECTANGLES(nr)%VertNo(iv)
        xp(1:3) = POINTS(np)%Rdata(1:3)
        s = xp(2)**2 + xp(3)**2
        s = sqrt(s)
! ......if vertex is not on the cylindrical surface        
        if (abs(s - 0.35d0) .gt. GEOM_TOL) then
          idec = 0
          exit
        endif
      enddo
! ....if triangle is on the surface
      if (idec .eq. 1) then
#if I_PRINT  >= 2              
        write(*,*)'customize_geometry_NEW: attaching rectangle = ',nr
#endif
! ......update rectangle type          
        RECTANGLES(nr)%Type = 'PTIRec'
        allocate(RECTANGLES(nr)%Idata(1), STAT = status)
        if (status .ne. 0 ) then
          write(*,*)'customize_geometry_NEW: Idata not allocated for nr =',nr
          stop
        endif
! ......set surface number to inner cylinder of ear canal (1)         
        RECTANGLES(nr)%Idata(1) = 1
! ....end if rectangle is on the surface        
      endif
! ..end of loop through rectangles      
    enddo
#if I_PRINT >= 2
      write(*,*)'customize_geometry_NEW: rectangles attached to cylinder.'
#endif      
!
!-----------------------------------------------------------------------------------
!  STEP 4: add membranes
!
! ..generate the ear drum
#if I_PRINT >= 2
    write(*,*)'customize_geometry_NEW: NRDOMAIN before EAR DRUM = ',NRDOMAIN
!   REMARK: NRDOMAIN = 11    
#endif
    nsplit = 11
    nr_bound = 1;  ns_bound(1) = 10
    nr_confm = 2;  ns_confm(1) = 1; ns_confm(2) = 8
    dh1 = 0.05d0;  dh2 = 0.05d0 
#if I_PRINT >= 2
    write(*,*)'**************************************************************************'
#endif    
    call split_surface(nsplit,nr_bound,ns_bound,nr_confm,ns_confm,dh1,dh2)
#if I_PRINT >= 2    
    write(*,*)'customize_geometry_NEW: EAR DRUM added.'
    write(*,*)'customize_geometry_NEW: NRDOMAIN after EAR DRUM = ',NRDOMAIN
!   REMARK: NRDOMAIN = 12    
#endif      
!
! ..generate the oval window
#if I_PRINT >= 2
    write(*,*)'customize_geometry_NEW: NRDOMAIN before OVAL WINDOW = ',NRDOMAIN
!   REMARK: NRDOMAIN = 12    
#endif
    nsplit = 2
    nr_bound = 1;  ns_bound(1) = 3
    nr_confm = 3;  ns_confm(1) = 4; ns_confm(2) = 6; ns_confm(3) = 8
    dh1 = 0.d0; dh2 = 0.02d0 
#if I_PRINT >= 2
    write(*,*)'**************************************************************************'
#endif    
    call split_surface(nsplit,nr_bound,ns_bound,nr_confm,ns_confm,dh1,dh2)
#if I_PRINT >= 2     
    write(*,*)'customize_geometry_NEW: OVAL WINDOW added.'
    write(*,*)'customize_geometry_NEW: NRDOMAIN after OVAL WINDOW = ',NRDOMAIN
!   REMARK: NRDOMAIN = 13    
#endif      
!
! ..generate the basilar membrane
#if I_PRINT >= 2
    write(*,*)'customize_geometry_NEW: NRDOMAIN before BASILAR MEMBRANE = ',NRDOMAIN
!   REMARK: NRDOMAIN = 13    
#endif
    nsplit = 6
    nr_bound = 3;  ns_bound(1) = 7; ns_bound(2) = 3; ns_bound(3) = 13
    nr_confm = 3;  ns_confm(1) = 5; ns_confm(2) = 2; ns_confm(3) = 4
    dh1 = 0.005d0; dh2 = 0.005d0 
#if I_PRINT >= 2
    write(*,*)'**************************************************************************'
#endif    
    call split_surface(nsplit,nr_bound,ns_bound,nr_confm,ns_confm,dh1,dh2)
#if I_PRINT >= 2    
    write(*,*)'customize_geometry_NEW: BASILAR MEMBRANE added.'
    write(*,*)'customize_geometry_NEW: NRDOMAIN after BASILAR MEMBRANE = ',NRDOMAIN
!   REMARK: NRDOMAIN = 14    
#endif
!
#if I_PRINT >= 1
      write(*,*) 'customize_geometry_NEW: done!'
#endif      
!      
end subroutine customize_geometry_head
