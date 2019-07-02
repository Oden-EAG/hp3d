!----------------------------------------------------------------------------
!> Purporse : reconstuction of surfaces surrounding a domain
!!
!! @param[in] Ndom - domain number
!!
!! @revision Aug 11
!----------------------------------------------------------------------------
!
subroutine smooth_domain(Ndom)
  use GMP
  implicit none  
  integer,intent(in) :: Ndom
!  
! workspaces
  integer,dimension(NRPOINT)   :: work_pt
  integer,dimension(NRCURVE)   :: work_cu
  integer,dimension(NRTRIAN,3) :: work_tr
! other local variables
  real*8,dimension(3) :: rdata_save
  integer :: i,j,k,np,nc,nt,nick,nb,or,iflag,lab,ie,ifig,is,ns,iprint,jflag, &
             ntet,npyr,npri,nflag,nr
!
!----------------------------------------------------------------------------
! STEP 0 : check that Ndom is valid
!
  if ((Ndom.le.0).or.(Ndom.gt.NRDOMAIN)) then
    write(*,*)'smooth_domain: invalid domain number!'
    stop
  endif
!
!----------------------------------------------------------------------------
! STEP 1 : record point, curve, triangle types pertaining to reconstructed 
!          geometry. Reset all of those types to linear geometry.
!
! P O I N T S       : 0 - Regular ; 1 - CoorNrm ; 2 - SharpPt  
  work_pt = 0
  do i=1,NRPOINT
    select case(POINTS(i)%Type)
    case('CoorNrm') ; work_pt(i)=1 ; POINTS(i)%Type='Regular'
    case('SharpPt') ; work_pt(i)=2 ; POINTS(i)%Type='Regular'
    endselect
  enddo
!
! C U R V E S       : 0 - Seglin ; 1 - HermCur
  work_cu = 0
  do i=1,NRCURVE
    select case(CURVES(i)%Type)
    case('HermCur') ; work_cu(i)=1 ; CURVES(i)%Type='Seglin'
    endselect
  enddo
!  
! T R I A N G L E S : 1 - G1RecTri
  work_tr = 0
  do i=1,NRTRIAN
    select case(TRIANGLES(i)%Type)
    case('G1RecTri') 
      work_tr(i,1)=1 ; work_tr(i,2:3)=TRIANGLES(i)%BlockNo(1:2)
      TRIANGLES(i)%Type   ='PlaneTri'
      TRIANGLES(i)%BlockNo=0
    endselect
  enddo
!
!----------------------------------------------------------------------------
! STEP 2 : for each tet in the domain of interest reset point, curve,
!          triangle types to the appropriate ones.
!
! loop over tets
  do i=1,NRTETRA
    if (TETRAS(i)%Domain.ne.Ndom) cycle
!   loop over faces
    do j=1,4
      nick=TETRAS(i)%FigNo(j)
      call decode(nick, nt,or)
      if (.not.associated(TRIANGLES(nt)%Idata)) cycle
      ns=TRIANGLES(nt)%Idata(1)
!     if triangle is on reconstructed surface 
      if (SURFACES(ns)%Type.eq.'RecSurf') then
!
!  .....update face............................................................
        TRIANGLES(nt)%Type='G1RecTri'
        work_tr(nt,1)=1 ; work_tr(nt,2:3)=TRIANGLES(nt)%BlockNo(1:2)
        TRIANGLES(nt)%BlockNo(1:2)=0 ; jflag=0
!       loop over adjacent blocks        
        do k=2,3
          call decode(work_tr(nt,k), nb,lab)
          if (lab.ne.3) cycle
          if (TETRAS(nb)%Domain.eq.Ndom) then
            jflag=k ; exit
          endif                  
        enddo
        if (jflag.eq.0) then
          write(*,*)'smooth_domain: neighboring block not found!'
          stop
        endif
        TRIANGLES(nt)%BlockNo(1)=work_tr(nt,jflag)
!        
!  .....update vertices........................................................        
        do k=1,3
          np = TRIANGLES(nt)%VertNo(k)
!         if vertex has already been visited, cycle
          if (POINTS(np)%Type.eq.'CoorNrm') cycle
          POINTS(np)%Type='CoorNrm' ; rdata_save(1:3)=POINTS(np)%Rdata(1:3)
          deallocate(POINTS(np)%Rdata, STAT=is)
          if (is.ne.0) then
            write(*,*)'smooth_domain: Rdata not deallocated for np = ',np
            stop
          endif
          allocate(POINTS(np)%Rdata(6), STAT=is)
          if (is.ne.0) then
            write(*,*)'smooth_domain: Rdata not allocated for np = ',np
            stop
          endif
          POINTS(np)%Rdata(1:3) = rdata_save(1:3)
        enddo
!        
!  .....update edges...........................................................        
        do k=1,3
          nc=iabs(TRIANGLES(nt)%EdgeNo(k))
          CURVES(nc)%Type='HermCur'
        enddo
!
      else
        write(*,*)'smooth_domain: have found a tet with a face on an algebraic surf'
        call pause
      endif
    enddo
!
! loop over tets
  enddo
!
!----------------------------------------------------------------------------
! STEP 3: perfomr reconstruction
!
  call reconstruct
!
!----------------------------------------------------------------------------
! STEP 4: restore original types for geometrical entities
!
  do i=1,NRPOINT
    select case(work_pt(i))
    case(1) ; POINTS(i)%Type='CoorNrm'
    case(2) ; POINTS(i)%Type='SharpPt'
    endselect
  enddo
!
  do i=1,NRCURVE
    select case(work_cu(i))
    case(1) 
      if(CURVES(i)%Type.eq.'5Bezier') cycle
      if(CURVES(i)%Type.eq.'7Bezier') cycle
      CURVES(i)%Type='HermCur'
    endselect
  enddo
!  
  do i=1,NRTRIAN
    select case(work_tr(i,1))
    case(1) 
      TRIANGLES(i)%Type         ='G1RecTri'
      TRIANGLES(i)%BlockNo(1:2) = work_tr(i,2:3)
    endselect
  enddo
!
!----------------------------------------------------------------------------
! STEP 5: redefine curvilinear triangles and rectangles
!       
  do i=1,NRTRIAN
    if (TRIANGLES(i)%Type.eq.'PlaneTri') then
      iflag=0
      do j=1,3
        nc=iabs(TRIANGLES(i)%EdgeNo(j))
        if (CURVES(nc)%Type.ne.'Seglin') iflag=1
      enddo
      if (iflag.eq.1) TRIANGLES(i)%Type = 'TransTri'
    endif
  enddo
!
  do i=1,NRRECTA
    if (RECTANGLES(i)%Type.eq.'BilQua') then
      iflag=0
      do j=1,4
        nc=iabs(RECTANGLES(i)%EdgeNo(j))
        if (CURVES(nc)%Type.ne.'Seglin') iflag=1
      enddo
      if (iflag.eq.1) RECTANGLES(i)%Type = 'TraQua'
    endif
  enddo
!
!----------------------------------------------------------------------
! STEP 6: redefine curvilinear tets, prisms and pyramids
!       
  do ntet=1,NRTETRA
    nflag=0
    do ie=1,6
      nc = iabs(TETRAS(ntet)%EdgeNo(ie))
      if (CURVES(nc)%Type.ne.'Seglin') nflag=1
    enddo
    do ifig=1,4
      if (TETRAS(ntet)%FigNo(ifig).eq.0) cycle
      call decode(TETRAS(ntet)%FigNo(ifig), nt,lab)
      if (TRIANGLES(nt)%Type.ne.'PlaneTri') nflag=1
    enddo
    if (nflag.eq.1) then
      if (iprint.eq.1) then
        write(*,*) 'input_RECONSTRUCT: REDEFINING ntet = ',ntet
      endif
      TETRAS(ntet)%Type = 'TraTet'
    endif
  enddo  
!
  do npri=1,NRPRISM
    nflag=0
    do ie=1,9
      nc = iabs(PRISMS(npri)%EdgeNo(ie))
      if (CURVES(nc)%Type.ne.'Seglin') nflag=1
    enddo
    do ifig=1,2
      if (PRISMS(npri)%FigNo(ifig).eq.0) cycle
      call decode(PRISMS(npri)%FigNo(ifig), nt,lab)
      if (TRIANGLES(nt)%Type.ne.'PlaneTri') nflag=1
    enddo
    do ifig=3,5
      if (PRISMS(npri)%FigNo(ifig).eq.0) cycle
      call decode(PRISMS(npri)%FigNo(ifig), nr,lab)
      if (RECTANGLES(nr)%Type.ne.'BilQua') nflag=1
    enddo
    if (nflag.eq.1) then
      if (iprint.eq.1) then
        write(*,*) 'input_RECONSTRUCT: REDEFINING npri = ',npri
      endif
      PRISMS(npri)%Type = 'TIprism'
    endif
  enddo  
!
  do npyr=1,NRPYRAM
    nflag=0
    do ie=1,8
      nc = iabs(PYRAMIDS(npyr)%EdgeNo(ie))
      if (CURVES(nc)%Type.ne.'Seglin') nflag=1
    enddo
    do ifig=1,1
      if (PYRAMIDS(npyr)%FigNo(ifig).eq.0) cycle
      call decode(PYRAMIDS(npyr)%FigNo(ifig), nr,lab)
      if (RECTANGLES(nr)%Type.ne.'BilQua') nflag=1
    enddo
    do ifig=2,5
      if (PYRAMIDS(npyr)%FigNo(ifig).eq.0) cycle
      call decode(PYRAMIDS(npyr)%FigNo(ifig), nt,lab)
      if (TRIANGLES(nt)%Type.ne.'PlaneTri') nflag=1
    enddo
    if (nflag.eq.1) then
      if (iprint.eq.1) then
        write(*,*) 'input_RECONSTRUCT: REDEFINING npyr = ',npyr
      endif
      PYRAMIDS(npyr)%Type = 'TIpyram'
    endif
  enddo  
!
!
endsubroutine smooth_domain
