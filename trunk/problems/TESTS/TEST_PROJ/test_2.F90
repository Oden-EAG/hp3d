!----------------------------------------------------------------------
!> Purpose : check orientations
!!
!!   LOOP over blocks
!!
!!     LOOP over edges 
!!
!!       change orientation
!!       solve projections
!!
!!     END
!!
!!     LOOP over faces
!!
!!       LOOP over orientations
!!
!!         solve projections
!!
!!       END
!!     END
!!   END
!!
!> @date : Nov 14
!----------------------------------------------------------------------
!
subroutine test_2(Isol)
!      
      use environment       , only : FILE_PHYS,QUIET_MODE
      use parameters
      use data_structure3D
      use GMP
      use PROJ              , only : FILE_TAGS,ITAG,ITEST
!
      implicit none
      integer, intent(in) :: Isol
!      
      character(len=4) :: btype , ftype
      integer,dimension(4) :: nv
      integer :: i,j,k,nc,mdle,nbln,lab,nor,iv,nfig,inod,nodf,nick,nvoid,iblock
      integer, parameter :: nin = 24
      integer :: l,ic,iface,iedge,idec
!      
      integer, parameter :: maxvis = 2000
      integer         , dimension(maxvis,4), save :: iwork
      character(len=4), dimension(maxvis,1), save :: swork
!
!----------------------------------------------------------------------
!
!     -- INTERACTIVE MODE --
      if (ITEST == 0) then
!
        idec=1
        do while(idec /= 0)
!        
!         display blocks info
          write(*,*) ''
          write(*,*) 'Blocks:'
!
!         loop over INITIAL mesh elements
          do iblock=1,NRELIS
!
!           middle nodes enumeration coincides with block enumeration
            mdle=iblock
!
            write(*,8000) mdle,NODES(mdle)%Type
 8000       format(' mdle = 'i4,' ; type = ',a4 )         
          enddo
          write(*,*) ''
!
!         select block
          write(*,*)'Set: mdle ='
          read( *,*) mdle
!          
          btype=NODES(mdle)%Type
          call decode(ELEMS(mdle)%GMPblock, nbln,lab)
!
!         select edge
          write(*,*)'Set: edge (0 - skip to faces) ='
          read( *,*) iedge
!
!         rotate edge
          if (iedge > 0) then
!                  
            select case(lab)
            case(1) ; nc=PRISMS(  nbln)%EdgeNo(iedge)
            case(2) ; nc=HEXAS(   nbln)%EdgeNo(iedge)
            case(3) ; nc=TETRAS(  nbln)%EdgeNo(iedge)
            case(4) ; nc=PYRAMIDS(nbln)%EdgeNo(iedge)
            endselect
!
!           edge curve number
            nc=iabs(nc)
!            
!           swap vertices
            iv=EDGE_L2G(1,1) ; nv(1)=CURVES(nc)%EndPoNo(iv)
            iv=EDGE_L2G(2,1) ; nv(2)=CURVES(nc)%EndPoNo(iv)
            CURVES(nc)%EndPoNo(1:2)=nv(1:2)
!
!           just for later printing
            j=1 ; iface=0 ; ftype = 'void'
          endif
!
!         select face
          if (iedge == 0) then
            write(*,*)'Set: face (0 - skip to solve) ='
            read( *,*) iface
!
!           rotate face
            if (iface > 0) then
!
              select case(lab)
              case(1) ; nick=PRISMS(  nbln)%FigNo(iface)
              case(2) ; nick=HEXAS(   nbln)%FigNo(iface)
              case(3) ; nick=TETRAS(  nbln)%FigNo(iface)
              case(4) ; nick=PYRAMIDS(nbln)%FigNo(iface)
              endselect
!              
!             face figure number        
              call decode(nick, nfig,nvoid)
!              
!             face type
              inod=nvert(btype)+nedge(btype)+iface
              nodf=ELEMS(mdle)%nodes(inod)
              ftype=NODES(nodf)%type
!      
!             number of orientations for the face
              select case(ftype)
              case('mdlt') ; nor=5
              case('mdlq') ; nor=7
              endselect
!
!             select orientation
              write(*,8001) nor
 8001         format(' Set: orientation (1 - ',i1,') =')             
              read(*,*) j
!
!             rotate vertices
              select case(ftype)
!
!             triangular face          
              case('mdlt')
                do k=1,3
                  iv=TRIAN_L2G(k,j) ; nv(k)=TRIANGLES(nfig)%VertNo(iv)
                enddo
                TRIANGLES(nfig)%VertNo(1:3)=nv(1:3)
!
!             rectangular face          
              case('mdlq')
                do k=1,4
                  iv=QUADR_L2G(k,j) ; nv(k)=RECTANGLES(nfig)%VertNo(iv)
                enddo
                RECTANGLES(nfig)%VertNo(1:4)=nv(1:4)
              endselect
            endif
          endif
!
!         complete connectivities
          call clean_GMP
          call connect
!          
!         generate initial mesh
          call deallocds
          call hp3gen(trim(FILE_PHYS))
!
!         increment tag
          ITAG=ITAG+1
!
!         prevent funcky printing
          if ((iedge == 0) .and. (iface == 0)) then
            j=0 ; ftype='void'
          endif
!          
!         print
          write(*,*) ''
          write(*,*) '-- Test 2 --'
          write(*,8002) ITAG,btype,iedge,ftype,iface,j
 8002     format(' tag = ',i2,' --> type = ',a4,' ; edge = ',i2,' ; face (',a4,') = ',i1,' ; orientation = ',i1)         
          write(*,*) ''
!
!         solve projection problems
          call test_1
!
          write(*,*)'Test 2: 0 - Return ; 1 - Continue testing'
          read( *,*) idec
        enddo
!        
!       exit here
        return
      endif        
!
!     -- TESTING MODE --
!
!     loop over initial mesh blocks
!
! P. Gatto, Jan 15
!!!      do iblock=1,NRELIS
      do iblock=1,1
        mdle=iblock
!  
        btype=NODES(mdle)%Type
        call decode(ELEMS(mdle)%GMPblock, nbln,lab)
!      
!======================================================================
!  E D G E S                                                          |
!======================================================================
!
!     loop over edges     
      do i=1,nedge(btype)
        select case(lab)
        case(1) ; nc=PRISMS(  nbln)%EdgeNo(i)
        case(2) ; nc=HEXAS(   nbln)%EdgeNo(i)
        case(3) ; nc=TETRAS(  nbln)%EdgeNo(i)
        case(4) ; nc=PYRAMIDS(nbln)%EdgeNo(i)
        endselect
!
!       edge curve number
        nc=iabs(nc)
!        
!       rotate edge
        iv=EDGE_L2G(1,1) ; nv(1)=CURVES(nc)%EndPoNo(iv)
        iv=EDGE_L2G(2,1) ; nv(2)=CURVES(nc)%EndPoNo(iv)
        CURVES(nc)%EndPoNo(1:2)=nv(1:2)
!
!       complete connectivities
        call clean_GMP
        call connect
!        
!       generate initial mesh
        call deallocds
        call hp3gen(trim(FILE_PHYS))
!
!       increment tag
        ITAG=ITAG+1
!
!       file for printing
!
!       -- 1st visit -- 
        if (ITAG == 1) then

!         open file for printing      
          open(unit=nin,file=trim(FILE_TAGS),form='formatted',access='sequential',status='replace',iostat=ic)
          if (ic /= 0) then
            write(*,*)'test_2: COULD NOT OPEN FILE! [0]'
            stop
          endif
!
!         print header to file
          write(nin,6000)
 6000     format(' Tag    // Block Type // Edge // Face // Orientation')
!
!       -- subsequent visits
        else 
!                
!         append to file
          open(unit=nin,file=trim(FILE_TAGS),form='formatted',access='sequential',status='old',position='append',iostat=ic)
          if (ic /= 0) then
            write(*,*)'test_2: COULD NOT OPEN FILE! [1]'
            stop
          endif
        endif
!
!       print to file
        write(nin,7000) ITAG,btype,i,0,1
 7000   format(1x,i6,' ; ',1x,a4,' ; ',7x,i2,' ; ',3x,i1,' ; ',14x,i1)               
!
!       print to screen
IF (.NOT. QUIET_MODE) THEN              
!       
!       check
        if (ITAG > maxvis) then
          write(*,*)'test_2: INCREASE maxvis!'
          stop
        endif
!
!       store
        iwork(ITAG,1) = ITAG 
        iwork(ITAG,2) = i
        iwork(ITAG,3) = 0 
        iwork(ITAG,4) = 1
        swork(ITAG,1) = btype 
!
!       print
        write(*,*)'' 
        write(*,6000)
        do l=1,ITAG
          write(*,7000) iwork(l,1),swork(l,1),iwork(l,2:4)
        enddo
        write(*,*)'' 
ENDIF
!        
!       close file
        close(unit=nin,iostat=ic)
        if (ic /= 0) then
           write(*,*)'test_2: COULD NOT CLOSE FILE! [0]'
           stop
        endif
!
!       solve projection problems
        call test_1
!
!     end loop over edges
      enddo
!
!======================================================================
!  F A C E S                                                          |
!======================================================================
!
!     loop over faces
      do i=1,nface(btype)
        select case(lab)
        case(1) ; nick=PRISMS(  nbln)%FigNo(i)
        case(2) ; nick=HEXAS(   nbln)%FigNo(i)
        case(3) ; nick=TETRAS(  nbln)%FigNo(i)
        case(4) ; nick=PYRAMIDS(nbln)%FigNo(i)
        endselect
!        
!       face figure number        
        call decode(nick, nfig,nvoid)
!        
!       face type
        inod=nvert(btype)+nedge(btype)+i
        nodf=ELEMS(mdle)%nodes(inod)
        ftype=NODES(nodf)%type
!
!       number or orientations for the face
        select case(ftype)
        case('mdlt') ; nor=5
        case('mdlq') ; nor=7
        endselect
!
!       loop over orientations (skip 0-th one)
        do j=1,nor
!
!         rotate vertices
          select case(ftype)
!
!         triangular face          
          case('mdlt')
            do k=1,3
              iv=TRIAN_L2G(k,j) ; nv(k)=TRIANGLES(nfig)%VertNo(iv)
            enddo
            TRIANGLES(nfig)%VertNo(1:3)=nv(1:3)
!
!         rectangular face          
          case('mdlq')
            do k=1,4
              iv=QUADR_L2G(k,j) ; nv(k)=RECTANGLES(nfig)%VertNo(iv)
            enddo
            RECTANGLES(nfig)%VertNo(1:4)=nv(1:4)
          endselect
!
!         complete connectivities
          call clean_GMP
          call connect
!          
!         generate initial mesh
          call deallocds
          call hp3gen(trim(FILE_PHYS))
!          
!         increment tag
          ITAG=ITAG+1
!          
!         file for printing
!
!         -- 1st visit -- (should not happen, but maybe one day you want to skip edges)
          if (ITAG == 1) then

!           open file for printing      
            open(unit=nin,file=trim(FILE_TAGS),form='formatted',access='sequential',status='replace',iostat=ic)
            if (ic /= 0) then
              write(*,*)'test_2: COULD NOT OPEN FILE! [2]'
              stop
            endif
!
!           print header to file
            write(nin,6000)
!
!         -- subsequent visits
          else 
!                  
!           append to file
            open(unit=nin,file=trim(FILE_TAGS),form='formatted',access='sequential',status='old',position='append',iostat=ic)
            if (ic /= 0) then
              write(*,*)'test_2: COULD NOT OPEN FILE! [3]'
              stop
            endif
          endif
!          
!         print to file
          write(nin,7000) ITAG,btype,0,i,j
!          
!         print to screen
IF (.NOT. QUIET_MODE) THEN              
!       
!         check
          if (ITAG > maxvis) then
            write(*,*)'test_2: INCREASE maxvis!'
            stop
          endif
!
!         store
          iwork(ITAG,1) = ITAG 
          iwork(ITAG,2) = 0
          iwork(ITAG,3) = i 
          iwork(ITAG,4) = j
          swork(ITAG,1) = btype 
!
!         print
          write(*,*)'' 
          write(*,6000)
          do l=1,ITAG
            write(*,7000) iwork(l,1),swork(l,1),iwork(l,2:4)
          enddo
          write(*,*)'' 
ENDIF
!
!         close file
          close(unit=nin,iostat=ic)
          if (ic /= 0) then
             write(*,*)'test_2: COULD NOT CLOSE FILE! [1]'
             stop
          endif
!
!         solve projection problems
          call test_1
!
!       end loop over orientations
        enddo
!        
!     end loop over faces        
      enddo
!
!     end loop over blocks      
      enddo
!
!
endsubroutine test_2
