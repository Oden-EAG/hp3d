!----------------------------------------------------------------------
!> Purpose : check costrained approximation
!!
!!   LOOP over blocks
!!
!!     LOOP over face orientations
!!
!!       LOOP over refinement kinds
!!
!!         solve projections
!!
!!       END
!!     END
!!   END
!!
!> @date : Nov 14
!----------------------------------------------------------------------
! Jan 15
!
! This test is really intended to by used on the hybrid mesh (prism +
! hexa + tet [ + pyramid ] ). For the non interactive mode, it has been
! simplified as follows:
!
!  LOOP over refinements of 1st block in DS (i.e., the prism)
!    solve projections
!  END
!  
!----------------------------------------------------------------------
!
subroutine test_3
!      
      use environment , only : FILE_PHYS,QUIET_MODE
      use refinements
      use data_structure3D
      use GMP
      use PROJ , only : ITAG,FILE_TAGS,ITEST
!
      implicit none
!      
      character(len=4) :: btype , ftype
      integer,dimension(4) :: nv
      integer :: i,j,k,nc,mdle,nbln,lab,nor,iv,nfig,inod,nodf, &
                 nick,nvoid,iref,kref,iblock,ic,l,idec,iface
!                 
      integer, parameter :: nin = 23
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
!         header
          write(*,*) ''
          write(*,*) 'Blocks:'
!
!         loop over initial mesh blocks
          do iblock=1,NRELIS
!          
            mdle=iblock
!
!           display block info
            write(*,8000) mdle,NODES(mdle)%Type
 8000       format(' mdle = 'i4,' ; type = ',a4 )         
          enddo
          write(*,*) ''
!
!         display refinement kind
          call display_ref_kind
!
!         select block
          write(*,*)'Set: mdle ='
          read( *,*) mdle
!          
          btype=NODES(mdle)%Type
          call decode(ELEMS(mdle)%GMPblock, nbln,lab)
!
!         select face
          write(*,*)'Set: face (0 - skip to refine) ='
          read( *,*) iface
!
!         set to prevent funcky printing
          j=0 ; ftype='void'
!
!         rotate face
          if (iface > 0) then
!
            select case(lab)
            case(1) ; nick=PRISMS(  nbln)%FigNo(iface)
            case(2) ; nick=HEXAS(   nbln)%FigNo(iface)
            case(3) ; nick=TETRAS(  nbln)%FigNo(iface)
            case(4) ; nick=PYRAMIDS(nbln)%FigNo(iface)
            endselect
!            
!           face figure number        
            call decode(nick, nfig,nvoid)
!            
!           face type
            inod=nvert(btype)+nedge(btype)+iface
            nodf=ELEMS(mdle)%nodes(inod)
            ftype=NODES(nodf)%type
!      
!           number of orientations for the face
            select case(ftype)
            case('mdlt') ; nor=5
            case('mdlq') ; nor=7
            endselect
!
!           select orientation
            write(*,8001) nor
 8001       format(' Set: orientation (1 - ',i1,') =')             
            read(*,*) j
!
!           rotate vertices
            select case(ftype)
!
!           triangular face          
            case('mdlt')
              do k=1,3
                iv=TRIAN_L2G(k,j) ; nv(k)=TRIANGLES(nfig)%VertNo(iv)
              enddo
              TRIANGLES(nfig)%VertNo(1:3)=nv(1:3)
!
!           rectangular face          
            case('mdlq')
              do k=1,4
                iv=QUADR_L2G(k,j) ; nv(k)=RECTANGLES(nfig)%VertNo(iv)
              enddo
              RECTANGLES(nfig)%VertNo(1:4)=nv(1:4)
            endselect
          endif
!
!         complete connectivities
          call clean_GMP
          call connect
!          
!         deallocate ds, generate initial mesh
          call cleanup_href
!
!         refinement kind
          write(*,*)'Set: refinement (0 - skip to solve) ='
          read( *,*) kref
!
!         refine element
          if (kref > 0) then
            call refine(mdle,kref)
            call close
            call update_gdof
            call update_ddof
            call verify_orient
            call verify_neig
          endif
!
!         increment tag
          ITAG=ITAG+1
!          
!         print
          write(*,*) ''
          write(*,*) '-- Test 3 --'
          write(*,8002) ITAG,btype,ftype,iface,j,kref
 8002     format(' tag = ',i2,' --> type = ',a4,' ; face (',a4,') = ',i1,' ; orientation = ',i1, &
                                                                         ' ; refinement = ',i3)         
          write(*,*) ''
!
!         solve projection problems
          call test_1
!
          write(*,*)'0 - Return ; 1 - Continue testing'
          read( *,*) idec
!
        enddo
!
!       clean up refinement
        call cleanup_href
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
!       loop over faces
!
! P. Gatto, Jan 15
!!!        do i=1,nface(btype)
        do i=1,1
!        
          select case(lab)
          case(1) ; nick=PRISMS(  nbln)%FigNo(i)
          case(2) ; nick=HEXAS(   nbln)%FigNo(i)
          case(3) ; nick=TETRAS(  nbln)%FigNo(i)
          case(4) ; nick=PYRAMIDS(nbln)%FigNo(i)
          endselect
!        
!         face figure number        
          call decode(nick, nfig,nvoid)
!          
!         face type
          inod=nvert(btype)+nedge(btype)+i
          nodf=ELEMS(mdle)%nodes(inod)
          ftype=NODES(nodf)%type
!
!         number of orientations for the face
          select case(ftype)
          case('mdlt') ; nor=5
          case('mdlq') ; nor=7
          endselect
!
!         loop over orientations (skip 0-th one)
!
! P. Gatto, Jan 15
!!!          do j=1,nor
          do j=1,1
!
!           rotate vertices
            select case(ftype)
!
!           triangular face          
            case('mdlt')
              do k=1,3
                iv=TRIAN_L2G(k,j) ; nv(k)=TRIANGLES(nfig)%VertNo(iv)
              enddo
              TRIANGLES(nfig)%VertNo(1:3)=nv(1:3)
!              
!           rectangular face          
            case('mdlq')
              do k=1,4
                iv=QUADR_L2G(k,j) ; nv(k)=RECTANGLES(nfig)%VertNo(iv)
              enddo
              RECTANGLES(nfig)%VertNo(1:4)=nv(1:4)
            endselect
!
!           complete connectivities
            call clean_GMP
            call connect
!            
!           loop over refinement kinds
            do iref=1,nr_ref(btype)
!
!             deallocate ds, generate initial mesh
              call cleanup_href
!
!             determine refinement kind
              kref=kref_kind(iref, btype)
!
!             refine element
              call refine(mdle,kref)
              call close
              call update_gdof
              call update_ddof
              call verify_orient
              call verify_neig
!              
!             increment tag
              ITAG=ITAG+1
!
!             file for printing
!
!             -- 1st visit --
              if (ITAG == 1) then
!
!               open file for printing      
                open(unit=nin,file=trim(FILE_TAGS),form='formatted',access='sequential', & 
                                                   status='replace',iostat=ic)
                if (ic /= 0) then
                  write(*,*)'test_3: COULD NOT OPEN FILE! [0]'
                  stop
                endif
!
!               print header to file
                write(nin,6000)
 6000           format(' Tag    // Block Type // Face // Orientation // Refinement')
!
!             -- subsequent visits --             
              else
!                      
!               append to file
                open(unit=nin,file=trim(FILE_TAGS),form='formatted',access='sequential', &
                                                   status='old',position='append',iostat=ic)
                if (ic /= 0) then
                  write(*,*)'test_3: COULD NOT OPEN FILE! [1]'
                  stop
                endif
              endif
!              
!             print to file
              write(nin,7000) ITAG,btype,i,j,kref
!
!             close file
              close(unit=nin,iostat=ic)
              if (ic /= 0) then
                 write(*,*)'test_3: COULD NOT CLOSE FILE!'
                 stop
              endif
!
!             print to screen
IF (.NOT. QUIET_MODE) THEN              
!        
!             check
              if (ITAG > maxvis) then
                write(*,*)'test_3: INCREASE maxvis!'
                stop
              endif
!
!             store
              iwork(ITAG,1) = ITAG 
              iwork(ITAG,2) = i
              iwork(ITAG,3) = j 
              iwork(ITAG,4) = kref 
              swork(ITAG,1) = btype 
!
!             print
              write(*,*)'' 
              write(*,6000)
              do l=1,ITAG
                write(*,7000) iwork(l,1),swork(l,1),iwork(l,2:4)
 7000           format(1x,i6,' ; ',1x,a4,' ; ',7x,i1,' ; ',4x,i1,' ; ',11x,i3)               
              enddo
              write(*,*)'' 
ENDIF
!
!             solve projection problems
              call test_1
!
!           end loop over refinement kind
            enddo
!
!         end loop over orientations
          enddo
!        
!       end loop over faces        
        enddo
!        
!     end loop over blocks
      enddo
!
!
endsubroutine test_3
!
!
!
subroutine cleanup_href
!
        use environment      , only : FILE_PHYS
        use data_structure3D , only : deallocds
!
        implicit none

!       clean up refinements
        call deallocds
        call hp3gen(trim(FILE_PHYS))
!
!
endsubroutine cleanup_href
