
!--------------------------------------------------------------------------
!> Purpose : verify routines find_neig and neig_face againts each other
!            
!! @revision May 20
!--------------------------------------------------------------------------
!
      subroutine verify_neig
!
      use element_data
      use data_structure3D
      implicit none
!
!  ...locals
      character(len=4) :: type
      integer, dimension(27) :: nodesl, norientl
!
!  ...work space for find_neig
      integer                 :: mdle,iel,nve,nrf,i,j,mdlen,k,nflag,loc
      integer, dimension(4,6) :: neig_list, neign_list
!
!  ...work space for neig_face
      integer               :: mface,nrneig,nsidn,nrfn
      integer, dimension(2) :: neig,nsid_list, norient_list
!
      integer :: iprint
!
!-------------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7010) 
 7010   format(' verify_neig: VERIFYING ROUTINES find_neig AND neig_face')
      endif
      nflag=0
!
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
!
!  .....determine the element nodes
        call elem_nodes(Mdle, nodesl,norientl)
!
!  .....determine neighbors of the element
        call find_neig(mdle, neig_list)
        type = NODES(mdle)%type
        nve = nvert(type)+nedge(type)
        nrf = nface(type)
        if (iprint.eq.1) then
          write(*,7020) neig_list(1:4,1:nrf)
 7020     format(' neig_list = ',6(4i5,3x))
        endif
        do i=1,nrf
          mface = nodesl(nve+i)
          call neig_face(mface, nrneig,neig,nsid_list,norient_list)
          if (nrneig.eq.1) cycle
          if (neig(1).eq.mdle) then
            nsidn = nsid_list(2)
          elseif (neig(2).eq.mdle) then
            nsidn = nsid_list(1)
          else
            write(*,*) 'verify_neig: INCONSISTENVY 1'
            call pause
            nflag=1
          endif
          do j=1,4
            mdlen = neig_list(j,i)
!
!  .........determine neighbors for the neighboring element
            call find_neig(mdlen, neign_list)
!
!  .........check if 'mdle' is on the list of neighbors of 'mdlen'
            call locate(mdle,neign_list(1:4,nsidn),4,loc)
            if (loc.eq.0) then
              write(*,*) 'verify_neig: INCONSISTENCY 2'
              write(*,7030) mdle,i,j,mdlen,nsidn
 7030         format(' verify_neig:  mdle,i,j,mdlen,nsidn = ',i6,3x,2i2,3x,i6,i3)
              nrfn = nface(NODES(mdlen)%type)
              write(*,7040) neign_list(1:4,1:nrfn)
 7040         format(' neign_list = ',6(4i5,3x))
              call pause
              nflag=1
            endif
          enddo
        enddo
!
!  ...end of loop through elements
      enddo
      select case(nflag)
      case(0)
        write(*,*) 'verify_neig: PASSED THE CONSISTENCY TEST'
      case(1)
        write(*,*) 'verify_neig: FAILED THE CONSISTENCY TEST'
      end select
!
      end subroutine verify_neig




!--------------------------------------------------------------------------
!> Purpose : verify routines elem_nodes and neig_edge against each other
!            
!! @revision May 20
!--------------------------------------------------------------------------
!
      subroutine verify_neig_edge
!
      use element_data
      use data_structure3D
      implicit none
      common /cneig_edge/ iprint_neig_edge
      integer :: iprint_neig_edge
!
!  ...locals
      character(len=4) :: type
      integer, dimension(27) :: nodesl, norientl
!
!  ...work space for neig_edge
      integer            :: medge
      integer, parameter :: maxn=20
      integer            :: nrneig
      integer, dimension(maxn) :: neig, nedg_list, norient_list
!
      integer :: iprint, iel, mdle, nrv, nre, ie, loc, nflag
!
!-------------------------------------------------------------------------
!
      iprint_neig_edge=0
      iprint=0
      if (iprint.eq.1) then
        write(*,7010) 
 7010   format(' verify_neig_edge: VERIFYING ROUTINES elem_nodes AND neig_edge')
      endif
      nflag=0
!
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
!
!  .....determine the element nodes
        call elem_nodes(mdle, nodesl,norientl)
        type = NODES(mdle)%type
        nrv = nvert(type); nre = nedge(type)
        if (iprint.eq.1) then
          write(*,7015) iel,mdle
 7015     format(' verify_neig_edge: iel,mdle  = ',i4,i10)
          write(*,7020) mdle, nodesl(nrv+1:nrv+nre)
        endif
        do ie=1,nre
   20     continue
          medge = nodesl(nrv+ie)
          if (iprint.eq.1) then
            write(*,7016) ie,medge
 7016       format('                   ie, medge = ',i4,i10)
          endif
          call neig_edge(medge,maxn, nrneig,neig,nedg_list,norient_list)
!
!  .......look for the middle on the list of edge neighbors
          call locate(mdle, neig(1:nrneig),nrneig, loc)
          if (loc.eq.0) then
            write(*,*) 'verify_neig_edge: INCONSISTENCY 1'
            write(*,7020) mdle, nodesl(nrv+1:nrv+nre)
 7020       format(' mdle = ',i6,' edges = ',12i6)
            write(*,7030) ie,medge,neig(1:nrneig)
 7030       format(' ie = ',i2,' medge = ',i6,' neig = ',10i6)
            call result
            nflag=1; iprint_neig_edge=1
            go to 20
          endif
          if ((nedg_list(loc).ne.ie).or.(norient_list(loc).ne.norientl(nrv+ie))) then
            write(*,*) 'verify_neig_edge: INCONSISTENCY 2'
            write(*,7020) mdle, nodesl(nrv+1:nrv+nre)
            write(*,7030) ie,medge,neig(1:nrneig)
            write(*,7040) nedg_list(1:nrneig)
 7040       format(' nedg_list    = ',12i6)
            write(*,7050) norient_list(1:nrneig)
 7050       format(' norient_list = ',12i6)
            call result
            nflag=1
          endif
        enddo
!
!  ...end of loop through elements
      enddo
      select case(nflag)
      case(0)
        write(*,*) 'verify_neig: PASSED THE CONSISTENCY TEST'
      case(1)
        write(*,*) 'verify_neig: FAILED THE CONSISTENCY TEST'
      end select
!
      end subroutine verify_neig_edge







