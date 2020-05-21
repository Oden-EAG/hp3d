
!--------------------------------------------------------------------------
!> Purpose : verify routines find_neig and nface_neig againts each other
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
 7010   format(' verify_neig: VERIFYING ROUTINES find_neig AND nface_neig')
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






