!-------------------------------------------------------------------------------------
!> Purpose : check if all faces attached to an edge share a particular BC flag
!!
!!
!! @param[in]  Medge        - an edge node
!! @param[in]  NBCflag      - a BC flag  
!! @param[in]  Icomp        - component number
!! @param[out] if_pure_BC__edge
!!                          = .true.  if all adjacent faces have the same BC flag
!!                            for the 'Icomp' component
!!                          - .false. elsewhere     
!!
!! @revision Jun 21
!-------------------------------------------------------------------------------------
!
      logical function if_pure_BC_edge(Medge,NBCflag,Icomp)
!
      use data_structure3D
      use element_data
      implicit none
!
!  ...BC flags
      integer :: ibc(6,NRINDEX)
!
      integer :: Medge,NBCflag,Icomp
!
      integer, parameter :: maxn=10
      integer :: nrneig,i,mdle,nofaces(2),j
      integer :: neig(maxn), nedg_list(maxn), norient_list(maxn)
!
      if_pure_BC_edge = .true.
!
!  ...determine mdle node neighbors of the edge 
      call neig_edge(Medge,maxn, nrneig,neig,nedg_list,norient_list)
!
!  ...loop through the neighbors
      do i=1,nrneig
        mdle = neig(i)
!
!  .....determine BC flags for the neighbor
        call find_bc(mdle, ibc)
!
!  .....determine faces attached to the edge
        call edge_to_faces(NODES(mdle)%type,nedg_list(i), nofaces)
        do j=1,2
          if (ibc(nofaces(j),Icomp).ne.NBCflag) then
            if_pure_BC_edge = .false.
            return
          endif
        enddo
      enddo
!
!
      end function if_pure_BC_edge

      
