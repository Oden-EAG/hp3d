!--------------------------------------------------------------------
!                                                                     
!     routine name      - set_initial_mesh
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - May 20
!                                                                     
!     purpose:          - runs interactively simple tests 
!                                                                    
!     arguments:                                                     
!                                                                     
!
!---------------------------------------------------------------------
!    
      subroutine my_tests
!   
      use GMP
      use element_data
      use data_structure3D
!   
      implicit none
      common /common_recta_TraQua/ iprint_recta_TraQua
      common /common_curve_SegCir/ iprint_curve_SegCir
!

!----------------------------------------------------------------------
!  
      integer :: iprint,iselect,iprint_recta_TraQua,iprint_curve_SegCir
!
!  ...work space for recta_TraQua, curve_SegCir  
      integer               :: no
      real*8,dimension(2  ) :: eta
      real*8,dimension(3  ) :: x
      real*8,dimension(3,2) :: dxdeta
!
!  ...work space for find_neig
      integer                 :: mdle,iel,nrf,i
      integer, dimension(4,6) :: neig_list
!
!  ...work space for neig_edge
      integer            :: medge
      integer, parameter :: maxn=20
      integer            :: nrneig
      integer, dimension(maxn) :: neig, nedg_list, norient_list
      integer                  :: nrv, ie
      integer, dimension(27)   :: nodesl, norientl
!
!------------------------------------------------------------------------------------
!
      iprint=0
   10 write(*,*) 'my_tests: SELECT:'
      write(*,*) 'EXIT...................................0'
      write(*,*) 'TEST recta_TraQua......................1'
      write(*,*) 'TEST curve_SegCir......................2'
      write(*,*) 'TEST find_neig.........................3'
      write(*,*) 'verify_neig............................4'
      write(*,*) 'TEST neig_edge.........................5'
      write(*,*) 'verify_neig_edge.......................6'
      read(*,*) iselect
!
      select case(iselect) 
      case(0)
        return
      case(1)
        write(*,*) 'my_tests: INPUT No, Eta(1:2)'
        read(*,*) no, eta(1:2)
        iprint_recta_TraQua=1
        call recta_TraQua(no,eta, x,dxdeta)
      case(2)
        write(*,*) 'my_tests: INPUT No, Eta'
        read(*,*) no, eta(1)
        iprint_curve_SegCir=1
        call curve_SegCir(no,eta(1), x,dxdeta)
      case(3)
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          call find_neig(mdle, neig_list)
          nrf = nface(NODES(mdle)%type)
          write(*,7010) mdle, (neig_list(1:4,i),i=1,nrf)
 7010     format('my_tests: mdle = ',i6,' neig_list = ',6(4i6,3x))
        enddo
      case(4)
        call verify_neig
      case(5)
        write(*,*) 'input mdle, ie'
        read(*,*)  mdle,ie
        call elem_nodes(mdle, nodesl,norientl)
        nrv = nvert(NODES(mdle)%type)
        medge = nodesl(nrv+ie)
        call neig_edge(medge,maxn, nrneig,neig,nedg_list,norient_list)
      case(6)
        call verify_neig_edge
      end select
      go to 10
!
!
      end subroutine my_tests
