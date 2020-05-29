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
      common /cneig_edge/ iprint_neig_edge
!

!----------------------------------------------------------------------
!  
      integer :: iprint,iselect,iprint_recta_TraQua,iprint_curve_SegCir,iprint_neig_edge
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
      integer, dimension(maxn) :: neig, nedg_list, norient_list, nface_list
      integer                  :: nrv, ie
      integer, dimension(27)   :: nodesl, norientl
!
!  ...work space for random_refine
      integer   :: niter
      real*8    :: percent
!
!  ...work space for nodmod
      integer :: nod,newp,ndofH,ndofE,ndofV,ndofQ,oldp,ndofHn,ndofEn,ndofVn,ndofQn
!
!------------------------------------------------------------------------------------
!
      iprint=0
   10 write(*,*) 'my_tests: SELECT:'
      write(*,*) 'EXIT...........................................0'
      write(*,*) 'TEST recta_TraQua..............................1'
      write(*,*) 'TEST curve_SegCir..............................2'
      write(*,*) 'TEST find_neig.................................3'
      write(*,*) 'verify_neig....................................4'
      write(*,*) 'TEST neig_edge.................................5'
      write(*,*) 'verify_neig_edge...............................6'
      write(*,*) 'call result....................................7'
      write(*,*) 'call random_refine and verify neig routines....8'
      write(*,*) 'test nodmod....................................9'
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
        iprint_neig_edge=1
        write(*,*) 'input mdle, ie'
        read(*,*)  mdle,ie
        call elem_nodes(mdle, nodesl,norientl)
        nrv = nvert(NODES(mdle)%type)
        medge = nodesl(nrv+ie)
        write(*,*) 'my_tests: medge = ',medge
        call neig_edge(medge,maxn, nrneig,neig,nedg_list,norient_list,nface_list)
      case(6)
        call verify_neig_edge
      case(7)
        call result
     case(8)
        write(*,*) 'my_tests: SET percent AND niter'
        read(*,*) percent, niter
        write(*,7020) percent, niter
 7020   format(' my_tests: percent = ',f7.2,' niter = ',i3)
        call random_refine(percent,niter)
      case(9)
        write(*,*) 'my_tests: SET nod'
        read(*,*) nod
        write(*,*) 'NODES(nod)%type, NODES(nod)%order = ', &
                    NODES(nod)%type, NODES(nod)%order
        write(*,*) '          SET newp'
        read(*,*) newp
        call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
        do i=1,ndofH
          NODES(nod)%dof%zdofH(1,i) = float(i)
        enddo
        do i=1,ndofE
          NODES(nod)%dof%zdofE(1,i) = float(i)
        enddo
        do i=1,ndofV
          NODES(nod)%dof%zdofV(1,i) = float(i)
        enddo
        do i=1,ndofQ
          NODES(nod)%dof%zdofQ(1,i) = float(i)
        enddo
        oldp = NODES(nod)%order
        call nodmod(nod,newp)
        call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofHn,ndofEn,ndofVn,ndofQn)
        write(*,7030) NODES(nod)%dof%zdofH(1,1:ndofHn)
 7030   format(20f5.1)
        write(*,7030) NODES(nod)%dof%zdofE(1,1:ndofEn)
        write(*,7030) NODES(nod)%dof%zdofV(1,1:ndofVn)
        write(*,7030) NODES(nod)%dof%zdofQ(1,1:ndofQn)
        call nodmod(nod,oldp)
        write(*,7030) NODES(nod)%dof%zdofH(1,1:ndofH)
        write(*,7030) NODES(nod)%dof%zdofE(1,1:ndofE)
        write(*,7030) NODES(nod)%dof%zdofV(1,1:ndofV)
        write(*,7030) NODES(nod)%dof%zdofQ(1,1:ndofQ)
!
      end select
      go to 10
!
!
      end subroutine my_tests
