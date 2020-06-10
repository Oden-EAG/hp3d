#include "typedefs.h"
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
      use refinements
      use refinements_history
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
!  ...work space for solelm
      VTYPE ::  zdofH(MAXEQNH,MAXbrickH)
      VTYPE ::  zdofE(MAXEQNE,MAXbrickE)
      VTYPE ::  zdofV(MAXEQNV,MAXbrickV)
      VTYPE ::  zdofQ(MAXEQNQ,MAXbrickQ)
!
!  ...workspace for compute_error
      integer :: flag(NR_PHYSA), itag
!
!  ...workspace for nodcor
      integer :: kref,nrs
      real*8 :: xnod(3,MAXbrickH)
!
!  ...workspace for initiate_dof
      integer :: nel, nrH, nrE, nrQ
!
!------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...set MAX_ELEMS_REF in module refinements_history
      MAX_ELEMS_REF = NRELES
!
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
      write(*,*) 'test solelm...................................10'
      write(*,*) 'compute_error.................................11'
      write(*,*) 'initiate dof..................................12'
      write(*,*) 'initiate vertex coordinates...................13'
      write(*,*) 'execute initiate_dof..........................14'
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
      case(10)
        do nod=1,NRNODS
          call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
          select case(NODES(nod)%case)
          case(2)
            NODES(nod)%dof%zdofH(1,1:ndofH) = float(nod)
          case(1)
            NODES(nod)%dof%zdofH(1,1:ndofH) = float(nod)*10
          case(3)
            NODES(nod)%dof%zdofH(1,1:ndofH) = float(nod)
            NODES(nod)%dof%zdofH(2,1:ndofH) = float(nod)*10
          end select
        enddo
        write(*,*) 'my_tests: SET mdle'
        read(*,*) mdle
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
      case(11)
        flag=1
        call compute_error(flag,itag)
      case(12)
        do nod=1,NRNODS
          call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
          if (ndofH.ne.0) then
            NODES(nod)%dof%zdofH=float(nod)
          endif
          if (ndofE.ne.0) then
            NODES(nod)%dof%zdofE=float(nod)
          endif
          if (ndofV.ne.0) then
            NODES(nod)%dof%zdofV=float(nod)
          endif
          if (ndofQ.ne.0) then
            NODES(nod)%dof%zdofQ=float(nod)
          endif
        enddo
      case(13)
        write(*,*) 'DID U REMEMBER TO MODIFY break ?'
        write(*,*) 'SELECT kref'
        read(*,*) kref
        call break(1,kref)
        call nr_mdle_sons(NODES(1)%type,kref, nrs)
        do i=1,nrs
          mdle = son(1,i)
          call nodcor(mdle, xnod)
          write(*,7040) xnod(1:3,1:nvert(NODES(mdle)%type))
 7040     format(10(6x,3(3(e8.3,','),1x)/))
        enddo     
!
      case(14)
        do nel=1,NR_ELEMS_REF
          mdle = ELEMS_REF(nel)%mdle
          nrH = ubound(ELEMS_REF(nel)%xnod,2)
          nrE = ubound(ELEMS_REF(nel)%zdofE,2)
          nrV = ubound(ELEMS_REF(nel)%zdofV,2)
          nrQ = ubound(ELEMS_REF(nel)%zdofQ,2)
          call initiate_dof(mdle,ELEMS_REF(nel)%xnod,ELEMS_REF(nel)%zdofH,  &
                                 ELEMS_REF(nel)%zdofE,ELEMS_REF(nel)%zdofV, &
                                 ELEMS_REF(nel)%zdofQ,nrH,nrE,nrV,nrQ)
        enddo
        call deallocref
      end select
      go to 10
!
!
      end subroutine my_tests
