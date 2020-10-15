
    subroutine all_tests()

        call test1

        write(*,*)
        call test2

        write(*,*)
        call test3

        write(*,*)
        call test4

        write(*,*)
        call test5

        write(*,*)
        call test6

    end subroutine

    ! https://stackoverflow.com/questions/20006253/using-metis-libraries-in-fortran-code-the-basics
    ! http://glaros.dtc.umn.edu/gkhome/node/852
    subroutine test1()

        use metis_interface, only: METIS_SetDefaultOptions, METIS_PartMeshNodal, METIS_PartMeshDual, &
            METIS_MeshToNodal, METIS_Free, METIS_NOPTIONS
        use metis_enum, only: METIS_OK, METIS_OPTION_NUMBERING
        use, intrinsic :: iso_c_binding, only : c_long, c_ptr, c_f_pointer

        integer(c_long), parameter :: ne0 = 3 ! number of elements
        integer(c_long), parameter :: nn0 = 8 ! number of nodes
        integer(c_long), parameter :: npel0 = 4 ! nodes per element

        integer(c_long):: eptr0(ne0+1)
        integer(c_long):: eind0(ne0*npel0)
        integer(c_long):: epart0(ne0), npart0(nn0)
        integer(c_long):: options0(0:METIS_NOPTIONS-1)
        integer(c_long):: ios, objval0

        ! integer(c_long),bind(C) :: ne,nn
        ! integer(c_long), parameter :: npel = 4 ! nodes per element

        ! integer(c_long),bind(C),allocatable :: eptr(:),eind(:),epart(:),npart(:),options(:)

        type(c_ptr) :: c_xadj, c_adjncy
        integer(c_long), pointer :: xadj(:) => null(), adjncy(:) => null()

        write(*,'(A)') "TEST 1"

        ! 0---1---4---6
        ! | 0 | 1 | 2 |
        ! 3---2---5---7
        ! allocate(ne,nn)
        ! ne=3
        ! nn=8
        ! allocate(eptr(ne+1))
        ! allocate(eind(ne*npel))
        ! allocate(epart(ne),npart(nn))
        eptr0 = [0,4,8,12]
        eind0 = [0,1,2,3,1,4,5,2,4,6,7,5]
        ! allocate(options(0:METIS_NOPTIONS-1))
        options0 = -1
        ! ios = METIS_SetDefaultOptions(options)
        ! if (ios /= METIS_OK) then
        !     write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
        !     error stop 1
        ! end if
        options0(METIS_OPTION_NUMBERING) = 0 ! C-style numbering

        ios = METIS_PartMeshNodal(ne=ne0,nn=nn0,eptr=eptr0,eind=eind0,nparts=2, options=options0(0:39),&
            objval=objval0,epart=epart0,npart=npart0)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartMeshNodal failed with error: ", ios
            error stop 1
        end if
        write(*,*) "NODAL"    
        write(*,'(A,I0)') "objval = ", objval0
        write(*,'(A,*(I1,:,1X))') "epart = ", epart0
        write(*,'(A,*(I1,:,1X))') "npart = ", npart0

        options0(METIS_OPTION_NUMBERING) = 0 ! C-style numbering

        ios = METIS_PartMeshDual(ne=ne0,nn=nn0,eptr=eptr0,eind=eind0,ncommon=2,nparts=2, options=options0(0:39),&
            objval=objval0,epart=epart0,npart=npart0)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartMeshDual failed with error: ", ios
            error stop 1
        end if
        write(*,*) "DUAL" 
        write(*,'(A,I0)') "objval = ", objval0
        write(*,'(A,*(I1,:,1X))') "epart = ", epart0
        write(*,'(A,*(I1,:,1X))') "npart = ", npart0

        ios = METIS_MeshToNodal(ne=ne0,nn=nn0,eptr=eptr0,eind=eind0,numflag=0,xadj=c_xadj,adjncy=c_adjncy)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_MeshToNodal failed with error: ", ios
            error stop 1
        end if

        call c_f_pointer(c_xadj,xadj,shape=[nn0+1]) ! size of adjacency list is one more than number of nodes
        call c_f_pointer(c_adjncy,adjncy,shape=[xadj(nn0+1)]) ! size of edge list is in the last element of xadj

        write(*,'(A,*(I0,:,1X))') "xadj = ", xadj
        write(*,'(A,*(I1,:,1X))') "adjncy = ", adjncy

        ! call write_graph("test1.graph",xadj,adjncy,1)

        ios = METIS_Free(c_xadj); xadj => null()
        if (ios /= METIS_OK) then
            write(*,*) "METIS_Free failed with error: ", ios
            error stop 1
        end if
        ios = METIS_Free(c_adjncy); adjncy => null()
        if (ios /= METIS_OK) then
            write(*,*) "METIS_Free failed with error: ", ios
            error stop 1
        end if

    end subroutine

    ! https://stackoverflow.com/questions/8155160/metis-with-fortran
    ! http://glaros.dtc.umn.edu/gkhome/node/799
    subroutine test2()

        use metis_interface, only: METIS_PartGraphRecursive
        use metis_enum, only: METIS_OK
        use iso_c_binding, only: c_long

        integer(c_long), parameter :: nvtxs = 15     ! number of vertices
        integer(c_long), parameter :: nedgs = 22     ! number of edges

        integer(c_long) :: xadj(nvtxs+1),adjncy(2*nedgs) ! adjacency arrays
        integer(c_long) :: part(nvtxs)                   ! partiotion vector
        integer(c_long) :: objval, ios

        write(*,'(A)') "TEST 2"

        !  0---1---2---3---4
        !  |   |   |   |   |
        !  5---6---7---8---9
        !  |   |   |   |   |
        ! 10--11--12--13--14

        ! Note that we are using C-style numbering!

        xadj = [0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44]
        adjncy = [1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13]

        ios = METIS_PartGraphRecursive(nvtxs,ncon=1,xadj=xadj,adjncy=adjncy,nparts=2,objval=objval,part=part)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartGraphRecursive failed with error: ", ios
            error stop 1
        end if

        write(*,'(A,I0)') "objval = ", objval
        write(*,'(A,*(I1,:,1X))') "part = ", part

    end subroutine


!>  Example of computing a fill-reducing ordering of a sparse matrix
!
!   Source: http://comp.lang.fortran.narkive.com/uFmDM7Bo/how-to-call-a-metis-subroutine-from-my-fortran-code
!
    subroutine test3()

        use metis_interface, only: METIS_SetDefaultOptions, METIS_NodeND, METIS_NOPTIONS
        use metis_enum, only: METIS_OK, METIS_OPTION_NUMBERING
        use iso_c_binding, only: c_long

        integer(c_long), parameter :: n = 15 ! number of vertices
        integer(c_long), parameter :: m = 22 ! number of edges

        integer(c_long) :: xadj(n+1), adjncy(2*m) ! graph adjacency structure
        integer(c_long) :: perm(n), iperm(n) ! fill-reducing permutation and inverse permutation

        integer(c_long) :: options(0:METIS_NOPTIONS-1), ios

        write(*,'(A)') "TEST 3"

        !  1---2---3---4---5
        !  |   |   |   |   |
        !  6---7---8---9---10
        !  |   |   |   |   |
        !  11--12--13--14--15

        xadj = [1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45]
        adjncy = [2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6, &
                  8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,7,11,13, &
                  8,12,14,9,13,15,10,14]

        ios = METIS_SetDefaultOptions(options)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
            error stop 1
        end if
        options(METIS_OPTION_NUMBERING) = 1 ! Fortran-style numbering

        ios = METIS_NodeND(n,xadj,adjncy,options=options,perm=perm,iperm=iperm)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_NodeND failed with error: ", ios
            error stop 1
        end if

        write(*,'(A,*(I2,:,1X))') "perm  = ", perm
        write(*,'(A,*(I2,:,1X))') "iperm = ", iperm

    end subroutine


!>  Example of partitioning a mesh composed of 4 quadrilaterals and 9 nodes
!   based upon its nodal graph.
!
!   Source: https://www.cfd-online.com/Forums/main/112366-using-metis-functions-fortran.html#post404734
!
    subroutine test4

        use metis_interface, only: METIS_SetDefaultOptions, METIS_PartMeshDual, METIS_PartMeshNodal, METIS_NOPTIONS
        use metis_enum, only: METIS_OPTION_NUMBERING, METIS_OPTION_CONTIG, METIS_OK
        use iso_c_binding, only: c_long

        integer(c_long), parameter :: ne = 4!651    ! number of elements
        integer(c_long), parameter :: nn = 9!163    ! number of nodes

        integer(c_long) :: eptr(ne+1), eind(4*ne)   ! arrays storing mesh structure
        integer(c_long) :: epart(ne), npart(nn)     ! element and node partition vectors

        integer(c_long) :: opts(0:METIS_NOPTIONS-1), ios, objval , i

        write(*,'(A)') "TEST 4"

         ! 1---2---5
         ! | 1 | 2 |
         ! 4---3---6
         ! | 4 | 3 |
         ! 9---8---7

        eptr = [1,5,9,13,17]
        eind = [1,2,3,4,&
                2,5,6,3,&
                3,6,7,8,&
                4,3,8,9]
        ! open(unit=1,file='eptr.txt', &
        !       form='formatted',access='sequential',status='unknown')
        ! do i=1,ne+1
        !     read(1,*) eptr(i)
        ! enddo
        ! close(1)
        ! open(unit=2,file='eind.txt', &
        !       form='formatted',access='sequential',status='unknown')
        ! do i=1,4*ne
        !     read(2,*) eind(i)
        ! enddo
        ! close(2)

        ios = METIS_SetDefaultOptions(opts)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
            error stop 1
        end if
        opts(METIS_OPTION_NUMBERING) = 1    ! Fortran-style numbering
        opts(METIS_OPTION_CONTIG) = 1       ! Force contigous partitions
        ! call print_metis_options(opts)
        
        ios = METIS_PartMeshDual(ne,nn,eptr,eind,ncommon=2,nparts=2,options=opts, &
                objval=objval,epart=epart,npart=npart)
        ! ios = METIS_PartMeshNodal(ne,nn,eptr,eind,nparts=4,options=opts, &
        !         objval=objval,epart=epart,npart=npart)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartMeshNodal failed with error: ", ios
            error stop 1
        end if

        write(*,'(A,I0)') "objval = ", objval
        write(*,'(A,*(I1,:,1X))') "epart = ", epart
        write(*,'(A,*(I1,:,1X))') "npart = ", npart

 !        open(unit=3,file='output/epart', &
 !              form='formatted',access='sequential',status='unknown')
 !        write(3,1111) epart
 ! 1111   format(i3.3)
 !        close(3)
 !        open(unit=4,file='output/npart', &
 !              form='formatted',access='sequential',status='unknown')
 !        write(4,1111) npart
 !        close(4)

    end subroutine

!>  Partition a graph specified by a matrix
!
!   Source: http://people.eecs.berkeley.edu/~demmel/cs267/lecture18/lecture18.html
!
    subroutine test5

        use metis_interface, only: METIS_PartGraphKway, METIS_SetDefaultOptions, METIS_NOPTIONS
        use metis_enum, only: METIS_OPTION_NUMBERING, METIS_OK
        use iso_c_binding, only: c_long

        integer(c_long), parameter :: npart = 3  ! number of partitions
        integer(c_long), parameter :: n = 8      ! number of nodes
        integer(c_long), parameter :: m = 10     ! number of edges

        integer(c_long) :: xadj(n+1), adjncy(2*m) ! graph adjacency structure
        integer(c_long) :: perm(n), iperm(n) ! fill-reducing permutation and inverse permutation
        integer(c_long) :: part(n)
        integer(c_long) :: options(0:METIS_NOPTIONS-1), ios, objval

        write(*,'(A)') "TEST 5"

        !    1  2  3  4  5  6  7  8
        !  | a  a  a     a          | 1
        !  | a  a     a  a          | 2
        !  | b     b        b       | 3
        !  |    b     b     b       | 4
        !  | a  a           a     a | 5
        !  |       b  b  b          | 6
        !  |                   c  c | 7
        !  |             c     c  c | 8

        xadj = [1,4,7,9,11,15,18,19,21]
        adjncy = [2,3,5,&
                  1,4,5,&
                  1,6,&
                  2,6,&
                  1,2,6,8,&
                  3,4,5,&
                  8,&
                  5,7]

        ios = METIS_SetDefaultOptions(options)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
            error stop 1
        end if
        options(METIS_OPTION_NUMBERING) = 1 ! Fortran-style numbering

        ios = METIS_PartGraphKway(n,ncon=1,xadj=xadj,adjncy=adjncy,nparts=npart,options=options,objval=objval,part=part)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartGraphKway failed with error: ", ios
            error stop 1
        end if

        write(*,'(A,I0)') "objval = ", objval
        write(*,'(A,*(I1,:,1X))') "part = ", part

    end subroutine

!>  Partitioning a graph with non-equal number of edges per node.
!
    subroutine test6

        use metis_interface, only: METIS_PartGraphRecursive, METIS_SetDefaultOptions, METIS_NOPTIONS
        use metis_enum, only: METIS_OPTION_NUMBERING, METIS_OK
        use iso_c_binding, only: c_long

        integer(c_long) :: n, m
        integer(c_long), allocatable :: xadj(:), adjncy(:), part(:)
        integer(c_long) :: options(0:METIS_NOPTIONS-1), ios, ncon, objval

        write(*,'(A)') "TEST 6"

        !   1---5
        !   |\ / \ 
        !   | 3   6      
        !   |/ \ / \    
        !   2---4---7 

        xadj = [1,4,7,11,15,18,21,23]
        adjncy = [5,3,2,1,3,4,5,4,2,1,2,3,6,7,1,3,6,5,4,7,6,4]

        n = size(xadj) - 1
        m = size(adjncy)/2
        write(*,'(A,I0)') "n = ", n
        write(*,'(A,I0)') "m = ", m

        ios = METIS_SetDefaultOptions(options)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_SetDefaultOptions failed with error: ", ios
            error stop 1
        end if
        options(METIS_OPTION_NUMBERING) = 1     ! Fortran-style numbering

        ncon = 1
        allocate(part(n))
        ios = METIS_PartGraphRecursive(n,ncon,xadj,adjncy,&
                    nparts=2,objval=objval,part=part,options=options)
        if (ios /= METIS_OK) then
            write(*,*) "METIS_PartGraphKway failed with error: ", ios
            error stop 1
        end if

        write(*,'(A,I0)') "objval = ", objval
        write(*,'(A,*(I1,:,1X))') "part = ", part

    end subroutine

    subroutine test_load_and_write()

        ! type(MetisGraph) :: mgraph

        ! call load_graph(mgraph,"../graphs/4elt.graph",numflag=0)

        ! print *, mgraph%nvtxs
        ! print *, mgraph%nedgs
        ! print *, mgraph%xadj
        ! print *, mgraph%adjncy
        ! call write_graph("print.graph",)
    end subroutine
