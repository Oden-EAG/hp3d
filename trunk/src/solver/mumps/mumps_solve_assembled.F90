SUBROUTINE mumps_solve_assembled(numrhs)
  ! ONLY SINGLE PHYSICS (1 H1 variable) AT THE MOMENT...
  USE data_structure3D
  USE element_data
  USE assembly
  USE control, ONLY: ISYM_FLAG
  USE MPI
  IMPLICIT NONE
#if C_MODE
  INCLUDE 'zmumps_struc.h'
  TYPE(zmumps_struc) :: mumps_par
#else
  INCLUDE 'dmumps_struc.h'
  TYPE(dmumps_struc) :: mumps_par
#endif
  INTEGER, INTENT(IN) :: numrhs
  INTEGER :: info, mdle, i_elem, nrdofs(NR_PHYSA), nrdofm, nrdofc, nodm(MAXNODM), &
       ndofmH(MAXNODM), ndofmE(MAXNODM), ndofmV(MAXNODM), ndofmQ(MAXNODM), &
       nrnodm, num_entries, i, i_node, i_global, j, i_local, k
  DOUBLE COMPLEX :: zvoid
  DOUBLE COMPLEX, ALLOCATABLE :: A_loc(:), b_loc(:)
  INTEGER, ALLOCATABLE :: ELTPTR(:), ELTVAR(:), node_global_index(:)

  NR_RHS = numrhs
  IF (ISYM_FLAG == 2) THEN
     PRINT *, "mumps_solve_assembled: running in UNSYMMETRIC mode"
  ELSE
     PRINT *, "mumps_solve_assembled: running in SYMMETRIC mode"
  ENDIF
!  CALL mpi_init(info)
  mumps_par%COMM = MPI_COMM_WORLD
  mumps_par%SYM = 0 ! unsymmetric
  mumps_par%PAR = 1
  mumps_par%JOB = -1
  CALL zmumps(mumps_par)
  mumps_par%ICNTL(22) = 1 ! out-of-core
  mumps_par%ICNTL(14) = 70
  ALLOCATE(ELTPTR(NRELES+1))
  ELTPTR(1) = 1
  MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
  ALLOCATE(NEXTRACT(MAXDOFM))
  ALLOCATE(IDBC(MAXDOFM))
  ALLOCATE(ZDOFD(MAXDOFM))
  ALLOCATE(MAXDOFS(NR_PHYSA))
  MAXDOFS = 0
  MAXDOFM = 0
  MAXDOFC = 0
  num_entries = 0
  mdle = 0
  DO i_elem = 1, NRELES
     CALL nelcon(mdle, mdle)
     CALL celem(mdle, 1, &
          nrdofs, nrdofm, nrdofc, &
          nodm, ndofmH, ndofmE, ndofmV, ndofmQ, nrnodm, &
          zvoid, zvoid)
     !PRINT *, "nrdofm, nrdofc, SUM(ndofmH(1:nrnodm)) =", nrdofm, nrdofc, SUM(ndofmH(1:nrnodm))
     ELTPTR(i_elem+1) = ELTPTR(i_elem) + nrdofc
     !ELTPTR(i_elem+1) = ELTPTR(i_elem) + nrdofm
     IF (ISYM_FLAG == 2) THEN
        num_entries = num_entries + nrdofc**2
        !num_entries = num_entries + nrdofm**2
     ELSE
        num_entries = num_entries + (nrdofc+1)*nrdofc/2
        !num_entries = num_entries + (nrdofm+1)*nrdofm/2
     ENDIF
     DO i = 1, NR_PHYSA
        MAXDOFS(i) = MAX(MAXDOFS(i), nrdofs(i))
     ENDDO
     MAXDOFM = MAX(MAXDOFM, nrdofm)
     MAXDOFC = MAX(MAXDOFC, nrdofc)
  ENDDO

  mumps_par%NZ = num_entries
  PRINT *, "mumps_solve_assembled: number of non-zero entries =", num_entries

  ALLOCATE(ELTVAR(ELTPTR(NRELES+1)-1))
  ALLOCATE(node_global_index(NRNODS))
  i_global = 0
  node_global_index = 0
  mdle = 0
  DO i_elem = 1, NRELES
     CALL nelcon(mdle, mdle)
     CALL celem(mdle, 1, &
          nrdofs, nrdofm, nrdofc, &
          nodm, ndofmH, ndofmE, ndofmV, ndofmQ, nrnodm, &
          zvoid, zvoid)
     i_local = 0
     DO i_node = nrnodm, 1, -1
        IF (ndofmH(i_node) > 0) THEN
           IF (node_global_index(nodm(i_node)) == 0) THEN
              node_global_index(nodm(i_node)) = i_global + 1
              FORALL (i = 1:ndofmH(i_node))
                 ELTVAR(ELTPTR(i_elem)-1+i_local+i) = node_global_index(nodm(i_node))-1+i
              END FORALL
              i_local = i_local + ndofmH(i_node)
              i_global = i_global + ndofmH(i_node)
           ELSE
              FORALL (i = 1:ndofmH(i_node))
                 ELTVAR(ELTPTR(i_elem)-1+i_local+i) = node_global_index(nodm(i_node))-1+i
              END FORALL
              i_local = i_local + ndofmH(i_node)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  mumps_par%N = i_global
  PRINT *, "mumps_solve_assembled: N =", mumps_par%N
  ALLOCATE(mumps_par%IRN(mumps_par%NZ))
  ALLOCATE(mumps_par%JCN(mumps_par%NZ))
  ALLOCATE(mumps_par%A(mumps_par%NZ))
  ALLOCATE(mumps_par%RHS(mumps_par%N))
  mumps_par%RHS = 0
  ALLOCATE(A_loc(MAXDOFC*MAXDOFC))
  !ALLOCATE(A_loc(MAXDOFM, MAXDOFM))
  ALLOCATE(b_loc(MAXDOFC))
  !ALLOCATE(b_loc(MAXDOFM))

  ALLOCATE(BLOC(NR_PHYSA))
  ALLOCATE(AAUX(NR_PHYSA))
  ALLOCATE(ALOC(NR_PHYSA, NR_PHYSA))
  DO i = 1, NR_PHYSA
     BLOC(i)%nrow = MAXDOFS(i)
     BLOC(i)%ncol = numrhs
     ALLOCATE(BLOC(i)%array(MAXDOFS(i), numrhs))
     DO j = 1, NR_PHYSA
        ALOC(i,j)%nrow = MAXDOFS(i)
        ALOC(i,j)%ncol = MAXDOFS(j)
        ALLOCATE(ALOC(i,j)%array(MAXDOFS(i), MAXDOFS(j)))
     ENDDO
     AAUX(i)%nrow = MAXDOFM
     AAUX(i)%ncol = MAXDOFS(i)
     ALLOCATE(AAUX(i)%array(MAXDOFM, MAXDOFS(i)))
  ENDDO
  ALLOCATE(ZBMOD(MAXDOFM, numrhs))
  ALLOCATE(ZAMOD(MAXDOFM, MAXDOFM))

  k = 0
  mdle = 0
  DO i_elem = 1, NRELES
     IF (MOD(i_elem, MAX(1,NRELES/10)) == 0) THEN
        PRINT *, "mumps_solve_assembled: computing element matrix", i_elem, "of", NRELES
     ENDIF
     CALL nelcon(mdle, mdle)
     CALL celem(mdle, 2, &
          nrdofs, nrdofm, nrdofc, &
          nodm, ndofmH, ndofmE, ndofmV, ndofmQ, nrnodm, &
          b_loc, A_loc)
     FORALL (i = 1:nrdofc)
        mumps_par%RHS(ELTVAR(ELTPTR(i_elem)-1+i)) = &
        mumps_par%RHS(ELTVAR(ELTPTR(i_elem)-1+i)) + b_loc(i)
     END FORALL
     DO j = 1, nrdofc
        DO i = 1, nrdofc
           k = k + 1
           !mumps_par%A(k) = A_loc(i,j) ! <-- THIS IS STUPID
           mumps_par%A(k) = A_loc(i + (j-1)*nrdofc) ! <-- THIS IS STUPID
           mumps_par%IRN(k) = ELTVAR(ELTPTR(i_elem)-1+i)
           mumps_par%JCN(k) = ELTVAR(ELTPTR(i_elem)-1+j)
        ENDDO
     ENDDO
  ENDDO

  DEALLOCATE(ELTPTR)
  DEALLOCATE(ELTVAR)
  DEALLOCATE(A_loc)
  DEALLOCATE(b_loc)
  DO i = 1, NR_PHYSA
     DEALLOCATE(BLOC(i)%array)
     DO j = 1, NR_PHYSA
        DEALLOCATE(ALOC(i,j)%array)
     ENDDO
     DEALLOCATE(AAUX(i)%array)
  ENDDO
  DEALLOCATE(BLOC)
  DEALLOCATE(AAUX)
  DEALLOCATE(ALOC)
  DEALLOCATE(ZBMOD)
  DEALLOCATE(ZAMOD)

  !mumps_par%WRITE_PROBLEM = "/workspace/jzitelli/hp3d/mumps_matrix"

  mumps_par%JOB = 1
#if C_MODE
  CALL zmumps(mumps_par)
#else
  CALL dmumps(mumps_par)
#endif

  mumps_par%JOB = 2
#if C_MODE
  CALL zmumps(mumps_par)
#else
  CALL dmumps(mumps_par)
#endif

  mumps_par%JOB = 3
#if C_MODE
  CALL zmumps(mumps_par)
#else
  CALL dmumps(mumps_par)
#endif

  mumps_par%JOB = -2
#if C_MODE
  CALL zmumps(mumps_par)
#else
  CALL dmumps(mumps_par)
#endif

  DEALLOCATE(mumps_par%A)
  DEALLOCATE(mumps_par%IRN)
  DEALLOCATE(mumps_par%JCN)

  PRINT *, "mumps_solve_assembled: storing results in data structure..."
  mdle = 0
  DO i_elem = 1, NRELES
     CALL nelcon(mdle, mdle)
     CALL celem(mdle, 1, &
          nrdofs, nrdofm, nrdofc, &
          nodm, ndofmH, ndofmE, ndofmV, ndofmQ, nrnodm, &
          zvoid, zvoid)
     DO i_node = nrnodm, 1, -1
        IF (ndofmH(i_node) > 0) THEN
           IF (ASSOCIATED(NODES(nodm(i_node))%zdofH)) THEN
              FORALL (i = 1:ndofmH(i_node))
                 NODES(nodm(i_node))%zdofH(1,i) = mumps_par%RHS(node_global_index(nodm(i_node))-1+i)
              END FORALL
              i = NODES(nodm(i_node))%father
              IF (i > 0) THEN
                 i = NODES(i)%father
                 IF (i > 0) THEN
                    PRINT *, "ASSOCIATED, NODES(grandfather)%type =", NODES(i)%type
                 ENDIF
              ENDIF
           ELSE
              PRINT *, "NOT ASSOCIATED!"
              PRINT *, "  node =", nodm(i_node)
              PRINT *, "  NODES(nodm(i_node))%type =", NODES(nodm(i_node))%type
              PRINT *, "  NODES(nodm(i_node))%act =", NODES(nodm(i_node))%act
              PRINT *, "  NODES(nodm(i_node))%father =", NODES(nodm(i_node))%father
              i = NODES(nodm(i_node))%father
              IF (i > 0) THEN
                 PRINT *, "    NODES(father)%type =", NODES(i)%type
                 i = NODES(i)%father
                 IF (i > 0) THEN
                    PRINT *, "      NODES(grandfather)%type =", NODES(i)%type
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  DEALLOCATE(mumps_par%RHS)
  DEALLOCATE(node_global_index)
  DEALLOCATE(IDBC)
  DEALLOCATE(ZDOFD)
  DEALLOCATE(NEXTRACT)
  DEALLOCATE(MAXDOFS)

END SUBROUTINE mumps_solve_assembled
