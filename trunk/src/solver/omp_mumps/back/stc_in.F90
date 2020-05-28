
  subroutine stc_in(iel,nc,nb,xc,xb)
!
   use m_assembly, only : CLOC, STORE_STC
!
   IMPLICIT NONE
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
   integer           :: nc,nb, iel
#if C_MODE
   complex*16        :: xc(nc), xb(nb), Abc(nb,nc)
#else
   real*8            :: xc(nc), xb(nb), Abc(nb,nc)
#endif
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!

   if (STORE_STC) then
      xb = CLOC(iel)%vect
      Abc = CLOC(iel)%array
      deallocate(CLOC(iel)%vect)
      deallocate(CLOC(iel)%array)
#if C_MODE
      call ZGEMV('N',nb,nc,(-1.0d0,0.0d0),Abc,nb,xc,1,(1.0d0,0.0d0),xb,1)
#else
      call DGEMV('N',nb,nc,-1.0d0,Abc,nb,xc,1,1.0d0,xb,1)
#endif
   else
      call compute_schur_complement(iel,nc,nb,xc,xb)
   endif
!
!
   end subroutine stc_in



   subroutine compute_schur_complement(iel,nc,nb,xc,xb)
!
   use m_assembly, only: CLOC, ZLOAD, ZTEMP
   use assembly,   only: NR_RHS,NR_PHYSA, MAXNODM, MAXDOFS, MAXDOFM, &
                         BLOC, AAUX,ALOC, ZBMOD, ZAMOD
!
   IMPLICIT NONE

   integer           :: iel, mdle,i,j, nrdofm, nrdofc, nrnodm
   integer           :: nc,nb,n, nHc, nEc, nVc, nHb, nEb, nVb, nQ
#if C_MODE
   complex*16        :: xc(nc), xb(nb)
#else
   real*8            :: xc(nc), xb(nb)
#endif
!..workspace for celem
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!..integer counters

!..work space for celem
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM), &
              ndofmE(MAXNODM),ndofmV(MAXNODM),ndofmQ(MAXNODM)



!...find the mdle number
   mdle = 0
   do i = 1, iel
      call nelcon(mdle,mdle)
   enddo

!..get the matrix and the load vector from celem
   allocate(BLOC(NR_PHYSA))
   allocate(AAUX(NR_PHYSA))
   allocate(ALOC(NR_PHYSA,NR_PHYSA))
   do i=1,NR_PHYSA
      BLOC(i)%nrow = MAXDOFS(i)
      BLOC(i)%ncol = NR_RHS
      allocate(BLOC(i)%array(MAXDOFS(i),NR_RHS))
      do j=1,NR_PHYSA
         ALOC(i,j)%nrow = MAXDOFS(i)
         ALOC(i,j)%ncol = MAXDOFS(j)
         allocate(ALOC(i,j)%array(MAXDOFS(i),MAXDOFS(j)))
      enddo
      AAUX(i)%nrow = MAXDOFM
      AAUX(i)%ncol = MAXDOFS(i)
      allocate(AAUX(i)%array(MAXDOFM,MAXDOFS(i)))
   enddo
   allocate(ZBMOD(MAXDOFM,NR_RHS))
   allocate(ZAMOD(MAXDOFM,MAXDOFM))
   allocate(ZLOAD(MAXDOFM),ZTEMP(MAXDOFM**2))

   call celem(mdle,2,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
              ndofmE,ndofmV,ndofmQ,nrnodm, ZLOAD,ZTEMP)

   deallocate(BLOC,AAUX,ALOC,ZBMOD,ZAMOD)


   nHc =  CLOC(iel)%nHc
   nEc =  CLOC(iel)%nEc
   nVc =  CLOC(iel)%nVc
!
   nHb =  CLOC(iel)%nHb
   nEb =  CLOC(iel)%nEb
   nVb =  CLOC(iel)%nVb
   nQ =   CLOC(iel)%nQ
!
   nc = nHc + nEc + nVc
   nb = nHb + nEb + nVb + nQ
   n = nc + nb
!
   call stc_herm_in(iel,ZTEMP(1:n**2),ZLOAD(1:n),n,nc,nb,nHc,nHb,nEc,nEb,nVc,nVb,nQ,xc,xb)
!
   deallocate(ZLOAD,ZTEMP)
!
!
   end subroutine compute_schur_complement


! -----------------------------------------------------------------------
!
!    routine name       - stc_herm
!
! -----------------------------------------------------------------------
!
!    latest revision    - July 17
!
!    purpose            - routine performs static condensation for the
!                         hermitian case
!                         (eliminates, H^1,H(curl),H(div) and L^2 bubbles)
!
!    arguments :
!      in:
!              A        - stiffness matrix
!              b        - load vector
!              n        - dimension of matrix A
!              nc       - total number of condensed dof
!              nb       - total number of bubbles
!              nHc      - number of condensed H1 dof
!              nHb      - number of bubble H1 dof
!              nEc      - number of condensed H(curl) dof
!              nEb      - number of bubble H(curl) dof
!              nVc      - number of condensed H(div) dof
!              nVb      - number of bubble H(div) dof
!              nQ       - number of bubble L^2 dof
!      out:
!              Ac       - the resulting Schur complement , condensed
!              bc       - modified right-hand side
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
   subroutine stc_herm_in(iel,Atemp,b,n,nc,nb,nHc,nHb,nEc,nEb,nVc,nVb,nQ,xc,xb)
!
   use m_assembly, only : CLOC
!
   IMPLICIT NONE
!
!----------------------------------------------------------------------
!
   integer           :: n, nc,nb,nHc,nHb,nEc,nEb,nVc,nVb,nQ, k1, k2,k
#if C_MODE
   complex*16        :: Atemp(n**2), b(n), Ac(nc**2), bc(nc), A(n,n)
   complex*16        :: Acc(nc,nc), Acb(nc,nb),Abc(nb,nc) , Abb(nb,nb)
   complex*16        :: AP(nb*(nb+1)/2) ,bb(nb), xc(nc), xb(nb)
#else
   real*8            :: Atemp(n**2), b(n), Ac(nc**2), bc(nc), A(n,n)
   real*8            :: Acc(nc,nc), Acb(nc,nb),Abc(nb,nc) , Abb(nb,nb)
   real*8            :: AP(nb*(nb+1)/2) ,bb(nb) ,xc(nc), xb(nb)
#endif
!
   integer           :: nk, nH,nE,nV,i,j,iel
   character*1       :: uplo
   integer           :: info,info1,info2
!..function for vector storage of a symmetric matrix
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!..set up the indices.
!..check if the indices are compatible
   if(nHc+nEc+nVc .ne. nc) then
      write(*,*) 'stc: incompatible indices for the condensed system '
      write(*,*) 'iel, nHc+nEc+nVc, nc = ', iel, nHc+nEc+nVc, nc
      stop 1
   endif
   if(nHb+nEb+nVb+nQ .ne. nb) then
      write(*,*) 'stc: incompatible indices for the bubbles'
      stop 2
   endif
   if(nc+nb .ne. n) then
      write(*,*) 'stc: incompatible indices for the total system'
      stop 3
   endif

!..total dof for each space
   nH = nHc + nHb
   nE = nEc + nEb
   nV = nVc + nVb

!..extract matrices
!..change format of the matrix(for now)

   do i = 1, n
      do j = 1, n
         k = (i-1)*n + j
         A(i,j) = Atemp(k)
      enddo
   enddo
!
!
!..first Acc and bc
!..H1 dof
   Acc(1:nHc,1:nHc)                 = A(nHb+1:nH,nHb+1:nH)
   Acc(1:nHc,nHc+1:nHc+nEc)         = A(nHb+1:nH,nH+nEb+1:nH+nE)
   Acc(1:nHc,nHc+nEc+1:nHc+nEc+nVc) = A(nHb+1:nH,nH+nE+nVb+1:nH+nE+nV)
   bc(1:nHc) = b(nHb+1:nH)
!..H(curl) dof
   Acc(nHc+1:nHc+nEc,1:nHc)                 = A(nH+nEb+1:nH+nE,nHb+1:nH)
   Acc(nHc+1:nHc+nEc,nHc+1:nHc+nEc)         = A(nH+nEb+1:nH+nE,nH+nEb+1:nH+nE)
   Acc(nHc+1:nHc+nEc,nHc+nEc+1:nHc+nEc+nVc) = A(nH+nEb+1:nH+nE,nH+nE+nVb+1:nH+nE+nV)
   bc(nHc+1:nHc+nEc) = b(nH+nEb+1:nH+nE)
!..H(div) dof
   Acc(nHc+nEc+1:nHc+nEc+nVc,1:nHc)         = A(nH+nE+nVb+1:nH+nE+nV,nHb+1:nH)
   Acc(nHc+nEc+1:nHc+nEc+nVc,nHc+1:nHc+nEc) = A(nH+nE+nVb+1:nH+nE+nV,nH+nEb+1:nH+nE)
   Acc(nHc+nEc+1:nHc+nEc+nVc,nHc+nEc+1:nHc+nEc+nVc) =    &
                                              A(nH+nE+nVb+1:nH+nE+nV,nH+nE+nVb+1:nH+nE+nV)
   bc(nHc+nEc+1:nHc+nEc+nVc) = b(nH+nE+nVb+1:nH+nE+nV)
!
!..Acb
!..H1 dof
   Acb(1:nHc,1:nHb) = A(nHb+1:nH,1:nHb)
   Acb(1:nHc,nHb+1:nHb+nEb) = A(nHb+1:nH,nH+1:nH+nEb)
   Acb(1:nHc,nHb+nEb+1:nHb+nEb+nVb) = A(nHb+1:nH,nH+nE+1:nH+nE+nVb)
   Acb(1:nHc,nHb+nEb+nVb+1:nHb+nEb+nVb+nQ) = A(nHb+1:nH,nH+nE+nV+1:nH+nE+nV+nQ)
!..H(curl) dof
   Acb(nHc+1:nHc+nEc,1:nHb) = A(nH+nEb+1:nH+nE,1:nHb)
   Acb(nHc+1:nHc+nEc,nHb+1:nHb+nEb) = A(nH+nEb+1:nH+nE,nH+1:nH+nEb)
   Acb(nHc+1:nHc+nEc,nHb+nEb+1:nHb+nEb+nVb) = A(nH+nEb+1:nH+nE,nH+nE+1:nH+nE+nVb)
   Acb(nHc+1:nHc+nEc,nHb+nEb+nVb+1:nHb+nEb+nVb+nQ) = A(nH+nEb+1:nH+nE,nH+nE+nV+1:nH+nE+nV+nQ)
!..H(div) dof

   Acb(nHc+nEc+1:nHc+nEc+nVc,1:nHb) = A(nH+nE+nVb+1:nH+nE+nV,1:nHb)
   Acb(nHc+nEc+1:nHc+nEc+nVc,nHb+1:nHb+nEb) = A(nH+nE+nVb+1:nH+nE+nV,nH+1:nH+nEb)
   Acb(nHc+nEc+1:nHc+nEc+nVc,nHb+nEb+1:nHb+nEb+nVb) =    &
                                        A(nH+nE+nVb+1:nH+nE+nV,nH+nE+1:nH+nE+nVb)
   Acb(nHc+nEc+1:nHc+nEc+nVc,nHb+nEb+nVb+1:nHb+nEb+nVb+nQ) = &
                                      A(nH+nE+nVb+1:nH+nE+nV,nH+nE+nV+1:nH+nE+nV+nQ)
!..Abc
#if C_MODE
   Abc = conjg(transpose(Acb))
#else
   Abc = transpose(Acb)
#endif
!
!..Abb
!..H1 dof
   Abb(1:nHb,1:nHb) = A(1:nHb,1:nHb)
   Abb(1:nHb,nHb+1:nHb+nEb) = A(1:nHb,nH+1:nH+nEb)
   Abb(1:nHb,nHb+nEb+1:nHb+nEb+nVb) = A(1:nHb,nH+nE+1:nH+nE+nVb)
   Abb(1:nHb,nHb+nEb+nVb+1:nb) = A(1:nHb,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(1:nHb) = b(1:nHb)
!..H(curl) dof
   Abb(nHb+1:nHb+nEb,1:nHb) = A(nH+1:nH+nEb,1:nHb)
   Abb(nHb+1:nHb+nEb,nHb+1:nHb+nEb) = A(nH+1:nH+nEb,nH+1:nH+nEb)
   Abb(nHb+1:nHb+nEb,nHb+nEb+1:nHb+nEb+nVb) = A(nH+1:nH+nEb,nH+nE+1:nH+nE+nVb)
   Abb(nHb+1:nHb+nEb,nHb+nEb+nVb+1:nb) = A(nH+1:nH+nEb,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(nHb+1:nHb+nEb) = b(nH+1:nH+nEb)
!..H(div) dof
   Abb(nHb+nEb+1:nHb+nEb+nVb,1:nHb) = A(nH+nE+1:nH+nE+nVb,1:nHb)
   Abb(nHb+nEb+1:nHb+nEb+nVb,nHb+1:nHb+nEb) = A(nH+nE+1:nH+nE+nVb,nH+1:nH+nEb)
   Abb(nHb+nEb+1:nHb+nEb+nVb,nHb+nEb+1:nHb+nEb+nVb) =         &
                                       A(nH+nE+1:nH+nE+nVb,nH+nE+1:nH+nE+nVb)
   Abb(nHb+nEb+1:nHb+nEb+nVb,nHb+nEb+nVb+1:nb) =              &
                                       A(nH+nE+1:nH+nE+nVb,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(nHb+nEb+1:nHb+nEb+nVb) = b(nH+nE+1:nH+nE+nVb)
!
!..L2 dof
   Abb(nHb+nEb+nVb+1:nb,1:nHb) = A(nH+nE+nV+1:nH+nE+nV+nQ,1:nHb)
   Abb(nHb+nEb+nVb+1:nb,nHb+1:nHb+nEb) = A(nH+nE+nV+1:nH+nE+nV+nQ,nH+1:nH+nEb)
   Abb(nHb+nEb+nVb+1:nb,nHb+nEb+1:nHb+nEb+nVb) = A(nH+nE+nV+1:nH+nE+nV+nQ,nH+nE+1:nH+nE+nVb)
   Abb(nHb+nEb+nVb+1:nb,nHb+nEb+nVb+1:nb) = A(nH+nE+nV+1:nH+nE+nV+nQ,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(nHb+nEb+nVb+1:nb) = b(nH+nE+nV+1:nH+nE+nV+nQ)
!
!
!..factorize Abb
!..first rewrite in packed form
   do i=1,nb
      do j=i,nb
         k = nk(i,j)
         AP(k) = Abb(i,j)
      enddo
   enddo
!
!..factorize AP
   uplo = 'U'
#if C_MODE
   call ZPPTRF(uplo, nb, AP, info)
!..store to module the LU factorization of A22
#else
   call DPPTRF(uplo, nb, AP, info)
#endif
   if (info .ne. 0) then
      write(*,*) 'stc: info = ',info
      write(*,*) 'stc: nb = ',nb
      stop 1
   endif
!
! ...compute the products of inverted AP with Abc and bb
#if C_MODE
   call ZPPTRS(uplo, nb, 1, AP, bb, nb, info1 )
   call ZGEMV('N',nb,nc,(1.0d0,0.0d0),Abc,nb,xc,1,(0.0d0,0.0d0),xb,1)
   call ZPPTRS(uplo, nb, 1, AP, xb, nb, info2 )
#else
   call DPPTRS(uplo, nb, 1, AP, bb, nb, info1 )
   call DGEMV('N',nb,nc,1.0d0,Abc,nb,xc,1,0.0d0,xb,1)
   call DPPTRS(uplo, nb, 1, AP, xb, nb, info2 )
#endif
   if (info1.ne.0) then
      write(*,*) '_herm_in: info1 = ',info1
      stop 1
   endif
!
   if (info2.ne.0) then
      write(*,*) '_herm_in: info2 = ',info2
      stop 2
   endif
   xb = bb - xb
!
!
   end subroutine stc_herm_in
!
!
! -----------------------------------------------------------------------
!
!    routine name       - stc_gen
!
! -----------------------------------------------------------------------
!
!    latest revision    - Jan 17
!
!    purpose            - routine performs static condensation for the
!                         general case
!                         (eliminates, H^1,H(curl),H(div) and L^2 bubbles)
!
!    arguments :
!      in:
!              A        - stiffness matrix
!              b        - load vector
!              n        - dimension of matrix A
!              nc       - total number of condensed dof
!              nb       - total number of bubbles
!              nHc      - number of condensed H1 dof
!              nHb      - number of bubble H1 dof
!              nEc      - number of condensed H(curl) dof
!              nEb      - number of bubble H(curl) dof
!              nVc      - number of condensed H(div) dof
!              nVb      - number of bubble H(div) dof
!              nQ       - number of bubble L^2 dof
!      out:
!              Ac       - the resulting Schur complement , condensed
!              bc       - modified right-hand side
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
   subroutine stc_gen_in(iel,Atemp,b,n,nc,nb,nHc,nHb,nEc,nEb,nVc,nVb,nQ,Abc,bb)
!
   use m_assembly, only: CLOC, STORE_STC
!
   IMPLICIT NONE
!
!----------------------------------------------------------------------
!
   integer           :: n, nc,nb,nHc,nHb,nEc,nEb,nVc,nVb,nQ, k1, k2,k
#if C_MODE
   complex*16        :: Atemp(n**2), b(n), Ac(nc**2), bc(nc), A(n,n)
   complex*16        :: Acc(nc,nc), Acb(nc,nb),Abc(nb,nc) , Abb(nb,nb)
   complex*16        :: bb(nb)
#else
   real*8            :: Atemp(n**2), b(n), Ac(nc**2), bc(nc), A(n,n)
   real*8            :: Acc(nc,nc), Acb(nc,nb),Abc(nb,nc) , Abb(nb,nb)
   real*8            :: bb(nb)
#endif
!
   integer           :: nk, nH,nE,nV,i,j,iel
   character*1       :: uplo
   integer           :: info,info1,info2
   integer, allocatable :: IPIV(:)
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
!..set up the indices.
!..check if the indices are compatible
   if(nHc+nEc+nVc .ne. nc) then
      write(*,*) 'stc: incompatible indices for the condensed system '
      write(*,*) 'iel, nHc+nEc+nVc, nc = ', iel, nHc+nEc+nVc, nc
      stop 1
   endif
   if(nHb+nEb+nVb+nQ .ne. nb) then
      write(*,*) 'stc: incompatible indices for the bubbles'
      stop 2
   endif
   if(nc+nb .ne. n) then
      write(*,*) 'stc: incompatible indices for the total system'
      stop 3
   endif

!..total dof for each space
   nH = nHc + nHb
   nE = nEc + nEb
   nV = nVc + nVb

!..extract matrices
!..change format of the matrix(for now)

   do i = 1, n
      do j = 1, n
         k = (i-1)*n + j
         A(i,j) = Atemp(k)
      enddo
   enddo
!
!
!..first Acc and bc
!..H1 dof
   Acc(1:nHc,1:nHc)                 = A(nHb+1:nH,nHb+1:nH)
   Acc(1:nHc,nHc+1:nHc+nEc)         = A(nHb+1:nH,nH+nEb+1:nH+nE)
   Acc(1:nHc,nHc+nEc+1:nHc+nEc+nVc) = A(nHb+1:nH,nH+nE+nVb+1:nH+nE+nV)
   bc(1:nHc) = b(nHb+1:nH)
!..H(curl) dof
   Acc(nHc+1:nHc+nEc,1:nHc)                 = A(nH+nEb+1:nH+nE,nHb+1:nH)
   Acc(nHc+1:nHc+nEc,nHc+1:nHc+nEc)         = A(nH+nEb+1:nH+nE,nH+nEb+1:nH+nE)
   Acc(nHc+1:nHc+nEc,nHc+nEc+1:nHc+nEc+nVc) = A(nH+nEb+1:nH+nE,nH+nE+nVb+1:nH+nE+nV)
   bc(nHc+1:nHc+nEc) = b(nH+nEb+1:nH+nE)
!..H(div) dof
   Acc(nHc+nEc+1:nHc+nEc+nVc,1:nHc)         = A(nH+nE+nVb+1:nH+nE+nV,nHb+1:nH)
   Acc(nHc+nEc+1:nHc+nEc+nVc,nHc+1:nHc+nEc) = A(nH+nE+nVb+1:nH+nE+nV,nH+nEb+1:nH+nE)
   Acc(nHc+nEc+1:nHc+nEc+nVc,nHc+nEc+1:nHc+nEc+nVc) =    &
                                              A(nH+nE+nVb+1:nH+nE+nV,nH+nE+nVb+1:nH+nE+nV)
   bc(nHc+nEc+1:nHc+nEc+nVc) = b(nH+nE+nVb+1:nH+nE+nV)
!
!..Acb
!..H1 dof
   Acb(1:nHc,1:nHb) = A(nHb+1:nH,1:nHb)
   Acb(1:nHc,nHb+1:nHb+nEb) = A(nHb+1:nH,nH+1:nH+nEb)
   Acb(1:nHc,nHb+nEb+1:nHb+nEb+nVb) = A(nHb+1:nH,nH+nE+1:nH+nE+nVb)
   Acb(1:nHc,nHb+nEb+nVb+1:nHb+nEb+nVb+nQ) = A(nHb+1:nH,nH+nE+nV+1:nH+nE+nV+nQ)
!..H(curl) dof
   Acb(nHc+1:nHc+nEc,1:nHb) = A(nH+nEb+1:nH+nE,1:nHb)
   Acb(nHc+1:nHc+nEc,nHb+1:nHb+nEb) = A(nH+nEb+1:nH+nE,nH+1:nH+nEb)
   Acb(nHc+1:nHc+nEc,nHb+nEb+1:nHb+nEb+nVb) = A(nH+nEb+1:nH+nE,nH+nE+1:nH+nE+nVb)
   Acb(nHc+1:nHc+nEc,nHb+nEb+nVb+1:nHb+nEb+nVb+nQ) = A(nH+nEb+1:nH+nE,nH+nE+nV+1:nH+nE+nV+nQ)
!..H(div) dof

   Acb(nHc+nEc+1:nHc+nEc+nVc,1:nHb) = A(nH+nE+nVb+1:nH+nE+nV,1:nHb)
   Acb(nHc+nEc+1:nHc+nEc+nVc,nHb+1:nHb+nEb) = A(nH+nE+nVb+1:nH+nE+nV,nH+1:nH+nEb)
   Acb(nHc+nEc+1:nHc+nEc+nVc,nHb+nEb+1:nHb+nEb+nVb) =    &
                                        A(nH+nE+nVb+1:nH+nE+nV,nH+nE+1:nH+nE+nVb)
   Acb(nHc+nEc+1:nHc+nEc+nVc,nHb+nEb+nVb+1:nHb+nEb+nVb+nQ) = &
                                      A(nH+nE+nVb+1:nH+nE+nV,nH+nE+nV+1:nH+nE+nV+nQ)
!..Abc
#if C_MODE
   Abc = conjg(transpose(Acb))
#else
   Abc = transpose(Acb)
#endif
!
!..Abb
!..H1 dof
   Abb(1:nHb,1:nHb) = A(1:nHb,1:nHb)
   Abb(1:nHb,nHb+1:nHb+nEb) = A(1:nHb,nH+1:nH+nEb)
   Abb(1:nHb,nHb+nEb+1:nHb+nEb+nVb) = A(1:nHb,nH+nE+1:nH+nE+nVb)
   Abb(1:nHb,nHb+nEb+nVb+1:nb) = A(1:nHb,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(1:nHb) = b(1:nHb)
!..H(curl) dof
   Abb(nHb+1:nHb+nEb,1:nHb) = A(nH+1:nH+nEb,1:nHb)
   Abb(nHb+1:nHb+nEb,nHb+1:nHb+nEb) = A(nH+1:nH+nEb,nH+1:nH+nEb)
   Abb(nHb+1:nHb+nEb,nHb+nEb+1:nHb+nEb+nVb) = A(nH+1:nH+nEb,nH+nE+1:nH+nE+nVb)
   Abb(nHb+1:nHb+nEb,nHb+nEb+nVb+1:nb) = A(nH+1:nH+nEb,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(nHb+1:nHb+nEb) = b(nH+1:nH+nEb)
!..H(div) dof
   Abb(nHb+nEb+1:nHb+nEb+nVb,1:nHb) = A(nH+nE+1:nH+nE+nVb,1:nHb)
   Abb(nHb+nEb+1:nHb+nEb+nVb,nHb+1:nHb+nEb) = A(nH+nE+1:nH+nE+nVb,nH+1:nH+nEb)
   Abb(nHb+nEb+1:nHb+nEb+nVb,nHb+nEb+1:nHb+nEb+nVb) =         &
                                       A(nH+nE+1:nH+nE+nVb,nH+nE+1:nH+nE+nVb)
   Abb(nHb+nEb+1:nHb+nEb+nVb,nHb+nEb+nVb+1:nb) =              &
                                       A(nH+nE+1:nH+nE+nVb,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(nHb+nEb+1:nHb+nEb+nVb) = b(nH+nE+1:nH+nE+nVb)
!
!..L2 dof
   Abb(nHb+nEb+nVb+1:nb,1:nHb) = A(nH+nE+nV+1:nH+nE+nV+nQ,1:nHb)
   Abb(nHb+nEb+nVb+1:nb,nHb+1:nHb+nEb) = A(nH+nE+nV+1:nH+nE+nV+nQ,nH+1:nH+nEb)
   Abb(nHb+nEb+nVb+1:nb,nHb+nEb+1:nHb+nEb+nVb) = A(nH+nE+nV+1:nH+nE+nV+nQ,nH+nE+1:nH+nE+nVb)
   Abb(nHb+nEb+nVb+1:nb,nHb+nEb+nVb+1:nb) = A(nH+nE+nV+1:nH+nE+nV+nQ,nH+nE+nV+1:nH+nE+nV+nQ)
   bb(nHb+nEb+nVb+1:nb) = b(nH+nE+nV+1:nH+nE+nV+nQ)
!
!
!..factorize Abb
   allocate(IPIV(nb))
!
#if C_MODE
      call ZGETRF(nb,nb,Abb,nb,IPIV,info)
!..store to module the LU factorization of A22
#else
      call DGETRF(nb,nb,Abb,nb,IPIV,info)
#endif
      if (info .ne. 0) then
         write(*,*) 'stc_gen: info = ',info
         stop 1
      endif
! ...compute the products of inverted AP with Abc and bb
#if C_MODE
      call ZGETRS('N',nb,1,Abb,nb,IPIV,bb,nb,info1)
#else
      call DGETRS('N',nb,1,Abb,nb,IPIV,bb,nb,info1)
#endif
      if (info1.ne.0) then
         write(*,*) 'stc: info1 = ',info1
         stop
      endif
!
#if C_MODE
      call ZGETRS('N',nb,nc,Abb,nb,IPIV,Abc,nb,info2)
#else
      call DGETRS('N',nb,nc,Abb,nb,IPIV,Abc,nb,info2)
#endif
      if (info2.ne.0) then
         write(*,*) 'stc: info2 = ',info2
         stop
      endif
      deallocate(IPIV)

   end subroutine stc_gen_in
!
!

