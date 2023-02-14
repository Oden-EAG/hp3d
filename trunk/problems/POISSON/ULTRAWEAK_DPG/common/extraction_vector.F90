!-----------------------------------------------------------------------
! routine name - extraction_vector
!---------------------------------------------------------------------

!last revision: 25 Nov 2022

!Purpose: extracts vector for location of smaller symmetric submatrix in larger symmetric matrix
!         Here larger matrix corresponds to a L2 gram matrix of larger space and smaller submatrix
!         is the L2 gram matrix of the smaller space
!arguments
! in:
!     Nord_current: order for the current step 
!     Nord_glob : order which resulted in the largest system.
!     nrdofmQ : L2 dofs for the current step
!     nrdofgQ : L2 dofs for the global matrix

    
! out: 

!     Nextract: array containing the indices of rows or columns corresponding to Nord_current
!------------------------------------------------------------------------------------





subroutine extraction_vector(Nord_current,Nord_glob,nrdofmQ,nrdofgQ,Nextract)


    implicit none

    integer, intent(in) ::  Nord_current
    integer, intent(in) ::  Nord_glob
    integer,    intent(in) ::  nrdofmQ
    integer,    intent(in) ::  nrdofgQ

    integer,dimension(nrdofmQ),   intent(out) :: Nextract


    integer :: nordx,nordy,nordz
    integer :: nordxg,nordyg,nordzg

    integer :: i,j,k,m
    
    call ddecode(Nord_current,nordx,nordy,nordz)

    call ddecode(Nord_glob,nordxg,nordyg,nordzg)


    m = 0

    ! L2 order is nordx-1, nordy-1 and nordz - 1 due to exact sequence.
    do  k = 0,nordz-1
        do j = 0,nordy-1
            do i = 0,nordx-1
            m = m + 1
            Nextract(m) = k * (nordxg * nordyg) + j * nordxg + i + 1
            enddo
        enddo
    enddo





end subroutine extraction_vector


subroutine extraction_vector_new(Nord_prev,Nord_current,Nord_glob,nrdofmQ,nrdofgQ,Nextract_prev,Nextract)


    implicit none

    integer, intent(in) ::  Nord_prev
    integer, intent(in) ::  Nord_current
    integer, intent(in) ::  Nord_glob
    integer,    intent(in) ::  nrdofmQ
    integer,    intent(in) ::  nrdofgQ

    integer,dimension(nrdofmQ),   intent(out) :: Nextract
    integer,    intent(in)  :: Nextract_prev(*)

    integer :: nordxp,nordyp,nordzp
    integer :: nordx,nordy,nordz
    integer :: nordxg,nordyg,nordzg

    integer :: i,j,k,m
    
    call ddecode(Nord_current,nordx,nordy,nordz)
    call ddecode(Nord_prev,nordxp,nordyp,nordzp)
    call ddecode(Nord_glob,nordxg,nordyg,nordzg)


    m = 0

    ! L2 order is nordx-1, nordy-1 and nordz - 1 due to exact sequence.

    !adding old indexes
    ! do  k = 0,nordzp-1
    !     do j = 0,nordyp-1
    !         do i = 0,nordxp-1
    !         m = m + 1
    !         Nextract(m) = k * (nordxg * nordyg) + j * nordxg + i + 1
    !         enddo
    !     enddo
    ! enddo

    do i = 1,nordxp*nordyp*nordzp
        m = m + 1
        Nextract(m) = Nextract_prev(m)
    enddo



    !adding new indexes
    do k = 0,nordzp-1
        do j = 0,nordyp-1
            do i = nordxp,nordx-1
                m = m + 1
                Nextract(m) =  k * (nordxg * nordyg) + j * nordxg + i + 1
            enddo
        enddo
    enddo
    
    do k = 0,nordzp-1
        do j = nordyp,nordy-1
            do i = 0,nordx-1
                m = m + 1
                Nextract(m) =  k * (nordxg * nordyg) + j * nordxg + i + 1

            enddo
        enddo
    enddo





    do k = nordzp,nordz-1
        do j = 0,nordy-1
            do i = 0,nordx-1
                m = m + 1
                Nextract(m) =  k * (nordxg * nordyg) + j * nordxg + i + 1

            enddo
        enddo
    enddo



end subroutine extraction_vector_new