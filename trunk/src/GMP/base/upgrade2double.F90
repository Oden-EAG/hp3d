!========================================================================================
!  subroutine - upgrade2double                                                          |
!========================================================================================
!  LATEST REVISION: Jul 09                                                              |
!
!  PURPOSE: routine upgrades to real(8) geometry using Newton
!           method iterations
!
!  REMARKS: new general routine!
!========================================================================================
subroutine upgrade2double
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!----------------------------------------------------------------------------------------
! MODULES
  use U2D
  use control
!----------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------
! VARIABLES
  integer                   :: np,nr_surf
  integer, dimension(MAXSU) :: SURFS
#if I_PRINT >= 2
  real(8), dimension(3)     :: aux
  real(8)                   :: norm_new,norm_old
  integer                   :: np_MAX
!----------------------------------------------------------------------------------------
! FUNCTIONS
  real(8), external          :: norm
#endif
!----------------------------------------------------------------------------------------
!
#if I_PRINT >= 1
    write(*,*)'upgrade2double: upgrading to real(8) geometry...'
#endif
! ..allocate
    call generate_POINT_TYPE
    call allocate_EXTRA_PLANES
#if I_PRINT >= 2
    write(*,*)'upgrade2double: POINT_TYPE generated.'
    write(*,*)'upgrade2double: EXTRA_PLANES allocated.'
    norm_old = 0.d0
    np_MAX = 0
#endif
! ..loop through points
    do np = 1, NRPOINT
      nr_surf = POINT_TYPE(1,np)
! ....if point is not on any surface cycle
      if (nr_surf .eq. 0 )  cycle
      if (nr_surf .gt. MAXSU) then
        write(*,*)'upgrade2double: increase MAX_surf.'
        stop
      endif
      SURFS(1:nr_surf) = POINT_TYPE(2:nr_surf+1,np)
#if I_PRINT >= 2
      aux(1:3) = POINTS(np)%Rdata(1:3)
#endif
      call reproject_point(np,nr_surf,SURFS(1:nr_surf))
#if I_PRINT >= 2
      aux(1:3) = aux(1:3) - POINTS(np)%Rdata(1:3)
      norm_new = norm(aux)
      if (norm_new .gt. norm_old) then
        norm_old = norm_new
        np_MAX = np
      endif
      if (norm_new .gt. GEOM_TOL) then
        write(*,*)'upgrade2double: exceeding GEOM_TOL, np    = ',np
        write(*,*)'upgrade2double:                     shift = ',norm_old
      endif
#endif
    enddo
! ..deallocate
    call deallocate_EXTRA_PLANES
    call deallocate_POINT_TYPE
#if I_PRINT >= 2
    write(*,*)'upgrade2double: np_MAX    = ',np_MAX
    write(*,*)'upgrade2double: max shift = ',norm_old
    call print_GMP
#endif
#if I_PRINT >= 1
write(*,*)'upgrade2double: done!'
#endif
!
end subroutine upgrade2double
!
!
!-----------------------------------------------------------------------------------
subroutine reproject_point(Np,Nr_surf,Surfs)
!-----------------------------------------------------------------------------------
  use GMP
!-----------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------
! DUMMY ARGUMENTS
! point number
  integer, intent(in)                   :: Np
! number of surfaces it will be reprojected on
  integer, intent(in)                   :: Nr_surf
! list of surfaces
  integer, dimension(MAXSU), intent(in) :: Surfs
!-----------------------------------------------------------------------------------
! VARIABLES
  integer, dimension(3) :: N_SURF
  real(8), dimension(3)  :: VOID_1
  real(8), dimension(4)  :: VOID_2
!-----------------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'reproject_point: reprojecting point Np = ',Np
#endif
    if (Nr_surf .lt. 1) then
      write(*,*)'reproject_point: WARNING, Nr_surf = 0'
    else
! ....select appropriate surfaces
      call select_surfaces(Np,Nr_surf,Surfs, N_SURF)
! ....apply Newton method
      VOID_1 = 0.d0;  VOID_2 = 0.d0
      call mnewt(1,N_SURF(1:3),VOID_1,VOID_2,POINTS(Np)%Rdata(1:3),VOID_2, POINTS(Np)%Rdata(1:3))
    endif
#if I_PRINT >= 1
    write(*,*)'reproject_point: done!'
#endif
!
end subroutine reproject_point
!
!
!---------------------------------------------------------------------------------------
subroutine select_surfaces(Np,Nr_surf,Surfs, SEL_Surfs)
!---------------------------------------------------------------------------------------
! MODULES
  use GMP
  use control
!---------------------------------------------------------------------------------------
IMPLICIT NONE
!---------------------------------------------------------------------------------------
! DUMMY ARGUMENTS
! point number
  integer, intent(in)                   :: Np
! number of surfaces the point lies on
  integer, intent(in)                   :: Nr_surf
! list of surfaces
  integer, dimension(MAXSU), intent(in) :: Surfs
! best choice of 3 surfaces
  integer, dimension(3), intent(out)    :: SEL_Surfs
!----------------------------------------------------------------------------------------
! VARIABLES
  real(8)                 :: prod_new,prod_old
  integer                 :: i,j,k,l,i1,i2
  real(8), dimension(3,3) :: GRAD_old, GRAD_new
  real(8), dimension(3)   :: VOID_1
!----------------------------------------------------------------------------------------
! FUNCTIONS
  integer, external :: combination,my_mod
!----------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'select_surfaces: selecting 3 surfaces for Np = ',Np
#endif
#if I_PRINT >= 2
    write(*,*)'select_surfaces: Nr_surf  = ',Nr_surf
    write(*,*)'select_surfaces: surfaces = ',Surfs(1:Nr_surf)
#endif
! ..select number of surfaces
    select case(Nr_surf)
! ....if no surface, do nothing
      case(0)
        write(*,*)'select_surfaces: WARNING, Nr_surf = 0'
! ....1 surface
      case(1)
        call select_1surf
! ....2 surfaces
      case(2)
        call select_2surf
! ....3 or more surfaces
      case default
        call select_3surf
    endselect
#if I_PRINT >= 1
    write(*,*)'select_surfaces: done!'
#endif
    contains
!
!
!
!------------------------------------------------------------------------------------------------
subroutine select_1surf
!------------------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT_ONE 0
!
#if I_PRINT_ONE >= 1
    write(*,*)'select_1surf: 1 surf point found.'
#endif
#if I_PRINT_ONE >= 2
    write(*,*)'select_1surf: surface = ',Surfs(1:Nr_surf)
#endif
    SEL_Surfs(1) = Surfs(1)
    SEL_Surfs(2) = NRSURFS - 1
    SEL_Surfs(3) = NRSURFS
#if I_PRINT_ONE >= 2
    write(*,*)'select_1surf: selected surfaces = ',SEL_Surfs
#endif
! ..compute normal
    call surf(SEL_Surfs(1),POINTS(Np)%Rdata(1:3), VOID_1,GRAD_new(1:3,1))
    call normalize(GRAD_new(1:3,1))
! ..rotate normal by 90 deg
    do i = 1, 3
      i1 = my_mod(i + 1,3)
      i2 = my_mod(i + 2,3)
! ....choose direction of rotation based on NORMAL(i)
      if (GRAD_new(i,1) .ne. 0d0) then
        GRAD_new(i1,2) =  GRAD_new(i,1)
        GRAD_new(i,2)  = -GRAD_new(i1,1)
        GRAD_new(i2,2) = 0.d0
        exit
      endif
    enddo
    call normalize(GRAD_new(1:3,2))
! ..determine 3rd vector
    call cross_product(GRAD_new(1:3,1),GRAD_new(1:3,2), GRAD_new(1:3,3))
! ..create 2 extra planes
    SURFACES(SEL_Surfs(2))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
    SURFACES(SEL_Surfs(2))%Rdata(4:6) = GRAD_new(1:3,2)
!
    SURFACES(SEL_Surfs(3))%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
    SURFACES(SEL_Surfs(3))%Rdata(4:6) = GRAD_new(1:3,3)
#if I_PRINT_ONE >= 1
    write(*,*)'select_1surf: done!'
#endif
!
end subroutine select_1surf
!------------------------------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------------------------------
subroutine select_2surf
!------------------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT_TWO 0
!
#if I_PRINT_TWO >= 1
    write(*,*)'select_2surf: 2 surf point found.'
#endif
#if I_PRINT_TWO >= 2
    write(*,*)'select_2surf: surfaces = ',Surfs(1:Nr_surf)
#endif
    prod_old = 1.d10
! ..loop through all possible choices of 2 surfaces
    l = 0
    do i = 1, (Nr_surf - 1)
      do j = 1, Nr_surf
! ......increment counter
        l = l + 1
! ......compute gradients of 2 selected surfaces and normalize them
        call surf(Surfs(i),POINTS(Np)%Rdata(1:3), VOID_1,GRAD_new(1:3,1))
        call surf(Surfs(j),POINTS(Np)%Rdata(1:3), VOID_1,GRAD_new(1:3,2))
        call normalize(GRAD_new(1:3,1))
        call normalize(GRAD_new(1:3,2))
! ......compute absolute value of scalar product and store pertainig surfaces
        call scalar_product(GRAD_new(1:3,1),GRAD_new(1:3,2), prod_new)
        prod_new = abs(prod_new)
! ......if scalar product is smaller, update
        if (prod_new .lt. prod_old) then
          SEL_Surfs(1)    = Surfs(i)
          SEL_Surfs(2)    = Surfs(j)
          prod_old        = prod_new
          GRAD_old(1:3,1) = GRAD_new(1:3,1)
          GRAD_old(1:3,2) = GRAD_new(1:3,2)
        endif
      enddo
    enddo
! ..if the two normals are almost parallel
    if (abs(prod_old - 1.d0) .lt. GEOM_TOL) then
! ....treat as 1 surf point
#if I_PRINT_TWO >= 2
      write(*,*)'select_2surf: point redefined as 1 surf.'
#endif
      call select_1surf
    else
! ....set up 3rd surface
      SEL_Surfs(3) = NRSURFS
      call cross_product(GRAD_old(1:3,1),GRAD_old(1:3,2), GRAD_old(1:3,3))
      call normalize(GRAD_old(1:3,3))
      SURFACES(NRSURFS)%Rdata(1:3) = POINTS(Np)%Rdata(1:3)
      SURFACES(NRSURFS)%Rdata(4:6) = GRAD_old(1:3,3)
    endif
#if I_PRINT_TWO >= 2
    write(*,*)'select_2surf: selected surfaces = ',SEL_Surfs
#endif
#if I_PRINT_TWO >= 1
    write(*,*)'select_2surf: done!'
#endif
!
end subroutine select_2surf
!-----------------------------------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------------------------------
subroutine select_3surf
!-----------------------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT_THREE 0
!
#if I_PRINT_THREE >= 1
    write(*,*)'select_3surf: 3 surf point found.'
#endif
#if I_PRINT_THREE >= 2
    write(*,*)'select_3surf: surfaces = ',Surfs(1:Nr_surf)
#endif
! ..loop through all possible choices of 3 surfaces
    prod_old = 0.d0
    l = 0
    do i = 1, (Nr_surf - 2)
      do j = (i + 1), (Nr_surf - 1)
        do k = (j + 1), Nr_surf
! ........increment counter
          l = l + 1
! ........compute gradients of 3 selected surfaces and normalize them
          call surf(Surfs(i),POINTS(Np)%Rdata(1:3), VOID_1,GRAD_new(1:3,1))
          call surf(Surfs(j),POINTS(Np)%Rdata(1:3), VOID_1,GRAD_new(1:3,2))
          call surf(Surfs(k),POINTS(Np)%Rdata(1:3), VOID_1,GRAD_new(1:3,3))
          call normalize(GRAD_new(1:3,1))
          call normalize(GRAD_new(1:3,2))
          call normalize(GRAD_new(1:3,3))
! ........compute absolute value of mixed product and store pertaining surfaces
          call mixed_product(GRAD_new(1:3,1),GRAD_new(1:3,2),GRAD_new(1:3,3), prod_new)
          prod_new = abs(prod_new)
          if (prod_new .gt. prod_old) then
            SEL_Surfs(1) = Surfs(i)
            SEL_Surfs(2) = Surfs(j)
            SEL_Surfs(3) = Surfs(k)
            prod_old     = prod_new
          endif
        enddo
      enddo
    enddo
#if I_PRINT_THREE >= 2
! ..run additional check
    if (combination(Nr_surf,3) .ne. l) then
      write(*,*)'select_3surf: something is wrong...'
      write(*,*) 'l = ',l
      write(*,*) 'combinations = ',combination(Nr_surf,3)
      stop
    endif
#endif
! ..if less than GEOM_TOL redefines as a 2 surfs point
    if (prod_old .lt. GEOM_TOL) then
#if I_PRINT_THREE >= 2
      write(*,*)'select_3surf: redefining point as 2 surf.'
#endif
      call select_2surf
    else
#if I_PRINT_THREE >= 2
      write(*,*)'select_3surf: found 3 appropriate surfaces = ',SEL_Surfs
#endif
    endif
!
end subroutine select_3surf
!-----------------------------------------------------------------------------------------------
!
!
end subroutine select_surfaces
!-----------------------------------------------------------------------------------------------
