
! input 
! Mdle : element number
! flag_pref : flag stating that the element is p-refined or not during fine mesh formation
! error_org : projection error at the original polynomial and h configuration.
! rate_p : highest error rate obtained from project_p
! Polyflag: polynomial refinement flag associated with rate_p

! output
! elem_grate: guranteed rate for the element
! Ref_indicator_flag: array containing information for type of refinement and polynomial order of the sons generated if h-ref is preferred.


subroutine project_h(Mdle,flag_pref_loc,error_org,rate_p,Poly_flag,elem_grate,Ref_indicator_flag)    


use control
use data_structure3D
use element_data
use parametersDPG

!
implicit none
integer,    intent(in)  :: Mdle !middle node number 
integer,    intent(in)  :: flag_pref_loc    !flag for refinement during fine hp mesh generation
real(8),    intent(in)  :: error_org    ! projection error for the father element
real(8),    intent(in)  :: rate_p   ! error rate for optimal p-refinement
integer,    intent(in)  :: Poly_flag    ! polynomial flag or order corresponding to opt p-refinement
real(8),    intent(out) :: elem_grate   ! element guranteed rate
integer,    intent(out) :: Ref_indicator_flag(10)   ! refinement indicator flag array


integer :: ref_opts !number of h-refinement options for an elemenet type
integer ::  k,j,m,l  
character(len=4) :: etype

real(8), allocatable   ::  error_opt(:)     !projection error for each refinement options
real(8), allocatable   ::  g_rate_ref(:)    !guranteed rate for each refinement option
integer, allocatable   ::  kref_opts(:)     !refinement markers ffor each refinement option
integer, allocatable   ::  Nord_max_href(:,:)   !array containing the poly distributio for each ref opt for which we obtain guranteed rate

real(8), allocatable   :: rate_hcomp(:) !array to store error rate for competitive h-refinements
real(8) :: max_rate_element  !maximum rate over all refinement options
real(8) :: max_rate_hcomp, error_max_rate  ! maximum rate for comp refs and error for correspoding comp-href
! real(8) :: rate_p
integer :: max_rate_element_loc(1)

integer :: nrdofgQ, nrdof_org, nrdof_polyref
integer :: Nord_glob, Nord_org ,px,py,pz
integer :: norder(19)

etype = NODES(Mdle)%type


if(etype .eq. 'mdlb') then

    call find_order(Mdle, norder)
    Nord_glob = norder(19) 
    call ddecode(Nord_glob,px,py,pz)
    nrdofgQ = px * py *pz
    
    if(flag_pref_loc .eq. 1) then
        Nord_org = norder(19) - 111
    else
        Nord_org = norder(19)
    endif

    call ddecode(Nord_org,px,py,pz)
    nrdof_org = px * py * pz

    ref_opts = 8
    
    allocate(error_opt(ref_opts))
    error_opt = ZERO

    allocate(kref_opts(ref_opts))
    kref_opts = ZERO

    allocate(g_rate_ref(ref_opts))
    g_rate_ref = ZERO

    allocate(rate_hcomp(ref_opts))
    rate_hcomp = ZERO

    allocate(Nord_max_href(ref_opts,8))
    Nord_max_href = ZERO

    m = 1
    do k = 0,1
    
        do j = 0,1
            
            do l = 0,1
                ! kref_opts(m,3) = l
                ! kref_opts(m,2) = j
                ! kref_opts(m,1) = k
                kref_opts(m) = k * 10**(2) + j * 10 **1 + l * 10**0
                m = m + 1
            enddo

        enddo
    
    enddo

    m = m - 1

    ! j should start from 2 as firts eleent of Kref_opt is 000 which means no refinement which is useless.
    do j = 2,m
        
        call elem_proj_h(Mdle,flag_pref_loc,kref_opts(j),error_org,error_opt(j),g_rate_ref(j),rate_hcomp(j),Nord_max_href(j,1:8))
        ! write(*,*) Nord_max_href(j,1:8)
        ! write(*,*) kref_opts(j),error_opt(j),g_rate_ref(j)
    enddo


    max_rate_element = maxval(g_rate_ref(2:m))
    error_max_rate  = maxval(error_opt(2:m))
    ! max_rate_element_href(1) = maxloc(g_rate_ref)
    max_rate_element_loc =  maxloc(g_rate_ref(2:m))
    max_rate_hcomp = maxval(rate_hcomp(2:m))
    ! write(*,*) "Max elemental rate is = ", max_rate_element
    ! write(*,*) "Max elemental comp rate is = ", max_rate_hcomp

    !selection between href or pref
    ! rate_p = (log(error_org**2) - log(error_p**2))/(log(real(nrdofgQ,8)) - log(real(nrdof_org,8)))
    ! call ddecode(Poly_flag,px,py,pz)
    ! nrdof_polyref = px * py * pz

    if(flag_pref_loc .eq. 0) then
        ! rate_p = 0.0

        Ref_indicator_flag(1) = 1
        !since 1 element of kref_opts means 000, hence we need to offset it by one. Thus,we have max_rate_element(1) + 1
        Ref_indicator_flag(2) = kref_opts(max_rate_element_loc(1)+1)
        Ref_indicator_flag(3:10) = Nord_max_href(max_rate_element_loc(1) + 1,1:8)
        elem_grate = max_rate_element

    else
        ! rate_p = (error_org - error_p)/(real(nrdof_polyref,8) - real(nrdof_org,8))

        if(max_rate_hcomp .ge. rate_p) then !href

            Ref_indicator_flag(1) = 1
            !since 1 element of kref_opts means 000, hence we need to offset it by one. Thus,we have max_rate_element(1) + 1
            Ref_indicator_flag(2) = kref_opts(max_rate_element_loc(1)+1)
            Ref_indicator_flag(3:10) = Nord_max_href(max_rate_element_loc(1) + 1,1:8)
            elem_grate = max_rate_element
    
        else
    
            Ref_indicator_flag(1) = 2
            Ref_indicator_flag(2) = 000
            Ref_indicator_flag(3) = Poly_flag !p-refinement by degree 1.
            elem_grate = rate_p
    
        endif

    endif

    write(*,*) "Rate for elemental p-ref = ",rate_p,max_rate_hcomp,Mdle,sqrt(error_max_rate),error_org

    ! if(max_rate_hcomp .ge. rate_p) then !href

    !     Ref_indicator_flag(1) = 1
    !     !since 1 element of kref_opts means 000, hence we need to offset it by one. Thus,we have max_rate_element(1) + 1
    !     Ref_indicator_flag(2) = kref_opts(max_rate_element_loc(1)+1)
    !     Ref_indicator_flag(3:10) = Nord_max_href(max_rate_element_loc(1) + 1,1:8)
    !     elem_grate = max_rate_element

    ! else

    !     Ref_indicator_flag(1) = 2
    !     Ref_indicator_flag(2) = 000
    !     Ref_indicator_flag(3) = Poly_flag !p-refinement by degree 1.
    !     elem_grate = rate_p

    ! endif

    ! write(*,*) "The elemental Refine_indicator_flag is  = ", Ref_indicator_flag

endif





end subroutine project_h