




subroutine project_h(Mdle,error_org,error_hopt,kref_out)    


use control
use data_structure3D
use element_data
use parametersDPG

!
implicit none
integer,    intent(in)  :: Mdle
real(8),    intent(in)  :: error_org
real(8),    intent(out) :: error_hopt
integer,    intent(out) :: kref_out


integer :: ref_opts
integer ::  k,j,m,l  
character(len=4) :: etype

real(8), allocatable   ::  error_opt(:)
real(8), allocatable   ::  g_rate_ref(:)
integer, allocatable   ::  kref_opts(:)

real(8) :: max_rate_element

etype = NODES(Mdle)%type

if(etype .eq. 'mdlb') then
    ref_opts = 8
    
    allocate(error_opt(ref_opts))
    allocate(kref_opts(ref_opts))
    allocate(g_rate_ref(ref_opts))
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
        
        call elem_proj_h(Mdle,kref_opts(j),error_org,error_opt(j),g_rate_ref(j))
        write(*,*) kref_opts(j),error_opt(j),g_rate_ref(j)
    enddo


    max_rate_element = maxval(g_rate_ref)
    write(*,*) "Max rate is = ", max_rate_element

endif





end subroutine project_h