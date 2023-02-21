!---------------------------------------------------------------------
! routine name : poly_adap_subson
!---------------------------------------------------------------------


!last revision: 17 Jan 2023

!purpose: computes optimal p-refinement for subson

! arguments
! in:
! Mdle: index for the father

!out:

!loc_poly_flag: the flag for the optimal

subroutine poly_adap_subson(kref_loc,Mdle,nrdofgQ,nrdof_org,Nord_old,Nord_glob,subson_overlap,base_error, &
                            Ap,zbload,Awork_prev,Bwork_prev,Nextract_prev,Polyflag)


    use control
    use data_structure3D
    use element_data
    use parametersDPG

    !
    implicit none
    integer,    intent(in)  :: kref_loc
    integer,    intent(in)  :: Mdle
    integer,    intent(in)  :: nrdofgQ
    integer,    intent(in)  :: nrdof_org
    integer,    intent(in)  :: Nord_old
    integer,    intent(in)  :: Nord_glob

    real(8),    intent(in)  ::  base_error !error at order at which we start subson telescopic solves
    real(8), dimension(nrdofgQ,nrdofgQ),    intent(in) :: Ap
    real(8), dimension(nrdofgQ),    intent(in) ::zbload
    real(8), dimension(nrdofgQ,nrdofgQ),    intent(in)   :: Awork_prev
    real(8), dimension(nrdofgQ),    intent(in)   :: Bwork_prev
    integer, intent(in) :: Nextract_prev(*)
    integer, intent(in) :: subson_overlap(*)
    integer, intent(out) :: Polyflag
    ! integer, dimension()

    real(8),    allocatable     ::   Ap_subson(:,:)
    real(8),    allocatable     ::   zbload_subson(:)
    real(8),    allocatable     ::   Awork_subson(:,:)
    real(8),    allocatable     ::   Bwork_subson(:)
    integer,    allocatable     ::   Nextract_prev_subson(:)
    integer,    allocatable     ::   Nextract_subson(:)

    integer :: nrdofmQ,pxm,pym,pzm
    real(8) :: coeff_sum(3)
    real(8) :: coeff_max
    real(8) :: error_rate_aniso(3)
    integer :: poly_order_added(3)
    integer :: poly_flag_store(3)
    integer :: max_loc(1), Nrhs, Ldwork, Ldglob
    integer,dimension(50)  :: Mblock
    integer :: mstep
    character(len=4) :: etype
!---------------------------------------------------------
!       aux variables
    integer ::  k,k1,k2,k3,j
    integer ::  Nord_mod, order_add,Nord_prev, nrdof_old
    real(8) ::  proj_error,rate_p
 

    allocate(Ap_subson(nrdofgQ,nrdofgQ))
    allocate(zbload_subson(nrdofgQ))
    allocate(Awork_subson(nrdofgQ,nrdofgQ))
    allocate(Bwork_subson(nrdofgQ))

    Ap_subson(1:nrdofgQ,1:nrdofgQ) = Ap(1:nrdofgQ,1:nrdofgQ)
    zbload_subson(1:nrdofgQ) = zbload(1:nrdofgQ)
    Awork_subson(1:nrdofgQ,1:nrdofgQ) = Awork_prev(1:nrdofgQ,1:nrdofgQ)
    Bwork_subson(1:nrdofgQ) = Bwork_prev(1:nrdofgQ)

    coeff_sum  = ZERO
    Ldglob = nrdofgQ    


    Nrhs = 1
    Mblock = ZERO
    etype = NODES(Mdle)%type

    if(etype .eq. 'mdlb') then

        call ddecode(Nord_old,pxm,pym,pzm)
        nrdofmQ = pxm * pym * pzm
        nrdof_old = nrdofmQ

        allocate(Nextract_prev_subson(nrdofgQ))
        Nextract_prev_subson(1:nrdofmQ) = Nextract_prev(1:nrdofmQ)


        call coeff_polynomial_aniso(Nord_old,Nextract_prev_subson(1:nrdofmQ),etype,Bwork_subson(1:nrdofmQ),coeff_sum)


        mstep = 1
        Mblock(mstep) = 0
        Mblock(mstep + 1) = nrdofmQ


        error_rate_aniso = ZERO
        poly_order_added = ZERO
        poly_flag_store = ZERO

        coeff_max = 0.d0
        k2 = 0
        Nord_mod = Nord_old
        mstep = 2
        

        do k = 1,3

            do k1 = 1,3
                if((coeff_max .le. coeff_sum(k1)) .and. (poly_order_added(k1) .eq. 0)) then
                    coeff_max = coeff_sum(k1)
                    k2 = k1
                endif
            enddo

            poly_order_added(k2) = 1
            order_add = int(10**(3-k2))
            Nord_prev = Nord_mod
            Nord_mod =  Nord_mod + order_add
            call ddecode(Nord_mod,pxm,pym,pzm)

            nrdofmQ = pxm * pym * pzm
            Mblock(mstep + 1) = nrdofmQ
            Ldwork = nrdofmQ

            allocate(Nextract_subson(nrdofmQ))
            Nextract_subson = ZERO
            call extraction_vector_new(Nord_prev,Nord_mod,Nord_glob,nrdofmQ,nrdofgQ,Nextract_prev_subson,Nextract_subson)
            call pbisolver3(mstep,Mblock,Ap_subson,Ldglob,Awork_subson,Ldwork,Nextract_subson,zbload_subson,Bwork_subson,Nrhs)
            call fine_to_subson_projection_error(kref_loc,Bwork_subson(1:nrdofmQ),Nextract_subson,subson_overlap, &
                                                    nrdofmQ,nrdofgQ,Mdle,proj_error)

            Nextract_prev_subson = ZERO
            Nextract_prev_subson(1:nrdofmQ) = Nextract_subson(1:nrdofmQ)
            deallocate(Nextract_subson)

            error_rate_aniso(k) = (base_error - proj_error)/(abs(real(nrdofmQ,8) - real(nrdof_org,8)))
            poly_flag_store(k) = Nord_mod

            mstep = mstep + 1
            coeff_sum = ZERO

            call coeff_polynomial_aniso(Nord_mod,Nextract_prev_subson(1:nrdofmQ),etype,Bwork_subson(1:nrdofmQ),coeff_sum)
            coeff_max = ZERO

        enddo

        rate_p = maxval(error_rate_aniso)
        max_loc = maxloc(error_rate_aniso)
        Polyflag = poly_flag_store(max_loc(1))




    endif
    









end subroutine poly_adap_subson




!-----------------------new poly_adap_subson -------------------------------------------------


subroutine poly_adap_subson_new(kref_loc,Mdle,nrdofgQ,nrdof_org,Nord_old,Nord_glob,subson_overlap,base_error, &
                                Ap,zbload,Awork_prev,Bwork_prev,Nextract_prev,Polyflag,&
                                weights_fine_store,quad_point_store,nint_pp_store,&
                                shap3DQ_fine_store,shap3DQ_coarse_store)


use control
use data_structure3D
use element_data
use parametersDPG

!
implicit none
integer,    intent(in)  :: kref_loc
integer,    intent(in)  :: Mdle
integer,    intent(in)  :: nrdofgQ
integer,    intent(in)  :: nrdof_org
integer,    intent(in)  :: Nord_old
integer,    intent(in)  :: Nord_glob

real(8),    intent(in)  ::  base_error !error at order at which we start subson telescopic solves
real(8), dimension(nrdofgQ,nrdofgQ),    intent(in) :: Ap
real(8), dimension(nrdofgQ),    intent(in) ::zbload
real(8), dimension(nrdofgQ,nrdofgQ),    intent(in)   :: Awork_prev
real(8), dimension(nrdofgQ),    intent(in)   :: Bwork_prev
integer, intent(in) :: Nextract_prev(*)
integer, intent(in) :: subson_overlap(*)
integer, intent(out) :: Polyflag


real(8),dimension(MAXNINT3ADD,2,*),intent(in) :: weights_fine_store
real(8),dimension(3,MAXNINT3ADD,*),intent(in) :: quad_point_store
real(8),dimension(MAXbrickQQ,MAXNINT3ADD,*), intent(in) :: shap3DQ_fine_store
real(8),dimension(MAXbrickQQ,MAXNINT3ADD,*), intent(in) :: shap3DQ_coarse_store
integer,dimension(*),intent(in) :: nint_pp_store
! integer, dimension()

real(8),    allocatable     ::   Ap_subson(:,:)
real(8),    allocatable     ::   zbload_subson(:)
real(8),    allocatable     ::   Awork_subson(:,:)
real(8),    allocatable     ::   Bwork_subson(:)
integer,    allocatable     ::   Nextract_prev_subson(:)
integer,    allocatable     ::   Nextract_subson(:)

integer :: nrdofmQ,pxm,pym,pzm
real(8) :: coeff_sum(3)
real(8) :: coeff_max
real(8) :: error_rate_aniso(3)
integer :: poly_order_added(3)
integer :: poly_flag_store(3)
integer :: max_loc(1), Nrhs, Ldwork, Ldglob
integer,dimension(50)  :: Mblock
integer :: mstep
character(len=4) :: etype
!---------------------------------------------------------
!       aux variables
integer ::  k,k1,k2,k3,j,l
integer ::  Nord_mod, order_add,Nord_prev, nrdof_old
real(8) ::  proj_error,rate_p


allocate(Ap_subson(nrdofgQ,nrdofgQ))
allocate(zbload_subson(nrdofgQ))
allocate(Awork_subson(nrdofgQ,nrdofgQ))
allocate(Bwork_subson(nrdofgQ))

Ap_subson(1:nrdofgQ,1:nrdofgQ) = Ap(1:nrdofgQ,1:nrdofgQ)
zbload_subson(1:nrdofgQ) = zbload(1:nrdofgQ)
Awork_subson(1:nrdofgQ,1:nrdofgQ) = Awork_prev(1:nrdofgQ,1:nrdofgQ)
Bwork_subson(1:nrdofgQ) = Bwork_prev(1:nrdofgQ)

coeff_sum  = ZERO
Ldglob = nrdofgQ    


Nrhs = 1
Mblock = ZERO
etype = NODES(Mdle)%type

if(etype .eq. 'mdlb') then

    call ddecode(Nord_old,pxm,pym,pzm)
    nrdofmQ = pxm * pym * pzm
    nrdof_old = nrdofmQ

    allocate(Nextract_prev_subson(nrdofgQ))
    Nextract_prev_subson(1:nrdofmQ) = Nextract_prev(1:nrdofmQ)


    call coeff_polynomial_aniso(Nord_old,Nextract_prev_subson(1:nrdofmQ),etype,Bwork_subson(1:nrdofmQ),coeff_sum)


    mstep = 1
    Mblock(mstep) = 0
    Mblock(mstep + 1) = nrdofmQ


    error_rate_aniso = ZERO
    poly_order_added = ZERO
    poly_flag_store = ZERO

    coeff_max = 0.d0
    k2 = 0
    Nord_mod = Nord_old
    mstep = 2
    

    do k = 1,3

        do k1 = 1,3
            if((coeff_max .le. coeff_sum(k1)) .and. (poly_order_added(k1) .eq. 0)) then
                coeff_max = coeff_sum(k1)
                k2 = k1
            endif
        enddo

        poly_order_added(k2) = 1
        order_add = int(10**(3-k2))
        Nord_prev = Nord_mod
        Nord_mod =  Nord_mod + order_add
        call ddecode(Nord_mod,pxm,pym,pzm)

        nrdofmQ = pxm * pym * pzm
        Mblock(mstep + 1) = nrdofmQ
        Ldwork = nrdofmQ

        allocate(Nextract_subson(nrdofmQ))
        Nextract_subson = ZERO
        call extraction_vector_new(Nord_prev,Nord_mod,Nord_glob,nrdofmQ,nrdofgQ,Nextract_prev_subson,Nextract_subson)
       
        ! call pbisolver3(mstep,Mblock,Ap_subson,Ldglob,Awork_subson,Ldwork,Nextract_subson,zbload_subson,Bwork_subson,Nrhs)
        
        ! call fine_to_subson_projection_error(kref_loc,Bwork_subson(1:nrdofmQ),Nextract_subson,subson_overlap, &
        !                                         nrdofmQ,nrdofgQ,Mdle,proj_error)

        do l = 1,nrdofmQ
        
            Bwork_subson(l) = zbload_subson(Nextract_subson(l))/Ap_subson(Nextract_subson(l),Nextract_subson(l))
    
        enddo

        call fine_to_subson_projection_error_new(kref_loc,Bwork_subson(1:nrdofmQ),Nextract_subson,subson_overlap,&
                                                    nrdofmQ,nrdofgQ,Mdle,proj_error,&
                                                    weights_fine_store,quad_point_store,nint_pp_store,&
                                                    shap3DQ_fine_store,shap3DQ_coarse_store)


        Nextract_prev_subson = ZERO
        Nextract_prev_subson(1:nrdofmQ) = Nextract_subson(1:nrdofmQ)
        deallocate(Nextract_subson)

        error_rate_aniso(k) = (base_error - proj_error)/(abs(real(nrdofmQ,8) - real(nrdof_org,8)))
        poly_flag_store(k) = Nord_mod

        mstep = mstep + 1
        coeff_sum = ZERO

        call coeff_polynomial_aniso(Nord_mod,Nextract_prev_subson(1:nrdofmQ),etype,Bwork_subson(1:nrdofmQ),coeff_sum)
        coeff_max = ZERO

    enddo

    rate_p = maxval(error_rate_aniso)
    max_loc = maxloc(error_rate_aniso)
    Polyflag = poly_flag_store(max_loc(1))




endif







end subroutine poly_adap_subson_new