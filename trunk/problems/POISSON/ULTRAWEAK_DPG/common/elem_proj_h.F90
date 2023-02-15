!-----------------------------------------------------------------------
! routine name - elem_proj_h
!---------------------------------------------------------------------

!last revision: 7 Dec 2022

!Purpose: projects the fine mesh to coarse element for h-refinement option

!arguments
! in:
!     Mdle : index for element selected for refinement on the coarse mesh.
!     Kref_loc : h-refinement option

    
! out: 

!     error_opt: will compute the final error for the given h-refinement
!                after following the largest subelement error path.
!------------------------------------------------------------------------------------



subroutine elem_proj_h(Mdle,flag_pref_loc,kref_loc,error_org,error_opt,g_rate_max,rate_hcomp,Nord_href)

    use control
    use data_structure3D
    use element_data
    use parametersDPG
    use MPI 
    !
    implicit none

    integer,    intent(in)  :: Mdle
    integer,    intent(in)  :: flag_pref_loc
    integer,    intent(in)  :: kref_loc
    real(8),    intent(in)  :: error_org
    real(8),    intent(out) :: g_rate_max
    real(8),    intent(out) :: rate_hcomp
    real(8),    intent(out) :: error_opt
    integer,    dimension(8), intent(out)   :: Nord_href
    character(len=4) :: etype
!..element order, orientation for edges and faces of the original element
    
    integer :: norder(19)
    integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
!.. number of shape functions in trial space for enriched sons (its same for all of them)
    integer :: nrdofH_pp, nrdofE_pp, nrdofV_pp, nrdofQ_pp  
!.. L2 shape functions
    real(8) :: shapQ(MAXbrickQQ)
!..H1 shape functions
    real(8) :: shapH(MAXbrickHH), gradH(3,MAXbrickHH)
!..3D quadrature data
    real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)

!..geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3),xis(3)
   real(8) :: dxidt(3,2), dxdt(3,2), t(2)

!..geometry dof of sons
   real(8) :: xnod_pp(3,MAXbrickHH)


!..to store extracted coefficeints from 
   real(8) :: zdofQ_pp(MAXEQNQ,MAXbrickQQ), zdofH_pp(MAXEQNH,MAXbrickHH),zdofE_pp(MAXEQNE,MAXbrickEE), &
   zdofV_pp(MAXEQNV,MAXbrickVV)

!..to store the values of the solution for the original and sons respectively               
   real(8) :: zvalQ(MAXEQNQ), zvalQpp(MAXEQNQ)

    integer :: nr_subsons,hx,hy,hz
    integer :: pxm,pym,pzm ! max poly order in each direction
    real(8),    allocatable :: subsons_Ap(:,:,:) !projection matrix for the child of the coarse elements 
                                               !produced by anisotropic h ref
    real(8),    allocatable :: subsons_zbload(:,:) !load vector for the projection problem for the child of the coarse element
                                                    !produced by anisotropic h ref
    integer,    allocatable :: subsons_overlap(:,:)! stores the overlap of the fine grid with subsons
    integer,    allocatable :: Mdle_sons(:)
    integer,    allocatable :: Nord_org_subsons(:)
    integer,    allocatable :: Nextract(:)


    real(8),    allocatable :: Awork(:,:),Bwork(:),Ap(:,:),zbload(:)
    real(8),    allocatable :: subsons_Awork(:,:,:),subsons_Bwork(:,:)
    integer :: Ldglob,Ldwork,Nrhs 
    real(8) :: proj_error_subson, proj_error
    integer,dimension(50)  :: Mblock
    integer,dimension(50)  :: dof_diff
    real(8),dimension(50)  :: error_rate


    integer :: nrdofgQ,nrdofmQ,mstep, nrdof_org
    integer :: nr_mdle_sons,mdle_fine,first_son
    integer :: nint_pp,nrdof,Nord_glob, Nord_org, Polyflag
    integer, allocatable    :: ipi(:),Nord_mep(:,:),Nord_old(:)
    real(8), allocatable    :: error_subsons(:) 
    
    integer, allocatable    :: subsons_Nextract(:,:)
    integer, allocatable    :: subsons_Mblock(:,:)
    integer, allocatable    :: subsons_mstep(:),Nord_max(:)

    real(8) :: q1,q2 !shape functions for computing the projection matrix
    real(8) :: q !shape functions for the load vector computation in the projection problem
    real(8) :: max_error_subson, ratio_mep

    !Auxiliary variable
    integer :: j,iss,is,l,k,iflag,k1,k2,iso_p,iss_max,nrdof_tmp
    integer :: Nref, local_order_check
    integer :: pxc,pyc,pzc, pxg,pyg,pzg
    real(8) :: wa, weight,rjac,g_rate_tmp, timer_a,timer_b,timer_c
    real(8) :: var_a, var_b, var_c !multipurpose temporary variable
    integer,    allocatable :: subsons_Nextract_prev(:,:)

    etype = NODES(Mdle)%type
    call find_order(Mdle, norder)

    Nord_glob = norder(19) !interior order for L2 variables
    

    if(etype .eq. 'mdlb') then

        !dofs on the coarse element before hp refinement

        ! call ddecode(norder(19),pxm,pym,pzm)
        timer_a = MPI_Wtime()
        if(flag_pref_loc .eq. 1) then
            Nord_org = Nord_glob - 111
        else
            Nord_org = Nord_glob
        endif
    
        call ddecode(Nord_org,pxm,pym,pzm)
        nrdof_org = pxm * pym * pzm

        nr_mdle_sons = 8 !number of fine childs

        allocate(Mdle_sons(nr_mdle_sons))

        call ddecode(Nord_glob,pxg,pyg,pzg)
        nrdofgQ = pxg * pyg *pzg
        iso_p = MAX(pxg,pyg,pzg)

        call ddecode(kref_loc,hx,hy,hz)
        nr_subsons = 2**(hx+hy+hz)

        allocate(subsons_Ap(nrdofgQ,nrdofgQ,nr_subsons))
        subsons_Ap = ZERO

        allocate(subsons_zbload(nrdofgQ,nr_subsons))
        subsons_zbload = ZERO

        ! write(*,*) "number of subsons = ",nr_subsons,kref_loc

        allocate(subsons_overlap(nr_subsons,8/nr_subsons))

        subsons_overlap = ZERO

        select case(kref_loc)
            case(100)

                subsons_overlap(1,:) = (/1,4,5,8/)
                subsons_overlap(2,:) = (/2,3,6,7/)

            case(010)
                
                subsons_overlap(1,:) = (/1,2,5,6/)
                subsons_overlap(2,:) = (/3,4,7,8/)

            case(001)

                subsons_overlap(1,:) = (/1,2,3,4/)
                subsons_overlap(2,:) = (/5,6,7,8/)

            case(110)

                subsons_overlap(1,:) = (/1,5/)
                subsons_overlap(2,:) = (/2,6/)
                subsons_overlap(3,:) = (/3,7/)
                subsons_overlap(4,:) = (/4,8/) 

            case(101)

                subsons_overlap(1,:) = (/1,4/)
                subsons_overlap(2,:) = (/2,3/)
                subsons_overlap(3,:) = (/5,8/)
                subsons_overlap(4,:) = (/6,7/) 

            case(011)

                subsons_overlap(1,:) = (/1,2/)
                subsons_overlap(2,:) = (/3,4/)
                subsons_overlap(3,:) = (/5,6/)
                subsons_overlap(4,:) = (/7,8/) 
            
            case(111)
                do j = 1,nr_subsons

                    subsons_overlap(j,1) = j

                enddo

            case default
                write(*,*) " Not a valid kref"

        end select

        first_son = NODES(Mdle)%first_son
        do is = 1,nr_mdle_sons
           Mdle_sons(is) = first_son + is - 1
        enddo
        
        do iss = 1,nr_subsons
            
            do is = 1,8/nr_subsons

                mdle_fine = Mdle_sons(subsons_overlap(iss,is))
                etype = NODES(mdle_fine)%type

                call find_order(mdle_fine, norder_pp)
                call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
                ! write(*,*) nrdofH_pp
                call find_orient(mdle_fine, norient_edge_pp,norient_face_pp)
                call nodcor(mdle_fine, xnod_pp)

                INTEGRATION = 1
                ! xiloc = ZERO
                ! waloc = ZERO
                call set_3D_int_DPG(etype,norder_pp,norient_face_pp, nint_pp,xiloc,waloc)
                !extract the coefficeints of the fine grid solution for the is^th son
                ! call solelm(mdle_fine,zdofH_pp,zdofE_pp,zdofV_pp,zdofQ_pp)
                call solelm_L2(mdle_fine,zdofQ_pp)
            
                do l = 1,nint_pp

                    xi(1:3) = xiloc(1:3,l)
                    wa = waloc(l)
                                !  ...H1 shape functions (for geometry)
                    call shape3DH(etype,xi,norder_pp,norient_edge_pp,norient_face_pp, nrdof,shapH,gradH)
                    !  ...L2 shape function calls
                    call shape3DQ(etype,xi,norder_pp, nrdof,shapQ)
                    !  ...geometry map
                    call geom3D(mdle_fine,xi,xnod_pp,shapH,gradH,nrdofH_pp, x,dxdxi,dxidx,rjac,iflag)
 
                    weight = rjac*wa

                    zvalQpp = ZERO

                    !modify this by adding a loop over all L2 variables for cumlative adaptation
                    do k = 1,nrdofQ_pp
                        q = shapQ(k)/rjac
                        zvalQpp(1) = zvalQpp(1) + zdofQ_pp(1,k) * q
                    enddo

                    !calling the map between son's master element and coarse element master element
                    ! call fine_to_coarse_gp_map(subsons_overlap(iss,is),xi,xis,etype)
                    call fine_to_subson_gp_map(etype,kref_loc,subsons_overlap(iss,is),xi,xis)
                    ! xis = xi
                    shapQ = ZERO

                    call shape3DQ(etype,xis,norder_pp,nrdof,shapQ)


                    !scaling the jacobian for isotropic refinement of coarse element
                    rjac = rjac * real(nr_mdle_sons/nr_subsons,8)
                    
                    do k1 = 1,nrdofQ_pp
                        q1 = shapQ(k1)/rjac
                        do k2 = 1,nrdofQ_pp

                            q2 = shapq(k2)/rjac
                            
                            subsons_Ap(k1,k2,iss) = subsons_Ap(k1,k2,iss) + weight * q1 * q2

                        enddo
                    enddo

                    do k = 1,nrdofQ_pp

                        q = shapQ(k)/rjac
                        subsons_zbload(k,iss) = subsons_zbload(k,iss) + weight * q * zvalQpp(1)

                    enddo


                enddo

            enddo


        enddo

        timer_b = MPI_Wtime()
        !now solving the projection problem starting from base order on each subson by looping over them and calling
        ! calling telescopic solver.
        Nrhs = 1
        Ldglob = nrdofgQ  

        allocate(Ap(nrdofgQ,nrdofgQ))
        allocate(zbload(nrdofgQ))


        proj_error = 0.d0
        proj_error_subson = 0.d0

        allocate(Nord_org_subsons(nr_subsons))
        Nord_org_subsons = ZERO

        allocate(Nord_mep(nr_subsons,nr_subsons * iso_p))
        Nord_mep = ZERO
        Nord_mep(1:nr_subsons,1) = 222

        allocate(Nord_old(nr_subsons))
        Nord_old = ZERO

        allocate(Nord_max(nr_subsons))
        Nord_max = ZERO

        allocate(error_subsons(nr_subsons))
        error_subsons   = ZERO

        allocate(subsons_Mblock(50,nr_subsons))
        subsons_Mblock  = ZERO

        allocate(subsons_mstep(nr_subsons))
        subsons_mstep   = 1
        
        allocate(subsons_Nextract_prev(nr_subsons,nrdofgQ))
        subsons_Nextract_prev = ZERO

        allocate(subsons_Awork(nrdofgQ,nrdofgQ,nr_subsons))
        subsons_Awork = ZERO

        allocate(subsons_Bwork(nrdofgQ,nr_subsons))
        subsons_Bwork = ZERO
        
        dof_diff = ZERO
        error_rate = ZERO
        ! write(*,*) " I am here too "
 
            
        Nref = 1
        ratio_mep = 0.7d0 !ratio for selecting subsons to p-refine  
        g_rate_max = 0.d0 !maximum guranteed rate

        ! call ddecode(Nord_glob,pxg,pyg,pzg)

        do

            ! if(sum(Nord_mep(1:nr_subsons,Nref)) .gt. nr_subsons * Nord_glob) then
            !     exit
            ! endif

            ! write(*,*) "sum = ",sum(Nord_mep(1:nr_subsons,Nref)),nr_subsons * Nord_glob,nr_subsons, Nord_mep(1:nr_subsons,Nref)
            !solve over all subsons
            do iss = 1,nr_subsons
                if(Nord_old(iss) .ne. Nord_mep(iss,Nref)) then !only solve if the subson has order changed

                    allocate(Awork(nrdofgQ,nrdofgQ))
                    allocate(Bwork(nrdofgQ))
        
                    Awork = ZERO
                    Bwork = ZERO

                    !copying the projection matrix
                    Awork(1:nrdofgQ,1:nrdofgQ) = subsons_Awork(1:nrdofgQ,1:nrdofgQ,iss)
                    Bwork(1:nrdofgQ)           = subsons_Bwork(1:nrdofgQ,iss)

        
                    Nord_org_subsons(iss) = Nord_mep(iss,Nref)

                    
                    call ddecode(Nord_org_subsons(iss),pxm,pym,pzm)
                    nrdofmQ = pxm * pym * pzm

                    mstep = subsons_mstep(iss)

                    if(mstep .eq. 1) then
                        subsons_Mblock(mstep,iss) = 0
                    endif

                    Ldwork = nrdofmQ
                    subsons_Mblock(mstep + 1,iss) = nrdofmQ
                    
                    Mblock(1:50) = subsons_Mblock(1:50,iss)


                    allocate(Nextract(nrdofmQ))
                    Nextract = ZERO
        
                    Ap(1:nrdofgQ,1:nrdofgQ) = subsons_Ap(1:nrdofgQ,1:nrdofgQ,iss)
                    zbload(1:nrdofgQ) = subsons_zbload(1:nrdofgQ,iss)
        
                    ! call extraction_vector(Nord_org_subsons(iss),Nord_glob,nrdofmQ,nrdofgQ,Nextract)
                    call extraction_vector_new(Nord_old(iss),Nord_org_subsons(iss),Nord_glob,nrdofmQ,nrdofgQ,subsons_Nextract_prev(iss,:),Nextract)

                    call pbisolver3(mstep,Mblock,Ap,Ldglob,Awork,Ldwork,Nextract,zbload,Bwork,Nrhs)

                    subsons_Nextract_prev(iss,1:nrdofmQ) = Nextract(1:nrdofmQ)

                    proj_error_subson = 0.d0

                    call fine_to_subson_projection_error(kref_loc,Bwork(1:nrdofmQ),Nextract,subsons_overlap(iss,:), &
                                                        nrdofmQ,nrdofgQ,Mdle,proj_error_subson)
        
                    ! proj_error =  proj_error + proj_error_subson

                    error_subsons(iss) = proj_error_subson
                    
                    subsons_Awork(1:nrdofgQ,1:nrdofgQ,iss) = Awork(1:nrdofgQ,1:nrdofgQ)
                    subsons_Bwork(1:nrdofgQ,iss) = Bwork(1:nrdofgQ)

                    deallocate(Awork)
                    deallocate(Bwork)
                    deallocate(Nextract)

                    subsons_mstep(iss) = subsons_mstep(iss) + 1

                endif
            enddo

            proj_error = sum(error_subsons)
            
            !conditional max
            max_error_subson = 0.d0
            iss_max = 0
            Nord_old(1:nr_subsons) = Nord_mep(1:nr_subsons,Nref)
            

            nrdof_tmp = 0
            do iss = 1,nr_subsons

                if((error_subsons(iss) .ge. max_error_subson)   .and. (Nord_mep(iss,Nref) .lt. Nord_glob)) then

                    max_error_subson = error_subsons(iss)

                endif
                
                call ddecode(Nord_old(iss),pxm,pym,pzm)
                nrdof_tmp = nrdof_tmp + pxm*pym*pzm

            enddo

            dof_diff(Nref) = nrdofgQ - nrdof_tmp

            !increasing the order for max error in subsons with order  .le. Nord_glob

            do iss = 1,nr_subsons

                Nord_mep(iss,Nref + 1) = Nord_mep(iss,Nref)
                if((error_subsons(iss) .ge. ratio_mep * max_error_subson) .and. (Nord_mep(iss,Nref) .lt. Nord_glob)) then

                    ! Nord_mep(iss, Nref + 1) = Nord_mep(iss,Nref) + 111

                    call ddecode(Nord_mep(iss,Nref),pxc,pyc,pzc)
                    if(((pxg - pxc) .ge. 1) .and. ((pyg - pyc) .ge. 1) .and. ((pzg - pzc) .ge. 1)) then

                        allocate(Awork(nrdofgQ,nrdofgQ))
                        allocate(Bwork(nrdofgQ))
                        nrdofmQ = pxc * pyc * pzc

                        Awork = ZERO
                        Bwork = ZERO

                        Awork(1:nrdofgQ,1:nrdofgQ) = subsons_Awork(1:nrdofgQ,1:nrdofgQ,iss)
                        Bwork(1:nrdofgQ)           = subsons_Bwork(1:nrdofgQ,iss)

                        Ap(1:nrdofgQ,1:nrdofgQ) = subsons_Ap(1:nrdofgQ,1:nrdofgQ,iss)
                        zbload(1:nrdofgQ) = subsons_zbload(1:nrdofgQ,iss)

                        call poly_adap_subson(kref_loc,Mdle,nrdofgQ,nrdof_org,Nord_mep(iss,Nref),Nord_glob,subsons_overlap(iss,:),error_org, &
                                              Ap,zbload,Awork,Bwork,subsons_Nextract_prev(iss,1:nrdofmQ),Polyflag)
                                              
                        Nord_mep(iss, Nref + 1) = Polyflag

                        deallocate(Awork)
                        deallocate(Bwork)

                    else
                        
                        Nord_mep(iss, Nref + 1) = Nord_glob
                    
                    endif

                endif

            enddo

            ! g_rate_tmp = (log(error_org**2) - log(proj_error))/(log(real(nrdof_tmp,8)) - log(real(nrdof_org,8)))
            if(nrdof_tmp .ne. nrdof_org) then
                g_rate_tmp = (error_org - proj_error)/abs(real(nrdof_tmp,8) - real(nrdof_org,8))
                error_rate(Nref) = g_rate_tmp         
                ! write(*,*) " The rate is ",g_rate_tmp,Nref,Mdle,real(nrdof_tmp,8),real(nrdof_org,8),error_org,proj_error
                ! write(*,*) Nord_old(1:nr_subsons)
            endif

            if(g_rate_tmp .ge. g_rate_max) then

                g_rate_max = g_rate_tmp
                error_opt = proj_error
                Nord_max(1:nr_subsons) = Nord_old(1:nr_subsons) 

            endif

            ! error_opt = proj_error
            ! g_rate_ref = g_rate_tmp
            local_order_check = 0

            do iss = 1,nr_subsons

                call ddecode(Nord_mep(iss,Nref+1),pxc,pyc,pzc)
                
                if((pxc .gt. pxg) .or. (pyc .gt. pyg) .or. (pzc .gt. pzg)) then
                    local_order_check = local_order_check + 1
                endif

            enddo

            ! write(*,*) "local_order_check", local_order_check

            if(local_order_check .gt. 0) then
                exit
            endif


            local_order_check = 0

            do iss = 1,nr_subsons
                call ddecode(Nord_mep(iss,Nref),pxc,pyc,pzc)
                
                if((pxc .eq. pxg) .and. (pyc .eq. pyg) .and. (pzc .eq. pzg)) then
                    local_order_check = local_order_check + 1
                endif

            enddo

            ! write(*,*) "local_order_check", local_order_check

            if(local_order_check .eq. nr_subsons) then
                exit
            endif



            ! if(sum(Nord_mep(1:nr_subsons,Nref)) .eq. nr_subsons * Nord_glob) then !exit at the last step when every one reaches at the highest polynomial order
            !     exit
            ! endif

           

            
            Nref = Nref + 1

        enddo
        timer_c = MPI_Wtime()
        ! write(*,*) Nord_max

        Nord_href(1:nr_subsons) = Nord_max(1:nr_subsons)

        rate_hcomp = 0.d0

        do k = 1,Nref

            if((dof_diff(k)) .ge. 0) then

                    if(rate_hcomp .le. error_rate(k)) then
                        rate_hcomp = error_rate(k)
                    endif
            
            endif
        enddo

    endif

    if(kref_loc .eq. 111) then
        write(*,*) timer_b - timer_a,timer_c-timer_a

    endif



end subroutine elem_proj_h