!-----------------------------------------------------------------------
! routine name - project_p
!---------------------------------------------------------------------

!last revision: 18 Nov 2022

!Purpose: projects the fine mesh to coarse element and p enriched choices

!arguments
! in:
!     Mdle : index for element selected for refinement on the coarse mesh.
!     flag_pref : flag stating that the element is p-refined or not during fine mesh formation
    
! out: 

!     Error_org : will contain the projection error for the coarse element
!     rate_p : will contain the max rate obtained via the anisotropic p-refinement.
!     Poly_flag : will contain the associated polynomial refinement flag
!------------------------------------------------------------------------------------


subroutine project_p_linear(Mdle,flag_pref_loc, Error_org,rate_p,Poly_flag)

    use control
    use data_structure3D
    use element_data
    use parametersDPG
    use physics
    use mpi_param, only: RANK
 !
    implicit none

    integer,    intent(in)  :: Mdle
    integer,    intent(in)  :: flag_pref_loc ! indicator whether element is p-refined while producing the hp fine mesh
    real(8),    intent(out) :: Error_org !will store the projection error for original element
    real(8),    intent(out) :: rate_p   !will store the error drop rate for optimal p-refined element
    integer,    intent(out) :: Poly_flag !will store the flag for suggested p-refinement


    real(8), allocatable    :: Ap(:)
    real(8), allocatable    :: zbload(:,:)
    integer, allocatable    :: Mdle_sons(:)
    ! real(8), dimension(MAXbrickQQ*(MAXbrickQQ+1)/2) :: Ap_copy

    
!..to store the coefficients of the basis functions of interpolant on original element after projection from fine to coarse element.
    real(8) :: zdofQ(MAXEQNQ,MAXbrickQQ)
!..to store extracted coefficeints from 
    real(8) :: zdofQ_pp(MAXEQNQ,MAXbrickQQ), zdofH_pp(MAXEQNH,MAXbrickHH),zdofE_pp(MAXEQNE,MAXbrickEE), &
               zdofV_pp(MAXEQNV,MAXbrickVV)

!..to store the values of the solution for the original and sons respectively               
    real(8) :: zvalQ(MAXEQNQ), zvalQpp(MAXEQNQ)

! variables required for coarse element    
!.. number of shape functions in trial space for original element (in case of ultra weak formulation: only nrdofQ is useful)
    integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ  
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

! .. LAPACK SOLVER
    character uplo
    character(len=4) :: etype,ftype

!..element order, orientation for edges and faces of the original element
    integer :: norder(19), norient_edge(12), norient_face(6)
    integer :: nint
!..element order, orientation for edges and faces for enriched sons
    integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
    integer :: nint_pp, nrdof
!.. variables needed for telescopic solver
    integer :: px,py,pz
    integer :: nrdofmQ, nrdofgQ ! L2 shape functions for mth step and global system with order p+1 in x,y,z repectively
    integer :: Nord_org,Nord_glob,Nord_prev !ndof corresponding to nrdofmQ and nrdofgQ respectively
    integer,dimension(50)  :: Mblock
    integer :: mstep
    real(8), allocatable    :: Awork(:,:), Bwork(:,:)
    real(8), allocatable    :: proj_error_lvl(:) 
    real(8), allocatable    :: proj_coeff(:,:)
    integer, allocatable    :: Nextract(:),Nord_lvl(:)
    integer, allocatable    :: Nextract_prev(:) !will hold prevvious Nextract during telescopic solves
    integer, allocatable    :: Nextract_save(:,:)
    integer :: Ldglob,Ldwork,Nrhs
    real ::  check 
!.. variables needed for polynomial anisotropy computations
    real(8) :: coeff_sum(3),coeff_sum_net(3)
    real(8) :: coeff_max
    real(8) :: error_rate_aniso(3)
    integer :: poly_order_added(3)
    integer :: poly_flag_store(3)
    integer :: max_loc(1)
!..auxiliary variables
    real(8) :: wa,weight,rjac,bjac,fval,proj_error, proj_error_base,proj_error_net
    real(8) :: ordselect,error_tmp,error_tmp_prev
    integer :: k1,k2,l,iflag,k,is,mdle_fine,j,order_add,Nord_mod, first_son,nrdofmQ_base,nrdof_org 
    integer :: idec,iattr,icomp,ibeg
    integer :: var_loc !location of L2 variable in L2 solution array
    integer :: nr_mdle_sons ! number of sons of a refined element
    real(8) :: q1,q2 !shape functions values at gauss points
    real(8) :: q !shape functions while computing load vector
    real(8) :: MPI_Wtime,start_time,end_time, timer_a, timer_b,timer_c,timer_d,timer_e,timer_f, net_time_a, net_time_b,net_time_c
    real(8), allocatable :: weights_fine_store(:,:,:), quad_point_store(:,:,:)
    real(8), allocatable :: shap3DQ_fine_store(:,:,:),shap3DQ_coarse_store(:,:,:)
    integer, allocatable :: nint_pp_store(:)
    !..function for vector storage for symmertic L2 gram matrix in lower triangular format
    !integer :: nk
    !nk(k1,k2,nrdofQ) = nrdofQ * (k2-1) + k1 - k2*(k2-1)/2

    !allocating memory to the p+1 order projection matrix
    etype = NODES(Mdle)%type
    call find_order(Mdle, norder)
    call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)

    allocate(Ap(nrdofQ))
    allocate(zbload(nrdofQ,MAXEQNQ))


    if(etype .eq. 'mdlb') then
        nr_mdle_sons = 8
    endif

    allocate(Mdle_sons(nr_mdle_sons))

    first_son = NODES(Mdle)%first_son

    do is = 1,nr_mdle_sons
       Mdle_sons(is) = first_son + is - 1
    enddo

!.. looping over sons to compute the load vector of projection from fine to coarse
    Ap = ZERO
    zbload = ZERO

    allocate(weights_fine_store(MAXNINT3ADD,2,nr_mdle_sons))
    allocate(quad_point_store(3,MAXNINT3ADD,nr_mdle_sons))
    allocate(nint_pp_store(nr_mdle_sons))
    allocate(shap3DQ_fine_store(MAXbrickQQ,MAXNINT3ADD,nr_mdle_sons))
    allocate(shap3DQ_coarse_store(MAXbrickQQ,MAXNINT3ADD,nr_mdle_sons))

    weights_fine_store = ZERO
    quad_point_store = ZERO
    nint_pp_store  =ZERO
    shap3DQ_coarse_store = ZERO
    shap3DQ_fine_store = ZERO

    timer_e = MPI_Wtime()

    do is = 1,nr_mdle_sons
        mdle_fine = Mdle_sons(is)
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
        nint_pp_store(is) = nint_pp
        ! weights_fine_store(1:nint_pp,is) = waloc
        ! quad_point_store(1:3,1:nint_pp,is) = xiloc(1:3,1:nint_pp)

        !extract the coefficeints of the fine grid solution for the is^th son
        ! call solelm(mdle_fine,zdofH_pp,zdofE_pp,zdofV_pp,zdofQ_pp)
        call solelm_L2(mdle_fine,zdofQ_pp)

        do l = 1,nint_pp
            !..coordinates and weight of this integration point
            xi(1:3)=xiloc(1:3,l)
            quad_point_store(1:3,l,is) = xi(1:3)
            wa=waloc(l)
            !  ...H1 shape functions (for geometry)
            call shape3DH(etype,xi,norder_pp,norient_edge_pp,norient_face_pp, nrdof,shapH,gradH)
            !  ...L2 shape function calls
            call shape3DQ(etype,xi,norder_pp, nrdof,shapQ)

            shap3DQ_fine_store(1:nrdofQ_pp,l,is) = shapQ(1:nrdofQ_pp)
            !  ...geometry map
            call geom3D(mdle_fine,xi,xnod_pp,shapH,gradH,nrdofH_pp, x,dxdxi,dxidx,rjac,iflag)
 
            weight = rjac*wa
            weights_fine_store(l,1,is) = weight
            weights_fine_store(l,2,is) = rjac
            ! check = check + weight
            zvalQpp = ZERO

            !modify this by adding a loop over all L2 variables for cumlative adaptation

            do iattr = 1,NR_PHYSA
                if(DTYPE(iattr) .eq. 'discon') then
                    ibeg = ADRES(iattr)
                    do icomp = 1, NR_COMP(iattr)
                        var_loc = ibeg + icomp
                        do k = 1,nrdofQ_pp
                            q = shapQ(k)/rjac
                            zvalQpp(var_loc) = zvalQpp(var_loc) + zdofQ_pp(var_loc,k) * q
                        enddo
                    enddo
                endif
            enddo

            !calling the map between son's master element and coarse element master element
            call fine_to_coarse_gp_map(is,xi,xis,etype)
            
            shapQ = ZERO
            call shape3DQ(etype,xis,norder_pp, nrdof,shapQ)

            shap3DQ_coarse_store(1:nrdofQ_pp,l,is) = shapQ(1:nrdofQ_pp)
            !scaling the jacobian for isotropic refinement of coarse element
            rjac = rjac * real(nr_mdle_sons,8)

            !contribution to p+1 coarse projection matrix using subelement integration
            do k1 = 1,nrdofQ_pp
                q1 = shapQ(k1)/rjac
                q2 = shapQ(k1)/rjac
                    ! k = nk(k1,k2,nrdofQ_pp)
                Ap(k1) = Ap(k1) + weight * q1 * q2 
            
            enddo

            !computing the load vector contributions for each son
            do iattr = 1,NR_PHYSA
                if(DTYPE(iattr) .eq. 'discon') then
                    ibeg = ADRES(iattr)
                    do icomp = 1,NR_COMP(iattr)
                        var_loc = ibeg + icomp
                        
                        do k = 1,nrdofQ_pp
                            q = shapQ(k)/rjac
                            zbload(k,var_loc) = zbload(k,var_loc) + weight *  q * zvalQpp(var_loc)
                        enddo
                        
                    enddo
                endif
            enddo
            

        enddo

    enddo
    timer_f = MPI_Wtime()


    ! dont need telescopic solver in case of linear element
    net_time_a = 0.d0
    net_time_b = 0.d0
    net_time_c = 0.d0

    Nord_glob = norder_pp(19)
    call ddecode(Nord_glob,px,py,pz)

    if(etype .eq. 'mdlb') then
        nrdofgQ = px * py * pz
    endif


    allocate(Nextract_prev(nrdofgQ))
    allocate(Bwork(nrdofgQ,MAXEQNQ))
    Bwork = ZERO
    Nrhs = 1

    if(flag_pref_loc .eq. 1) then !if p-refinement done during hp-fine mesh generation
        Nord_org = Nord_glob - 111
    else
        Nord_org = Nord_glob
    endif


    call ddecode(Nord_org,px,py,pz)
    ! write(*,*) "p selec = ",px,py,pz

    if(etype .eq. 'mdlb') then
        nrdofmQ = px * py * pz
    endif

    nrdof_org = nrdofmQ ! number of dofs for original order before iso hp-refinement
    allocate(Nextract(nrdofmQ))

    Nextract = ZERO
    Nextract_prev = ZERO

    timer_a = MPI_Wtime()
    ! call extraction_vector(Nord_org,Nord_glob,nrdofmQ,nrdofgQ,Nextract)
    call extraction_vector_new(0,Nord_org,Nord_glob,nrdofmQ,nrdofgQ,Nextract_prev,Nextract)
    
    timer_b = MPI_Wtime()
    !diagonal solve due to L2 orthogonality of basis functions on linear elements
    do iattr = 1,NRQVAR
        do l = 1,nrdofmQ
            
            Bwork(l,iattr) = zbload(Nextract(l),iattr)/Ap(Nextract(l))

        enddo
    
    enddo

    
    Nextract_prev(1:nrdofmQ) = Nextract(1:nrdofmQ)

    timer_c = MPI_Wtime()
    ! call fine_to_coarse_projection_error(nr_mdle_sons,Bwork(1:nrdofmQ),Nextract,nrdofmQ,nrdofgQ,Mdle,Mdle_sons,proj_error)
    
    proj_error_net = 0.d0
    do iattr = 1,NRQVAR
        
        call fine_to_coarse_projection_error_new(nr_mdle_sons,Bwork(1:nrdofmQ,iattr),Nextract,nrdofmQ, &
                                                nrdofgQ,Mdle,iattr,Mdle_sons,proj_error,&
                                                weights_fine_store,quad_point_store,nint_pp_store,&
                                                shap3DQ_fine_store,shap3DQ_coarse_store)
       
        proj_error_net = proj_error_net + proj_error

    enddo

    timer_d = MPI_Wtime()
    ! write(*,*) "time is  = ", timer_b - timer_a, timer_c - timer_b,timer_d - timer_c
    net_time_a = net_time_a + timer_b - timer_a
    net_time_b = net_time_b + timer_c - timer_b
    net_time_c = net_time_c + timer_d - timer_c


    deallocate(Nextract)
    Error_org = proj_error_net
    ! storing the coefficients of polynomial along x,y, and z


    coeff_sum_net = ZERO

    do iattr = 1,NRQVAR

        coeff_sum = ZERO

        call coeff_polynomial_aniso(Nord_org,Nextract_prev(1:nrdofmQ),etype,Bwork(1:nrdofmQ,iattr),coeff_sum)

        coeff_sum_net(1:3) = coeff_sum_net(1:3) + coeff_sum(1:3)
    enddo

    !Anisotropic polynomial order increase for linear element

    if(etype .eq. 'mdlb') then

        error_rate_aniso = ZERO
        poly_order_added = ZERO

        coeff_max = 0.d0
        k2 = 0
        Nord_mod = Nord_org
        mstep = 2
        poly_flag_store = ZERO

        if(flag_pref_loc .eq. 1) then

            do k = 1,3

                do k1 = 1,3
                    if((coeff_max .le. coeff_sum_net(k1)) .and. (poly_order_added(k1) .eq. 0)) then
                        coeff_max = coeff_sum_net(k1)
                        k2 = k1
                    endif
                enddo

                poly_order_added(k2) = 1

                order_add = int(10**(3-k2))
                Nord_prev = Nord_mod
                Nord_mod =  Nord_mod + order_add
                call ddecode(Nord_mod,px,py,pz)

                nrdofmQ = px * py * pz
                allocate(Nextract(nrdofmQ))
                Nextract = ZERO

                timer_a = MPI_Wtime()
                call extraction_vector_new(Nord_prev,Nord_mod,Nord_glob,nrdofmQ,nrdofgQ,Nextract_prev,Nextract)

                !diagonal solve
                timer_b = MPI_Wtime()

                do iattr = 1,NRQVAR
                    do l = 1,nrdofmQ
                        
                        Bwork(l,iattr) = zbload(Nextract(l),iattr)/Ap(Nextract(l))
            
                    enddo
                
                enddo
                timer_c = MPI_Wtime()

                ! call fine_to_coarse_projection_error(nr_mdle_sons,Bwork(1:nrdofmQ),Nextract,nrdofmQ,nrdofgQ,Mdle,Mdle_sons,proj_error)
                proj_error_net = 0.d0
                do iattr = 1,NRQVAR
                    
                    call fine_to_coarse_projection_error_new(nr_mdle_sons,Bwork(1:nrdofmQ,iattr),Nextract,nrdofmQ, &
                                                            nrdofgQ,Mdle,iattr,Mdle_sons,proj_error,&
                                                            weights_fine_store,quad_point_store,nint_pp_store,&
                                                            shap3DQ_fine_store,shap3DQ_coarse_store)
                   
                    proj_error_net = proj_error_net + proj_error
            
                enddo

                timer_d = MPI_Wtime()

                net_time_a = net_time_a + timer_b - timer_a
                net_time_b = net_time_b + timer_c - timer_b
                net_time_c = net_time_c + timer_d - timer_c

                Nextract_prev = ZERO
                Nextract_prev(1:nrdofmQ) = Nextract(1:nrdofmQ)

                deallocate(Nextract)

                error_rate_aniso(k) = (Error_org - proj_error_net)/abs(real(NRQVAR * nrdofmQ,8) - real(NRQVAR * nrdof_org,8))

                poly_flag_store(k) = Nord_mod

                mstep = mstep + 1


                coeff_sum_net = ZERO
                do iattr = 1,NRQVAR
            
                    coeff_sum = ZERO
            
                    call coeff_polynomial_aniso(Nord_org,Nextract_prev(1:nrdofmQ),etype,Bwork(1:nrdofmQ,iattr),coeff_sum)
            
                    coeff_sum_net(1:3) = coeff_sum_net(1:3) + coeff_sum(1:3)
                enddo



                coeff_max = ZERO

            enddo    

            rate_p = maxval(error_rate_aniso)
            max_loc = maxloc(error_rate_aniso)
            Poly_flag = poly_flag_store(max_loc(1))
        
        else

            rate_p = 0.0
            Poly_flag = Nord_org

        endif
    endif
        !write(*,*) "time is  = ",net_time_a,net_time_b,net_time_c,timer_f - timer_e
    
    deallocate(weights_fine_store)
    deallocate(quad_point_store)
    deallocate(nint_pp_store)
    deallocate(shap3DQ_fine_store)
    deallocate(shap3DQ_coarse_store)

end subroutine project_p_linear
