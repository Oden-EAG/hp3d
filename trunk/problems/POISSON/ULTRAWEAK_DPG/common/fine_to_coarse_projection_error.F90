

!-----------------------------------------------------------------------
! routine name - fine_to_coarse_projection_error
!---------------------------------------------------------------------

! last revision: 7 Dec 2022

! arguments:

!in :
!     nr_mdle_sons: number of sons
!     coeff   : coefficeints of the basis functions on coarse element after projection
!     Nextract    : extraction vector from extraction_vector.F90
!     nrdofmQ     : number of Dofs on the coarse mesh for the correspoding order
!     nrdofgQ     : number of dofs on fine grid for the sons.
!     Mdle        : index of the coarse element.
!     Mdle_sons   : sons of the coarse element on the fine grid
! out:
!     proj_error  : projection error

!-----------------------------------------------------------------------------








subroutine fine_to_coarse_projection_error(nr_mdle_sons,coeff,Nextract,nrdofmQ,nrdofgQ,Mdle,Mdle_sons,proj_error)

    use control
    use data_structure3D
    use element_data
    use parametersDPG

    implicit none
    
    integer,    intent(in)  :: nrdofmQ
    integer,    intent(in)  :: nrdofgQ
    integer,    intent(in)  :: nr_mdle_sons
    integer,    intent(in)  :: Mdle
    integer, dimension(nr_mdle_sons),       intent(in) :: Mdle_sons
    real(8), dimension(nrdofmQ), intent(in) :: coeff
    integer, dimension(nrdofmQ), intent(in) :: Nextract

    real(8),    intent(out) :: proj_error

!..to store the coefficients of the basis functions of interpolant on original element after projection from fine to coarse element.
    real(8) :: zdofQ(MAXbrickQQ)
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

!..geometry dof of original element
   real(8) :: xnod(3,MAXbrickHH)
!..geometry dof of sons
   real(8) :: xnod_pp(3,MAXbrickHH)

    character(len=4) :: etype,ftype

!..element order, orientation for edges and faces of the original element
    integer :: norder(19), norient_edge(12), norient_face(6)
    integer :: nint, nrdof
!..element order, orientation for edges and faces for enriched sons
    integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
    integer :: nint_pp

    !..auxiliary variables
    real(8) :: wa,weight,rjac,bjac,fval
    integer :: k1,k2,l,iflag,k,is,mdle_fine,j,order_add,Nord_mod    
    integer :: idec
    real(8) :: q !shape functions while computing load vector

    etype = NODES(Mdle)%type
    call find_order(Mdle, norder)
    call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)

    ! if(etype .eq. 'mdlb') then
    !     nr_mdle_sons = 8
    ! endif

    proj_error = 0.d0

    do is = 1,nr_mdle_sons

        mdle_fine = Mdle_sons(is)
        etype = NODES(mdle_fine)%type
        call find_order(mdle_fine, norder_pp)
        call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
        ! write(*,*) nrdofH_pp
        call find_orient(mdle_fine, norient_edge_pp,norient_face_pp)
        call nodcor(mdle_fine, xnod_pp)
        INTEGRATION = 1
        call set_3D_int_DPG(etype,norder_pp,norient_face_pp, nint_pp,xiloc,waloc)
        !extract the coefficeints of the fine grid solution for the is^th son
        call solelm(mdle_fine,zdofH_pp,zdofE_pp,zdofV_pp,zdofQ_pp)

        do l = 1,nint_pp
            !..coordinates and weight of this integration point
            xi(1:3)=xiloc(1:3,l)
            wa=waloc(l)
            !  ...H1 shape functions (for geometry)
            call shape3DH(etype,xi,norder_pp,norient_edge_pp,norient_face_pp, nrdof,shapH,gradH)
            !  ...L2 shape function calls
            call shape3DQ(etype,xi,norder_pp, nrdof,shapQ)
            !  ...geometry map
            call geom3D(mdle_fine,xi,xnod_pp,shapH,gradH,nrdofH_pp, x,dxdxi,dxidx,rjac,iflag)

            weight = rjac*wa
            ! check = check + weight
            zvalQpp = ZERO

            do k = 1,nrdofQ_pp
                q = shapQ(k)/rjac
                zvalQpp(1) = zvalQpp(1) + zdofQ_pp(1,k) * q
            enddo

                        !calling the map between son's master element and coarse element master element
            call fine_to_coarse_gp_map(is,xi,xis,etype)

            call shape3DQ(etype,xis,norder_pp, nrdof,shapQ)
                        !scaling the jacobian for isotropic refinement of coarse element
            rjac = rjac * real(nr_mdle_sons,8)

            !recontructing the projection at the transformed gauss point for
            !the transformed coarse elements.

            zvalQ = ZERO
            do k = 1,nrdofmQ
                q = shapQ(Nextract(k))/rjac
                zvalQ(1) = zvalQ(1) + coeff(k) * q
            enddo

            proj_error = proj_error + weight * (zvalQ(1) - zvalQpp(1))**2

        enddo

    enddo

    ! proj_error = sqrt(proj_error)

    
end subroutine fine_to_coarse_projection_error





!-----------------------------------------------------------------------
! routine name - fine_to_subson_projection_error
!---------------------------------------------------------------------

! last revision: 11 Dec 2022

! arguments:

!in :
!     coeff   : coefficeints of the basis functions on coarse element after projection
!     Nextract    : extraction vector from extraction_vector.F90
!     overlap     : small array with overlapping fine_sons 
!      
!     nrdofmQ     : number of Dofs on the coarse mesh for the correspoding order
!     nrdofgQ     : number of dofs on fine grid for the sons.
!     Mdle        : index of the coarse element.
!     Mdle_sons   : sons of the coarse element on the fine grid
! out:
!     proj_error  : projection error

!-----------------------------------------------------------------------------

subroutine fine_to_subson_projection_error(kref,coeff,Nextract,overlap,nrdofmQ,nrdofgQ, &
                                            Mdle,proj_error)
    use control
    use data_structure3D
    use element_data
    use parametersDPG

    implicit none

    integer,    intent(in)  ::  kref
    integer,    intent(in)  ::  nrdofmQ
    integer,    intent(in)  ::  nrdofgQ
    real(8), dimension(nrdofmQ), intent(in) :: coeff
    integer, dimension(nrdofmQ), intent(in) :: Nextract
    
    integer,    intent(in)  ::  Mdle
    integer,    intent(in)  ::  overlap(*)

    real(8),    intent(out)  ::  proj_error

    integer :: norder_pp(19), norient_edge_pp(12), norient_face_pp(6)
    integer :: nrdofH_pp, nrdofE_pp, nrdofV_pp, nrdofQ_pp  
    !..to store extracted coefficeints from 
    real(8) :: zdofQ_pp(MAXEQNQ,MAXbrickQQ), zdofH_pp(MAXEQNH,MAXbrickHH),zdofE_pp(MAXEQNE,MAXbrickEE), &
                zdofV_pp(MAXEQNV,MAXbrickVV)

    !..to store the values of the solution for the original and sons respectively               
    real(8) :: zvalQ(MAXEQNQ), zvalQpp(MAXEQNQ)
!.. L2 shape functions
    real(8) :: shapQ(MAXbrickQQ)
    !..H1 shape functions
    real(8) :: shapH(MAXbrickHH), gradH(3,MAXbrickHH)
    
    !..geometry dof of sons
    real(8) :: xnod_pp(3,MAXbrickHH)
    !..3D quadrature data
    real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
    
    !..geometry
    real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3),xis(3)
    real(8) :: dxidt(3,2), dxdt(3,2), t(2)

    integer :: hx,hy,hz !refinement indicators in x,y,z directions

    character(len=4) :: etype
    character(len=4) :: etype_coarse
    !auxiliary variables
    real(8) :: rjac,wa,weight,q,check
    integer :: first_son,nr_subsons,nr_mdle_sons,iss,iflag,nrdof
    integer :: is, mdle_fine,overlap_count,nint_pp
    integer :: k,l
    integer,allocatable :: Mdle_olp_fine(:)

    etype_coarse = NODES(Mdle)%type

    if(etype_coarse .eq. 'mdlb') then

        nr_mdle_sons = 8

        call ddecode(kref,hx,hy,hz)
        nr_subsons = 2**(hx+hy+hz)

        overlap_count = nr_mdle_sons/nr_subsons

        allocate(Mdle_olp_fine(overlap_count))
        
        first_son = NODES(Mdle)%first_son

        do is = 1,overlap_count

            Mdle_olp_fine(is) = first_son + overlap(is) - 1

        enddo

        proj_error = 0.d0

        ! if(kref .eq. 100) then 
        !     write(*,*) Mdle_olp_fine
        ! endif

        check = 0.d0

        do is = 1,overlap_count

            mdle_fine = Mdle_olp_fine(is)
            etype = NODES(mdle_fine)%type
            call find_order(mdle_fine, norder_pp)
            call celndof(etype,norder_pp, nrdofH_pp,nrdofE_pp,nrdofV_pp,nrdofQ_pp)
            ! write(*,*) nrdofH_pp
            call find_orient(mdle_fine, norient_edge_pp,norient_face_pp)
            call nodcor(mdle_fine, xnod_pp)
    
            INTEGRATION = 1

            call set_3D_int_DPG(etype,norder_pp,norient_face_pp, nint_pp,xiloc,waloc)
            !extract the coefficeints of the fine grid solution for the is^th son
            call solelm(mdle_fine,zdofH_pp,zdofE_pp,zdofV_pp,zdofQ_pp)

            do l = 1,nint_pp

            !..coordinates and weight of this integration point
                xi(1:3)=xiloc(1:3,l)
                wa=waloc(l)

                call shape3DH(etype,xi,norder_pp,norient_edge_pp,norient_face_pp, nrdof,shapH,gradH)
                !  ...L2 shape function calls
                call shape3DQ(etype,xi,norder_pp, nrdof,shapQ)
                !  ...geometry map
                call geom3D(mdle_fine,xi,xnod_pp,shapH,gradH,nrdofH_pp, x,dxdxi,dxidx,rjac,iflag)
                
                weight = rjac*wa

                zvalQpp = ZERO

                do k = 1,nrdofQ_pp
                    q = shapQ(k)/rjac
                    zvalQpp(1) = zvalQpp(1) + zdofQ_pp(1,k) * q
                enddo

                call fine_to_subson_gp_map(etype_coarse,kref,overlap(is),xi,xis)

                shapQ = ZERO

                call shape3DQ(etype,xis,norder_pp,nrdof,shapQ)

                !scaling the jacobian for isotropic refinement of coarse element
                rjac = rjac * real(overlap_count,8)

                !recontructing the projection at the transformed gauss point for
                !the transformed coarse elements.

                zvalQ = ZERO
                do k = 1,nrdofmQ
                    q = shapQ(Nextract(k))/rjac
                    zvalQ(1) = zvalQ(1) + coeff(k) * q
                enddo

                proj_error = proj_error + weight * (zvalQ(1) - zvalQpp(1))**2
                ! proj_error = proj_error + weight * (1.d0)**2
                ! check = check + weight * (1.d0)**2

            enddo

        enddo
        ! if(kref .eq. 100) then 
            ! write(*,*) "check = ", proj_error
        ! endif

    endif

end subroutine fine_to_subson_projection_error


