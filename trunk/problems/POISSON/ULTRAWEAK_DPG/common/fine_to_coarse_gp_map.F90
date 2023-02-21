
!-----------------------------------------------------------------------
! routine name - fine_to_coarse_gp_map
!---------------------------------------------------------------------

! last revision: 7 Dec 2022

! arguments:

!in :
!   iel: index of the son
!   xi: gauss point correspoding to the son
!   xis: gauss point corresponding to the coarse element.
!   etype: element type

!-----------------------------------------------------------------------------

subroutine fine_to_coarse_gp_map(iel,xi,xis,etype)

implicit none

real(8), intent(in) :: xi(3)
integer, intent(in) :: iel
character(len=4), intent(in) :: etype

real(8), intent(out) :: xis(3)


select case(etype)
    case('mdlb')
        select case(iel)
            case(1)
                xis  = 0.5 * xi
            case(2)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 * xi(2)
                xis(3) = 0.5 * xi(3)
            case(3)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 * xi(3)
            case(4)
                xis(1) = 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 * xi(3)
            case(5)
                xis(1) = 0.5 * xi(1)
                xis(2) = 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case(6)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case(7)
                xis(1) = 0.5 + 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case(8)
                xis(1) = 0.5 * xi(1)
                xis(2) = 0.5 + 0.5 * xi(2)
                xis(3) = 0.5 + 0.5 * xi(3)
            case default
                write(*,*) "Wrong son"
        end select
    case default
        write(*,*) "not implemented yet"


end select






end subroutine fine_to_coarse_gp_map


subroutine  fine_to_subson_gp_map(etype,kref,sel,xi,xis)

implicit none
integer,    intent(in)  :: kref
integer,    intent(in)  :: sel
character(len=4), intent(in) :: etype

real(8), intent(in) :: xi(3)
real(8), intent(out) :: xis(3)


if(etype .eq. 'mdlb') then

    select case(kref)

        case(100)
            select case(sel)

                case(1,2)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,0.5d0,0.5d0,xi,xis)
                case(4,3)
                    call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,1.d0,0.5d0,0.5d0,xi,xis)
                case(5,6)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,1.d0,0.5d0,0.5d0,xi,xis)
                case(8,7)
                    call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.5d0,1.d0,0.5d0,0.5d0,xi,xis)

                case default
                    write(*,*) "wrong fine son"

            end select
        case(010)
            select case(sel)
                case(1,4)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,1.d0,0.5d0,xi,xis)
                case(2,3)
                    call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.d0,0.5d0,1.d0,0.5d0,xi,xis)
                case(5,8)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,0.5d0,1.d0,0.5d0,xi,xis)
                case(6,7)
                    call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.5d0,0.5d0,1.d0,0.5d0,xi,xis)
                case default
                    write(*,*) "wrong fine son"

            end select

        case(001)
            select case(sel)

                case(1,5)
                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case(2,6)
                    call singlefineson_to_subson_gp_map(0.5d0,0.d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case(3,7)
                    call singlefineson_to_subson_gp_map(0.5d0,0.5d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case(4,8)
                    call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,0.5d0,0.5d0,1.d0,xi,xis)
                case default
                    write(*,*) "wrong fine son"

            end select

        case(110)
            select case(sel)
                case(1,2,3,4)

                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,1.d0,0.5d0,xi,xis)

                case(5,6,7,8)

                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.5d0,1.d0,1.d0,0.5d0,xi,xis) 

            end select


        case(101)
            select case(sel)
                case(1,2,5,6)

                call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,1.d0,0.5d0,1.d0,xi,xis)

                case(4,3,8,7)

                call singlefineson_to_subson_gp_map(0.d0,0.5d0,0.d0,1.d0,0.5d0,1.d0,xi,xis) 

            end select

        case(011)
            select case(sel)
                case(1,4,5,8)

                    call singlefineson_to_subson_gp_map(0.d0,0.d0,0.d0,0.5d0,1.d0,1.d0,xi,xis)

                case(2,3,6,7)

                    call singlefineson_to_subson_gp_map(0.5d0,0.0d0,0.d0,0.5d0,1.d0,1.d0,xi,xis) 

            end select


        case(111)

             xis = xi



    end select


endif



end subroutine fine_to_subson_gp_map





subroutine singlefineson_to_subson_gp_map(ofx,ofy,ofz,fx,fy,fz,xi,xis)

    implicit none
    real(8), intent(in) :: ofx
    real(8), intent(in) :: ofy
    real(8), intent(in) :: ofz

    real(8), intent(in) :: fx
    real(8), intent(in) :: fy
    real(8), intent(in) :: fz

    real(8),dimension(3),   intent(in)    :: xi
    real(8),dimension(3),   intent(out)   :: xis

    xis(1) = ofx + fx * xi(1)
    xis(2) = ofy + fy * xi(2)
    xis(3) = ofz + fz * xi(3)


end subroutine singlefineson_to_subson_gp_map


subroutine subson_one_irregularity_map(etype,kref_intent,kref_appl,nr_sons_intent,pref_intent,nr_sons_appl,pref_appl)


    implicit none

    character(len=4), intent(in) :: etype
    integer,    intent(in)  :: kref_intent !intented h-ref
    integer,    intent(in)  :: kref_appl  !isoref subson index when iso ref is implemented in place of intended h-refinement
    integer,    intent(in)  :: nr_sons_intent 
    integer,    intent(in)  :: pref_intent(nr_sons_intent)
    integer,    intent(in)  :: nr_sons_appl
    integer,    intent(out) :: pref_appl(nr_sons_appl)

    
    integer :: i
    integer, allocatable :: el_pmap(:)
    


    if(etype .eq. 'mdlb') then


        allocate(el_pmap(8))

            !making an transfer map using an isotropic refinement
            select case(kref_intent)

                case(100)
                    do i = 1,8
                        select case(i)
                        
                            case(1,4,5,8)
                                el_pmap(i) = pref_intent(1)
                            case(2,3,6,7)
                                el_pmap(i) = pref_intent(2)

                        end select
                    enddo
                case(010)

                    do i = 1,8
                        select case(i)
                        
                            case(1,2,5,6)
                                el_pmap(i) = pref_intent(1)
                            case(3,4,7,8)
                                el_pmap(i) = pref_intent(2)

                        end select
                    enddo

                case(001)

                    do i = 1,8
                        select case(i)
                        
                            case(1,2,3,4)
                                el_pmap(i) = pref_intent(1)
                            case(5,6,7,8)
                                el_pmap(i) = pref_intent(2)

                        end select
                    enddo

                case(110)

                    do i = 1,8
                        select case(i)
                        
                            case(1,5)
                                el_pmap(i) = pref_intent(1)
                            case(2,6)
                                el_pmap(i) = pref_intent(2)
                            case(3,7)
                                el_pmap(i) = pref_intent(3)
                            case(4,8)
                                el_pmap(i) = pref_intent(4)

                        end select
                    enddo

                case(101)

                    do i = 1,8

                        select case(i)
                        
                            case(1,4)
                                el_pmap(i) = pref_intent(1)
                            case(2,3)
                                el_pmap(i) = pref_intent(2)
                            case(5,8)
                                el_pmap(i) = pref_intent(3)
                            case(6,7)
                                el_pmap(i) = pref_intent(4)

                        end select


                    enddo

                case(011)

                    do i = 1,8

                        select case(i)
                        
                            case(1,2)
                                el_pmap(i) = pref_intent(1)
                            case(3,4)
                                el_pmap(i) = pref_intent(2)
                            case(5,6)
                                el_pmap(i) = pref_intent(3)
                            case(7,8)
                                el_pmap(i) = pref_intent(4)

                        end select


                    enddo

                case(111)

                    do i = 1,8

                        el_pmap(i) = pref_intent(i)

                    enddo

            end select


            !transferring the p-order to kref_close childs


            select case(kref_appl)

                case(100)
                     
                    pref_appl(1) = (el_pmap(1)+el_pmap(4)+el_pmap(5)+el_pmap(8))/4
                    pref_appl(2) = (el_pmap(2)+el_pmap(3)+el_pmap(6)+el_pmap(7))/4


                case(010)

                    pref_appl(1) = (el_pmap(1)+el_pmap(2)+el_pmap(5)+el_pmap(6))/4
                    pref_appl(2) = (el_pmap(3)+el_pmap(4)+el_pmap(7)+el_pmap(8))/4

                case(001)

                    pref_appl(1) = (el_pmap(1)+el_pmap(2)+el_pmap(3)+el_pmap(4))/4
                    pref_appl(2) = (el_pmap(5)+el_pmap(6)+el_pmap(7)+el_pmap(8))/4

                case(110)

                    pref_appl(1) = (el_pmap(1) + el_pmap(5))/2
                    pref_appl(2) = (el_pmap(2) + el_pmap(6))/2
                    pref_appl(3) = (el_pmap(3) + el_pmap(7))/2
                    pref_appl(4) = (el_pmap(4) + el_pmap(8))/2

                case(101)


                    pref_appl(1) = (el_pmap(1) + el_pmap(4))/2
                    pref_appl(2) = (el_pmap(2) + el_pmap(3))/2
                    pref_appl(3) = (el_pmap(5) + el_pmap(8))/2
                    pref_appl(4) = (el_pmap(6) + el_pmap(7))/2

                case(011)


                    pref_appl(1) = (el_pmap(1) + el_pmap(2))/2
                    pref_appl(2) = (el_pmap(3) + el_pmap(4))/2
                    pref_appl(3) = (el_pmap(5) + el_pmap(6))/2
                    pref_appl(4) = (el_pmap(7) + el_pmap(8))/2


                case(111)

                    do i = 1,8
                        pref_appl(i) = el_pmap(i)
                    enddo

            end select
        

    endif
    




end subroutine subson_one_irregularity_map