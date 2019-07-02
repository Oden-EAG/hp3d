!--------------------------------------------------------------------------------------------------------------
subroutine generate_remaining_CUR_TRI_REC
!--------------------------------------------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPORSE: routine generates TWIN CURVES, CONNECTING FIGURES, TWIN FIGURES
!            and updates ORIGINAL FIGURES
!  
! STRUCTURE:
!
!    LOOP through figs to split
!    | LOOP through fig edges
!    | | IF 1st visit to curve
!    | | | 1. generate TWIN CURVE
!    | | | 2. generate CONN FIGURE
!    | | ENDIF
!    | ENDDO
!    | 3. modify original FIGURE
!    | 4. generate TWIN FIGURE
!    ENDDO
!
! CONTAINS: update_conforming_surface
!--------------------------------------------------------------------------------------------------------------
! MODULES
  use SPLIT_SURF
  use GMP
  use control
!--------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------
! VARIABLES
  integer :: ic_new,ifig,nt,nr,nedg,np,np1,iv,iv1,ie,ifound,ns,is,idec
  integer :: status
!--------------------------------------------------------------------------------------------------------------
! FUNCTIONS
  integer :: my_mod
!--------------------------------------------------------------------------------------------------------------
! printint flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'generate_remaining_CUR_TRI_REC: generating curves, triangles, rectangles'
#endif
! ..allocate new_curves and initialize to 0
    call allocate_NEW_CURVES
    new_curves = 0
! ..loop through figures to be split      
    ic_new = 0
    do ifig = 1, nr_figs_to_split
! ....select figure type      
      select case(figs_split(1,ifig))
        case(1);  nt = figs_split(2,ifig);  nedg = 3
        case(2);  nr = figs_split(2,ifig);  nedg = 4
      end select
! ***********************************  TWIN CURVES & CONNECTING FIGURES  ************************************  |
! ....loop through figure edges
      do ie = 1, nedg
! ......local edge endpoints   
        iv = ie;  iv1 = my_mod(iv + 1,nedg)
        select case(figs_split(1,ifig))
          case(1)
            np  = TRIANGLES(nt)%VertNo(iv)
            np1 = TRIANGLES(nt)%VertNo(iv1)
          case(2)
            np  = RECTANGLES(nr)%VertNo(iv)
            np1 = RECTANGLES(nr)%VertNo(iv1)
        end select
! ......if both vertices are on the bounding surface don't duplicate curve
        if ((new_point(np) .eq. np) .and. (new_point(np1) .eq. np1))  cycle  ! <---------------- (*)
! ......determine if curve has already been visited (curves are visited twice)
        call locate_curve(np,np1,new_curves,ic_new, ifound)
! ......if 1st visit to curve
        if (ifound .eq. 0) then
          NRCURVE = NRCURVE + 1
          if (NRCURVE .gt. MAXNC) then
            write(*,*)'generate_remaining_CUR_TRI_REC: increase MAXNC.'
            stop
          endif
          CURVES(NRCURVE)%Type = 'Seglin'                                    ! <-- TWIN CURVE
! ........endpoints of new curve are twin points of old curve            
          CURVES(NRCURVE)%EndPoNo(1) = new_point(np)                         ! <-- TWIN CURVE
          CURVES(NRCURVE)%EndPoNo(2) = new_point(np1)                        ! <-- TWIN CURVE
          ic_new = ic_new + 1
          if (ic_new .gt. MAX_NEW_CURVES) then
            write(*,*)'generate_remaining_CUR_TRI_REC: increase MAX_NEW_CURVES.'
            stop
          endif
! ........update list of new curves endpoints            
          new_curves(1,ic_new) = np
          new_curves(2,ic_new) = np1
!==============================================================================================================
!  REMARK: recall that we are looping over curves that need to be duplicated, see (*) above, hence at least   |
!    one curve endpoint is not on the bounding surface!                                                       |
!==============================================================================================================   
! ........if 1st endpoint lays on bounding surface (hence 2nd doesn't)
          if (new_point(np) .eq. np) then
! ..........generate a new triangle
            NRTRIAN = NRTRIAN + 1
            if (NRTRIAN .gt. MAXTR) then
              write(*,*)'generate_remaining_CUR_TRI_REC: increase MAXTR.'
              stop
            endif
            TRIANGLES(NRTRIAN)%Type      = 'PlaneTri'                        ! <-- CONNECTING TRIANGLE
            TRIANGLES(NRTRIAN)%VertNo(1) = np                                ! <-- CONNECTING TRIANGLE
            TRIANGLES(NRTRIAN)%VertNo(2) = np1                               ! <-- CONNECTING TRIANGLE
            TRIANGLES(NRTRIAN)%VertNo(3) = new_point(np1)                    ! <-- CONNECTING TRIANGLE
            call update_conforming_surface(1)                                ! <-- CONNECTING TRIANGLE
          elseif (new_point(np1) .eq. np1) then
! ..........generate a new triangle
            NRTRIAN = NRTRIAN + 1
            if (NRTRIAN .gt. MAXTR) then
              write(*,*)'generate_remaining_CUR_TRI_REC: increase MAXTR.'
              stop
            endif
            TRIANGLES(NRTRIAN)%Type      = 'PlaneTri'                        ! <-- CONNECTING TRIANGLE
            TRIANGLES(NRTRIAN)%VertNo(1) = np                                ! <-- CONNECTING TRIANGLE
            TRIANGLES(NRTRIAN)%VertNo(2) = np1                               ! <-- CONNECTING TRIANGLE
            TRIANGLES(NRTRIAN)%VertNo(3) = new_point(np)                     ! <-- CONNECTING TRIANGLE
            call update_conforming_surface(1)                                ! <-- CONNECTING TRIANGLE
! ........neither endpoint lays on bounding surface, then generate a rectangle
          else
! ..........generate new rectangle
            NRRECTA = NRRECTA + 1
            if (NRRECTA .gt. MAXRE) then
              write(*,*) 'generate_remaining_CUR_TRI_REC: increase MAXRE.'
              stop
            endif
            RECTANGLES(NRRECTA)%Type      = 'BilQua'                         ! <-- CONNECTING RECTANGLE 
            RECTANGLES(NRRECTA)%VertNo(1) = np                               ! <-- CONNECTING RECTANGLE
            RECTANGLES(NRRECTA)%VertNo(2) = np1                              ! <-- CONNECTING RECTANGLE
            RECTANGLES(NRRECTA)%VertNo(3) = new_point(np1)                   ! <-- CONNECTING RECTANGLE
            RECTANGLES(NRRECTA)%VertNo(4) = new_point(np)                    ! <-- CONNECTING RECTANGLE
            call update_conforming_surface(2)                                ! <-- CONNECTING RECTANGLE
         endif
! ......end if 1st visit to curve
        endif
! ....end of loop through figure edges          
      enddo
! ***********************************************  TWIN FIGURE *********************************************** |
! ....select figure type        
      select case(figs_split(1,ifig))
! ......figure is a triangle        
        case(1)
! ........add triangle          
          NRTRIAN = NRTRIAN + 1
          if (NRTRIAN .gt. MAXTR) then
            write(*,*)'generate_remaining_CUR_TRI_REC: increase MAXTR.'
            stop
          endif
! ........duplicate vertices          
          do iv = 1, 3
            np = TRIANGLES(nt)%VertNo(iv)
            TRIANGLES(NRTRIAN)%VertNo(iv) = new_point(np)                    ! <-- TWIN TRIANGLE
          enddo
! ........determine whether the old and new triangle need to conform to plane on POS or NEG side
          idec = 0
          do iv = 1, 3
            np = TRIANGLES(nt)%VertNo(iv)
! ..........if np is different from its twin point idec++           
            if (new_point(np) .ne. np) idec = idec + 1
          enddo
!===========================================================================
!  SELECT number of duplicated vertices:
!    0,1,2  detach original triangle from plane, if Dh1 > 0
!        3  attach original triangle to plane on NEG side
!==========================================================================
          select case(idec)
! ..........up to 2 vertices have been duplicated
            case(0,1,2)
! ............detach old triangle from the original surface if dh1>0 and not a plane triangle (presence of a plane triangle may result from successive splittings...)
!!!              if ((Dh1_MOD .gt. 0.d0) .and. (TRIANGLES(nt)%Type .ne. 'PlaneTri')) then  
              if (Dh1_MOD .gt. 0.d0) then
                deallocate(TRIANGLES(nt)%Idata, STAT = status)
                if (status .ne. 0) then
                  write(*,*)'generate_remaining_CUR_TRI_REC: Idata not deallocated for nt = ',nt
                  stop
                endif  
                TRIANGLES(nt)%Type = 'PlaneTri'                              ! <-- ORIG TRIANGLE
              endif 
              TRIANGLES(NRTRIAN)%Type = 'PlaneTri'                           ! <-- TWIN TRIANGLE
! ..........all 3 vertices have been duplicated            
            case(3)
! ............attach old triangle to the plane on NEG side
              if (Dh1_MOD .gt. 0.d0)  TRIANGLES(nt)%Idata(1) = NRSURFS - 1   ! <-- ORIG TRIANGLE
              TRIANGLES(NRTRIAN)%Type = 'PTITri'                             ! <-- TWIN TRIANGLE
              allocate(TRIANGLES(NRTRIAN)%Idata(1))                
! ............attach new triangle to the plane on POS side   
              TRIANGLES(NRTRIAN)%Idata(1) = NRSURFS                          ! <-- TWIN TRIANGLE
! ........end select number of duplicated vertices              
          end select
! ......figure is a rectangle          
        case(2)
! ........add rectangle          
          NRRECTA = NRRECTA + 1
          if (NRRECTA .gt. MAXRE) then
            write(*,*)'generate_remaining_CUR_TRI_REC: increase MAXRE.'
            stop
          endif
! ........duplicate vertices 
          do iv = 1, 4
            np = RECTANGLES(nr)%VertNo(iv)
            RECTANGLES(NRRECTA)%VertNo(iv) = new_point(np)                   ! <-- TWIN RECTANGLE
          enddo
! ........determine whether the old and new rectangles need to conform to new planes
          idec = 0
          do iv = 1, 4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np) .ne. np)  idec = idec + 1
          enddo
! ........select number of duplicated vertices         
          select case(idec)
! ..........up to 3 duplicated vertices          
            case(0,1,2,3)
! ............detach old rectangle from the original surface if dh1>0
              if (Dh1_MOD .gt. 0.d0) then   
                deallocate(RECTANGLES(nr)%Idata)  
                RECTANGLES(nr)%Type = 'BilQua'                               ! <-- ORIG RECTANGLE
              endif 
              RECTANGLES(NRRECTA)%Type = 'BilQua'                            ! <-- TWIN RECTANGLE
! ..........all vertices have been duplicated              
            case(4)
              if (Dh1_MOD .gt. 0.d0) RECTANGLES(nr)%Idata(1) = NRSURFS - 1   ! <-- ORIG RECTANGLE
              RECTANGLES(NRRECTA)%Type = 'PTIRec'                            ! <-- TWIN RECTANGLE
              allocate(RECTANGLES(NRRECTA)%Idata(1), STAT = status)
              if (status .ne. 0) then
                write(*,*)'generate_remaining_CUR_TRI_REC: Idata not allocated for nr = ',nr
                stop
              endif
              RECTANGLES(NRRECTA)%Idata(1) = NRSURFS                         ! <-- TWIN RECTANGLE
          end select
#if I_PRINT >= 2          
            write(*,7024) nr,NRRECTA
 7024       format('generate_remaining_CUR_TRI_REC: have duplicated rectangle ',i5,' ; twin rectangle = ',i5)
!!!            call print_GMP
#endif
! ....end select figure type          
      end select
! ..end of loop through figures to split
    enddo
#if I_PRINT >= 1
    write(*,*)'generate_remaining_CUR_TRI_REC: done!'   
#endif    
!
    contains
!
!
!
!------------------------------------------------------------------------------------------------------
subroutine update_conforming_surface(fig_type)
!------------------------------------------------------------------------------------------------------
  use control
  use SPLIT_SURF
!------------------------------------------------------------------------------------------------------  
! DUMMY ARGUMENTS
  integer, intent(in)    :: fig_type
!------------------------------------------------------------------------------------------------------  
! VARIABLES
  real*8               :: fval
  real*8, dimension(3) :: xp,dfdx
!------------------------------------------------------------------------------------------------------  
! printing flag (0,1)
#define I_PRINT 0
!
! ..select figure type
    select case(fig_type)
! **********************************************  TRIANGLE  ***************************************** |    
      case(1)
! ......loop through surfaces
        do is = 1, Nr_confm_MOD
          ns = Ns_confm_MOD(is)
          idec = 0
          do iv = 1, 3
            np = TRIANGLES(NRTRIAN)%VertNo(iv)
            xp(1:3) = POINTS(np)%Rdata(1:3)
            call surf(ns,xp, fval,dfdx)
            if (abs(fval) .lt. GEOM_TOL)  idec = idec + 1
          enddo
! ........if a surface was found         
          if (idec .eq. 3) then
            TRIANGLES(NRTRIAN)%Type = 'PTITri'                          ! <-- CONNECTING TRIANGLE
            allocate(TRIANGLES(NRTRIAN)%Idata(1), STAT = status)
            if (status .ne. 0) then
              write(*,*)'update_conforming_surface: Idata not allocated for nt = ',NRTRIAN
            endif        
            TRIANGLES(NRTRIAN)%Idata(1) = ns                            ! <-- CONNECTING TRIANGLE
#if I_PRINT >= 1        
            write(*,1) NRTRIAN,ns
1           format('update_conforming_surface: attaching triangle ',I5,' to surface ', I3)
#endif
! ..........since triangle can conform only to 1 surface, exit when done
            exit
          endif
! ......end of loop through surfaces to conform to
        enddo   
! **********************************************  RECTANGLE  **************************************** |    
      case(2)
! ......loop through surfaces to conform to
        do is = 1, Nr_confm_MOD
          ns = Ns_confm_MOD(is)
          idec = 0
          do iv = 1, 4
            np = RECTANGLES(NRRECTA)%VertNo(iv)
            xp(1:3) = POINTS(np)%Rdata(1:3)
            call surf(ns,xp, fval,dfdx)
            if (abs(fval) .lt. GEOM_TOL)  idec = idec + 1
          enddo
! ........if a surface was found         
          if (idec .eq. 4) then
            RECTANGLES(NRRECTA)%Type = 'PTIRec'                          ! <-- CONNECTING RECTANGLE
            allocate(RECTANGLES(NRRECTA)%Idata(1), STAT = status)
            if (status .ne. 0) then
              write(*,*)'update_conforming_surface: Idata not allocated for nr = ',NRRECTA
            endif        
            RECTANGLES(NRRECTA)%Idata(1) = ns                            ! <-- CONNECTING RECTANGLE
#if I_PRINT >= 1            
            write(*,2) NRRECTA,ns
2           format('update_conforming_surface: attaching rectangle ',I5,' to surface ', I3)
#endif
! ..........since rectangle can conform only to 1 surface, exit when done
            exit
          endif
! ......end of loop through surfaces to conform to  
        enddo   
      case default
        write(*,*)'update_conforming_surface: unkwon figure type.'
        stop        
    end select
! ..end select figure type  
!
end subroutine update_conforming_surface
!------------------------------------------------------------------------------------------------------  
!
!
!
end subroutine generate_remaining_CUR_TRI_REC
!------------------------------------------------------------------------------------------------------  
