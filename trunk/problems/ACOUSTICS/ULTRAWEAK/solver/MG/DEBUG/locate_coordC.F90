!----------------------------------------------------------------------
!
!   routine name       - locate_coordC
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2018
!
!   purpose            - routine calculates coordinates of an element
!                        vertices wrt its coarse mesh ancestor ;
!
!   arguments :
!     in:
!             MdleC    - coarse grid ancestor
!             Mdle     - an element middle node, identified with
!                        the element
!     out:
!             Xsubel   - coordinates of the element
!                        vertices wrt to the ancestor
!
!----------------------------------------------------------------------
!
   subroutine locate_coordC(MdleC,Mdle, Xsubel)
!
   use parameters
   use refinements
   use element_data
   use data_structure3D
!   
   IMPLICIT NONE
!   
   integer, intent(in ) :: MdleC, mdle
   real*8,  intent(out) :: Xsubel(3,8)
!
!..locals
   real*8  ::  x(3,8), x_new(3,8)
!..history information: father, its type, refinement kind, son number
   character(len=4) :: father_type(MAXGEN),son_type(MAXGEN),ftype
   integer :: nfather(MAXGEN), nfather_ref_kind(MAXGEN), noson(MAXGEN)
   integer :: nfath, nson, igen, nrgen, j, jp, is, if, iv,ie, nrsons
   integer :: kref1, kref2, kref3, nv1, nv2, nv3, nv4
!
!----------------------------------------------------------------------
!
   Xsubel = 0.d0
!
!..reference coordinates of ancestor      
   select case(NODES(mdleC)%type)
   case('mdlp'); x(1:3,1:6) = PRISM_COORD(1:3,1:6)
   case('mdlb'); x(1:3,1:8) = BRICK_COORD(1:3,1:8)
   case('mdln'); x(1:3,1:4) = TETRA_COORD(1:3,1:4)
   case('mdld'); x(1:3,1:5) = PYRAM_COORD(1:3,1:5)
   end select
!
!..go up the tree to find the initial mesh ancestor
   nson = Mdle; nfath = NODES(nson)%father;   igen=0
   do 
      if (nson.eq.MdleC) exit
      igen = igen + 1
      nfather(igen) = nfath
      father_type(igen) = NODES(nfath)%type
      son_type(igen) = NODES(nson)%type
      nfather_ref_kind(igen) = NODES(nfath)%ref_kind
      call nr_mdle_sons(NODES(nfath)%type, NODES(nfath)%ref_kind,nrsons)
      call locate(nson,NODES(nfath)%sons,nrsons, noson(igen))
!        
!  ...update        
      nson = nfath
      nfath = NODES(nson)%father
   enddo
   nrgen = igen
!
!..go down the tree reconstructing reference coordinates
!..loop through the generations
   do igen=nrgen,1,-1
      ftype = father_type(igen)
!
!  ...loop through son's vertex nodes
      do j=1,Nvert(son_type(igen))
!
!     ...decode father refinement flag
         call decode_ref(ftype,nfather_ref_kind(igen),kref1,kref2,kref3)
!
!     ...find the parent node of the vertex
         jp = Npar_ref(ftype,j,noson(igen),kref1,kref2,kref3)
!
!     ...find the son number of the vertex
         is = Nson_ref(ftype,j,noson(igen),kref1,kref2,kref3)
!
!     ...node shared with the father
         if (is.eq.0) then
            if (jp.eq.0) jp=j
            x_new(1:3,j) = x(1:3,jp)
!
!     ...node generated through refinement
         else
!
!        ...parent edge node
            if (jp.le.nvert(ftype)+nedge(ftype)) then
               ie = jp-nvert(ftype)
               call edge_to_vert(ftype,ie, nv1,nv2)
               x_new(1:3,j) = (x(1:3,nv1) + x(1:3,nv2))/2.d0
!              
!     ...parent rectangular face node
            elseif (jp.le.nvert(ftype)+nedge(ftype)+nface(ftype)) then
               if = jp-nvert(ftype)-nedge(ftype)
               call face_to_vert(ftype,if, nv1,nv2,nv3,nv4)
               x_new(1:3,j) = (x(1:3,nv1) + x(1:3,nv2)          &
                            +  x(1:3,nv3) + x(1:3,nv4))/4.d0
!
!        ...parent middle node (brick element only)
            else
               x_new(1:3,j) = 0.d0
               do iv=1,8
                  x_new(1:3,j) = x_new(1:3,j) + x(1:3,iv)
               enddo
               x_new(1:3,j) = x_new(1:3,j)/8.d0
!
            endif
         endif
!   
!  ...end of loop through vertex nodes
      enddo
!
!  ...update        
      x = x_new
!
!..end of loop through generations
   enddo
!
   Xsubel(1:3,1:nvert(NODES(Mdle)%type)) = x(1:3,1:nvert(NODES(Mdle)%type))
! 
!
   end subroutine locate_coordC