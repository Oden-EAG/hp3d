
      nint = 0
      !     loop over faces
      do jf = 1,nrf
        mdlf = Nfaces(jf)
        call face_vert_list(mdlf,nrv_f,nverts_f)
        allocate(xverts_f(3,nrv_f),xverts_f_a(3,nrv_f))
        do jv=1,Nrv_f
          xverts_f(1:3,jv) = NODES(nverts_f(jv))%coord(1:3,1)
        enddo
        ! transform vertices coordinates from physical to element's
        call poly_x_to_affine_3d(xverts_f,nrv_f,a_aff,xverts_f_a)
        ! compute normal in affine coordinates
        call face_normal(xverts_f_a(:,1:3),fn)
        x_aux(:) = (xverts_f_a(:,1)+xverts_f_a(:,2)+xverts_f_a(:,3))/3.d0
        ! compute distance from origin to plane (dotp)
        call scalar_product(x_aux,fn,dotp)
! 
        if (abs(dotp).le.1e-6) cycle
        !
        ! get face affine coordinates points b0,b1,b2
        call poly_affine_2d(xverts_f_a,nrv_f,b0,b1,b2)
        b_aff(:,1) = b0
        b_aff(:,2) = b1
        b_aff(:,3) = b2
        ! 
        !       set up face order in usual structure
        nordf = 0
        nordf(1:4) = NODES(Mdlf)%order
        norie = 0
        !  .....set up the face quadrature
        allocate(tloc(2,(nrv_f-2)*MAXNINT2ADD),wt((nrv_f-2)*MAXNINT2ADD))
        tloc = 0.d0; wt = 0.d0
        ! get quadrature points for polygon using the affine coordinates
        INTEGRATION = max(0, nord_add_local - 1)
        call set_quadrature_affine_face_DPG(mdlf,NODES(Mdlf)%order,nrv_f,xverts_f,b_aff,nint_f,tloc,wt)
        INTEGRATION = 0
        ! transform integration points' face affine coordinates to element affine coordinates
        xiloc_f = 0.d0
        ! dxi_eta = 0.d0
        call poly_affine_2d_to_x(tloc,nint_f,b_aff,xiloc(:,nint+1:nint+nint_f))
        ! 
!       check if face normal locally goes outward (Norientf = 0). Also correct sign of integral (dir)
        dir = 1.d0
        if(Norientf(jf).ne.0) then
        !       if normal fn goes inward, correct it
          ! fn = -1.d0*fn
          dir = -1.d0
        endif
! 
        wxi(nint+1:nint+nint_f) = wt(1:nint_f) * dotp * dir
! 
        nint = nint + nint_f
! 
        deallocate(xverts_f,xverts_f_a,tloc,wt)
      enddo
 

