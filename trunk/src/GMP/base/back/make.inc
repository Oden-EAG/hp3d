#-----------------------------------------------------------------
#  GMP data structure module                                     |
#-----------------------------------------------------------------
#
SRC0 += \
	$(SRC_DIR)/GMP/module/GMP.F \
	$(SRC_DIR)/GMP/module/PML.F \
	$(SRC_DIR)/GMP/module/SPLIT_SURF.F90 \
	$(SRC_DIR)/GMP/module/U2D.F90
#
GMP =   $(SRC_DIR)/GMP/add_layer2plane.F90 \
        $(SRC_DIR)/GMP/add_prism_layer.F \
        $(SRC_DIR)/GMP/check_jacobian.F \
        $(SRC_DIR)/GMP/cone.F \
        $(SRC_DIR)/GMP/connect.F \
        $(SRC_DIR)/GMP/curve.F \
        $(SRC_DIR)/GMP/cylinder1.F \
        $(SRC_DIR)/GMP/def_impl.F \
        $(SRC_DIR)/GMP/determine_plane.F \
        $(SRC_DIR)/GMP/dip.F \
        $(SRC_DIR)/GMP/ellipsoid.F \
        $(SRC_DIR)/GMP/gen_Cyl_Geod.F \
        $(SRC_DIR)/GMP/GMP_kernels.F \
        $(SRC_DIR)/GMP/find_curve_to_block.F \
        $(SRC_DIR)/GMP/find_point_to_block.F \
	$(SRC_DIR)/GMP/Hbasis.F \
        $(SRC_DIR)/GMP/hexa.F \
        $(SRC_DIR)/GMP/input_geometry.F \
        $(SRC_DIR)/GMP/input_GMPdata.F \
        $(SRC_DIR)/GMP/input_NEW.F \
        $(SRC_DIR)/GMP/input_NETGEN.F \
        $(SRC_DIR)/GMP/input_POINTS.F \
        $(SRC_DIR)/GMP/input_CURVES.F \
        $(SRC_DIR)/GMP/input_reconstruct.F \
        $(SRC_DIR)/GMP/input_SURFACES.F \
        $(SRC_DIR)/GMP/local_para.F \
        $(SRC_DIR)/GMP/mnewt.F \
        $(SRC_DIR)/GMP/modify_NETGEN_input.F \
        $(SRC_DIR)/GMP/new_old.F \
	$(SRC_DIR)/GMP/orientations.F \
        $(SRC_DIR)/GMP/plane1.F \
        $(SRC_DIR)/GMP/plane2.F \
        $(SRC_DIR)/GMP/planeANDline.F \
        $(SRC_DIR)/GMP/pointr.F \
        $(SRC_DIR)/GMP/prism.F90 \
        $(SRC_DIR)/GMP/prism_PML.F \
        $(SRC_DIR)/GMP/pyram.F \
        $(SRC_DIR)/GMP/read_imprec.F \
        $(SRC_DIR)/GMP/reconstruct.F \
        $(SRC_DIR)/GMP/reconstruct_curvature.F \
        $(SRC_DIR)/GMP/reconstruct_curve.F \
        $(SRC_DIR)/GMP/reconstruct_normal.F \
        $(SRC_DIR)/GMP/reconstruct_patch.F \
        $(SRC_DIR)/GMP/reconstruct_sharp_edge_curve.F \
        $(SRC_DIR)/GMP/recta.F \
        $(SRC_DIR)/GMP/rotation.F \
        $(SRC_DIR)/GMP/saruss.F \
        $(SRC_DIR)/GMP/set_sharp_edge_bundle.F \
        $(SRC_DIR)/GMP/set_tangents.F \
        $(SRC_DIR)/GMP/set_bundle.F \
        $(SRC_DIR)/GMP/sphere1.F \
        $(SRC_DIR)/GMP/split_surface.F90 \
        $(SRC_DIR)/GMP/surf.F \
        $(SRC_DIR)/GMP/tetra.F \
        $(SRC_DIR)/GMP/tetra_TraTet.F \
        $(SRC_DIR)/GMP/transform.F \
        $(SRC_DIR)/GMP/trian.F \
        $(SRC_DIR)/GMP/udefsu.F \
        $(SRC_DIR)/GMP/upgrade2double.F90 \
        $(SRC_DIR)/GMP/usrfun.F \
#
# GMP/curves
GMP += \
        $(SRC_DIR)/GMP/curves/circular_segment.F \
        $(SRC_DIR)/GMP/curves/coneANDplane.F \
        $(SRC_DIR)/GMP/curves/cone_geodesic.F \
        $(SRC_DIR)/GMP/curves/curve_CylGeod.F \
        $(SRC_DIR)/GMP/curves/curve_1SurfrCur.F \
        $(SRC_DIR)/GMP/curves/curve_1SurfsCur.F \
        $(SRC_DIR)/GMP/curves/curve_2SurfrCur.F \
        $(SRC_DIR)/GMP/curves/curve_2SurfsCur.F \
        $(SRC_DIR)/GMP/curves/curve_3SurfsCur.F \
        $(SRC_DIR)/GMP/curves/curve.F \
        $(SRC_DIR)/GMP/curves/curve_HermCur.F \
        $(SRC_DIR)/GMP/curves/curve_ImpCir.F \
        $(SRC_DIR)/GMP/curves/curve_QuaCir.F \
        $(SRC_DIR)/GMP/curves/curve_SegCir.F \
        $(SRC_DIR)/GMP/curves/cylinderANDplane.F \
        $(SRC_DIR)/GMP/curves/cylinderANDsphere.F \
        $(SRC_DIR)/GMP/curves/cylinder_geodesic.F \
        $(SRC_DIR)/GMP/curves/sphereANDplane.F
#
# GMP/geometry_customization
#GMP += \
#        $(SRC_DIR)/GMP/geometry_customization/customize_sphere.F90 \
#        $(SRC_DIR)/GMP/geometry_customization/customize_shells.F90 \
#        $(SRC_DIR)/GMP/geometry_customization/customize_head.F90
# 
# GMP/DP_geometry
GMP += \
        $(SRC_DIR)/GMP/DP_geometry/reproject_point.F90 \
        $(SRC_DIR)/GMP/DP_geometry/select_surfaces.F90
# 
# GMP/prisms
GMP += \
        $(SRC_DIR)/GMP/prisms/face_bubble_prism.F90 \
        $(SRC_DIR)/GMP/prisms/face_kernel_prism.F90 \
        $(SRC_DIR)/GMP/prisms/prism_TI.F90 \
        $(SRC_DIR)/GMP/prisms/prism_PML.F \
        $(SRC_DIR)/GMP/prisms/proj_prism2face.F90
# 
# GMP/pyramids
GMP += \
        $(SRC_DIR)/GMP/pyramids/pyram_TI.F90
# 
# GMP/GMP_utilities
GMP += \
        $(SRC_DIR)/GMP/GMP_utilities/cartesian2cylindrical.F \
        $(SRC_DIR)/GMP/GMP_utilities/cartesian2spherical.F \
        $(SRC_DIR)/GMP/GMP_utilities/cart_to_polar.F \
        $(SRC_DIR)/GMP/GMP_utilities/check_orientation.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/combination.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/cone2polar.F \
        $(SRC_DIR)/GMP/GMP_utilities/cross_product.F \
        $(SRC_DIR)/GMP/GMP_utilities/cylindrical2cartesian.F \
        $(SRC_DIR)/GMP/GMP_utilities/Dcartesian2polar.F \
        $(SRC_DIR)/GMP/GMP_utilities/Dcylindrical2cartesian.F \
        $(SRC_DIR)/GMP/GMP_utilities/Dpolar2cone.F \
        $(SRC_DIR)/GMP/GMP_utilities/Dspherical2cartesian.F \
        $(SRC_DIR)/GMP/GMP_utilities/factorial.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/find_GMPtriangle.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/give_surf.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/if_bound.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/lambda.F \
        $(SRC_DIR)/GMP/GMP_utilities/locate_curve.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/mod3.F \
        $(SRC_DIR)/GMP/GMP_utilities/mod4.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/modulus.F90 \
        $(SRC_DIR)/GMP/GMP_utilities/polar2cone.F \
        $(SRC_DIR)/GMP/GMP_utilities/spherical2cartesian.F
#
# GMP/rectangles
GMP += \
        $(SRC_DIR)/GMP/rectangles/quadB.F90 \
        $(SRC_DIR)/GMP/rectangles/proj_quad2edge.F90 \
        $(SRC_DIR)/GMP/rectangles/recta_linear.F \
        $(SRC_DIR)/GMP/rectangles/recta_ConfRec.F \
        $(SRC_DIR)/GMP/rectangles/recta_CylRec.F \
        $(SRC_DIR)/GMP/rectangles/recta_HermRec.F \
        $(SRC_DIR)/GMP/rectangles/recta_ImpRec.F \
        $(SRC_DIR)/GMP/rectangles/recta_TraQua.F 
#
# GMP/surface_splitting
GMP += \
        $(SRC_DIR)/GMP/surface_splitting/allocate_FIGS_SPLIT.F90 \
        $(SRC_DIR)/GMP/surface_splitting/allocate_NEW_POINT.F90 \
        $(SRC_DIR)/GMP/surface_splitting/create_TW_POINT_AND_CONN_CURVE.F90 \
        $(SRC_DIR)/GMP/surface_splitting/duplicate_POINTS_AND_CURVES.F90 \
        $(SRC_DIR)/GMP/surface_splitting/generate_PRI_TET_PYR.F90 \
        $(SRC_DIR)/GMP/surface_splitting/generate_remaining_CUR_TRI_REC.F90 \
        $(SRC_DIR)/GMP/surface_splitting/split_interface.F90 \
        $(SRC_DIR)/GMP/surface_splitting/update_connectivities_ORIG_ENTITIES.F90
#
# GMP/OLD_routines
GMP += \
        $(SRC_DIR)/GMP/OLD_routines/conic_triangle_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/cylindrical_triangle_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/prism_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/prism_TI_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/spherical_triangle_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/split_surface_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/tetra_TraTet_OLD.F \
        $(SRC_DIR)/GMP/OLD_routines/upgrade2double_OLD.F
#
# GMP/triangles	
GMP += \
        $(SRC_DIR)/GMP/triangles/conic_triangle.F \
        $(SRC_DIR)/GMP/triangles/cylindrical_triangle.F \
        $(SRC_DIR)/GMP/triangles/plane_triangle.F \
        $(SRC_DIR)/GMP/triangles/spherical_triangle.F \
        $(SRC_DIR)/GMP/triangles/TI_triangle.F \
        $(SRC_DIR)/GMP/triangles/trian_G1RecTri.F \
        $(SRC_DIR)/GMP/triangles/trian_ImpliTri.F \
        $(SRC_DIR)/GMP/triangles/trian_PaSphTri.F \
        $(SRC_DIR)/GMP/triangles/trian_PTITri.F \
        $(SRC_DIR)/GMP/triangles/trian_QtCirTri.F \
        $(SRC_DIR)/GMP/triangles/trian_SpherTri.F \
        $(SRC_DIR)/GMP/triangles/trian_TransTri.F
