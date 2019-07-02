SRC_MODULE = \
	module/thermo_elast.F90

SRC_CORE = \
         celem.F90 \
         data.F90 \
         dirichlet.F90 \
         elem.F90 \
         exact.F90 \
         initialize.F90 \
         set_initial_mesh.F90 \
         set_mesh_ref_filter.F90 \
         soldis.F90 

SRC_MAIN = \
         main_thermo_elast.F90


SRC_FILE = $(SRC_MODULE) $(SRC_CORE)