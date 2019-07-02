
# OPTIONAL USER-DEFINED INPUT FILE (IT CAN BE USED TO
# SPECIFY DIFFERENT PATHS FOR VTK_LIB_DIR AND VTK_INC_DIR
# ON YOUR LAPTOP)
-include m_options

# DEFAULT FLAG FOR OPTIMIZED OR DEBUG VERSION
ifndef OPT
OPT = 1
endif

# DEFAULT FLAG FOR C++ COMPILER
CXX = g++

# COMPILER FLAGS
FLAGS = -DIPRINT=0 -Wno-deprecated
OPT_FLAGS = -O3
DBG_FLAGS = -g
ifeq ($(OPT),1)
  FLAGS += $(OPT_FLAGS)
  EXEC = vis
else
  FLAGS += $(DBG_FLAGS)
  EXEC = vis_debug
endif

# SHARED LIBRARIES
ifndef VTK_LIB_DIR
VTK_LIB_DIR = /opt/local/lib/vtk-5.6
endif
LIB = $(VTK_LIB_DIR)/libvtkCommon.dylib \
      $(VTK_LIB_DIR)/libvtkDICOMParser.dylib \
      $(VTK_LIB_DIR)/libvtkexoIIc.dylib \
      $(VTK_LIB_DIR)/libvtkexpat.dylib \
      $(VTK_LIB_DIR)/libvtkFiltering.dylib \
      $(VTK_LIB_DIR)/libvtkfreetype.dylib \
      $(VTK_LIB_DIR)/libvtkftgl.dylib \
      $(VTK_LIB_DIR)/libvtkGenericFiltering.dylib \
      $(VTK_LIB_DIR)/libvtkGraphics.dylib \
      $(VTK_LIB_DIR)/libvtkHybrid.dylib \
      $(VTK_LIB_DIR)/libvtkImaging.dylib \
      $(VTK_LIB_DIR)/libvtkInfovis.dylib \
      $(VTK_LIB_DIR)/libvtkIO.dylib \
      $(VTK_LIB_DIR)/libvtkjpeg.dylib \
      $(VTK_LIB_DIR)/libvtklibxml2.dylib \
      $(VTK_LIB_DIR)/libvtkmetaio.dylib \
      $(VTK_LIB_DIR)/libvtkNetCDF.dylib \
      $(VTK_LIB_DIR)/libvtkpng.dylib \
      $(VTK_LIB_DIR)/libvtkRendering.dylib \
      $(VTK_LIB_DIR)/libvtksqlite.dylib \
      $(VTK_LIB_DIR)/libvtksys.dylib \
      $(VTK_LIB_DIR)/libvtktiff.dylib \
      $(VTK_LIB_DIR)/libvtkverdict.dylib \
      $(VTK_LIB_DIR)/libvtkViews.dylib \
      $(VTK_LIB_DIR)/libvtkVolumeRendering.dylib \
      $(VTK_LIB_DIR)/libvtkWidgets.dylib \
      $(VTK_LIB_DIR)/libvtkzlib.dylib

# INCLUDE PATHS
ifndef VTK_INC_DIR
VTK_INC_DIR = /opt/local/include/vtk-5.6
endif
INC = -I$(VTK_LIB_DIR) \
      -I$(VTK_INC_DIR)

# OBJECT FILE DIRECTORY
ifeq ($(OPT),1)
  OBJ_PATH = _RELEASE
else
  OBJ_PATH = _DEBUG
endif

# LIST OF OBJECT FILES
OBJ = $(OBJ_PATH)/main.o \
      $(OBJ_PATH)/util.o \
      $(OBJ_PATH)/geometry.o \
      $(OBJ_PATH)/geometry_hexa.o \
      $(OBJ_PATH)/geometry_tetra.o \
      $(OBJ_PATH)/geometry_prism.o \
      $(OBJ_PATH)/geometry_pyramid.o \
      $(OBJ_PATH)/solution.o

# LINK (DEFAUlT TARGET)
$(EXEC): $(OBJ_PATH)/.dummy $(OBJ)
	$(CXX) $(FLAGS) -o $(EXEC) $(OBJ) $(LIB)

# OBJECT FILE DEPENDENCIES

$(OBJ_PATH)/main.o: main.cxx \
                    geometry.h \
                    solution.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/util.o: util.cxx \
                    util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/geometry.o: geometry.cxx \
                        geometry.h \
                        util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/geometry_hexa.o: geometry_hexa.cxx \
                             geometry.h \
                             util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/geometry_tetra.o: geometry_tetra.cxx \
                              geometry.h \
                              util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/geometry_prism.o: geometry_prism.cxx \
                              geometry.h \
                              util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/geometry_pyramid.o: geometry_pyramid.cxx \
                                geometry.h \
                                util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/solution.o: solution.cxx \
                        solution.h \
                        geometry.h \
                        util.h
	$(CXX) $(FLAGS) $(INC) -o $@ -c $<

$(OBJ_PATH)/.dummy:
	@if [ -d  $(OBJ_PATH) ] ; then \
		touch $@; \
	else \
		mkdir $(OBJ_PATH); touch $@ ; \
	fi
