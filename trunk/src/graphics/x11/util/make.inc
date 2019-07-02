#
# 3D graphics utilities 
#
#
SRC0 += \
	$(SRC_DIR)/graph_util_3d/module/graphmod.F
#
graph_util_3d = $(SRC_DIR)/graph_util_3d/arctan.F \
	$(SRC_DIR)/graph_util_3d/cartobs.F \
	$(SRC_DIR)/graph_util_3d/clip2.F \
	$(SRC_DIR)/graph_util_3d/colorset.F \
	$(SRC_DIR)/graph_util_3d/draw_line.F \
	$(SRC_DIR)/graph_util_3d/dpborder.F \
	$(SRC_DIR)/graph_util_3d/dpvisid.F \
	$(SRC_DIR)/graph_util_3d/finstr.F \
	$(SRC_DIR)/graph_util_3d/initwin.F \
	$(SRC_DIR)/graph_util_3d/mat3prod.F \
	$(SRC_DIR)/graph_util_3d/selwin.F \
	$(SRC_DIR)/graph_util_3d/sortz.F \
	$(SRC_DIR)/graph_util_3d/trobs.F

